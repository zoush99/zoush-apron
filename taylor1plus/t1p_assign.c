/*
   APRON Library / Taylor1+ Domain (beta version)
   Copyright (C) 2009-2011 Khalil Ghorbal

*/


#include <stdlib.h>
#include <stdio.h>

#include "t1p_internal.h"
#include "t1p_representation.h"
#include "t1p_fun.h"
#include "t1p_itv_utils.h"
#include "t1p_meetjoin.h"

#include "itv_linearize.h"

#include "ap_dimension.h"
#include "ap_texpr0.h"
#include "ap_generic.h"

#include "t1p_assign.h"

/* Assign and Substitute */
/*************************/
t1p_t* t1p_assign_linexpr_array(ap_manager_t* man,
		bool destructive,
		t1p_t* a,
		ap_dim_t* tdim, ap_linexpr0_t** lexpr,
		size_t size,
		t1p_t* dest)
{
    CALL();
    t1p_internal_t* pr = t1p_init_from_manager(man, AP_FUNID_ASSIGN_LINEXPR_ARRAY);
    ap_texpr0_t** texpr;
    checked_malloc(texpr, ap_texpr0_t*, size, abort(););
    size_t i = 0;
    for (i=0; i<size;i++) {
	texpr[i] = ap_texpr0_from_linexpr0(lexpr[i]);
    }
    t1p_t* res = t1p_assign_texpr_array(man, destructive, a, tdim, texpr, size, dest);
    for (i=0; i<size;i++) {
	ap_texpr0_free(texpr[i]);
    }
    free(texpr);
    return res;
}

t1p_t* t1p_substitute_linexpr_array(ap_manager_t* man,
		bool destructive,
		t1p_t* org,
		ap_dim_t* tdim, ap_linexpr0_t** lexpr,
		size_t size,
		t1p_t* dest)
{
    CALL();
	t1p_internal_t* pr = t1p_init_from_manager(man, AP_FUNID_SUBSTITUTE_LINEXPR_ARRAY);
	not_implemented();
}

/* evaluate an apron tree expression */
/* TODO: destructive not used */
t1p_t* t1p_assign_texpr_array(ap_manager_t* man,
			      bool destructive,
			      t1p_t* a,
			      ap_dim_t* tdim,
			      ap_texpr0_t** texpr,
			      size_t size,
			      t1p_t* dest)
{
    CALL();
    /*
     * NOTA: si on initialise "pr" avant la copie, le pr contiendra la fonction "copy".
     * en mettant la copie avant, on s'assure que le pr contiendra reellement le bon funid.
     */
    t1p_internal_t* pr = t1p_init_from_manager(man, AP_FUNID_ASSIGN_TEXPR_ARRAY);
    size_t i = 0;
    // By zoush99 { @}
    bool real_dim_flag = false;
    real_dim_flag = true;
    ap_interval_t* ap_relative_err = ap_interval_alloc();
    ap_interval_t* ap_absolute_err = ap_interval_alloc();
    ap_interval_set_double(ap_relative_err, 1-pow(2, -23), 1+pow(2, -23));
    ap_interval_set_double(ap_absolute_err, -pow(2, -23), pow(2, -149));

    itv_t itv_relative_err; itv_init(itv_relative_err);
    itv_t itv_absolute_err; itv_init(itv_absolute_err);
    itv_set_ap_interval(pr->itv, itv_relative_err, ap_relative_err);
    itv_set_ap_interval(pr->itv, itv_absolute_err, ap_absolute_err);
    // @}

#ifdef _T1P_DEBUG
    fprintf(stdout, "### ASSIGN TEXPR ARRAY (des %d) %tx ###\n", destructive,(intptr_t)a);
    //t1p_fprint(stdout, man, a, NULL);
    for (i=0; i<size; i++) {
        fprintf(stdout, "(%d) = ", tdim[i]);
        ap_texpr0_fprint(stdout, texpr[i], NULL);
    }
    fprintf(stdout, "\n### ### ###\n");
#endif
    t1p_t* res = t1p_copy(man, a);
    for (i=0; i<res->dims; i++) itv_set(res->paf[i]->itv, res->box[i]);
    for (i=0; i<size; i++) {
        t1p_aff_check_free(pr, res->paf[tdim[i]]);
        res->paf[tdim[i]] = t1p_aff_eval_ap_texpr0(pr, texpr[i], a);

        // By zoush99: we need to consider rounding error for real dimensions
        if(real_dim_flag){
            // relative error
            t1p_aff_mul_itv_inplace(pr, res->paf[tdim[i]], itv_relative_err);
            // absolute error
            // construct the absolute error noise symbol (only need one)
            t1p_aff_t* abs_err = t1p_aff_alloc_init(pr);
            t1p_aff_mul_itv(pr, abs_err, itv_absolute_err);
            // add the absolute error noise symbol
            t1p_aff_add_aff(pr, res->paf[tdim[i]],res->paf[tdim[i]], abs_err);

            // \todo By zoush99
        }

        t1p_aff_reduce(pr, res->paf[tdim[i]]);

        if (t1p_aff_is_top(pr, res->paf[tdim[i]])) {
            t1p_aff_check_free(pr, res->paf[tdim[i]]);
            res->paf[tdim[i]] = pr->top;
        } else if (t1p_aff_is_bottom(pr, res->paf[tdim[i]])) {
            t1p_aff_check_free(pr, res->paf[tdim[i]]);
            res->paf[tdim[i]] = pr->bot;
        }

        itv_set(res->box[tdim[i]],res->paf[tdim[i]]->itv);
        res->paf[tdim[i]]->pby++;
    }
    /* TODO: mettre top pour le moment */
    man->result.flag_best = tbool_top;
    man->result.flag_exact = tbool_top;
#ifdef _T1P_DEBUG
    fprintf(stdout, "### RESULT OF ASSIGN TEXPR ARRAY (des %d) [%tx] ###\n", destructive, (intptr_t)res);
    t1p_fprint(stdout, man, res, NULL);
    fprintf(stdout, "### ### ###\n");
#endif

    // By zoush99 { @}
    ap_interval_free(ap_relative_err);
    ap_interval_free(ap_absolute_err);
    itv_clear(itv_relative_err);
    itv_clear(itv_absolute_err);
    // @}

    /* intersect the result with dest */ 
    if (dest != NULL) {
        t1p_t* ress = t1p_meet(man, false, res, dest);
        free(res);
        return ress;
    } else return res; 
}

t1p_t* t1p_substitute_texpr_array(ap_manager_t* man,
		bool destructive,
		t1p_t* org,
		ap_dim_t* tdim, ap_texpr0_t** texpr,
		size_t size,
		t1p_t* dest)
{
    CALL();
    t1p_internal_t* pr = t1p_init_from_manager(man, AP_FUNID_SUBSTITUTE_TEXPR_ARRAY);
    not_implemented();
}
