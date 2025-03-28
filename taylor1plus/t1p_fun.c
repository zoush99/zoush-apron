/*
   APRON Library / Taylor1+ Domain (beta version)
   Copyright (C) 2009-2011 Khalil Ghorbal

*/



#include "itv_linearize.h"

#include "t1p_internal.h"
#include "t1p_representation.h"
#include "t1p_itv_utils.h"

#ifdef _USE_SDP
#include "declarations.h"
#include "vSDP.h"
#endif

#include "ap_texpr0.h"
#include "t1p_fun.h"


	/* pseudo-taken from http://kiwi.emse.fr/POLE/SDA/corr-iset.html */
	/* we initialize p,q with the head of each list :
	 * - if p and q are null pointers, we reach the halt condition;
	 * - if one of p or q is a null pointer, we copy the non null list;
	 * - else, p and q are non null pointers :
	 *   - if nsym_index p = nsym_index q then 
	 *    * add coeff, 
	 *    * copy the index, 
	 *    * increment p and q to p->next and q->next respectively.
	 *   - if nsym_index p < nsymindex q then
	 *    * copy coeff and nsym_index of p,
	 *    * increment p.
	 *   - else :
	 *    * copy coeff and nsym_index of q,
	 *    * increment q.
	 */
t1p_aff_t* t1p_aff_add(t1p_internal_t* pr, t1p_aff_t* exprA, t1p_aff_t* exprB, t1p_t* abs)
{
    itv_t box; itv_init(box);
    itv_t tmp; itv_init(tmp);
    t1p_aff_t* res = t1p_aff_alloc_init(pr);
    t1p_aaterm_t *p, *q, *ptr;
    itv_add(res->c, exprA->c, exprB->c);
    itv_set(box, res->c);
    if (exprA->q || exprB->q) {
	ptr = t1p_aaterm_alloc_init();
	for(p = exprA->q, q = exprB->q; p || q;) {
	    if (p && q) {
		if (p->pnsym->index == q->pnsym->index) {
		    itv_add(ptr->coeff, p->coeff, q->coeff);
		    ptr->pnsym = p->pnsym;
		    p = p->n ;
		    q = q->n ;
		} else if (p->pnsym->index < q->pnsym->index) {
		    itv_set(ptr->coeff, p->coeff);
		    ptr->pnsym = p->pnsym;
		    p = p->n ;
		} else {
		    itv_set(ptr->coeff, q->coeff);
		    ptr->pnsym = q->pnsym;
		    q = q->n ;
		}
	    } else if (p) {
		itv_set(ptr->coeff, p->coeff);
		ptr->pnsym = p->pnsym;
		p = p->n ;
	    } else {
		itv_set(ptr->coeff, q->coeff);
		ptr->pnsym = q->pnsym;
		q = q->n ;
	    }
	    if (itv_is_zero(ptr->coeff)) {
		if (!(p||q)) {
		    /* the last iteration */
		    t1p_aaterm_free(pr, ptr);
		    if (res->end) res->end->n = NULL;
		}
	    } else {
		/* keep this term */
		if (!res->q) res->q = ptr;
		res->end = ptr;
		res->l++;
		t1p_nsymcons_get_gamma(pr, tmp, ptr->pnsym->index, abs);
		itv_mul(pr->itv, tmp, tmp, ptr->coeff);
		itv_add(box, box, tmp);
		if (p||q) {
		    /* continuing */
		    ptr->n = t1p_aaterm_alloc_init();
		    ptr=ptr->n;
		}
	    }
	}
    }
    itv_add(res->itv, exprA->itv, exprB->itv);
    itv_meet(pr->itv, res->itv, res->itv, box);
    itv_clear(box);
    itv_clear(tmp);
    return res;
}

t1p_aff_t* t1p_aff_sub(t1p_internal_t* pr, t1p_aff_t* exprA, t1p_aff_t* exprB, t1p_t* abs)
{
    itv_t box; itv_init(box);
    itv_t tmp; itv_init(tmp);
    t1p_aff_t* res = t1p_aff_alloc_init(pr);
    t1p_aaterm_t *p, *q, *ptr;
    if (!itv_is_eq(exprA->c, exprB->c)) itv_sub(res->c, exprA->c, exprB->c);
    itv_set(box, res->c);
    if (exprA->q || exprB->q) {
	ptr = t1p_aaterm_alloc_init();
	for(p = exprA->q, q = exprB->q; p || q;) {
	    if (p && q) {
		if (p->pnsym->index == q->pnsym->index) {
		    if (!itv_is_eq(p->coeff, q->coeff)) itv_sub(ptr->coeff, p->coeff, q->coeff);
		    ptr->pnsym = p->pnsym;
		    p = p->n ;
		    q = q->n ;
		} else if (p->pnsym->index < q->pnsym->index) {
		    itv_set(ptr->coeff, p->coeff);
		    ptr->pnsym = p->pnsym;
		    p = p->n ;
		} else {
		    itv_neg(ptr->coeff, q->coeff);
		    ptr->pnsym = q->pnsym;
		    q = q->n ;
		}
	    } else if (p) {
		itv_set(ptr->coeff, p->coeff);
		ptr->pnsym = p->pnsym;
		p = p->n ;
	    } else {
		itv_neg(ptr->coeff, q->coeff);
		ptr->pnsym = q->pnsym;
		q = q->n ;
	    }
	    if (itv_is_zero(ptr->coeff)) {
		if (!(p||q)) {
		    /* the last iteration */
		    t1p_aaterm_free(pr, ptr);
		    if (res->end) res->end->n = NULL;
		}
	    } else {
		/* keep this term */
		if (!res->q) res->q = ptr;
		res->end = ptr;
		res->l++;
		t1p_nsymcons_get_gamma(pr, tmp, ptr->pnsym->index, abs);
		itv_mul(pr->itv, tmp, tmp, ptr->coeff);
		itv_add(box, box, tmp);
		if (p||q) {
		    /* continuing */
		    ptr->n = t1p_aaterm_alloc_init();
		    ptr=ptr->n;
		} else {
		    /* the last iteration */
		    ptr->n = NULL;
		}
	    }
	}
    }
    itv_sub(res->itv, exprA->itv, exprB->itv);
    itv_meet(pr->itv, res->itv, res->itv, box);
    itv_clear(box);
    itv_clear(tmp);
    return res;
}

t1p_aff_t* t1p_aff_mul_itv(t1p_internal_t* pr, t1p_aff_t* src, itv_t lambda)
{
    if (itv_is_zero(lambda) || t1p_aff_is_known_to_be_zero(pr, src)) {
	return t1p_aff_alloc_init(pr);
    } else if (t1p_aff_is_bottom(pr, src) || itv_is_bottom(pr->itv, lambda)) {
	return t1p_aff_bottom(pr);
    } else if (t1p_aff_is_top(pr, src) || itv_is_top(lambda)) {
	return t1p_aff_top(pr);
    } else {
	t1p_aff_t* dst = NULL;
	t1p_aaterm_t *p,*q;
	q = NULL;
	dst = t1p_aff_alloc_init(pr);
	itv_mul(pr->itv, dst->c, lambda, src->c);
	if (src->q) {
	    dst->q = q = t1p_aaterm_alloc_init();
	    for (p=src->q; p; p=p->n) {
		/*if (p == src->lastu) {
		    dst->lastu = q;
		}*/
		itv_mul(pr->itv, q->coeff, lambda, p->coeff);
		q->pnsym = p->pnsym;
		if (p->n) {
		    /* continue */
		    q->n = t1p_aaterm_alloc_init();
		    q = q->n;
		} else {
		    /* the last iteration */
		    dst->end = q; 
		}
	    }
	}
	dst->l = src->l;
	itv_mul(pr->itv, dst->itv, src->itv, lambda);
	return dst;
    }
}

void t1p_aff_mul_itv_inplace(t1p_internal_t* pr, t1p_aff_t* exprA, itv_t lambda)
{
    if (itv_is_zero(lambda) || t1p_aff_is_known_to_be_zero(pr, exprA)) {
	t1p_aff_clear(pr, exprA);
    } else if (t1p_aff_is_bottom(pr, exprA) || itv_is_bottom(pr->itv, lambda)) {
	t1p_aff_check_free(pr, exprA);
	exprA = t1p_aff_bottom(pr);
    } else if (t1p_aff_is_top(pr, exprA) || itv_is_top(lambda)) {
	t1p_aff_check_free(pr, exprA);
	exprA = t1p_aff_top(pr);
    } else {
	t1p_aaterm_t *p;
	itv_mul(pr->itv, exprA->c, exprA->c, lambda);
	for (p=exprA->q; p; p=p->n) {
	    itv_mul(pr->itv, p->coeff, p->coeff, lambda);
	}
	itv_mul(pr->itv, exprA->itv, exprA->itv, lambda);
    }
}

// Commented by zoush99: t1p_aff_mul is a function that multiplies two t1p_aff_t structs together.
// Choose different algothrims according to the hypercube flag in the t1p_t struct, which is the flag of whether it is constrained affine sets
t1p_aff_t* t1p_aff_mul(t1p_internal_t* pr, t1p_aff_t* exprA, t1p_aff_t* exprB, t1p_t* env)
{
    t1p_aff_t* res;
    if (t1p_aff_is_known_to_be_zero(pr, exprA) || t1p_aff_is_known_to_be_zero(pr, exprB)) {
	res = t1p_aff_alloc_init(pr);
    } else if (t1p_aff_is_bottom(pr, exprA) || t1p_aff_is_bottom(pr, exprB)) {
	res = t1p_aff_bottom(pr);
    } else if (t1p_aff_is_top(pr, exprA) || t1p_aff_is_top(pr, exprB)) {
	res = t1p_aff_top(pr);
    } else {
	itv_t box; itv_init(box);
	if (env->hypercube) res = t1p_aff_mul_non_constrained(pr, exprA, exprB, env);
	else res = t1p_aff_mul_constrained(pr, exprA, exprB, env);
	itv_mul(pr->itv, box, exprA->itv, exprB->itv);
	itv_meet(pr->itv, res->itv, res->itv, box);
	itv_clear(box);
    }
    return res;
}

t1p_aff_t* t1p_aff_mul_constrained_backup(t1p_internal_t* pr, t1p_aff_t* exprA, t1p_aff_t* exprB, t1p_t* env)
{
    t1p_aff_t* res = t1p_aff_alloc_init(pr);
    t1p_aaterm_t *p, *q, *ptr;
    itv_t itv_nlin; itv_init(itv_nlin);
    int i = 0;
    int SDP_dim = 0;
    unsigned int dim = 0;	/* compute the number of different eps */
    int cond_SDP;		/* use SDP or not ? */

    itv_t* tmp = itv_array_alloc(5);	/* to store temporary variables */

    itv_t midi, devi, mid1, mid2;
    itv_init(midi); itv_init(devi); itv_init(mid1); itv_init(mid2);

    itv_t* gammabis = itv_array_alloc(pr->dim);

    itv_middev(pr->itv, mid1, tmp[2], exprA->itv);
    itv_middev(pr->itv, mid2, tmp[3], exprB->itv);

    itv_mul(pr->itv, tmp[0], mid1, mid2);
    itv_mul(pr->itv, tmp[1], mid1, exprB->c);
    itv_mul(pr->itv, tmp[2], mid2, exprA->c);

    itv_add(res->c, tmp[1], tmp[2]);
    itv_sub(res->c, res->c, tmp[0]);
    itv_set(res->itv, res->c);

    itv_t eps_itv; itv_init(eps_itv);

    /* hash table to sort noise symbols "de facon contigue". 
     * If we have eps1, eps2, eps4, eps5 shared between the two forms, they are seen as eps1, eps2, eps3, eps4 which optimizes the size of the SDP problem.  
     */
    int* hash = (int *)malloc((2+pr->dim)*sizeof(int));

    int* nsym_shared; 
    checked_malloc(nsym_shared,int,pr->dim,abort(););

    int* nsym_shared_hash;
    /* si on accumule, pour chaque expression, les coeff des symboles ind锟絧endants, on rajoute 2 symboles, un pour chaque expression.
       le hash prend en compte ces deux eventuels ajouts.
     */
    checked_malloc(nsym_shared_hash,int,2+pr->dim,abort(););
    ap_dim_t abs_dim = 0;

    if (exprA->q || exprB->q) {
	res->q = ptr = t1p_aaterm_alloc_init();
	res->l++;
	/* linear part */
	for(p = exprA->q, q = exprB->q; p || q;) {
	    dim += 1;	/* computes the number of different noise symbols */
	    if (p && q) {
		if (p->pnsym->index == q->pnsym->index) {
		    itv_mul(pr->itv, tmp[3], p->coeff, mid2);
		    itv_mul(pr->itv, tmp[4], q->coeff, mid1);
		    itv_add(ptr->coeff, tmp[3], tmp[4]);
		    ptr->pnsym = p->pnsym;
		    hash[p->pnsym->index] = dim;
		    nsym_shared[p->pnsym->index] = true;
		    SDP_dim++;
		    nsym_shared_hash[p->pnsym->index] = SDP_dim;
		    if (env->hypercube) itv_set(gammabis[p->pnsym->index], pr->muu);
		    else {
			if (t1p_nsymcons_get_dimpos(pr, &abs_dim, p->pnsym->index, env)) {
			    itv_set_ap_interval(pr->itv, eps_itv, env->gamma[abs_dim]);
			    itv_middev(pr->itv, midi, devi, eps_itv);
			    itv_sub(gammabis[p->pnsym->index], eps_itv, midi);
			} else {
			    itv_set(gammabis[p->pnsym->index], pr->muu);
			}
		    }
		    p = p->n ;
		    q = q->n ;
		} else if (p->pnsym->index < q->pnsym->index) {
		    itv_mul(pr->itv, ptr->coeff, p->coeff, mid2);
		    ptr->pnsym = p->pnsym;
		    hash[p->pnsym->index] = dim;
		    if (env->hypercube) itv_set(gammabis[p->pnsym->index], pr->muu);
		    else {
			if (t1p_nsymcons_get_dimpos(pr, &abs_dim, p->pnsym->index, env)) {
			    itv_set_ap_interval(pr->itv, eps_itv, env->gamma[abs_dim]);
			    itv_middev(pr->itv, midi, devi, eps_itv);
			    itv_sub(gammabis[p->pnsym->index], eps_itv, midi);
			} else {
			    itv_set(gammabis[p->pnsym->index], pr->muu);
			}
		    }
		    p = p->n ;
		} else {
		    itv_mul(pr->itv, ptr->coeff, q->coeff, mid1);
		    ptr->pnsym = q->pnsym;
		    hash[q->pnsym->index] = dim;
		    if (env->hypercube) itv_set(gammabis[q->pnsym->index], pr->muu);
		    else {
			if (t1p_nsymcons_get_dimpos(pr, &abs_dim, q->pnsym->index, env)) {
			    itv_set_ap_interval(pr->itv, eps_itv, env->gamma[abs_dim]);
			    itv_middev(pr->itv, midi, devi, eps_itv);
			    itv_sub(gammabis[q->pnsym->index], eps_itv, midi);
			} else {
			    itv_set(gammabis[q->pnsym->index], pr->muu);
			}
		    }
		    q = q->n ;
		}
	    } else if (p) {
		itv_mul(pr->itv, ptr->coeff, p->coeff, exprB->c);
		ptr->pnsym = p->pnsym;
		hash[p->pnsym->index] = dim;
		if (env->hypercube) itv_set(gammabis[p->pnsym->index], pr->muu);
		else {
		    if (t1p_nsymcons_get_dimpos(pr, &abs_dim, p->pnsym->index, env)) {
			itv_set_ap_interval(pr->itv, eps_itv, env->gamma[abs_dim]);
			itv_middev(pr->itv, midi, devi, eps_itv);
			itv_sub(gammabis[p->pnsym->index], eps_itv, midi);
		    } else {
			itv_set(gammabis[p->pnsym->index], pr->muu);
		    }
		}
		p = p->n ;
	    } else {
		itv_mul(pr->itv, ptr->coeff, q->coeff, exprA->c);
		ptr->pnsym = q->pnsym;
		hash[q->pnsym->index] = dim;
		if (env->hypercube) itv_set(gammabis[q->pnsym->index], pr->muu);
		else {
		    if (t1p_nsymcons_get_dimpos(pr, &abs_dim, q->pnsym->index, env)) {
			itv_set_ap_interval(pr->itv, eps_itv, env->gamma[abs_dim]);
			itv_middev(pr->itv, midi, devi, eps_itv);
			itv_sub(gammabis[q->pnsym->index], eps_itv, midi);
		    } else {
			itv_set(gammabis[q->pnsym->index], pr->muu);
		    }
		}
		q = q->n ;
	    }
	    /* we don't create a new aaterm if p and q are null */
	    if (p||q) {
		ptr->n = t1p_aaterm_alloc_init();
		res->l++;
		ptr=ptr->n;
	    }
	}
	ptr->n = NULL;
	res->end = ptr;
    }
    /* non linear part */
    /* Calling SDP solver (non guaranteed):
     * - compute different noise symboles, say "n", 
     * - we have "n" constraints, each is a diagonal matrix (can also be seen as 2 block matrix) of dimension 2n, containing 1 in (i,i) and (i+n,i+n)  positions,
     * - objectif matrix (in primal form) is a 2 block matrix of dimension 2n. First block is deduced from the non linear (quadric) part of multiplication, we zero all values of the second block,
     * - objective vector (in dual form) is a one vector of dimension "n".
     */
    if (exprA->q && exprB->q) {
	ap_funopt_t option = ap_manager_get_funopt(pr->man, pr->funid);
	cond_SDP = option.algorithm == INT_MAX ? true : false;
	itv_t itv1; itv_init(itv1);
	itv_t itv2; itv_init(itv2);
	if (cond_SDP) { 
	    call_sdp(pr->itv, itv1, exprA->q, exprB->q, dim, hash, true);
	}
	else 
	    square_dep(pr, itv2, exprA->q, exprB->q, hash, dim, gammabis);
	if (exprA == exprB) {
	    if (bound_sgn(itv1->inf) > 0) {bound_set_int(itv1->inf,0);}
	    if (bound_sgn(itv2->inf) > 0) {bound_set_int(itv2->inf,0);}
	}
	if (cond_SDP) itv_set(itv_nlin,itv1);
	else itv_set(itv_nlin,itv2);
	itv_add(res->c, res->c, itv_nlin);
	itv_clear(itv1);
	itv_clear(itv2);
    }
    /* hash may be NULL if pr->index is zero */
    if (hash) {free(hash);free(nsym_shared);free(nsym_shared_hash);}
    itv_array_free(gammabis, pr->dim);
    itv_array_free(tmp, 5);
    itv_clear(mid1);
    itv_clear(mid2);
    itv_clear(midi);
    itv_clear(devi);
    itv_clear(itv_nlin);
    itv_clear(eps_itv);
    return res;
}

t1p_aff_t* t1p_aff_mul_constrained(t1p_internal_t* pr, t1p_aff_t* exprA, t1p_aff_t* exprB, t1p_t* env)
{
    t1p_aff_t* res = t1p_aff_alloc_init(pr);
    t1p_aaterm_t *p, *q, *ptr;
    itv_t itv_nlin; itv_init(itv_nlin);
    int i = 0;
    int SDP_dim = 0;
    unsigned int dim = 0;	/* compute the number of different eps */
    int cond_SDP;		/* use SDP or not ? */

    itv_t* tmp = itv_array_alloc(5);	/* to store temporary variables */

    itv_t midi, devi, mid1, mid2;
    itv_init(midi); itv_init(devi); itv_init(mid1); itv_init(mid2);

    itv_t* gammabis = itv_array_alloc(pr->dim);

    itv_t itv1; itv_init(itv1);
    itv_t itv2; itv_init(itv2);
    t1p_aff_bound(pr, itv1, exprA, env);
    t1p_aff_bound(pr, itv2, exprB, env);
    itv_middev(pr->itv, mid1, tmp[2], itv1);
    itv_middev(pr->itv, mid2, tmp[3], itv2);
    //itv_middev(pr->itv, mid1, tmp[2], exprA->itv);
    //itv_middev(pr->itv, mid2, tmp[3], exprB->itv);

    itv_mul(pr->itv, tmp[0], mid1, mid2);
    itv_mul(pr->itv, tmp[1], mid1, exprB->c);
    itv_mul(pr->itv, tmp[2], mid2, exprA->c);

    itv_add(res->c, tmp[1], tmp[2]);
    itv_sub(res->c, res->c, tmp[0]);
    itv_set(res->itv, res->c);

    itv_clear(itv1);
    itv_clear(itv2);
    itv_t eps_itv; itv_init(eps_itv);

    /* hash table to sort noise symbols "de facon contigue". 
     * If we have eps1, eps2, eps4, eps5 shared between the two forms, they are seen as eps1, eps2, eps3, eps4 which optimizes the size of the SDP problem.  
     */
    int* hash = (int *)malloc((2+pr->dim)*sizeof(int));

    int* nsym_shared; 
    checked_malloc(nsym_shared,int,pr->dim,abort(););

    int* nsym_shared_hash;
    /* si on accumule, pour chaque expression, les coeff des symboles independants, on rajoute 2 symboles, un pour chaque expression.
       le hash prend en compte ces deux eventuels ajouts.
     */
    checked_malloc(nsym_shared_hash,int,2+pr->dim,abort(););
    ap_dim_t abs_dim = 0;

    if (exprA->q || exprB->q) {
	ptr = t1p_aaterm_alloc_init();
	//res->l++;
	/* linear part */
	for(p = exprA->q, q = exprB->q; p || q;) {
	    dim += 1;	/* computes the number of different noise symbols */
	    if (p && q) {
		if (p->pnsym->index == q->pnsym->index) {
		    itv_mul(pr->itv, tmp[3], p->coeff, mid2);
		    itv_mul(pr->itv, tmp[4], q->coeff, mid1);
		    itv_add(ptr->coeff, tmp[3], tmp[4]);
		    ptr->pnsym = p->pnsym;
		    hash[p->pnsym->index] = dim;
		    nsym_shared[p->pnsym->index] = true;
		    SDP_dim++;
		    nsym_shared_hash[p->pnsym->index] = SDP_dim;
		    if (env->hypercube) itv_set(gammabis[p->pnsym->index], pr->muu);
		    else {
			if (t1p_nsymcons_get_dimpos(pr, &abs_dim, p->pnsym->index, env)) {
			    itv_set_ap_interval(pr->itv, eps_itv, env->gamma[abs_dim]);
			    itv_middev(pr->itv, midi, devi, eps_itv);
			    itv_sub(gammabis[p->pnsym->index], eps_itv, midi);
			} else {
			    itv_set(gammabis[p->pnsym->index], pr->muu);
			    itv_set(eps_itv,pr->muu);
			}
		    }
		    p = p->n ;
		    q = q->n ;
		} else if (p->pnsym->index < q->pnsym->index) {
		    itv_mul(pr->itv, ptr->coeff, p->coeff, mid2);
		    ptr->pnsym = p->pnsym;
		    hash[p->pnsym->index] = dim;
		    if (env->hypercube) itv_set(gammabis[p->pnsym->index], pr->muu);
		    else {
			if (t1p_nsymcons_get_dimpos(pr, &abs_dim, p->pnsym->index, env)) {
			    itv_set_ap_interval(pr->itv, eps_itv, env->gamma[abs_dim]);
			    itv_middev(pr->itv, midi, devi, eps_itv);
			    itv_sub(gammabis[p->pnsym->index], eps_itv, midi);
			} else {
			    itv_set(gammabis[p->pnsym->index], pr->muu);
			    itv_set(eps_itv,pr->muu);
			}
		    }
		    p = p->n ;
		} else {
		    itv_mul(pr->itv, ptr->coeff, q->coeff, mid1);
		    ptr->pnsym = q->pnsym;
		    hash[q->pnsym->index] = dim;
		    if (env->hypercube) itv_set(gammabis[q->pnsym->index], pr->muu);
		    else {
			if (t1p_nsymcons_get_dimpos(pr, &abs_dim, q->pnsym->index, env)) {
			    itv_set_ap_interval(pr->itv, eps_itv, env->gamma[abs_dim]);
			    itv_middev(pr->itv, midi, devi, eps_itv);
			    itv_sub(gammabis[q->pnsym->index], eps_itv, midi);
			} else {
			    itv_set(gammabis[q->pnsym->index], pr->muu);
			    itv_set(eps_itv,pr->muu);
			}
		    }
		    q = q->n ;
		}
	    } else if (p) {
		itv_mul(pr->itv, ptr->coeff, p->coeff, mid2);
		ptr->pnsym = p->pnsym;
		hash[p->pnsym->index] = dim;
		if (env->hypercube) itv_set(gammabis[p->pnsym->index], pr->muu);
		else {
		    if (t1p_nsymcons_get_dimpos(pr, &abs_dim, p->pnsym->index, env)) {
			itv_set_ap_interval(pr->itv, eps_itv, env->gamma[abs_dim]);
			itv_middev(pr->itv, midi, devi, eps_itv);
			itv_sub(gammabis[p->pnsym->index], eps_itv, midi);
		    } else {
			itv_set(gammabis[p->pnsym->index], pr->muu);
			itv_set(eps_itv,pr->muu);
		    }
		}
		p = p->n ;
	    } else {
		itv_mul(pr->itv, ptr->coeff, q->coeff, mid1);
		ptr->pnsym = q->pnsym;
		hash[q->pnsym->index] = dim;
		if (env->hypercube) itv_set(gammabis[q->pnsym->index], pr->muu);
		else {
		    if (t1p_nsymcons_get_dimpos(pr, &abs_dim, q->pnsym->index, env)) {
			itv_set_ap_interval(pr->itv, eps_itv, env->gamma[abs_dim]);
			itv_middev(pr->itv, midi, devi, eps_itv);
			itv_sub(gammabis[q->pnsym->index], eps_itv, midi);
		    } else {
			itv_set(gammabis[q->pnsym->index], pr->muu);
			itv_set(eps_itv,pr->muu);
		    }
		}
		q = q->n ;
	    }
	    if (itv_is_zero(ptr->coeff)) {
		if (!(p||q)) {
		    /* the last iteration */
		    t1p_aaterm_free(pr, ptr);
		    if (res->end) res->end->n = NULL;
		}
	    } else {
		/* keep this term */
		if (!res->q) res->q = ptr;
		res->end = ptr;
		res->l++;
		itv_mul(pr->itv, eps_itv, eps_itv, ptr->coeff);
		itv_add(res->itv, res->itv, eps_itv);
		if (p||q) {
		    /* continuing */
		    ptr->n = t1p_aaterm_alloc_init();
		    ptr=ptr->n;
		} else {
		    /* the last iteration */
		    ptr->n = NULL;
		}
	    }
	}
    }
    /* non linear part */
    /* Calling SDP solver (non guaranteed):
     * - compute different noise symboles, say "n", 
     * - we have "n" constraints, each is a diagonal matrix (can also be seen as 2 block matrix) of dimension 2n, containing 1 in (i,i) and (i+n,i+n)  positions,
     * - objectif matrix (in primal form) is a 2 block matrix of dimension 2n. First block is deduced from the non linear (quadric) part of multiplication, we zero all values of the second block,
     * - objective vector (in dual form) is a one vector of dimension "n".
     */
    if (exprA->q && exprB->q) {
	ap_funopt_t option = ap_manager_get_funopt(pr->man, pr->funid);
	cond_SDP = option.algorithm == INT_MAX ? true : false;
	itv_t itv1; itv_init(itv1);
	itv_t itv2; itv_init(itv2);
	if (cond_SDP) { 
	    call_sdp(pr->itv, itv1, exprA->q, exprB->q, dim, hash, true);
	}
	else 
	    square_dep(pr, itv2, exprA->q, exprB->q, hash, dim, gammabis);
	if (exprA == exprB) {
	    if (bound_sgn(itv1->inf) > 0) {bound_set_int(itv1->inf,0);}
	    if (bound_sgn(itv2->inf) > 0) {bound_set_int(itv2->inf,0);}
	}
	if (cond_SDP) itv_set(itv_nlin,itv1);
	else itv_set(itv_nlin,itv2);
	itv_add(res->c, res->c, itv_nlin);
	itv_add(res->itv, res->itv, itv_nlin);
	itv_clear(itv1);
	itv_clear(itv2);
    }
    /* hash may be NULL if pr->index is zero */
    if (hash) {free(hash);free(nsym_shared);free(nsym_shared_hash);}
    itv_array_free(gammabis, pr->dim);
    itv_array_free(tmp, 5);
    itv_clear(mid1);
    itv_clear(mid2);
    itv_clear(midi);
    itv_clear(devi);
    itv_clear(eps_itv);
    itv_clear(itv_nlin);
    return res;
}

t1p_aff_t* t1p_aff_mul_non_constrained(t1p_internal_t* pr, t1p_aff_t* exprA, t1p_aff_t* exprB, t1p_t* env)
{
    t1p_aaterm_t *p, *q, *ptr;
    itv_t itv_nlin; itv_init(itv_nlin);
    int i = 0;
    int SDP_dim = 0;
    unsigned int dim = 0;	/* compute the number of different eps */
    int cond_SDP;		/* use SDP or not ? */
    itv_t tmp; itv_init(tmp);
    /* hash table to sort noise symbols "de facon contigue". 
     * If we have eps1, eps2, eps4, eps5 shared between the two forms, they are seen as eps1, eps2, eps3, eps4 which optimizes the size of the SDP problem.  
     */
    int* hash = (int *)malloc((2+pr->dim)*sizeof(int));

    int* nsym_shared; 
    checked_malloc(nsym_shared,int,pr->dim,abort(););

    int* nsym_shared_hash;
    /* si on accumule, pour chaque expression, les coeff des symboles ind锟絧endants, on rajoute 2 symboles, un pour chaque expression.
       le hash prend en compte ces deux eventuels ajouts.
     */
    checked_malloc(nsym_shared_hash,int,2+pr->dim,abort(););

    t1p_aff_t* res = t1p_aff_alloc_init(pr);
    itv_mul(pr->itv, res->c, exprA->c, exprB->c);
    itv_set(res->itv, res->c);
    if (exprA->q || exprB->q) {
	ptr = t1p_aaterm_alloc_init();
	for(p = exprA->q, q = exprB->q; p || q;) {
	    dim += 1;	/* compute the number of different noise symbols */
	    if (p && q) {
		if (p->pnsym->index == q->pnsym->index) {
		    itv_mul(pr->itv, tmp, p->coeff, exprB->c);
		    itv_mul(pr->itv, ptr->coeff, q->coeff, exprA->c);
		    itv_add(ptr->coeff, ptr->coeff, tmp);
		    ptr->pnsym = p->pnsym;

		    hash[p->pnsym->index] = dim;
		    nsym_shared[p->pnsym->index] = true;
		    SDP_dim++;
		    nsym_shared_hash[p->pnsym->index] = SDP_dim;

		    p = p->n ;
		    q = q->n ;
		} else if (p->pnsym->index < q->pnsym->index) {
		    itv_mul(pr->itv, ptr->coeff, p->coeff, exprB->c);
		    ptr->pnsym = p->pnsym;
		    hash[p->pnsym->index] = dim;
		    p = p->n ;
		} else {
		    itv_mul(pr->itv, ptr->coeff, q->coeff, exprA->c);
		    ptr->pnsym = q->pnsym;
		    hash[q->pnsym->index] = dim;
		    q = q->n ;
		}
	    } else if (p) {
		itv_mul(pr->itv, ptr->coeff, p->coeff, exprB->c);
		ptr->pnsym = p->pnsym;
		hash[p->pnsym->index] = dim;
		p = p->n ;
	    } else {
		itv_mul(pr->itv, ptr->coeff, q->coeff, exprA->c);
		ptr->pnsym = q->pnsym;
		hash[q->pnsym->index] = dim;
		q = q->n ;
	    }
	    if (itv_is_zero(ptr->coeff)) {
		if (!(p||q)) {
		    /* the last iteration */
		    t1p_aaterm_free(pr, ptr);
		    if (res->end) res->end->n = NULL;
		}
	    } else {
		/* keep this term */
		if (!res->q) res->q = ptr;
		res->end = ptr;
		res->l++;
		t1p_nsymcons_get_gamma(pr, tmp, ptr->pnsym->index, env);
		itv_mul(pr->itv, tmp, tmp, ptr->coeff);
		itv_add(res->itv, res->itv, tmp);
		if (p||q) {
		    /* continuing */
		    ptr->n = t1p_aaterm_alloc_init();
		    ptr=ptr->n;
		} else {
		    /* the last iteration */
		    ptr->n = NULL;
		}
	    }
	}
    }
    /* non linear part */
    /* Calling SDP solver (non guaranteed):
     * - compute different noise symboles, say "n", 
     * - we have "n" constraints, each is a diagonal matrix (can also be seen as 2 block matrix) of dimension 2n, containing 1 in (i,i) and (i+n,i+n)  positions,
     * - objectif matrix (in primal form) is a 2 block matrix of dimension 2n. First block is deduced from the non linear (quadric) part of multiplication, we zero all values of the second block,
     * - objective vector (in dual form) is a one vector of dimension "n".
     */
    if (exprA->q && exprB->q) {
	ap_funopt_t option = ap_manager_get_funopt(pr->man, pr->funid);
	cond_SDP = option.algorithm == INT_MAX ? true : false;
	itv_t itv1; itv_init(itv1);
	itv_t itv2; itv_init(itv2);
	if (cond_SDP) { 
	    call_sdp(pr->itv, itv1, exprA->q, exprB->q, dim, hash, true);
	    /* tentative de mise en place d'un solveur garantit, to be continued ... */
	    //		    itv_t* array = itv_array_alloc(dim*(dim+1)/2);
	    //		    buildIntervalUpperTriangle(pr, array, exprA->q, exprB->q, dim, hash);
	    //		    double* phi = buildDoubleUpperTriangle(pr, array, dim);
	    //		    call_sound_sdp(dim, phi);
	    //		    free(phi);
	    //		    itv_array_free(array,(dim*(dim+1)/2));
	}
	else 
	    square_dep(pr, itv2, exprA->q, exprB->q, hash, dim, NULL);

	if (exprA == exprB) {
	    if (bound_sgn(itv1->inf) > 0) {bound_set_int(itv1->inf,0);}
	    if (bound_sgn(itv2->inf) > 0) {bound_set_int(itv2->inf,0);}
	}
	if (cond_SDP) itv_set(itv_nlin,itv1);
	else itv_set(itv_nlin,itv2);
	itv_add(res->c, res->c, itv_nlin);
	itv_add(res->itv, res->itv, itv_nlin);
	itv_clear(itv1);
	itv_clear(itv2);
    }
    /* hash may be NULL if pr->index is zero */
    if (hash) {free(hash);free(nsym_shared);free(nsym_shared_hash);}
    itv_clear(tmp);
    itv_clear(itv_nlin);
    return res;
}

void refine_affine_div_bounds(t1p_internal_t* pr, t1p_aff_t* x, t1p_aff_t* y, itv_t lower_bound, itv_t upper_bound, size_t counter) {
    // 初始化变量
    itv_t itv_temp, itv_mid, w_lower, w_upper;
    itv_t itv_two;
    itv_init(itv_temp);
    itv_init(itv_mid);
    itv_init(w_lower);
    itv_init(w_upper);
    itv_init(itv_two);

    itv_set_int(itv_two, 2); // itv_two = 2

    size_t current_counter = 0;

    // temp = (lower_bound + upper_bound) / 2
    itv_add(itv_mid, lower_bound, upper_bound); // mid = lower_bound + upper_bound
    itv_div(pr->itv, itv_temp, itv_mid, itv_two); // temp = (lower_bound + upper_bound) / 2

    // 开始区间搜索
    while ((itv_is_pos(lower_bound) || itv_is_neg(upper_bound)) && current_counter < counter) {
        // 计算 w = x - temp * y
        t1p_aff_t* temp_y = t1p_aff_mul_itv(pr, y, itv_temp); // temp_y = temp * y
        t1p_aff_t* w = t1p_aff_sub(pr, x, temp_y, NULL);      // w = x - temp * y

        // 获取 w 的上下界
        itv_set(w_lower, w->itv); // w_lower = lower bound of w
        itv_set(w_upper, w->itv); // w_upper = upper bound of w

        // 更新计数器
        current_counter++;

        // 根据 w 的上下界更新 temp 和 lower_bound/upper_bound
        if (itv_is_neg(w_upper)) {
            // 如果 w 的上界小于等于 0，更新 upper_bound
            itv_set(upper_bound, itv_temp); // upper_bound = temp
            itv_add(itv_mid, lower_bound, upper_bound);
            itv_div(pr->itv, itv_temp, itv_mid, itv_two); // temp = (lower_bound + upper_bound) / 2
        } else if (itv_is_pos(w_lower)) {
            // 如果 w 的下界大于等于 0，更新 lower_bound
            itv_set(lower_bound, itv_temp); // lower_bound = temp
            itv_add(itv_mid, lower_bound, upper_bound);
            itv_div(pr->itv, itv_temp, itv_mid, itv_two); // temp = (lower_bound + upper_bound) / 2
        } else {
            // 如果 w 的上下界跨越 0，根据 lower_bound 和 upper_bound 的关系调整 temp
            if (itv_is_neg(lower_bound) && itv_is_pos(upper_bound)) {
                itv_add(itv_mid, upper_bound, itv_temp);
                itv_div(pr->itv, itv_temp, itv_mid, itv_two); // temp = (upper_bound + temp) / 2
            } else {
                itv_add(itv_mid, lower_bound, itv_temp);
                itv_div(pr->itv, itv_temp, itv_mid, itv_two); // temp = (lower_bound + temp) / 2
            }
        }

        // 清理临时变量
        t1p_aff_clear(pr, temp_y);
        t1p_aff_clear(pr, w);
    }

    // 清理资源
    itv_clear(itv_temp);
    itv_clear(itv_mid);
    itv_clear(w_lower);
    itv_clear(w_upper);
    itv_clear(itv_two);
}

/* computes exprA/exprB.
   If 0 \in exprB->itv, return top + caveat
   else if exprB->itv has an infity bound return (1/exprB->itv)*exprA
   else return (1/exprB)*exprA
 */
// \todo By zoush99, I need to write the code here which is my novelty.
t1p_aff_t* t1p_aff_div(t1p_internal_t* pr, t1p_aff_t* exprA, t1p_aff_t* exprB, t1p_t* env)
{
    t1p_aff_t* res;
    itv_t box;
    itv_t mid, dev;
    itv_t a,b;
    itv_t zeta, alpha, beta;
    itv_t one;
    itv_t tmp;

    itv_init(box);
    itv_init(mid); itv_init(dev);
    itv_init(a); itv_init(b);
    itv_init(zeta); itv_init(alpha); itv_init(beta);
    itv_init(one);
    itv_init(tmp);

    itv_set_int(one, (long int)1);
    if (itv_is_leq(tmp, exprB->itv)) {
	/// \todo bugs here!!! By zoush99
	fprintf(stderr,"(caveat) division by zero.\n");
	res = t1p_aff_alloc_init(pr);
	itv_set_top(res->c);
	itv_set_top(res->itv);
    } else if (itv_has_infty_bound(exprB->itv)) {
	itv_div(pr->itv, tmp, one, exprB->itv);
	res = t1p_aff_mul_itv(pr, exprA, tmp);
    } else {
	if (exprB->q == NULL) {
	    itv_div(pr->itv, tmp, one, exprB->c);
	    res = t1p_aff_mul_itv(pr,exprA,tmp);
	} else {
	    if (itv_is_neg(exprB->itv)) itv_neg(box,exprB->itv);
	    else itv_set(box, exprB->itv);

	    itv_set(a,box);bound_neg(a->sup,a->inf);	// a <- [box->inf,box->inf]
	    itv_set(b,box);bound_neg(b->inf,b->sup);	// b <- [box->sup,box->sup]
	    itv_middev(pr->itv, mid, dev, box);

	    /* alpha = -1/mid^2 */
	    itv_mul(pr->itv, alpha, mid, mid);
	    itv_div(pr->itv, alpha, one, alpha);
	    itv_neg(alpha,alpha);
	    /* beta = 1/2a - 1/(a^2/2b + b/2 + a) */
	    itv_mul(pr->itv, tmp, a, a);	// a*a
	    itv_div(pr->itv, tmp, tmp, b);	// a*a/b
	    itv_mul_2exp(beta, tmp, -1);	// a*a/2b
	    itv_mul_2exp(tmp, b, -1);		// b/2
	    itv_add(beta, beta, tmp);
	    itv_add(beta, beta, a);
	    itv_div(pr->itv, beta, one, beta);
	    itv_neg(beta, beta);
	    itv_div(pr->itv, tmp, one, a);	// 1/a
	    itv_mul_2exp(tmp, tmp, -1);		// 1/2a
	    itv_add(beta, beta, tmp);
	    /* zeta = 2/mid + beta  = 1/2a + 1/(b^2/2a + a/2 + b) + 2/(a+b) */
	    itv_mul(pr->itv, tmp, b, b);	// b*b
	    itv_div(pr->itv, tmp, tmp, a);	// b*b/a
	    itv_mul_2exp(zeta, tmp, -1);	// b*b/2a
	    itv_mul_2exp(tmp, a, -1);		// a/2
	    itv_add(zeta, zeta, tmp);
	    itv_add(zeta, zeta, b);
	    itv_div(pr->itv, zeta, one, zeta);
	    itv_div(pr->itv, tmp, one, a);	// 1/a
	    itv_mul_2exp(tmp, tmp, -1);         // 1/2a
	    itv_add(zeta, zeta, tmp);
	    itv_div(pr->itv, tmp, one, mid);	// 2/(a+b)
	    itv_add(zeta, zeta, tmp);

	    t1p_aff_t* inv = t1p_aff_mul_itv(pr, exprB, alpha);
	    if (itv_is_neg(exprB->itv)) itv_sub(inv->c, inv->c, zeta);
	    else itv_add(inv->c, inv->c, zeta);
	    bound_set(beta->inf, beta->sup);
	    itv_add(inv->c, inv->c, beta);
	    //t1p_aff_nsym_create(pr, inv, beta);
	    itv_div(pr->itv, inv->itv, one, exprB->itv);	/* update inv->itv */
	    /* x/y = x*y^{-1} */
	    res = t1p_aff_mul(pr,exprA,inv,env);	/* inplace , boxes are needed */
	    //if (itv_is_neg(exprB->itv)) t1p_aff_neg_inplace(pr,res);
	    t1p_aff_free(pr, inv);
	    itv_div(pr->itv, box, exprA->itv, exprB->itv);
	    itv_meet(pr->itv, res->itv, res->itv, box);
	}
    }

    
	// begin use interval search to get the result
		// 初始化上下界和搜索次数
	itv_t lower_bound, upper_bound;
	itv_init(lower_bound);
	itv_init(upper_bound);
	
	// 设置初始上下界
	itv_set_num(lower_bound, res->itv->inf); // 下界为结果区间的下界
	itv_set_num(upper_bound, res->itv->sup); // 上界为结果区间的上界
	
	size_t max_counter = 10; // 设置最大搜索次数
	
	// 调用 refine_affine_div_bounds 函数进行精化
	refine_affine_div_bounds(pr, exprA, exprB, lower_bound, upper_bound, max_counter);
	
	// 更新结果区间
	itv_set_num2(res->itv, lower_bound->inf, upper_bound->sup);
	
	// \todo By zoush99, have not use the refined interval to update the coefficient of the new noise symbol.
	
	// 清理资源
	itv_clear(lower_bound);
	itv_clear(upper_bound);
	
	itv_clear(box);
    itv_clear(one);
    itv_clear(a); itv_clear(b);
    itv_clear(zeta); itv_clear(alpha); itv_clear(beta);
    itv_clear(mid); itv_clear(dev);
    itv_clear(tmp);
    return res;
}

t1p_aff_t* t1p_aff_mod (t1p_internal_t* pr, t1p_aff_t* exprA, t1p_aff_t* exprB)
{
	not_implemented();
}

/* a = sqrt(b) */
t1p_aff_t* t1p_aff_sqrt( t1p_internal_t* pr, t1p_aff_t* expr, t1p_t* env)
{
    t1p_aff_t* res;
    itv_t zero; itv_init(zero);
    itv_t box; itv_init(box);
    bool eval_aff = false;

    if (t1p_aff_is_bottom(pr, expr)) res = t1p_aff_bottom(pr);
    else if (itv_is_neg(expr->itv)) {
	if (itv_is_leq(zero,expr->itv)) res = t1p_aff_alloc_init(pr);
	else {
	    fprintf(stderr," (error) computing the square root of something negative.\n");
	    res = t1p_aff_bottom(pr);
	}
    } else if (itv_is_pos(expr->itv)) {
	if (t1p_aff_is_zero(pr, expr)) {
	    res = t1p_aff_alloc_init(pr);
	} else {
	    itv_set(box,expr->itv);
	    eval_aff = true;
	}
    } else {
	fprintf(stderr," (caveat) computing the square root of something negative.\n");
	itv_set(box,expr->itv);
	bound_set_int(box->inf, 0);
	eval_aff = true;
    }
    if (eval_aff) {
	itv_t alpha, beta, zeta;
	itv_t mid, dev;
	itv_t one, tmp, tmp1;
	itv_t a, b;

	itv_init(alpha); itv_init(beta); itv_init(zeta);
	itv_init(mid); itv_init(dev);
	itv_init(one); itv_init(tmp); itv_init(tmp1);
	itv_init(a); itv_init(b);

	itv_set_int(one, 1);
	itv_set(a,box);bound_neg(a->sup,a->inf);	// a <- [box->inf,box->inf]
	itv_set(b,box);bound_neg(b->inf,b->sup);	// b <- [box->sup,box->sup]

	itv_middev(pr->itv, mid, dev, box);
	itv_sqrt(pr->itv, tmp, mid);	// sqrt(mid)
	itv_div(pr->itv, alpha, one, tmp);	// 1/sqrt(mid)
	itv_mul_2exp(alpha, alpha, -1);

	itv_mul_2exp(beta, tmp, -1);
	if (!itv_is_zero(a)) {
	    itv_sqrt(pr->itv, tmp1, a);		// sqrt(a)
	    itv_sub(beta, beta, tmp1);
	    itv_mul(pr->itv, tmp1, a, a);	// a*a
	    itv_div(pr->itv, tmp1, b, tmp1);	// b/a*a
	    itv_div(pr->itv, tmp, one, a);	// 1/a
	    itv_add(tmp, tmp1, tmp);		// 1/a + b/a*a
	    itv_mul_2exp(tmp, tmp, 1);		// 2(1/a + b/a*a)
	    itv_sqrt(pr->itv, tmp, tmp);	// sqrt(2(1/a + b/a*a))
	    itv_div(pr->itv, tmp, one, tmp);	// 1/sqrt(2(1/a + b/a*a))
	    itv_add(beta, beta, tmp);
	}
	itv_mul_2exp(beta, beta, -1);

	itv_sqrt(pr->itv, tmp, mid);	// sqrt(mid)
	itv_mul_2exp(zeta, tmp, -1);
	itv_sub(zeta, zeta, beta);

	res = t1p_aff_mul_itv(pr, expr, alpha);
	itv_add(res->c, res->c, zeta);
	bound_set(beta->inf, beta->sup);
	itv_add(res->c, res->c, beta);

	itv_sqrt(pr->itv, res->itv, box);	/* update a->itv */

	itv_clear(a); itv_clear(b);
	itv_clear(alpha); itv_clear(beta); itv_clear(zeta);
	itv_clear(mid); itv_clear(dev);
	itv_clear(one);
	itv_clear(tmp); itv_clear(tmp1);

    }
    itv_clear(zero); itv_clear(box);
    return res;
}

/* a = -b */
t1p_aff_t* t1p_aff_neg(t1p_internal_t* pr, t1p_aff_t* b)
{
    t1p_aff_t* res = t1p_aff_alloc_init(pr);
    t1p_aaterm_t *p, *q;
    itv_neg(res->c, b->c);
    res->l = b->l;
    if (b->q) {
	res->q = t1p_aaterm_alloc_init();

	for (p=b->q, q=res->q; p->n; p=p->n) {
	    itv_neg(q->coeff, p->coeff);
	    q->pnsym = p->pnsym;
	    q->n = t1p_aaterm_alloc_init();
	    q = q->n;
	}
	itv_neg(q->coeff, p->coeff);
	q->pnsym = p->pnsym;
	res->end = q;
    }
    itv_neg(res->itv, b->itv);
    return res;
}

void t1p_aff_neg_inplace(t1p_internal_t* pr, t1p_aff_t* b)
{
    t1p_aaterm_t *p;
    itv_neg(b->c, b->c);
    for (p=b->q; p->n; p=p->n) {
	itv_neg(p->coeff, p->coeff);
    }
    itv_neg(b->itv, b->itv);
}

t1p_aff_t* t1p_aff_eval_node_unary (t1p_internal_t* pr, ap_texpr0_node_t* node, t1p_t* env)
{
    t1p_aff_t* res;
    t1p_aff_t* pB = t1p_aff_eval_ap_texpr0(pr, node->exprA, env);
#ifdef _T1P_DEBUG_FUN
    fprintf(stdout, "### UNARY OPERANDES ###\n");
    t1p_aff_fprint(pr, stdout, pB);
#endif
    switch (node->op) {
	case AP_TEXPR_NEG:
	      {
		res = t1p_aff_neg(pr, pB);
		break;
	      }
	case AP_TEXPR_SQRT:
	      {
		res = t1p_aff_sqrt(pr, pB, env);
		break;
	      }
	case AP_TEXPR_CAST:
	      {
		not_implemented();
		break;
	      }
	default:
	    fatal("Unknown unary operation");
    }
    /* if pB->pby is non zero pB pointes to an affine form of a dimension, we shall not free it here */
    //if (node->exprA->discr != AP_TEXPR_DIM) t1p_aff_free(pr, pB);
    t1p_aff_check_free(pr, pB);
    return res;
}

// \todo By zoush99, I need to write the code here which is my novelty.
t1p_aff_t* t1p_aff_eval_node_binary (t1p_internal_t* pr, ap_texpr0_node_t* node, t1p_t* env)
{
    t1p_aff_t* res;
    t1p_aff_t* exprAB[2];
    exprAB[0] = t1p_aff_eval_ap_texpr0(pr, node->exprA, env);
    exprAB[1] = t1p_aff_eval_ap_texpr0(pr, node->exprB, env);
	// By zoush99 {@
	bool real_dim_flag = false;
	if (node->type == AP_RTYPE_REAL)
		real_dim_flag = true;
	ap_interval_t* ap_relative_err = ap_interval_alloc();
    ap_interval_t* ap_absolute_err = ap_interval_alloc();
    ap_interval_set_double(ap_relative_err, pow(2, -23), pow(2, -23));
    ap_interval_set_double(ap_absolute_err, pow(2, -149), pow(2, -149));

    itv_t delta_r; itv_init(delta_r);
    itv_t delta_a; itv_init(delta_a);
	itv_set_ap_interval(pr->itv, delta_r, ap_relative_err);
	itv_set_ap_interval(pr->itv, delta_a, ap_absolute_err);
	ap_interval_free(ap_relative_err);
	ap_interval_free(ap_absolute_err);
	
	// @}

#ifdef _T1P_DEBUG_FUN
    fprintf(stdout, "### BINARY OPERANDES ###\n");
    t1p_aff_fprint(pr, stdout, exprAB[0]);
    t1p_aff_fprint(pr, stdout, exprAB[1]);
#endif

    switch (node->op) 
      {
	case AP_TEXPR_ADD:
	      {
		res = t1p_aff_add(pr, exprAB[0], exprAB[1], env);
		if(real_dim_flag){
			itv_t itv_temp;itv_init(itv_temp);
			itv_t itv_two_delta; itv_init(itv_two_delta);
			itv_t itv_one_plus_two_delta; itv_init(itv_one_plus_two_delta);
			t1p_aff_t* temp_res;
			itv_set_int(itv_temp, 2);
			itv_mul(pr->itv, itv_two_delta, delta_r, itv_temp);
			itv_set_int(itv_temp, 1);
			itv_add(itv_one_plus_two_delta, itv_temp, delta_r);	// 1+2*delta_r
			// relative error
			temp_res=t1p_aff_mul_itv(pr, res, itv_one_plus_two_delta);	// temp_res = res*(1+2*delta_r)
			
			// absolute error
			itv_t r0_x_plus_y;
			itv_t temp1, temp2, result_interval;
			
			num_t lower_bound,upper_bound;

			itv_init(r0_x_plus_y);
			itv_init(temp1);
			itv_init(temp2);
			itv_init(result_interval);

			itv_set(r0_x_plus_y, res->c); // r0_x_plus_y = res->c

			// compute lower_bound = -4 * delta_r * (r0_x_plus_y) - 3 * delta_a
			itv_set_int(itv_temp, -4);
			itv_mul(pr->itv, temp1, r0_x_plus_y, itv_temp);   // temp1 = -4 * r0_x_plus_y
			itv_mul(pr->itv, temp2, delta_r, temp1);     // temp2 = -4 * delta_r * r0_x_plus_y
			itv_set_int(itv_temp, -3);
			itv_mul(pr->itv, temp1, delta_a, itv_temp); // temp1 = -3 * delta_a
			itv_add(temp1, temp2, temp1);          // temp1 = -4 * delta_r * r0_x_plus_y - 3 * delta_a

			num_set(lower_bound, temp1->inf); // compute lower_bound = -4 * delta_r * (r0_x_plus_y) - 3 * delta_a
			// compute upper_bound = 3 * delta_a
			itv_set_int(itv_temp, 3);
			itv_mul(pr->itv, temp2, delta_a, itv_temp); // upper_bound = 3 * delta_a
			num_set(upper_bound, temp2->sup); // compute upper_bound = 3 * delta_a

			// result_interval = [lower_bound, upper_bound]
			itv_set_num2(result_interval, lower_bound, upper_bound);

			// construct the absolute error noise symbol (only need one)
			t1p_aff_t* new_noise = create_affine_form_with_interval(pr, result_interval);

			t1p_aff_add_aff(pr, res, temp_res, new_noise);

			t1p_aff_free(pr,temp_res);
			t1p_aff_free(pr,new_noise);
			itv_clear(itv_two_delta);
			itv_clear(itv_one_plus_two_delta);
			num_clear(lower_bound);
			num_clear(upper_bound);
			itv_clear(itv_temp);
			itv_clear(temp1);
			itv_clear(temp2);
			itv_clear(result_interval);
			itv_clear(r0_x_plus_y);
		}

		break;
	      }
	case AP_TEXPR_SUB:
	      {
		res = t1p_aff_sub(pr, exprAB[0], exprAB[1], env);
		if(real_dim_flag){
			itv_t itv_temp;itv_init(itv_temp);
			itv_t itv_two_delta; itv_init(itv_two_delta);
			itv_t itv_one_plus_two_delta; itv_init(itv_one_plus_two_delta);
			t1p_aff_t* temp_res;
			itv_set_int(itv_temp, 2);
			itv_mul(pr->itv, itv_two_delta, delta_r, itv_temp);
			itv_set_int(itv_temp, 1);
			itv_add(itv_one_plus_two_delta, itv_temp, delta_r);	// 1+2*delta_r
			// relative error
			temp_res=t1p_aff_mul_itv(pr, res, itv_one_plus_two_delta);	// temp_res = res*(1+2*delta_r)
			
			// absolute error
			itv_t r0_x_plus_y;
			itv_t temp1, temp2, result_interval;
			
			num_t lower_bound,upper_bound;

			itv_init(r0_x_plus_y);
			itv_init(temp1);
			itv_init(temp2);
			itv_init(result_interval);

			itv_set(r0_x_plus_y, res->c); // r0_x_plus_y = res->c

			// compute lower_bound = -4 * delta_r * (r0_x_plus_y) - 3 * delta_a
			itv_set_int(itv_temp, -4);
			itv_mul(pr->itv, temp1, r0_x_plus_y, itv_temp);   // temp1 = -4 * r0_x_plus_y
			itv_mul(pr->itv, temp2, delta_r, temp1);     // temp2 = -4 * delta_r * r0_x_plus_y
			itv_set_int(itv_temp, -3);
			itv_mul(pr->itv, temp1, delta_a, itv_temp); // temp1 = -3 * delta_a
			itv_add(temp1, temp2, temp1);          // temp1 = -4 * delta_r * r0_x_plus_y - 3 * delta_a

			num_set(lower_bound, temp1->inf); // compute lower_bound = -4 * delta_r * (r0_x_plus_y) - 3 * delta_a
			// compute upper_bound = 3 * delta_a
			itv_set_int(itv_temp, 3);
			itv_mul(pr->itv, temp2, delta_a, itv_temp); // upper_bound = 3 * delta_a
			num_set(upper_bound, temp2->sup); // compute upper_bound = 3 * delta_a

			// result_interval = [lower_bound, upper_bound]
			itv_set_num2(result_interval, lower_bound, upper_bound);

			// construct the absolute error noise symbol (only need one)
			t1p_aff_t* new_noise = create_affine_form_with_interval(pr, result_interval);

			t1p_aff_add_aff(pr, res, temp_res, new_noise);

			t1p_aff_free(pr,temp_res);
			t1p_aff_free(pr,new_noise);
			itv_clear(itv_two_delta);
			itv_clear(itv_one_plus_two_delta);
			num_clear(lower_bound);
			num_clear(upper_bound);
			itv_clear(itv_temp);
			itv_clear(temp1);
			itv_clear(temp2);
			itv_clear(result_interval);
			itv_clear(r0_x_plus_y);
		}

		break;
	      }
	case AP_TEXPR_MUL:
	      {
		res = t1p_aff_mul(pr, exprAB[0], exprAB[1], env);
		if (real_dim_flag) {
			// 初始化变量
			itv_t itv_temp; itv_init(itv_temp);
			itv_t itv_three; itv_init(itv_three);
			itv_t itv_one; itv_init(itv_one);
			itv_t itv_delta_r_term; itv_init(itv_delta_r_term);
			itv_t itv_delta_a_term; itv_init(itv_delta_a_term);
			itv_t itv_relative_error; itv_init(itv_relative_error);
			itv_t itv_absolute_error; itv_init(itv_absolute_error);
			itv_t result_interval; itv_init(result_interval);
		
			num_t lower_bound, upper_bound;
			num_init(lower_bound);
			num_init(upper_bound);
			t1p_aff_t * r_x_times_y;
			t1p_aff_t * r_x_plus_y;
		
			// compute r_x_times_y = r_x * r_y
			r_x_times_y=res;
		
			// compute r_x_plus_y = r_x + r_y
			r_x_plus_y=t1p_aff_add(pr, exprAB[0], exprAB[1],env); 
		
			// compute delta_r_term = delta_r * (3 + delta_a)
			itv_set_int(itv_temp, 3);
			itv_add(itv_temp, itv_temp, delta_a); // itv_temp = 3 + delta_a
			itv_mul(pr->itv, itv_delta_r_term, delta_r, itv_temp); // itv_delta_r_term = delta_r * (3 + delta_a)
		
			// compute relative_error = (1 + delta_r_term) * r_x_times_y
			itv_set_int(itv_temp, 1);
			// itv_relative_error = 1 + delta_r * (3 + delta_a)
			itv_add(itv_relative_error, itv_temp, itv_delta_r_term); 
			t1p_aff_mul_itv_inplace(pr, r_x_times_y, itv_relative_error);
		
			// compute delta_a_term = delta_a * (1 + delta_a)
			itv_set_int(itv_temp, 1);
			itv_add(itv_temp, itv_temp, delta_a); // itv_temp = 1 + delta_a
			itv_mul(pr->itv, itv_delta_a_term, delta_a, itv_temp); // itv_delta_a_term = delta_a * (1 + delta_a)
		
			// compute absolute_error = delta_a_term * r_x_plus_y
			t1p_aff_mul_itv_inplace(pr, r_x_plus_y, itv_absolute_error);

			// compute result_interval = relative_error + absolute_error
			// itv_add(result_interval, itv_relative_error, itv_absolute_error);
			t1p_aff_add_aff(pr, res, r_x_times_y, r_x_plus_y);
		
			// \todo By zoush99
			// construct new noise symbol for the result
			itv_t itv_temp1; itv_init(itv_temp1);
			// 计算 delta_r_term = 3 * delta_r
			itv_set_int(itv_temp, 3);							  // itv_temp = 3
			itv_mul(pr->itv, itv_delta_r_term, delta_r, itv_temp); // itv_delta_r_term = 3 * delta_r
			
			// 计算 delta_a_term = delta_a * (1 + delta_r + delta_a)
			itv_add(itv_temp1, delta_r, delta_a);                  // itv_temp1 = delta_r + delta_a
			itv_set_int(itv_temp,1);
			itv_add(itv_temp, itv_temp, itv_temp1);                 // itv_temp = 1 + delta_r + delta_a
			itv_mul(pr->itv, itv_delta_a_term, delta_a, itv_temp); // itv_delta_a_term = delta_a * (1 + delta_r + delta_a)
			
			// 计算 lower_bound = -1 - itv_delta_r_term - itv_delta_a_term
			itv_set_int(itv_temp, -1);                         // itv_temp = -1
			itv_sub(itv_temp, itv_temp, itv_delta_r_term);         // itv_temp = -1 - 3 * delta_r
			itv_sub(itv_temp, itv_temp, itv_delta_a_term);         // itv_temp = -1 - 3 * delta_r - delta_a * (1 + delta_r + delta_a)
			
			// 设置最终区间 result_interval = [lower_bound, 0]
			num_set(lower_bound, itv_temp->inf);
			num_set_int(upper_bound, 0);
			itv_set_num2(result_interval, lower_bound, upper_bound);	// result_interval = [-1 - 3 * delta_r - delta_a * (1 + delta_r + delta_a), 0]
			
			// 清理资源
			
			t1p_aff_t* new_noise = create_affine_form_with_interval(pr, result_interval);
		
			// 将绝对误差噪声符号添加到结果中
			t1p_aff_add_aff(pr, res, res, new_noise);
		
			// 清理资源
			itv_clear(itv_temp);
			itv_clear(itv_temp1);
			itv_clear(itv_delta_r_term);
			itv_clear(itv_delta_a_term);
			itv_clear(itv_relative_error);
			itv_clear(itv_absolute_error);
			itv_clear(result_interval);
			num_clear(lower_bound);
			num_clear(upper_bound);
			t1p_aff_clear(pr, r_x_times_y);
			t1p_aff_clear(pr, r_x_plus_y);
			t1p_aff_clear(pr,new_noise);
		}
		break;
	      }
	case AP_TEXPR_DIV:
		{
		// Compute res = exprAB[0] / exprAB[1]
		res = t1p_aff_div(pr, exprAB[0], exprAB[1], env);
	
		if (real_dim_flag) {
			// Initialize variables
			itv_t itv_temp, itv_temp1,itv_delta_r_term, itv_delta_a_term;
			itv_t result_interval1, result_interval2, result_interval3;
			t1p_aff_t *r_x_div_y, *new_noise1, *new_noise2, *new_noise3;
	
			itv_init(itv_temp);
			itv_init(itv_temp1);
			itv_init(itv_delta_r_term);
			itv_init(itv_delta_a_term);
			itv_init(result_interval1);
			itv_init(result_interval2);
			itv_init(result_interval3);
	
			// Compute r_x_div_y = r_x / r_y
			r_x_div_y = t1p_aff_copy(pr, res);
	
			// --- Compute Δ1 = [-2δ_r * r_0^x - δ_a * (r_0^z + 1), 0] ---
			itv_t r0_x, r0_z;
			itv_init(r0_x);
			itv_init(r0_z);
	
			itv_set(r0_x, exprAB[0]->c); // r0_x = midpoint of r_x
			itv_set(r0_z, r_x_div_y->c); // r0_z = midpoint of r_x / r_y
			
			itv_set_int(itv_temp, -2);
			itv_mul(pr->itv,itv_temp, r0_x, itv_temp); // itv_temp = -2 * δ_r * r0_x
			itv_add(itv_temp1, itv_temp, r0_z);           // itv_temp = -2 * δ_r * r0_x - δ_a * r0_z
			itv_set_int(itv_temp, 1);
			itv_add(itv_temp, itv_temp, itv_temp1);        // itv_temp = -2 * δ_r * r0_x - δ_a * (r0_z + 1)
	
			num_t lower_bound, upper_bound;
			num_init(lower_bound);
			num_init(upper_bound);
	
			num_set(lower_bound, itv_temp->inf); // lower_bound = -2 * δ_r * r0_x - δ_a * (r0_z + 1)
			num_set_int(upper_bound, 0);         // upper_bound = 0
			itv_set_num2(result_interval1, lower_bound, upper_bound);
	
			new_noise1 = create_affine_form_with_interval(pr, result_interval1);
	
			// --- Compute Δ2 = [-(1 + δ_r)δ_r * r_0^y, 0] ---
			itv_t r0_y;
			itv_init(r0_y);
			itv_set(r0_y, exprAB[1]->c); // r0_y = midpoint of r_y
	
			itv_set_int(itv_temp, 1);
			itv_add(itv_temp, itv_temp, delta_r);         // itv_temp = 1 + δ_r
			itv_mul(pr->itv, itv_temp, itv_temp, delta_r); // itv_temp = (1 + δ_r) * δ_r
			itv_mul(pr->itv, itv_temp, itv_temp, r0_y);    // itv_temp = (1 + δ_r) * δ_r * r0_y
	
			num_set(lower_bound, itv_temp->inf); // lower_bound = -(1 + δ_r) * δ_r * r0_y
			num_set_int(upper_bound, 0);         // upper_bound = 0
			itv_set_num2(result_interval2, lower_bound, upper_bound);
	
			new_noise2 = create_affine_form_with_interval(pr, result_interval2);
	
			// --- Compute Δ3 = [-(1 + δ_r) * r_0^z - δ_a, δ_a] ---
			itv_set_int(itv_temp, 1);
			itv_add(itv_temp, itv_temp, delta_r);         // itv_temp = 1 + δ_r
			itv_mul(pr->itv, itv_temp, itv_temp, r0_z);    // itv_temp = (1 + δ_r) * r0_z
			itv_add(itv_temp, itv_temp, delta_a);         // itv_temp = (1 + δ_r) * r0_z + δ_a
	
			num_set(lower_bound, itv_temp->inf); // lower_bound = -(1 + δ_r) * r0_z - δ_a
			num_set(upper_bound, delta_a->sup);  // upper_bound = δ_a
			itv_set_num2(result_interval3, lower_bound, upper_bound);
	
			new_noise3 = create_affine_form_with_interval(pr, result_interval3);

			// 初始化变量
			itv_t term1, term2;
			itv_init(term1);
			itv_init(term2);

			// 计算 term1 = 2 * δ_r * r_x
			itv_set_int(itv_temp, 2); // itv_temp = 2
			itv_mul(pr->itv, term1, delta_r, itv_temp); // term1 = 2 * δ_r
			t1p_aff_t* aff_temp1 = t1p_aff_mul_itv(pr, exprAB[0], term1);

			// 计算 term2 = δ_a * (r_x / r_y + 1)
			// 计算 r_x / r_y + 1
			itv_set_int(itv_temp, 1); // itv_temp = 1
			t1p_aff_t* aff_temp2=t1p_aff_copy(pr, r_x_div_y);
			itv_add(aff_temp2->c,r_x_div_y->c,itv_temp); // itv_temp = r_x / r_y + 1
			// 计算 δ_a * (r_x / r_y + 1)
			aff_temp2=t1p_aff_mul_itv(pr,aff_temp2,delta_a);

			// 计算最终结果 result = term1 + term2
			t1p_aff_add_aff(pr,aff_temp1,aff_temp1,aff_temp2);

			t1p_aff_add_aff(pr, aff_temp1, aff_temp1, new_noise1);

			//{@ 
			// 计算 (1 + δ_r)
			itv_set_int(itv_temp, 1); // itv_temp = 1
			itv_add(itv_temp, itv_temp, delta_r); // itv_temp = 1 + δ_r

			// 计算 (1 + δ_r) * δ_r
			itv_mul(pr->itv, itv_temp, itv_temp, delta_r); // itv_temp = (1 + δ_r) * δ_r

			// 计算 (1 + δ_r) * δ_r * r_y
			aff_temp2=t1p_aff_mul_itv(pr,exprAB[1],itv_temp);

			t1p_aff_add_aff(pr, aff_temp1, aff_temp2, new_noise2);

			// @}

			aff_temp1=t1p_aff_div(pr, aff_temp1, aff_temp2, env);

			// {@
			itv_set_int(itv_temp, 1); // itv_temp = 1
			itv_add(itv_temp, itv_temp, delta_r); // itv_temp = 1 + δ_r
			t1p_aff_t* aff_temp3 = t1p_aff_mul_itv(pr, r_x_div_y, itv_temp);
			// @}

			t1p_aff_add_aff(pr, aff_temp1, aff_temp1, aff_temp3);
			t1p_aff_add_aff(pr, res, aff_temp1, new_noise3);

			// 清理资源
			itv_clear(itv_temp);
			itv_clear(itv_temp1);
			itv_clear(itv_delta_r_term);
			itv_clear(itv_delta_a_term);
			itv_clear(result_interval1);
			itv_clear(result_interval2);
			itv_clear(result_interval3);
			itv_clear(r0_x);
			itv_clear(r0_y);
			itv_clear(r0_z);
			t1p_aff_clear(pr, r_x_div_y);
			t1p_aff_clear(pr, new_noise1);
			t1p_aff_clear(pr, new_noise2);
			t1p_aff_clear(pr, new_noise3);
			num_clear(lower_bound);
			num_clear(upper_bound);
		}
		break;
	}
	case AP_TEXPR_MOD:
	      {
		res = t1p_aff_mod(pr, exprAB[0], exprAB[1]);
		break;
	      }
	default: fatal("Unknown binary operation");
      }
    //if (node->exprA->discr != AP_TEXPR_DIM) t1p_aff_free(pr, exprAB[0]);
    //if (node->exprB->discr != AP_TEXPR_DIM) t1p_aff_free(pr, exprAB[1]);

    t1p_aff_check_free(pr, exprAB[0]);
    t1p_aff_check_free(pr, exprAB[1]);
#ifdef _T1P_DEBUG_FUN
    fprintf(stdout, "### BINARY RESULT ###\n");
    t1p_aff_fprint(pr, stdout, res);
#endif
    // By zoush99 { @}
    itv_clear(delta_r);
    itv_clear(delta_a);
    // @}
    return res;
}

t1p_aff_t* t1p_aff_eval_node( t1p_internal_t* pr, ap_texpr0_node_t* node, t1p_t* env)
{
    arg_assert((node && env), abort(););
    if (node->exprB) return t1p_aff_eval_node_binary(pr, node, env);
    else return t1p_aff_eval_node_unary(pr, node, env);
}

t1p_aff_t* t1p_aff_eval_ap_texpr0(t1p_internal_t* pr, ap_texpr0_t* expr, t1p_t* env)
{
    t1p_aff_t* res;
    arg_assert((env && expr), abort(););
#ifdef _T1P_DEBUG_FUN
    fprintf(stdout, "### eval TEXPR ###\n");
    ap_texpr0_fprint(stdout, expr, NULL);
    fprintf(stdout, "\n");
#endif
    switch(expr->discr){
	case AP_TEXPR_CST:
	      {
		res = t1p_aff_alloc_init(pr);
		itv_set_ap_coeff(pr->itv, res->c, &(expr->val.cst));
		itv_set(res->itv, res->c);
		break;
	      }
	case AP_TEXPR_DIM:
	      {
		/* involved affine forms should have ->pby > 0 */
		res = env->paf[expr->val.dim];
		if (res->pby == 0) fatal("partage des formes affines foireux...");
		res->pby++;
		break;
	      }
	case AP_TEXPR_NODE:
	      {
		res = t1p_aff_eval_node(pr, expr->val.node, env);
		break;
	      }
	default:
	    fatal("Unknown texpr type");
    }
#ifdef _T1P_DEBUG_FUN
    fprintf(stdout, "### RESULT of eval TEXPR ###\n");
    t1p_aff_fprint(pr, stdout, res);
    fprintf(stdout, "\n");
#endif
    return res;
}

/* [a,b]*[-1,1] = [min(a,-b), max(-a,b)] 
 * => sum{[ai,bi]*[-1,1]} = [sum(min(ai,-bi)), sum(max(-ai,bi))] (itv1) 
 * [a,b]*[0,1] = [min(a,0),max(0,b)]
 * => sum{[ai,bi]*[0,1]} = [sum(min(ai,0)), sum(max(0,bi))] (itv2)
 * => non linear part \in (itv1 + itv2)
 * add_itv(itv1 + itv2)
*/
void square_dep(t1p_internal_t* pr, itv_t res, t1p_aaterm_t* p, t1p_aaterm_t* q, int* hash, unsigned int dim, itv_t* gamma)
{
    t1p_aaterm_t* ptr_p; 
    t1p_aaterm_t* ptr_q;
    itv_t tmp, tmp1, tmp2;
    itv_init(tmp); itv_init(tmp1); itv_init(tmp2);
    itv_t ** itv_matrix;
    checked_malloc(itv_matrix, itv_t*, dim, abort(););
    size_t i = 0; size_t j = 0;

    if (gamma == NULL) {
	itv_t zeroone, moneone;
	itv_init(zeroone);
	itv_init(moneone);
	itv_set_int2(zeroone,(long int)0,(long int)1);
	itv_set_int2(moneone,(long int)(-1),(long int)1);
	for (i=0; i<dim; i++) itv_matrix[i] = itv_array_alloc(dim);
	for (ptr_p=p; ptr_p; ptr_p = ptr_p->n) {
	    for (ptr_q=q; ptr_q; ptr_q = ptr_q->n) {
		itv_mul(pr->itv, tmp, ptr_p->coeff, ptr_q->coeff);
		itv_set(itv_matrix[hash[ptr_p->pnsym->index]-1][hash[ptr_q->pnsym->index]-1], tmp);
	    }
	}
	for (i=1; i<=dim; i++) {
	    for (j=1; j<=i; j++) {
		if (i==j) {
		    itv_mul(pr->itv, tmp, itv_matrix[i-1][j-1], zeroone);
		    itv_add(res,res,tmp);
		} else {
		    itv_neg(tmp1, itv_matrix[i-1][j-1]);
		    if (!(itv_is_eq(tmp1,itv_matrix[j-1][i-1]))) {
			itv_add(tmp1, itv_matrix[i-1][j-1], itv_matrix[j-1][i-1]);
			itv_mul(pr->itv, tmp, tmp1, moneone);
			itv_add(res,res,tmp);
		    }
		}
	    }
	}
	itv_clear(zeroone); itv_clear(moneone);
	for (i=0; i<dim; i++) {
	    itv_array_free(itv_matrix[i],dim);
	}
	free(itv_matrix);
    } else {
	ap_dim_t dim1 = 0;
	ap_dim_t dim2 = 0;
	itv_t eps_itv1, eps_itv2;
	itv_init(eps_itv1);
	itv_init(eps_itv2);
	itv_t ** itv_matrix2;
	checked_malloc(itv_matrix2, itv_t*, dim, abort(););
	for (i=0; i<dim; i++) {
	    itv_matrix[i] = itv_array_alloc(dim);
	    itv_matrix2[i] = itv_array_alloc(i+1);
	}
	for (ptr_p=p; ptr_p; ptr_p = ptr_p->n) {
	    for (ptr_q=q; ptr_q; ptr_q = ptr_q->n) {
		itv_mul(pr->itv, tmp, ptr_p->coeff, ptr_q->coeff);
		itv_set(itv_matrix[hash[ptr_p->pnsym->index]-1][hash[ptr_q->pnsym->index]-1], tmp);
		itv_mul(pr->itv, tmp, gamma[ptr_p->pnsym->index], gamma[ptr_q->pnsym->index]);
		if (hash[ptr_p->pnsym->index]-1 < hash[ptr_q->pnsym->index]-1) {
		    itv_set(itv_matrix2[hash[ptr_q->pnsym->index] -1][hash[ptr_p->pnsym->index] -1], tmp);
		} else {
		    itv_set(itv_matrix2[hash[ptr_p->pnsym->index]-1][hash[ptr_q->pnsym->index]-1], tmp);
		}
	    }
	}
	for (i=1; i<=dim; i++) {
	    for (j=1; j<=i; j++) {
		if (i==j) {
		    if (bound_sgn(itv_matrix2[i-1][j-1]->inf) > 0) bound_set_int(itv_matrix2[i-1][j-1]->inf,0);
		    itv_mul(pr->itv, tmp, itv_matrix[i-1][j-1], itv_matrix2[i-1][j-1]);
		    itv_add(res,res,tmp);
		} else {
		    itv_neg(tmp1, itv_matrix[i-1][j-1]);
		    if (!(itv_is_eq(tmp1,itv_matrix[j-1][i-1]))) {
			itv_add(tmp1, itv_matrix[i-1][j-1], itv_matrix[j-1][i-1]);
			itv_mul(pr->itv, tmp, tmp1, itv_matrix2[i-1][j-1]);
			itv_add(res,res,tmp);
		    }
		}
	    }
	}
	for (i=0; i<dim; i++) {
	    itv_array_free(itv_matrix[i],dim);
	    itv_array_free(itv_matrix2[i],i+1);
	}
	free(itv_matrix);
	free(itv_matrix2);
	itv_clear(eps_itv1);
	itv_clear(eps_itv2);
    }
    itv_clear(tmp); itv_clear(tmp1); itv_clear(tmp2);
}

    void buildIntervalUpperTriangle(t1p_internal_t* pr, itv_t* array, t1p_aaterm_t* p, t1p_aaterm_t* q, size_t dim, int* hash) 
{
}

double* buildDoubleUpperTriangle(t1p_internal_t* pr, itv_t* array, size_t dim) {
    itv_t mid, dev;
    double dbl;
    size_t i;
    double* res = (double*)calloc((dim*(dim+1)/2),sizeof(double));
    for (i=0; i<(dim*(dim+1)/2); i++) {
    	itv_middev(pr->itv, mid, dev, array[i]);
    	double_set_num(&dbl, mid->sup);
	res[i] = dbl;
    }
    return res;
}

#ifndef _USE_SDP
bool call_sdp(itv_internal_t* itv, itv_t res, t1p_aaterm_t* p, t1p_aaterm_t* q, size_t dim, int* hash, bool square) 
{
	printf("SDP not supported, you have to enable it and rebuild the library\n");
	abort();
}
#else
bool call_sdp(itv_internal_t* itv, itv_t res, t1p_aaterm_t* p, t1p_aaterm_t* q, size_t dim, int* hash, bool square)
{
    //	printf("Calling SDP %d\n",dim);
    t1p_aaterm_t* ptr_p; 
    t1p_aaterm_t* ptr_q;
    itv_t tmp;
    itv_t mid, dev;
    itv_init(tmp);
    itv_init(mid);
    itv_init(dev);

    //printf("dim of SDP problem: %d\n",dim);fflush(stdout);

    /************* SDP data ************/
    unsigned int i=1,j=1;
    /*
     * The problem and solution data.
     */

    struct blockmatrix C;
    double *b;
    struct constraintmatrix *constraints;

    /*
     * Storage for the initial and final solutions.
     */

    struct blockmatrix X,Z;
    double *y;
    double pobj,dobj;

    num_t sol1, sol2;

    /*
     * blockptr will be used to point to blocks in constraint matrices.
     */

    struct sparseblock *blockptr;

    /*
     * A return code for the call to easy_sdp().
     */

    int ret;

    double dbl = 0;
    num_t temp;
    size_t dim_of_X, nb_of_cons;
    /************* end SDP data ****************/

    /*
     * The first major task is to setup the C matrix and right hand side b.
     */

    /*
     * First, allocate storage for the C matrix.  We have always two blocks, but
     * because C starts arrays with index 0, we have to allocate space for
     * three blocks- we'll waste the 0th block.  Notice that we check to 
     * make sure that the malloc succeeded.
     */

    C.nblocks=2;
    C.blocks=(struct blockrec *)malloc((2+1)*sizeof(struct blockrec));
    if (C.blocks == NULL)
      {
	printf("Couldn't allocate storage for C!\n");
	exit(1);
      };

    /*
     * Setup the first block.
     */

    C.blocks[1].blockcategory=MATRIX;
    C.blocks[1].blocksize=dim;
    C.blocks[1].data.mat=(double *)malloc(dim*dim*sizeof(double));
    if (C.blocks[1].data.mat == NULL)
      {
	printf("Couldn't allocate storage for C!\n");
	exit(1);
      };

    /*
     * Put the entries into the first block.
     */

    /* Initialize */
    for (i=1; i<=dim; i++) {
	for (j=1; j<=dim; j++) {
	    C.blocks[1].data.mat[ijtok(i,j,dim)] = 0.0;
	}
    }

    num_t eps; num_init(eps);
    bound_t err; bound_init(err);
    num_set_double(eps,(double)30*1.11022302462515654042e-16);
    for (ptr_p=p; ptr_p; ptr_p = ptr_p->n) {
	for (ptr_q=q; ptr_q; ptr_q = ptr_q->n) {
	    /* TODO: We shall resolve an interval SDP problem and not reduce the interval to tmp->sup */
	    itv_mul(itv, tmp, ptr_p->coeff, ptr_q->coeff);
	    if (itv_is_point(itv, tmp)) double_set_num(&dbl, tmp->sup); 
	    else {
		itv_range_rel(itv,err,tmp);
		itv_middev(itv, mid, dev, tmp);
		bound_print(err);
		printf("\n");
		//	printf("\ndev:");
		//	itv_print(dev);
		//	printf("\n");
		//	num_print(eps);
		if (bound_cmp_num(err,eps) <= 0) {
		    double_set_num(&dbl, mid->sup);
		} else {
		    printf("Warning: (SDP) "); itv_print(mid); printf("is taken instead of "); itv_print(tmp);printf(" continue ...\n");abort();
		}
	    }
	    if (ptr_p->pnsym->index == ptr_q->pnsym->index) {
		C.blocks[1].data.mat[ijtok(hash[ptr_p->pnsym->index],hash[ptr_q->pnsym->index],dim)]=dbl;
	    } else {
		C.blocks[1].data.mat[ijtok(hash[ptr_p->pnsym->index],hash[ptr_q->pnsym->index],dim)] += dbl/2;
		C.blocks[1].data.mat[ijtok(hash[ptr_q->pnsym->index],hash[ptr_p->pnsym->index],dim)] += dbl/2;
	    }
	}
    }

    /*
     * Setup the second block.
     */

    C.blocks[2].blockcategory=MATRIX;
    C.blocks[2].blocksize=dim;
    C.blocks[2].data.mat=(double *)malloc(dim*dim*sizeof(double));
    if (C.blocks[2].data.mat == NULL)
      {
	printf("Couldn't allocate storage for C!\n");
	exit(1);
      };
    /*
     * Put the entries into the second block.
     */
    for (i=1; i<=dim; i++) {
	for (j=1; j<=dim; j++) {
	    C.blocks[2].data.mat[ijtok(i,j,dim)] = 0.0;
	}
    }

    /*
     * Allocate storage for the right hand side, b.
     */

    b=(double *)malloc((dim+1)*sizeof(double));
    if (b==NULL)
      {
	printf("Failed to allocate storage for a!\n");
	exit(1);
      };

    /*
     * Fill in the entries in b.
     */

    for (i=1; i<=dim; i++) {
	b[i]=1.0;
    }

    /*
     * The next major step is to setup the constraint matrices Ai.
     * Again, because C indexing starts with 0, we have to allocate space for
     * one more constraint.  constraints[0] is not used.
     */

    constraints=(struct constraintmatrix *)malloc((dim+1)*sizeof(struct constraintmatrix));
    if (constraints==NULL)
      {
	printf("Failed to allocate storage for constraints!\n");
	exit(1);
      };

    /*
     * Setup the A1 matrix.  Note that we start with block 3 of A1 and then
     * do block 1 of A1.  We do this in this order because the blocks will
     * be inserted into the linked list of A1 blocks in reverse order.  
     */

    /*
     * Terminate the linked list with a NULL pointer.
     */

    for (i=1; i<=dim; i++) {
	constraints[i].blocks=NULL;

	/*
	 * Now, we handle block 2 of Ai.
	 */

	/*
	 * Allocate space for block 2 of Ai.
	 */

	blockptr=(struct sparseblock *)malloc(sizeof(struct sparseblock));
	if (blockptr==NULL)
	  {
	    printf("Allocation of constraint block failed!\n");
	    exit(1);
	  };

	/*
	 * Initialize block 2.
	 */

	blockptr->blocknum=2;
	blockptr->blocksize=dim;
	blockptr->constraintnum=i;
	blockptr->next=NULL;
	blockptr->nextbyblock=NULL;
	blockptr->entries=(double *) malloc((1+1)*sizeof(double));
	if (blockptr->entries==NULL)
	  {
	    printf("Allocation of constraint block failed!\n");
	    exit(1);
	  };
	blockptr->iindices=(int *) malloc((1+1)*sizeof(int));
	if (blockptr->iindices==NULL)
	  {
	    printf("Allocation of constraint block failed!\n");
	    exit(1);
	  };
	blockptr->jindices=(int *) malloc((1+1)*sizeof(int));
	if (blockptr->jindices==NULL)
	  {
	    printf("Allocation of constraint block failed!\n");
	    exit(1);
	  };

	/*
	 * We have 1 nonzero entry in the upper triangle of block 2 of Ai.
	 */

	blockptr->numentries=1;

	/*
	 * The entry in the i,i position of block 2 of Ai is 1.0
	 */

	blockptr->iindices[1]=i;
	blockptr->jindices[1]=i;
	blockptr->entries[1]=1.0;

	/*
	 * Note that the entry in the other positions of block 2 of Ai is 0, 
	 * So we don't store anything for it.
	 */

	/*
	 * Insert block 2 into the linked list of Ai blocks.  
	 */

	blockptr->next=constraints[i].blocks;
	constraints[i].blocks=blockptr;

	/*
	 * Now, we handle block 1.  
	 */

	/*
	 * Allocate space for block 1.
	 */

	blockptr=(struct sparseblock *)malloc(sizeof(struct sparseblock));
	if (blockptr==NULL)
	  {
	    printf("Allocation of constraint block failed!\n");
	    exit(1);
	  };

	/*
	 * Initialize block 1.
	 */

	blockptr->blocknum=1;
	blockptr->blocksize=dim;
	blockptr->constraintnum=i;
	blockptr->next=NULL;
	blockptr->nextbyblock=NULL;
	blockptr->entries=(double *) malloc((1+1)*sizeof(double));
	if (blockptr->entries==NULL)
	  {
	    printf("Allocation of constraint block failed!\n");
	    exit(1);
	  };
	blockptr->iindices=(int *) malloc((1+1)*sizeof(int));
	if (blockptr->iindices==NULL)
	  {
	    printf("Allocation of constraint block failed!\n");
	    exit(1);
	  };
	blockptr->jindices=(int *) malloc((1+1)*sizeof(int));
	if (blockptr->jindices==NULL)
	  {
	    printf("Allocation of constraint block failed!\n");
	    exit(1);
	  };

	/*
	 * We have 3 nonzero entries in the upper triangle of block 1 of A1.
	 */

	blockptr->numentries=1;

	/*
	 * The entry in the i,i position of block 1 of Ai is 1.0
	 */

	blockptr->iindices[1]=i;
	blockptr->jindices[1]=i;
	blockptr->entries[1]=1.0;

	/*
	 * Insert block 1 into the linked list of A1 blocks.  
	 */

	blockptr->next=constraints[i].blocks;
	constraints[i].blocks=blockptr;
    }


    /*
     * At this point, we have all of the problem data setup.
     */

    /*
     * Write the problem out in SDPA sparse format.
     */
    dim_of_X = 2*dim;
    nb_of_cons = dim;

    write_prob("prob.dat-s", dim_of_X, nb_of_cons, C, b, constraints);

    /*
     * Create an initial solution.  This allocates space for X, y, and Z,
     * and sets initial values.
     */

    initsoln(dim_of_X, nb_of_cons, C, b, constraints, &X, &y, &Z);

    /*
     * Solve the problem.
     */

    ret=easy_sdp(dim_of_X, nb_of_cons,C,b,constraints,0.0,&X,&y,&Z,&pobj,&dobj);

    /*
     * Write out the problem solution.
     */

    //num_set_double(sol1, (dobj+pobj)/(2));
    //	printf("primal : %.16e \t dual : %.16e\n",pobj,dobj);
    num_set_double(sol1, dobj);
    //	if (square) {
    //		num_set_double(sol2, (double)0.0);
    //		itv_set_num2(res, sol1, sol2);
    //		free_prob(dim_of_X, nb_of_cons, C,b,constraints,X,y,Z);
    //		return ret;
    //	}
    //	write_sol("prob.sol", dim_of_X, nb_of_cons, X, y, Z);


    /*
     * solve the same problem with C = -C
     */

    for (i=1; i<=dim; i++) {
	for (j=1; j<=dim; j++) {
	    C.blocks[1].data.mat[ijtok(i,j,dim)] = -1*C.blocks[1].data.mat[ijtok(i,j,dim)];
	}
    }
    //write_prob("prob.dat-s", dim_of_X, nb_of_cons, C, b, constraints);
    initsoln(dim_of_X, nb_of_cons, C, b, constraints, &X, &y, &Z);
    ret+=easy_sdp(dim_of_X, nb_of_cons,C,b,constraints,0.0,&X,&y,&Z,&pobj,&dobj);

    //num_set_double(sol2, (dobj+pobj)/(-2));
    //	printf("primal : %.16e \t dual : %.16e\n",-pobj,-dobj);
    num_set_double(sol2, -dobj);
    itv_set_num2(res, sol1, sol2);

    /*
     * Free storage allocated for the problem and return.
     */

    free_prob(dim_of_X, nb_of_cons, C,b,constraints,X,y,Z);
    itv_clear(tmp);
    itv_clear(mid);
    itv_clear(dev);
    num_clear(eps);
    bound_clear(err);
    return ret;
}

#endif
