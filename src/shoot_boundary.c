/* $Id: shoot_boundary.c 4875 2012-08-30 20:04:30Z john $ */
/* (c) John Ashburner (2011) */

#include "mex.h"
#include "shoot_boundary.h"

/* Neumann/Dirichlet boundary condition -- index */
static mwSignedIndex neumann_dirichlet_boundary(mwSignedIndex i, mwSize m)
{
    mwSignedIndex m2 = m*2;
    i = (i<0) ? m2-((-i-1)%m2)-1 : (i%m2);
    return((m<=i)? m2-i-1: i);
}

/* Circulant boundary condition -- index */
static mwSignedIndex circulant_boundary(mwSignedIndex i, mwSize m)
{
    return((i>=0) ? i%((signed)m) : (((signed)m)+i%((signed)m))%(signed)m);
}

/* Neuman/Circulant boundary condition -- factor */
static int neumann_circulant_factor(mwSignedIndex i, mwSize m, int diag)
{
    return(1);
}

/* Dirichlet boundary condition -- factor */
static int dirichlet_factor(mwSignedIndex i, mwSize m, int diag)
{
    return(bperiod(i,m) % 2 ? -1 : 1);
}

/* Sliding (Neumann & Dirichlet) boundary condition -- factor */
static int sliding_factor(mwSignedIndex i, mwSize m, int diag)
{
    if(diag)
        return(bperiod(i,m) % 2 ? -1 : 1);
    else
        return(1);
}

mwSignedIndex (*bound)()        = circulant_boundary;
int           (*bound_factor)() = neumann_circulant_factor;
static int bound_type = BOUND_CIRCULANT;

void set_bound(int t)
{
    bound_type = t;
    if (t==BOUND_CIRCULANT)
    {
        bound = circulant_boundary;
        bound_factor = neumann_circulant_factor;
    }
    else if (t==BOUND_NEUMANN)
    {
        bound = neumann_dirichlet_boundary;
        bound_factor = neumann_circulant_factor;
    }
    else if (t==BOUND_DIRICHLET)
    {
        bound = neumann_dirichlet_boundary;
        bound_factor = dirichlet_factor;
    }
    else if (t==BOUND_SLIDING)
    {
        bound = neumann_dirichlet_boundary;
        bound_factor = neumann_circulant_factor;
    }
    else
        mexErrMsgTxt("Undefined boundary condition.");
}

int get_bound()
{
    return(bound_type);
}

mwSignedIndex bperiod(mwSignedIndex i, mwSize m)
{
    return(i >= 0 ? i/m : (i+1)/m - 1);
}


