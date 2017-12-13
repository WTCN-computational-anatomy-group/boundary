/* $Id: shoot_boundary.h 4875 2012-08-30 20:04:30Z john $ */
/* (c) John Ashburner (2011) */

#define BOUND_CIRCULANT 0
#define BOUND_NEUMANN   1
#define BOUND_DIRICHLET 2
#define BOUND_SLIDING   3

extern mwSignedIndex (*bound)();
extern float (*bound_factor)();
extern void set_bound(int t);
extern int  get_bound();
extern mwSignedIndex bperiod(mwSignedIndex i, mwSize m);

#define BOUND(a,b) bound(a,b)
#define BOUND_FACTOR(a,b,c) bound_factor(a,b,c)

