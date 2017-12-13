/* $Id: shoot_regularisers.c 6799 2016-05-20 16:50:25Z john $ */
/* (c) John Ashburner (2011) */

#include<mex.h>
#include<math.h>
extern double log(double x);

#include "shoot_boundary.h"
/*
% MATLAB code (requiring Symbolic Toolbox) for computing the
% regulariser for linear elasticity.

syms mu lam
syms v0 v1 v2
N  = 3;
O  = sym(zeros(N));
OO = kron(kron(O,O),O);
I  = sym(eye(N));
G1 = sym(spdiags(repmat([-1 1],N,1),[ 0 1],N,N)); G1(N,1) =  1;
G2 = sym(spdiags(repmat([-1 1],N,1),[-1 0],N,N)); G2(1,N) = -1;
G  = {G1 G2};
LL = sym(zeros(3*N^3,3*N^3));
for i=1:2,
    Di = kron(I,kron(I,G{i}))/v0;
    for j=1:2,
        Dj = kron(I,kron(G{j},I))/v1;
        for k=1:2,
            Dk = kron(G{k},kron(I,I))/v2;
            L1 = [    Di     OO     OO
                      OO     Dj     OO
                      OO     OO     Dk
                  0.5*Dj 0.5*Di     OO
                  0.5*Dk     OO 0.5*Di
                      OO 0.5*Dk 0.5*Dj
                  0.5*Dj 0.5*Di     OO
                  0.5*Dk     OO 0.5*Di
                      OO 0.5*Dk 0.5*Dj]; % Elements of symmetric part of Jacobian

            L2 = [    Di     Dj  Dk];    % Divergence
            LL = LL + mu*(L1.'*L1) + lam*(L2.'*L2);
        end
    end
end
LL = LL/8;
o  = sym(ones(1,27));
LL = LL.*kron([v0*o v1*o v2*o],[o o o].').*kron([v0*o v1*o v2*o].',[o o o]);


L1 = simplify(reshape(LL(:,2+3+9   ),[3 3 3 3]))
L2 = simplify(reshape(LL(:,2+3+9+27),[3 3 3 3]))
L3 = simplify(reshape(LL(:,2+3+9+54),[3 3 3 3]))

*/
/*
% MATLAB code (requiring Symbolic Toolbox) for computing the
% regulariser for bending energy. Note that an additional multiplication
% by the square of the voxel size is needed.

syms s1 s2 s3
syms l1 l2 l3
zz = sym(zeros(3,3));
K1 = cat(3,zz,[0 -1/(v0*v0) 0; 0 2/(v0*v0) 0; 0 -1/(v0*v0) 0],zz);
K2 = cat(3,zz,[0 0 0; -1/(v1*v1) 2/(v1*v1) -1/(v1*v1); 0 0 0],zz);
K3 = sym(zeros(3,3,3));
K3(2,2,1) = -1/(v2*v2);
K3(2,2,2) = 2/(v2*v2);
K3(2,2,3) = -1/(v2*v2);

K  = K1+K2+K3;
K1 = K*l1; K1(2,2,2) = K1(2,2,2)+l2;
K2 = K*l1; K2(2,2,2) = K2(2,2,2)+l3;

% L  = convn(K,K)
L  = sym(zeros(5,5,5));
for i=1:3,
    for j=1:3,
        for k=1:3,
            L(i-1+1:i+1+1,j-1+1:j+1+1,k-1+1:k+1+1) = L(i-1+1:i+1+1,j-1+1:j+1+1,k-1+1:k+1+1) + K1(i,j,k)*K2;
        end;
    end;
end;
disp(simplify(L(3:end,3:end,3:end)))
*/


double trapprox(mwSize dm[], float a[], double s[])
{
    double v0 = s[0]*s[0], v1 = s[1]*s[1], v2 = s[2]*s[2];
    double lam0 = s[3], lam1 = s[4], lam2 = s[5], mu = s[6], lam = s[7];
    double w000, wx000, wy000, wz000;
    double tr = 0.0;
    mwSignedIndex i, j, k;

    w000  =  lam2*(6*(v0*v0+v1*v1+v2*v2) +8*(v0*v1+v0*v2+v1*v2)) +lam1*2*(v0+v1+v2) + lam0;
    wx000 =  2*mu*(2*v0+v1+v2)/v0+2*lam + w000/v0;
    wy000 =  2*mu*(v0+2*v1+v2)/v1+2*lam + w000/v1;
    wz000 =  2*mu*(v0+v1+2*v2)/v2+2*lam + w000/v2;

    for(k=0; k<dm[2]; k++)
    {
        for(j=0; j<dm[1]; j++)
        {
           float *paxx, *payy, *pazz, *paxy, *paxz, *payz;

            paxx = a+dm[0]*(j+dm[1]* k);
            payy = a+dm[0]*(j+dm[1]*(k+dm[2]));
            pazz = a+dm[0]*(j+dm[1]*(k+dm[2]*2));
            paxy = a+dm[0]*(j+dm[1]*(k+dm[2]*3));
            paxz = a+dm[0]*(j+dm[1]*(k+dm[2]*4));
            payz = a+dm[0]*(j+dm[1]*(k+dm[2]*5));

            for(i=0; i<dm[0]; i++)
            {
                double axx, ayy, azz, axy, axz, ayz, dt;
                axx  = paxx[i] + wx000;
                ayy  = payy[i] + wy000;
                azz  = pazz[i] + wz000;
                axy  = paxy[i];
                axz  = paxz[i];
                ayz  = payz[i];
                dt   = axx*ayy*azz -axx*ayz*ayz-ayy*axz*axz-azz*axy*axy +2*axy*axz*ayz;
                tr  += (wx000*(ayy*azz - ayz*ayz) 
                       +wy000*(axx*azz - axz*axz)
                       +wz000*(axx*ayy - axy*axy))/dt;

            }
        }
    }
    return(tr);
}

/* Returns the velocity-to-momentum convolution kernel.
 *
 * \param[in]  dm     Dimensions of the velocity lattice (voxels)
 * \param[in]  s[0:2] Voxel size of the lattice (voxels/mm)
 * \param[in]  s[3]   Parameter of the absolute displacement penalty
 * \param[in]  s[4]   Parameter of the membrane energy (penalizes 
 *                    elements of the Jacobian matrix -> 1st order 
 *                    smoothness)
 * \param[in]  s[5]   Parameter of the bending energy (penalizes 
 *                    elements of the Hessian matrix -> 2nd order 
 *                    smoothness)
 * \param[in]  s[6]   Parameter of the linear elastic energy (penalizes 
 *                    elements of the symmetric part of the Jacobian matrix 
 *                    -> penalizes scaling and shearing)
 * \param[in]  s[7]   Parameter of the linear elastic energy (penalizes 
 *                    the divergence of the Jacobian matrix -> preserves 
 *                    volumes)
 * \param[out] f      Allocated array of size dm in which to store the  
 *                    kernel
 */
void kernel(mwSize dm[], double s[], float f[])
{
    double w000,w100,w200,
           w010,w110,
           w020,
           w001,w101,
           w011,
           w002;
    double v0 = s[0]*s[0], v1 = s[1]*s[1], v2 = s[2]*s[2];
    double lam0 = s[3], lam1 = s[4], lam2 = s[5], mu = s[6], lam = s[7];
    mwSignedIndex im2,im1,ip1,ip2, jm2,jm1,jp1,jp2, km2,km1,kp1,kp2;

    km2 = (bound((mwSignedIndex)-2,dm[2]))*dm[0]*dm[1];
    km1 = (bound((mwSignedIndex)-1,dm[2]))*dm[0]*dm[1];
    kp1 = (bound((mwSignedIndex) 1,dm[2]))*dm[0]*dm[1];
    kp2 = (bound((mwSignedIndex) 2,dm[2]))*dm[0]*dm[1];

    jm2 = (bound((mwSignedIndex)-2,dm[1]))*dm[0];
    jm1 = (bound((mwSignedIndex)-1,dm[1]))*dm[0];
    jp1 = (bound((mwSignedIndex) 1,dm[1]))*dm[0];
    jp2 = (bound((mwSignedIndex) 2,dm[1]))*dm[0];

    im2 = bound((mwSignedIndex)-2,dm[0]);
    im1 = bound((mwSignedIndex)-1,dm[0]);
    ip1 = bound((mwSignedIndex) 1,dm[0]);
    ip2 = bound((mwSignedIndex) 2,dm[0]);

    w000 = lam2*(6*(v0*v0+v1*v1+v2*v2) +8*(v0*v1+v0*v2+v1*v2)) +lam1*2*(v0+v1+v2) + lam0;
    w100 = lam2*(-4*v0*(v0+v1+v2)) -lam1*v0;
    w010 = lam2*(-4*v1*(v0+v1+v2)) -lam1*v1;
    w001 = lam2*(-4*v2*(v0+v1+v2)) -lam1*v2;
    w200 = lam2*v0*v0;
    w020 = lam2*v1*v1;
    w002 = lam2*v2*v2;
    w110 = lam2*2*v0*v1;
    w101 = lam2*2*v0*v2;
    w011 = lam2*2*v1*v2;

    if (mu==0 && lam==0)
    {
        f[0]           += w000;
        f[im1        ] += w100; f[ip1        ] += w100;
        f[    jm1    ] += w010; f[    jp1    ] += w010;
        f[        km1] += w001; f[        kp1] += w001;
        f[im2        ] += w200; f[ip2        ] += w200;
        f[    jm2    ] += w020; f[    jp2    ] += w020;
        f[        km2] += w002; f[        kp2] += w002;
        f[im1+jm1    ] += w110; f[ip1+jm1    ] += w110; f[im1+jp1    ] += w110; f[ip1+jp1    ] += w110;
        f[im1    +km1] += w101; f[ip1    +km1] += w101; f[im1    +kp1] += w101; f[ip1    +kp1] += w101;
        f[    jm1+km1] += w011; f[    jp1+km1] += w011; f[    jm1+kp1] += w011; f[    jp1+kp1] += w011;
    }
    else
    {
        double wx000, wx100, wx010, wx001, wy000, wy100, wy010, wy001, wz000, wz100, wz010, wz001, w2;
        float *pxx, *pxy, *pxz, *pyx, *pyy, *pyz, *pzx, *pzy, *pzz;
        mwSignedIndex m = dm[0]*dm[1]*dm[2];

        wx000 =  2*mu*(2*v0+v1+v2)/v0+2*lam + w000/v0;
        wx100 = -2*mu-lam + w100/v0;
        wx010 = -mu*v1/v0 + w010/v0;
        wx001 = -mu*v2/v0 + w001/v0;
        wy000 =  2*mu*(v0+2*v1+v2)/v1+2*lam + w000/v1;
        wy100 = -mu*v0/v1 + w100/v1;
        wy010 = -2*mu-lam + w010/v1;
        wy001 = -mu*v2/v1 + w001/v1;
        wz000 =  2*mu*(v0+v1+2*v2)/v2+2*lam + w000/v2;
        wz100 = -mu*v0/v2 + w100/v2;
        wz010 = -mu*v1/v2 + w010/v2;
        wz001 = -2*mu-lam + w001/v2;
        w2    = 0.25*mu+0.25*lam;

        pxx = f;
        pyx = f + m;
        pzx = f + m*2;
        pxy = f + m*3;
        pyy = f + m*4;
        pzy = f + m*5;
        pxz = f + m*6;
        pyz = f + m*7;
        pzz = f + m*8;

        pxx[0          ] += wx000;
        pxx[im1        ] += wx100;   pxx[ip1        ] += wx100;
        pxx[    jm1    ] += wx010;   pxx[    jp1    ] += wx010;
        pxx[        km1] += wx001;   pxx[        kp1] += wx001;
        pxy[ip1+jm1    ] += w2;      pxy[ip1+jp1    ] -= w2;      pxy[im1+jm1    ] -= w2;      pxy[im1+jp1    ] += w2;
        pxz[ip1    +km1] += w2;      pxz[ip1    +kp1] -= w2;      pxz[im1    +km1] -= w2;      pxz[im1    +kp1] += w2;
        pxx[im1+jm1    ] += w110/v0; pxx[ip1+jm1    ] += w110/v0; pxx[im1+jp1    ] += w110/v0; pxx[ip1+jp1    ] += w110/v0;
        pxx[im1    +km1] += w101/v0; pxx[ip1    +km1] += w101/v0; pxx[im1    +kp1] += w101/v0; pxx[ip1    +kp1] += w101/v0;
        pxx[    jm1+km1] += w011/v0; pxx[    jp1+km1] += w011/v0; pxx[    jm1+kp1] += w011/v0; pxx[    jp1+kp1] += w011/v0;
        pxx[im2        ] += w200/v0; pxx[ip2        ] += w200/v0;
        pxx[    jm2    ] += w020/v0; pxx[    jp2    ] += w020/v0;
        pxx[        km2] += w002/v0; pxx[        kp2] += w002/v0;

        pyy[0          ] += wy000;
        pyy[im1        ] += wy100;   pyy[ip1        ] += wy100;
        pyy[    jm1    ] += wy010;   pyy[    jp1    ] += wy010;
        pyy[        km1] += wy001;   pyy[        kp1] += wy001;
        pyx[ip1+jm1    ] += w2;      pyx[ip1+jp1    ] -= w2;      pyx[im1+jm1    ] -= w2;      pyx[im1+jp1    ] += w2;
        pyz[    jp1+km1] += w2;      pyz[    jp1+kp1] -= w2;      pyz[    jm1+km1] -= w2;      pyz[    jm1+kp1] += w2;
        pyy[im1+jm1    ] += w110/v1; pyy[ip1+jm1    ] += w110/v1; pyy[im1+jp1    ] += w110/v1; pyy[ip1+jp1    ] += w110/v1;
        pyy[im1    +km1] += w101/v1; pyy[ip1    +km1] += w101/v1; pyy[im1    +kp1] += w101/v1; pyy[ip1    +kp1] += w101/v1;
        pyy[    jm1+km1] += w011/v1; pyy[    jp1+km1] += w011/v1; pyy[    jm1+kp1] += w011/v1; pyy[    jp1+kp1] += w011/v1;
        pyy[im2        ] += w200/v1; pyy[ip2        ] += w200/v1;
        pyy[    jm2    ] += w020/v1; pyy[    jp2    ] += w020/v1;
        pyy[        km2] += w002/v1; pyy[        kp2] += w002/v1;

        pzz[0          ] += wz000;
        pzz[im1        ] += wz100;   pzz[ip1        ] += wz100;
        pzz[    jm1    ] += wz010;   pzz[    jp1    ] += wz010;
        pzz[        km1] += wz001;   pzz[        kp1] += wz001;
        pzx[im1    +kp1] += w2;      pzx[ip1    +kp1] -= w2;      pzx[im1    +km1] -= w2;      pzx[ip1    +km1] += w2;
        pzy[jm1    +kp1] += w2;      pzy[    jp1+kp1] -= w2;      pzy[    jm1+km1] -= w2;      pzy[    jp1+km1] += w2;
        pzz[im1+jm1    ] += w110/v2; pzz[ip1+jm1    ] += w110/v2; pzz[im1+jp1    ] += w110/v2; pzz[ip1+jp1    ] += w110/v2;
        pzz[im1    +km1] += w101/v2; pzz[ip1    +km1] += w101/v2; pzz[im1    +kp1] += w101/v2; pzz[ip1    +kp1] += w101/v2;
        pzz[    jm1+km1] += w011/v2; pzz[    jp1+km1] += w011/v2; pzz[    jm1+kp1] += w011/v2; pzz[    jp1+kp1] += w011/v2;
        pzz[im2        ] += w200/v2; pzz[ip2        ] += w200/v2;
        pzz[    jm2    ] += w020/v2; pzz[    jp2    ] += w020/v2;
        pzz[        km2] += w002/v2; pzz[        kp2] += w002/v2;
    }
}

/* Returns the sum of squares of (b - (L + H) * u)
 *
 * \param[in]    a      (Optional) Symmetric tensor field H (i.e. a 3x3   
 *                      symmetric matrix at each point of the lattice.).
 * \param[in]    b      Point at which to solve the system.
 * \param[in]    s[0:2] Voxel size of the lattice (voxels/mm)
 * \param[in]    s[3]   Parameter of the absolute displacement penalty
 * \param[in]    s[4]   Parameter of the membrane energy (penalizes 
 *                      elements of the Jacobian matrix -> 1st order 
 *                      smoothness)
 * \param[in]    s[5]   Parameter of the bending energy (penalizes 
 *                      elements of the Hessian matrix -> 2nd order 
 *                      smoothness)
 * \param[in]    s[6]   Parameter of the linear elastic energy (penalizes 
 *                      elements of the symmetric part of the Jacobian  
 *                      matrix -> penalizes scaling and shearing)
 * \param[in]    s[7]   Parameter of the linear elastic energy (penalizes 
 *                      the divergence of the Jacobian matrix -> preserves 
 *                      volumes)
 * \param[in]    u      Approximate solution.
 */
double sumsq(mwSize dm[], float a[], float b[], double s[], float u[])
{
    double w000,w100,w200,
           w010,w110,
           w020,
           w001,w101,
           w011,
           w002;
    double ss = 0.0;
    mwSignedIndex k;
    double v0 = s[0]*s[0], v1 = s[1]*s[1], v2 = s[2]*s[2];
    double lam0 = s[3], lam1 = s[4], lam2 = s[5], mu = s[6], lam = s[7];
    double wx000, wx100, wx010, wx001, wy000, wy100, wy010, wy001, wz000, wz100, wz010, wz001, w2;

    w000 = lam2*(6*(v0*v0+v1*v1+v2*v2) +8*(v0*v1+v0*v2+v1*v2)) +lam1*2*(v0+v1+v2) + lam0;
    w100 = lam2*(-4*v0*(v0+v1+v2)) -lam1*v0;
    w010 = lam2*(-4*v1*(v0+v1+v2)) -lam1*v1;
    w001 = lam2*(-4*v2*(v0+v1+v2)) -lam1*v2;
    w200 = lam2*v0*v0;
    w020 = lam2*v1*v1;
    w002 = lam2*v2*v2;
    w110 = lam2*2*v0*v1;
    w101 = lam2*2*v0*v2;
    w011 = lam2*2*v1*v2;

    wx000 =  2*mu*(2*v0+v1+v2)/v0+2*lam + w000/v0;
    wx100 = -2*mu-lam + w100/v0;
    wx010 = -mu*v1/v0 + w010/v0;
    wx001 = -mu*v2/v0 + w001/v0;
    wy000 =  2*mu*(v0+2*v1+v2)/v1+2*lam + w000/v1;
    wy100 = -mu*v0/v1 + w100/v1;
    wy010 = -2*mu-lam + w010/v1;
    wy001 = -mu*v2/v1 + w001/v1;
    wz000 =  2*mu*(v0+v1+2*v2)/v2+2*lam + w000/v2;
    wz100 = -mu*v0/v2 + w100/v2;
    wz010 = -mu*v1/v2 + w010/v2;
    wz001 = -2*mu-lam + w001/v2;
    w2    = 0.25*mu+0.25*lam;

    for(k=0; k<dm[2]; k++)
    {
        mwSignedIndex j, km2,km1,kp1,kp2;
        int fkm2, fkm1, fkp2, fkp1, gkm2, gkm1, gkp2, gkp1;
        
        km2 = (bound(k-2,dm[2])-k)*dm[0]*dm[1];
        km1 = (bound(k-1,dm[2])-k)*dm[0]*dm[1];
        kp1 = (bound(k+1,dm[2])-k)*dm[0]*dm[1];
        kp2 = (bound(k+2,dm[2])-k)*dm[0]*dm[1];
        fkm1 = bound_factor(k-1, dm[2], 1);
        fkm2 = bound_factor(k-2, dm[2], 1);
        fkp1 = bound_factor(k+1, dm[2], 1);
        fkp2 = bound_factor(k+2, dm[2], 1);
        gkm1 = bound_factor(k-1, dm[2], 0);
        gkm2 = bound_factor(k-2, dm[2], 0);
        gkp1 = bound_factor(k+1, dm[2], 0);
        gkp2 = bound_factor(k+2, dm[2], 0);

        for(j=0; j<dm[1]; j++)
        {
            float *pux, *puy, *puz, *pbx, *pby, *pbz, *paxx, *payy, *pazz, *paxy, *paxz, *payz;
            mwSignedIndex i, jm2,jm1,jp1,jp2;
            int fjm2, fjm1, fjp2, fjp1, gjm2, gjm1, gjp2, gjp1;

            pux  = u+dm[0]*(j+dm[1]*k);
            puy  = u+dm[0]*(j+dm[1]*(k+dm[2]));
            puz  = u+dm[0]*(j+dm[1]*(k+dm[2]*2));
            pbx  = b+dm[0]*(j+dm[1]*k);
            pby  = b+dm[0]*(j+dm[1]*(k+dm[2]));
            pbz  = b+dm[0]*(j+dm[1]*(k+dm[2]*2));

            if (a)
            {
                paxx = a+dm[0]*(j+dm[1]*k);
                payy = a+dm[0]*(j+dm[1]*(k+dm[2]));
                pazz = a+dm[0]*(j+dm[1]*(k+dm[2]*2));
                paxy = a+dm[0]*(j+dm[1]*(k+dm[2]*3));
                paxz = a+dm[0]*(j+dm[1]*(k+dm[2]*4));
                payz = a+dm[0]*(j+dm[1]*(k+dm[2]*5));
            }

            // Wrap/mirror outside domain
            jm2 = (bound(j-2,dm[1])-j)*dm[0];
            jm1 = (bound(j-1,dm[1])-j)*dm[0];
            jp1 = (bound(j+1,dm[1])-j)*dm[0];
            jp2 = (bound(j+2,dm[1])-j)*dm[0];
            // Factor if dirichlet/sliding
            // 1) Diagonal convolution
            fjm1 = bound_factor(j-1, dm[1], 1);
            fjm2 = bound_factor(j-2, dm[1], 1);
            fjp1 = bound_factor(j+1, dm[1], 1);
            fjp2 = bound_factor(j+2, dm[1], 1);
            // 2) Off-Diagonal convolution
            gjm1 = bound_factor(j-1, dm[1], 0);
            gjm2 = bound_factor(j-2, dm[1], 0);
            gjp1 = bound_factor(j+1, dm[1], 0);
            gjp2 = bound_factor(j+2, dm[1], 0);

            for(i=0; i<dm[0]; i++)
            {
                mwSignedIndex im2,im1,ip1,ip2;
                int fim1, fim2, fip1, fip2, gim1, gim2, gip1, gip2;
                float *px = pux+i, *py = puy+i, *pz = puz+i;
                double tmp, abx, aby, abz, c;

                im2 = bound(i-2,dm[0])-i;
                im1 = bound(i-1,dm[0])-i;
                ip1 = bound(i+1,dm[0])-i;
                ip2 = bound(i+2,dm[0])-i;
                fim1 = bound_factor(i-1, dm[0], 1);
                fim2 = bound_factor(i-2, dm[0], 1);
                fip1 = bound_factor(i+1, dm[0], 1);
                fip2 = bound_factor(i+2, dm[0], 1);
                gim1 = bound_factor(i-1, dm[0], 0);
                gim2 = bound_factor(i-2, dm[0], 0);
                gip1 = bound_factor(i+1, dm[0], 0);
                gip2 = bound_factor(i+2, dm[0], 0);

                if (a)
                {
                    abx = paxx[i]*px[0] + paxy[i]*py[0] + paxz[i]*pz[0];
                    aby = paxy[i]*px[0] + payy[i]*py[0] + payz[i]*pz[0];
                    abz = paxz[i]*px[0] + payz[i]*py[0] + pazz[i]*pz[0];
                }
                else
                {
                    abx = aby = abz = 0.0;
                }

                /* Note that a few things have been done here to reduce rounding errors.
                   This may slow things down, but it does lead to more accuracy. */
                c   = px[0];
                tmp = abx - pbx[i]
                        + wx100*((fim1*px[im1        ]-c) + (fip1*px[ip1        ]-c))
                        + wx010*((gjm1*px[    jm1    ]-c) + (gjp1*px[    jp1    ]-c))
                        + wx001*((gkm1*px[        km1]-c) + (gkp1*px[        kp1]-c))
                        + w2   *( gip1*fjm1*py[ip1+jm1] - gip1*fjp1*py[ip1+jp1] +
                                  gim1*fjp1*py[im1+jp1] - gim1*fjm1*py[im1+jm1] + 
                                  gip1*fkm1*pz[ip1+km1] - gip1*fkp1*pz[ip1+kp1] + 
                                  gim1*fkp1*pz[im1+kp1] - gim1*fkm1*pz[im1+km1] )
                        + (lam0*c
                        +  w110*((fim1*gjm1*px[im1+jm1    ]-c) + (fip1*gjm1*px[ip1+jm1    ]-c) + (fim1*gjp1*px[im1+jp1    ]-c) + (fip1*gjp1*px[ip1+jp1    ]-c))
                        +  w101*((fim1*gkm1*px[im1    +km1]-c) + (fip1*gkm1*px[ip1    +km1]-c) + (fim1*gkp1*px[im1    +kp1]-c) + (fip1*gkp1*px[ip1    +kp1]-c))
                        +  w011*((gjm1*gkm1*px[    jm1+km1]-c) + (gjp1*gkm1*px[    jp1+km1]-c) + (gjm1*gkp1*px[    jm1+kp1]-c) + (gjp1*gkp1*px[    jp1+kp1]-c))
                        +  w200*((fim2*px[im2        ]-c) + (fip2*px[ip2        ]-c))
                        +  w020*((gjm2*px[    jm2    ]-c) + (gjp2*px[    jp2    ]-c))
                        +  w002*((gkm2*px[        km2]-c) + (gkp2*px[        kp2]-c)))/v0;
                ss += tmp*tmp;

                c   = py[0];
                tmp = aby - pby[i]
                        + wy100*((gim1*py[im1        ]-c) + (gip1*py[ip1        ]-c))
                        + wy010*((fjm1*py[    jm1    ]-c) + (fjp1*py[    jp1    ]-c))   
                        + wy001*((gkm1*py[        km1]-c) + (gkp1*py[        kp1]-c))   
                        + w2   *( gjp1*fim1*px[jp1+im1] - gjp1*fip1*px[jp1+ip1] + 
                                  gjm1*fip1*px[jm1+ip1] - gjm1*fim1*px[jm1+im1] + 
                                  gjp1*fkm1*pz[jp1+km1] - gjp1*fkp1*pz[jp1+kp1] + 
                                  gjm1*fkp1*pz[jm1+kp1] - gjm1*fkm1*pz[jm1+km1] )
                        + (lam0*c
                        +  w110*((gim1*fjm1*py[im1+jm1    ]-c) + (gip1*fjm1*py[ip1+jm1    ]-c) + (gim1*fjp1*py[im1+jp1    ]-c) + (gip1*fjp1*py[ip1+jp1    ]-c))
                        +  w101*((gim1*gkm1*py[im1    +km1]-c) + (gip1*gkm1*py[ip1    +km1]-c) + (gim1*gkp1*py[im1    +kp1]-c) + (gip1*gkp1*py[ip1    +kp1]-c))
                        +  w011*((fjm1*gkm1*py[    jm1+km1]-c) + (fjp1*gkm1*py[    jp1+km1]-c) + (fjm1*gkp1*py[    jm1+kp1]-c) + (fjp1*gkp1*py[    jp1+kp1]-c))
                        +  w200*((gim2*py[im2        ]-c) + (gip2*py[ip2        ]-c))
                        +  w020*((fjm2*py[    jm2    ]-c) + (fjp2*py[    jp2    ]-c))
                        +  w002*((gkm2*py[        km2]-c) + (gkp2*py[        kp2]-c)))/v1;
                ss += tmp*tmp;

                c   = pz[0];
                tmp = abz - pbz[i]
                        + wz100*((gim1*pz[im1        ]-c) + (gip1*pz[ip1        ]-c))
                        + wz010*((gjm1*pz[    jm1    ]-c) + (gjp1*pz[    jp1    ]-c))   
                        + wz001*((fkm1*pz[        km1]-c) + (fkp1*pz[        kp1]-c))   
                        + w2   *( gkp1*fim1*px[kp1+im1] - gkp1*fip1*px[kp1+ip1] + 
                                  gkm1*fip1*px[km1+ip1] - gkm1*fim1*px[km1+im1] + 
                                  gkp1*fjm1*py[kp1+jm1] - gkp1*fjp1*py[kp1+jp1] + 
                                  gkm1*fjp1*py[km1+jp1] - gkm1*fjm1*py[km1+jm1] )
                        + (lam0*c
                        +  w110*((gim1*gjm1*pz[im1+jm1    ]-c) + (gip1*gjm1*pz[ip1+jm1    ]-c) + (gim1*gjp1*pz[im1+jp1    ]-c) + (gip1*gjp1*pz[ip1+jp1    ]-c))
                        +  w101*((gim1*fkm1*pz[im1    +km1]-c) + (gip1*fkm1*pz[ip1    +km1]-c) + (gim1*fkp1*pz[im1    +kp1]-c) + (gip1*fkp1*pz[ip1    +kp1]-c))
                        +  w011*((gjm1*fkm1*pz[    jm1+km1]-c) + (gjp1*fkm1*pz[    jp1+km1]-c) + (gjm1*fkp1*pz[    jm1+kp1]-c) + (gjp1*fkp1*pz[    jp1+kp1]-c))
                        +  w200*((gim2*pz[im2        ]-c) + (gip1*pz[ip2        ]-c))
                        +  w020*((gjm2*pz[    jm2    ]-c) + (gjp2*pz[    jp2    ]-c))
                        +  w002*((fkm2*pz[        km2]-c) + (fkp2*pz[        kp2]-c)))/v2;
                ss += tmp*tmp;
            }
        }
    }
    return(ss);
}

/* Pointwise matrix multiplication.
 * Returns A[i] * p[i], where A contains M symmetric 3*3 matrix 
 * (tensor field) and p contains M 3*1 vectors (velocity field).
 * 
 * \param[in]  dm  Dimension of the lattice (M = prod(dm))
 * \param[in]  A   Symmetric tensor field
 * \param[in]  p   Vector field
 * \param[out] Ap  Result of the pointwise matrix multiplication (A * p)
 */
static void Atimesp1(mwSize dm[], float A[], float p[], float Ap[])
{
    mwSignedIndex i, m = dm[0]*dm[1]*dm[2];
    float *pa11 = A ,     *pa22 = A +m,   *pa33 = A +2*m,
          *pa12 = A +3*m, *pa13 = A +4*m, *pa23 = A +5*m;
    float *pap1 = Ap,     *pap2 = Ap+m,   *pap3 = Ap+2*m;
    float *pp1  = p,      *pp2  = p +m,   *pp3  = p +2*m;

    if (A==0) return;

    for(i=0; i<m; i++)
    {
        pap1[i] += pa11[i]*pp1[i] + pa12[i]*pp2[i] + pa13[i]*pp3[i];
        pap2[i] += pa12[i]*pp1[i] + pa22[i]*pp2[i] + pa23[i]*pp3[i];
        pap3[i] += pa13[i]*pp1[i] + pa23[i]*pp2[i] + pa33[i]*pp3[i];
    }
}

/* Returns u = L * v in the case where L is based on the linear 
 * elastic energy model
 *
 * \param dm     Dimensions of the velocity lattice (voxels)
 * \param f      Velocities computed on the above lattice
 * \param s[0:3] Voxel sizes (voxels/mm)
 * \param s[3]   Parameter of the absolute displacement penalty (lam0)
 * \param s[6]   Parameter of the linear elastic energy (mu) (penalizes 
 *               elements of the symmetric part of the Jacobian matrix 
 *               -> penalizes scaling and shearing)
 * \param s[7]   Parameter of the linear elastic energy (lam) (penalizes 
 *               the divergence of the Jacobian matrix -> preserves 
 *               volumes)
 * \param g      Allocated array in which to store the output momentum 
 */
void vel2mom_le(mwSize dm[], float f[], double s[], float g[])
{
    mwSignedIndex k;
    double wx000, wx100, wx010, wx001, wy000, wy100, wy010, wy001, wz000, wz100, wz010, wz001, w2;
    double v0 = s[0]*s[0], v1 = s[1]*s[1], v2 = s[2]*s[2];
    double lam0 = s[3], mu = s[6], lam = s[7];

    // Convolution kernel
    wx000 =  2*mu*(2*v0+v1+v2)/v0+2*lam + lam0/v0;
    wx100 = -2*mu-lam;
    wx010 = -mu*v1/v0;
    wx001 = -mu*v2/v0;
    wy000 =  2*mu*(v0+2*v1+v2)/v1+2*lam + lam0/v1;
    wy100 = -mu*v0/v1;
    wy010 = -2*mu-lam;
    wy001 = -mu*v2/v1;
    wz000 =  2*mu*(v0+v1+2*v2)/v2+2*lam + lam0/v2;
    wz100 = -mu*v0/v2;
    wz010 = -mu*v1/v2;
    wz001 = -2*mu-lam;
    w2    = 0.25*mu+0.25*lam;

    for(k=0; k<dm[2]; k++)
    {
        mwSignedIndex j, km1, kp1;
        int fkm1, fkp1, gkm1, gkp1;
        km1 = (bound(k-1,dm[2])-k)*dm[0]*dm[1]; // Wrap/mirror outside domain
        kp1 = (bound(k+1,dm[2])-k)*dm[0]*dm[1];
        fkm1 = bound_factor(k-1, dm[2], 1); // Factor if sliding/dircihlet
        fkp1 = bound_factor(k+1, dm[2], 1);
        gkm1 = bound_factor(k-1, dm[2], 0);
        gkp1 = bound_factor(k+1, dm[2], 0);

        for(j=0; j<dm[1]; j++)
        {
            mwSignedIndex i, jm1, jp1;
            int fjm1, fjp1, gjm1, gjp1;
            float *pgx, *pgy, *pgz, *pfx, *pfy, *pfz;

            // Pointers to output momentum volume
            pgx = g+dm[0]*(j+dm[1]*k);              //< ux(:,j,k)
            pgy = g+dm[0]*(j+dm[1]*(k+dm[2]));      //< uy(:,j,k)
            pgz = g+dm[0]*(j+dm[1]*(k+dm[2]*2));    //< uz(:,j,k)

            // Pointers to input velocity volume
            pfx = f+dm[0]*(j+dm[1]*k);              //< vx(:,j,k)
            pfy = f+dm[0]*(j+dm[1]*(k+dm[2]));      //< vy(:,j,k)
            pfz = f+dm[0]*(j+dm[1]*(k+dm[2]*2));    //< vz(:,j,k)

            jm1 = (bound(j-1,dm[1])-j)*dm[0];
            jp1 = (bound(j+1,dm[1])-j)*dm[0];
            fjm1 = bound_factor(j-1, dm[1], 1);
            fjp1 = bound_factor(j+1, dm[1], 1);
            gjm1 = bound_factor(j-1, dm[1], 0);
            gjp1 = bound_factor(j+1, dm[1], 0);

            for(i=0; i<dm[0]; i++)
            {
                mwSignedIndex im1, ip1;
                int fim1, fip1, gim1, gip1;
                float *px = &pfx[i],    //< vx(i,j,k)
                      *py = &pfy[i],    //< vy(i,j,k)
                      *pz = &pfz[i];    //< vz(i,j,k)
                double c;

                im1 = bound(i-1,dm[0])-i;
                ip1 = bound(i+1,dm[0])-i;
                fim1 = bound_factor(i-1, dm[0], 1);
                fip1 = bound_factor(i+1, dm[0], 1);
                gim1 = bound_factor(i-1, dm[0], 0);
                gip1 = bound_factor(i+1, dm[0], 0);

                /* Note that a few things have been done here to reduce rounding errors.
                   This may slow things down, but it does lead to more accuracy. */
                
                // ux = lxx (*) vx + lxy (*) vy + lxz (*) vz
                c      = px[0];
                pgx[i] =  lam0/v0*c
                       + wx100*((fip1*px[ip1]-c) + (fim1*px[im1]-c))
                       + wx010*((gjp1*px[jp1]-c) + (gjm1*px[jm1]-c))
                       + wx001*((gkp1*px[kp1]-c) + (gkm1*px[km1]-c))
                       + w2   *(gip1*fjm1*py[ip1+jm1] - gip1*fjp1*py[ip1+jp1] + 
                                gim1*fjp1*py[im1+jp1] - gim1*fjm1*py[im1+jm1] + 
                                gip1*fkm1*pz[ip1+km1] - gip1*fkp1*pz[ip1+kp1] + 
                                gim1*fkp1*pz[im1+kp1] - gim1*fkm1*pz[im1+km1]);

                // uy = lyx (*) vx + lyy (*) vy + lyz (*) vz
                c      = py[0];
                pgy[i] =  lam0/v1*c
                       + wy100*((gip1*py[ip1]-c) + (gim1*py[im1]-c))
                       + wy010*((fjp1*py[jp1]-c) + (fjm1*py[jm1]-c))
                       + wy001*((gkp1*py[kp1]-c) + (gkm1*py[km1]-c))
                       + w2   *(gjp1*fim1*px[jp1+im1] - gjp1*fip1*px[jp1+ip1] + 
                                gjm1*fip1*px[jm1+ip1] - gjm1*fim1*px[jm1+im1] + 
                                gjp1*fkm1*pz[jp1+km1] - gjp1*fkp1*pz[jp1+kp1] + 
                                gjm1*fkp1*pz[jm1+kp1] - gjm1*fkm1*pz[jm1+km1]);

                // uz = lzx (*) vx + lzy (*) vy + lzz (*) vz
                c      = pz[0];
                pgz[i] = lam0/v2*c
                       + wz100*((gip1*pz[ip1]-c) + (gim1*pz[im1]-c))
                       + wz010*((gjp1*pz[jp1]-c) + (gjm1*pz[jm1]-c))
                       + wz001*((fkp1*pz[kp1]-c) + (fkm1*pz[km1]-c))
                       + w2   *(gkp1*fim1*px[kp1+im1] - gkp1*fip1*px[kp1+ip1] + 
                                gkm1*fip1*px[km1+ip1] - gkm1*fim1*px[km1+im1] + 
                                gkp1*fjm1*py[kp1+jm1] - gkp1*fjp1*py[kp1+jp1] + 
                                gkm1*fjp1*py[km1+jp1] - gkm1*fjm1*py[km1+jm1]);
            }
        }
    }
}

/* Relaxation iterations for linear energy penalty.
 *
 * . We solve for (L + H) * u = b
 * . F = nondiag(L) (i.e., L without its diagonal elements)
 * . E = H + diag(L) + sI (to insure diagonal dominance)
 * . u = E^{-1} * ( b - F * u )
 *
 * \param[in]    dm     Dimension of the lattice.
 * \param[in]    a      Symmetric tensor field H (i.e. a 3x3 symmetric  
 *                      matrix at each point of the lattice.).
 * \param[in]    b      Point at which to solve the system.
 * \param[in]    s[0:2] Voxel size of the lattice (voxels/mm)
 * \param[in]    s[3]   Parameter of the absolute displacement penalty
 * \param[in]    s[4]   Parameter of the membrane energy (penalizes 
 *                      elements of the Jacobian matrix -> 1st order 
 *                      smoothness)
 * \param[in]    s[5]   Parameter of the bending energy (penalizes 
 *                      elements of the Hessian matrix -> 2nd order 
 *                      smoothness)
 * \param[in]    s[6]   Parameter of the linear elastic energy (penalizes 
 *                      elements of the symmetric part of the Jacobian  
 *                      matrix -> penalizes scaling and shearing)
 * \param[in]    s[7]   Parameter of the linear elastic energy (penalizes 
 *                      the divergence of the Jacobian matrix -> preserves 
 *                      volumes)
 * \param[in]    nit    Number of relaxation iterations.
 * \param[inout] u      [in] Initial guess for the solution.
 *                      [out] Output relaxation solution.
 */
void relax_le(mwSize dm[], float a[], float b[], double s[], int nit, float u[])
{
    int it;
    double wx000, wx100, wx010, wx001, wy000, wy100, wy010, wy001, wz000, wz100, wz010, wz001, w2;
    double v0 = s[0]*s[0], v1 = s[1]*s[1], v2 = s[2]*s[2];
    double lam0 = s[3], mu = s[6], lam = s[7];

    // Direct convolution kernel (i.e., local L)
    wx000 =  2*mu*(2*v0+v1+v2)/v0+2*lam + lam0/v0;
    wx100 = -2*mu-lam;
    wx010 = -mu*v1/v0;
    wx001 = -mu*v2/v0;
    wy000 =  2*mu*(v0+2*v1+v2)/v1+2*lam + lam0/v1;
    wy100 = -mu*v0/v1;
    wy010 = -2*mu-lam;
    wy001 = -mu*v2/v1;
    wz000 =  2*mu*(v0+v1+2*v2)/v2+2*lam + lam0/v2;
    wz100 = -mu*v0/v2;
    wz010 = -mu*v1/v2;
    wz001 = -2*mu-lam;
    w2    = 0.25*mu+0.25*lam;

    // Insure diagonal dominance
    wx000 = wx000*1.00001 + 1e-6;
    wy000 = wy000*1.00001 + 1e-6;
    wz000 = wz000*1.00001 + 1e-6;

    // Deal with 2D cases
    if (dm[0]==1)
    {
        wx000 += 2*wx100 + 2*wy100 + 2*wz100;
        wx100  = 0.0;
        wy100  = 0.0;
        wz100  = 0.0;
    }
    if (dm[1]==1)
    {
        wx000 += 2*wx010 + 2*wy010 + 2*wz010;
        wx010  = 0.0;
        wy010  = 0.0;
        wz010  = 0.0;
    }
    if (dm[2]==1)
    {
        wx000 += 2*wx001 + 2*wy001 + 2*wz001;
        wx001  = 0.0;
        wy001  = 0.0;
        wz001  = 0.0;
    }

#   ifdef VERBOSE
        for(it=0; it< 10-(int)ceil(1.44269504088896*log((double)dm[0])); it++) printf("  ");
        printf("L%dx%dx%d: %g ", dm[0],dm[1],dm[2], sumsq(dm, a, b, s, u));
#   endif

    for(it=0; it<8*nit; it++) // Relaxation loop
    {
        mwSignedIndex k;
        for(k=it&1; k<dm[2]; k+=2) // Z loop
        {
            mwSignedIndex j, km1, kp1;
            int fkm1, fkp1, gkm1, gkp1;
            km1 = (bound(k-1,dm[2])-k)*dm[0]*dm[1];
            kp1 = (bound(k+1,dm[2])-k)*dm[0]*dm[1];
            fkm1 = bound_factor(k-1, dm[2], 1); // Factor if sliding
            fkp1 = bound_factor(k+1, dm[2], 1);
            gkm1 = bound_factor(k-1, dm[2], 0); // Factor if sliding
            gkp1 = bound_factor(k+1, dm[2], 0);

            for(j=(it>>1)&1; j<dm[1]; j+=2) // Y loop
            {
                float *pux, *puy, *puz, *pbx, *pby, *pbz, *paxx, *payy, *pazz, *paxy, *paxz, *payz;
                mwSignedIndex i, jm1,jp1;
                int fjm1, fjp1, gjm1, gjp1;

                // pu <- Get slice of previous guess of the solution (:,j,k)
                pux  = u+dm[0]*(j+dm[1]* k);
                puy  = u+dm[0]*(j+dm[1]*(k+dm[2]));
                puz  = u+dm[0]*(j+dm[1]*(k+dm[2]*2));
                // pb <- Get slice of point to solve for (:,j,k)
                pbx  = b+dm[0]*(j+dm[1]* k);
                pby  = b+dm[0]*(j+dm[1]*(k+dm[2]));
                pbz  = b+dm[0]*(j+dm[1]*(k+dm[2]*2));

                if (a)
                {
                    // pa <- Get slice of symmetric tensor field (:,j,k)
                    paxx = a+dm[0]*(j+dm[1]* k);
                    payy = a+dm[0]*(j+dm[1]*(k+dm[2]));
                    pazz = a+dm[0]*(j+dm[1]*(k+dm[2]*2));
                    paxy = a+dm[0]*(j+dm[1]*(k+dm[2]*3));
                    paxz = a+dm[0]*(j+dm[1]*(k+dm[2]*4));
                    payz = a+dm[0]*(j+dm[1]*(k+dm[2]*5));
                }

                jm1 = (bound(j-1,dm[1])-j)*dm[0];
                jp1 = (bound(j+1,dm[1])-j)*dm[0];
                fjm1 = bound_factor(j-1, dm[1], 1); // Factor if sliding
                fjp1 = bound_factor(j+1, dm[1], 1);
                fjm1 = bound_factor(j-1, dm[1], 0);
                fjp1 = bound_factor(j+1, dm[1], 0);

                for(i=(it>>2)&1; i<dm[0]; i+=2) // X loop
                {
                    mwSignedIndex im1,ip1;
                    int fim1, fip1, gim1, gip1;
                    double sux, suy, suz;
                    // p <- Get previous guess of the solution at (i,j,k)
                    float *px = pux+i, *py = puy+i, *pz = puz+i;

                    im1 = bound(i-1,dm[0])-i;
                    ip1 = bound(i+1,dm[0])-i;
                    fim1 = bound_factor(i-1, dm[0], 1); // Factor if sliding
                    fip1 = bound_factor(i+1, dm[0], 1);
                    fim1 = bound_factor(i-1, dm[0], 0);
                    fip1 = bound_factor(i+1, dm[0], 0);

                    // su <- Get relaxation point to solve for (b - F * u)
                    sux = pbx[i] - ( wx100*(fim1*px[im1] + fip1*px[ip1])
                                   + wx010*(gjm1*px[jm1] + gjp1*px[jp1])
                                   + wx001*(gkm1*px[km1] + gkp1*px[kp1])
                                   + w2   *(gip1*fjm1*py[ip1+jm1] - gip1*fjp1*py[ip1+jp1] + 
                                            gim1*fjp1*py[im1+jp1] - gim1*fjm1*py[im1+jm1] + 
                                            gip1*fkm1*pz[ip1+km1] - gip1*fkp1*pz[ip1+kp1] + 
                                            gim1*fkp1*pz[im1+kp1] - gim1*fkm1*pz[im1+km1]));

                    suy = pby[i] - ( wy100*(gim1*py[im1] + gip1*py[ip1])
                                   + wy010*(fjm1*py[jm1] + fjp1*py[jp1])
                                   + wy001*(gkm1*py[km1] + gkp1*py[kp1])
                                   + w2   *(gjp1*fim1*px[jp1+im1] - gjp1*fip1*px[jp1+ip1] + 
                                            gjm1*fip1*px[jm1+ip1] - gjm1*fim1*px[jm1+im1] + 
                                            gjp1*fkm1*pz[jp1+km1] - gjp1*fkp1*pz[jp1+kp1] + 
                                            gjm1*fkp1*pz[jm1+kp1] - gjm1*fkm1*pz[jm1+km1]));

                    suz = pbz[i] - ( wz100*(gim1*pz[im1] + gip1*pz[ip1])
                                   + wz010*(gjm1*pz[jm1] + gjp1*pz[jp1])
                                   + wz001*(fkm1*pz[km1] + fkp1*pz[kp1])
                                   + w2   *(gkp1*fim1*px[kp1+im1] - gkp1*fip1*px[kp1+ip1] + 
                                            gkm1*fip1*px[km1+ip1] - gkm1*fim1*px[km1+im1] + 
                                            gkp1*fjm1*py[kp1+jm1] - gkp1*fjp1*py[kp1+jp1] + 
                                            gkm1*fjp1*py[km1+jp1] - gkm1*fjm1*py[km1+jm1]));

                    if (a)
                    {
                        double axx, ayy, azz, axy, axz, ayz, idt;

                        // Get 3x3 symmetric matrix to invert
                        // a <- (E = H + diag(L) + sI)
                        axx  = paxx[i] + wx000;
                        ayy  = payy[i] + wy000;
                        azz  = pazz[i] + wz000;
                        axy  = paxy[i];
                        axz  = paxz[i];
                        ayz  = payz[i];
                        idt  = 1.0/(axx*ayy*azz -axx*ayz*ayz-ayy*axz*axz-azz*axy*axy +2*axy*axz*ayz);

                        // Update solution:
                        // p <- E^{-1} * su = E^{-1} * (b - F * u)
                        // with E = H + diag(L) + sI
                        *px = idt*(sux*(ayy*azz-ayz*ayz)+suy*(axz*ayz-axy*azz)+suz*(axy*ayz-axz*ayy));
                        *py = idt*(sux*(axz*ayz-axy*azz)+suy*(axx*azz-axz*axz)+suz*(axy*axz-axx*ayz));
                        *pz = idt*(sux*(axy*ayz-axz*ayy)+suy*(axy*axz-axx*ayz)+suz*(axx*ayy-axy*axy));
                    }
                    else
                    {
                        // Update solution:
                        // p <- E^{-1} * su = E^{-1} * (b - F * u)
                        // with E = diag(L) + sI
                        *px = sux/wx000;
                        *py = suy/wy000;
                        *pz = suz/wz000;
                    }
                }
            }
        } 
#       ifdef VERBOSE
            if ((it%8)==7) printf(" %g", sumsq(dm, a, b, s, u));
#       endif
    }
#   ifdef VERBOSE
        printf("\n");
#   endif
}


void Atimesp_le(mwSize dm[], float A[], double param[], float p[], float Ap[])
{
    vel2mom_le(dm, p, param, Ap);
    Atimesp1(dm, A, p, Ap);
}

/* Returns m = L * v in the case where L is based on the membrane 
 * energy model.
 *
 * \param dm     Dimensions of the velocity lattice (voxels)
 * \param f      Velocities computed on the above lattice
 * \param s[0:3] Voxel sizes (voxels/mm)
 * \param s[3]   Parameter of the absolute displacement penalty (lam0)
 * \param s[4]   Parameter of the membrane energy (lam1) (penalizes 
 *               elements of the Jacobian matrix -> 1st order smoothness)
 * \param g      Allocated array in which to store the output momentum 
 */
void vel2mom_me(mwSize dm[], float f[], double s[], float g[])
{
    mwSignedIndex i, j, k, km1,kp1, jm1,jp1, im1,ip1;
    int fkm1, fkp1, fjm1, fjp1, fim1, fip1, 
        gkm1, gkp1, gjm1, gjp1, gim1, gip1;
    float *pgx, *pgy, *pgz, *pfx, *pfy, *pfz;
    double w000,w001,w010,w100;
    double lam0 = s[3], lam1 = s[4];

    w000 = lam1*(2*s[0]*s[0]+2*s[1]*s[1]+2*s[2]*s[2]) + lam0;
    w001 = lam1*(-s[2]*s[2]);
    w010 = lam1*(-s[1]*s[1]);
    w100 = lam1*(-s[0]*s[0]);

    for(k=0; k<dm[2]; k++)
    {
        km1 = (bound(k-1,dm[2])-k)*dm[0]*dm[1];
        kp1 = (bound(k+1,dm[2])-k)*dm[0]*dm[1];
        fkm1 = bound_factor(k-1, dm[2], 1); // Factor if sliding
        fkp1 = bound_factor(k+1, dm[2], 1);
        gkm1 = bound_factor(k-1, dm[2], 0);
        gkp1 = bound_factor(k+1, dm[2], 0);

        for(j=0; j<dm[1]; j++)
        {
            pgx = g+dm[0]*(j+dm[1]*k);
            pgy = g+dm[0]*(j+dm[1]*(k+dm[2]));
            pgz = g+dm[0]*(j+dm[1]*(k+dm[2]*2));

            pfx = f+dm[0]*(j+dm[1]*k);
            pfy = f+dm[0]*(j+dm[1]*(k+dm[2]));
            pfz = f+dm[0]*(j+dm[1]*(k+dm[2]*2));

            jm1 = (bound(j-1,dm[1])-j)*dm[0];
            jp1 = (bound(j+1,dm[1])-j)*dm[0];
            fjm1 = bound_factor(j-1, dm[1], 1); // Factor if sliding
            fjp1 = bound_factor(j+1, dm[1], 1);
            gjm1 = bound_factor(j-1, dm[1], 0);
            gjp1 = bound_factor(j+1, dm[1], 0);

            for(i=0; i<dm[0]; i++)
            {
                float *px = &pfx[i], *py = &pfy[i], *pz = &pfz[i];

                im1 = bound(i-1,dm[0])-i;
                ip1 = bound(i+1,dm[0])-i;
                fim1 = bound_factor(i-1, dm[0], 1); // Factor if sliding
                fip1 = bound_factor(i+1, dm[0], 1);
                gim1 = bound_factor(i-1, dm[0], 0);
                gip1 = bound_factor(i+1, dm[0], 0);

                pgx[i] = (w000*px[0] + w001*(gkm1*px[km1] + gkp1*px[kp1]) +
                                       w010*(gjm1*px[jm1] + gjp1*px[jp1]) + 
                                       w100*(fim1*px[im1] + fip1*px[ip1]) )/(s[0]*s[0]);
                pgy[i] = (w000*py[0] + w001*(gkm1*py[km1] + gkp1*py[kp1]) + 
                                       w010*(fjm1*py[jm1] + fjp1*py[jp1]) + 
                                       w100*(gim1*py[im1] + gip1*py[ip1]) )/(s[1]*s[1]);
                pgz[i] = (w000*pz[0] + w001*(fkm1*pz[km1] + fkp1*pz[kp1]) + 
                                       w010*(gjm1*pz[jm1] + gjp1*pz[jp1]) + 
                                       w100*(gim1*pz[im1] + gip1*pz[ip1]) )/(s[2]*s[2]);
            }
        }
    }
}

void relax_me(mwSize dm[], float a[], float b[], double s[], int nit, float u[])
{
    int it;
    double w000,w001,w010,w100;
    double lam0 = s[3], lam1 = s[4];

    w000 = lam1*(2*s[0]*s[0]+2*s[1]*s[1]+2*s[2]*s[2]) + lam0;
    w001 = lam1*(-s[2]*s[2]);
    w010 = lam1*(-s[1]*s[1]);
    w100 = lam1*(-s[0]*s[0]);

    w000 = w000*1.00001 + 1e-6;

    if (dm[0]==1)
    {
        w000 += 2*w100;
        w100  = 0.0;
    }
    if (dm[1]==1)
    {
        w000 += 2*w010;
        w010  = 0.0;
    }
    if (dm[2]==1)
    {
        w000 += 2*w001;
        w001  = 0.0;
    }

#   ifdef VERBOSE
        for(it=0; it< 10-(int)ceil(1.44269504088896*log((double)dm[0])); it++) printf("  ");
        printf("M%dx%dx%d: %g ", dm[0],dm[1],dm[2], sumsq(dm, a, b, s, u));
#   endif

    for(it=0; it<2*nit; it++)
    {
        mwSignedIndex k, kstart;
        mwSignedIndex j, jstart;
        mwSignedIndex i, istart;

        kstart = it%2;
        for(k=0; k<dm[2]; k++)
        {
            mwSignedIndex km1, kp1;
            int fkm1, fkp1, gkm1, gkp1;
            km1 = (bound(k-1,dm[2])-k)*dm[0]*dm[1];
            kp1 = (bound(k+1,dm[2])-k)*dm[0]*dm[1];
            fkm1 = bound_factor(k-1, dm[2], 1); // Factor if sliding
            fkp1 = bound_factor(k+1, dm[2], 1);
            gkm1 = bound_factor(k-1, dm[2], 0);
            gkp1 = bound_factor(k+1, dm[2], 0);

            jstart = (kstart == (k%2));
            for(j=0; j<dm[1]; j++)
            {
                float *pux, *puy, *puz, *pbx, *pby, *pbz, *paxx, *paxy, *payy, *paxz, *payz, *pazz;
                mwSignedIndex jm1,jp1, im1,ip1;
                int fjm1, fjp1, fim1, fip1, gjm1, gjp1, gim1, gip1;

                pux  = u+dm[0]*(j+dm[1]*k);
                puy  = u+dm[0]*(j+dm[1]*(k+dm[2]));
                puz  = u+dm[0]*(j+dm[1]*(k+dm[2]*2));
                pbx  = b+dm[0]*(j+dm[1]*k);
                pby  = b+dm[0]*(j+dm[1]*(k+dm[2]));
                pbz  = b+dm[0]*(j+dm[1]*(k+dm[2]*2));

                if (a)
                {
                    paxx = a+dm[0]*(j+dm[1]*k);
                    payy = a+dm[0]*(j+dm[1]*(k+dm[2]));
                    pazz = a+dm[0]*(j+dm[1]*(k+dm[2]*2));
                    paxy = a+dm[0]*(j+dm[1]*(k+dm[2]*3));
                    paxz = a+dm[0]*(j+dm[1]*(k+dm[2]*4));
                    payz = a+dm[0]*(j+dm[1]*(k+dm[2]*5));
                }

                jm1 = (bound(j-1,dm[1])-j)*dm[0];
                jp1 = (bound(j+1,dm[1])-j)*dm[0];
                fjm1 = bound_factor(j-1, dm[1], 1); // Factor if sliding
                fjp1 = bound_factor(j+1, dm[1], 1);
                gjm1 = bound_factor(j-1, dm[1], 0);
                gjp1 = bound_factor(j+1, dm[1], 0);

                istart = (jstart == (j%2));

                for(i=istart; i<dm[0]; i+=2)
                {
                    double sux, suy, suz;
                    float *px = pux+i, *py = puy+i, *pz = puz+i;

                    im1 = bound(i-1,dm[0])-i;
                    ip1 = bound(i+1,dm[0])-i;
                    fim1 = bound_factor(i-1, dm[0], 1); // Factor if sliding
                    fip1 = bound_factor(i+1, dm[0], 1);
                    gim1 = bound_factor(i-1, dm[0], 0);
                    gip1 = bound_factor(i+1, dm[0], 0);

                    sux = pbx[i]-(w001*(gkm1*px[km1] + gkp1*px[kp1]) + 
                                  w010*(gjm1*px[jm1] + gjp1*px[jp1]) + 
                                  w100*(fim1*px[im1] + fip1*px[ip1]) )/(s[0]*s[0]);
                    suy = pby[i]-(w001*(gkm1*py[km1] + gkp1*py[kp1]) + 
                                  w010*(fjm1*py[jm1] + fjp1*py[jp1]) + 
                                  w100*(gim1*py[im1] + gip1*py[ip1]) )/(s[1]*s[1]);
                    suz = pbz[i]-(w001*(fkm1*pz[km1] + fkp1*pz[kp1]) + 
                                  w010*(gjm1*pz[jm1] + gjp1*pz[jp1]) + 
                                  w100*(gim1*pz[im1] + gip1*pz[ip1]) )/(s[2]*s[2]);

                    if (a)
                    {
                        double axx, ayy, azz, axy, axz, ayz, idt;
                        /*
                           syms axx ayy azz axy axz ayz sux suy suz
                           A = [axx axy axz; axy ayy ayz; axz ayz azz];
                           su = [sux ; suy; suz]
                           simplify(inv(A)*su)
                        */
                        axx = paxx[i] + w000/(s[0]*s[0]);
                        ayy = payy[i] + w000/(s[1]*s[1]);
                        azz = pazz[i] + w000/(s[2]*s[2]);
                        axy = paxy[i];
                        axz = paxz[i];
                        ayz = payz[i];
                        idt = 1.0/(axx*ayy*azz -axx*ayz*ayz-ayy*axz*axz-azz*axy*axy +2*axy*axz*ayz);
                        *px = idt*(sux*(ayy*azz-ayz*ayz)+suy*(axz*ayz-axy*azz)+suz*(axy*ayz-axz*ayy));
                        *py = idt*(sux*(axz*ayz-axy*azz)+suy*(axx*azz-axz*axz)+suz*(axy*axz-axx*ayz));
                        *pz = idt*(sux*(axy*ayz-axz*ayy)+suy*(axy*axz-axx*ayz)+suz*(axx*ayy-axy*axy));
                    }
                    else
                    {
                        *px = (s[0]*s[0])*sux/w000;
                        *py = (s[1]*s[1])*suy/w000;
                        *pz = (s[2]*s[2])*suz/w000;
                    }
                }
            }
        }
#   ifdef VERBOSE
        if ((it%2)==1) printf(" %g", sumsq(dm, a, b, s, u));
#   endif
    }
#ifdef VERBOSE
    printf("\n");
#endif

}

void Atimesp_me(mwSize dm[], float A[], double param[], float p[], float Ap[])
{
    vel2mom_me(dm, p, param, Ap);
    Atimesp1(dm, A, p, Ap);
}

/* Returns m = L * v in the case where L is made of a mixture of 
 * membrane and bending energies.
 *
 * \param dm     Dimensions of the velocity lattice (voxels)
 * \param f      Velocities computed on the above lattice
 * \param s[0:3] Voxel sizes (voxels/mm)
 * \param s[3]   Parameter of the absolute displacement penalty (lam0)
 * \param s[4]   Parameter of the membrane energy (lam1) (penalizes 
 *               elements of the Jacobian matrix -> 1st order smoothness)
 * \param s[5]   Parameter of the bending energy (lam2) (penalizes 
 *               elements of the Hessian matrix -> 2nd order smoothness)
 * \param g      Allocated array in which to store the output momentum 
 */
void vel2mom_be(mwSize dm[], float f[], double s[], float g[])
{
    mwSignedIndex k;
    double w000,w100,w200,
           w010,w110,
           w020,
           w001,w101,
           w011,
           w002;
    double lam0 = s[3], lam1 = s[4], lam2 = s[5];
    double v0 = s[0]*s[0], v1 = s[1]*s[1], v2 = s[2]*s[2];

    w000 = lam2*(6*(v0*v0+v1*v1+v2*v2) +8*(v0*v1+v0*v2+v1*v2)) +lam1*2*(v0+v1+v2) + lam0;
    w100 = lam2*(-4*v0*(v0+v1+v2)) -lam1*v0;
    w010 = lam2*(-4*v1*(v0+v1+v2)) -lam1*v1;
    w001 = lam2*(-4*v2*(v0+v1+v2)) -lam1*v2;
    w200 = lam2*v0*v0;
    w020 = lam2*v1*v1;
    w002 = lam2*v2*v2;
    w110 = lam2*2*v0*v1;
    w101 = lam2*2*v0*v2;
    w011 = lam2*2*v1*v2;

    /*
        syms s1 s2 s3
        syms l1 l2 l3
        zz = sym(zeros(3,3));
        K1 = cat(3,zz,[0 -s1*s1 0; 0 2*s1*s1 0; 0 -s1*s1 0],zz);
        K2 = cat(3,zz,[0 0 0; -s2*s2 2*s2*s2 -s2*s2; 0 0 0],zz);
        K3 = sym(zeros(3,3,3));
        K3(2,2,1) = -s3*s3;
        K3(2,2,2) = 2*s3*s3;
        K3(2,2,3) = -s3*s3;

        K  = K1+K2+K3;
        K1 = K*l1; K1(2,2,2) = K1(2,2,2)+l2;
        K2 = K*l1; K2(2,2,2) = K2(2,2,2)+l3;

        % L  = convn(K,K)
        L  = sym(zeros(5,5,5));
        for i=1:3,
            for j=1:3,
                for k=1:3,
                    L(i-1+1:i+1+1,j-1+1:j+1+1,k-1+1:k+1+1) = L(i-1+1:i+1+1,j-1+1:j+1+1,k-1+1:k+1+1) + K1(i,j,k)*K2;
                end;
            end;
        end;
        disp(simplify(L(3:end,3:end,3:end)))
    */

    for(k=0; k<dm[2]; k++)
    {
        mwSignedIndex j, km2,km1,kp1,kp2;
        int fkm1, fkm2, fkp1, fkp2, gkm1, gkm2, gkp1, gkp2;
        km2 = (bound(k-2,dm[2])-k)*dm[0]*dm[1];
        km1 = (bound(k-1,dm[2])-k)*dm[0]*dm[1];
        kp1 = (bound(k+1,dm[2])-k)*dm[0]*dm[1];
        kp2 = (bound(k+2,dm[2])-k)*dm[0]*dm[1];
        fkm1 = bound_factor(k-1, dm[2], 1); // Factor if sliding
        fkm2 = bound_factor(k-2, dm[2], 1);
        fkp1 = bound_factor(k+1, dm[2], 1);
        fkp2 = bound_factor(k+2, dm[2], 1);
        gkm1 = bound_factor(k-1, dm[2], 0);
        gkm2 = bound_factor(k-2, dm[2], 0);
        gkp1 = bound_factor(k+1, dm[2], 0);
        gkp2 = bound_factor(k+2, dm[2], 0);

        for(j=0; j<dm[1]; j++)
        {
            mwSignedIndex i, jm2,jm1,jp1,jp2;
            int fjm1, fjm2, fjp1, fjp2, gjm1, gjm2, gjp1, gjp2;
            float *pgx, *pgy, *pgz, *pfx, *pfy, *pfz;

            pgx = g+dm[0]*(j+dm[1]*k);
            pgy = g+dm[0]*(j+dm[1]*(k+dm[2]));
            pgz = g+dm[0]*(j+dm[1]*(k+dm[2]*2));

            pfx = f+dm[0]*(j+dm[1]*k);
            pfy = f+dm[0]*(j+dm[1]*(k+dm[2]));
            pfz = f+dm[0]*(j+dm[1]*(k+dm[2]*2));

            jm2 = (bound(j-2,dm[1])-j)*dm[0];
            jm1 = (bound(j-1,dm[1])-j)*dm[0];
            jp1 = (bound(j+1,dm[1])-j)*dm[0];
            jp2 = (bound(j+2,dm[1])-j)*dm[0];
            fjm1 = bound_factor(j-1, dm[1], 1);
            fjm2 = bound_factor(j-2, dm[1], 1);
            fjp1 = bound_factor(j+1, dm[1], 1);
            fjp2 = bound_factor(j+2, dm[1], 1);
            gjm1 = bound_factor(j-1, dm[1], 0);
            gjm2 = bound_factor(j-2, dm[1], 0);
            gjp1 = bound_factor(j+1, dm[1], 0);
            gjp2 = bound_factor(j+2, dm[1], 0);

            for(i=0; i<dm[0]; i++)
            {
                mwSignedIndex im2,im1,ip1,ip2;
                int fim1, fim2, fip1, fip2, gim1, gim2, gip1, gip2;
                float *px = &pfx[i], *py = &pfy[i], *pz = &pfz[i];
                double c;

                im2 = bound(i-2,dm[0])-i;
                im1 = bound(i-1,dm[0])-i;
                ip1 = bound(i+1,dm[0])-i;
                ip2 = bound(i+2,dm[0])-i;
                fim1 = bound_factor(i-1, dm[0], 1); // Factor if sliding
                fim2 = bound_factor(i-2, dm[0], 1);
                fip1 = bound_factor(i+1, dm[0], 1);
                fip2 = bound_factor(i+2, dm[0], 1);
                gim1 = bound_factor(i-1, dm[0], 0);
                gim2 = bound_factor(i-2, dm[0], 0);
                gip1 = bound_factor(i+1, dm[0], 0);
                gip2 = bound_factor(i+2, dm[0], 0);

                /* Note that a few things have been done here to reduce rounding errors.
                   This may slow things down, but it does lead to more accuracy. */
                c      = px[0];
                pgx[i] =(lam0*c
                       + w100*((fim1*px[im1        ]-c) + (fip1*px[ip1        ]-c))
                       + w010*((gjm1*px[    jm1    ]-c) + (gjp2*px[    jp1    ]-c))
                       + w001*((gkm1*px[        km1]-c) + (gkp2*px[        kp1]-c))
                       + w200*((fim2*px[im2        ]-c) + (fip2*px[ip2        ]-c))
                       + w020*((gjm2*px[    jm2    ]-c) + (gjp2*px[    jp2    ]-c))
                       + w002*((gkm2*px[        km2]-c) + (gkp2*px[        kp2]-c))
                       + w110*((fim1*gjm1*px[im1+jm1    ]-c) + (fip1*gjm1*px[ip1+jm1    ]-c) + (fim1*gjp1*px[im1+jp1    ]-c) + (fip1*gjp1*px[ip1+jp1    ]-c))
                       + w101*((fim1*gkm1*px[im1    +km1]-c) + (fip1*gkm1*px[ip1    +km1]-c) + (fim1*gkp1*px[im1    +kp1]-c) + (fip1*gkp1*px[ip1    +kp1]-c))
                       + w011*((gjm1*gkm1*px[    jm1+km1]-c) + (gjp1*gkm1*px[    jp1+km1]-c) + (gjm1*gkp1*px[    jm1+kp1]-c) + (gjp1*gkp1*px[    jp1+kp1]-c)))/v0;

                c      = py[0];
                pgy[i] =(lam0*c
                       + w100*((gim1*py[im1        ]-c) + (gip1*py[ip1        ]-c))
                       + w010*((fjm1*py[    jm1    ]-c) + (fjp1*py[    jp1    ]-c))
                       + w001*((gkm1*py[        km1]-c) + (gkp1*py[        kp1]-c))
                       + w200*((gim1*py[im2        ]-c) + (gip2*py[ip2        ]-c))
                       + w020*((fjm2*py[    jm2    ]-c) + (fjp2*py[    jp2    ]-c))
                       + w002*((gkm2*py[        km2]-c) + (gkp2*py[        kp2]-c))
                       + w110*((gim1*fjm1*py[im1+jm1    ]-c) + (gip1*fjm1*py[ip1+jm1    ]-c) + (gim1*fjp1*py[im1+jp1    ]-c) + (gip1*fjp1*py[ip1+jp1    ]-c))
                       + w101*((gim1*gkm1*py[im1    +km1]-c) + (gip1*gkm1*py[ip1    +km1]-c) + (gim1*gkp1*py[im1    +kp1]-c) + (gip1*gkp1*py[ip1    +kp1]-c))
                       + w011*((fjm1*gkm1*py[    jm1+km1]-c) + (fjp1*gkm1*py[    jp1+km1]-c) + (fjm1*gkp1*py[    jm1+kp1]-c) + (fjp1*gkp1*py[    jp1+kp1]-c)))/v1;

                c      = pz[0];
                pgz[i] =(lam0*c
                       + w100*((gim1*pz[im1        ]-c) + (gip1*pz[ip1        ]-c))
                       + w010*((gjm1*pz[    jm1    ]-c) + (gjp1*pz[    jp1    ]-c))
                       + w001*((fkm1*pz[        km1]-c) + (fkp1*pz[        kp1]-c))
                       + w200*((gim2*pz[im2        ]-c) + (gip2*pz[ip2        ]-c))
                       + w020*((gjm2*pz[    jm2    ]-c) + (gjp2*pz[    jp2    ]-c))
                       + w002*((fkm2*pz[        km2]-c) + (fkp2*pz[        kp2]-c))
                       + w200*((gim2*pz[im2        ]-c) + (gip2*pz[ip2        ]-c))
                       + w110*((gim1*gjm1*pz[im1+jm1    ]-c) + (gip1*gjm1*pz[ip1+jm1    ]-c) + (gim1*gjp1*pz[im1+jp1    ]-c) + (gip1*gjp1*pz[ip1+jp1    ]-c))
                       + w101*((gim1*fkm1*pz[im1    +km1]-c) + (gip1*fkm1*pz[ip1    +km1]-c) + (gim1*fkp1*pz[im1    +kp1]-c) + (gip1*fkp1*pz[ip1    +kp1]-c))
                       + w011*((gjm1*fkm1*pz[    jm1+km1]-c) + (gjp1*fkm1*pz[    jp1+km1]-c) + (gjm1*fkp1*pz[    jm1+kp1]-c) + (fkp1*pz[    jp1+kp1]-c)))/v2;
            }
        }
    }
}

void relax_be(mwSize dm[], float a[], float b[], double s[], int nit, float u[])
{
    int it;
    double w000,w100,w200,
           w010,w110,
           w020,
           w001,w101,
           w011,
           w002;
    double lam0 = s[3], lam1 = s[4], lam2 = s[5];
    double v0 = s[0]*s[0], v1 = s[1]*s[1], v2 = s[2]*s[2];

    w000 = lam2*(6*(v0*v0+v1*v1+v2*v2) +8*(v0*v1+v0*v2+v1*v2)) +lam1*2*(v0+v1+v2) + lam0;
    w100 = lam2*(-4*v0*(v0+v1+v2)) -lam1*v0;
    w010 = lam2*(-4*v1*(v0+v1+v2)) -lam1*v1;
    w001 = lam2*(-4*v2*(v0+v1+v2)) -lam1*v2;
    w200 = lam2*v0*v0;
    w020 = lam2*v1*v1;
    w002 = lam2*v2*v2;
    w110 = lam2*2*v0*v1;
    w101 = lam2*2*v0*v2;
    w011 = lam2*2*v1*v2;

    w000 = w000*1.00001 + 1e-6;

    if (dm[0]<=2)
    {
        w000 += 2*w200;
        w200  = 0.0;
    }
    if (dm[1]<=2)
    {
        w000 += 2*w020;
        w020  = 0.0;
    }
    if (dm[2]<=2)
    {
        w000 += 2*w002;
        w002  = 0.0;
    }

    if (dm[0]==1)
    {
        w000 += 2*w100;
        w100  = 0.0;
        if (dm[1]==1)
        {
            w000 += 4*w110;
            w110  = 0.0;
        }
        if (dm[2]==1)
        {
            w000 += 4*w101;
            w101  = 0.0;
        }
    }
    if (dm[1]==1)
    {
        w000 += 2*w010;
        w010  = 0.0;
        if (dm[2]==1)
        {
            w000 += 4*w011;
            w011  = 0.0;
        }
    }
    if (dm[2]==1)
    {
        w000 += 2*w001;
        w001  = 0.0;
    }

#   ifdef VERBOSE
        for(it=0; it< 10-(int)ceil(1.44269504088896*log((double)dm[0])); it++) printf("  ");
        printf("B%dx%dx%d: %g ", dm[0],dm[1],dm[2],sumsq(dm, a, b, s, u));
#   endif

    for(it=0; it<27*nit; it++)
    {
        mwSignedIndex i, j, k;
        for(k=(it/9)%3; k<dm[2]; k+=3)
        {
            mwSignedIndex km2, km1, kp1, kp2;
            int fkm1, fkm2, fkp1, fkp2, gkm1, gkm2, gkp1, gkp2;
            km2 = (bound(k-2,dm[2])-k)*dm[0]*dm[1];
            km1 = (bound(k-1,dm[2])-k)*dm[0]*dm[1];
            kp1 = (bound(k+1,dm[2])-k)*dm[0]*dm[1];
            kp2 = (bound(k+2,dm[2])-k)*dm[0]*dm[1];
            fkm1 = bound_factor(k-1, dm[2], 1);
            fkm2 = bound_factor(k-2, dm[2], 1);
            fkp1 = bound_factor(k+1, dm[2], 1);
            fkp2 = bound_factor(k+2, dm[2], 1);
            gkm1 = bound_factor(k-1, dm[2], 0);
            gkm2 = bound_factor(k-2, dm[2], 0);
            gkp1 = bound_factor(k+1, dm[2], 0);
            gkp2 = bound_factor(k+2, dm[2], 0);

            for(j=(it/3)%3; j<dm[1]; j+=3)
            {
                float *pux, *puy, *puz, *pbx, *pby, *pbz, *paxx, *payy, *pazz, *paxy, *paxz, *payz;
                mwSignedIndex jm2,jm1,jp1,jp2;
                int fjm1, fjm2, fjp1, fjp2, gjm1, gjm2, gjp1, gjp2;

                pux  = u+dm[0]*(j+dm[1]* k);
                puy  = u+dm[0]*(j+dm[1]*(k+dm[2]));
                puz  = u+dm[0]*(j+dm[1]*(k+dm[2]*2));
                pbx  = b+dm[0]*(j+dm[1]* k);
                pby  = b+dm[0]*(j+dm[1]*(k+dm[2]));
                pbz  = b+dm[0]*(j+dm[1]*(k+dm[2]*2));

                if (a)
                {
                    paxx = a+dm[0]*(j+dm[1]* k);
                    payy = a+dm[0]*(j+dm[1]*(k+dm[2]));
                    pazz = a+dm[0]*(j+dm[1]*(k+dm[2]*2));
                    paxy = a+dm[0]*(j+dm[1]*(k+dm[2]*3));
                    paxz = a+dm[0]*(j+dm[1]*(k+dm[2]*4));
                    payz = a+dm[0]*(j+dm[1]*(k+dm[2]*5));
                }

                jm2 = (bound(j-2,dm[1])-j)*dm[0];
                jm1 = (bound(j-1,dm[1])-j)*dm[0];
                jp1 = (bound(j+1,dm[1])-j)*dm[0];
                jp2 = (bound(j+2,dm[1])-j)*dm[0];
                fjm1 = bound_factor(j-1, dm[1], 1); // Factor if sliding
                fjm2 = bound_factor(j-2, dm[1], 1);
                fjp1 = bound_factor(j+1, dm[1], 1);
                fjp2 = bound_factor(j+2, dm[1], 1);
                gjm1 = bound_factor(j-1, dm[1], 0);
                gjm2 = bound_factor(j-2, dm[1], 0);
                gjp1 = bound_factor(j+1, dm[1], 0);
                gjp2 = bound_factor(j+2, dm[1], 0);

                for(i=it%3; i<dm[0]; i+=3)
                {
                    mwSignedIndex im2,im1,ip1,ip2;
                    int fim1, fim2, fip1, fip2, gim1, gim2, gip1, gip2;
                    double sux, suy, suz, c;
                    float *px = pux+i, *py = puy+i, *pz = puz+i;

                    im2 = bound(i-2,dm[0])-i;
                    im1 = bound(i-1,dm[0])-i;
                    ip1 = bound(i+1,dm[0])-i;
                    ip2 = bound(i+2,dm[0])-i;
                    fim1 = bound_factor(i-1, dm[0], 1); // Factor if sliding
                    fim2 = bound_factor(i-2, dm[0], 1);
                    fip1 = bound_factor(i+1, dm[0], 1);
                    fip2 = bound_factor(i+2, dm[0], 1);
                    gim1 = bound_factor(i-1, dm[0], 0);
                    gim2 = bound_factor(i-2, dm[0], 0);
                    gip1 = bound_factor(i+1, dm[0], 0);
                    gip2 = bound_factor(i+2, dm[0], 0);

                    /* Note that a few things have been done here to reduce rounding errors.
                       This may slow things down, but it does lead to more accuracy. */
                    c   = px[0];
                    sux = pbx[i] - (lam0*c
                                  + w100*((fim1*px[im1        ]-c) + (fip1*px[ip1        ]-c))
                                  + w010*((gjm1*px[    jm1    ]-c) + (gjp1*px[    jp1    ]-c))
                                  + w001*((gkm1*px[        km1]-c) + (gkp1*px[        kp1]-c))
                                  + w200*((fim2*px[im2        ]-c) + (fip2*px[ip2        ]-c))
                                  + w020*((gjm2*px[    jm2    ]-c) + (gjp2*px[    jp2    ]-c))
                                  + w002*((gkm2*px[        km2]-c) + (gkp2*px[        kp2]-c))
                                  + w110*((fim1*gjm1*px[im1+jm1    ]-c) + (fip1*gjm1*px[ip1+jm1    ]-c) + (fim1*gjp1*px[im1+jp1    ]-c) + (fip1*gjp1*px[ip1+jp1    ]-c))
                                  + w101*((fim1*gkm1*px[im1    +km1]-c) + (fip1*gkm1*px[ip1    +km1]-c) + (fim1*gkp1*px[im1    +kp1]-c) + (fip1*gkp1*px[ip1    +kp1]-c))
                                  + w011*((gjm1*gkm1*px[    jm1+km1]-c) + (gjp1*gkm1*px[    jp1+km1]-c) + (gjm1*gkp1*px[    jm1+kp1]-c) + (gjp1*gkp1*px[    jp1+kp1]-c)))/v0;

                    c   = py[0];
                    suy = pby[i] - (lam0*c
                                  + w100*((gim1*py[im1        ]-c) + (gip1*py[ip1        ]-c))
                                  + w010*((fjm1*py[    jm1    ]-c) + (fjp1*py[    jp1    ]-c))
                                  + w001*((gkm1*py[        km1]-c) + (gkp1*py[        kp1]-c))
                                  + w200*((gim2*py[im2        ]-c) + (gip2*py[ip2        ]-c))
                                  + w020*((fjm2*py[    jm2    ]-c) + (fjp2*py[    jp2    ]-c))
                                  + w002*((gkm2*py[        km2]-c) + (gkp2*py[        kp2]-c))
                                  + w110*((gim1*fjm1*py[im1+jm1    ]-c) + (gip1*fjm1*py[ip1+jm1    ]-c) + (gim1*fjp1*py[im1+jp1    ]-c) + (gip1*fjp1*py[ip1+jp1    ]-c))
                                  + w101*((gim1*gkm1*py[im1    +km1]-c) + (gip1*gkm1*py[ip1    +km1]-c) + (gim1*gkp1*py[im1    +kp1]-c) + (gip1*gkp1*py[ip1    +kp1]-c))
                                  + w011*((fjm1*gkm1*py[    jm1+km1]-c) + (fjp1*gkm1*py[    jp1+km1]-c) + (fjm1*gkp1*py[    jm1+kp1]-c) + (fjp1*gkp1*py[    jp1+kp1]-c)))/v1;

                    c   = pz[0];
                    suz = pbz[i] - (lam0*c
                                  + w100*((gim1*pz[im1        ]-c) + (gip1*pz[ip1        ]-c))
                                  + w010*((gjm1*pz[    jm1    ]-c) + (gjp1*pz[    jp1    ]-c))
                                  + w001*((fkm1*pz[        km1]-c) + (fkp1*pz[        kp1]-c))
                                  + w200*((gim2*pz[im2        ]-c) + (gip2*pz[ip2        ]-c))
                                  + w020*((gjm2*pz[    jm2    ]-c) + (gjp2*pz[    jp2    ]-c))
                                  + w002*((fkm2*pz[        km2]-c) + (fkp2*pz[        kp2]-c))
                                  + w110*((gim1*gjm1*pz[im1+jm1    ]-c) + (gip1*gjm1*pz[ip1+jm1    ]-c) + (gim1*gjp1*pz[im1+jp1    ]-c) + (gip1*gjp1*pz[ip1+jp1    ]-c))
                                  + w101*((gim1*fkm1*pz[im1    +km1]-c) + (gip1*fkm1*pz[ip1    +km1]-c) + (gim1*fkp1*pz[im1    +kp1]-c) + (gip1*fkp1*pz[ip1    +kp1]-c))
                                  + w011*((gjm1*fkm1*pz[    jm1+km1]-c) + (gjp1*fkm1*pz[    jp1+km1]-c) + (gjm1*fkp1*pz[    jm1+kp1]-c) + (gjp1*fkp1*pz[    jp1+kp1]-c)))/v2;

                    if (a)
                    {
                        double axx, ayy, azz, axy, axz, ayz, idt;

                        sux -= (paxx[i]*px[0] + paxy[i]*py[0] + paxz[i]*pz[0]);
                        suy -= (paxy[i]*px[0] + payy[i]*py[0] + payz[i]*pz[0]);
                        suz -= (paxz[i]*px[0] + payz[i]*py[0] + pazz[i]*pz[0]);

                        axx  = paxx[i] + w000/v0;
                        ayy  = payy[i] + w000/v1;
                        azz  = pazz[i] + w000/v2;
                        axy  = paxy[i];
                        axz  = paxz[i];
                        ayz  = payz[i];
                        idt  = 1.0/(axx*ayy*azz -axx*ayz*ayz-ayy*axz*axz-azz*axy*axy +2*axy*axz*ayz);
                        *px += idt*(sux*(ayy*azz-ayz*ayz)+suy*(axz*ayz-axy*azz)+suz*(axy*ayz-axz*ayy));
                        *py += idt*(sux*(axz*ayz-axy*azz)+suy*(axx*azz-axz*axz)+suz*(axy*axz-axx*ayz));
                        *pz += idt*(sux*(axy*ayz-axz*ayy)+suy*(axy*axz-axx*ayz)+suz*(axx*ayy-axy*axy));
                    }
                    else
                    {
                        *px += v0*sux/w000;
                        *py += v1*suy/w000;
                        *pz += v2*suz/w000;
                    }
                }
            }
        }
#       ifdef VERBOSE
        if ((it%27) == 26)
            printf(" %g", sumsq(dm, a, b, s, u));
#       endif
    }
#   ifdef VERBOSE
        printf("\n");
#   endif
}


void Atimesp_be(mwSize dm[], float A[], double param[], float p[], float Ap[])
{
    vel2mom_be(dm, p, param, Ap);
    Atimesp1(dm, A, p, Ap);
}






/************************************************************************************************/

/* Returns m = L * v in the case where L is made of a mixture of 
 * regularization energies.
 *
 * \param dm     Dimensions of the velocity lattice (voxels)
 * \param f      Velocities computed on the above lattice
 * \param s[0:3] Voxel sizes (voxels/mm)
 * \param s[3]   Parameter of the absolute displacement penalty (lam0)
 * \param s[4]   Parameter of the membrane energy (lam1) (penalizes 
 *               elements of the Jacobian matrix -> 1st order smoothness)
 * \param s[5]   Parameter of the bending energy (lam2) (penalizes 
 *               elements of the Hessian matrix -> 2nd order smoothness)
 * \param s[6]   Parameter of the linear elastic energy (mu) (penalizes 
 *               elements of the symmetric part of the Jacobian matrix 
 *               -> penalizes scaling and shearing)
 * \param s[7]   Parameter of the linear elastic energy (lam) (penalizes 
 *               the divergence of the Jacobian matrix -> preserves 
 *               volumes)
 * \param g      Allocated array in which to store the output momentum 
 */
void vel2mom_all(mwSize dm[], float f[], double s[], float g[])
{
    mwSignedIndex k;
    double w000,w100,w200,
           w010,w110,
           w020,
           w001,w101,
           w011,
           w002;
    double v0 = s[0]*s[0], v1 = s[1]*s[1], v2 = s[2]*s[2];
    double lam0 = s[3], lam1 = s[4], lam2 = s[5], mu = s[6], lam = s[7];
    double wx000, wx100, wx010, wx001, wy000, wy100, wy010, wy001, wz000, wz100, wz010, wz001, w2;

    // Convolution kernel W for the absolute displacement, membrane and 
    // bending energies
    w000 = lam2*(6*(v0*v0+v1*v1+v2*v2) +8*(v0*v1+v0*v2+v1*v2)) +lam1*2*(v0+v1+v2) + lam0;
    w100 = lam2*(-4*v0*(v0+v1+v2)) -lam1*v0;
    w010 = lam2*(-4*v1*(v0+v1+v2)) -lam1*v1;
    w001 = lam2*(-4*v2*(v0+v1+v2)) -lam1*v2;
    w200 = lam2*v0*v0;
    w020 = lam2*v1*v1;
    w002 = lam2*v2*v2;
    w110 = lam2*2*v0*v1;
    w101 = lam2*2*v0*v2;
    w011 = lam2*2*v1*v2;

    // Convolution kernels W_x, W_y, W_z (in the x, y and z directions) 
    // for the linear elastic energy
    wx000 =  2*mu*(2*v0+v1+v2)/v0+2*lam + w000/v0;
    wx100 = -2*mu-lam + w100/v0;
    wx010 = -mu*v1/v0 + w010/v0;
    wx001 = -mu*v2/v0 + w001/v0;
    wy000 =  2*mu*(v0+2*v1+v2)/v1+2*lam + w000/v1;
    wy100 = -mu*v0/v1 + w100/v1;
    wy010 = -2*mu-lam + w010/v1;
    wy001 = -mu*v2/v1 + w001/v1;
    wz000 =  2*mu*(v0+v1+2*v2)/v2+2*lam + w000/v2;
    wz100 = -mu*v0/v2 + w100/v2;
    wz010 = -mu*v1/v2 + w010/v2;
    wz001 = -2*mu-lam + w001/v2;
    w2    = 0.25*mu+0.25*lam;

    // Now, let's perform the convolution. (*) will denote the convolution 
    // operator:
    // m_x = W (*) v_x + W_x (*) v_x
    // m_y = W (*) v_y + W_y (*) v_y
    // m_z = W (*) v_z + W_z (*) v_z
    // (Plus millimetric corrections induced by multiplications 
    // with voxel sizes)
    for(k=0; k<dm[2]; k++)
    {
        mwSignedIndex j, km2,km1,kp1,kp2;
        int fkm1, fkm2, fkp1, fkp2, gkm1, gkm2, gkp1, gkp2;
        km2 = (bound(k-2,dm[2])-k)*dm[0]*dm[1];
        km1 = (bound(k-1,dm[2])-k)*dm[0]*dm[1];
        kp1 = (bound(k+1,dm[2])-k)*dm[0]*dm[1];
        kp2 = (bound(k+2,dm[2])-k)*dm[0]*dm[1];
        fkm1 = bound_factor(k-1, dm[2], 1); // Factor if sliding
        fkm2 = bound_factor(k-2, dm[2], 1);
        fkp1 = bound_factor(k+1, dm[2], 1);
        fkp2 = bound_factor(k+2, dm[2], 1);
        gkm1 = bound_factor(k-1, dm[2], 0);
        gkm2 = bound_factor(k-2, dm[2], 0);
        gkp1 = bound_factor(k+1, dm[2], 0);
        gkp2 = bound_factor(k+2, dm[2], 0);

        for(j=0; j<dm[1]; j++)
        {
            mwSignedIndex i, jm2,jm1,jp1,jp2;
            int fjm1, fjm2, fjp1, fjp2, gjm1, gjm2, gjp1, gjp2;
            float *pgx, *pgy, *pgz, *pfx, *pfy, *pfz;

            pgx = g+dm[0]*(j+dm[1]*k);
            pgy = g+dm[0]*(j+dm[1]*(k+dm[2]));
            pgz = g+dm[0]*(j+dm[1]*(k+dm[2]*2));

            pfx = f+dm[0]*(j+dm[1]*k);
            pfy = f+dm[0]*(j+dm[1]*(k+dm[2]));
            pfz = f+dm[0]*(j+dm[1]*(k+dm[2]*2));

            jm2 = (bound(j-2,dm[1])-j)*dm[0];
            jm1 = (bound(j-1,dm[1])-j)*dm[0];
            jp1 = (bound(j+1,dm[1])-j)*dm[0];
            jp2 = (bound(j+2,dm[1])-j)*dm[0];
            fjm1 = bound_factor(j-1, dm[1], 1); // Factor if sliding
            fjm2 = bound_factor(j-2, dm[1], 1);
            fjp1 = bound_factor(j+1, dm[1], 1);
            fjp2 = bound_factor(j+2, dm[1], 1);
            gjm1 = bound_factor(j-1, dm[1], 0);
            gjm2 = bound_factor(j-2, dm[1], 0);
            gjp1 = bound_factor(j+1, dm[1], 0);
            gjp2 = bound_factor(j+2, dm[1], 0);

            for(i=0; i<dm[0]; i++)
            {
                mwSignedIndex im2,im1,ip1,ip2;
                int fim1, fim2, fip1, fip2, gim1, gim2, gip1, gip2;
                float *px = &pfx[i], *py = &pfy[i], *pz = &pfz[i];
                double c;

                im2 = bound(i-2,dm[0])-i;
                im1 = bound(i-1,dm[0])-i;
                ip1 = bound(i+1,dm[0])-i;
                ip2 = bound(i+2,dm[0])-i;
                fim1 = bound_factor(i-1, dm[0], 1); // Factor if sliding
                fim2 = bound_factor(i-2, dm[0], 1);
                fip1 = bound_factor(i+1, dm[0], 1);
                fip2 = bound_factor(i+2, dm[0], 1);
                gim1 = bound_factor(i-1, dm[0], 0);
                gim2 = bound_factor(i-2, dm[0], 0);
                gip1 = bound_factor(i+1, dm[0], 0);
                gip2 = bound_factor(i+2, dm[0], 0);

                /* Note that a few things have been done here to reduce rounding errors.
                   This may slow things down, but it does lead to more accuracy. */
                c      = px[0];
                pgx[i] = (wx100*((fim1*px[im1        ]-c) + (fip1*px[ip1        ]-c))
                        + wx010*((gjm1*px[    jm1    ]-c) + (gjp1*px[    jp1    ]-c))
                        + wx001*((gkm1*px[        km1]-c) + (gkp1*px[        kp1]-c))
                        + w2   *( gip1*fjm1*py[ip1+jm1] - gip1*fjp1*py[ip1+jp1] + 
                                  gim1*fjp1*py[im1+jp1] - gim1*fjm1*py[im1+jm1] + 
                                  gip1*fkm1*pz[ip1+km1] - gip1*fkp1*pz[ip1+kp1] + 
                                  gim1*fkp1*pz[im1+kp1] - gim1*fkm1*pz[im1+km1])
                        + (lam0*c
                        +  w110*((fim1*gjm1*px[im1+jm1    ]-c) + (fip1*gjm1*px[ip1+jm1    ]-c) + (fim1*gjp1*px[im1+jp1    ]-c) + (fip1*gjp1*px[ip1+jp1    ]-c))
                        +  w101*((fim1*gkm1*px[im1    +km1]-c) + (fip1*gkm1*px[ip1    +km1]-c) + (fim1*gkp1*px[im1    +kp1]-c) + (fip1*gkp1*px[ip1    +kp1]-c))
                        +  w011*((gjm1*gkm1*px[    jm1+km1]-c) + (gjp1*gkm1*px[    jp1+km1]-c) + (gjm1*gkp1*px[    jm1+kp1]-c) + (gjp1*gkp1*px[    jp1+kp1]-c))
                        +  w200*((fim2*px[im2        ]-c) + (fip2*px[ip2        ]-c))
                        +  w020*((gjm2*px[    jm2    ]-c) + (gjp2*px[    jp2    ]-c))
                        +  w002*((gkm2*px[        km2]-c) + (gkp2*px[        kp2]-c)))/v0);

                c      = py[0];
                pgy[i] = (wy100*((gim1*py[im1        ]-c) + (gip1*py[ip1        ]-c))
                        + wy010*((fjm1*py[    jm1    ]-c) + (fjp1*py[    jp1    ]-c))
                        + wy001*((gkm1*py[        km1]-c) + (gkp1*py[        kp1]-c))
                        + w2   *( gjp1*fim1*px[jp1+im1] - gjp1*fip1*px[jp1+ip1] + 
                                  gjm1*fip1*px[jm1+ip1] - gjm1*fim1*px[jm1+im1] + 
                                  gjp1*fkm1*pz[jp1+km1] - gjp1*fkp1*pz[jp1+kp1] + 
                                  gjm1*fkp1*pz[jm1+kp1] - gjm1*fkm1*pz[jm1+km1])
                        + (lam0*c
                        +  w110*((gim1*fjm1*py[im1+jm1    ]-c) + (gip1*fjm1*py[ip1+jm1    ]-c) + (gim1*fjp1*py[im1+jp1    ]-c) + (gip1*fjp1*py[ip1+jp1    ]-c))
                        +  w101*((gim1*gkm1*py[im1    +km1]-c) + (gip1*gkm1*py[ip1    +km1]-c) + (gim1*gkp1*py[im1    +kp1]-c) + (gip1*gkp1*py[ip1    +kp1]-c))
                        +  w011*((fjm1*gkm1*py[    jm1+km1]-c) + (fjp1*gkm1*py[    jp1+km1]-c) + (fjm1*gkp1*py[    jm1+kp1]-c) + (fjp1*gkp1*py[    jp1+kp1]-c))
                        +  w200*((gim2*py[im2        ]-c) + (gip2*py[ip2        ]-c))
                        +  w020*((fjm2*py[    jm2    ]-c) + (fjp2*py[    jp2    ]-c))
                        +  w002*((gkm2*py[        km2]-c) + (gkp2*py[        kp2]-c)))/v1);

                c      = pz[0];
                pgz[i] = (wz100*((gim1*pz[im1        ]-c) + (gip1*pz[ip1        ]-c))
                        + wz010*((gjm1*pz[    jm1    ]-c) + (gjp1*pz[    jp1    ]-c))
                        + wz001*((fkm1*pz[        km1]-c) + (fkp1*pz[        kp1]-c))
                        + w2   *( gkp1*fim1*px[kp1+im1] - gkp1*fip1*px[kp1+ip1] + 
                                  gkm1*fip1*px[km1+ip1] - gkm1*fim1*px[km1+im1] + 
                                  gkp1*fjm1*py[kp1+jm1] - gkp1*fjp1*py[kp1+jp1] + 
                                  gkm1*fjp1*py[km1+jp1] - gkm1*fjm1*py[km1+jm1])
                        + (lam0*c
                        +  w110*((gim1*gjm1*pz[im1+jm1    ]-c) + (gip1*gjm1*pz[ip1+jm1    ]-c) + (gim1*gjp1*pz[im1+jp1    ]-c) + (gip1*gjp1*pz[ip1+jp1    ]-c))
                        +  w101*((gim1*fkm1*pz[im1    +km1]-c) + (gip1*fkm1*pz[ip1    +km1]-c) + (gim1*fkp1*pz[im1    +kp1]-c) + (gip1*fkp1*pz[ip1    +kp1]-c))
                        +  w011*((gjm1*fkm1*pz[    jm1+km1]-c) + (gjp1*fkm1*pz[    jp1+km1]-c) + (gjm1*fkp1*pz[    jm1+kp1]-c) + (gjp1*fkp1*pz[    jp1+kp1]-c))
                        +  w200*((gim2*pz[im2        ]-c) + (gip2*pz[ip2        ]-c))
                        +  w020*((gjm2*pz[    jm2    ]-c) + (gjp2*pz[    jp2    ]-c))
                        +  w002*((fkm2*pz[        km2]-c) + (fkp2*pz[        kp2]-c)))/v2);
            }
        }
    }
}

void relax_all(mwSize dm[], float a[], float b[], double s[], int nit, float u[])
{
    int it;
    double w000,w100,w200,
           w010,w110,
           w020,
           w001,w101,
           w011,
           w002;
    double v0 = s[0]*s[0], v1 = s[1]*s[1], v2 = s[2]*s[2];
    double lam0 = s[3], lam1 = s[4], lam2 = s[5], mu = s[6], lam = s[7];
    double wx000, wx100, wx010, wx001, wy000, wy100, wy010, wy001, wz000, wz100, wz010, wz001, w2;

    w000 = lam2*(6*(v0*v0+v1*v1+v2*v2) +8*(v0*v1+v0*v2+v1*v2)) +lam1*2*(v0+v1+v2) + lam0;
    w100 = lam2*(-4*v0*(v0+v1+v2)) -lam1*v0;
    w010 = lam2*(-4*v1*(v0+v1+v2)) -lam1*v1;
    w001 = lam2*(-4*v2*(v0+v1+v2)) -lam1*v2;
    w200 = lam2*v0*v0;
    w020 = lam2*v1*v1;
    w002 = lam2*v2*v2;
    w110 = lam2*2*v0*v1;
    w101 = lam2*2*v0*v2;
    w011 = lam2*2*v1*v2;

    wx000 =  2*mu*(2*v0+v1+v2)/v0+2*lam + w000/v0;
    wx100 = -2*mu-lam + w100/v0;
    wx010 = -mu*v1/v0 + w010/v0;
    wx001 = -mu*v2/v0 + w001/v0;
    wy000 =  2*mu*(v0+2*v1+v2)/v1+2*lam + w000/v1;
    wy100 = -mu*v0/v1 + w100/v1;
    wy010 = -2*mu-lam + w010/v1;
    wy001 = -mu*v2/v1 + w001/v1;
    wz000 =  2*mu*(v0+v1+2*v2)/v2+2*lam + w000/v2;
    wz100 = -mu*v0/v2 + w100/v2;
    wz010 = -mu*v1/v2 + w010/v2;
    wz001 = -2*mu-lam + w001/v2;
    w2    = 0.25*mu+0.25*lam;

    wx000 = wx000*1.00001 + 1e-6;
    wy000 = wy000*1.00001 + 1e-6;
    wz000 = wz000*1.00001 + 1e-6;

    if (dm[0]<=2)
    {
        wx000 += 2*w200/v0;
        wy000 += 2*w200/v1;
        wz000 += 2*w200/v2;
        w200   = 0.0;
    }
    if (dm[1]<=2)
    {
        wx000 += 2*w020/v0;
        wy000 += 2*w020/v1;
        wz000 += 2*w020/v2;
        w020   = 0.0;
    }
    if (dm[2]<=2)
    {
        wx000 += 2*w002/v0;
        wy000 += 2*w002/v1;
        wz000 += 2*w002/v2; 
        w002   = 0.0;
    }

    if (dm[0]==1)
    {
        wx000 += 2*wx100; wx100  = 0.0;
        wy000 += 2*wy100; wy100  = 0.0;
        wz000 += 2*wz100; wz100  = 0.0;
        if (dm[1]==1)
        {
            wx000 += 4*w110/v0;
            wy000 += 4*w110/v1;
            wz000 += 4*w110/v2;
            w110   = 0.0;
        }
        if (dm[2]==1)
        {
            wx000 += 4*w101/v0;
            wy000 += 4*w101/v1;
            wz000 += 4*w101/v2;
            w101   = 0.0;
        }

    }
    if (dm[1]==1)
    {
        wx000 += 2*wx010; wx010  = 0.0;
        wy000 += 2*wy010; wy010  = 0.0;
        wz000 += 2*wz010; wz010  = 0.0;
        if (dm[2]==1)
        {
            wx000 += 4*w011/v0;
            wy000 += 4*w011/v1;
            wz000 += 4*w011/v2;
            w011   = 0.0;
        }
    }
    if (dm[2]==1)
    {
        wx000 += 2*wx001; wx001  = 0.0;
        wy000 += 2*wy001; wy001  = 0.0;
        wz000 += 2*wz001; wz001  = 0.0;
    }
    if (wx000<0.0) wx000;
    if (wy000<0.0) wy000;
    if (wz000<0.0) wz000;

#   ifdef VERBOSE
        for(it=0; it< 10-(int)ceil(1.44269504088896*log((double)dm[0])); it++) printf("  ");
        printf("A%dx%dx%d: %g ", dm[0],dm[1],dm[2],sumsq(dm, a, b, s, u));
#   endif

    for(it=0; it<27*nit; it++)
    {
        mwSignedIndex i, j, k;
        for(k=(it/9)%3; k<dm[2]; k+=3)
        {
            mwSignedIndex km2, km1, kp1, kp2;
            int fkm1, fkm2, fkp1, fkp2, gkm1, gkm2, gkp1, gkp2;
            km2 = (bound(k-2,dm[2])-k)*dm[0]*dm[1];
            km1 = (bound(k-1,dm[2])-k)*dm[0]*dm[1];
            kp1 = (bound(k+1,dm[2])-k)*dm[0]*dm[1];
            kp2 = (bound(k+2,dm[2])-k)*dm[0]*dm[1];
            fkm1 = bound_factor(k-1, dm[2], 1); // Factor if sliding
            fkm2 = bound_factor(k-2, dm[2], 1);
            fkp1 = bound_factor(k+1, dm[2], 1);
            fkp2 = bound_factor(k+2, dm[2], 1);
            gkm1 = bound_factor(k-1, dm[2], 0);
            gkm2 = bound_factor(k-2, dm[2], 0);
            gkp1 = bound_factor(k+1, dm[2], 0);
            gkp2 = bound_factor(k+2, dm[2], 0);

            for(j=(it/3)%3; j<dm[1]; j+=3)
            {
                float *pux, *puy, *puz, *pbx, *pby, *pbz, *paxx, *payy, *pazz, *paxy, *paxz, *payz;
                mwSignedIndex jm2,jm1,jp1,jp2;
                int fjm1, fjm2, fjp1, fjp2, gjm1, gjm2, gjp1, gjp2;

                pux  = u+dm[0]*(j+dm[1]* k);
                puy  = u+dm[0]*(j+dm[1]*(k+dm[2]));
                puz  = u+dm[0]*(j+dm[1]*(k+dm[2]*2));
                pbx  = b+dm[0]*(j+dm[1]* k);
                pby  = b+dm[0]*(j+dm[1]*(k+dm[2]));
                pbz  = b+dm[0]*(j+dm[1]*(k+dm[2]*2));

                if (a)
                {
                    paxx = a+dm[0]*(j+dm[1]* k);
                    payy = a+dm[0]*(j+dm[1]*(k+dm[2]));
                    pazz = a+dm[0]*(j+dm[1]*(k+dm[2]*2));
                    paxy = a+dm[0]*(j+dm[1]*(k+dm[2]*3));
                    paxz = a+dm[0]*(j+dm[1]*(k+dm[2]*4));
                    payz = a+dm[0]*(j+dm[1]*(k+dm[2]*5));
                }

                jm2 = (bound(j-2,dm[1])-j)*dm[0];
                jm1 = (bound(j-1,dm[1])-j)*dm[0];
                jp1 = (bound(j+1,dm[1])-j)*dm[0];
                jp2 = (bound(j+2,dm[1])-j)*dm[0];
                fjm1 = bound_factor(j-1, dm[1], 1); // Factor if sliding
                fjm2 = bound_factor(j-2, dm[1], 1);
                fjp1 = bound_factor(j+1, dm[1], 1);
                fjp2 = bound_factor(j+2, dm[1], 1);
                gjm1 = bound_factor(j-1, dm[1], 0);
                gjm2 = bound_factor(j-2, dm[1], 0);
                gjp1 = bound_factor(j+1, dm[1], 0);
                gjp2 = bound_factor(j+2, dm[1], 0);

                for(i=it%3; i<dm[0]; i+=3)
                {
                    mwSignedIndex im2,im1,ip1,ip2;
                    int fim1, fim2, fip1, fip2, gim1, gim2, gip1, gip2;
                    double sux, suy, suz, c;
                    float *px = pux+i, *py = puy+i, *pz = puz+i;

                    im2 = bound(i-2,dm[0])-i;
                    im1 = bound(i-1,dm[0])-i;
                    ip1 = bound(i+1,dm[0])-i;
                    ip2 = bound(i+2,dm[0])-i;
                    fim1 = bound_factor(i-1, dm[0], 1);
                    fim2 = bound_factor(i-2, dm[0], 1);
                    fip1 = bound_factor(i+1, dm[0], 1);
                    fip2 = bound_factor(i+2, dm[0], 1);

                    /* Note that a few things have been done here to reduce rounding errors.
                       This may slow things down, but it does lead to more accuracy. */
                    c   = px[0];
                    sux = pbx[i]
                          - ( wx100*((fim1*px[im1        ]-c) + (fip1*px[ip1        ]-c))
                            + wx010*((gjm1*px[    jm1    ]-c) + (gjp1*px[    jp1    ]-c))
                            + wx001*((gkm1*px[        km1]-c) + (gkp1*px[        kp1]-c))
                            + w2   *( gip1*fjm1*py[ip1+jm1] - gip1*fjp1*py[ip1+jp1] + 
                                      gim1*fjp1*py[im1+jp1] - gim1*fjm1*py[im1+jm1] + 
                                      gip1*fkm1*pz[ip1+km1] - gip1*fkp1*pz[ip1+kp1] + 
                                      gim1*fkp1*pz[im1+kp1] - gim1*fkm1*pz[im1+km1])
                            + (lam0*c
                            +  w110*((fim1*gjm1*px[im1+jm1    ]-c) + (fip1*gjm1*px[ip1+jm1    ]-c) + (fim1*gjp1*px[im1+jp1    ]-c) + (fip1*gjp1*px[ip1+jp1    ]-c))
                            +  w101*((fim1*gkm1*px[im1    +km1]-c) + (fip1*gkm1*px[ip1    +km1]-c) + (fim1*gkp1*px[im1    +kp1]-c) + (fip1*gkp1*px[ip1    +kp1]-c))
                            +  w011*((gjm1*gkm1*px[    jm1+km1]-c) + (gjp1*gkm1*px[    jp1+km1]-c) + (gjm1*gkp1*px[    jm1+kp1]-c) + (gjp1*gkp1*px[    jp1+kp1]-c))
                            +  w200*((fim2*px[im2        ]-c) + (fip2*px[ip2        ]-c))
                            +  w020*((gjm2*px[    jm2    ]-c) + (gjp2*px[    jp2    ]-c))
                            +  w002*((gkm2*px[        km2]-c) + (gkp2*px[        kp2]-c)))/v0);

                    c   = py[0];
                    suy = pby[i]
                          - ( wy100*((gim1*py[im1        ]-c) + (gip1*py[ip1        ]-c))
                            + wy010*((fjm1*py[    jm1    ]-c) + (fjp1*py[    jp1    ]-c))
                            + wy001*((gkm1*py[        km1]-c) + (gkp1*py[        kp1]-c))
                            + w2   *( gjp1*fim1*px[jp1+im1] - gjp1*fip1*px[jp1+ip1] + 
                                      gjm1*fip1*px[jm1+ip1] - gjm1*fim1*px[jm1+im1] + 
                                      gjp1*fkm1*pz[jp1+km1] - gjp1*fkp1*pz[jp1+kp1] + 
                                      gjm1*fkp1*pz[jm1+kp1] - gjm1*fkm1*pz[jm1+km1])
                            + (lam0*c
                            +  w110*((gim1*fjm1*py[im1+jm1    ]-c) + (gip1*fjm1*py[ip1+jm1    ]-c) + (gim1*fjp1*py[im1+jp1    ]-c) + (gip1*fjp1*py[ip1+jp1    ]-c))
                            +  w101*((gim1*gkm1*py[im1    +km1]-c) + (gip1*gkm1*py[ip1    +km1]-c) + (gim1*gkp1*py[im1    +kp1]-c) + (gip1*gkp1*py[ip1    +kp1]-c))
                            +  w011*((fjm1*gkm1*py[    jm1+km1]-c) + (fjp1*gkm1*py[    jp1+km1]-c) + (fjm1*gkp1*py[    jm1+kp1]-c) + (fjp1*gkp1*py[    jp1+kp1]-c))
                            +  w200*((gim2*py[im2        ]-c) + (gip2*py[ip2        ]-c))
                            +  w020*((fjm2*py[    jm2    ]-c) + (fjp2*py[    jp2    ]-c))
                            +  w002*((gkm2*py[        km2]-c) + (gkp2*py[        kp2]-c)))/v1);

                    c   = pz[0];
                    suz = pbz[i]
                          - ( wz100*((gim1*pz[im1        ]-c) + (gip1*pz[ip1        ]-c))
                            + wz010*((gjm1*pz[    jm1    ]-c) + (gjp1*pz[    jp1    ]-c))
                            + wz001*((fkm1*pz[        km1]-c) + (fkp1*pz[        kp1]-c))
                            + w2   *(gkp1*fim1*px[kp1+im1] - gkp1*fip1*px[kp1+ip1] + 
                                     gkm1*fip1*px[km1+ip1] - gkm1*fim1*px[km1+im1] + 
                                     gkp1*fjm1*py[kp1+jm1] - gkp1*fjp1*py[kp1+jp1] + 
                                     gkm1*fjp1*py[km1+jp1] - gkm1*fjm1*py[km1+jm1])
                            + (lam0*c
                            +  w110*((gim1*gjm1*pz[im1+jm1    ]-c) + (gip1*gjm1*pz[ip1+jm1    ]-c) + (gim1*gjp1*pz[im1+jp1    ]-c) + (gip1*gjp1*pz[ip1+jp1    ]-c))
                            +  w101*((gim1*fkm1*pz[im1    +km1]-c) + (gip1*fkm1*pz[ip1    +km1]-c) + (gim1*fkp1*pz[im1    +kp1]-c) + (gip1*fkp1*pz[ip1    +kp1]-c))
                            +  w011*((gjm1*fkm1*pz[    jm1+km1]-c) + (gjp1*fkm1*pz[    jp1+km1]-c) + (gjm1*fkp1*pz[    jm1+kp1]-c) + (gjp1*fkp1*pz[    jp1+kp1]-c))
                            +  w200*((gim2*pz[im2        ]-c) + (gip2*pz[ip2        ]-c))
                            +  w020*((gjm2*pz[    jm2    ]-c) + (gjp2*pz[    jp2    ]-c))
                            +  w002*((fkm2*pz[        km2]-c) + (fkp2*pz[        kp2]-c)))/v2);

                    if (a)
                    {
                        double axx, ayy, azz, axy, axz, ayz, idt;

                        sux -= (paxx[i]*px[0] + paxy[i]*py[0] + paxz[i]*pz[0]);
                        suy -= (paxy[i]*px[0] + payy[i]*py[0] + payz[i]*pz[0]);
                        suz -= (paxz[i]*px[0] + payz[i]*py[0] + pazz[i]*pz[0]);

                        axx  = paxx[i] + wx000;
                        ayy  = payy[i] + wy000;
                        azz  = pazz[i] + wz000;
                        axy  = paxy[i];
                        axz  = paxz[i];
                        ayz  = payz[i];
                        idt  = 1.0/(axx*ayy*azz -axx*ayz*ayz-ayy*axz*axz-azz*axy*axy +2*axy*axz*ayz);
                        *px += idt*(sux*(ayy*azz-ayz*ayz)+suy*(axz*ayz-axy*azz)+suz*(axy*ayz-axz*ayy));
                        *py += idt*(sux*(axz*ayz-axy*azz)+suy*(axx*azz-axz*axz)+suz*(axy*axz-axx*ayz));
                        *pz += idt*(sux*(axy*ayz-axz*ayy)+suy*(axy*axz-axx*ayz)+suz*(axx*ayy-axy*axy));
                    }
                    else
                    {
                        *px += sux/wx000;
                        *py += suy/wy000;
                        *pz += suz/wz000;
                    }
                }
            }
        }
#       ifdef VERBOSE
        if ((it%27) == 26)
            printf(" %g", sumsq(dm, a, b, s, u));
#       endif
    }
#   ifdef VERBOSE
        printf("\n");
#   endif
}

void relax(mwSize dm[], float a[], float b[], double s[], int nit, float u[])
{
    if (s[5]==0 && s[6]==0 && s[7]==0)
        relax_me(dm, a, b, s, nit, u);
    else if (s[6]==0 && s[7]==0)
        relax_be(dm, a, b, s, nit, u);
    else if (s[4]==0 && s[5]==0)
        relax_le(dm, a, b, s, nit, u);
    else
        relax_all(dm, a, b, s, nit, u);
}

/* Compute momentum from velocity.
 * It is equivalent to m = L * v (where L is the penalty operator), 
 * or to m = f (*) v (where f is the equivalent convolution kernel).
 *
 * \param dm     Dimensions of the velocity lattice (voxels)
 * \param f      Velocities computed on the above lattice
 * \param s[0:3] Voxel sizes (voxels/mm)
 * \param s[3]   Parameter of the absolute displacement penalty (lam0)
 * \param s[4]   Parameter of the membrane energy (lam1) (penalizes 
 *               elements of the Jacobian matrix -> 1st order smoothness)
 * \param s[5]   Parameter of the bending energy (lam2) (penalizes 
 *               elements of the Hessian matrix -> 2nd order smoothness)
 * \param s[6]   Parameter of the linear elastic energy (mu) (penalizes 
 *               elements of the symmetric part of the Jacobian matrix 
 *               -> penalizes scaling and shearing)
 * \param s[7]   Parameter of the linear elastic energy (lam) (penalizes 
 *               the divergence of the Jacobian matrix -> preserves 
 *               volumes)
 * \param g      Allocated array in which to store the output momentum 
 */
void vel2mom(mwSize dm[], float f[], double s[], float g[])
{
    if (s[5]==0 && s[6]==0 && s[7]==0)
        // Membrane energy only
        vel2mom_me(dm, f, s, g);
    else if (s[6]==0 && s[7]==0)
        // Membrane + Bending energy only
        vel2mom_be(dm, f, s, g);
    else if (s[4]==0 && s[5]==0)
        // Linear elastic energy only
        vel2mom_le(dm, f, s, g);
    else
        // Mixture of regularization energies
        vel2mom_all(dm, f, s, g);
}

/* Pointwise matrix multiplication in the momentum space.
 * Returns A[i] * vel2mom(p[i]), where A contains M symmetric 3*3 matrix 
 * (tensor field) and p contains M 3*1 vectors (velocity field).
 * 
 * \param[in]  dm  Dimension of the lattice (M = prod(dm))
 * \param[in]  A   Symmetric tensor field
 * \param[in]  p   Vector field
 * \param[out] Ap  Result of the pointwise matrix multiplication (A * p)
 */
void Atimesp(mwSize dm[], float A[], double param[], float p[], float Ap[])
{
    vel2mom(dm, p, param, Ap);
    Atimesp1(dm, A, p, Ap);
}



