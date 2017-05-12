#include "CALPHADConcSolverTernary.h"
#include "CALPHADFunctions.h"

#include <iostream>
#include <cmath>
#include <cassert>

using namespace std;

//=======================================================================

CALPHADConcentrationSolverTernary::CALPHADConcentrationSolverTernary()
{
}

//=======================================================================

void CALPHADConcentrationSolverTernary::computeXi(const double* const cL, const double* const cA, 
                                                  double xiL[2], double xiA[2])const
{
   //std::cout<<"CALPHADConcentrationSolverTernary::computeXi()"<<endl;

   xiL[0]=0.;
   xiL[1]=0.;
   CALPHADcomputeFMix_derivTernary(d_L_AB_L, d_L_AC_L, d_L_BC_L, cL[0], cL[1], &xiL[0]);
   xiL[0] += (d_fA[0] - d_fC[0]);
   xiL[1] += (d_fB[0] - d_fC[0]);
   
   xiL[0] *= d_RTinv;
   xiL[1] *= d_RTinv;

   xiA[0]=0.;
   xiA[1]=0.;
   CALPHADcomputeFMix_derivTernary(d_L_AB_A, d_L_AC_A, d_L_BC_A, cA[0], cA[1], &xiA[0]);
   xiA[0] += (d_fA[1] - d_fC[1]);
   xiA[1] += (d_fB[1] - d_fC[1]);
   
   xiA[0] *= d_RTinv;
   xiA[1] *= d_RTinv;
}
//=======================================================================

// solve for c=(c_L, c_A)
void CALPHADConcentrationSolverTernary::RHS(
   const double* const c,
   double* const fvec )
{
   const double* const cL=&c[0];
   const double* const cA=&c[2];
   
   double xiL[2]={0.,0.};
   double xiA[2]={0.,0.};

   computeXi(cL,cA,xiL,xiA);

   //equation for 1st species
   fvec[0] =
      -d_c0[0] + ( 1.0 - d_hphi ) * cL[0] + 
      d_hphi * cA[0];

   //equation for 2nd species
   fvec[1] =
      -d_c0[1] + ( 1.0 - d_hphi ) * cL[1] + 
      d_hphi * cA[1];

   //equal chemical potential equation for 1st species
   fvec[2] = xlogx_deriv(cL[0])-xlogx_deriv(1.-cL[0]-cL[1])
            -xlogx_deriv(cA[0])+xlogx_deriv(1.-cA[0]-cA[1])
            +(xiL[0] - xiA[0] );

   //equal chemical potential equation for 2nd species
   fvec[3] = xlogx_deriv(cL[1])-xlogx_deriv(1.-cL[0]-cL[1])
            -xlogx_deriv(cA[1])+xlogx_deriv(1.-cA[0]-cA[1])
            +(xiL[1] - xiA[1] );

}

//=======================================================================

void CALPHADConcentrationSolverTernary::computeDxiDc(const double* const cL, const double* const cA,
                                                     double xiL[2], double xiA[2], 
                                                     double dxidcL[2], double dxidcA[2])const
{
   //std::cout<<"CALPHADConcentrationSolverTernary::computeDxiDc()"<<endl;
   computeXi(cL,cA,xiL,xiA);
   
   CALPHADcomputeFMix_deriv2Ternary(d_L_AB_L, d_L_AC_L, d_L_BC_L,cL[0],cL[1],dxidcL);
   CALPHADcomputeFMix_deriv2Ternary(d_L_AB_A, d_L_AC_A, d_L_BC_A,cA[0],cA[1],dxidcA);


   for(short ii=0;ii<2;ii++){
      dxidcL[ii] *= d_RTinv;
      dxidcA[ii] *= d_RTinv;
   }
}

//=======================================================================

void CALPHADConcentrationSolverTernary::Jacobian(
   const double* const c,
   double** const fjac )
{
   const double* const cL=&c[0];
   const double* const cA=&c[2];

   double xiL[2];
   double xiA[2];
   double dxidcL[2];
   double dxidcA[2];
   double expxi[2];
   
   computeDxiDc(cL,cA,xiL,xiA,dxidcL,dxidcA);

   fjac[0][0] = ( 1.0 - d_hphi );
   fjac[0][1] = d_hphi;
   fjac[0][2] = 0.;
   fjac[0][3] = 0.;

   fjac[1][0] = 0.;
   fjac[1][1] = 0.;
   fjac[1][2] = ( 1.0 - d_hphi );
   fjac[1][3] = d_hphi;

   fjac[2][0] =  dxidcL[0]+xlogx_deriv2(cL[0])+xlogx_deriv2(1.-cL[0]-cL[1]);
   fjac[2][1] = -dxidcA[0]+xlogx_deriv2(cA[0])+xlogx_deriv2(1.-cA[0]-cA[1]);
   fjac[2][2] =            xlogx_deriv2(cL[1])+xlogx_deriv2(1.-cL[0]-cL[1]);
   fjac[2][3] =            xlogx_deriv2(cA[1])+xlogx_deriv2(1.-cA[0]-cA[1]);

   fjac[3][0] =            xlogx_deriv2(cL[0])+xlogx_deriv2(1.-cL[0]-cL[1]);
   fjac[3][1] =            xlogx_deriv2(cA[0])+xlogx_deriv2(1.-cA[0]-cA[1]);
   fjac[3][0] =  dxidcL[1]+xlogx_deriv2(cL[1])+xlogx_deriv2(1.-cL[0]-cL[1]);
   fjac[3][1] = -dxidcA[1]+xlogx_deriv2(cA[1])+xlogx_deriv2(1.-cA[0]-cA[1]);
}

//=======================================================================

int CALPHADConcentrationSolverTernary::ComputeConcentration(
   double* const conc,
   const double c0,
   const double c1,
   const double hphi,
   const double RTinv,
   const double* const L_AB_L,
   const double* const L_AC_L,
   const double* const L_BC_L,
   const double* const L_AB_A,
   const double* const L_AC_A,
   const double* const L_BC_A,
   const double* const fA,
   const double* const fB,
   const double* const fC )
{
   //std::cout<<"CALPHADConcentrationSolverTernary::ComputeConcentration()"<<endl;
   d_c0[0] = c0;
   d_c0[1] = c1;
   d_hphi = hphi;
   d_RTinv = RTinv;

   for ( int ii = 0; ii < 4; ii++ ) {
      d_L_AB_L[ii] = L_AB_L[ii];
      d_L_AC_L[ii] = L_AC_L[ii];
      d_L_BC_L[ii] = L_BC_L[ii];
   }
   for ( int ii = 0; ii < 4; ii++ ) {
      d_L_AB_A[ii] = L_AB_A[ii];
      d_L_AC_A[ii] = L_AC_A[ii];
      d_L_BC_A[ii] = L_BC_A[ii];
   }
   for ( int ii = 0; ii < 2; ii++ ) {
      d_fA[ii] = fA[ii];
      d_fB[ii] = fB[ii];
      d_fC[ii] = fC[ii];
   }

   int ret= NewtonSolver::ComputeSolution( conc, 4 );
   return ret;
}
