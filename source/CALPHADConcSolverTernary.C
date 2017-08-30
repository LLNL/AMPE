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

void CALPHADConcentrationSolverTernary::computeXi(const double* const cL, const double* const cS, 
                                                  double xiL[2], double xiS[2])const
{
   //std::cout<<"CALPHADConcentrationSolverTernary::computeXi()"<<endl;

   xiL[0]=0.;
   xiL[1]=0.;
   CALPHADcomputeFMix_derivTernary(d_L_AB_L, d_L_AC_L, d_L_BC_L, cL[0], cL[1], &xiL[0]);
   xiL[0] += (d_fA[0] - d_fC[0]);
   xiL[1] += (d_fB[0] - d_fC[0]);
   
   xiL[0] *= d_RTinv;
   xiL[1] *= d_RTinv;

   xiS[0]=0.;
   xiS[1]=0.;
   CALPHADcomputeFMix_derivTernary(d_L_AB_S, d_L_AC_S, d_L_BC_S, cS[0], cS[1], &xiS[0]);
   xiS[0] += (d_fA[1] - d_fC[1]);
   xiS[1] += (d_fB[1] - d_fC[1]);
   
   xiS[0] *= d_RTinv;
   xiS[1] *= d_RTinv;
}
//=======================================================================

// solve for c=(c_L, c_A)
void CALPHADConcentrationSolverTernary::RHS(
   const double* const c,
   double* const fvec )
{
   const double* const cL=&c[0];
   const double* const cS=&c[2];
   
   double xiL[2]={0.,0.};
   double xiS[2]={0.,0.};

   computeXi(cL,cS,xiL,xiS);

   //equation for 1st species
   fvec[0] =
      -d_c0[0] + ( 1.0 - d_hphi ) * cL[0] + 
      d_hphi * cS[0];

   //equation for 2nd species
   fvec[1] =
      -d_c0[1] + ( 1.0 - d_hphi ) * cL[1] + 
      d_hphi * cS[1];

   //equal chemical potential equation for 1st species
   fvec[2] = xlogx_deriv(cL[0])-xlogx_deriv(1.-cL[0]-cL[1])
            -xlogx_deriv(cS[0])+xlogx_deriv(1.-cS[0]-cS[1])
            +(xiL[0] - xiS[0] );

   //equal chemical potential equation for 2nd species
   fvec[3] = xlogx_deriv(cL[1])-xlogx_deriv(1.-cL[0]-cL[1])
            -xlogx_deriv(cS[1])+xlogx_deriv(1.-cS[0]-cS[1])
            +(xiL[1] - xiS[1] );

}

//=======================================================================

void CALPHADConcentrationSolverTernary::computeDxiDc(const double* const cL, const double* const cS,
                                                     double dxidcL[4], double dxidcS[4])const
{
   //std::cout<<"CALPHADConcentrationSolverTernary::computeDxiDc()"<<endl;
   
   CALPHADcomputeFMix_deriv2Ternary(d_L_AB_L, d_L_AC_L, d_L_BC_L,cL[0],cL[1],dxidcL);
   CALPHADcomputeFMix_deriv2Ternary(d_L_AB_S, d_L_AC_S, d_L_BC_S,cS[0],cS[1],dxidcS);

   for(short ii=0;ii<4;ii++){
      dxidcL[ii] *= d_RTinv;
      dxidcS[ii] *= d_RTinv;
   }
}

//=======================================================================

void CALPHADConcentrationSolverTernary::Jacobian(
   const double* const c,
   double** const fjac )
{
   const double* const cL=&c[0];
   const double* const cS=&c[2];

   double dxidcL[4];
   double dxidcS[4];
   double expxi[2];
   
   computeDxiDc(cL,cS,dxidcL,dxidcS);

   fjac[0][0] = ( 1.0 - d_hphi );
   fjac[0][1] = 0.;
   fjac[0][2] = d_hphi;
   fjac[0][3] = 0.;

   fjac[1][0] = 0.;
   fjac[1][1] = ( 1.0 - d_hphi );
   fjac[1][2] = 0.;
   fjac[1][3] = d_hphi;

   fjac[2][0] =  dxidcL[0]+xlogx_deriv2(cL[0])+xlogx_deriv2(1.-cL[0]-cL[1]);
   fjac[2][1] =  dxidcL[1]+                    xlogx_deriv2(1.-cL[0]-cL[1]);
   fjac[2][2] = -dxidcS[0]+xlogx_deriv2(cS[0])+xlogx_deriv2(1.-cS[0]-cS[1]);
   fjac[2][3] = -dxidcS[1]+                    xlogx_deriv2(1.-cS[0]-cS[1]);

   fjac[3][0] =  dxidcL[1]                    +xlogx_deriv2(1.-cL[0]-cL[1]);
   fjac[3][1] =  dxidcL[3]+xlogx_deriv2(cL[1])+xlogx_deriv2(1.-cL[0]-cL[1]);
   fjac[3][2] = -dxidcS[1]                    +xlogx_deriv2(1.-cS[0]-cS[1]);
   fjac[3][3] = -dxidcS[3]+xlogx_deriv2(cS[1])+xlogx_deriv2(1.-cS[0]-cS[1]);
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
   const double* const L_AB_S,
   const double* const L_AC_S,
   const double* const L_BC_S,
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
      d_L_AB_S[ii] = L_AB_S[ii];
      d_L_AC_S[ii] = L_AC_S[ii];
      d_L_BC_S[ii] = L_BC_S[ii];
   }
   for ( int ii = 0; ii < 2; ii++ ) {
      d_fA[ii] = fA[ii];
      d_fB[ii] = fB[ii];
      d_fC[ii] = fC[ii];
   }

   int ret= NewtonSolver::ComputeSolution( conc, 4 );
   return ret;
}
