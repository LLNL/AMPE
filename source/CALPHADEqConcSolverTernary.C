#include "CALPHADEqConcSolverTernary.h"
#include "CALPHADFunctions.h"

#include <iostream>
#include <cmath>
#include <cassert>

using namespace std;

//=======================================================================

void CALPHADEqConcentrationSolverTernary::RHS(
   const double* const c,
   double* const fvec )
{
   const double* const cL=&c[0];
   const double* const cA=&c[2];
   //tbox::pout<<"Compute RHS for CALPHAD..."<<endl;
   
   double derivIdealMixL[2];
   CALPHADcomputeFIdealMix_derivTernary( d_RT, cL[0], cL[1], derivIdealMixL );
   
   double derivFMixL[2];
   CALPHADcomputeFMix_derivTernary( d_L_AB_L, d_L_AC_L, d_L_BC_L, cL[0], cL[1], derivFMixL );
   
   double dfLdciL[2];
   //1st species
   dfLdciL[0] = d_fA[0] - d_fC[0]
       + derivFMixL[0]
       + derivIdealMixL[0];

   //2nd species
   dfLdciL[1] = d_fB[0] - d_fC[0]
       + derivFMixL[1]
       + derivIdealMixL[1];
   
   double derivIdealMixA[2];
   CALPHADcomputeFIdealMix_derivTernary( d_RT, cA[0], cA[1], derivIdealMixA );
   
   double derivFMixA[2];
   CALPHADcomputeFMix_derivTernary( d_L_AB_A, d_L_AC_A, d_L_BC_A, cA[0], cA[1], derivFMixA );
   
   double dfAdciA[2];
   //1st species
   dfAdciA[0] = d_fA[1] - d_fC[1]
       + derivFMixA[0]
       + derivIdealMixA[0];

   //2nd species
   dfAdciA[1] = d_fB[1] - d_fC[1]
       + derivFMixA[1]
       + derivIdealMixA[1];
   
   // equation fL-fA-(cL-cA)*dfL/dcL=0
   fvec[0] =
        cL[0] * d_fA[0] + cL[1] * d_fB[0] + ( 1.0 - cL[0]-cL[1] ) * d_fC[0]
      + CALPHADcomputeFIdealMixTernary( d_RT, cL[0], cL[1] )
      + CALPHADcomputeFMixTernary( d_L_AB_L, d_L_AC_L, d_L_BC_L, cL[0], cL[1] )
      - cA[0] * d_fA[1] - cA[1] * d_fB[1] - ( 1.0 - cA[0] - cA[1] ) * d_fC[1]
      - CALPHADcomputeFIdealMixTernary( d_RT, cA[0], cA[1] )
      - CALPHADcomputeFMixTernary( d_L_AB_A, d_L_AC_A, d_L_BC_A, cA[0], cA[1] )
      - (cL[0]-cA[0])* dfLdciL[0]
      - (cL[1]-cA[1])* dfLdciL[1];

   // equation: slope 0 in tangent plane for vector orthogonal to cL-CA
   fvec[1] = dfLdciL[0]*cL[1] - dfLdciL[1]*cA[1]
            -dfLdciL[0]*cL[0] + dfLdciL[1]*cA[0];


   fvec[2] = dfLdciL[0] - dfAdciA[0];
   fvec[3] = dfLdciL[1] - dfAdciA[1];
}

//=======================================================================

void CALPHADEqConcentrationSolverTernary::Jacobian(
   const double* const c,
   double** const fjac )
{
   //tbox::pout<<"Compute Jacobian for CALPHAD..."<<endl;
   const double* const cL=&c[0];
   const double* const cA=&c[2];
   //tbox::pout<<"Compute RHS for CALPHAD..."<<endl;
   
   double derivIdealMixL[2];
   CALPHADcomputeFIdealMix_derivTernary( d_RT, cL[0], cL[1], derivIdealMixL );
   
   double derivFMixL[2];
   CALPHADcomputeFMix_derivTernary( d_L_AB_L, d_L_AC_L, d_L_BC_L, cL[0], cL[1], derivFMixL );
   
   double deriv2IdealMixL[2];
   CALPHADcomputeFIdealMix_deriv2Ternary( d_RT, cL[0], cL[1], deriv2IdealMixL );
   
   double deriv2FMixL[2];
   CALPHADcomputeFMix_deriv2Ternary( d_L_AB_L, d_L_AC_L, d_L_BC_L, cL[0], cL[1], deriv2FMixL );
   
   double dfLdciL[2];
   double d2fLdciL2[2];
   //1st species
   dfLdciL[0] = d_fA[0] - d_fC[0]
       + derivFMixL[0]
       + derivIdealMixL[0];
   d2fLdciL2[0] =
         deriv2FMixL[0]
       + deriv2IdealMixL[0];

   //2nd species
   dfLdciL[1] = d_fB[0] - d_fC[0]
       + derivFMixL[1]
       + derivIdealMixL[1];
   d2fLdciL2[1] =
         deriv2FMixL[1]
       + deriv2IdealMixL[1];
   
   double derivIdealMixA[2];
   CALPHADcomputeFIdealMix_derivTernary( d_RT, cA[0], cA[1], derivIdealMixA );
   
   double derivFMixA[2];
   CALPHADcomputeFMix_derivTernary( d_L_AB_A, d_L_AC_A, d_L_BC_A, cA[0], cA[1], derivFMixA );
   
   double deriv2IdealMixA[2];
   CALPHADcomputeFIdealMix_deriv2Ternary( d_RT, cA[0], cA[1], deriv2IdealMixA );
   
   double deriv2FMixA[2];
   CALPHADcomputeFMix_deriv2Ternary( d_L_AB_A, d_L_AC_A, d_L_BC_A, cA[0], cA[1], deriv2FMixA );
   
   double dfAdciA[2];
   double d2fAdciA2[2];
   //1st species
   dfAdciA[0] = d_fA[1] - d_fC[1]
       + derivFMixA[0]
       + derivIdealMixA[0];
   d2fAdciA2[0] =
         deriv2FMixA[0]
       + deriv2IdealMixA[0];

   //2nd species
   dfAdciA[1] = d_fB[1] - d_fC[1]
       + derivFMixA[1]
       + derivIdealMixA[1];
   d2fAdciA2[1] =
         deriv2FMixA[1]
       + deriv2IdealMixA[1];


   // f[i][j]=df[i]/dc[j]
   fjac[0][0] = 
        d_fA[0] - d_fC[0]
      + derivIdealMixL[0]
      + derivFMixL[0]
      - dfLdciL[0];

   fjac[0][1] = 
        d_fB[0] - d_fC[0]
      + derivIdealMixL[1]
      + derivFMixL[1]
      -  dfLdciL[1];
   
   fjac[0][2] = 
      - d_fA[1] + d_fC[1]
      - derivIdealMixA[0]
      - derivFMixA[0]
      + dfLdciL[0];
 
   fjac[0][3] = 
      - d_fB[1] + d_fC[1]
      - derivIdealMixA[1]
      - derivFMixA[1]
      + dfLdciL[1];

   fjac[1][0] = d2fLdciL2[0]*cL[1]
               -d2fLdciL2[0]*cL[0] -dfLdciL[0];

   fjac[1][1] = dfLdciL[0] - d2fLdciL2[1]*cA[1]
              + d2fLdciL2[1]*cA[0];

   fjac[1][2] = dfLdciL[1];
   
   fjac[1][3] = - dfLdciL[1];
   
   
   fjac[2][0] = d2fLdciL2[0];
   fjac[2][1] = 0.;
   fjac[2][2] = - d2fAdciA2[0];
   fjac[2][3] = 0.;
   
   fjac[3][0] = 0.;
   fjac[3][1] = d2fLdciL2[1];
   fjac[3][2] = 0.;
   fjac[3][3] = - d2fAdciA2[1];
}

//=======================================================================

int CALPHADEqConcentrationSolverTernary::ComputeConcentration(
   double* const conc,
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
   d_RTinv = RTinv;
   d_RT    = 1./RTinv;

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

   return DampedNewtonSolver::ComputeSolution( conc, 4 );
}
