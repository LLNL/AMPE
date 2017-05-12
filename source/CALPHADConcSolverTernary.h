#ifndef included_CALPHADConcSolverTernary
#define included_CALPHADConcSolverTernary

#include "DampedNewtonSolver.h"

class CALPHADConcentrationSolverTernary :
   public DampedNewtonSolver
{
public :

   CALPHADConcentrationSolverTernary();
      
   virtual ~CALPHADConcentrationSolverTernary() {};
      
   int ComputeConcentration(
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
      const double* const fC );
   
protected:

   //energies of 3 species, in two phase each
   double d_fA[2];
   double d_fB[2];
   double d_fC[2];

   // L coefficients for phase L
   double d_L_AB_L[4];
   double d_L_AC_L[4];
   double d_L_BC_L[4];
   
   // L coefficients for phase A
   double d_L_AB_A[4];
   double d_L_AC_A[4];
   double d_L_BC_A[4];
   
   
   double d_RTinv;

   void computeXi(const double* const cL, const double* const cA, 
                  double xiL[2], double xiA[2])const;

   void computeDxiDc(const double* const cL, const double* const cA,
                     double xiL[2], double xiA[2], 
                     double dxidcL[2], double dxidcA[2])const;
   
private :

   void RHS(
      const double* const c,
      double* const fvec );

   void Jacobian(
      const double* const c,
      double** const fjac );
   
   double d_c0[2];
   double d_hphi;
};

#endif
