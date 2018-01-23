#include "CALPHADEqPhaseConcSolverTernary.h"
#include "CALPHADConcSolverTernary.h"

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;


int main( int argc, char *argv[] )
{
   cout<<"====================================="<<endl;
   cout<<"Test CALPHADEqPhaseConcSolverTernary."<<endl;

   double epsilon=1.e-8;
   double tol=1.e-6;

   double cA=0.1;
   double cB=0.2;
   CALPHADEqPhaseConcentrationSolverTernary solver(cA,cB);

   //energies of 3 species, in two phase each
   double fA[2]={2.3,4.5};
   double fB[2]={0.5,2.6};
   double fC[2]={0.9,3.1};

   // L coefficients for 2 possible phases (L and S)
   double L_AB_L[4]={2.,3.,4.,5.};
   double L_AC_L[4]={8.,5.,1.,3.};
   double L_BC_L[4]={6.,4.,9.,2.};

   double L_AB_S[4]={2.,2.,3.,4.};
   double L_AC_S[4]={6.,8.,5.,1.};
   double L_BC_S[4]={2.,6.,4.,9.};

   double L_ABC_L[3]={3.1,4.1,5.1};
   double L_ABC_S[3]={2.1,3.1,4.1};

   double RTinv=10.;

   solver.setup(RTinv,
                L_AB_L,L_AC_L,L_BC_L,
                L_AB_S,L_AC_S,L_BC_S,
                L_ABC_L,L_ABC_S,
                fA,fB,fC);

   double fvec1[5];
   double fvec2[5];
   double* fjac[5];
   for(int i=0;i<5;i++)fjac[i]=new double[5];

   double x[5]={0.1,0.2,0.3,0.4,0.5};

   solver.RHS(x,fvec1);

   solver.Jacobian(x,fjac);

   cerr<<setprecision(12);

   //loop over variables
   for(int j=0;j<5;j++)
   {
      cout<<"----------------------------"<<endl;
      cout<<"Test variations of variable "<<j<<endl;
      x[j]+=epsilon;
      solver.RHS(x,fvec2);

      double fd[5];
      //loop over equations
      for(int i=0;i<5;i++)fd[i]=(fvec2[i]-fvec1[i])/epsilon;

      for(int i=0;i<5;i++)
      {
         if( fabs(fd[i]-fjac[i][j])>tol )
         {
            cerr<<"ERROR: Equation "<<i<<", FD="<<fd[i]<<", fjac="<<fjac[i][j]<<endl;
         }
      }

      x[j]-=epsilon;
   }

   cout<<"====================================="<<endl;
   cout<<"Test CALPHADConcentrationSolverTernary"<<endl;
   CALPHADConcentrationSolverTernary solver2;

   solver2.setup(cA,cB,0.5,
                RTinv,
                L_AB_L,L_AC_L,L_BC_L,
                L_AB_S,L_AC_S,L_BC_S,
                L_ABC_L,L_ABC_S,
                fA,fB,fC);

   solver2.RHS(x,fvec1);

   solver2.Jacobian(x,fjac);

   cerr<<setprecision(12);

   //loop over variables
   for(int j=0;j<4;j++)
   {
      cout<<"----------------------------"<<endl;
      cout<<"Test variations of variable "<<j<<endl;
      x[j]+=epsilon;
      solver2.RHS(x,fvec2);

      double fd[4];
      //loop over equations
      for(int i=0;i<4;i++)fd[i]=(fvec2[i]-fvec1[i])/epsilon;

      for(int i=0;i<4;i++)
      {
         if( fabs(fd[i]-fjac[i][j])>tol )
         {
            cerr<<"ERROR: Equation "<<i<<", FD="<<fd[i]<<", fjac="<<fjac[i][j]<<endl;
         }
      }

      x[j]-=epsilon;
   }

   return(0);
}

