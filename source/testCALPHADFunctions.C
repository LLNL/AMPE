#include "CALPHADFunctions.h"

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;


int main( int argc, char *argv[] )
{
   cout<<"Test CALPHAD functions."<<endl;

   const double epsilon=1.e-8;
   const double tol=1.e-6;
   cout<<setprecision(12);
   cerr<<setprecision(12);

   {
   double l0=2.3;
   double l1=5.1;
   double l2=3.2;
   double l3=-2.5;
   double c=0.33;

   double f0=CALPHADcomputeFMixBinary(l0,l1,l2,l3,c);
   double f1=CALPHADcomputeFMixBinary(l0,l1,l2,l3,c+epsilon);

   double fd=(f1-f0)/epsilon;
   cout<<"Numerical derivative   = "<<fd<<endl;
   double ad=CALPHADcomputeFMix_derivBinary(l0,l1,l2,l3,c);
   cout<<"Analytical derivative = "<<ad<<endl;

   
   if( fabs(fd-ad)>tol )
   {
      cerr<<"TEST: CALPHADcomputeFMix_derivBinary failed!!!"<<endl;
      cerr<<"Difference = "<<fd-ad<<endl;
   }else{
      cout<<"TEST successful!"<<endl;
   }
   }

   {
   double lAB[4]={1.2,3.5,6.1,-1.2};
   double lAC[4]={3.2,4.5,-3.1,7.1};
   double lBC[4]={0.2,0.7,4.1,-8.2};

   double cA=0.33;
   double cB=0.21;

   double f0=CALPHADcomputeFMixTernary(lAB,lAC,lBC,cA,cB);
   double f1=CALPHADcomputeFMixTernary(lAB,lAC,lBC,cA+epsilon,cB);
   double f2=CALPHADcomputeFMixTernary(lAB,lAC,lBC,cA,cB+epsilon);

   double deriv[2];
   CALPHADcomputeFMix_derivTernary(lAB,lAC,lBC,cA,cB,deriv);

   double fd1=(f1-f0)/epsilon;
   if( fabs(fd1-deriv[0])>tol )
   {
      cerr<<"TEST: CALPHADcomputeFMix_derivTernary failed!!!"<<endl;
      cerr<<"Difference = "<<fd1-deriv[0]<<endl;
   }else{
      cout<<"TEST CALPHADcomputeFMix_derivTernary, c0, successful!"<<endl;
   }
   double fd2=(f2-f0)/epsilon;
   if( fabs(fd2-deriv[1])>tol )
   {
      cerr<<"TEST: CALPHADcomputeFMix_derivTernary failed!!!"<<endl;
      cerr<<"Difference = "<<fd2-deriv[1]<<endl;
   }else{
      cout<<"TEST CALPHADcomputeFMix_derivTernary, c1, successful!"<<endl;
   }

   double deriv2[4];
   CALPHADcomputeFMix_deriv2Ternary(lAB,lAC,lBC,cA,cB,deriv2);

   double fderiv0[2];
   CALPHADcomputeFMix_derivTernary(lAB,lAC,lBC,cA+epsilon,cB,fderiv0);
   double fderiv1[2];
   CALPHADcomputeFMix_derivTernary(lAB,lAC,lBC,cA,cB+epsilon,fderiv1);

   double fd=(fderiv0[0]-deriv[0])/epsilon;
   cout<<"FD="<<fd<<", exact="<<deriv2[0]<<endl;
   if( fabs(fd-deriv2[0])>tol )
   {
      cerr<<"TEST: CALPHADcomputeFMix_deriv2Ternary failed!!!"<<endl;
      cerr<<"Difference = "<<fd-deriv2[0]<<endl;
   }else{
      cout<<"TEST CALPHADcomputeFMix_deriv2Ternary, c0, c0, successful!"<<endl;
   }

   fd=(fderiv1[1]-deriv[1])/epsilon;
   cout<<"FD="<<fd<<", exact="<<deriv2[3]<<endl;
   if( fabs(fd-deriv2[3])>tol )
   {
      cerr<<"TEST: CALPHADcomputeFMix_deriv2Ternary failed!!!"<<endl;
      cerr<<"Difference = "<<fd-deriv2[3]<<endl;
   }else{
      cout<<"TEST CALPHADcomputeFMix_deriv2Ternary, c1, c1, successful!"<<endl;
   }

   //cross derivatives
   fd=(fderiv0[1]-deriv[1])/epsilon;
   cout<<"FD="<<fd<<", exact="<<deriv2[1]<<endl;
   if( fabs(fd-deriv2[1])>tol )
   {
      cerr<<"TEST: CALPHADcomputeFMix_deriv2Ternary failed!!!"<<endl;
      cerr<<"Difference = "<<fd-deriv2[0]<<endl;
   }else{
      cout<<"TEST CALPHADcomputeFMix_deriv2Ternary, c0, c1, successful!"<<endl;
   }

   fd=(fderiv1[0]-deriv[0])/epsilon;
   cout<<"FD="<<fd<<", exact="<<deriv2[2]<<endl;
   if( fabs(fd-deriv2[2])>tol )
   {
      cerr<<"TEST: CALPHADcomputeFMix_deriv2Ternary failed!!!"<<endl;
      cerr<<"Difference = "<<fd-deriv2[0]<<endl;
   }else{
      cout<<"TEST CALPHADcomputeFMix_deriv2Ternary, c1, c0, successful!"<<endl;
   }




   

   }

   return(0);
}
