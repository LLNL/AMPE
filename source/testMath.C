#include "math_utilities.h"

#include <iostream>
#include <math.h>

using namespace std;


int main( int argc, char *argv[] )
{
   cout<<"Test determinant computation."<<endl;

   double* mat[4];
   double work[16]={11,2,3,4,5,16,7,8,9,10,11,12,13,14,15,16};
   for(short i=0;i<4;++i)mat[i]=&work[4*i];

   cout<<"Test function Determinant4..."<<endl;
   double d=Determinant4(mat); 

   const double tol=1.e-8;
   if( fabs(d+400.)>tol )
   {
      cerr<<"TEST: Determinant of 4x4 matrix failed!!!"<<endl;
   }else{
      cout<<"TEST successful!"<<endl;
   }
   
   cout<<"Test function DeterminantN..."<<endl;
   d=DeterminantN(mat,4);
   if( fabs(d+400.)>tol )
   {
      cerr<<"TEST: Determinant of 4x4 matrix failed!!!"<<endl;
      cerr<<"computed d="<<d<<endl;
      cerr<<"exact d="<<400.<<endl;
   }else{
      cout<<"TEST successful!"<<endl;
   }

   return(0);
}
