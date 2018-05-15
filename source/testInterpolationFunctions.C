#include "FuncFort.h"

#include <iostream>
#include <math.h>

using namespace std;


int main( int argc, char *argv[] )
{
   cout<<"Test interpolation functions..."<<endl;

   const double tol = 1.e-8;

   {
      std::string interp_func_type = "pbg";

      double val = FORT_INTERP_FUNC( 0.5,  interp_func_type.c_str() );
      if( fabs( val-0.5 )>tol ){
         cerr<<"Test failed!"<<endl;
         return 1;
      }

      val = FORT_INTERP_FUNC( -0.5,  interp_func_type.c_str() );
      if( fabs( val-0. )>tol ){
         cerr<<"Test failed!"<<endl;
         return 1;
      }
   }

   /*
    * Ratio pbg/lin
    */
   {
      std::string interp_func_type1 = "pbg";
      std::string interp_func_type2 = "lin";
      double phi=0.05;
      double val1 = FORT_INTERP_RATIO_FUNC( phi, interp_func_type1.c_str(),
                                                 interp_func_type2.c_str() );
      double val2 = phi*phi*(10.-15.*phi+6*phi*phi);
      if( fabs( val1-val2 )>tol ){
         cerr<<"Test failed!"<<endl;
         return 1;
      }

   }

   {
      std::string interp_func_type1 = "pbg";
      std::string interp_func_type2 = "lin";
      double phi=0.05;
      double val1 = FORT_COMPL_INTERP_RATIO_FUNC( phi,
                                                  interp_func_type1.c_str(),
                                                  interp_func_type2.c_str() );
      double a = 1.-FORT_INTERP_FUNC( phi, interp_func_type1.c_str() );
      double b = 1.-FORT_INTERP_FUNC( phi, interp_func_type2.c_str() );

      if( fabs( val1-(a/b) )>tol ){
         cerr<<"Test failed!"<<endl;
         return 1;
      }

   }

   cout<<"TEST successful!"<<endl;

   return(0);
}
