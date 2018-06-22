// Copyright (c) 2018, Lawrence Livermore National Security, LLC.
// Produced at the Lawrence Livermore National Laboratory
// Written by M.R. Dorr, J.-L. Fattebert and M.E. Wickett
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE. 
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
// - Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the disclaimer below.
// - Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the disclaimer (as noted below) in the
//   documentation and/or other materials provided with the distribution.
// - Neither the name of the LLNS/LLNL nor the names of its contributors may be
//   used to endorse or promote products derived from this software without
//   specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
// LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
// OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
// 
#include "CALPHADSpeciesPhaseGibbsEnergy.h"
#include "SAMRAI/tbox/Utilities.h"

#include <cassert>

using namespace std;

void CALPHADSpeciesPhaseGibbsEnergy::initialize(const string& name, 
   boost::shared_ptr<tbox::Database> db)
{
   d_name=name;
   size_t ntc=db->getArraySize("Tc");
   assert( ntc>1 );   
   d_tc.resize(ntc);
   db->getDoubleArray("Tc",&d_tc[0],ntc);

   size_t nintervals=db->getArraySize("a");
   assert( nintervals==ntc-1 );
   d_a.resize(nintervals);
   db->getDoubleArray("a",&d_a[0],nintervals);
   
   nintervals=db->getArraySize("b");
   d_b.resize(nintervals);
   db->getDoubleArray("b",&d_b[0],nintervals);
   
   nintervals=db->getArraySize("c");
   d_c.resize(nintervals);
   db->getDoubleArray("c",&d_c[0],nintervals);
   
   nintervals=db->getArraySize("d2");
   d_d2.resize(nintervals);
   db->getDoubleArray("d2",&d_d2[0],nintervals);
   
   if ( db->keyExists( "d3" ) ) {
      nintervals=db->getArraySize("d3");
      d_d3.resize(nintervals);
      db->getDoubleArray("d3",&d_d3[0],nintervals);
   }else{
      d_d3.resize(nintervals);
      for(unsigned i=0;i<nintervals;i++)d_d3[i]=0.;
   }
   
   if ( db->keyExists( "d4" ) ) {
      nintervals=db->getArraySize("d4");
      d_d4.resize(nintervals);
      db->getDoubleArray("d4",&d_d4[0],nintervals);
   }
   
   if ( db->keyExists( "d7" ) ) {
      nintervals=db->getArraySize("d7");
      d_d7.resize(nintervals);
      db->getDoubleArray("d7",&d_d7[0],nintervals);
   }

   if ( db->keyExists( "dm1" ) ) {
      nintervals=db->getArraySize("dm1");
      d_dm1.resize(nintervals);
      db->getDoubleArray("dm1",&d_dm1[0],nintervals);
   }else{
      d_dm1.resize(nintervals);
      for(unsigned i=0;i<nintervals;i++)d_dm1[i]=0.;
   }
   
   if ( db->keyExists( "dm9" ) ) {
      nintervals=db->getArraySize("dm9");
      d_dm9.resize(nintervals);
      db->getDoubleArray("dm9",&d_dm9[0],nintervals);
   }
   
   if ( db->keyExists( "d5" ) ) {
      TBOX_ERROR( "CALPHADSpeciesPhaseGibbsEnergy: T**5 not implemented!!!" << endl );
   }
   if ( db->keyExists( "d6" ) ) {
      TBOX_ERROR( "CALPHADSpeciesPhaseGibbsEnergy: T**6 not implemented!!!" << endl );
   }
   if ( db->keyExists( "dm2" ) ) {
      TBOX_ERROR( "CALPHADSpeciesPhaseGibbsEnergy: T**-2 not implemented!!!" << endl );
   }

   d_test_fenergy_done= false;
}

/////////////////////////////////////////////////////////////////////
// Free energy function
// parameters are in J/mol
// returned values are in J/mol
double CALPHADSpeciesPhaseGibbsEnergy::fenergy(const double T) // expect T in Kelvin
{
   int index=-1;
   const int n=(int)d_tc.size();
   //tbox::pout<<"n="<<n<<endl;
   assert( n>1 );
   
   for(int i=0;i<n-1;i++)
   if( T>=d_tc[i] && T<d_tc[i+1] ){
      index=i;
   }
   
   if( index<0 ){
      cout<<"T="<<T<<", Tmin="<<d_tc[0]
                   <<", Tmax="<<d_tc[n-1]<<endl;
      TBOX_ERROR( "T out of range for fenergy" << endl );
   }

#if 0
   cout<<"a="<<d_a[index]<<endl;
   cout<<"b="<<d_b[index]<<endl;
   cout<<"c="<<d_c[index]<<endl;
   cout<<"d2="<<d_d2[index]<<endl;
   cout<<"d3="<<d_d3[index]<<endl;
   cout<<"dm1="<<d_dm1[index]<<endl;
#endif
   
   const double t2=T*T;
   double f= (d_a[index]+d_b[index]*T+d_c[index]*T*log(T)
             +d_d2[index]*t2
             +d_d3[index]*t2*T
             +d_dm1[index]/T);
   if( d_d4.size()>0 )
      f += d_d4[index]*t2*t2;
   if( d_d7.size()>0 )
      f += d_d7[index]*t2*t2*t2*T;
   if( d_dm9.size()>0 )
      f += d_dm9[index]/pow(T,9.);

   return f;
}

void CALPHADSpeciesPhaseGibbsEnergy::plotFofT(std::ostream& os, 
   const double T0, const double T1)
{
   const double dT=10.;
   const int npts=(int)trunc((T1-T0)/dT);
   os<<"# fenergy(J/mol) vs. T(K) for species "<<d_name<<endl;
   for(int i=0;i<npts;i++)
   {
       double testT=T0+dT*i;
       os<<testT<<'\t';
       os<<fenergy(testT)<<endl;
   }
   os<<endl;
}
