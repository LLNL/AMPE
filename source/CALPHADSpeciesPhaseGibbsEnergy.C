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

   const size_t nintervals=db->getArraySize("a");
   assert( nintervals==ntc-1 );
   vector<double> a(nintervals);
   db->getDoubleArray("a",&a[0],nintervals);
   
   assert( nintervals==db->getArraySize("b") );
   vector<double> b(nintervals);
   db->getDoubleArray("b",&b[0],nintervals);
   
   assert( nintervals==db->getArraySize("c") );
   vector<double> c(nintervals);
   db->getDoubleArray("c",&c[0],nintervals);
   
   assert( nintervals==db->getArraySize("d2") );
   vector<double> d2(nintervals);
   db->getDoubleArray("d2",&d2[0],nintervals);

   vector<double> d3(nintervals);   
   if ( db->keyExists( "d3" ) ) {
      db->getDoubleArray("d3",&d3[0],nintervals);
   }else{
      for(unsigned i=0;i<nintervals;i++)d3[i]=0.;
   }
   
   vector<double> d4(nintervals);
   if ( db->keyExists( "d4" ) ) {
      db->getDoubleArray("d4",&d4[0],nintervals);
   }else{
      for(unsigned i=0;i<nintervals;i++)d4[i]=0.;
   }

   vector<double> d7(nintervals); 
   if ( db->keyExists( "d7" ) ) {
      db->getDoubleArray("d7",&d7[0],nintervals);
   }else{
      for(unsigned i=0;i<nintervals;i++)d7[i]=0.;
   }

   vector<double> dm1(nintervals);
   if ( db->keyExists( "dm1" ) ) {
      db->getDoubleArray("dm1",&dm1[0],nintervals);
   }else{
      for(unsigned i=0;i<nintervals;i++)dm1[i]=0.;
   }

   vector<double> dm9(nintervals); 
   if ( db->keyExists( "dm9" ) ) {
      db->getDoubleArray("dm9",&dm9[0],nintervals);
   }else{
      for(unsigned i=0;i<nintervals;i++)dm9[i]=0.;
   }
   
   if ( db->keyExists( "d5" ) ) {
      TBOX_ERROR( "CALPHADSpeciesPhaseGibbsEnergy: T**5 not implemented!!!"
                  << endl );
   }
   if ( db->keyExists( "d6" ) ) {
      TBOX_ERROR( "CALPHADSpeciesPhaseGibbsEnergy: T**6 not implemented!!!"
                  << endl );
   }
   if ( db->keyExists( "dm2" ) ) {
      TBOX_ERROR( "CALPHADSpeciesPhaseGibbsEnergy: T**-2 not implemented!!!"
                  << endl );
   }

   for(unsigned i=0;i<nintervals;i++){
      CALPHADSpeciesPhaseGibbsEnergyExpansion expan(a[i],b[i],c[i],
         d2[i],d3[i],d4[i],d7[i],dm1[i],dm9[i]);
      d_expansion.push_back(expan);
   }
}

/////////////////////////////////////////////////////////////////////
// Free energy function
// parameters are in J/mol
// returned values are in J/mol
double CALPHADSpeciesPhaseGibbsEnergy::fenergy(
   const double T) // expect T in Kelvin
{
   const int n=(int)d_tc.size();
   //tbox::pout<<"n="<<n<<endl;
   assert( n>1 );
   
   for(int i=0;i<n-1;i++)
   if( T>=d_tc[i] && T<d_tc[i+1] ){
      return d_expansion[i].value(T);
   }
   
   cerr<<"T="<<T<<", Tmin="<<d_tc[0]
             <<", Tmax="<<d_tc[n-1]<<endl;
   TBOX_ERROR( "T out of range for fenergy" << endl );

   return 0.;   
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
