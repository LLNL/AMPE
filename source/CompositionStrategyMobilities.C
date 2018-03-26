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
#include "CompositionStrategyMobilities.h"
#include "QuatParams.h"
#include "ConcFort.h"
#include "QuatFort.h"
#include "CALPHADFunctions.h"

#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"

#include <cassert>

using namespace std;
const double m2toum2 = 1.e12;

CompositionStrategyMobilities::CompositionStrategyMobilities(
   boost::shared_ptr<tbox::Database> input_db,
   const int phase_scratch_id,
   const int eta_scratch_id,
   const unsigned short ncompositions,
   const int temperature_scratch_id,
   const int Mq_id,
   const vector<double>& Q_heat_transport,
   FreeEnergyStrategy* free_energy_strategy
   ):d_free_energy_strategy(free_energy_strategy)
{
   assert( ncompositions==1 );
   assert( temperature_scratch_id>=0 );
   assert( d_free_energy_strategy!=NULL );
   
   if( Mq_id>=0 ){
      assert( Q_heat_transport.size()>0 );
      tbox::plog<<"CompositionStrategyMobilities with thermal diffusion"<<endl;
   }else{
      tbox::plog<<"CompositionStrategyMobilities without thermal diffusion"<<endl;
   }
   
   d_ncompositions=ncompositions;
   
   d_phase_scratch_id = phase_scratch_id;
   d_eta_scratch_id = phase_scratch_id;
   
   d_Mq_id=Mq_id;
   
   d_temperature_scratch_id = temperature_scratch_id;

   d_with_third_phase = ( eta_scratch_id>=0 );
   
   string calphad_filename = input_db->getString( "filename" );

   boost::shared_ptr<tbox::MemoryDatabase> calphad_db(new tbox::MemoryDatabase( "calphad_db" ));
   tbox::InputManager::getManager()->parseInputFile(
      calphad_filename, calphad_db );

   boost::shared_ptr<tbox::Database> mobility_db
      = calphad_db->getDatabase( "MobilityParameters" );
   boost::shared_ptr<tbox::Database> species0_db
      = mobility_db->getDatabase( "Species0" );
   boost::shared_ptr<tbox::Database> species1_db
      = mobility_db->getDatabase( "Species1" );

   CALPHADMobility calphad_mobility0_phaseL("MobilitySpecies0");   
   calphad_mobility0_phaseL.initialize(species0_db->getDatabase( "PhaseL" ));

   CALPHADMobility calphad_mobility1_phaseL("MobilitySpecies1");
   calphad_mobility1_phaseL.initialize(species1_db->getDatabase( "PhaseL" ));

   d_calphad_mobilities_phaseL.push_back(calphad_mobility0_phaseL);
   d_calphad_mobilities_phaseL.push_back(calphad_mobility1_phaseL);

   CALPHADMobility calphad_mobility0_phaseA("MobilitySpecies0");
   calphad_mobility0_phaseA.initialize(species0_db->getDatabase( "PhaseA" ));

   CALPHADMobility calphad_mobility1_phaseA("MobilitySpecies1");
   calphad_mobility1_phaseA.initialize(species1_db->getDatabase( "PhaseA" ));

   d_calphad_mobilities_phaseA.push_back(calphad_mobility0_phaseA);
   d_calphad_mobilities_phaseA.push_back(calphad_mobility1_phaseA);

   CALPHADMobility calphad_mobility0_phaseB("MobilitySpecies0");
   CALPHADMobility calphad_mobility1_phaseB("MobilitySpecies1");
   if( eta_scratch_id>-1 ){
      calphad_mobility0_phaseB.initialize(species0_db->getDatabase( "PhaseB" ));
      calphad_mobility1_phaseB.initialize(species1_db->getDatabase( "PhaseB" ));

      d_calphad_mobilities_phaseB.push_back(calphad_mobility0_phaseB);
      d_calphad_mobilities_phaseB.push_back(calphad_mobility1_phaseB);

   }
}

//-----------------------------------------------------------------------

void CompositionStrategyMobilities::printDiagnostics(const boost::shared_ptr<hier::PatchHierarchy > hierarchy)
{
   tbox::plog<<"CompositionRHSStrategy::printDiagnostics()"<<endl;
   
   math::HierarchyCellDataOpsReal<double> cell_ops(hierarchy, 0, 0);
   const double Tmax=cell_ops.max(d_temperature_scratch_id);
   const double Tmin=cell_ops.min(d_temperature_scratch_id);
   
   assert( Tmin>0. );
   assert( Tmax>0. );
   assert( Tmin<100000. );
   assert( Tmax<100000. );
   
   int nT=0;
   double dT=0.;
   if( Tmax-Tmin > 0. )
   {
      nT=10;
      dT=(Tmax-Tmin)/nT;
   }
   
   const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
   if( mpi.getRank()==0 )
   for(int iT=0;iT<=nT;iT++)
   {
      string filename ("DvsC_T");
      double temperature = Tmin+iT*dT;
      stringstream ss (stringstream::in | stringstream::out);
      ss << temperature;
      string temps = ss.str();
      filename += temps;
      filename += "K.dat";

      ofstream tfile(&filename[0], ios::out);
      tfile<<"#Diffusion vs. composition"<<endl;
      printDiffusionVsComposition( temperature, tfile );
      //printDiffusionVsComposition( temperature, tbox::plog );
   }
}

//=======================================================================

double CompositionStrategyMobilities::computeDiffusionMobilityBinaryPhaseL(
   const double c0,
   const double temp)
{
   return computeDiffusionMobilityBinaryPhase(c0, temp, d_calphad_mobilities_phaseL);
}

double CompositionStrategyMobilities::computeDiffusionMobilityBinaryPhaseA(
   const double c0,
   const double temp)
{
   return computeDiffusionMobilityBinaryPhase(c0, temp, d_calphad_mobilities_phaseA);
}

double CompositionStrategyMobilities::computeDiffusionMobilityBinaryPhaseB(
   const double c0,
   const double temp)
{
   return computeDiffusionMobilityBinaryPhase(c0, temp, d_calphad_mobilities_phaseB);
}

//=======================================================================

void CompositionStrategyMobilities::computeDiffusionMobilityTernaryPhaseL(
   const double c0,
   const double c1,
   const double temp,
   vector<double>& mobility)
{
   computeDiffusionMobilityTernaryPhase(c0, c1, temp, d_calphad_mobilities_phaseL, mobility);
}

void CompositionStrategyMobilities::computeDiffusionMobilityTernaryPhaseA(
   const double c0,
   const double c1,
   const double temp,
   vector<double>& mobility)
{
   computeDiffusionMobilityTernaryPhase(c0, c1, temp, d_calphad_mobilities_phaseA, mobility);
}

void CompositionStrategyMobilities::computeDiffusionMobilityTernaryPhaseB(
   const double c0,
   const double c1,
   const double temp,
   vector<double>& mobility)
{
   computeDiffusionMobilityTernaryPhase(c0, c1, temp, d_calphad_mobilities_phaseB, mobility);
}

//=======================================================================

void CompositionStrategyMobilities::computeDiffusionMobilityPhaseL(
   const vector<double>& c,
   const double temp,
   vector<double>& mobility)
{
   switch( c.size() ){
      case 1:
         mobility[0]=computeDiffusionMobilityBinaryPhaseL(c[0], temp);
         break;
         
      case 2:
         computeDiffusionMobilityTernaryPhaseL(c[0], c[1], temp, mobility);
         break;

      default:
         tbox::pout<<"size of c vector="<<c.size()<<endl;
         tbox::pout<<"Error: diffusion mobility implemented for up to 3 species only!!!"<<endl;
         tbox::SAMRAI_MPI::abort();
   }
}

void CompositionStrategyMobilities::computeDiffusionMobilityPhaseA(
   const vector<double>& c,
   const double temp,
   vector<double>& mobility)
{
   switch( c.size() ){
      case 1:
         mobility[0]=computeDiffusionMobilityBinaryPhaseA(c[0], temp);
         break;
         
      case 2:
         computeDiffusionMobilityTernaryPhaseA(c[0], c[1], temp, mobility);
         break;

      default:
         tbox::pout<<"Error: diffusion mobility implemented for up to 3 species only!!!"<<endl;
         tbox::SAMRAI_MPI::abort();
   }
   
   assert( mobility[0]>=0. );
}

void CompositionStrategyMobilities::computeDiffusionMobilityPhaseB(
   const vector<double>& c,
   const double temp,
   vector<double>& mobility)
{
   switch( c.size() ){
      case 1:
         mobility[0]=computeDiffusionMobilityBinaryPhaseB(c[0], temp);
         break;
         
      case 2:
         computeDiffusionMobilityTernaryPhaseB(c[0], c[1], temp, mobility);
         break;

      default:
         tbox::pout<<"Error: diffusion mobility implemented for up to 3 species only!!!"<<endl;
         tbox::SAMRAI_MPI::abort();
   }
}


//=======================================================================

void CompositionStrategyMobilities::printMobilitiesVsComposition( const double temperature, std::ostream &os )
{
   assert( temperature>0. );
   assert( temperature<100000. );
   assert( d_free_energy_strategy != NULL );

   int nc=20;
   double eps=1.e-7;
   double dc=(1.-2*eps)/nc;
   vector<double> amob(1);
   vector<double> conc(d_ncompositions);
   const double rt=gas_constant_R_JpKpmol*temperature;

   os<<fixed;
   os<<"#Species 0: log10(atomic mobility*RT[m2/s]) vs. composition in phase L at T="<<temperature<<endl;
   for(int i=0;i<=nc;i++)
   {
      conc[0]=eps+i*dc;
      const double m0=d_calphad_mobilities_phaseL[0].getAtomicMobilityBinary(conc[0], temperature);
      assert( m0>0. );
            
      os << conc[0] <<"  "<< log10(m0*rt) << endl;
   }

   os<<"#Species 0: log10(atomic mobility*RT[m2/s]) vs. composition in phase A at T="<<temperature<<endl;
   for(int i=0;i<=nc;i++)
   {
      conc[0]=eps+i*dc;
            
      const double m0=d_calphad_mobilities_phaseA[0].getAtomicMobilityBinary(conc[0], temperature);
      assert( m0>0. );

      os << conc[0] <<"  "<< log10(m0*rt) << endl;
   }

   if( d_with_third_phase )
   {
      os<<"#Species 0: log10(atomic mobility*RT[m2/s]) vs. composition in phase B at T="<<temperature<<endl;
      for(int i=0;i<=nc;i++)
      {
         conc[0]=eps+i*dc;
            
         const double m0=d_calphad_mobilities_phaseB[0].getAtomicMobilityBinary(conc[0], temperature);
         assert( m0>0. );

         os << conc[0] <<"  "<< log10(m0*rt) << endl;
      }
   }
}

//=======================================================================

void CompositionStrategyMobilities::printDiffusionVsComposition( const double temperature, std::ostream &os )
{
   assert( temperature>0. );
   assert( temperature<100000. );
   assert( d_free_energy_strategy != NULL );

   int nc=50;
   double eps=1.e-7;
   double dc=(1.-2*eps)/nc;
   vector<double> amob(1);
   vector<double> conc(d_ncompositions);
   vector<double> d2f(d_ncompositions*d_ncompositions);
   const double rt=gas_constant_R_JpKpmol*temperature;

   os<<"#Interdiffusion[m2/s] vs. composition in phase L at T="<<temperature<<endl;
   for(int i=0;i<=nc;i++)
   {
      conc[0]=eps+i*dc;
      d_free_energy_strategy->computeSecondDerivativeEnergyPhaseL(
         temperature, conc, d2f, false);
        
      computeDiffusionMobilityPhaseL(conc, temperature, amob );

      assert( amob[0]>0. );
      os << fixed << conc[0] <<"  " <<scientific<< amob[0]*d2f[0]/m2toum2<< endl;
   }

   os << endl;
   os<<"#Interdiffusion[m2/s]  vs. composition in phase A at T="<<temperature<<endl;
   for(int i=0;i<=nc;i++)
   {
      conc[0]=eps+i*dc;
      d_free_energy_strategy->computeSecondDerivativeEnergyPhaseA(
         temperature, conc, d2f, false);
            
      computeDiffusionMobilityPhaseA(conc, temperature, amob );

      assert( amob[0]>=0. );
      os << fixed << conc[0] <<"  " <<scientific<< amob[0]*d2f[0]/m2toum2<< endl;
   }
   
   if( d_with_third_phase )
   {
      os << endl;
      os<<"#Interdiffusion[m2/s]  vs. composition in phase B at T="<<temperature<<endl;
      for(int i=0;i<=nc;i++)
      {
         conc[0]=eps+i*dc;
         d_free_energy_strategy->computeSecondDerivativeEnergyPhaseB(
            temperature, conc, d2f, false);
            
         computeDiffusionMobilityPhaseB(conc, temperature, amob );

         assert( amob[0]>0. );
         os << fixed << conc[0] <<"  " <<scientific<< amob[0]*d2f[0]/m2toum2<< endl;
      }
   }
}

