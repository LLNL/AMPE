// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
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
// LLC, UT BATTELLE, LLC, 
// THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
// OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
// 
#include "CALPHADFunctions.h"
#include "CALPHADFreeEnergyFunctionsBinary.h"
#include "PhysicalConstants.h"
#include "CALPHADEqConcSolverBinary.h"
#include "FuncFort.h"


#include <string>
using namespace std;

#ifdef HAVE_TLOT
void readLcoefficients(boost::shared_ptr<tbox::Database> db,
                       double (&LmixPhase)[4][3]){
#else
void readLcoefficients(boost::shared_ptr<tbox::Database> db,
                       double (&LmixPhase)[4][2]){
#endif
   size_t lsize = db->getArraySize("L0");
   db->getDoubleArray( "L0", LmixPhase[0], lsize );
#ifdef HAVE_TLOGT
   if(lsize<3)LmixPhase[0][2]=0.;
#endif

   lsize = db->getArraySize("L1");
   db->getDoubleArray( "L1", LmixPhase[1], lsize );
#ifdef HAVE_TLOGT
   if(lsize<3)LmixPhase[1][2]=0.;
#endif

   if ( db->keyExists( "L2" ) ) {
      lsize = db->getArraySize("L2");
      db->getDoubleArray( "L2", LmixPhase[2], lsize );
#ifdef HAVE_TLOGT
      if(lsize<3)LmixPhase[2][2]=0.;
#endif
   }else {
      LmixPhase[2][0] = 0.0;
      LmixPhase[2][1] = 0.0;
#ifdef HAVE_TLOGT
      LmixPhase[2][2] = 0.0;
#endif
   }

   if ( db->keyExists( "L3" ) ) {
      lsize = db->getArraySize("L3");
      db->getDoubleArray( "L3", LmixPhase[3], lsize );
#ifdef HAVE_TLOGT
      if(lsize<3)LmixPhase[3][2]=0.;
#endif
   }else {
      LmixPhase[3][0] = 0.0;
      LmixPhase[3][1] = 0.0;
#ifdef HAVE_TLOGT
      LmixPhase[3][2] = 0.0;
#endif
   }
}

CALPHADFreeEnergyFunctionsBinary::CALPHADFreeEnergyFunctionsBinary(
   boost::shared_ptr<SAMRAI::tbox::Database> calphad_db,
   boost::shared_ptr<SAMRAI::tbox::Database> newton_db,
   const std::string& energy_interp_func_type,
   const std::string& conc_interp_func_type,
   const bool with_third_phase):
      d_energy_interp_func_type(energy_interp_func_type),
      d_conc_interp_func_type(conc_interp_func_type),
      d_with_third_phase(with_third_phase)
{
   d_fenergy_diag_filename = "energy.vtk";

   int N = 2;
   if ( d_with_third_phase ) {
      N = 3;
   }

   d_fA = new double[N];
   d_fB = new double[N];
   d_L0 = new double[N];
   d_L1 = new double[N];
   d_L2 = new double[N];
   d_L3 = new double[N];

   d_ceq_l=-1;
   d_ceq_a=-1;
   d_ceq_b=-1;
   
   readParameters(calphad_db);

   setupSolver(newton_db);
}

//=======================================================================

void CALPHADFreeEnergyFunctionsBinary::setupSolver(
   boost::shared_ptr<tbox::Database> newton_db)
{
   tbox::plog << "CALPHADFreeEnergyFunctionsBinary::setupSolver()..." << endl;
   d_solver = new CALPHADConcentrationSolverBinary( d_with_third_phase );

   readNewtonparameters(newton_db);
}

//=======================================================================

void CALPHADFreeEnergyFunctionsBinary::readNewtonparameters(boost::shared_ptr<tbox::Database> newton_db)
{
   if( newton_db!=NULL ){
      double tol =newton_db->getDoubleWithDefault( "tol", 1.e-8 );
      double alpha = newton_db->getDoubleWithDefault( "alpha", 1. );
      int maxits = newton_db->getIntegerWithDefault( "max_its", 20 );
      const bool verbose = newton_db->getBoolWithDefault( "verbose", false );
   
      d_solver->SetTolerance( tol );
      d_solver->SetMaxIterations( maxits );
      d_solver->SetDamping( alpha );
      d_solver->SetVerbose( verbose );
   }
}

//=======================================================================

void CALPHADFreeEnergyFunctionsBinary::readParameters(
   boost::shared_ptr<tbox::Database> calphad_db)
{   
   boost::shared_ptr<tbox::Database> species0_db =
      calphad_db->getDatabase( "SpeciesA" );
   string name = species0_db->getStringWithDefault( "name", "unknown" );
   string dbnameL("PhaseL");
   if( species0_db->keyExists( "Phase0" ) )
   {
      tbox::pout << "Input Phase0 is deprecated.  Use PhaseL."<< endl;
      dbnameL="Phase0";
   }
   d_g_species_phaseL[0].initialize( name, species0_db->getDatabase(dbnameL) );
   string dbnameA("PhaseA");
   if( species0_db->keyExists( "Phase1" ) )
   {
      tbox::pout << "Input Phase1 is deprecated.  Use PhaseA."<< endl;
      dbnameA="Phase1";
   }
   d_g_species_phaseA[0].initialize( name, species0_db->getDatabase(dbnameA) );
   string dbnameB("PhaseB");
   if( species0_db->keyExists( "Phase2" ) )
   {
      tbox::pout << "Input Phase2 is deprecated.  Use PhaseB."<< endl;
      dbnameB="PhaseB";
   }
   if ( d_with_third_phase ) {
      d_g_species_phaseB[0].initialize( name,
                                        species0_db->getDatabase(dbnameB) );
   }      

   boost::shared_ptr<tbox::Database> speciesB_db =
      calphad_db->getDatabase( "SpeciesB" );
   name = speciesB_db->getStringWithDefault( "name", "unknown" );
   d_g_species_phaseL[1].initialize( name, speciesB_db->getDatabase(dbnameL) );
   d_g_species_phaseA[1].initialize( name, speciesB_db->getDatabase(dbnameA) );
   if ( d_with_third_phase ) {
      d_g_species_phaseB[1].initialize( name,
                                        speciesB_db->getDatabase(dbnameB) );
   }      

   // read Lmix coefficients
   string dbnamemixL("LmixPhaseL");
   if( calphad_db->keyExists( "LmixPhase0" ) )
   {
      tbox::pout << "Input LmixPhase0 is deprecated.  Use LmixPhaseL."<< endl;
      dbnamemixL="LmixPhase0";
   }
   string dbnamemixA("LmixPhaseA");
   if( calphad_db->keyExists( "LmixPhase1" ) )
   {
      tbox::pout << "Input LmixPhase1 is deprecated.  Use LmixPhaseA."<< endl;
      dbnamemixA="LmixPhase1";
   }
   string dbnamemixB("LmixPhaseB");
   if( calphad_db->keyExists( "LmixPhase2" ) )
   {
      tbox::pout << "Input LmixPhase2 is deprecated.  Use LmixPhaseB."<< endl;
      dbnamemixB="LmixPhase2";
   }
   boost::shared_ptr<tbox::Database> Lmix0_db =
      calphad_db->getDatabase( dbnamemixL );
   readLcoefficients(Lmix0_db, d_LmixPhaseL);

   boost::shared_ptr<tbox::Database> Lmix1_db =
      calphad_db->getDatabase(dbnamemixA);
   readLcoefficients(Lmix1_db, d_LmixPhaseA);

   if ( d_with_third_phase ) {
      boost::shared_ptr<tbox::Database> Lmix2_db =
         calphad_db->getDatabase(dbnamemixB);
      readLcoefficients(Lmix2_db, d_LmixPhaseB);
   }

   // print database just read
   tbox::plog << "CALPHAD database..." << endl;
   calphad_db->printClassData( tbox::plog );
}

//-----------------------------------------------------------------------

double CALPHADFreeEnergyFunctionsBinary::computeFreeEnergy(
   const double temperature,
   const double* const conc,
   const PHASE_INDEX pi,
   const bool gp )
{
   const double l0 = lmixPhase( 0, pi, temperature );
   const double l1 = lmixPhase( 1, pi, temperature );
   const double l2 = lmixPhase( 2, pi, temperature );
   const double l3 = lmixPhase( 3, pi, temperature );
   
   CALPHADSpeciesPhaseGibbsEnergy* g_species;
   
   switch( pi ){
      case phaseL:
         g_species=&d_g_species_phaseL[0];
         break;
      case phaseA:
         g_species=&d_g_species_phaseA[0];
         break;
      case phaseB:
         g_species=&d_g_species_phaseB[0];
         break;
      default:
         SAMRAI::tbox::pout<<
            "CALPHADFreeEnergyFunctionsBinary::computeFreeEnergy(), undefined phase="
            <<pi<<"!!!"<<std::endl;
         SAMRAI::tbox::SAMRAI_MPI::abort();
      return 0.;
   }
      
   double fe =
      conc[0] * g_species[0].fenergy( temperature ) +
      ( 1. - conc[0] ) * g_species[1].fenergy( temperature ) + 
      CALPHADcomputeFMixBinary( l0, l1, l2, l3, conc[0] ) +
      CALPHADcomputeFIdealMixBinary(
         gas_constant_R_JpKpmol * temperature,
         conc[0] );

   // subtract -mu*c to get grand potential
   if( gp ){
      double deriv;
      computeDerivFreeEnergy(temperature,conc,pi,&deriv);
      fe -= deriv*conc[0];
   }
   
   return fe;
}

//=======================================================================

void CALPHADFreeEnergyFunctionsBinary::computeDerivFreeEnergy(
   const double temperature,
   const double* const conc,
   const PHASE_INDEX pi,
   double* deriv )
{
   const double l0 = lmixPhase( 0, pi, temperature );
   const double l1 = lmixPhase( 1, pi, temperature );
   const double l2 = lmixPhase( 2, pi, temperature );
   const double l3 = lmixPhase( 3, pi, temperature );
   
   CALPHADSpeciesPhaseGibbsEnergy* g_species;
   
   switch( pi ){
      case phaseL:
         g_species=&d_g_species_phaseL[0];
         break;
      case phaseA:
         g_species=&d_g_species_phaseA[0];
         break;
      case phaseB:
         g_species=&d_g_species_phaseB[0];
         break;
      default:
         SAMRAI::tbox::pout<<"CALPHADFreeEnergyFunctionsBinary::computeFreeEnergy(), undefined phase="<<pi<<"!!!"<<std::endl;
         SAMRAI::tbox::SAMRAI_MPI::abort();
      return;
   }

   double mu =
      ( g_species[0].fenergy( temperature ) - g_species[1].fenergy( temperature ) ) 
      +
      CALPHADcomputeFMix_derivBinary( l0, l1, l2, l3, conc[0] )
      +
      CALPHADcomputeFIdealMix_derivBinary(
         gas_constant_R_JpKpmol * temperature,
         conc[0] );

   deriv[0] = mu;
}

//=======================================================================

void CALPHADFreeEnergyFunctionsBinary::computeSecondDerivativeFreeEnergy(
   const double temp,
   const double* const conc,
   const PHASE_INDEX pi,
   std::vector<double>& d2fdc2)
{
   assert( conc[0]>=0. );
   assert( conc[0]<=1. );
   
   const double l0_l=lmixPhase( 0, pi, temp );
   const double l1_l=lmixPhase( 1, pi, temp );
   const double l2_l=lmixPhase( 2, pi, temp );
   const double l3_l=lmixPhase( 3, pi, temp );
   const double rt = gas_constant_R_JpKpmol * temp;

   d2fdc2[0] = 
          ( CALPHADcomputeFMix_deriv2Binary(l0_l,l1_l,l2_l,l3_l,conc[0])
           +CALPHADcomputeFIdealMix_deriv2Binary(rt,conc[0]) );   
}

//=======================================================================

void CALPHADFreeEnergyFunctionsBinary::setupValuesForTwoPhasesSolver(
   const double temperature,
   double* L0, double* L1, double* L2, double* L3,
   double* fA, double* fB,
   const PHASE_INDEX pi0, const PHASE_INDEX pi1)
{
   PHASE_INDEX pis[2]={pi0,pi1};
   
   //loop over two phases
   for(short i=0;i<2;i++)
   {
      L0[i] = lmixPhase( 0, pis[i], temperature );
      L1[i] = lmixPhase( 1, pis[i], temperature );
      L2[i] = lmixPhase( 2, pis[i], temperature );
      L3[i] = lmixPhase( 3, pis[i], temperature );

      switch ( pis[i] ) {
 
         case phaseL:
            fA[i] = d_g_species_phaseL[0].fenergy( temperature );
            fB[i] = d_g_species_phaseL[1].fenergy( temperature );
            break;
   
         case phaseA:
            fA[i] = d_g_species_phaseA[0].fenergy( temperature );
            fB[i] = d_g_species_phaseA[1].fenergy( temperature );
            break;
   
         case phaseB:
            fA[i] = d_g_species_phaseB[0].fenergy( temperature );
            fB[i] = d_g_species_phaseB[1].fenergy( temperature );
            break;
      
         default:
            cerr<<"CALPHADFreeEnergyFunctionsBinary::setupValuesForTwoPhasesSolver: Undefined phase"<<endl;
      }
   
   }
}

//=======================================================================

void CALPHADFreeEnergyFunctionsBinary::setup(const double temperature)
{
   d_fA[0] = d_g_species_phaseL[0].fenergy( temperature );
   d_fA[1] = d_g_species_phaseA[0].fenergy( temperature );

   d_fB[0] = d_g_species_phaseL[1].fenergy( temperature );
   d_fB[1] = d_g_species_phaseA[1].fenergy( temperature );

   d_L0[0] = lmixPhase( 0, phaseL, temperature );
   d_L1[0] = lmixPhase( 1, phaseL, temperature );
   d_L2[0] = lmixPhase( 2, phaseL, temperature );
   d_L3[0] = lmixPhase( 3, phaseL, temperature );

   d_L0[1] = lmixPhase( 0, phaseA, temperature );
   d_L1[1] = lmixPhase( 1, phaseA, temperature );
   d_L2[1] = lmixPhase( 2, phaseA, temperature );
   d_L3[1] = lmixPhase( 3, phaseA, temperature );

   if ( d_with_third_phase ) {
      d_fA[2] = d_g_species_phaseB[0].fenergy( temperature );
      d_fB[2] = d_g_species_phaseB[1].fenergy( temperature );

      d_L0[2] = lmixPhase( 0, phaseB, temperature );
      d_L1[2] = lmixPhase( 1, phaseB, temperature );
      d_L2[2] = lmixPhase( 2, phaseB, temperature );
      d_L3[2] = lmixPhase( 3, phaseB, temperature );
   }
}

//=======================================================================

// compute equilibrium concentrations in various phases for given temperature
bool CALPHADFreeEnergyFunctionsBinary::computeCeqT(
   const double temperature,
   const PHASE_INDEX pi0, const PHASE_INDEX pi1,
   double* ceq,
   const int maxits,
   const bool verbose )
{
   if(verbose)tbox::pout<<"CALPHADFreeEnergyFunctionsBinary::computeCeqT()"<<endl;
   assert( temperature>0. );

   setupValuesForTwoPhasesSolver(temperature, d_L0, d_L1, d_L2, d_L3,
                                 d_fA, d_fB, pi0, pi1);
   double RTinv = 1.0 / ( gas_constant_R_JpKpmol * temperature );
   CALPHADEqConcentrationSolverBinary eq_solver;
   eq_solver.SetMaxIterations(maxits);

   int ret = eq_solver.ComputeConcentration(
      ceq,
      RTinv,
      d_L0, d_L1, d_L2, d_L3,
      d_fA, d_fB );

   if( ret>=0 )
   {
      if(verbose){
         tbox::pout<<"CALPHAD, c_eq phase0="<<ceq[0]<<endl;
         tbox::pout<<"CALPHAD, c_eq phase1="<<ceq[1]<<endl;
      }
      
      if( pi1==phaseB ){
         d_ceq_b=ceq[1];
         d_ceq_l=0.5*(ceq[0]+d_ceq_l);
      }else{
         d_ceq_l=ceq[0];
         d_ceq_a=ceq[1];
      }
   }else{
      tbox::pout<<"CALPHADFreeEnergyFunctionsBinary, WARNING: ceq computation did not converge"<<endl;
   }
   
   return (ret>=0);
}

//=======================================================================

void CALPHADFreeEnergyFunctionsBinary::computePhasesFreeEnergies(
   const double temperature,
   const double hphi,
   const double heta,
   const double conc,
   double& fl,
   double& fa,
   double& fb)
{
   //tbox::pout<<"CALPHADFreeEnergyFunctionsBinary::computePhasesFreeEnergies()"<<endl;
   
   const int N = d_with_third_phase ? 3 : 2;

   double* c  = new double[N];
   for ( int ii = 0; ii < N; ii++ ) {
      c[ii] = conc;
   }
   //tbox::pout<<"d_ceq_l="<<d_ceq_l<<endl;
   //tbox::pout<<"d_ceq_a="<<d_ceq_a<<endl;
   //tbox::pout<<"d_ceq_b="<<d_ceq_b<<endl;
   if( d_ceq_l>=0. )c[0]=d_ceq_l;
   if( d_ceq_a>=0. )c[1]=d_ceq_a;
   if( d_ceq_b>=0. )c[2]=d_ceq_b;
   
   setup(temperature);

   double RTinv = 1.0 / ( gas_constant_R_JpKpmol * temperature );
   int ret = d_solver->ComputeConcentration(
      c,
      conc,
      hphi, heta,
      RTinv,
      d_L0, d_L1, d_L2, d_L3,
      d_fA, d_fB );

   if( ret<0 )
   {
      if ( d_with_third_phase ) {
         cerr<<"d_ceq_l="<<d_ceq_l<<endl;
         cerr<<"d_ceq_a="<<d_ceq_a<<endl;
         cerr<<"d_ceq_b="<<d_ceq_b<<endl;
         cerr<<"ERROR in CALPHADFreeEnergyFunctionsBinary::computePhasesFreeEnergies() ---"
             <<"conc="<<conc<<", hphi="<<hphi<<", heta="<<heta<<endl;
      }else{
         cerr<<"ERROR in CALPHADFreeEnergyFunctionsBinary::computePhasesFreeEnergies() ---"
             <<"conc="<<conc<<", hphi="<<hphi<<endl;
      }
      tbox::SAMRAI_MPI::abort();
   }

   assert( c[0]>=0. );
   fl = computeFreeEnergy(temperature,&c[0],phaseL,false);

   assert( c[1]>=0. );
   fa = computeFreeEnergy(temperature,&c[1],phaseA,false);

   fb = 0.;
   if ( d_with_third_phase ){
      assert( c[2]>=0. );
      assert( c[2]<=1. );
      fb = computeFreeEnergy(temperature,&c[2],phaseB,false);
   }
   delete[] c;
}

//-----------------------------------------------------------------------

int CALPHADFreeEnergyFunctionsBinary::computePhaseConcentrations(
   const double temperature, const double* const conc, 
   const double phi, const double eta,
   double* x)

{
   assert( x[0]>=0. );
   assert( x[1]>=0. );
   assert( x[0]<=1. );
   assert( x[1]<=1. );
   
   const double conc0 = conc[0];
   
   const double RTinv = 1.0 / ( gas_constant_R_JpKpmol * temperature );
   
   d_fA[0] = getFenergyPhaseL( 0, temperature );
   d_fA[1] = getFenergyPhaseA( 0, temperature );

   d_fB[0] = getFenergyPhaseL( 1, temperature );
   d_fB[1] = getFenergyPhaseA( 1, temperature );

   d_L0[0] = lmixPhase( 0, phaseL, temperature );
   d_L1[0] = lmixPhase( 1, phaseL, temperature );
   d_L2[0] = lmixPhase( 2, phaseL, temperature );
   d_L3[0] = lmixPhase( 3, phaseL, temperature );

   d_L0[1] = lmixPhase( 0, phaseA, temperature );
   d_L1[1] = lmixPhase( 1, phaseA, temperature );
   d_L2[1] = lmixPhase( 2, phaseA, temperature );
   d_L3[1] = lmixPhase( 3, phaseA, temperature );

   const double hphi =
      FORT_INTERP_FUNC(
         phi,
         d_conc_interp_func_type.c_str() );

   double heta = 0.0;
   //tbox::pout<<"d_ceq_a="<<d_ceq_a<<endl;
   //x[0] = ( d_ceq_l>=0. ) ? d_ceq_l : 0.5;
   //x[1] = ( d_ceq_a>=0. ) ? d_ceq_a : 0.5;

   if ( d_with_third_phase ) {
      assert( x[2]>=0. );
      assert( x[2]<=1. );
      //x[2] = ( d_ceq_b>=0. ) ? d_ceq_b : 0.5;
      
      d_fA[2] = getFenergyPhaseB( 0, temperature );
      d_fB[2] = getFenergyPhaseB( 1, temperature );

      d_L0[2] = lmixPhase( 0, phaseB, temperature );
      d_L1[2] = lmixPhase( 1, phaseB, temperature );
      d_L2[2] = lmixPhase( 2, phaseB, temperature );
      d_L3[2] = lmixPhase( 3, phaseB, temperature );

      heta =
         FORT_INTERP_FUNC(
            eta,
            d_conc_interp_func_type.c_str() );
   }
   
   // conc could be outside of [0.,1.] in a trial step
   double c0 = conc[0]>=0. ? conc[0] : 0.;
   c0 = c0<=1. ? c0 : 1.;
   int ret = d_solver->ComputeConcentration(
      x,
      c0,
      hphi, heta,
      RTinv,
      d_L0, d_L1, d_L2, d_L3,
      d_fA, d_fB );
   if( ret==-1 )
   {
      cerr<<"ERROR, CALPHADFreeEnergyFunctionsBinary::computePhaseConcentrations() failed for conc="<<conc0
          <<", hphi="<<hphi
          <<", heta="<<heta<<endl;
      sleep(5);
      tbox::SAMRAI_MPI::abort();
   }

   return ret;
}

//-----------------------------------------------------------------------

void CALPHADFreeEnergyFunctionsBinary::energyVsPhiAndC(
   const double temperature, 
   const double* const ceq,
   const bool found_ceq,
   const double phi_well_scale,
   const std::string& phi_well_type,
   const int npts_phi,
   const int npts_c)
{
   tbox::plog<<"CALPHADFreeEnergyFunctionsBinary::energyVsPhiAndC()..."<<endl;

   const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());

   double slopec = 0.;
   double fc0=0.;
   double fc1=0.;
   if( found_ceq )
   if( mpi.getRank()==0){
      // compute slope of f between equilibrium concentrations
      // to add slopec*conc to energy later on
   
      fc0 = computeFreeEnergy(temperature,&ceq[0],phaseL);
      fc1 = computeFreeEnergy(temperature,&ceq[1],phaseA);
      slopec = -(fc1-fc0)/(ceq[1]-ceq[0]);
   }
   tbox::plog<<setprecision(8)<<"fc0: "<<fc0<<"..."<<", fc1: "<<fc1<<"..."<<endl;
   tbox::plog<<"CALPHADFreeEnergyFunctionsBinary: Use slope: "<<slopec<<"..."<<endl;
   mpi.Barrier();
      
   if( mpi.getRank()==0 ){
   
      // reset cmin, cmax, deltac
      double cmin = min(ceq[0],ceq[1]);
      double cmax = max(ceq[0],ceq[1]);
      double dc=cmax-cmin;
      cmin = max( 0.25*cmin,         cmin-0.25*dc );
      cmax = min( 1.-0.25*(1.-cmax), cmax+0.25*dc );
      cmax = max( cmax, cmin+dc );
      double deltac = (cmax - cmin) / (npts_c-1);

      ofstream tfile(d_fenergy_diag_filename.data(), ios::out);
      
      printEnergyVsPhiHeader(temperature, npts_phi, 
                             npts_c, cmin, cmax, slopec,
                             tfile);

      for ( int i = 0; i < npts_c; i++ ) {
         double conc=cmin+deltac*i;
         printEnergyVsPhi( &conc, temperature, phi_well_scale, phi_well_type, 
                           npts_phi, slopec,
                           tfile );
      }
   }
   
}

// Print out free energy as a function of phase
// for given composition and temperature
// File format: ASCII VTK, readble with Visit
void CALPHADFreeEnergyFunctionsBinary::printEnergyVsPhiHeader(
   const double temperature,
   const int nphi,
   const int nc,
   const double cmin,
   const double cmax,
   const double slopec,
   std::ostream& os )const
{
   os << "# vtk DataFile Version 2.0"<<endl;
   os << "Free energy + "<<slopec<<"*c [J/mol] at T="<<temperature<<endl;
   os << "ASCII"<<endl;
   os << "DATASET STRUCTURED_POINTS"<<endl;

   os << "DIMENSIONS   "<<nphi<<" "<<nc<<" 1"<<endl;
   double asp_ratio_c = ( nc>1 ) ? (cmax-cmin)/(nc-1) : 1.;
   os << "ASPECT_RATIO "<<1./(nphi-1)<<" "<<asp_ratio_c<<" 1."
      <<endl;
   os << "ORIGIN        0. "<<cmin<<" 0."<<endl;
   os << "POINT_DATA   "<<nphi*nc<<endl;
   os << "SCALARS energy float 1"<<endl;
   os << "LOOKUP_TABLE default"<<endl;
}

//=======================================================================

void CALPHADFreeEnergyFunctionsBinary::printEnergyVsPhi(
   const double* const conc,
   const double temperature,
   const double phi_well_scale,
   const string& phi_well_type,
   const int npts,
   const double slopec,
   std::ostream& os )
{
   //tbox::pout << "CALPHADFreeEnergyFunctionsBinary::printEnergyVsPhi()..." << endl;
   const double dphi = 1.0 / (double)(npts-1);
   const double eta = 0.0;
   
   //os << "# phi     f(phi)     for c=" << conc
   //           << " eta=" << eta
   //           << " and T=" << temperature << endl;
   for ( int i = 0; i < npts; i++ ) {
      const double phi = i*dphi;

      double e = fchem( phi, eta, conc, temperature );
      const double w =
         phi_well_scale *
         FORT_WELL_FUNC( phi, phi_well_type.c_str() );

      os << e + w + slopec*conc[0] << endl;
   }
   //os << endl;
}

//=======================================================================

void CALPHADFreeEnergyFunctionsBinary::printEnergyVsEta(
   const double* const conc,
   const double temperature,
   const double eta_well_scale,
   const string& eta_well_type,
   const int npts,
   const double slopec,
   std::ostream& os )
{
   //tbox::pout << "CALPHADFreeEnergyFunctionsBinary::printEnergyVsEta()..." << endl;
   const double deta = 1.0 / (double)(npts-1);
   const double phi = 1.0;
   
   for ( int i = 0; i < npts; i++ ) {
      const double eta = i*deta;

      double e = fchem( phi, eta, conc, temperature );

      double w =
         eta_well_scale *
         FORT_WELL_FUNC( eta, eta_well_type.c_str() );

      os << e + w + slopec*conc[0] << endl;
   }
   //os << endl;
}

//=======================================================================
// compute free energy in [J/mol]
double CALPHADFreeEnergyFunctionsBinary::fchem(
   const double phi,
   const double eta,
   const double* const conc,
   const double temperature )
{
   const double hcphi =
      FORT_INTERP_FUNC( phi, d_conc_interp_func_type.c_str() );
   double heta = 0.0;
   if ( d_with_third_phase ) {
      heta = FORT_INTERP_FUNC( eta, d_conc_interp_func_type.c_str() );
   }      

   const double tol=1.e-8;
   double fl=0.;
   double fa=0.; 
   double fb=0.;
   if( (phi>tol) & (phi<(1.-tol)) )
   {
      computePhasesFreeEnergies(
         temperature, hcphi, heta, conc[0],
         fl, fa, fb);
   }else{
      if( phi<=tol )
      {
         fl=computeFreeEnergy(temperature,conc,phaseL);
      }else{
         fa=computeFreeEnergy(temperature,conc,phaseA);
         if ( d_with_third_phase )
            fb=computeFreeEnergy(temperature,conc,phaseB);
      }
   }

   const double hfphi =
      FORT_INTERP_FUNC( phi, d_energy_interp_func_type.c_str() );
   double e =
      ( 1.0 - hfphi ) * fl +
      hfphi * ( ( 1.0 - heta ) * fa + heta * fb );

   return e;
}

//=======================================================================

void CALPHADFreeEnergyFunctionsBinary::printEnergyVsComposition(
   const double temperature,
   const int npts )
{
   ofstream os("FvsC.dat", ios::out);

   const double dc = 1.0 / (double)(npts-1);
   
   os << "#phi=0" << endl;
   for ( int i = 0; i < npts; i++ ) {
      const double conc = i*dc;

      double e = fchem( 0., 0., &conc, temperature );
      os << conc <<"\t"<< e << endl;
   }
   os << endl;

   os << "#phi=1" << endl;
   for ( int i = 0; i < npts; i++ ) {
      const double conc = i*dc;

      double e = fchem( 1., 0., &conc, temperature );
      os << conc <<"\t"<< e << endl;
   }
   
   if( d_with_third_phase )
   {
      os << endl;

      os << "#eta=1" << endl;
      for ( int i = 0; i < npts; i++ ) {
         const double conc = i*dc;

         double e = fchem( 1., 1., &conc, temperature );
         os << conc <<"\t"<< e << endl;
      }
   }
}

