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
#ifndef included_CALPHADFreeEnergyFunctionsBinary
#define included_CALPHADFreeEnergyFunctionsBinary 

#include "Phases.h"
#include "CALPHADSpeciesPhaseGibbsEnergy.h"
#include "CALPHADConcSolver.h"
#include "CALPHADEqConcSolver.h"
#include "FreeEnergyFunctions.h"

#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"

class CALPHADFreeEnergyFunctionsBinary:
   public FreeEnergyFunctions
{
public:
   CALPHADFreeEnergyFunctionsBinary(
      boost::shared_ptr<SAMRAI::tbox::Database> input_db,
      boost::shared_ptr<SAMRAI::tbox::Database> newton_db,
      const std::string& phase_interp_func_type,
      const std::string& eta_interp_func_type,
      const std::string& avg_func_type,
      const bool with_third_phase,
      const double  phase_well_scale,
      const double eta_well_scale,
      const std::string& phase_well_func_type,
      const std::string& eta_well_func_type );

   ~CALPHADFreeEnergyFunctionsBinary()
   {
      delete[] d_fA;
      delete[] d_fB;
      delete[] d_L0;
      delete[] d_L1;
      delete[] d_L2;
      delete[] d_L3;
      delete d_solver;
   };
   
   virtual double computeFreeEnergy(
      const double temperature,
      const double conc,
      const PHASE_INDEX pi,
      const bool gp=false  );   
   virtual double computeDerivFreeEnergy(
      const double temperature,
      const double conc,
      const PHASE_INDEX pi );
   virtual void computeSecondDerivativeFreeEnergy(
      const double temp,
      const std::vector<double>& conc,
      const PHASE_INDEX pi,
      std::vector<double>& d2fdc2);

   virtual bool computeCeqT(
      const double temperature,
      const PHASE_INDEX pi0, const PHASE_INDEX pi1,
      double* ceq,
      const bool verbose = false);

   void preRunDiagnostics(std::ostream& os)
   {
      os<<"#Species 0, Phase L"<<std::endl;
      d_g_species_phaseL[0].plotFofT(os);
      os<<"#Species 1, Phase L"<<std::endl;
      d_g_species_phaseL[1].plotFofT(os);
      os<<"#Species 0, Phase A"<<std::endl;
      d_g_species_phaseA[0].plotFofT(os);
      os<<"#Species 1, Phase A"<<std::endl;
      d_g_species_phaseA[1].plotFofT(os);
   }

   int computePhaseConcentrations(
      const double temperature, const double conc, const double phi, const double eta,
      double* x);
   void energyVsPhiAndC(const double temperature, 
                        const double* const ceq,
                        const bool found_ceq,
                        const bool third_phase,
                        const int npts_phi=51,
                        const int npts_c=50); // number of compositions to use (>1)
   void printEnergyVsComposition(const double temperature, const int npts=100 );
   double fenergy(
      const double phi,
      const double eta,
      const double conc,
      const double temperature );
   void printEnergyVsPhiHeader(
      const double temperature,
      const int nphi,
      const int nc,
      const double cmin,
      const double cmax,
      const double slopec,
      std::ostream& os )const;
   void printEnergyVsPhi(
      const double conc,
      const double temperature,
      const int npts,
      const double slopec,
      std::ostream& os );
   void printEnergyVsEta(
      const double conc,
      const double temperature,
      const int npts,
      const double slopec,
      std::ostream& os );

   // empty default implementation to avoid downcasting
   virtual double computePenalty(const PHASE_INDEX index, const double conc){return 0.;};
   virtual double computeDerivPenalty(const PHASE_INDEX index, const double conc){return 0.;};
   virtual double compute2ndDerivPenalty(const PHASE_INDEX index, const double conc){return 0.;};
   
protected:

   CALPHADConcentrationSolver* d_solver;

   double d_ceq_l;
   double d_ceq_a;
   double d_ceq_b;
   
   std::string d_phase_interp_func_type;
   std::string d_eta_interp_func_type;
   std::string d_avg_func_type;
   
   bool d_with_third_phase;

   void readNewtonparameters(boost::shared_ptr<tbox::Database> newton_db);

   void setupValuesForTwoPhasesSolver(const double temperature,
                                      double* L0, double* L1, double* L2, double* L3,
                                      double* fA, double* fB,
                                      const PHASE_INDEX pi0, const PHASE_INDEX pi1);

   void setupValuesForThreePhasesSolver(const double temperature,
                                        double* L0, double* L1, double* L2, double* L3,
                                        double* fA, double* fB);

private:

   std::string d_fenergy_diag_filename;
   
   // size 2 for species 0 and 1
   CALPHADSpeciesPhaseGibbsEnergy d_g_species_phaseL[2];
   CALPHADSpeciesPhaseGibbsEnergy d_g_species_phaseA[2];
   CALPHADSpeciesPhaseGibbsEnergy d_g_species_phaseB[2];
   
   // size 4 for L0, L1, L2, L3
   double d_LmixPhaseL[4][2];
   double d_LmixPhaseA[4][2];
   double d_LmixPhaseB[4][2];

   double* d_fA;
   double* d_fB;
   double* d_L0;
   double* d_L1;
   double* d_L2;
   double* d_L3;

   double d_phase_well_scale;
   double d_eta_well_scale;

   std::string d_phase_well_func_type;
   std::string d_eta_well_func_type;
   
   void readParameters(boost::shared_ptr<SAMRAI::tbox::Database> calphad_db);

   void setupSolver(boost::shared_ptr<tbox::Database> newton_db);

   // energy of species "is" in phase L,A,B
   double getFenergyPhaseL(const short is, const double temperature )
   {
      return d_g_species_phaseL[is].fenergy( temperature );
   }
   double getFenergyPhaseA(const short is, const double temperature )
   {
      return d_g_species_phaseA[is].fenergy( temperature );
   }
   double getFenergyPhaseB(const short is, const double temperature )
   {
      return d_g_species_phaseB[is].fenergy( temperature );
   }

   double lmix0Phase( const PHASE_INDEX pi, const double temperature )
   {
      switch( pi ){
         case phaseL:
            return lmix0PhaseL( temperature );
         case phaseA:
            return lmix0PhaseA( temperature );
         case phaseB:
            return lmix0PhaseB( temperature );
         default:
            SAMRAI::tbox::pout<<"CALPHADFreeEnergyStrategy::lmix0Phase(), undefined phase="<<pi<<"!!!"<<std::endl;
            SAMRAI::tbox::SAMRAI_MPI::abort();
         return 0.;
      }
   }
   
   double lmix1Phase( const PHASE_INDEX pi, const double temperature )
   {
      switch( pi ){
         case phaseL:
            return lmix1PhaseL( temperature );
         case phaseA:
            return lmix1PhaseA( temperature );
         case phaseB:
            return lmix1PhaseB( temperature );
         default:
            SAMRAI::tbox::pout<<"CALPHADFreeEnergyStrategy::lmix1Phase(), undefined phase="<<pi<<"!!!"<<std::endl;
            SAMRAI::tbox::SAMRAI_MPI::abort();
         return 0.;
      }
   }
   
   double lmix2Phase( const PHASE_INDEX pi, const double temperature )
   {
      switch( pi ){
         case phaseL:
            return lmix2PhaseL( temperature );
         case phaseA:
            return lmix2PhaseA( temperature );
         case phaseB:
            return lmix2PhaseB( temperature );
         default:
            SAMRAI::tbox::pout<<"CALPHADFreeEnergyStrategy::lmix2Phase(), undefined phase="<<pi<<"!!!"<<std::endl;
            SAMRAI::tbox::SAMRAI_MPI::abort();
         return 0.;
      }
   }
   
   double lmix3Phase( const PHASE_INDEX pi, const double temperature )
   {
      switch( pi ){
         case phaseL:
            return lmix3PhaseL( temperature );
         case phaseA:
            return lmix3PhaseA( temperature );
         case phaseB:
            return lmix3PhaseB( temperature );
         default:
            SAMRAI::tbox::pout<<"CALPHADFreeEnergyStrategy::lmix3Phase(), undefined phase="<<pi<<"!!!"<<std::endl;
            SAMRAI::tbox::SAMRAI_MPI::abort();
         return 0.;
      }
   }
   
   double lmix0PhaseL( const double temperature )
   {
      return d_LmixPhaseL[0][0] + d_LmixPhaseL[0][1] * temperature;
   }
   
   double lmix1PhaseL( const double temperature )
   {
      return d_LmixPhaseL[1][0] + d_LmixPhaseL[1][1] * temperature;
   }
   
   double lmix2PhaseL( const double temperature )
   {
      return d_LmixPhaseL[2][0] + d_LmixPhaseL[2][1] * temperature;
   }

   double lmix3PhaseL( const double temperature )
   {
      return d_LmixPhaseL[3][0] + d_LmixPhaseL[3][1] * temperature;
   }

   double lmix0PhaseA( const double temperature )
   {
      return d_LmixPhaseA[0][0] + d_LmixPhaseA[0][1] * temperature;
   }
   
   double lmix1PhaseA( const double temperature )
   {
      return d_LmixPhaseA[1][0] + d_LmixPhaseA[1][1] * temperature;
   }
   
   double lmix2PhaseA( const double temperature )
   {
      return d_LmixPhaseA[2][0] + d_LmixPhaseA[2][1] * temperature;
   }

   double lmix3PhaseA( const double temperature )
   {
      return d_LmixPhaseA[3][0] + d_LmixPhaseA[3][1] * temperature;
   }

   double lmix0PhaseB( const double temperature )
   {
      return d_LmixPhaseB[0][0] + d_LmixPhaseB[0][1] * temperature;
   }
   
   double lmix1PhaseB( const double temperature )
   {
      return d_LmixPhaseB[1][0] + d_LmixPhaseB[1][1] * temperature;
   }
   
   double lmix2PhaseB( const double temperature )
   {
      return d_LmixPhaseB[2][0] + d_LmixPhaseB[2][1] * temperature;
   }
   
   double lmix3PhaseB( const double temperature )
   {
      return d_LmixPhaseB[3][0] + d_LmixPhaseB[3][1] * temperature;
   }

   void computePhasesFreeEnergies(
      const double temperature,
      const double hphi,
      const double heta,
      const double conc,
      double& fl,
      double& fa,
      double& fb);
      
};

#endif
