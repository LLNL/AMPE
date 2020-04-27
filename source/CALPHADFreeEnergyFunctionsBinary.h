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
#ifndef included_CALPHADFreeEnergyFunctionsBinary
#define included_CALPHADFreeEnergyFunctionsBinary

#include "Phases.h"
#include "CALPHADSpeciesPhaseGibbsEnergy.h"
#include "CALPHADConcSolverBinary.h"
#include "CALPHADEqConcSolverBinary.h"
#include "CALPHADFreeEnergyFunctions.h"
#include "InterpolationType.h"

#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"

#ifdef HAVE_TLOGT
#define MAX_POL_T_INDEX 3
#else
#define MAX_POL_T_INDEX 2
#endif

class CALPHADFreeEnergyFunctionsBinary : public CALPHADFreeEnergyFunctions
{
 public:
   CALPHADFreeEnergyFunctionsBinary(
       std::shared_ptr<SAMRAI::tbox::Database> input_db,
       std::shared_ptr<SAMRAI::tbox::Database> newton_db,
       const EnergyInterpolationType energy_interp_func_type,
       const ConcInterpolationType conc_interp_func_type,
       const bool with_third_phase);

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

   virtual double computeFreeEnergy(const double temperature,
                                    const double* const conc,
                                    const PhaseIndex pi, const bool gp = false);
   virtual void computeDerivFreeEnergy(const double temperature,
                                       const double* const conc,
                                       const PhaseIndex pi, double*);
   virtual void computeSecondDerivativeFreeEnergy(const double temp,
                                                  const double* const conc,
                                                  const PhaseIndex pi,
                                                  std::vector<double>& d2fdc2);

   virtual bool computeCeqT(const double temperature, const PhaseIndex pi0,
                            const PhaseIndex pi1, double* ceq,
                            const int maxits = 20, const bool verbose = false);

   void preRunDiagnostics(const double T0 = 300., const double T1 = 3000.)
   {
      std::ofstream os1("FlC0vsT.dat", std::ios::out);
      os1 << "#Species 0, Phase L" << std::endl;
      d_g_species_phaseL[0].plotFofT(os1, T0, T1);

      std::ofstream os2("FlC1vsT.dat", std::ios::out);
      os2 << "#Species 1, Phase L" << std::endl;
      d_g_species_phaseL[1].plotFofT(os2, T0, T1);

      std::ofstream os3("FsC0vsT.dat", std::ios::out);
      os3 << "#Species 0, Phase A" << std::endl;
      d_g_species_phaseA[0].plotFofT(os3, T0, T1);

      std::ofstream os4("FsC1vsT.dat", std::ios::out);
      os4 << "#Species 1, Phase A" << std::endl;
      d_g_species_phaseA[1].plotFofT(os4, T0, T1);
   }

   int computePhaseConcentrations(const double temperature, const double* conc,
                                  const double phi, const double eta,
                                  double* x);
   void energyVsPhiAndC(
       const double temperature, const double* const ceq, const bool found_ceq,
       const double phi_well_scale, const std::string& phi_well_type,
       const int npts_phi = 51,
       const int npts_c = 50);  // # of compositions to use (>1)
   void printEnergyVsComposition(const double temperature,
                                 const int npts = 100);
   double fchem(const double phi, const double eta, const double* const conc,
                const double temperature);
   void printEnergyVsPhiHeader(const double temperature, const int nphi,
                               const int nc, const double cmin,
                               const double cmax, const double slopec,
                               std::ostream& os) const;
   void printEnergyVsPhi(const double* const conc, const double temperature,
                         const double phi_well_scale,
                         const std::string& phi_well_type, const int npts,
                         const double slopec, std::ostream& os);
   void printEnergyVsEta(const double* const conc, const double temperature,
                         const double eta_well_scale,
                         const std::string& eta_well_type, const int npts,
                         const double slopec, std::ostream& os);

   // empty default implementation to avoid downcasting
   virtual double computePenalty(const PhaseIndex, const double) { return 0.; };
   virtual double computeDerivPenalty(const PhaseIndex, const double)
   {
      return 0.;
   };
   virtual double compute2ndDerivPenalty(const PhaseIndex, const double)
   {
      return 0.;
   };

 protected:
   CALPHADConcentrationSolverBinary* d_solver;

   double d_ceq_l;
   double d_ceq_a;
   double d_ceq_b;

   EnergyInterpolationType d_energy_interp_func_type;
   ConcInterpolationType d_conc_interp_func_type;

   bool d_with_third_phase;

   void readNewtonparameters(std::shared_ptr<tbox::Database> newton_db);

   void setupValuesForTwoPhasesSolver(const double temperature, double* L0,
                                      double* L1, double* L2, double* L3,
                                      double* fA, double* fB,
                                      const PhaseIndex pi0,
                                      const PhaseIndex pi1);

   void setup(const double temperature);

 private:
   std::string d_fenergy_diag_filename;

   // size 2 for species 0 and 1
   CALPHADSpeciesPhaseGibbsEnergy d_g_species_phaseL[2];
   CALPHADSpeciesPhaseGibbsEnergy d_g_species_phaseA[2];
   CALPHADSpeciesPhaseGibbsEnergy d_g_species_phaseB[2];

   // size 4 for L0, L1, L2, L3,
   // can contain up to 3 coefficients a,b,c for a+b*T,
   // possibly +c*T*ln(T) if compiled with -DHAVE_TLOGT
   double d_LmixPhaseL[4][MAX_POL_T_INDEX];
   double d_LmixPhaseA[4][MAX_POL_T_INDEX];
   double d_LmixPhaseB[4][MAX_POL_T_INDEX];

   double* d_fA;
   double* d_fB;
   /*
    * L values evaluated at temperature T (index corresponding to phase)
    */
   double* d_L0;
   double* d_L1;
   double* d_L2;
   double* d_L3;

   void readParameters(std::shared_ptr<SAMRAI::tbox::Database> calphad_db);

   void setupSolver(std::shared_ptr<tbox::Database> newton_db);

   // energy of species "is" in phase L,A,B
   double getFenergyPhaseL(const short is, const double temperature)
   {
      return d_g_species_phaseL[is].fenergy(temperature);
   }
   double getFenergyPhaseA(const short is, const double temperature)
   {
      return d_g_species_phaseA[is].fenergy(temperature);
   }
   double getFenergyPhaseB(const short is, const double temperature)
   {
      return d_g_species_phaseB[is].fenergy(temperature);
   }

   double lmixPhase(const unsigned index, const PhaseIndex pi,
                    const double temperature)
   {
      TBOX_ASSERT(index < 4);

      switch (pi) {
         case PhaseIndex::phaseL:
            return d_LmixPhaseL[index][0] + d_LmixPhaseL[index][1] * temperature
#ifdef HAVE_TLOGT
                   + d_LmixPhaseL[index][2] * temperature * log(temperature)
#endif
                ;
         case PhaseIndex::phaseA:
            return d_LmixPhaseA[index][0] + d_LmixPhaseA[index][1] * temperature
#ifdef HAVE_TLOGT
                   + d_LmixPhaseA[index][2] * temperature * log(temperature)
#endif
                ;
         case PhaseIndex::phaseB:
            return d_LmixPhaseB[index][0] + d_LmixPhaseB[index][1] * temperature
#ifdef HAVE_TLOGT
                   + d_LmixPhaseB[index][2] * temperature * log(temperature)
#endif
                ;
         default:
            SAMRAI::tbox::pout << "CALPHADFreeEnergyStrategy::lmix0Phase(), "
                                  "undefined phase"
                               << "!!!" << std::endl;
            SAMRAI::tbox::SAMRAI_MPI::abort();
            return 0.;
      }
   }

   void computePhasesFreeEnergies(const double temperature, const double hphi,
                                  const double heta, const double conc,
                                  double& fl, double& fa, double& fb);
};

#endif
