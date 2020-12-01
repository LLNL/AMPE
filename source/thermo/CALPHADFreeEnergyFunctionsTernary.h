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
#ifndef included_CALPHADFreeEnergyFunctionsTernary
#define included_CALPHADFreeEnergyFunctionsTernary

#include "Phases.h"
#include "CALPHADSpeciesPhaseGibbsEnergy.h"
#include "CALPHADConcSolverTernary.h"
#include "CALPHADEqConcSolverTernary.h"
#include "CALPHADEqPhaseConcSolverTernary.h"
#include "CALPHADFreeEnergyFunctions.h"
#include "InterpolationType.h"

#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"

class CALPHADFreeEnergyFunctionsTernary : public CALPHADFreeEnergyFunctions
{
 public:
   CALPHADFreeEnergyFunctionsTernary(
       std::shared_ptr<SAMRAI::tbox::Database> input_db,
       std::shared_ptr<SAMRAI::tbox::Database> newton_db,
       const EnergyInterpolationType energy_interp_func_type,
       const ConcInterpolationType conc_interp_func_type);

   ~CALPHADFreeEnergyFunctionsTernary() { delete d_solver; };

   virtual double computeFreeEnergy(const double temperature,
                                    const double* const conc,
                                    const PhaseIndex pi, const bool gp = false);
   virtual void computeDerivFreeEnergy(const double temperature,
                                       const double* const conc,
                                       const PhaseIndex pi, double* deriv);
   virtual void computeSecondDerivativeFreeEnergy(const double temp,
                                                  const double* const conc,
                                                  const PhaseIndex pi,
                                                  std::vector<double>& d2fdc2);

   virtual bool computeCeqT(const double temperature, const PhaseIndex pi0,
                            const PhaseIndex pi1, double* ceq,
                            const int maxits = 20, const bool verbose = false);

   virtual bool computeCeqT(const double temperature, const PhaseIndex pi0,
                            const PhaseIndex pi1, const double c0,
                            const double c1, double* ceq, const int maxits = 20,
                            const bool verbose = false);

   void preRunDiagnostics(const double T0 = 300., const double T1 = 3000.)
   {
      std::ofstream os1("FlC0vsT.dat", std::ios::out);
      os1 << "#Species 0, Phase L" << std::endl;
      d_g_species_phaseL[0].plotFofT(os1, T0, T1);

      std::ofstream os2("FlC1vsT.dat", std::ios::out);
      os2 << "#Species 1, Phase L" << std::endl;
      d_g_species_phaseL[1].plotFofT(os2, T0, T1);

      std::ofstream os3("FlC2vsT.dat", std::ios::out);
      os3 << "#Species 2, Phase L" << std::endl;
      d_g_species_phaseL[2].plotFofT(os3, T0, T1);

      std::ofstream os4("FsC0vsT.dat", std::ios::out);
      os4 << "#Species 0, Phase A" << std::endl;
      d_g_species_phaseA[0].plotFofT(os4, T0, T1);

      std::ofstream os5("FsC1vsT.dat", std::ios::out);
      os5 << "#Species 1, Phase A" << std::endl;
      d_g_species_phaseA[1].plotFofT(os5, T0, T1);

      std::ofstream os6("FsC2vsT.dat", std::ios::out);
      os6 << "#Species 2, Phase A" << std::endl;
      d_g_species_phaseA[2].plotFofT(os6, T0, T1);
   }

   int computePhaseConcentrations(const double temperature,
                                  const double* const conc, const double phi,
                                  const double eta, double* x);
   void energyVsPhiAndC(
       const double temperature, const double* const ceq, const bool found_ceq,
       const double phi_well_scale, const std::string& phi_well_type,
       const int npts_phi = 51,
       const int npts_c = 50);  // number of compositions to use (>1)
   void printEnergyVsComposition(const double temperature,
                                 const int npts = 100);
   double fchem(const double phi, const double eta, const double* const conc,
                const double temperature);
   void printEnergyVsPhiHeader(const double temperature, const int nphi,
                               const int nc0, const int nc1, const double c0min,
                               const double c0max, const double c1min,
                               const double c1max, std::ostream& os) const;
   void printEnergyVsPhi(const double* const conc, const double temperature,
                         const double phi_well_scale,
                         const std::string& phi_well_type, const int npts,
                         std::ostream& os);

   // empty default implementation to avoid downcasting
   virtual double computePenalty(const PhaseIndex index, const double conc)
   {
      (void)index;
      (void)conc;
      return 0.;
   };
   virtual double computeDerivPenalty(const PhaseIndex index, const double conc)
   {
      (void)index;
      (void)conc;
      return 0.;
   };
   virtual double compute2ndDerivPenalty(const PhaseIndex index,
                                         const double conc)
   {
      (void)index;
      (void)conc;
      return 0.;
   };

 protected:
   CALPHADConcentrationSolverTernary* d_solver;

   double d_ceq_l[2];
   double d_ceq_s[2];

   EnergyInterpolationType d_energy_interp_func_type;
   ConcInterpolationType d_conc_interp_func_type;

   void readNewtonparameters(std::shared_ptr<tbox::Database> newton_db);

   void setupValuesForTwoPhasesSolver(const double temperature,
                                      const PhaseIndex pi0,
                                      const PhaseIndex pi1);

   void setup(const double temperature);

 private:
   std::string d_fenergy_diag_filename;

   // size 3 for species A, B, C
   CALPHADSpeciesPhaseGibbsEnergy d_g_species_phaseL[3];
   CALPHADSpeciesPhaseGibbsEnergy d_g_species_phaseA[3];

   // size 4 for L0, L1, L2, L3, with 2 coefficient for linear expansion in T
   double d_LmixABPhaseL[4][2];
   double d_LmixABPhaseA[4][2];

   double d_LmixACPhaseL[4][2];
   double d_LmixACPhaseA[4][2];

   double d_LmixBCPhaseL[4][2];
   double d_LmixBCPhaseA[4][2];

   double d_LmixABCPhaseL[3][2];
   double d_LmixABCPhaseA[3][2];

   double d_L_AB_L[4];
   double d_L_AC_L[4];
   double d_L_BC_L[4];

   double d_L_ABC_L[3];

   double d_L_AB_S[4];
   double d_L_AC_S[4];
   double d_L_BC_S[4];

   double d_L_ABC_S[3];

   double d_fA[2];
   double d_fB[2];
   double d_fC[2];

   void readParameters(std::shared_ptr<SAMRAI::tbox::Database> calphad_db);

   void setupSolver(std::shared_ptr<tbox::Database> newton_db);

   void setupValuesL(const double temperature);
   void setupValuesS(const double temperature);

   // energy of species "is" in phase L,A,B
   double getFenergyPhaseL(const short is, const double temperature)
   {
      return d_g_species_phaseL[is].fenergy(temperature);
   }
   double getFenergyPhaseA(const short is, const double temperature)
   {
      return d_g_species_phaseA[is].fenergy(temperature);
   }

   double lmix0ABPhase(const PhaseIndex pi, const double temperature)
   {
      switch (pi) {
         case PhaseIndex::phaseL: return lmix0ABPhaseL(temperature);
         case PhaseIndex::phaseA: return lmix0ABPhaseA(temperature);
         default:
            SAMRAI::tbox::pout << "CALPHADFreeEnergyStrategy::lmix0ABPhase(), "
                                  "undefined phase!!!"
                               << std::endl;
            SAMRAI::tbox::SAMRAI_MPI::abort();
            return 0.;
      }
   }

   double lmix1ABPhase(const PhaseIndex pi, const double temperature)
   {
      switch (pi) {
         case PhaseIndex::phaseL: return lmix1ABPhaseL(temperature);
         case PhaseIndex::phaseA: return lmix1ABPhaseA(temperature);
         default:
            SAMRAI::tbox::pout << "CALPHADFreeEnergyStrategy::lmix1ABPhase(), "
                                  "undefined phase!!!"
                               << std::endl;
            SAMRAI::tbox::SAMRAI_MPI::abort();
            return 0.;
      }
   }

   double lmix2ABPhase(const PhaseIndex pi, const double temperature)
   {
      switch (pi) {
         case PhaseIndex::phaseL: return lmix2ABPhaseL(temperature);
         case PhaseIndex::phaseA: return lmix2ABPhaseA(temperature);
         default:
            SAMRAI::tbox::pout << "CALPHADFreeEnergyStrategy::lmix2ABPhase(), "
                                  "undefined phase!!!"
                               << std::endl;
            SAMRAI::tbox::SAMRAI_MPI::abort();
            return 0.;
      }
   }

   double lmix3ABPhase(const PhaseIndex pi, const double temperature)
   {
      switch (pi) {
         case PhaseIndex::phaseL: return lmix3ABPhaseL(temperature);
         case PhaseIndex::phaseA: return lmix3ABPhaseA(temperature);
         default:
            SAMRAI::tbox::pout << "CALPHADFreeEnergyStrategy::lmix3ABPhase(), "
                                  "undefined phase!!!"
                               << std::endl;
            SAMRAI::tbox::SAMRAI_MPI::abort();
            return 0.;
      }
   }

   double lmix0ABPhaseL(const double temperature)
   {
      return d_LmixABPhaseL[0][0] + d_LmixABPhaseL[0][1] * temperature;
   }

   double lmix1ABPhaseL(const double temperature)
   {
      return d_LmixABPhaseL[1][0] + d_LmixABPhaseL[1][1] * temperature;
   }

   double lmix2ABPhaseL(const double temperature)
   {
      return d_LmixABPhaseL[2][0] + d_LmixABPhaseL[2][1] * temperature;
   }

   double lmix3ABPhaseL(const double temperature)
   {
      return d_LmixABPhaseL[3][0] + d_LmixABPhaseL[3][1] * temperature;
   }

   double lmix0ABPhaseA(const double temperature)
   {
      return d_LmixABPhaseA[0][0] + d_LmixABPhaseA[0][1] * temperature;
   }

   double lmix1ABPhaseA(const double temperature)
   {
      return d_LmixABPhaseA[1][0] + d_LmixABPhaseA[1][1] * temperature;
   }

   double lmix2ABPhaseA(const double temperature)
   {
      return d_LmixABPhaseA[2][0] + d_LmixABPhaseA[2][1] * temperature;
   }

   double lmix3ABPhaseA(const double temperature)
   {
      return d_LmixABPhaseA[3][0] + d_LmixABPhaseA[3][1] * temperature;
   }

   double lmix0ACPhase(const PhaseIndex pi, const double temperature)
   {
      switch (pi) {
         case PhaseIndex::phaseL: return lmix0ACPhaseL(temperature);
         case PhaseIndex::phaseA: return lmix0ACPhaseA(temperature);
         default:
            SAMRAI::tbox::pout << "CALPHADFreeEnergyStrategy::lmix0ACPhase(), "
                                  "undefined phase!!!"
                               << std::endl;
            SAMRAI::tbox::SAMRAI_MPI::abort();
            return 0.;
      }
   }

   double lmix1ACPhase(const PhaseIndex pi, const double temperature)
   {
      switch (pi) {
         case PhaseIndex::phaseL: return lmix1ACPhaseL(temperature);
         case PhaseIndex::phaseA: return lmix1ACPhaseA(temperature);
         default:
            SAMRAI::tbox::pout << "CALPHADFreeEnergyStrategy::lmix1ACPhase(), "
                                  "undefined phase!!!"
                               << std::endl;
            SAMRAI::tbox::SAMRAI_MPI::abort();
            return 0.;
      }
   }

   double lmix2ACPhase(const PhaseIndex pi, const double temperature)
   {
      switch (pi) {
         case PhaseIndex::phaseL: return lmix2ACPhaseL(temperature);
         case PhaseIndex::phaseA: return lmix2ACPhaseA(temperature);
         default:
            SAMRAI::tbox::pout << "CALPHADFreeEnergyStrategy::lmix2ACPhase(), "
                                  "undefined phase!!!"
                               << std::endl;
            SAMRAI::tbox::SAMRAI_MPI::abort();
            return 0.;
      }
   }

   double lmix3ACPhase(const PhaseIndex pi, const double temperature)
   {
      switch (pi) {
         case PhaseIndex::phaseL: return lmix3ACPhaseL(temperature);
         case PhaseIndex::phaseA: return lmix3ACPhaseA(temperature);
         default:
            SAMRAI::tbox::pout << "CALPHADFreeEnergyStrategy::lmix3ACPhase(), "
                                  "undefined phase!!!"
                               << std::endl;
            SAMRAI::tbox::SAMRAI_MPI::abort();
            return 0.;
      }
   }

   double lmix0ACPhaseL(const double temperature)
   {
      return d_LmixACPhaseL[0][0] + d_LmixACPhaseL[0][1] * temperature;
   }

   double lmix1ACPhaseL(const double temperature)
   {
      return d_LmixACPhaseL[1][0] + d_LmixACPhaseL[1][1] * temperature;
   }

   double lmix2ACPhaseL(const double temperature)
   {
      return d_LmixACPhaseL[2][0] + d_LmixACPhaseL[2][1] * temperature;
   }

   double lmix3ACPhaseL(const double temperature)
   {
      return d_LmixACPhaseL[3][0] + d_LmixACPhaseL[3][1] * temperature;
   }

   double lmix0ACPhaseA(const double temperature)
   {
      return d_LmixACPhaseA[0][0] + d_LmixACPhaseA[0][1] * temperature;
   }

   double lmix1ACPhaseA(const double temperature)
   {
      return d_LmixACPhaseA[1][0] + d_LmixACPhaseA[1][1] * temperature;
   }

   double lmix2ACPhaseA(const double temperature)
   {
      return d_LmixACPhaseA[2][0] + d_LmixACPhaseA[2][1] * temperature;
   }

   double lmix3ACPhaseA(const double temperature)
   {
      return d_LmixACPhaseA[3][0] + d_LmixACPhaseA[3][1] * temperature;
   }

   double lmix0BCPhase(const PhaseIndex pi, const double temperature)
   {
      switch (pi) {
         case PhaseIndex::phaseL: return lmix0BCPhaseL(temperature);
         case PhaseIndex::phaseA: return lmix0BCPhaseA(temperature);
         default:
            SAMRAI::tbox::pout << "CALPHADFreeEnergyStrategy::lmix0BCPhase(), "
                                  "undefined phase!!!"
                               << std::endl;
            SAMRAI::tbox::SAMRAI_MPI::abort();
            return 0.;
      }
   }

   double lmix1BCPhase(const PhaseIndex pi, const double temperature)
   {
      switch (pi) {
         case PhaseIndex::phaseL: return lmix1BCPhaseL(temperature);
         case PhaseIndex::phaseA: return lmix1BCPhaseA(temperature);
         default:
            SAMRAI::tbox::pout << "CALPHADFreeEnergyStrategy::lmix1BCPhase(), "
                                  "undefined phase!!!"
                               << std::endl;
            SAMRAI::tbox::SAMRAI_MPI::abort();
            return 0.;
      }
   }

   double lmix2BCPhase(const PhaseIndex pi, const double temperature)
   {
      switch (pi) {
         case PhaseIndex::phaseL: return lmix2BCPhaseL(temperature);
         case PhaseIndex::phaseA: return lmix2BCPhaseA(temperature);
         default:
            SAMRAI::tbox::pout << "CALPHADFreeEnergyStrategy::lmix2BCPhase(), "
                                  "undefined phase!!!"
                               << std::endl;
            SAMRAI::tbox::SAMRAI_MPI::abort();
            return 0.;
      }
   }

   double lmix3BCPhase(const PhaseIndex pi, const double temperature)
   {
      switch (pi) {
         case PhaseIndex::phaseL: return lmix3BCPhaseL(temperature);
         case PhaseIndex::phaseA: return lmix3BCPhaseA(temperature);
         default:
            SAMRAI::tbox::pout << "CALPHADFreeEnergyStrategy::lmix3BCPhase(), "
                                  "undefined phase!!!"
                               << std::endl;
            SAMRAI::tbox::SAMRAI_MPI::abort();
            return 0.;
      }
   }

   double lmix0BCPhaseL(const double temperature)
   {
      return d_LmixBCPhaseL[0][0] + d_LmixBCPhaseL[0][1] * temperature;
   }

   double lmix1BCPhaseL(const double temperature)
   {
      return d_LmixBCPhaseL[1][0] + d_LmixBCPhaseL[1][1] * temperature;
   }

   double lmix2BCPhaseL(const double temperature)
   {
      return d_LmixBCPhaseL[2][0] + d_LmixBCPhaseL[2][1] * temperature;
   }

   double lmix3BCPhaseL(const double temperature)
   {
      return d_LmixBCPhaseL[3][0] + d_LmixBCPhaseL[3][1] * temperature;
   }

   double lmix0BCPhaseA(const double temperature)
   {
      return d_LmixBCPhaseA[0][0] + d_LmixBCPhaseA[0][1] * temperature;
   }

   double lmix1BCPhaseA(const double temperature)
   {
      return d_LmixBCPhaseA[1][0] + d_LmixBCPhaseA[1][1] * temperature;
   }

   double lmix2BCPhaseA(const double temperature)
   {
      return d_LmixBCPhaseA[2][0] + d_LmixBCPhaseA[2][1] * temperature;
   }

   double lmix3BCPhaseA(const double temperature)
   {
      return d_LmixBCPhaseA[3][0] + d_LmixBCPhaseA[3][1] * temperature;
   }

   // ABC
   double lmix0ABCPhase(const PhaseIndex pi, const double temperature)
   {
      switch (pi) {
         case PhaseIndex::phaseL: return lmix0ABCPhaseL(temperature);
         case PhaseIndex::phaseA: return lmix0ABCPhaseA(temperature);
         default:
            SAMRAI::tbox::pout << "CALPHADFreeEnergyStrategy::lmix0ABCPhase(), "
                                  "undefined phase!!!"
                               << std::endl;
            SAMRAI::tbox::SAMRAI_MPI::abort();
            return 0.;
      }
   }

   double lmix1ABCPhase(const PhaseIndex pi, const double temperature)
   {
      switch (pi) {
         case PhaseIndex::phaseL: return lmix1ABCPhaseL(temperature);
         case PhaseIndex::phaseA: return lmix1ABCPhaseA(temperature);
         default:
            SAMRAI::tbox::pout << "CALPHADFreeEnergyStrategy::lmix1ABCPhase(), "
                                  "undefined phase!!!"
                               << std::endl;
            SAMRAI::tbox::SAMRAI_MPI::abort();
            return 0.;
      }
   }

   double lmix2ABCPhase(const PhaseIndex pi, const double temperature)
   {
      switch (pi) {
         case PhaseIndex::phaseL: return lmix2ABCPhaseL(temperature);
         case PhaseIndex::phaseA: return lmix2ABCPhaseA(temperature);
         default:
            SAMRAI::tbox::pout << "CALPHADFreeEnergyStrategy::lmix2ABCPhase(), "
                                  "undefined phase!!!"
                               << std::endl;
            SAMRAI::tbox::SAMRAI_MPI::abort();
            return 0.;
      }
   }

   // ABC liquid
   double lmix0ABCPhaseL(const double temperature)
   {
      assert(d_LmixABCPhaseL[0][0] == d_LmixABCPhaseL[0][0]);

      return d_LmixABCPhaseL[0][0] + d_LmixABCPhaseL[0][1] * temperature;
   }

   double lmix1ABCPhaseL(const double temperature)
   {
      return d_LmixABCPhaseL[1][0] + d_LmixABCPhaseL[1][1] * temperature;
   }

   double lmix2ABCPhaseL(const double temperature)
   {
      return d_LmixABCPhaseL[2][0] + d_LmixABCPhaseL[2][1] * temperature;
   }

   // ABC solid
   double lmix0ABCPhaseA(const double temperature)
   {
      return d_LmixABCPhaseA[0][0] + d_LmixABCPhaseA[0][1] * temperature;
   }

   double lmix1ABCPhaseA(const double temperature)
   {
      return d_LmixABCPhaseA[1][0] + d_LmixABCPhaseA[1][1] * temperature;
   }

   double lmix2ABCPhaseA(const double temperature)
   {
      return d_LmixABCPhaseA[2][0] + d_LmixABCPhaseA[2][1] * temperature;
   }

   void computePhasesFreeEnergies(const double temperature, const double hphi,
                                  const double conc0, const double conc1,
                                  double& fl, double& fa);
};

#endif
