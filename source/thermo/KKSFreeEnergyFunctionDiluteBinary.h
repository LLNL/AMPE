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
#ifndef included_KKSFreeEnergyFunctionDiluteBinary
#define included_KKSFreeEnergyFunctionDiluteBinary

#include "Phases.h"
#include "KKSdiluteBinaryConcentrationSolver.h"
#include "FreeEnergyFunctions.h"
#include "InterpolationType.h"

#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"

using namespace SAMRAI;

class KKSFreeEnergyFunctionDiluteBinary : public FreeEnergyFunctions
{
 public:
   KKSFreeEnergyFunctionDiluteBinary(
       std::shared_ptr<SAMRAI::tbox::Database> conc_db,
       const EnergyInterpolationType energy_interp_func_type,
       const ConcInterpolationType conc_interp_func_type);

   ~KKSFreeEnergyFunctionDiluteBinary() { delete d_solver; };

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

   void preRunDiagnostics(const double T0 = 300., const double T1 = 3000.) {}

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

 private:
   KKSdiluteBinaryConcentrationSolver* d_solver;

   double d_ceq_l;
   double d_ceq_a;

   EnergyInterpolationType d_energy_interp_func_type;
   ConcInterpolationType d_conc_interp_func_type;

   void readNewtonparameters(std::shared_ptr<tbox::Database> newton_db);

   void setupFB(const double temperature);

   std::string d_fenergy_diag_filename;

   double d_fA;
   double d_fB;

   double d_Tm;
   double d_me;
   double d_ke;

   void readParameters(std::shared_ptr<SAMRAI::tbox::Database> conc_db);

   void setupSolver(std::shared_ptr<tbox::Database> newton_db);

   void computePhasesFreeEnergies(const double temperature, const double hphi,
                                  const double conc, double& fl, double& fa);
};

#endif
