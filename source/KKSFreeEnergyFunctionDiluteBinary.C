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
#include "xlogx.h"
#include "KKSFreeEnergyFunctionDiluteBinary.h"
#include "PhysicalConstants.h"
#include "FuncFort.h"


#include <string>


KKSFreeEnergyFunctionDiluteBinary::KKSFreeEnergyFunctionDiluteBinary(
    boost::shared_ptr<SAMRAI::tbox::Database> conc_db,
    const EnergyInterpolationType energy_interp_func_type,
    const ConcInterpolationType conc_interp_func_type)
    : d_energy_interp_func_type(energy_interp_func_type),
      d_conc_interp_func_type(conc_interp_func_type)
{
   d_fenergy_diag_filename = "energy.vtk";

   d_ceq_l = -1;
   d_ceq_a = -1;

   readParameters(conc_db);

   d_fA = log(1. / d_ke);

   boost::shared_ptr<tbox::Database> newton_db;
   if (conc_db->isDatabase("NewtonSolver"))
      newton_db = conc_db->getDatabase("NewtonSolver");

   setupSolver(newton_db);
}

//=======================================================================

void KKSFreeEnergyFunctionDiluteBinary::setupSolver(
    boost::shared_ptr<tbox::Database> newton_db)
{
   tbox::plog << "KKSFreeEnergyFunctionDiluteBinary::setupSolver()..."
              << std::endl;
   d_solver = new KKSdiluteBinaryConcentrationSolver();

   readNewtonparameters(newton_db);
}

//=======================================================================

void KKSFreeEnergyFunctionDiluteBinary::readNewtonparameters(
    boost::shared_ptr<tbox::Database> newton_db)
{
   if (newton_db) {
      double tol = newton_db->getDoubleWithDefault("tol", 1.e-8);
      double alpha = newton_db->getDoubleWithDefault("alpha", 1.);
      int maxits = newton_db->getIntegerWithDefault("max_its", 20);
      const bool verbose = newton_db->getBoolWithDefault("verbose", false);

      d_solver->SetTolerance(tol);
      d_solver->SetMaxIterations(maxits);
      d_solver->SetDamping(alpha);
      d_solver->SetVerbose(verbose);
   }
}

//=======================================================================

void KKSFreeEnergyFunctionDiluteBinary::readParameters(
    boost::shared_ptr<tbox::Database> conc_db)
{
   d_me = conc_db->getDouble("liquidus_slope");
   d_Tm = conc_db->getDouble("meltingT");
   d_ke = conc_db->getDouble("keq");

   assert(d_me < 0.);
   assert(d_Tm > 0.);
   assert(d_ke <= 1.);
}

//-----------------------------------------------------------------------

double KKSFreeEnergyFunctionDiluteBinary::computeFreeEnergy(
    const double temperature, const double* const conc, const PhaseIndex pi,
    const bool gp)
{
   double fe = xlogx(conc[0]) + xlogx(1. - conc[0]);
   setupFB(temperature);

   switch (pi) {
      case PhaseIndex::phaseL: break;
      case PhaseIndex::phaseA:
         fe += conc[0] * d_fA + (1. - conc[0]) * d_fB;
         break;
      default:
         SAMRAI::tbox::pout << "KKSFreeEnergyFunctionDiluteBinary::"
                               "computeFreeEnergy(), undefined phase!!!"
                            << std::endl;
         SAMRAI::tbox::SAMRAI_MPI::abort();
         return 0.;
   }

   fe *= gas_constant_R_JpKpmol * temperature;

   // subtract -mu*c to get grand potential
   if (gp) {
      double deriv;
      computeDerivFreeEnergy(temperature, conc, pi, &deriv);
      fe -= deriv * conc[0];
   }

   return fe;
}

//=======================================================================

void KKSFreeEnergyFunctionDiluteBinary::computeDerivFreeEnergy(
    const double temperature, const double* const conc, const PhaseIndex pi,
    double* deriv)
{
   double mu = xlogx_deriv(conc[0]) - xlogx_deriv(1.0 - conc[0]);

   switch (pi) {
      case PhaseIndex::phaseL: break;
      case PhaseIndex::phaseA:
         setupFB(temperature);
         mu += (d_fA - d_fB);
         break;
      default:
         SAMRAI::tbox::pout << "KKSFreeEnergyFunctionDiluteBinary::"
                               "computeFreeEnergy(), undefined phase!!!"
                            << std::endl;
         SAMRAI::tbox::SAMRAI_MPI::abort();
         return;
   }

   deriv[0] = gas_constant_R_JpKpmol * temperature * mu;
}

//=======================================================================

void KKSFreeEnergyFunctionDiluteBinary::computeSecondDerivativeFreeEnergy(
    const double temp, const double* const conc, const PhaseIndex pi,
    std::vector<double>& d2fdc2)
{
   assert(conc[0] >= 0.);
   assert(conc[0] <= 1.);

   const double rt = gas_constant_R_JpKpmol * temp;

   d2fdc2[0] = rt * (xlogx_deriv2(conc[0]) + xlogx_deriv2(1.0 - conc[0]));
}

//=======================================================================

void KKSFreeEnergyFunctionDiluteBinary::setupFB(const double temperature)
{
   assert(d_ke > 0.);
   assert(temperature < d_Tm);
   assert(d_me < 0.);

   const double cLe = (temperature - d_Tm) / d_me;
   const double cSe = cLe * d_ke;

   assert(cLe < 1.);
   assert(cSe < 1.);

   d_fB = log(1. - cLe) - log(1. - cSe);
}

//=======================================================================

// compute equilibrium concentrations in various phases for given temperature
bool KKSFreeEnergyFunctionDiluteBinary::computeCeqT(
    const double temperature, const PhaseIndex pi0, const PhaseIndex pi1,
    double* ceq, const int maxits, const bool verbose)
{
   if (verbose)
      tbox::pout << "KKSFreeEnergyFunctionDiluteBinary::computeCeqT()"
                 << std::endl;
   assert(temperature > 0.);

   d_ceq_l = (temperature - d_Tm) / d_me;
   d_ceq_a = d_ceq_l * d_ke;

   ceq[0] = d_ceq_l;
   ceq[1] = d_ceq_a;

   return true;
}

//=======================================================================

void KKSFreeEnergyFunctionDiluteBinary::computePhasesFreeEnergies(
    const double temperature, const double hphi, const double conc, double& fl,
    double& fa)
{
   // tbox::pout<<"KKSFreeEnergyFunctionDiluteBinary::computePhasesFreeEnergies()"<<endl;

   double c[2] = {conc, conc};
   if (d_ceq_l >= 0.) c[0] = d_ceq_l;
   if (d_ceq_a >= 0.) c[1] = d_ceq_a;

   setupFB(temperature);

   double RTinv = 1.0 / (gas_constant_R_JpKpmol * temperature);
   int ret = d_solver->ComputeConcentration(c, conc, hphi, RTinv, d_fA, d_fB);

   if (ret < 0) {
      std::cerr << "ERROR in "
                   "KKSFreeEnergyFunctionDiluteBinary::"
                   "computePhasesFreeEnergies() "
                   "---"
                << "conc=" << conc << ", hphi=" << hphi << std::endl;
      tbox::SAMRAI_MPI::abort();
   }

   assert(c[0] >= 0.);
   fl = computeFreeEnergy(temperature, &c[0], PhaseIndex::phaseL, false);

   assert(c[1] >= 0.);
   fa = computeFreeEnergy(temperature, &c[1], PhaseIndex::phaseA, false);
}

//-----------------------------------------------------------------------

int KKSFreeEnergyFunctionDiluteBinary::computePhaseConcentrations(
    const double temperature, const double* const conc, const double phi,
    const double eta, double* x)

{
   assert(x[0] >= 0.);
   assert(x[1] >= 0.);
   assert(x[0] <= 1.);
   assert(x[1] <= 1.);

   const double conc0 = conc[0];

   const char interp_func_type = concInterpChar(d_conc_interp_func_type);
   const double hphi = INTERP_FUNC(phi, &interp_func_type);

   setupFB(temperature);

   // conc could be outside of [0.,1.] in a trial step
   double c0 = conc[0] >= 0. ? conc[0] : 0.;
   c0 = c0 <= 1. ? c0 : 1.;
   int ret = d_solver->ComputeConcentration(x, c0, hphi,
                                            -1.,  // unused parameter
                                            d_fA, d_fB);
   if (ret == -1) {
      std::cerr << "ERROR, "
                   "KKSFreeEnergyFunctionDiluteBinary::"
                   "computePhaseConcentrations() "
                   "failed for conc="
                << conc0 << ", hphi=" << hphi << std::endl;
      tbox::SAMRAI_MPI::abort();
   }

   return ret;
}

//-----------------------------------------------------------------------

void KKSFreeEnergyFunctionDiluteBinary::energyVsPhiAndC(
    const double temperature, const double* const ceq, const bool found_ceq,
    const double phi_well_scale, const std::string& phi_well_type,
    const int npts_phi, const int npts_c)
{
   tbox::plog << "KKSFreeEnergyFunctionDiluteBinary::energyVsPhiAndC()..."
              << std::endl;

   const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());

   double slopec = 0.;
   double fc0 = 0.;
   double fc1 = 0.;
   if (found_ceq)
      if (mpi.getRank() == 0) {
         // compute slope of f between equilibrium concentrations
         // to add slopec*conc to energy later on

         fc0 = computeFreeEnergy(temperature, &ceq[0], PhaseIndex::phaseL);
         fc1 = computeFreeEnergy(temperature, &ceq[1], PhaseIndex::phaseA);
         slopec = -(fc1 - fc0) / (ceq[1] - ceq[0]);
      }
   tbox::plog << std::setprecision(8) << "fc0: " << fc0 << "..."
              << ", fc1: " << fc1 << "..." << std::endl;
   tbox::plog << "KKSFreeEnergyFunctionDiluteBinary: Use slope: " << slopec
              << "..." << std::endl;
   mpi.Barrier();

   if (mpi.getRank() == 0) {

      // reset cmin, cmax, deltac
      double cmin = std::min(ceq[0], ceq[1]);
      double cmax = std::max(ceq[0], ceq[1]);
      double dc = cmax - cmin;
      cmin = std::max(0.25 * cmin, cmin - 0.25 * dc);
      cmax = std::min(1. - 0.25 * (1. - cmax), cmax + 0.25 * dc);
      cmax = std::max(cmax, cmin + dc);
      double deltac = (cmax - cmin) / (npts_c - 1);

      std::ofstream tfile(d_fenergy_diag_filename.data(), std::ios::out);

      printEnergyVsPhiHeader(temperature, npts_phi, npts_c, cmin, cmax, slopec,
                             tfile);

      for (int i = 0; i < npts_c; i++) {
         double conc = cmin + deltac * i;
         printEnergyVsPhi(&conc, temperature, phi_well_scale, phi_well_type,
                          npts_phi, slopec, tfile);
      }
   }
}

// Print out free energy as a function of phase
// for given composition and temperature
// File format: ASCII VTK, readble with Visit
void KKSFreeEnergyFunctionDiluteBinary::printEnergyVsPhiHeader(
    const double temperature, const int nphi, const int nc, const double cmin,
    const double cmax, const double slopec, std::ostream& os) const
{
   os << "# vtk DataFile Version 2.0" << std::endl;
   os << "Free energy + " << slopec << "*c [J/mol] at T=" << temperature
      << std::endl;
   os << "ASCII" << std::endl;
   os << "DATASET STRUCTURED_POINTS" << std::endl;

   os << "DIMENSIONS   " << nphi << " " << nc << " 1" << std::endl;
   double asp_ratio_c = (nc > 1) ? (cmax - cmin) / (nc - 1) : 1.;
   os << "ASPECT_RATIO " << 1. / (nphi - 1) << " " << asp_ratio_c << " 1."
      << std::endl;
   os << "ORIGIN        0. " << cmin << " 0." << std::endl;
   os << "POINT_DATA   " << nphi * nc << std::endl;
   os << "SCALARS energy float 1" << std::endl;
   os << "LOOKUP_TABLE default" << std::endl;
}

//=======================================================================

void KKSFreeEnergyFunctionDiluteBinary::printEnergyVsPhi(
    const double* const conc, const double temperature,
    const double phi_well_scale, const std::string& phi_well_type,
    const int npts, const double slopec, std::ostream& os)
{
   // tbox::pout << "KKSFreeEnergyFunctionDiluteBinary::printEnergyVsPhi()..."
   // << std::endl;
   const double dphi = 1.0 / (double)(npts - 1);
   const double eta = 0.0;

   for (int i = 0; i < npts; i++) {
      const double phi = i * dphi;

      double e = fchem(phi, eta, conc, temperature);
      const double w = phi_well_scale * WELL_FUNC(phi, phi_well_type.c_str());

      os << e + w + slopec * conc[0] << std::endl;
   }
}

//=======================================================================

void KKSFreeEnergyFunctionDiluteBinary::printEnergyVsEta(
    const double* const conc, const double temperature,
    const double eta_well_scale, const std::string& eta_well_type,
    const int npts, const double slopec, std::ostream& os)
{
   (void)conc;
   (void)temperature;
   (void)eta_well_scale;
   (void)eta_well_type;
   (void)npts;
   (void)slopec;
   (void)os;
}

//=======================================================================
// compute free energy in [J/mol]
double KKSFreeEnergyFunctionDiluteBinary::fchem(const double phi,
                                                const double eta,
                                                const double* const conc,
                                                const double temperature)
{
   (void)eta;

   const char interp_func_type = concInterpChar(d_conc_interp_func_type);
   const double hcphi = INTERP_FUNC(phi, &interp_func_type);

   const double tol = 1.e-8;
   double fl = 0.;
   double fa = 0.;
   if ((phi > tol) & (phi < (1. - tol))) {
      computePhasesFreeEnergies(temperature, hcphi, conc[0], fl, fa);
   } else {
      if (phi <= tol) {
         fl = computeFreeEnergy(temperature, conc, PhaseIndex::phaseL);
      } else {
         fa = computeFreeEnergy(temperature, conc, PhaseIndex::phaseA);
      }
   }

   const char interpf = energyInterpChar(d_energy_interp_func_type);
   const double hfphi = INTERP_FUNC(phi, &interpf);
   double e = (1.0 - hfphi) * fl + hfphi * fa;

   return e;
}

//=======================================================================

void KKSFreeEnergyFunctionDiluteBinary::printEnergyVsComposition(
    const double temperature, const int npts)
{
   std::ofstream os("FvsC.dat", std::ios::out);

   const double dc = 1.0 / (double)(npts - 1);

   os << "#phi=0" << std::endl;
   for (int i = 0; i < npts; i++) {
      const double conc = i * dc;

      double e = fchem(0., 0., &conc, temperature);
      os << conc << "\t" << e << std::endl;
   }
   os << std::endl << std::endl;

   os << "#phi=1" << std::endl;
   for (int i = 0; i < npts; i++) {
      const double conc = i * dc;

      double e = fchem(1., 0., &conc, temperature);
      os << conc << "\t" << e << std::endl;
   }
}
