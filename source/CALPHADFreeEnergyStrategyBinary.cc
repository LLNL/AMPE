// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.
//
#include "ConcFort.h"
#include "CALPHADFreeEnergyStrategyBinary.h"
#include "CALPHADFunctions.h"
#include "MolarVolumeStrategy.h"
#include "CALPHADFreeEnergyFunctionsBinary.h"
#include "CALPHADFreeEnergyFunctionsBinaryThreePhase.h"
#include "FuncFort.h"

#include "CALPHADFreeEnergyFunctionsBinary2Ph1Sl.h"
#include "Database2JSON.h"

namespace pt = boost::property_tree;

#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"

using namespace SAMRAI;


#include <cassert>

//=======================================================================

template <class FreeEnergyFunctionType>
CALPHADFreeEnergyStrategyBinary<FreeEnergyFunctionType>::
    CALPHADFreeEnergyStrategyBinary(
        pt::ptree calphad_db, std::shared_ptr<tbox::Database> newton_db,
        const Thermo4PFM::EnergyInterpolationType energy_interp_func_type,
        const Thermo4PFM::ConcInterpolationType conc_interp_func_type,
        MolarVolumeStrategy* mvstrategy, const int conc_l_id,
        const int conc_a_id, const int conc_b_id, const bool with_third_phase)
    : FreeEnergyStrategyBinary(energy_interp_func_type, conc_interp_func_type,
                               conc_l_id, conc_a_id, conc_b_id,
                               with_third_phase),
      d_mv_strategy(mvstrategy)
{
   // conversion factor from [J/mol] to [pJ/(mu m)^3]
   // vm^-1 [mol/m^3] * 10e-18 [m^3/(mu m^3)] * 10e12 [pJ/J]
   // d_jpmol2pjpmumcube = 1.e-6 / d_vm;

   // R = 8.314472 J · K-1 · mol-1
   // tbox::plog << "CALPHADFreeEnergyStrategyBinary:" << std::endl;
   // tbox::plog << "Molar volume L =" << vml << std::endl;
   // tbox::plog << "Molar volume A =" << vma << std::endl;
   // tbox::plog << "jpmol2pjpmumcube=" << d_jpmol2pjpmumcube << std::endl;

   setup(calphad_db, newton_db);
}

//=======================================================================

template <class FreeEnergyFunctionType>
void CALPHADFreeEnergyStrategyBinary<FreeEnergyFunctionType>::setup(
    pt::ptree calphad_pt, std::shared_ptr<tbox::Database> newton_db)
{
   pt::ptree newton_pt;
   if (newton_db) copyDatabase(newton_db, newton_pt);
   d_calphad_fenergy.reset(new FreeEnergyFunctionType(calphad_pt, newton_pt,
                                                      d_energy_interp_func_type,
                                                      d_conc_interp_func_type));
}

//=======================================================================

template <class FreeEnergyFunctionType>
bool CALPHADFreeEnergyStrategyBinary<FreeEnergyFunctionType>::computeCeqT(
    const double temperature, const Thermo4PFM::PhaseIndex pi0,
    const Thermo4PFM::PhaseIndex pi1, double* ceq)
{
   TBOX_ERROR("computeCeqT not implented for that class");
   return false;
}

template <>
bool CALPHADFreeEnergyStrategyBinary<
    Thermo4PFM::CALPHADFreeEnergyFunctionsBinary>::
    computeCeqT(const double temperature, const Thermo4PFM::PhaseIndex pi0,
                const Thermo4PFM::PhaseIndex pi1, double* ceq)
{
   return d_calphad_fenergy->computeCeqT(temperature, &ceq[0], 50, true);
}

//=======================================================================

template <class FreeEnergyFunctionType>
double CALPHADFreeEnergyStrategyBinary<
    FreeEnergyFunctionType>::computeFreeEnergy(const double temperature,
                                               double* c_i,
                                               const Thermo4PFM::PhaseIndex pi,
                                               const bool gp)
{
   double f = d_calphad_fenergy->computeFreeEnergy(temperature, c_i, pi, gp);
   f *= d_mv_strategy->computeInvMolarVolume(temperature, c_i, pi);
   return f;
}

//=======================================================================

template <class FreeEnergyFunctionType>
double CALPHADFreeEnergyStrategyBinary<FreeEnergyFunctionType>::
    computeDerivFreeEnergy(const double temperature, double* c_i,
                           const Thermo4PFM::PhaseIndex pi)
{
   double deriv;
   d_calphad_fenergy->computeDerivFreeEnergy(temperature, c_i, pi, &deriv);
   deriv *= d_mv_strategy->computeInvMolarVolume(temperature, c_i, pi);
   return deriv;
}

//=======================================================================

template <class FreeEnergyFunctionType>
double CALPHADFreeEnergyStrategyBinary<FreeEnergyFunctionType>::computeMuA(
    const double t, const double c)
{
   double mu;
   d_calphad_fenergy->computeDerivFreeEnergy(t, &c,
                                             Thermo4PFM::PhaseIndex::phaseA,
                                             &mu);
   mu *= d_mv_strategy->computeInvMolarVolume(t, &c,
                                              Thermo4PFM::PhaseIndex::phaseA);

   return mu;
}

//=======================================================================

template <class FreeEnergyFunctionType>
double CALPHADFreeEnergyStrategyBinary<FreeEnergyFunctionType>::computeMuL(
    const double t, const double c)
{
   double mu;
   d_calphad_fenergy->computeDerivFreeEnergy(t, &c,
                                             Thermo4PFM::PhaseIndex::phaseL,
                                             &mu);
   mu *= d_mv_strategy->computeInvMolarVolume(t, &c,
                                              Thermo4PFM::PhaseIndex::phaseL);

   return mu;
}

template <class FreeEnergyFunctionType>
double CALPHADFreeEnergyStrategyBinary<FreeEnergyFunctionType>::computeMuB(
    const double t, const double c)
{
   double mu;
   d_calphad_fenergy->computeDerivFreeEnergy(t, &c,
                                             Thermo4PFM::PhaseIndex::phaseB,
                                             &mu);
   mu *= d_mv_strategy->computeInvMolarVolume(t, &c,
                                              Thermo4PFM::PhaseIndex::phaseB);

   return mu;
}


//=======================================================================

template <class FreeEnergyFunctionType>
void CALPHADFreeEnergyStrategyBinary<FreeEnergyFunctionType>::
    computeSecondDerivativeEnergyPhaseL(const double temp,
                                        const std::vector<double>& c_l,
                                        std::vector<double>& d2fdc2,
                                        const bool use_internal_units)
{
   d_calphad_fenergy->computeSecondDerivativeFreeEnergy(
       temp, &c_l[0], Thermo4PFM::PhaseIndex::phaseL, d2fdc2.data());

   if (use_internal_units)
      d2fdc2[0] *=
          d_mv_strategy->computeInvMolarVolume(temp, &c_l[0],
                                               Thermo4PFM::PhaseIndex::phaseL);
}

//=======================================================================

template <class FreeEnergyFunctionType>
void CALPHADFreeEnergyStrategyBinary<FreeEnergyFunctionType>::
    computeSecondDerivativeEnergyPhaseA(const double temp,
                                        const std::vector<double>& c_a,
                                        std::vector<double>& d2fdc2,
                                        const bool use_internal_units)
{
   d_calphad_fenergy->computeSecondDerivativeFreeEnergy(
       temp, &c_a[0], Thermo4PFM::PhaseIndex::phaseA, d2fdc2.data());

   if (use_internal_units)
      d2fdc2[0] *=
          d_mv_strategy->computeInvMolarVolume(temp, &c_a[0],
                                               Thermo4PFM::PhaseIndex::phaseA);
}


template class CALPHADFreeEnergyStrategyBinary<
    Thermo4PFM::CALPHADFreeEnergyFunctionsBinary>;
template class CALPHADFreeEnergyStrategyBinary<
    Thermo4PFM::CALPHADFreeEnergyFunctionsBinary2Ph1Sl>;
template class CALPHADFreeEnergyStrategyBinary<
    Thermo4PFM::CALPHADFreeEnergyFunctionsBinaryThreePhase>;
