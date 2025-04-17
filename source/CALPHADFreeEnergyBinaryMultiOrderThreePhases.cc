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
#include "CALPHADFreeEnergyBinaryMultiOrderThreePhases.h"
#include "Database2JSON.h"

#include "CALPHADFreeEnergyFunctionsBinaryThreePhase.h"
#include "CALPHADFreeEnergyFunctionsBinary2Ph1Sl.h"

#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"

using namespace SAMRAI;
namespace pt = boost::property_tree;

#include <cassert>

//=======================================================================

template <class FreeEnergyFunctionType>
CALPHADFreeEnergyBinaryMultiOrderThreePhases<FreeEnergyFunctionType>::
    CALPHADFreeEnergyBinaryMultiOrderThreePhases(
        boost::property_tree::ptree calphad_pt,
        std::shared_ptr<tbox::Database> newton_db,
        const Thermo4PFM::ConcInterpolationType conc_interp_func_type,
        const short norderp_A, MolarVolumeStrategy* mvstrategy,
        const int conc_l_id, const int conc_a_id, const int conc_b_id)
    : FreeEnergyStrategyBinary(Thermo4PFM::EnergyInterpolationType::LINEAR,
                               conc_interp_func_type, conc_l_id, conc_a_id,
                               conc_b_id, false),
      d_mv_strategy(mvstrategy)
{
   setup(calphad_pt, newton_db);

   d_multiorder_driving_force.reset(
       new MultiOrderBinaryThreePhasesDrivingForce(this, norderp_A));
}

template <class FreeEnergyFunctionType>
void CALPHADFreeEnergyBinaryMultiOrderThreePhases<
    FreeEnergyFunctionType>::setup(pt::ptree calphad_pt,
                                   std::shared_ptr<tbox::Database> newton_db)
{
   pt::ptree newton_pt;
   if (newton_db) copyDatabase(newton_db, newton_pt);
   d_calphad_fenergy.reset(new FreeEnergyFunctionType(calphad_pt, newton_pt,
                                                      d_energy_interp_func_type,
                                                      d_conc_interp_func_type));
}

//=======================================================================

template <class FreeEnergyFunctionType>
bool CALPHADFreeEnergyBinaryMultiOrderThreePhases<
    FreeEnergyFunctionType>::computeCeqT(const double temperature,
                                         const Thermo4PFM::PhaseIndex pi0,
                                         const Thermo4PFM::PhaseIndex pi1,
                                         double* ceq)
{
   TBOX_ERROR(
       "CALPHADFreeEnergyBinaryMultiOrderThreePhases::computeCeqT() not "
       "implemented");
   return false;
   // return d_calphad_fenergy->computeCeqT(temperature, &ceq[0], 50, true);
}

//=======================================================================

template <class FreeEnergyFunctionType>
double CALPHADFreeEnergyBinaryMultiOrderThreePhases<
    FreeEnergyFunctionType>::computeFreeEnergy(const double temperature,
                                               double* c_i,
                                               const Thermo4PFM::PhaseIndex pi,
                                               const bool gp)
{
   double f = d_calphad_fenergy->computeFreeEnergy(temperature, c_i, pi, gp);
   return f * d_mv_strategy->computeInvMolarVolume(temperature, c_i, pi);
}

//=======================================================================

template <class FreeEnergyFunctionType>
double CALPHADFreeEnergyBinaryMultiOrderThreePhases<FreeEnergyFunctionType>::
    computeDerivFreeEnergy(const double temperature, double* c_i,
                           const Thermo4PFM::PhaseIndex pi)
{
   double deriv;
   d_calphad_fenergy->computeDerivFreeEnergy(temperature, c_i, pi, &deriv);
   return deriv * d_mv_strategy->computeInvMolarVolume(temperature, c_i, pi);
}

//=======================================================================

template <class FreeEnergyFunctionType>
void CALPHADFreeEnergyBinaryMultiOrderThreePhases<
    FreeEnergyFunctionType>::addDrivingForce(const double time,
                                             hier::Patch& patch,
                                             const int temperature_id,
                                             const int phase_id,
                                             const int conc_id,
                                             const int f_l_id, const int f_a_id,
                                             const int f_b_id, const int rhs_id)
{
   (void)time;

   d_multiorder_driving_force->addDrivingForce(patch, temperature_id, phase_id,
                                               d_conc_l_id, d_conc_a_id,
                                               d_conc_b_id, f_l_id, f_a_id,
                                               f_b_id, rhs_id);
}

//=======================================================================

template <class FreeEnergyFunctionType>
double CALPHADFreeEnergyBinaryMultiOrderThreePhases<
    FreeEnergyFunctionType>::computeMuL(const double t, const double c)
{
   double mu;
   d_calphad_fenergy->computeDerivFreeEnergy(t, &c,
                                             Thermo4PFM::PhaseIndex::phaseL,
                                             &mu);
   return mu *
          d_mv_strategy->computeInvMolarVolume(t, &c,
                                               Thermo4PFM::PhaseIndex::phaseL);
}

//=======================================================================

template <class FreeEnergyFunctionType>
double CALPHADFreeEnergyBinaryMultiOrderThreePhases<
    FreeEnergyFunctionType>::computeMuA(const double t, const double c)
{
   double mu;
   d_calphad_fenergy->computeDerivFreeEnergy(t, &c,
                                             Thermo4PFM::PhaseIndex::phaseA,
                                             &mu);
   return mu *
          d_mv_strategy->computeInvMolarVolume(t, &c,
                                               Thermo4PFM::PhaseIndex::phaseA);
}

template <class FreeEnergyFunctionType>
void CALPHADFreeEnergyBinaryMultiOrderThreePhases<FreeEnergyFunctionType>::
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

template <class FreeEnergyFunctionType>
void CALPHADFreeEnergyBinaryMultiOrderThreePhases<FreeEnergyFunctionType>::
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

template class CALPHADFreeEnergyBinaryMultiOrderThreePhases<
    Thermo4PFM::CALPHADFreeEnergyFunctionsBinaryThreePhase>;
template class CALPHADFreeEnergyBinaryMultiOrderThreePhases<
    Thermo4PFM::CALPHADFreeEnergyFunctionsBinary2Ph1Sl>;
