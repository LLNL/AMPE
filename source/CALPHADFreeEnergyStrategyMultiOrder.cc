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
#include "CALPHADFreeEnergyStrategyMultiOrder.h"
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
CALPHADFreeEnergyStrategyMultiOrder<FreeEnergyFunctionType>::
    CALPHADFreeEnergyStrategyMultiOrder(
        pt::ptree calphad_pt, std::shared_ptr<tbox::Database> newton_db,
        const Thermo4PFM::ConcInterpolationType conc_interp_func_type,
        const short norderp_A, MolarVolumeStrategy* mvstrategy,
        const int conc_l_id, const int conc_a_id, const int conc_b_id)
    : CALPHADFreeEnergyStrategyBinary<FreeEnergyFunctionType>(
          calphad_pt, newton_db, Thermo4PFM::EnergyInterpolationType::LINEAR,
          conc_interp_func_type, mvstrategy, conc_l_id, conc_a_id, conc_b_id,
          false),
      d_norderp_A(norderp_A)
{
}

//=======================================================================

template <class FreeEnergyFunctionType>
void CALPHADFreeEnergyStrategyMultiOrder<FreeEnergyFunctionType>::
    addDrivingForce(std::shared_ptr<pdat::CellData<double> > cd_rhs,
                    std::shared_ptr<pdat::CellData<double> > cd_temperature,
                    std::shared_ptr<pdat::CellData<double> > cd_phi,
                    std::shared_ptr<pdat::CellData<double> > cd_eta,
                    std::shared_ptr<pdat::CellData<double> > cd_f_l,
                    std::shared_ptr<pdat::CellData<double> > cd_f_a,
                    std::shared_ptr<pdat::CellData<double> > cd_f_b,
                    std::shared_ptr<pdat::CellData<double> > cd_c_l,
                    std::shared_ptr<pdat::CellData<double> > cd_c_a,
                    std::shared_ptr<pdat::CellData<double> > cd_c_b,
                    const hier::Box& pbox)
{
   (void)cd_eta;

   assert(cd_f_l->getGhostCellWidth()[0] == cd_f_a->getGhostCellWidth()[0]);
   assert(cd_c_l->getGhostCellWidth()[0] == cd_c_a->getGhostCellWidth()[0]);
   assert(cd_phi->getDepth() > 1);
   assert(cd_phi->getDepth() == cd_rhs->getDepth());
   const int norderp = cd_phi->getDepth();
   if (cd_c_b == nullptr) assert(d_norderp_A == norderp - 1);

   const bool with_phaseB = (cd_c_b != nullptr);
   // if(with_phaseB)std::cout<<"With phase B!"<<std::endl;
   // else std::cout<<"Without phase B!"<<std::endl;
   double* ptr_temp = cd_temperature->getPointer();
   double* ptr_f_l = cd_f_l->getPointer();
   double* ptr_f_a = cd_f_a->getPointer();
   double* ptr_f_b = with_phaseB ? cd_f_b->getPointer() : nullptr;

   double* ptr_c_l = cd_c_l->getPointer();
   double* ptr_c_a = cd_c_a->getPointer();
   double* ptr_c_b = with_phaseB ? cd_c_b->getPointer() : nullptr;

   const hier::Box& rhs_gbox = cd_rhs->getGhostBox();
   int imin_rhs = rhs_gbox.lower(0);
   int jmin_rhs = rhs_gbox.lower(1);
   int jp_rhs = rhs_gbox.numberCells(0);
   int kmin_rhs = 0;
   int kp_rhs = 0;
#if (NDIM == 3)
   kmin_rhs = rhs_gbox.lower(2);
   kp_rhs = jp_rhs * rhs_gbox.numberCells(1);
#endif

   const hier::Box& temp_gbox = cd_temperature->getGhostBox();
   int imin_temp = temp_gbox.lower(0);
   int jmin_temp = temp_gbox.lower(1);
   int jp_temp = temp_gbox.numberCells(0);
   int kmin_temp = 0;
   int kp_temp = 0;
#if (NDIM == 3)
   kmin_temp = temp_gbox.lower(2);
   kp_temp = jp_temp * temp_gbox.numberCells(1);
#endif

   const hier::Box& pf_gbox = cd_phi->getGhostBox();
   int imin_pf = pf_gbox.lower(0);
   int jmin_pf = pf_gbox.lower(1);
   int jp_pf = pf_gbox.numberCells(0);
   int kmin_pf = 0;
   int kp_pf = 0;
#if (NDIM == 3)
   kmin_pf = pf_gbox.lower(2);
   kp_pf = jp_pf * pf_gbox.numberCells(1);
#endif

   // Assuming f_l, f_a, all have same ghost box
   const hier::Box& f_i_gbox = cd_f_l->getGhostBox();
   int imin_f_i = f_i_gbox.lower(0);
   int jmin_f_i = f_i_gbox.lower(1);
   int jp_f_i = f_i_gbox.numberCells(0);
   int kmin_f_i = 0;
   int kp_f_i = 0;
#if (NDIM == 3)
   kmin_f_i = f_i_gbox.lower(2);
   kp_f_i = jp_f_i * f_i_gbox.numberCells(1);
#endif

   // Assuming c_l, c_a, all have same ghost box
   const hier::Box& c_i_gbox = cd_c_l->getGhostBox();
   int imin_c_i = c_i_gbox.lower(0);
   int jmin_c_i = c_i_gbox.lower(1);
   int jp_c_i = c_i_gbox.numberCells(0);
   int kmin_c_i = 0;
   int kp_c_i = 0;
#if (NDIM == 3)
   kmin_c_i = c_i_gbox.lower(2);
   kp_c_i = jp_c_i * c_i_gbox.numberCells(1);
#endif

   int imin = pbox.lower(0);
   int imax = pbox.upper(0);
   int jmin = pbox.lower(1);
   int jmax = pbox.upper(1);
   int kmin = 0;
   int kmax = 0;
#if (NDIM == 3)
   kmin = pbox.lower(2);
   kmax = pbox.upper(2);
#endif

   std::vector<double> rhs(norderp);

   std::vector<double*> ptr_rhs(norderp);
   for (short i = 0; i < norderp; i++)
      ptr_rhs[i] = cd_rhs->getPointer(i);

   std::vector<double*> ptr_phi(norderp);
   for (short i = 0; i < norderp; i++)
      ptr_phi[i] = cd_phi->getPointer(i);

   for (int kk = kmin; kk <= kmax; kk++) {
      for (int jj = jmin; jj <= jmax; jj++) {
         for (int ii = imin; ii <= imax; ii++) {

            const int idx_rhs = (ii - imin_rhs) + (jj - jmin_rhs) * jp_rhs +
                                (kk - kmin_rhs) * kp_rhs;

            const int idx_temp = (ii - imin_temp) + (jj - jmin_temp) * jp_temp +
                                 (kk - kmin_temp) * kp_temp;

            const int idx_pf = (ii - imin_pf) + (jj - jmin_pf) * jp_pf +
                               (kk - kmin_pf) * kp_pf;

            const int idx_f_i = (ii - imin_f_i) + (jj - jmin_f_i) * jp_f_i +
                                (kk - kmin_f_i) * kp_f_i;

            const int idx_c_i = (ii - imin_c_i) + (jj - jmin_c_i) * jp_c_i +
                                (kk - kmin_c_i) * kp_c_i;

            double t = ptr_temp[idx_temp];
            double fl = ptr_f_l[idx_f_i];
            double fa = ptr_f_a[idx_f_i];
            double cl = ptr_c_l[idx_c_i];
            double ca = ptr_c_a[idx_c_i];
            double mu = this->computeMuL(t, cl);

            //
            // see Moelans, Acta Mat 59 (2011)
            //

            // driving forces
            double dfa = (fa - mu * ca);
            double dfl = (fl - mu * cl);
            assert(!std::isnan(dfa));

            double dfb = 0.;
            if (with_phaseB) {
               double cb = ptr_c_b[idx_c_i];
               double fb = ptr_f_b[idx_f_i];
               dfb = (fb - mu * cb);
               // std::cout<<"dfb="<<dfb<<std::endl;
            }

            // interpolation polynomials
            double hphiA = 0.;
            for (short i = 0; i < d_norderp_A; i++)
               hphiA += ptr_phi[i][idx_pf] * ptr_phi[i][idx_pf];
            assert(!std::isnan(hphiA));

            double hphiB = 0.;
            for (short i = d_norderp_A; i < norderp - 1; i++)
               hphiB += ptr_phi[i][idx_pf] * ptr_phi[i][idx_pf];
            assert(!std::isnan(hphiB));

            double hphil =
                ptr_phi[norderp - 1][idx_pf] * ptr_phi[norderp - 1][idx_pf];

            const double sum2 = hphil + hphiA + hphiB;
            assert(sum2 > 0.);
            const double sum2inv = 1. / sum2;

            hphil *= sum2inv;
            hphiA *= sum2inv;
            hphiB *= sum2inv;

            assert(!std::isnan(hphiA));
            assert(!std::isnan(hphiB));

            // solid phase A order parameters
            for (short i = 0; i < d_norderp_A; i++)
               rhs[i] = 2. * ptr_phi[i][idx_pf] *
                        ((1. - hphiA) * dfa - hphil * dfl - hphiB * dfb) *
                        sum2inv;
            // solid phase B order parameters
            for (short i = d_norderp_A; i < norderp - 1; i++)
               rhs[i] = 2. * ptr_phi[i][idx_pf] *
                        ((1. - hphiB) * dfb - hphil * dfl - hphiA * dfa) *
                        sum2inv;
            // liquid phase order parameter
            rhs[norderp - 1] =
                2. * ptr_phi[norderp - 1][idx_pf] *
                ((1. - hphil) * dfl - hphiA * dfa - hphiB * dfb) * sum2inv;
            for (short i = 0; i < norderp; i++)
               assert(!std::isnan(rhs[i]));

            for (short i = 0; i < norderp; i++)
               ptr_rhs[i][idx_rhs] -= (rhs[i]);
         }
      }
   }
}

template class CALPHADFreeEnergyStrategyMultiOrder<
    Thermo4PFM::CALPHADFreeEnergyFunctionsBinaryThreePhase>;
template class CALPHADFreeEnergyStrategyMultiOrder<
    Thermo4PFM::CALPHADFreeEnergyFunctionsBinary2Ph1Sl>;
