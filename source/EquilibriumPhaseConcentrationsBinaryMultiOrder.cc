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
#include "EquilibriumPhaseConcentrationsBinaryMultiOrder.h"
#include "QuadraticFreeEnergyFunctionsBinary.h"
#include "FuncFort.h"

#include "SAMRAI/tbox/SAMRAI_MPI.h"

EquilibriumPhaseConcentrationsBinaryMultiOrder::
    EquilibriumPhaseConcentrationsBinaryMultiOrder(
        const int conc_l_id, const int conc_a_id,
        std::shared_ptr<tbox::Database> conc_db)
    : PhaseConcentrationsStrategy(conc_l_id, conc_a_id, -1)
{
}

int EquilibriumPhaseConcentrationsBinaryMultiOrder::
    computePhaseConcentrationsOnPatch(
        std::shared_ptr<pdat::CellData<double> > cd_temperature,
        std::shared_ptr<pdat::CellData<double> > cd_phi,
        std::shared_ptr<pdat::CellData<double> > cd_concentration,
        std::shared_ptr<pdat::CellData<double> > cd_c_l,
        std::shared_ptr<pdat::CellData<double> > cd_c_a,
        std::shared_ptr<pdat::CellData<double> > cd_c_b,
        std::shared_ptr<hier::Patch> patch)
{
   // tbox::plog << "computePhaseConcentrationsOnPatch..." << std::endl;
   assert(cd_temperature);
   assert(cd_phi);
   assert(cd_concentration);
   assert(cd_c_l);
   assert(cd_c_a);
   assert(cd_phi->getGhostCellWidth() == cd_concentration->getGhostCellWidth());
   assert(cd_phi->getDepth() > 1);

   (void)cd_c_b;

   const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());

   const int norderp = cd_phi->getDepth();
   assert(norderp > 1);

   double* ptr_temp = cd_temperature->getPointer();
   double* ptr_conc = cd_concentration->getPointer();
   double* ptr_c_l = cd_c_l->getPointer();
   double* ptr_c_a = cd_c_a->getPointer();

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

   // Assuming phi, and concentration all have same box
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

   // Assuming c_l, c_a, have same box
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

   int imin = c_i_gbox.lower(0);
   int imax = c_i_gbox.upper(0);
   int jmin = c_i_gbox.lower(1);
   int jmax = c_i_gbox.upper(1);
   int kmin = 0;
   int kmax = 0;
#if (NDIM == 3)
   kmin = c_i_gbox.lower(2);
   kmax = c_i_gbox.upper(2);
#endif

   std::vector<double*> ptr_phi(norderp);
   for (short i = 0; i < norderp; i++)
      ptr_phi[i] = cd_phi->getPointer(i);

   for (int kk = kmin; kk <= kmax; kk++) {
      for (int jj = jmin; jj <= jmax; jj++) {
         for (int ii = imin; ii <= imax; ii++) {

            const int idx_temp = (ii - imin_temp) + (jj - jmin_temp) * jp_temp +
                                 (kk - kmin_temp) * kp_temp;

            const int idx_pf = (ii - imin_pf) + (jj - jmin_pf) * jp_pf +
                               (kk - kmin_pf) * kp_pf;

            const int idx_c_i = (ii - imin_c_i) + (jj - jmin_c_i) * jp_c_i +
                                (kk - kmin_c_i) * kp_c_i;

            // last order parameter is assumed liquid
            double hphi[2] = {0., 0.};
            const double phiL = std::max(0., ptr_phi[norderp - 1][idx_pf]);
            hphi[0] = phiL * phiL;

            for (short i = 0; i < norderp - 1; i++) {
               const double phi = std::max(0., ptr_phi[i][idx_pf]);
               hphi[1] += phi * phi;
            }

            const double sum2 = hphi[0] + hphi[1];
            assert(sum2 > 0.);
            const double sum2inv = 1. / sum2;
            for (short i = 0; i < 2; i++)
               hphi[i] *= sum2inv;

            double c = ptr_conc[idx_pf];

            double t = ptr_temp[idx_temp];

            double x[2];
            int status = computePhaseConcentrations(t, &c, hphi, x);
            if (status < 0) {
               std::cerr << "computePhaseConcentrations failed for T=" << t
                         << ", " << norderp << " phases, hphi=";
               std::cerr << hphi << ", ";
               std::cerr << ", c=" << c << std::endl;
               std::cerr << "x=" << x[0] << "," << x[1] << std::endl;
               MPI_Abort(mpi.getCommunicator(), -1);
            }

            ptr_c_l[idx_c_i] = x[0];
            ptr_c_a[idx_c_i] = x[1];
         }
      }
   }

   return 0;
}
