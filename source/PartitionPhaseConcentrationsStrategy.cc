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
#include "PartitionPhaseConcentrationsStrategy.h"
#include "FuncFort.h"

int PartitionPhaseConcentrationsStrategy::computePhaseConcentrationsOnPatch(
    std::shared_ptr<pdat::CellData<double> > cd_temperature,
    std::shared_ptr<pdat::CellData<double> > cd_phi,
    std::shared_ptr<pdat::CellData<double> > cd_concentration,
    std::shared_ptr<pdat::CellData<double> > cd_c_l,
    std::shared_ptr<pdat::CellData<double> > cd_c_a,
    std::shared_ptr<pdat::CellData<double> > cd_c_b,
    std::shared_ptr<hier::Patch> patch)
{
   (void)cd_temperature;
   (void)cd_c_b;

   assert(cd_phi);
   assert(cd_concentration);
   assert(cd_c_l);
   assert(cd_c_a);
   assert(d_partition_coeff_id >= 0);

   const hier::Box& pbox = patch->getBox();

   std::shared_ptr<pdat::CellData<double> > cd_partition_coeff(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch->getPatchData(d_partition_coeff_id)));
   assert(cd_partition_coeff);

   const double* const ptr_partition = cd_partition_coeff->getPointer();
   const double* const ptr_phi = cd_phi->getPointer();

   const hier::Box& partition_gbox = cd_partition_coeff->getGhostBox();
   int imin_partition = partition_gbox.lower(0);
   int jmin_partition = partition_gbox.lower(1);
   int jp_partition = partition_gbox.numberCells(0);
   int kmin_partition = 0;
   int kp_partition = 0;
#if (NDIM == 3)
   kmin_partition = partition_gbox.lower(2);
   kp_partition = jp_partition * partition_gbox.numberCells(1);
#endif

   // Assuming phi and concentration all have same box
   assert(cd_phi->getGhostCellWidth()[0] ==
          cd_concentration->getGhostCellWidth()[0]);
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

   // Assuming c_l, c_a, and c_b all have same box
   assert(cd_c_l->getGhostCellWidth()[0] == cd_c_a->getGhostCellWidth()[0]);
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

   const char interp = Thermo4PFM::concInterpChar(d_phase_interp_func_type);

   for (int kk = kmin; kk <= kmax; kk++) {
      for (int jj = jmin; jj <= jmax; jj++) {
         for (int ii = imin; ii <= imax; ii++) {

            const int idx_partition = (ii - imin_partition) +
                                      (jj - jmin_partition) * jp_partition +
                                      (kk - kmin_partition) * kp_partition;

            const int idx_pf = (ii - imin_pf) + (jj - jmin_pf) * jp_pf +
                               (kk - kmin_pf) * kp_pf;

            const int idx_c_i = (ii - imin_c_i) + (jj - jmin_c_i) * jp_c_i +
                                (kk - kmin_c_i) * kp_c_i;

            const double k = ptr_partition[idx_partition];
            assert(k == k);
            assert(k >= 0.);
            assert(k <= 1.);

            const double phi = ptr_phi[idx_pf];
            assert(phi == phi);

            const double hphi = INTERP_FUNC(phi, &interp);
            assert(hphi >= 0.);
            // if( hphi>1. )std::cerr<<"phi="<<phi<<std::endl;
            // if( hphi>1. )std::cerr<<"hphi="<<hphi<<std::endl;
            assert(hphi <= 1.00000001);

            // loop over atomic species
            for (int ic = 0; ic < cd_concentration->getDepth(); ic++) {
               const double* const ptr_conc = cd_concentration->getPointer(ic);
               const double c = ptr_conc[idx_pf];

               //               assert( c>=0. );
               //               assert( c<=1. );
               assert(k >= 0.);

               double* ptr_c_l = cd_c_l->getPointer(ic);
               double* ptr_c_a = cd_c_a->getPointer(ic);

               double factor = 1. / (1. - hphi + k * hphi);

               ptr_c_l[idx_c_i] = c * factor;
               ptr_c_a[idx_c_i] = k * c * factor;

               //               assert( ptr_c_l[idx_c_i]>=0. );
               //               assert( ptr_c_l[idx_c_i]<=1. );
               //               assert( ptr_c_a[idx_c_i]>=0. );
               //               assert( ptr_c_a[idx_c_i]<=1. );
            }  // ic
         }
      }
   }

   return 0;
}
