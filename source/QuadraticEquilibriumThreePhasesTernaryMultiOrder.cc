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
#include "QuadraticEquilibriumThreePhasesTernaryMultiOrder.h"
#include "FuncFort.h"


QuadraticEquilibriumThreePhasesTernaryMultiOrder::
    QuadraticEquilibriumThreePhasesTernaryMultiOrder(
        const short norderp_A, const int conc_l_id, const int conc_a_id,
        const int conc_b_id, const QuatModelParameters& model_parameters,
        std::shared_ptr<tbox::Database> conc_db)
    : PhaseConcentrationsStrategy(conc_l_id, conc_a_id, conc_b_id),
      d_norderp_A(norderp_A),
      d_conc_interp_func_type(model_parameters.conc_interp_func_type())
{
   std::shared_ptr<tbox::Database> quad_db = conc_db->getDatabase("Quadratic");

   double A_liquid[2];
   quad_db->getDoubleArray("A_liquid", &A_liquid[0], 2);
   double Ceq_liquid[2];
   quad_db->getDoubleArray("Ceq_liquid", &Ceq_liquid[0], 2);

   double A_solid_A[2];
   quad_db->getDoubleArray("A_solid_A", &A_solid_A[0], 2);
   double Ceq_solid_A[2];
   quad_db->getDoubleArray("Ceq_solid_A", &Ceq_solid_A[0], 2);

   double A_solid_B[2];
   quad_db->getDoubleArray("A_solid_B", &A_solid_B[0], 2);
   double Ceq_solid_B[2];
   quad_db->getDoubleArray("Ceq_solid_B", &Ceq_solid_B[0], 2);

   d_fenergy.reset(
       new Thermo4PFM::QuadraticFreeEnergyFunctionsTernaryThreePhase(
           A_liquid, Ceq_liquid, A_solid_A, Ceq_solid_A, A_solid_B, Ceq_solid_B,
           model_parameters.energy_interp_func_type(),
           model_parameters.conc_interp_func_type()));
}

int QuadraticEquilibriumThreePhasesTernaryMultiOrder::
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
   assert(cd_c_b);
   assert(cd_phi->getGhostCellWidth() == cd_concentration->getGhostCellWidth());
   assert(cd_phi->getDepth() > 1);

   const int norderp = cd_phi->getDepth();

   double* ptr_temp = cd_temperature->getPointer();
   double* ptr_conc = cd_concentration->getPointer();
   double* ptr_c_l = cd_c_l->getPointer();
   double* ptr_c_a = cd_c_a->getPointer();
   double* ptr_c_b = cd_c_b->getPointer();

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
   const short norder = cd_phi->getDepth();

   std::vector<double*> ptr_phi(norderp);
   for (short i = 0; i < norderp; i++)
      ptr_phi[i] = cd_phi->getPointer(i);

   // number of cells for each field
   const size_t coffset = c_i_gbox.size();
   const size_t ncp = pf_gbox.size();

   for (int kk = kmin; kk <= kmax; kk++) {
      for (int jj = jmin; jj <= jmax; jj++) {
         for (int ii = imin; ii <= imax; ii++) {

            const int idx_temp = (ii - imin_temp) + (jj - jmin_temp) * jp_temp +
                                 (kk - kmin_temp) * kp_temp;

            const int idx_pf = (ii - imin_pf) + (jj - jmin_pf) * jp_pf +
                               (kk - kmin_pf) * kp_pf;

            const int idx_c_i = (ii - imin_c_i) + (jj - jmin_c_i) * jp_c_i +
                                (kk - kmin_c_i) * kp_c_i;

            // 3 phases, last order parameter is assumed liquid
            double hphi[3] = {0., 0., 0.};

            const double phiL = std::max(0., ptr_phi[norderp - 1][idx_pf]);
            hphi[0] = phiL * phiL;

            for (short i = 0; i < d_norderp_A; i++) {
               const double phi = std::max(0., ptr_phi[i][idx_pf]);
               hphi[1] += phi * phi;
            }
            for (short i = d_norderp_A; i < norderp - 1; i++) {
               const double phi = std::max(0., ptr_phi[i][idx_pf]);
               hphi[2] += phi * phi;
            }

            const double sum2 = hphi[0] + hphi[1] + hphi[2];
            const double sum2inv = 1. / sum2;
            for (short i = 0; i < 3; i++)
               hphi[i] *= sum2inv;

            // get local two components of ternary alloy
            double c[2] = {ptr_conc[idx_pf], ptr_conc[idx_pf + ncp]};

            double t = ptr_temp[idx_temp];

            // solve KKS equations
            double x[6];
            d_fenergy->computePhaseConcentrations(t, c, hphi, x);

            ptr_c_l[idx_c_i] = x[0];
            ptr_c_l[idx_c_i + coffset] = x[1];
            ptr_c_a[idx_c_i] = x[2];
            ptr_c_a[idx_c_i + coffset] = x[3];
            ptr_c_b[idx_c_i] = x[4];
            ptr_c_b[idx_c_i + coffset] = x[5];
         }
      }
   }

   return 0;
}
