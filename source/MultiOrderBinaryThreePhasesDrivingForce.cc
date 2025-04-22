#include "MultiOrderBinaryThreePhasesDrivingForce.h"

#include "SAMRAI/pdat/CellData.h"


using namespace SAMRAI;

MultiOrderBinaryThreePhasesDrivingForce::
    MultiOrderBinaryThreePhasesDrivingForce(
        FreeEnergyStrategyBinary* fenergy_strategy, const int norderp_A)
    : d_fenergy_strategy(fenergy_strategy), d_norderp_A(norderp_A)
{
}

void MultiOrderBinaryThreePhasesDrivingForce::addDrivingForce(
    hier::Patch& patch, const int temperature_id, const int phase_id,
    const int conc_l_id, const int conc_a_id, const int conc_b_id,
    const int f_l_id, const int f_a_id, const int f_b_id, const int rhs_id)
{
   assert(phase_id >= 0);
   assert(f_l_id >= 0);
   assert(f_a_id >= 0);
   assert(f_b_id >= 0);
   assert(rhs_id >= 0);
   assert(conc_l_id >= 0);
   assert(conc_a_id >= 0);
   assert(conc_b_id >= 0);
   assert(temperature_id >= 0);

   std::shared_ptr<pdat::CellData<double> > phase(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(phase_id)));
   assert(phase);
   assert(phase->getDepth() > 1);

   std::shared_ptr<pdat::CellData<double> > temperature(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(temperature_id)));
   assert(temperature);

   std::shared_ptr<pdat::CellData<double> > fl(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(f_l_id)));
   assert(fl);

   std::shared_ptr<pdat::CellData<double> > fa(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(f_a_id)));
   assert(fa);

   std::shared_ptr<pdat::CellData<double> > fb(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(f_b_id)));
   assert(fb);

   std::shared_ptr<pdat::CellData<double> > cl(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(conc_l_id)));
   assert(cl);

   std::shared_ptr<pdat::CellData<double> > ca(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(conc_a_id)));
   assert(ca);

   std::shared_ptr<pdat::CellData<double> > cb(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(conc_b_id)));
   assert(cb);

   std::shared_ptr<pdat::CellData<double> > rhs(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(rhs_id)));

   assert(fl->getGhostCellWidth()[0] == fa->getGhostCellWidth()[0]);
   assert(cl->getGhostCellWidth()[0] == ca->getGhostCellWidth()[0]);
   assert(phase->getDepth() > 1);
   assert(phase->getDepth() == rhs->getDepth());

   const hier::Box& pbox(patch.getBox());
   const int norderp = phase->getDepth();

   double* ptr_temp = temperature->getPointer();
   double* ptr_f_l = fl->getPointer();
   double* ptr_f_a = fa->getPointer();
   double* ptr_f_b = fb->getPointer();
   double* ptr_c_l = cl->getPointer();
   double* ptr_c_a = ca->getPointer();
   double* ptr_c_b = cb->getPointer();

   const hier::Box& rhs_gbox = rhs->getGhostBox();
   int imin_rhs = rhs_gbox.lower(0);
   int jmin_rhs = rhs_gbox.lower(1);
   int jp_rhs = rhs_gbox.numberCells(0);
   int kmin_rhs = 0;
   int kp_rhs = 0;
#if (NDIM == 3)
   kmin_rhs = rhs_gbox.lower(2);
   kp_rhs = jp_rhs * rhs_gbox.numberCells(1);
#endif

   const hier::Box& temp_gbox = temperature->getGhostBox();
   int imin_temp = temp_gbox.lower(0);
   int jmin_temp = temp_gbox.lower(1);
   int jp_temp = temp_gbox.numberCells(0);
   int kmin_temp = 0;
   int kp_temp = 0;
#if (NDIM == 3)
   kmin_temp = temp_gbox.lower(2);
   kp_temp = jp_temp * temp_gbox.numberCells(1);
#endif

   const hier::Box& pf_gbox = phase->getGhostBox();
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
   const hier::Box& f_i_gbox = fl->getGhostBox();
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
   const hier::Box& c_i_gbox = cl->getGhostBox();
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

   std::vector<double> rhs_local(norderp);

   std::vector<double*> ptr_rhs(norderp);
   for (short i = 0; i < norderp; i++)
      ptr_rhs[i] = rhs->getPointer(i);

   std::vector<double*> ptr_phi(norderp);
   for (short i = 0; i < norderp; i++)
      ptr_phi[i] = phase->getPointer(i);

   // tbox::plog<<"d_norderp_A = "<<d_norderp_A<<std::endl;

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
            double fb = ptr_f_b[idx_f_i];
            double cl = ptr_c_l[idx_c_i];
            double ca = ptr_c_a[idx_c_i];
            double cb = ptr_c_b[idx_c_i];

            double muL = d_fenergy_strategy->computeMuL(t, cl);

            //
            // see Moelans, Acta Mat 59 (2011)
            //

            // driving forces, assuming muL = muA = muB
            double dfa = (fa - muL * ca);
            double dfb = (fb - muL * cb);
            double dfl = (fl - muL * cl);
            assert(!std::isnan(dfa));

            // interpolation polynomials
            double hphiA = 0.;
            for (short i = 0; i < d_norderp_A; i++) {
               const double phi = std::max(0., ptr_phi[i][idx_pf]);
               hphiA += phi * phi;
            }
            assert(!std::isnan(hphiA));

            double hphiB = 0.;
            for (short i = d_norderp_A; i < norderp - 1; i++) {
               const double phi = std::max(0., ptr_phi[i][idx_pf]);
               hphiB += phi * phi;
            }
            assert(!std::isnan(hphiB));

            const double phiL = std::max(0., ptr_phi[norderp - 1][idx_pf]);
            double hphil = phiL * phiL;

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
               rhs_local[i] = 2. * std::max(0., ptr_phi[i][idx_pf]) *
                              ((1. - hphiA) * dfa - hphil * dfl - hphiB * dfb) *
                              sum2inv;
            // solid phase B order parameters
            for (short i = d_norderp_A; i < norderp - 1; i++)
               rhs_local[i] = 2. * std::max(0., ptr_phi[i][idx_pf]) *
                              ((1. - hphiB) * dfb - hphil * dfl - hphiA * dfa) *
                              sum2inv;
            // liquid phase order parameter
            rhs_local[norderp - 1] =
                2. * std::max(0., ptr_phi[norderp - 1][idx_pf]) *
                ((1. - hphil) * dfl - hphiA * dfa - hphiB * dfb) * sum2inv;
            for (short i = 0; i < norderp; i++)
               assert(!std::isnan(rhs_local[i]));

            for (short i = 0; i < norderp; i++)
               ptr_rhs[i][idx_rhs] -= (rhs_local[i]);
         }
      }
   }
}
