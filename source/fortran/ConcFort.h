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
#include "fc_internal_mangle.h"

// Function argument list interfaces
extern "C" {

void CONCENTRATIONFLUX(const int& ifirst0, const int& ilast0,
                       const int& ifirst1, const int& ilast1,
#if (NDIM == 3)
                       const int& ifirst2, const int& ilast2,
#endif
                       const double* dx, const double* conc, const int& ngconc,
                       const double* phase, const int& ngphase,
                       const double* diffconc0, const double* diffconc1,
#if (NDIM == 3)
                       const double* diffconc2,
#endif
                       const int& ngdiffconc, const double* dphicoupl0,
                       const double* dphicoupl1,
#if (NDIM == 3)
                       const double* dphicoupl2,
#endif
                       const int& ngdphicoupl, const double* flux0,
                       const double* flux1,
#if (NDIM == 3)
                       const double* flux2,
#endif
                       const int& ngflux);


void ADD_CAHNHILLIARDDOUBLEWELL_FLUX(const int& ifirst0, const int& ilast0,
                                     const int& ifirst1, const int& ilast1,
#if (NDIM == 3)
                                     const int& ifirst2, const int& ilast2,
#endif
                                     const double* dx, const double* conc,
                                     const int& ngconc, const double& m,
                                     const double& ca, const double& cb,
                                     const double& well_scale,
                                     const double& kappa, const double* flux0,
                                     const double* flux1,
#if (NDIM == 3)
                                     const double* flux2,
#endif
                                     const int& ngflux);

void ADD_WANG_SINTERING_FLUX(
    const int& ifirst0, const int& ilast0, const int& ifirst1,
    const int& ilast1,
#if (NDIM == 3)
    const int& ifirst2, const int& ilast2,
#endif
    const double* dx, const double* conc, const int& ngconc,
    const double* diffconc0, const double* diffconc1,
#if (NDIM == 3)
    const double* diffconc2,
#endif
    const int& ngdiffconc, const double* phi, const int& ngphi,
    const int& norder, const double& mobility, const double& a, const double& b,
    const double& kappa, const double* flux0, const double* flux1,
#if (NDIM == 3)
    const double* flux2,
#endif
    const int& ngflux, const double* tmp1, const double* tmp2);

void ADDCONCENTRATIONFLUXFROMGRADT(const int& ifirst0, const int& ilast0,
                                   const int& ifirst1, const int& ilast1,
#if (NDIM == 3)
                                   const int& ifirst2, const int& ilast2,
#endif
                                   const double* dx, const double* temperature,
                                   const int& ngtemperature, const double* mq0,
                                   const double* mq1,
#if (NDIM == 3)
                                   const double* mq2,
#endif
                                   const int& ngmq, const double* flux0,
                                   const double* flux1,
#if (NDIM == 3)
                                   const double* flux2,
#endif
                                   const int& ngflux,
                                   const char* const avg_type);

void ADDCONCENTRATIONFLUXFROMANTITRAPPING(
    const int& ifirst0, const int& ilast0, const int& ifirst1,
    const int& ilast1,
#if (NDIM == 3)
    const int& ifirst2, const int& ilast2,
#endif
    const double* dx, const double* phase, const int& ngphase, const double* cl,
    const double* ca, const int& ngc, const int& ncomp, const double* dphi,
    const int& ngdphi, const double& alpha, const double* flux0,
    const double* flux1,
#if (NDIM == 3)
    const double* flux2,
#endif
    const int& ngflux);

void ADDCONCENTRATIONFLUXFROMANTITRAPPING3PHASES(
    const int& ifirst0, const int& ilast0, const int& ifirst1,
    const int& ilast1,
#if (NDIM == 3)
    const int& ifirst2, const int& ilast2,
#endif
    const double* dx, const double* phase, const int& ngphase, const double* cl,
    const double* ca, const double* cb, const int& ngc, const double* dphi,
    const int& ngdphi, const double& alpha, const double* flux0,
    const double* flux1,
#if (NDIM == 3)
    const double* flux2,
#endif
    const int& ngflux);

void ADDCONCENTRATIONFLUXFROMANTITRAPPINGMULTIORDERP(
    const int& ifirst0, const int& ilast0, const int& ifirst1,
    const int& ilast1,
#if (NDIM == 3)
    const int& ifirst2, const int& ilast2,
#endif
    const double* dx, const double* phase, const int& ngphase,
    const int& nphases, const double* cl, const double* ca, const int& ngc,
    const double* dphi, const int& ngdphi, const double& alpha,
    const double* flux0, const double* flux1,
#if (NDIM == 3)
    const double* flux2,
#endif
    const int& ngflux);

void ADD_FLUX(const int& ifirst0, const int& ilast0, const int& ifirst1,
              const int& ilast1,
#if (NDIM == 3)
              const int& ifirst2, const int& ilast2,
#endif
              const double* dx, const double* conc, const int& ngconc,
              const int& ncomp, const double* d0, const double* d1,
#if (NDIM == 3)
              const double* d2,
#endif
              const int& ngd, const double* flux0, const double* flux1,
#if (NDIM == 3)
              const double* flux2,
#endif
              const int& ngflux);
void ADD_FLUX_ISO(const int& ifirst0, const int& ilast0, const int& ifirst1,
                  const int& ilast1,
#if (NDIM == 3)
                  const int& ifirst2, const int& ilast2,
#endif
                  const double* dx, const double* conc, const int& ngconc,
                  const int& ncomp, const double* d0, const double* d1,
#if (NDIM == 3)
                  const double* d2,
#endif
                  const int& ngd, const double* flux0, const double* flux1,
#if (NDIM == 3)
                  const double* flux2,
#endif
                  const int& ngflux);

void ADD_FLUX_4TH(const int& ifirst0, const int& ilast0, const int& ifirst1,
                  const int& ilast1,
#if (NDIM == 3)
                  const int& ifirst2, const int& ilast2,
#endif
                  const double* dx, const double* conc, const int& ng,
                  const int& depth, const double* d0, const double* d1,
#if (NDIM == 3)
                  const double* d2,
#endif
                  const int& ngd, const double* flux0, const double* flux1,
#if (NDIM == 3)
                  const double* flux2,
#endif
                  const int& ngflux, const int* physb);

void CONCENTRATION_FLUX_SPINODAL(
    const int& ifirst0, const int& ilast0, const int& ifirst1,
    const int& ilast1,
#if (NDIM == 3)
    const int& ifirst2, const int& ilast2,
#endif
    const double* dx, const double* conc, const int& ngconc, const int& ncomp,
    const double* conca, const int& ngconca, const double* concb,
    const int& ngconcb, const double* diffconc0, const double* diffconc1,
#if (NDIM == 3)
    const double* diffconc2,
#endif
    const int& ngdiffconcl, const double* eta, const int& ngeta,
    const double& kappa, const double* flux0, const double* flux1,
#if (NDIM == 3)
    const double* flux2,
#endif
    const int& ngflux);

void PFMDIFFUSION(const int& ifirst0, const int& ilast0, const int& ifirst1,
                  const int& ilast1,
#if (NDIM == 3)
                  const int& ifirst2, const int& ilast2,
#endif
                  const double* phi, const int& ngphi, const double* diff0,
                  const double* diff1,
#if (NDIM == 3)
                  const double* diff2,
#endif
                  const int& ngdiff, const double* t, const int& ngt,
                  const double& d_phase0, const double& q0_phase0,
                  const double& d_phase1, const double& q0_phase1,
                  const double& gas_constant_R, const char* phi_interp_type,
                  const char* avg_func_type);

void PFMDIFFUSION_OF_TEMPERATURE(
    const int& ifirst0, const int& ilast0, const int& ifirst1,
    const int& ilast1,
#if (NDIM == 3)
    const int& ifirst2, const int& ilast2,
#endif
    const double* phi, const int& ngphi, const double* diffl0,
    const double* diffl1,
#if (NDIM == 3)
    const double* diffl2,
#endif
    const double* diffa0, const double* diffa1,
#if (NDIM == 3)
    const double* diffa2,
#endif
    const int& ngdiff, const double* t, const int& ngt, const double& d_phase0,
    const double& q0_phase0, const double& d_phase1, const double& q0_phase1,
    const double& gas_constant_R, const char* phi_interp_type,
    const char* avg_func_type);

void PFMDIFFUSION_SCALAR_2PHASES(const int& ifirst0, const int& ilast0,
                                 const int& ifirst1, const int& ilast1,
#if (NDIM == 3)
                                 const int& ifirst2, const int& ilast2,
#endif
                                 const double* phiL, const double* phiA,
                                 const int& nphiA, const int& ngphi,
                                 const double* diff0, const double* diff1,
#if (NDIM == 3)
                                 const double* diff2,
#endif
                                 const int& ngdiff, const double& d0_phaseL,
                                 const double& d0_phaseA,
                                 const char* phi_interp_type,
                                 const char* avg_func_type);

void PFMDIFFUSION_SCALAR(const int& ifirst0, const int& ilast0,
                         const int& ifirst1, const int& ilast1,
#if (NDIM == 3)
                         const int& ifirst2, const int& ilast2,
#endif
                         const double* phi, const int& ngphi,
                         const double* diff0, const double* diff1,
#if (NDIM == 3)
                         const double* diff2,
#endif
                         const int& ngdiff, const double& d0_phaseL,
                         const double& d0_phaseA, const char* phi_interp_type,
                         const char* avg_func_type);

void PFMDIFFUSION_OF_TEMPERATURE_MULTIORDER(
    const int& ifirst0, const int& ilast0, const int& ifirst1,
    const int& ilast1,
#if (NDIM == 3)
    const int& ifirst2, const int& ilast2,
#endif
    const double* phi, const int& nphi, const int& ngphi, const double* diffl0,
    const double* diffl1,
#if (NDIM == 3)
    const double* diffl2,
#endif
    const double* diffa0, const double* diffa1,
#if (NDIM == 3)
    const double* diffa2,
#endif
    const int& ngdiff, const double* t, const int& ngt, const double& d_phase0,
    const double& q0_phase0, const double& d_phase1, const double& q0_phase1,
    const double& gas_constant_R, const char* avg_func_type);

void PFMDIFFUSION_OF_TEMPERATURE_MULTIORDER_THREEPHASES(
    const int& ifirst0, const int& ilast0, const int& ifirst1,
    const int& ilast1,
#if (NDIM == 3)
    const int& ifirst2, const int& ilast2,
#endif
    const double* phiL, const int& nphiL, const double* phiA, const int& nphiA,
    const double* phiB, const int& nphiB, const int& ngphi,
    const double* diffl0, const double* diffl1,
#if (NDIM == 3)
    const double* diffl2,
#endif
    const double* diffa0, const double* diffa1,
#if (NDIM == 3)
    const double* diffa2,
#endif
    const double* diffb0, const double* diffb1,
#if (NDIM == 3)
    const double* diffb2,
#endif
    const int& ngdiff, const double* t, const int& ngt, const double& d_phase0,
    const double& q0_phase0, const double& d_phase1, const double& q0_phase1,
    const double& d_phase2, const double& q0_phase2,
    const double& gas_constant_R, const char* avg_func_type);

void PFMDIFFUSION_OF_TEMPERATURE_FOLCHPLAPP(
    const int& ifirst0, const int& ilast0, const int& ifirst1,
    const int& ilast1,
#if (NDIM == 3)
    const int& ifirst2, const int& ilast2,
#endif
    const double* phi, const int& ngphi, const double* diffl0,
    const double* diffl1,
#if (NDIM == 3)
    const double* diffl2,
#endif
    const double* diffa0, const double* diffa1,
#if (NDIM == 3)
    const double* diffa2,
#endif
    const double* diffb0, const double* diffb1,
#if (NDIM == 3)
    const double* diffb2,
#endif
    const int& ngdiff, const double* t, const int& ngt, const double& d_phase0,
    const double& q0_phase0, const double& d_phase1, const double& q0_phase1,
    const double& d_phase2, const double& q0_phase2,
    const double& gas_constant_R, const char* phi_interp_type,
    const char* avg_func_type);


void PFMDIFFUSION_SCALAR_MULTIORDER_3PHASES(
    const int& ifirst0, const int& ilast0, const int& ifirst1,
    const int& ilast1,
#if (NDIM == 3)
    const int& ifirst2, const int& ilast2,
#endif
    const double* phiL, const int& nphiL, const double* phiA, const int& nphiA,
    const double* phiB, const int& nphiB, const int& ngphi, const double* diff0,
    const double* diff1,
#if (NDIM == 3)
    const double* diff2,
#endif
    const int& ngdiff, const double& d_phaseL, const double& d_phaseA,
    const double& d_phaseB, const char* phi_interp_type,
    const char* avg_func_type);


void AB_DIFFUSION_OF_TEMPERATURE(
    const int& ifirst0, const int& ilast0, const int& ifirst1,
    const int& ilast1,
#if (NDIM == 3)
    const int& ifirst2, const int& ilast2,
#endif
    const double* phia, const int& nphia, const double* phib, const int& nphib,
    const int& ngphi, const double* diffa0, const double* diffa1,
#if (NDIM == 3)
    const double* diffa2,
#endif
    const double* diffb0, const double* diffb1,
#if (NDIM == 3)
    const double* diffb2,
#endif
    const int& ngdiff, const double* t, const int& ngt, const double& d0,
    const double& q0, const double& gas_constant_R, const int& dupl);

void ADD_AB_DIFFUSION(const int& ifirst0, const int& ilast0, const int& ifirst1,
                      const int& ilast1,
#if (NDIM == 3)
                      const int& ifirst2, const int& ilast2,
#endif
                      const double* phia, const int& nphia, const double* phib,
                      const int& nphib, const int& ngphi, const double* diff0,
                      const double* diff1,
#if (NDIM == 3)
                      const double* diff2,
#endif
                      const int& ngdiff, const double& d0, const int& dupl);

void ADD_AB_DIFFUSION_SINGLE(const int& ifirst0, const int& ilast0,
                             const int& ifirst1, const int& ilast1,
#if (NDIM == 3)
                             const int& ifirst2, const int& ilast2,
#endif
                             const double* phi, const int& ngphi,
                             const double* diff0, const double* diff1,
#if (NDIM == 3)
                             const double* diff2,
#endif
                             const int& ngdiff, const double& d0);

void ADD_AB_DIFFUSION_OF_TEMPERATURE_SINGLE(
    const int& ifirst0, const int& ilast0, const int& ifirst1,
    const int& ilast1,
#if (NDIM == 3)
    const int& ifirst2, const int& ilast2,
#endif
    const double* phi, const int& ngphi, const double* temp, const int& ngtemp,
    const double* diff0, const double* diff1,
#if (NDIM == 3)
    const double* diff2,
#endif
    const int& ngdiff, const double& d0, const double& q0,
    const double& gas_constant_R);

void ADD_DIFFUSION_SCALAR_OF_TEMPERATURE(
    const int& ifirst0, const int& ilast0, const int& ifirst1,
    const int& ilast1,
#if (NDIM == 3)
    const int& ifirst2, const int& ilast2,
#endif
    const double* phi, const int& ngphi, const double* temp, const int& ngtemp,
    const double* diff0, const double* diff1,
#if (NDIM == 3)
    const double* diff2,
#endif
    const int& ngdiff, const double& d0l, const double& q0l, const double& d0s,
    const double& q0s, const double& gas_constant_R,
    const char* phi_interp_type);

void CONCENTRATION_DIFFCOEFF_OF_TEMPERATURE(
    const int& ifirst0, const int& ilast0, const int& ifirst1,
    const int& ilast1,
#if (NDIM == 3)
    const int& ifirst2, const int& ilast2,
#endif
    const double* diffl0, const double* diffl1,
#if (NDIM == 3)
    const double* diffl2,
#endif
    const double* diffa0, const double* diffa1,
#if (NDIM == 3)
    const double* diffa2,
#endif
    const double* diffb0, const double* diffb1,
#if (NDIM == 3)
    const double* diffb2,
#endif
    const int& ngdiff, const double* t, const int& ngt, const double& d_phase0,
    const double& q0_phase0, const double& d_phase1, const double& q0_phase1,
    const double& d_phase2, const double& q0_phase2,
    const double& gas_constant_R, const int& with_three_phases);

void CONCENTRATIONDIFFUSION_BECKERMANN(const int& ifirst0, const int& ilast0,
                                       const int& ifirst1, const int& ilast1,
#if (NDIM == 3)
                                       const int& ifirst2, const int& ilast2,
#endif
                                       const double* phi, const int& ngphi,
                                       const double* diff0, const double* diff1,
#if (NDIM == 3)
                                       const double* diff2,
#endif
                                       const int& ngdiff,
                                       const double* partition_coeff,
                                       const int& ngk, const double& d_phase0,
                                       const double& d_phase1,
                                       const char* phi_interp_type,
                                       const char* avg_func_type);

void COMPUTERHSCONCENTRATION(const int& ifirst0, const int& ilast0,
                             const int& ifirst1, const int& ilast1,
#if (NDIM == 3)
                             const int& ifirst2, const int& ilast2,
#endif
                             const double* dx, const double* flux0,
                             const double* flux1,
#if (NDIM == 3)
                             const double* flux2,
#endif
                             const int& ngflux, const double& conc_mobility,
                             const double* rhs, const int& ngrhs);

void PHASERHS_FENERGY(const int& ifirst0, const int& ilast0, const int& ifirst1,
                      const int& ilast1,
#if (NDIM == 3)
                      const int& ifirst2, const int& ilast2,
#endif
                      const double* fl, const double* fa, const double* phi,
                      const int& ngphi, const double* rhs, const int& ngrhs,
                      const char* phi_interp_type);

void PHASERHS_FENERGY_MULTIORDERP(const int& ifirst0, const int& ilast0,
                                  const int& ifirst1, const int& ilast1,
#if (NDIM == 3)
                                  const int& ifirst2, const int& ilast2,
#endif
                                  const double* fl, const double* fa,
                                  const double* phi, const int& ngphi,
                                  const int& nphi, const double* rhs,
                                  const int& ngrhs);

void ETARHS_FENERGY(const int& ifirst0, const int& ilast0, const int& ifirst1,
                    const int& ilast1,
#if (NDIM == 3)
                    const int& ifirst2, const int& ilast2,
#endif
                    const double* fl, const double* fa, const double* fb,
                    const double* phi, const int& ngphi, const double* eta,
                    const int& ngeta, const double* rhs, const int& ngrhs,
                    const char* phi_interp_type, const char* eta_interp_type);

void CALPHAD_CONC_SOLV_TWO(const double* x, const double& c, const double& hphi,
                           const double& d_gas_constant_R, const double* L0,
                           const double* L1, const double* L2, const double* fA,
                           const double* fB);

int CALPHAD_CONC_CEQ_TWO(const double* x, const double& cl, const double& cs,
                         const double& d_gas_constant_R, const double* L0,
                         const double* L1, const double* L2, const double* fA,
                         const double* fB);

void CALPHAD_CONC_SOLV_THREE(const double* x, const double& c,
                             const double& hphi, const double& heta,
                             const double& d_gas_constant_R, const double* L0,
                             const double* L1, const double* L2,
                             const double* fA, const double* fB);

void INITGAUSSIAN(const double*, const double*, const int& ifirst0,
                  const int& ilast0, const int& ifirst1, const int& ilast1,
#if (NDIM == 3)
                  const int& ifirst2, const int& ilast2,
#endif
                  const double*, const int&, const double*, const double*,
                  const double&, const double&, const double&);

void INITGRADIENT(const double*, const double*, const int& ifirst0,
                  const int& ilast0, const int& ifirst1, const int& ilast1,
#if (NDIM == 3)
                  const int& ifirst2, const int& ilast2,
#endif
                  const double*, const int&, const double*, const double&,
                  const double*);

void INITGAUSSIANSOURCE(const double*, const double*, const int& ifirst0,
                        const int& ilast0, const int& ifirst1,
                        const int& ilast1,
#if (NDIM == 3)
                        const int&, const int&,
#endif
                        const double*, const int&, const double*, const int&,
                        const double*, const double*, const double&,
                        const double&);

void LINEARMELTINGLINE(const int& ifirst0, const int& ilast0,
                       const int& ifirst1, const int& ilast1,
#if (NDIM == 3)
                       const int&, const int&,
#endif
                       const double*, const int&, const double*, const int&,
                       const double&, const double&, const double&);

void COMPUTE_CONCENTRATION_FROM_PHASE_CONCENTRATIONS(
    const int& ifirst0, const int& ilast0, const int& ifirst1,
    const int& ilast1,
#if (NDIM == 3)
    const int&, const int&,
#endif
    const double*, const int&, const double*, const double*, const int&,
    const double*, const int&, const char* const);
}
