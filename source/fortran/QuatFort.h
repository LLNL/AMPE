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
void SETTOZERO(const int&, const int&, const int&, const int&, const int&,
               const int&, const int&, const int&,
#if (NDIM == 3)
               const int&, const int&, const int&, const int&,
#endif
               const double*);

void SIDE2CELL(const int&, const int&, const int&, const int&,
#if (NDIM == 3)
               const int&, const int&,
#endif
               const double*, const double*,
#if (NDIM == 3)
               const double*,
#endif
               const int&, const double*, const int&);

void SOURCE_TEMPERATURE(const int&, const int&, const int&, const int&,
#if (NDIM == 3)
                        const int&, const int&,
#endif
                        double* conc, const int&, double* rhs, const int&,
                        const double* const cp, const int&,
                        const double* const s, const int&);

void HEAT_CAPACITY_NKR(const int&, const int&, const int&, const int&,
#if (NDIM == 3)
                       const int&, const int&,
#endif
                       const double* const conc, const int&, const int&,
                       const double* const temp, const int&, double* cp,
                       const int&, const int* const cp_powers, const int&,
                       const double* const cpcoeffs, const int&);

void GRADIENT_FLUX(const int& ifirst0, const int& ilast0, const int& ifirst1,
                   const int& ilast1,
#if (NDIM == 3)
                   const int& ifirst2, const int& ilast2,
#endif
                   const double* dx, const double& epsilon, double* phase,
                   const int& phaseng, double* flux0, double* flux1,
#if (NDIM == 3)
                   double* flux2,
#endif
                   const int& fluxng);

void COMPUTE_FLUX_ISOTROPIC(const int& ifirst0, const int& ilast0,
                            const int& ifirst1, const int& ilast1,
#if (NDIM == 3)
                            const int& ifirst2, const int& ilast2,
#endif
                            const double* dx, const double& epsilon,
                            double* phase, const int& phaseng, double* flux0,
                            double* flux1,
#if (NDIM == 3)
                            double* flux2,
#endif
                            const int& fluxng);

void ANISOTROPIC_GRADIENT_FLUX(const int& ifirst0, const int& ilast0,
                               const int& ifirst1, const int& ilast1,
#if (NDIM == 3)
                               const int& ifirst2, const int& ilast2,
#endif
                               const double* dx, const double& epsilon,
                               const double& nu, const int& knumber,
                               double* phase, const int& phaseng, double* quat,
                               const int& qng, const int& qlen, double* flux0,
                               double* flux1,
#if (NDIM == 3)
                               double* flux2,
#endif
                               const int& fluxng);

void COMPUTERHSPBG(const int& ifirst0, const int& ilast0, const int& ifirst1,
                   const int& ilast1,
#if (NDIM == 3)
                   const int& ifirst2, const int& ilast2,
#endif
                   const double* dx, const double& misorientation_factor,
                   const double& epsilonq, const double* phase_flux0,
                   const double* phase_flux1,
#if (NDIM == 3)
                   const double* phase_flux2,
#endif
                   const int& ngphaseflux, const double* temperature,
                   const int& ngtemp, const double& phi_well_scale,
                   const double& eta_well_scale, const double* phi,
                   const int& ngphi, const double* eta, const int& ngeta,
                   const double* orient_grad_mod, const int& ngogm, double* rhs,
                   const int& ngrhs, const char* phi_well_type,
                   const char* eta_well_type, const char* phi_interp_type,
                   const char* orient_interp_type1,
                   const char* orient_interp_type2, const int& with_orient,
                   const int& three_phase);

void COMPUTERHSTHREEPHASES(const int& ifirst0, const int& ilast0,
                           const int& ifirst1, const int& ilast1,
#if (NDIM == 3)
                           const int& ifirst2, const int& ilast2,
#endif
                           const double* dx, const double* phase_flux0,
                           const double* phase_flux1,
#if (NDIM == 3)
                           const double* phase_flux2,
#endif
                           const int& ngphaseflux, const double& phi_well_scale,
                           const double* phi, const int& ngphi, double* rhs,
                           const int& ngrhs);

void COMPUTERHSMULTIORDER(const int& ifirst0, const int& ilast0,
                          const int& ifirst1, const int& ilast1,
#if (NDIM == 3)
                          const int& ifirst2, const int& ilast2,
#endif
                          const int& norder, const double* dx,
                          const double* phase_flux0, const double* phase_flux1,
#if (NDIM == 3)
                          const double* phase_flux2,
#endif
                          const int& ngphaseflux, const double& gamma,
                          const double& m, const double* phi, const int& ngphi,
                          double* tmp, double* rhs, const int& ngrhs);

void COMPUTERHS_WANG_SINTERING(const int& ifirst0, const int& ilast0,
                               const int& ifirst1, const int& ilast1,
#if (NDIM == 3)
                               const int& ifirst2, const int& ilast2,
#endif
                               const int& norder, const double* dx,
                               const double* phase_flux0,
                               const double* phase_flux1,
#if (NDIM == 3)
                               const double* phase_flux2,
#endif
                               const int& ngphaseflux,
                               const double& parameter_b, const double* phi,
                               const int& ngphi, const double* conc,
                               const int& ngconc, double* tmp, double* rhs,
                               const int& ngrhs);

void ADDVDPHIDX(const int& ifirst0, const int& ilast0, const int& ifirst1,
                const int& ilast1,
#if (NDIM == 3)
                const int& ifirst2, const int& ilast2,
#endif
                const double* dx, const double* phase, const int& ngphase,
                const double& vel, const double* rhs, const int& ngrhs);

void ADDVDPHIDX_UPWIND3(const int& ifirst0, const int& ilast0,
                        const int& ifirst1, const int& ilast1,
#if (NDIM == 3)
                        const int& ifirst2, const int& ilast2,
#endif
                        const double* dx, const double* phase,
                        const int& ngphase, const double& vel,
                        const double* rhs, const int& ngrhs, const int& bc);

void COMPUTERHSTEMP(const int& ifirst0, const int& ilast0, const int& ifirst1,
                    const int& ilast1,
#if (NDIM == 3)
                    const int& ifirst2, const int& ilast2,
#endif
                    const double* dx, const double& thermal_diffusivity,
                    const double& latent_heat, const double* temperature,
                    const int& ngtemp, const double* cpl, const int& ngcpl,
                    const int& with_phase, const double* phirhs,
                    const int& ngphirhs, double* rhs, const int& ngrhs);

void COMPUTERHSBIASWELL(const int& ifirst0, const int& ilast0,
                        const int& ifirst1, const int& ilast1,
#if (NDIM == 3)
                        const int& ifirst2, const int& ilast2,
#endif
                        const double*, const int&, const double*, const int&,
                        const double&, const double&, double* teq, const int&,
                        double* rhs, const int&);

void COMPUTERHSBIASWELLBECKERMANN(const int& ifirst0, const int& ilast0,
                                  const int& ifirst1, const int& ilast1,
#if (NDIM == 3)
                                  const int& ifirst2, const int& ilast2,
#endif
                                  const double*, const int&, const double*,
                                  const int&, const double&, double* teq,
                                  const int&, double* rhs, const int&);

void COMPUTERHSDELTATEMPERATURE(const int& ifirst0, const int& ilast0,
                                const int& ifirst1, const int& ilast1,
#if (NDIM == 3)
                                const int& ifirst2, const int& ilast2,
#endif
                                const double*, const int&, const double*,
                                const int&, const double&, const double&,
                                double* rhs, const int&, const char* const);

void COMPUTERHSETA(const int& ifirst0, const int& ilast0, const int& ifirst1,
                   const int& ilast1,
#if (NDIM == 3)
                   const int& ifirst2, const int& ilast2,
#endif
                   const double* dx, const double* eta_flux0,
                   const double* eta_flux1,
#if (NDIM == 3)
                   const double* eta_flux2,
#endif
                   const int& ngetaflux, const double* temperature,
                   const int& ngtemp, const double& eta_well_scale,
                   const double* phi, const int& ngphi, const double* eta,
                   const int& ngeta, double* rhs, const int& ngrhs,
                   const char* eta_well_type, const char* phi_interp_type);

void SMOOTHQUAT(const int& ifirst0, const int& ilast0, const int& ifirst1,
                const int& ilast1,
#if (NDIM == 3)
                const int& ifirst2, const int& ilast2,
#endif
                const int& depth, const double*, const int&, const double*,
                const int&, const double*, const double*,
#if (NDIM == 3)
                const double*,
#endif
                const int&, const double*, const int&, const double&);

void DIFFS(const int&, const int&, const int&, const int&,
#if (NDIM == 3)
           const int&, const int&,
#endif
           const double*, const int&, const double*, const double*,
#if (NDIM == 3)
           const double*,
#endif
           const int&);

void GRAD_CELL(const int&, const int&, const int&, const int&,
#if (NDIM == 3)
               const int&, const int&,
#endif
               const double*, const double*,
#if (NDIM == 3)
               const double*,
#endif
               const int&, const double*, const double*, const double*,
#if (NDIM == 3)
               const double*,
#endif
               const int&);

void GRAD_SIDE(const int&, const int&, const int&, const int&,
#if (NDIM == 3)
               const int&, const int&,
#endif
               const double*, const double*,
#if (NDIM == 3)
               const double*,
#endif
               const int&, const int&, const int&, const int&,
#if (NDIM == 3)
               const int&, const int&,
#endif
               const double*, const double*, const double*,
#if (NDIM == 3)
               const double*,
#endif
               const double*, const double*,
#if (NDIM == 3)
               const double*, const double*, const double*, const double*,
#endif
               const int&, const int&, const int&, const int&
#if (NDIM == 3)
               ,
               const int&, const int&
#endif
);

void GRAD_SIDE_ISOTROPIC(const int&, const int&, const int&, const int&,
#if (NDIM == 3)
                         const int&, const int&,
#endif
                         const double*, const double*,
#if (NDIM == 3)
                         const double*,
#endif
                         const int&, const int&, const int&, const int&,
#if (NDIM == 3)
                         const int&, const int&,
#endif
                         const double*, const double*, const double*,
#if (NDIM == 3)
                         const double*,
#endif
                         const double*, const double*,
#if (NDIM == 3)
                         const double*, const double*, const double*,
                         const double*,
#endif
                         const int&, const int&, const int&, const int&
#if (NDIM == 3)
                         ,
                         const int&, const int&
#endif
);

void QUAT_SYMM_ROTATION(const int& lo0, const int& hi0, const int& lo1,
                        const int& hi1,
#if (NDIM == 3)
                        const int& lo2, const int& hi2,
#endif
                        const double* q, const int& ngq, const int& depth,
                        const int* rot_x, const int* rot_y,
#if (NDIM == 3)
                        const int* rot_z,
#endif
                        const int& ngrot);

void QUAT_FUNDAMENTAL(const int&, const int&, const int&, const int&,
#if (NDIM == 3)
                      const int&, const int&,
#endif
                      const double*, const int&, const int&, const int&,
                      const int&,
#if (NDIM == 3)
                      const int&, const int&,
#endif
                      const int&);

void QUATDIFFS(const int& lo0, const int& hi0, const int& lo1, const int& hi1,
#if (NDIM == 3)
               const int& lo2, const int& hi2,
#endif
               const int& depth, const double* q, const int& ngq,
               const double* diff_x, const double* diff_y,
#if (NDIM == 3)
               const double* diff_z,
#endif
               const int& ngdiff);

void QUATDIFFS_SYMM(const int& lo0, const int& hi0, const int& lo1,
                    const int& hi1,
#if (NDIM == 3)
                    const int& lo2, const int& hi2,
#endif
                    const int& depth, const double* q, const int& ngq,
                    const double* diff_x, const double* diff_y,
#if (NDIM == 3)
                    const double* diff_z,
#endif
                    const int& ngdiff, const int* iqrot_x, const int* iqrot_y,
#if (NDIM == 3)
                    const int* iqrot_z,
#endif
                    const int& ngiq);

void QUATGRAD_CELL(const int& lo0, const int& hi0, const int& lo1,
                   const int& hi1,
#if (NDIM == 3)
                   const int& lo2, const int& hi2,
#endif
                   const int& depth, const double* h, const double* diff_x,
                   const double* diff_y,
#if (NDIM == 3)
                   const double* diff_z,
#endif
                   const int& ngdiff, const double* grad_x,
                   const double* grad_y,
#if (NDIM == 3)
                   const double* grad_z,
#endif
                   const int& nggrad);

void QUATGRAD_CELL_SYMM(const int& lo0, const int& hi0, const int& lo1,
                        const int& hi1,
#if (NDIM == 3)
                        const int& lo2, const int& hi2,
#endif
                        const int& depth, const double* h, const double* diff_x,
                        const double* diff_y,
#if (NDIM == 3)
                        const double* diff_z,
#endif
                        const int& ngdiff, const double* grad_x,
                        const double* grad_y,
#if (NDIM == 3)
                        const double* grad_z,
#endif
                        const int& nggrad, const int* iqrot_x,
                        const int* iqrot_y,
#if (NDIM == 3)
                        const int* iqrot_z,
#endif
                        const int& ngiq);

void QUATGRAD_SIDE(const int& lo0, const int& hi0, const int& lo1,
                   const int& hi1,
#if (NDIM == 3)
                   const int& lo2, const int& hi2,
#endif
                   const int& depth, const double* h, const double* diff_x,
                   const double* diff_y,
#if (NDIM == 3)
                   const double* diff_z,
#endif
                   const int& ngdiff, const double* grad_x_x,
                   const double* grad_y_x,
#if (NDIM == 3)
                   const double* grad_z_x,
#endif
                   const double* grad_x_y, const double* grad_y_y,
#if (NDIM == 3)
                   const double* grad_z_y, const double* grad_x_z,
                   const double* grad_y_z, const double* grad_z_z,
#endif
                   const int& nggrad);

void QUATGRAD_SIDE_ISOTROPIC(const int& lo0, const int& hi0, const int& lo1,
                             const int& hi1,
#if (NDIM == 3)
                             const int& lo2, const int& hi2,
#endif
                             const int& depth, const double* h,
                             const double* diff_x, const double* diff_y,
#if (NDIM == 3)
                             const double* diff_z,
#endif
                             const int& ngdiff, const double* grad_x_x,
                             const double* grad_y_x,
#if (NDIM == 3)
                             const double* grad_z_x,
#endif
                             const double* grad_x_y, const double* grad_y_y,
#if (NDIM == 3)
                             const double* grad_z_y, const double* grad_x_z,
                             const double* grad_y_z, const double* grad_z_z,
#endif
                             const int& nggrad);

void QUATGRAD_SIDE_SYMM(const int& lo0, const int& hi0, const int& lo1,
                        const int& hi1,
#if (NDIM == 3)
                        const int& lo2, const int& hi2,
#endif
                        const int& depth, const double* h, const double* diff_x,
                        const double* diff_y,
#if (NDIM == 3)
                        const double* diff_z,
#endif
                        const int& ngdiff, const double* grad_x_x,
                        const double* grad_y_x,
#if (NDIM == 3)
                        const double* grad_z_x,
#endif
                        const double* grad_x_y, const double* grad_y_y,
#if (NDIM == 3)
                        const double* grad_z_y, const double* grad_x_z,
                        const double* grad_y_z, const double* grad_z_z,
#endif
                        const int& nggrad, const int* iqrot_x,
                        const int* iqrot_y,
#if (NDIM == 3)
                        const int* iqrot_z,
#endif
                        const int& ngiq);
void QUATGRAD_SIDE_SYMM_ISOTROPIC(
    const int& lo0, const int& hi0, const int& lo1, const int& hi1,
#if (NDIM == 3)
    const int& lo2, const int& hi2,
#endif
    const int& depth, const double* h, const double* diff_x,
    const double* diff_y,
#if (NDIM == 3)
    const double* diff_z,
#endif
    const int& ngdiff, const double* grad_x_x, const double* grad_y_x,
#if (NDIM == 3)
    const double* grad_z_x,
#endif
    const double* grad_x_y, const double* grad_y_y,
#if (NDIM == 3)
    const double* grad_z_y, const double* grad_x_z, const double* grad_y_z,
    const double* grad_z_z,
#endif
    const int& nggrad, const int* iqrot_x, const int* iqrot_y,
#if (NDIM == 3)
    const int* iqrot_z,
#endif
    const int& ngiq);

void QUATGRAD_MODULUS(const int&, const int&, const int&, const int&,
#if (NDIM == 3)
                      const int&, const int&,
#endif
                      const int&, const double*, const double*,
#if (NDIM == 3)
                      const double*,
#endif
                      const int&, const double*, const int&);

void QUATGRAD_MODULUS_FROM_SIDES(const int&, const int&, const int&, const int&,
#if (NDIM == 3)
                                 const int&, const int&,
#endif
                                 const int&, const double*, const double*,
#if (NDIM == 3)
                                 const double*,
#endif
                                 const double*, const int&);

void QUATGRAD_MODULUS_FROM_SIDES_COMPACT(const int&, const int&, const int&,
                                         const int&,
#if (NDIM == 3)
                                         const int&, const int&,
#endif
                                         const int&, const double*,
                                         const double*,
#if (NDIM == 3)
                                         const double*,
#endif
                                         const int&, const double*, const int&);

void QUATDIFFUSION(const int&, const int&, const int&, const int&,
#if (NDIM == 3)
                   const int&, const int&,
#endif
                   const double&, const double*, const int&, const double*,
                   const int&, const double*, const double*,
#if (NDIM == 3)
                   const double*,
#endif
                   const int&, const char*, const char*);

void QUATDIFFUSIONDERIV(const int&, const int&, const int&, const int&,
#if (NDIM == 3)
                        const int&, const int&,
#endif
                        const double&, const double*, const int&, const double*,
                        const int&, const int&, const double*, const double*,
#if (NDIM == 3)
                        const double*,
#endif
                        const int&, const double*, const double*,
#if (NDIM == 3)
                        const double*,
#endif
                        const int&, const double&, const char*, const char*,
                        const char*);

void CORRECTRHSQUATFORSYMMETRY(
    const int& lo0, const int& hi0, const int& lo1, const int& hi1,
#if (NDIM == 3)
    const int& lo2, const int& hi2,
#endif
    const int& depth, const double* dx, const double* nonsymm_diff_x,
    const double* nonsymm_diff_y,
#if (NDIM == 3)
    const double* nonsymm_diff_z,
#endif
    const double* symm_diff_x, const double* symm_diff_y,
#if (NDIM == 3)
    const double* symm_diff_z,
#endif
    const int& ngdiff, const double* rhs, const int& ngrhs, const double* quat,
    const int& ngq, const double* diffusion_x, const double* diffusion_y,
#if (NDIM == 3)
    const double* diffusion_z,
#endif
    const int& ngdiffusion, const double* mobility, const int& ngmob,
    const int* iqrot_x, const int* iqrot_y,
#if (NDIM == 3)
    const int* iqrot_z,
#endif
    const int& ngiq);

void CORRECTRHSQUATWITHCONSTRAINT(
    const int&, const int&, const int&, const int&,
#if (NDIM == 3)
    const int&, const int&,
#endif
    const int&, const double*, const double*, const double*, const double*,
#if (NDIM == 3)
    const double*, const double*,
#endif
    const int&, const double*, const int&, const int&, const int&, const int&,
#if (NDIM == 3)
    const int&, const int&,
#endif
    const double*, const int&, const int&, const int&, const int&,
#if (NDIM == 3)
    const int&, const int&,
#endif
    const double*, const double*,
#if (NDIM == 3)
    const double*,
#endif
    const double*, const double*);

void TAGFROMGRADS(const int&, const int&, const int&, const int&,
#if (NDIM == 3)
                  const int&, const int&,
#endif
                  const double*, const double*,
#if (NDIM == 3)
                  const double*,
#endif
                  const int&, const int&, const int&, const int&,
#if (NDIM == 3)
                  const int&, const int&,
#endif
                  const double*, const int*, const int&, const int&, const int&,
                  const int&,
#if (NDIM == 3)
                  const int&, const int&,
#endif
                  const int&, const double&);

void TAGFROMQUATGRADS(const int&, const int&, const int&, const int&,
#if (NDIM == 3)
                      const int&, const int&,
#endif
                      const double*, const double*,
#if (NDIM == 3)
                      const double*,
#endif
                      const int&, const int&, const int&, const int&,
#if (NDIM == 3)
                      const int&, const int&,
#endif
                      const int&, const double*, const int*, const int&,
                      const int&, const int&, const int&,
#if (NDIM == 3)
                      const int&, const int&,
#endif
                      const int&, const double&);

void QUATMOBILITY(const int& ifirst0, const int& ilast0, const int& ifirst1,
                  const int& ilast1,
#if (NDIM == 3)
                  const int& ifirst2, const int& ilast2,
#endif
                  const double* phase, const int& ngphase,
                  const double* mobility, const int& ngmobility,
                  const double& scale_mobility, const double& min_mobility,
                  const char* func_type, const double& alt_scale_factor);

void QUATMOBILITYDERIV(const int& ifirst0, const int& ilast0,
                       const int& ifirst1, const int& ilast1,
#if (NDIM == 3)
                       const int& ifirst2, const int& ilast2,
#endif
                       const double* phase, const int& ngphase,
                       const double* dmobility, const int& ngmobility,
                       const double& scale_mobility, const double& min_mobility,
                       const char* func_type, const double& alt_scale_factor);

void TEMPERATURE_ENERGY(const int& ifirst0, const int& ilast0,
                        const int& ifirst1, const int& ilast1,
#if (NDIM == 3)
                        const int& ifirst2, const int& ilast2,
#endif
                        const double* const temp, const int& ng, double* fs,
                        const double& tm, const double& la);

void QUATENERGY(const int& ifirst0, const int& ilast0, const int& ifirst1,
                const int& ilast1,
#if (NDIM == 3)
                const int& ifirst2, const int& ilast2,
#endif
                const int& depth, const double* const dx,
                const double* const gqx, const double* const gqy,
#if (NDIM == 3)
                const double* const gqz,
#endif
                const int& nggq, const double* const phi, const int& ngphi,
                const int& nphases, const double* const quat, const int& ngq,
                const double& epsilon_phi, const double& epsilon_q,
                const double& anisotropy, const int& knumber,
                const double& misorientation_factor,
                const double* const temperature, const int& ngtemp,
                const double& phi_well_scale, const double* const weight,
                double& total_energy, double& total_interface_e,
                double& total_orient_e, double& total_qint_e,
                const double* const energy, const int& eval_per_cell,
                const char* phi_interp_type, const char* phi_well_type,
                const char* orient_interp_type1,
                const char* orient_interp_type2, const char* avg_type,
                const char* floor_type, const double& floor2);

void PHIPHI_INTERFACIALENERGY(const int& ifirst0, const int& ilast0,
                              const int& ifirst1, const int& ilast1,
#if (NDIM == 3)
                              const int& ifirst2, const int& ilast2,
#endif
                              const double* const phi, const int& ngphi,
                              const int& nphases, const double& gamma,
                              const double& m, const double* const weight,
                              const double* const energy_density,
                              double* energy, const int& eval_per_cell);

void PHIPHI_FORCES(const int& ifirst0, const int& ilast0, const int& ifirst1,
                   const int& ilast1,
#if (NDIM == 3)
                   const int& ifirst2, const int& ilast2,
#endif
                   const double* const dx, const double* const phi,
                   const int& ngphi, const int& nphases, const double& gamma,
                   const double& m, const double* const weight, double* forces);

void WANG_FORCES(const int& ifirst0, const int& ilast0, const int& ifirst1,
                 const int& ilast1,
#if (NDIM == 3)
                 const int& ifirst2, const int& ilast2,
#endif
                 const double* const dx, const double* const phi,
                 const int& ngphi, const int& nphases, const double* const c,
                 const int& ngc, const int& nc, const double& rhoe,
                 const double& cthreshold, const double* const weight,
                 double* forces);

void BULKENERGY(const int& ifirst0, const int& ilast0, const int& ifirst1,
                const int& ilast1,
#if (NDIM == 3)
                const int& ifirst2, const int& ilast2,
#endif
                const double* const phi, const int& ngphi,
                const double* const fl, const double* const fa,
                const double* const weight, double& total_energy,
                double& total_bulk_e, const double* const energy,
                const int& eval_per_cell, const char* phi_interp_type);

void QUATAVG(const int&, const int&, const int&, const int&,
#if (NDIM == 3)
             const int&, const int&,
#endif
             const int&, const double*, const int&, const int&, const int&,
             const int&,
#if (NDIM == 3)
             const int&, const int&,
#endif
             const double*, const int&, const int&, const int&, const int&
#if (NDIM == 3)
             ,
             const int&, const int&
#endif
);

void VELOCITY(const int&, const int&, const int&, const int&,
#if (NDIM == 3)
              const int&, const int&,
#endif
              const double&, const double*, const double*,
#if (NDIM == 3)
              const double*,
#endif
              const double*, const double*);

void ADDRBMOTION(const int& ifirst0, const int& ilast0, const int& ifirst1,
                 const int& ilast1,
#if (NDIM == 3)
                 const int& ifirst2, const int& ilast2,
#endif
                 const double* const dx, const double* const phi,
                 const int& ngphi, const double* const vel, double* rhs,
                 const int& ngrhs);

#if NDIM == 2

// 2d prototypes

void COMPUTE_FACE_COEF2D(const int&, const int&, const int&, const int&,
                         const int&, const double&, const double*, const int&,
                         const double*, const int&, const double&,
                         const double*, const double*, const int&,
                         const double*, const double*, const int&,
                         const double&, const char*, const char*, const char*,
                         const char*);
void COMPUTE_DQUATDPHI_FACE_COEF2D(
    const int&, const int&, const int&, const int&, const int&, const double*,
    const int&, const int&, const int&, const int&, const double*, const int&,
    const int&, const int&, const int&, const double*, const int&, const int&,
    const int&, const int&, const double*, const int&, const int&, const int&,
    const int&, const double*, const int&, const int&, const int&, const int&);
void COMPUTE_FLUX2D(const int&, const int&, const int&, const int&, const int&,
                    const double*, const int&, const int&, const int&,
                    const int&, const double*, const int&, const int&,
                    const int&, const int&, const double*, const int&,
                    const int&, const int&, const int&, const double*,
                    const double*, const int&, const int&, const int&,
                    const int&, const double*, const int&, const int&,
                    const int&, const int&);
void COMPUTE_FLUX2D_FROM_GRADQ(
    const int&, const int&, const int&, const int&, const int&, const double*,
    const int&, const int&, const int&, const int&, const double*, const int&,
    const int&, const int&, const int&,
    const double*,  // const int&, const int&, const int&, const int&,
    const double*,  // const int&, const int&, const int&, const int&,
    const double*, const int&, const int&, const int&, const int&,
    const double*, const int&, const int&, const int&, const int&);
void COMPUTE_SYM_FLUX2D(const int&, const int&, const int&, const int&,
                        const int&, const double*, const int&, const int&,
                        const int&, const int&, const double*, const int&,
                        const int&, const int&, const int&, const double*,
                        const int&, const int&, const int&, const int&,
                        const double*, const int&, const int&, const int&,
                        const int&, const double*, const double*, const int&,
                        const int&, const int&, const int&, const double*,
                        const int&, const int&, const int&, const int&);
void COMPUTE_Q_RESIDUAL2D(const int&, const int&, const int&, const int&,
                          const int&, const double*, const int&, const int&,
                          const int&, const int&, const double*, const int&,
                          const int&, const int&, const int&, const double*,
                          const int&, const int&, const int&, const int&,
                          const double*, const int&, const int&, const int&,
                          const int&, const double*, const double&,
                          const double*, const int&, const int&, const int&,
                          const int&, const double*, const int&, const int&,
                          const int&, const int&);
void COMPUTE_Q_RESIDUAL2D_SYMM(
    const int&, const int&, const int&, const int&, const int&, const double*,
    const int&, const int&, const int&, const int&, const double*, const int&,
    const int&, const int&, const int&, const double*, const int&, const int&,
    const int&, const int&, const double*, const int&, const int&, const int&,
    const int&, const double*, const double&, const double*, const int&,
    const int&, const int&, const int&, const double*, const int&, const int&,
    const int&, const int&, const int*, const int*, const int&);
void ADD_QUAT_OP2D(const int&, const int&, const int&, const int&, const int&,
                   const double*, const int&, const int&, const int&,
                   const int&, const double*, const int&, const int&,
                   const int&, const int&, const double*, const int&,
                   const int&, const int&, const int&, const double*,
                   const double*, const int&, const int&, const int&,
                   const int&);
void ADD_QUAT_PROJ_OP2D(const int&, const int&, const int&, const int&,
                        const int&, const double*, const int&, const int&,
                        const int&, const int&, const double*, const int&,
                        const int&, const int&, const int&, const double*,
                        const int&, const int&, const int&, const int&,
                        const double*, const int&, const int&, const int&,
                        const int&, const double*, const int&, const int&,
                        const int&, const int&, const double*, const double*,
                        const int&, const int&, const int&, const int&);
void ADD_QUAT_PROJ_OP2D_SYMM(const int&, const int&, const int&, const int&,
                             const int&, const double*, const int&, const int&,
                             const int&, const int&, const double*, const int&,
                             const int&, const int&, const int&, const double*,
                             const int&, const int&, const int&, const int&,
                             const double*, const int&, const int&, const int&,
                             const int&, const double*, const int&, const int&,
                             const int&, const int&, const double*,
                             const double*, const int&, const int&, const int&,
                             const int&, const int*, const int*, const int&);
void COMPUTE_LAMBDA_FLUX2D(const int&, const int&, const int&, const int&,
                           const int&, const double*, const int&, const int&,
                           const int&, const int&, const double*, const int&,
                           const int&, const int&, const int&, const double*,
                           const int&, const int&, const int&, const int&,
                           const double*, const double*, const int&, const int&,
                           const int&, const int&);
void COMPUTE_LAMBDA_FLUX2D_SYMM(
    const int&, const int&, const int&, const int&, const int&, const double*,
    const int&, const int&, const int&, const int&, const double*, const int&,
    const int&, const int&, const int&, const double*, const int&, const int&,
    const int&, const int&, const double*, const double*, const int&,
    const int&, const int&, const int&, const int*, const int*, const int&);
void FIXFLUX2D(double*, double*, const int&, const int&, const int&,
               const double*, const double*, const int&, const int&, const int&,
               const double*, const int&, const int&, const int&, const int&,
               const int&, const int&, const int&, const int&, const int&,
               const int*, const int*, const int*, const double*);
void PROJECT2D(const int&, const int&, const int&, const int&, const int&,
               const double*, const int&, const int&, const int&, const int&,
               const double*, const int&, const int&, const int&, const int&,
               const double*, const int&, const int&, const int&, const int&);
void PROJECTPHI2D(const int&, const int&, const int&, const int&, const int&,
                  const double*, const int&, const double*, const int&,
                  const double*, const int&);
void TAKE_SQUARE_ROOT2D(const int&, const int&, const int&, const int&, double*,
                        const int&, const int&, const int&, const int&);
void MULTICOMPONENT_MULTIPLY2D(const int&, const int&, const int&, const int&,
                               const double*, const int&, const int&,
                               const int&, const int&, double*, const int&,
                               const int&, const int&, const int&, const int&);
void MULTICOMPONENT_DIVIDE2D(const int&, const int&, const int&, const int&,
                             const double*, const int&, const int&, const int&,
                             const int&, double*, const int&, const int&,
                             const int&, const int&, const int&);
#endif
#if NDIM == 3

// 3d prototypes

void COMPUTE_FACE_COEF3D(const int&, const int&, const int&, const int&,
                         const int&, const int&, const int&, const double&,
                         const double*, const int&, const double*, const int&,
                         const double&,
                         // grad_q_data
                         const double*, const double*, const double*,
                         const int&,
                         // face_coef_data
                         const double*, const double*, const double*,
                         const int&, const double&, const char*, const char*,
                         const char*, const char*);
void COMPUTE_DQUATDPHI_FACE_COEF3D(
    const int&, const int&, const int&, const int&, const int&, const int&,
    const int&, const double*, const int&, const int&, const int&, const int&,
    const int&, const int&, const double*, const int&, const int&, const int&,
    const int&, const int&, const int&, const double*, const int&, const int&,
    const int&, const int&, const int&, const int&, const double*, const int&,
    const int&, const int&, const int&, const int&, const int&, const double*,
    const int&, const int&, const int&, const int&, const int&, const int&,
    const double*, const int&, const int&, const int&, const int&, const int&,
    const int&, const double*, const int&, const int&, const int&, const int&,
    const int&, const int&);
void COMPUTE_FLUX3D(
    const int&, const int&, const int&, const int&, const int&, const int&,
    const int&, const double*, const int&, const int&, const int&, const int&,
    const int&, const int&, const double*, const int&, const int&, const int&,
    const int&, const int&, const int&, const double*, const int&, const int&,
    const int&, const int&, const int&, const int&, const double*, const int&,
    const int&, const int&, const int&, const int&, const int&, const double*,
    const double*, const int&, const int&, const int&, const int&, const int&,
    const int&, const double*, const int&, const int&, const int&, const int&,
    const int&, const int&, const double*, const int&, const int&, const int&,
    const int&, const int&, const int&);
void COMPUTE_FLUX3D_FROM_GRADQ(
    const int&, const int&, const int&, const int&, const int&, const int&,
    const int&, const double*, const int&, const int&, const int&, const int&,
    const int&, const int&, const double*, const int&, const int&, const int&,
    const int&, const int&, const int&, const double*, const int&, const int&,
    const int&, const int&, const int&, const int&, const double*,
    const double*, const double*, const double*, const int&, const int&,
    const int&, const int&, const int&, const int&, const double*, const int&,
    const int&, const int&, const int&, const int&, const int&, const double*,
    const int&, const int&, const int&, const int&, const int&, const int&);
void COMPUTE_SYM_FLUX3D(
    const int&, const int&, const int&, const int&, const int&, const int&,
    const int&, const double*, const int&, const int&, const int&, const int&,
    const int&, const int&, const double*, const int&, const int&, const int&,
    const int&, const int&, const int&, const double*, const int&, const int&,
    const int&, const int&, const int&, const int&, const double*, const int&,
    const int&, const int&, const int&, const int&, const int&, const double*,
    const int&, const int&, const int&, const int&, const int&, const int&,
    const double*, const double*, const int&, const int&, const int&,
    const int&, const int&, const int&, const double*, const int&, const int&,
    const int&, const int&, const int&, const int&, const double*, const int&,
    const int&, const int&, const int&, const int&, const int&);
void COMPUTE_Q_RESIDUAL3D(
    const int&, const int&, const int&, const int&, const int&, const int&,
    const int&, const double*, const int&, const int&, const int&, const int&,
    const int&, const int&, const double*, const int&, const int&, const int&,
    const int&, const int&, const int&, const double*, const int&, const int&,
    const int&, const int&, const int&, const int&, const double*, const int&,
    const int&, const int&, const int&, const int&, const int&, const double*,
    const int&, const int&, const int&, const int&, const int&, const int&,
    const double*, const double&, const double*, const int&, const int&,
    const int&, const int&, const int&, const int&, const double*, const int&,
    const int&, const int&, const int&, const int&, const int&);
void COMPUTE_Q_RESIDUAL3D_SYMM(
    const int&, const int&, const int&, const int&, const int&, const int&,
    const int&, const double*, const int&, const int&, const int&, const int&,
    const int&, const int&, const double*, const int&, const int&, const int&,
    const int&, const int&, const int&, const double*, const int&, const int&,
    const int&, const int&, const int&, const int&, const double*, const int&,
    const int&, const int&, const int&, const int&, const int&, const double*,
    const int&, const int&, const int&, const int&, const int&, const int&,
    const double*, const double&, const double*, const int&, const int&,
    const int&, const int&, const int&, const int&, const double*, const int&,
    const int&, const int&, const int&, const int&, const int&, const int*,
    const int*, const int*, const int&);
void ADD_QUAT_OP3D(const int&, const int&, const int&, const int&, const int&,
                   const int&, const int&, const double*, const int&,
                   const int&, const int&, const int&, const int&, const int&,
                   const double*, const int&, const int&, const int&,
                   const int&, const int&, const int&, const double*,
                   const int&, const int&, const int&, const int&, const int&,
                   const int&, const double*, const int&, const int&,
                   const int&, const int&, const int&, const int&,
                   const double*, const double*, const int&, const int&,
                   const int&, const int&, const int&, const int&);
void ADD_QUAT_PROJ_OP3D(
    const int&, const int&, const int&, const int&, const int&, const int&,
    const int&, const double*, const int&, const int&, const int&, const int&,
    const int&, const int&, const double*, const int&, const int&, const int&,
    const int&, const int&, const int&, const double*, const int&, const int&,
    const int&, const int&, const int&, const int&, const double*, const int&,
    const int&, const int&, const int&, const int&, const int&, const double*,
    const int&, const int&, const int&, const int&, const int&, const int&,
    const double*, const int&, const int&, const int&, const int&, const int&,
    const int&, const double*, const double*, const int&, const int&,
    const int&, const int&, const int&, const int&);
void COMPUTE_LAMBDA_FLUX3D(const int&, const int&, const int&, const int&,
                           const int&, const int&, const int&, const double*,
                           const int&, const int&, const int&, const int&,
                           const int&, const int&, const double*, const int&,
                           const int&, const int&, const int&, const int&,
                           const int&, const double*, const int&, const int&,
                           const int&, const int&, const int&, const int&,
                           const double*, const int&, const int&, const int&,
                           const int&, const int&, const int&, const double*,
                           const double*, const int&, const int&, const int&,
                           const int&, const int&, const int&);
void FIXFLUX3D(double*, double*, double*, const int&, const int&, const int&,
               const int&, const double*, const double*, const double*,
               const int&, const int&, const int&, const int&, const double*,
               const int&, const int&, const int&, const int&, const int&,
               const int&, const int&, const int&, const int&, const int&,
               const int&, const int&, const int*, const int*, const int*,
               const double*);
void PROJECT3D(const int&, const int&, const int&, const int&, const int&,
               const int&, const int&, const double*, const int&, const int&,
               const int&, const int&, const int&, const int&, const double*,
               const int&, const int&, const int&, const int&, const int&,
               const int&, const double*, const int&, const int&, const int&,
               const int&, const int&, const int&);
void PROJECTPHI3D(const int&, const int&, const int&, const int&, const int&,
                  const int&, const int&, const double*, const int&,
                  const double*, const int&, const double*, const int&);
void TAKE_SQUARE_ROOT3D(const int&, const int&, const int&, const int&,
                        const int&, const int&, double*, const int&, const int&,
                        const int&, const int&, const int&, const int&);
void MULTICOMPONENT_MULTIPLY3D(const int&, const int&, const int&, const int&,
                               const int&, const int&, const double*,
                               const int&, const int&, const int&, const int&,
                               const int&, const int&, double*, const int&,
                               const int&, const int&, const int&, const int&,
                               const int&, const int&);
void MULTICOMPONENT_DIVIDE3D(const int&, const int&, const int&, const int&,
                             const int&, const int&, const double*, const int&,
                             const int&, const int&, const int&, const int&,
                             const int&, double*, const int&, const int&,
                             const int&, const int&, const int&, const int&,
                             const int&);

#endif
}
