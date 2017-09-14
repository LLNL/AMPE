// Copyright (c) 2018, Lawrence Livermore National Security, LLC.
// Produced at the Lawrence Livermore National Laboratory
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
// LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
// OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
// 
// Link between C/C++ and Fortran files
//       name in             name in
//      C/C++ code            Fortran code
//      ----------            ------------
#define FORT_GRADIENT_FLUX gradient_flux_ 
#define FORT_COMPUTE_FLUX_ISOTROPIC compute_flux_isotropic_ 
#define FORT_ANISOTROPIC_GRADIENT_FLUX anisotropic_gradient_flux_ 
#define FORT_COMP_RHS_PBG computerhspbg_
#define FORT_COMP_RHS_TEMP computerhstemp_
#define FORT_COMP_RHS_BIASWELL computerhsbiaswell_
#define FORT_COMP_RHS_BIASWELL_BECKERMANN computerhsbiaswellbeckermann_
#define FORT_COMP_RHS_ETA computerhseta_
#define FORT_DIFFS diffs_
#define FORT_GRAD_CELL grad_cell_
#define FORT_GRAD_SIDE grad_side_
#define FORT_GRAD_SIDE_ISOTROPIC grad_side_isotropic_
#define FORT_QUATSYMMROTATION quat_symm_rotation_
#define FORT_QUATFUNDAMENTAL quat_fundamental_
#define FORT_QUATDIFFS quatdiffs_
#define FORT_QUATDIFFS_SYMM quatdiffs_symm_
#define FORT_QUATGRAD_CELL quatgrad_cell_
#define FORT_QUATGRAD_CELL_SYMM quatgrad_cell_symm_
#define FORT_QUATGRAD_SIDE quatgrad_side_
#define FORT_QUATGRAD_SIDE_SYMM quatgrad_side_symm_
#define FORT_QUATGRAD_SIDE_SYMM_ISOTROPIC quatgrad_side_symm_isotropic_
#define FORT_QUATGRAD_SIDE_ISOTROPIC quatgrad_side_isotropic_
#define FORT_QUATGRAD_MODULUS quatgrad_modulus_
#define FORT_QUATGRAD_MODULUS_FROM_SIDES quatgrad_modulus_from_sides_
#define FORT_QUATGRAD_MODULUS_FROM_SIDES_COMPACT quatgrad_modulus_from_sides_compact_
#define FORT_QUATDIFFUSION quatdiffusion_
#define FORT_QUATDIFFUSIONDERIV quatdiffusionderiv_
#define FORT_CORRECTRHSQUAT correctrhsquatforsymmetry_
#define FORT_TAG_GRAD tagfromgrads_
#define FORT_TAG_QUAT_GRAD tagfromquatgrads_
#define FORT_QUATMOBILITY quatmobility_
#define FORT_QUATMOBILITYDERIV quatmobilityderiv_
#define FORT_QUATENERGY quatenergy_
#define FORT_QUATAVG quatavg_
#define FORT_SMOOTHQUAT smoothquat_
#define FORT_SETTOZERO settozero_
#define FORT_SOURCE_TEMPERATURE source_temperature_
#define FORT_HEAT_CAPACITY_NKR heat_capacity_nkr_
#define FORT_VELOCITY velocity_

// Function argument list interfaces
extern "C" {
   void FORT_SETTOZERO(const int&, const int&,
                  const int&, const int&,
                  const int&, const int&,
                  const int&, const int&,
#if (NDIM == 3)
                  const int&, const int&,
                  const int&, const int&,
#endif
                  const double*);
   void FORT_SOURCE_TEMPERATURE(const int&, const int&,
                  const int&, const int&,
#if (NDIM == 3)
                  const int&, const int&,
#endif
                  double* conc, const int&,
                  double* rhs, const int&,
                  const double* const cp,  const int&,
                  const double* const s,
                  const int&);

   void FORT_HEAT_CAPACITY_NKR(const int&, const int&,
                  const int&, const int&,
#if (NDIM == 3)
                  const int&, const int&,
#endif
                  const double* const conc, const int&,
                  const double* const temp, const int&,
                  double* cp,  const int&,
                  const int* const cp_powers, const int&,
                  const double* const cpcoeffs, const int&);

   void FORT_GRADIENT_FLUX(
      const int& ifirst0, const int& ilast0,
      const int& ifirst1, const int& ilast1,
#if (NDIM == 3)
      const int& ifirst2, const int& ilast2,
#endif
      const double* dx,
      const double& epsilon,
      double* phase, const int& phaseng,
      double* flux0,
      double* flux1,
#if (NDIM == 3)
      double* flux2,
#endif
      const int& fluxng);

   void FORT_COMPUTE_FLUX_ISOTROPIC(
      const int& ifirst0, const int& ilast0,
      const int& ifirst1, const int& ilast1,
#if (NDIM == 3)
      const int& ifirst2, const int& ilast2,
#endif
      const double* dx,
      const double& epsilon,
      double* phase, const int& phaseng,
      double* flux0,
      double* flux1,
#if (NDIM == 3)
      double* flux2,
#endif
      const int& fluxng);

   void FORT_ANISOTROPIC_GRADIENT_FLUX(
      const int& ifirst0, const int& ilast0,
      const int& ifirst1, const int& ilast1,
#if (NDIM == 3)
      const int& ifirst2, const int& ilast2,
#endif
      const double* dx,
      const double& epsilon,
      const double& nu,
      const int& knumber,
      double* phase, const int& phaseng,
      double* quat, const int& qng, const int& qlen,
      double* flux0,
      double* flux1,
#if (NDIM == 3)
      double* flux2,
#endif
      const int& fluxng);

   void FORT_COMP_RHS_PBG(
      const int& ifirst0, const int& ilast0,
      const int& ifirst1, const int& ilast1,
#if (NDIM == 3)
      const int& ifirst2, const int& ilast2,
#endif
      const double* dx,
      const double& misorientation_factor,
      const double* phase_flux0, const double* phase_flux1,
#if (NDIM == 3)
      const double* phase_flux2,
#endif
      const int& ngphaseflux,
      const double* temperature, const int& ngtemp,
      const double& phi_well_scale,
      const double& eta_well_scale,
      const double* phi, const int& ngphi,
      const double* eta, const int& ngeta,
      const double* orient_grad_mod, const int& ngogm,
      double* rhs, const int& ngrhs,
      const char* phi_well_type,
      const char* eta_well_type,
      const char* phi_interp_type,
      const char* orient_interp_type,
      const int& with_orient,
      const int& three_phase
      );

   void FORT_COMP_RHS_TEMP(
      const int& ifirst0, const int& ilast0,
      const int& ifirst1, const int& ilast1,
#if (NDIM == 3)
      const int& ifirst2, const int& ilast2,
#endif
      const double* dx,
      const double& thermal_diffusivity,
      const double& latent_heat,
      const double* temperature, const int& ngtemp,
      const double* cpl, const int& ngcpl,
      const double* phirhs, const int& ngphirhs,
      double* rhs, const int& ngrhs
      );

   void FORT_COMP_RHS_BIASWELL(
      const int& ifirst0, const int& ilast0,
      const int& ifirst1, const int& ilast1,
#if (NDIM == 3)
      const int& ifirst2, const int& ilast2,
#endif
      const double*, const int&,
      const double*, const int&,
      const double&,const double&,
      double* teq, const int&,
      double* rhs, const int&
      );

   void FORT_COMP_RHS_BIASWELL_BECKERMANN(
      const int& ifirst0, const int& ilast0,
      const int& ifirst1, const int& ilast1,
#if (NDIM == 3)
      const int& ifirst2, const int& ilast2,
#endif
      const double*, const int&,
      const double*, const int&,
      const double&,
      double* teq, const int&,
      double* rhs, const int&
      );

   void FORT_COMP_RHS_ETA(
      const int& ifirst0, const int& ilast0,
      const int& ifirst1, const int& ilast1,
#if (NDIM == 3)
      const int& ifirst2, const int& ilast2,
#endif
      const double* dx,
      const double* eta_flux0, const double* eta_flux1,
#if (NDIM == 3)
      const double* eta_flux2,
#endif
      const int& ngetaflux,
      const double* temperature, const int& ngtemp,
      const double& eta_well_scale,
      const double* phi, const int& ngphi,
      const double* eta, const int& ngeta,
      double* rhs, const int& ngrhs,
      const char* eta_well_type,
      const char* phi_interp_type
      );

   void FORT_SMOOTHQUAT(
      const int& ifirst0, const int& ilast0,
      const int& ifirst1, const int& ilast1,
#if (NDIM == 3)
      const int& ifirst2, const int& ilast2,
#endif
      const int& depth,
      const double*, const int&,
      const double*, const int&,
      const double*, const double*,
#if (NDIM == 3)
      const double*,
#endif
      const int&,
      const double*, const int&,
      const double& );

   void FORT_DIFFS(
      const int&, const int&,
      const int&, const int&,
#if (NDIM == 3)
      const int&, const int&,
#endif
      const double*,
      const int&, const int&,
      const int&, const int&,
#if (NDIM == 3)
      const int&, const int&,
#endif
      const double*,
      const double*,
#if (NDIM == 3)
      const double*,
#endif
      const int&, const int&,
      const int&, const int&
#if (NDIM == 3)
      ,const int&, const int&
#endif
      );

   void FORT_GRAD_CELL(
      const int&, const int&,
      const int&, const int&,
#if (NDIM == 3)
      const int&, const int&,
#endif
      const double*,
      const double*,
#if (NDIM == 3)
      const double*,
#endif
      const int&, const int&,
      const int&, const int&,
#if (NDIM == 3)
      const int&, const int&,
#endif
      const double*,
      const double*,
      const double*,
#if (NDIM == 3)
      const double*,
#endif
      const int&, const int&,
      const int&, const int&
#if (NDIM == 3)
      , const int&, const int&
#endif
      );

   void FORT_GRAD_SIDE(
      const int&, const int&,
      const int&, const int&,
#if (NDIM == 3)
      const int&, const int&,
#endif
      const double*,
      const double*,
#if (NDIM == 3)
      const double*,
#endif
      const int&, const int&,
      const int&, const int&,
#if (NDIM == 3)
      const int&, const int&,
#endif
      const double*,
      const double*, const double*,
#if (NDIM == 3)
      const double*,
#endif
      const double*, const double*,
#if (NDIM == 3)
      const double*,
      const double*, const double*, const double*,
#endif
      const int&, const int&,
      const int&, const int&
#if (NDIM == 3)
      , const int&, const int&
#endif
      );

   void FORT_GRAD_SIDE_ISOTROPIC(
      const int&, const int&,
      const int&, const int&,
#if (NDIM == 3)
      const int&, const int&,
#endif
      const double*,
      const double*,
#if (NDIM == 3)
      const double*,
#endif
      const int&, const int&,
      const int&, const int&,
#if (NDIM == 3)
      const int&, const int&,
#endif
      const double*,
      const double*, const double*,
#if (NDIM == 3)
      const double*,
#endif
      const double*, const double*,
#if (NDIM == 3)
      const double*,
      const double*, const double*, const double*,
#endif
      const int&, const int&,
      const int&, const int&
#if (NDIM == 3)
      , const int&, const int&
#endif
      );

   void FORT_QUATSYMMROTATION(
      const int& lo0, const int& hi0,
      const int& lo1, const int& hi1,
#if (NDIM == 3)
      const int& lo2, const int& hi2,
#endif
      const double* q, const int& ngq,
      const int& depth,
      const int* rot_x,
      const int* rot_y,
#if (NDIM == 3)
      const int* rot_z,
#endif
      const int& ngrot
      );

   void FORT_QUATFUNDAMENTAL(
      const int&, const int&,
      const int&, const int&,
#if (NDIM == 3)
      const int&, const int&,
#endif
      const double *,
      const int &, const int &,
      const int &, const int &,
#if (NDIM == 3)
      const int &, const int &,
#endif
      const int &
      );

   void FORT_QUATDIFFS(
      const int& lo0, const int& hi0,
      const int& lo1, const int& hi1,
#if (NDIM == 3)
      const int& lo2, const int& hi2,
#endif
      const int& depth,
      const double* q,
      const int& ngq,
      const double* diff_x,
      const double* diff_y,
#if (NDIM == 3)
      const double* diff_z,
#endif
      const int& ngdiff
      );

   void FORT_QUATDIFFS_SYMM(
      const int& lo0, const int& hi0,
      const int& lo1, const int& hi1,
#if (NDIM == 3)
      const int& lo2, const int& hi2,
#endif
      const int& depth,
      const double* q,
      const int& ngq,
      const double* diff_x,
      const double* diff_y,
#if (NDIM == 3)
      const double* diff_z,
#endif
      const int& ngdiff,
      const int* iqrot_x,
      const int* iqrot_y,
#if (NDIM == 3)
      const int* iqrot_z,
#endif
      const int& ngiq
      );

   void FORT_QUATGRAD_CELL(
      const int& lo0, const int& hi0,
      const int& lo1, const int& hi1,
#if (NDIM == 3)
      const int& lo2, const int& hi2,
#endif
      const int& depth,
      const double* h,
      const double* diff_x,
      const double* diff_y,
#if (NDIM == 3)
      const double* diff_z,
#endif
      const int& ngdiff,
      const double* grad_x,
      const double* grad_y,
#if (NDIM == 3)
      const double* grad_z,
#endif
      const int& nggrad
      );

   void FORT_QUATGRAD_CELL_SYMM(
      const int& lo0, const int& hi0,
      const int& lo1, const int& hi1,
#if (NDIM == 3)
      const int& lo2, const int& hi2,
#endif
      const int& depth,
      const double* h,
      const double* diff_x,
      const double* diff_y,
#if (NDIM == 3)
      const double* diff_z,
#endif
      const int& ngdiff,
      const double* grad_x,
      const double* grad_y,
#if (NDIM == 3)
      const double* grad_z,
#endif
      const int& nggrad,
      const int* iqrot_x,
      const int* iqrot_y,
#if (NDIM == 3)
      const int* iqrot_z,
#endif
      const int& ngiq
      );

   void FORT_QUATGRAD_SIDE(
      const int& lo0, const int& hi0,
      const int& lo1, const int& hi1,
#if (NDIM == 3)
      const int& lo2, const int& hi2,
#endif
      const int& depth,
      const double* h,
      const double* diff_x,
      const double* diff_y,
#if (NDIM == 3)
      const double* diff_z,
#endif
      const int& ngdiff,
      const double* grad_x_x,
      const double* grad_y_x,
#if (NDIM == 3)
      const double* grad_z_x,
#endif
      const double* grad_x_y,
      const double* grad_y_y,
#if (NDIM == 3)
      const double* grad_z_y,
      const double* grad_x_z,
      const double* grad_y_z,
      const double* grad_z_z,
#endif
      const int& nggrad
      );

   void FORT_QUATGRAD_SIDE_ISOTROPIC(
      const int& lo0, const int& hi0,
      const int& lo1, const int& hi1,
#if (NDIM == 3)
      const int& lo2, const int& hi2,
#endif
      const int& depth,
      const double* h,
      const double* diff_x,
      const double* diff_y,
#if (NDIM == 3)
      const double* diff_z,
#endif
      const int& ngdiff,
      const double* grad_x_x,
      const double* grad_y_x,
#if (NDIM == 3)
      const double* grad_z_x,
#endif
      const double* grad_x_y,
      const double* grad_y_y,
#if (NDIM == 3)
      const double* grad_z_y,
      const double* grad_x_z,
      const double* grad_y_z,
      const double* grad_z_z,
#endif
      const int& nggrad
      );

   void FORT_QUATGRAD_SIDE_SYMM(
      const int& lo0, const int& hi0,
      const int& lo1, const int& hi1,
#if (NDIM == 3)
      const int& lo2, const int& hi2,
#endif
      const int& depth,
      const double* h,
      const double* diff_x,
      const double* diff_y,
#if (NDIM == 3)
      const double* diff_z,
#endif
      const int& ngdiff,
      const double* grad_x_x,
      const double* grad_y_x,
#if (NDIM == 3)
      const double* grad_z_x,
#endif
      const double* grad_x_y,
      const double* grad_y_y,
#if (NDIM == 3)
      const double* grad_z_y,
      const double* grad_x_z,
      const double* grad_y_z,
      const double* grad_z_z,
#endif
      const int& nggrad,
      const int* iqrot_x,
      const int* iqrot_y,
#if (NDIM == 3)
      const int* iqrot_z,
#endif
      const int& ngiq
      );
   void FORT_QUATGRAD_SIDE_SYMM_ISOTROPIC(
      const int& lo0, const int& hi0,
      const int& lo1, const int& hi1,
#if (NDIM == 3)
      const int& lo2, const int& hi2,
#endif
      const int& depth,
      const double* h,
      const double* diff_x,
      const double* diff_y,
#if (NDIM == 3)
      const double* diff_z,
#endif
      const int& ngdiff,
      const double* grad_x_x,
      const double* grad_y_x,
#if (NDIM == 3)
      const double* grad_z_x,
#endif
      const double* grad_x_y,
      const double* grad_y_y,
#if (NDIM == 3)
      const double* grad_z_y,
      const double* grad_x_z,
      const double* grad_y_z,
      const double* grad_z_z,
#endif
      const int& nggrad,
      const int* iqrot_x,
      const int* iqrot_y,
#if (NDIM == 3)
      const int* iqrot_z,
#endif
      const int& ngiq
      );

   void FORT_QUATGRAD_MODULUS(
      const int&, const int&,
      const int&, const int&,
#if (NDIM == 3)
      const int&, const int&,
#endif
      const int&,
      const double*,
      const double*,
#if (NDIM == 3)
      const double*,
#endif
      const int&, const int&,
      const int&, const int&,
#if (NDIM == 3)
      const int&, const int&,
#endif
      const double*,
      const int&, const int&,
      const int&, const int&,
#if (NDIM == 3)
      const int&, const int&,
#endif
      const char*,
      const double&
      );

   void FORT_QUATGRAD_MODULUS_FROM_SIDES(
      const int&, const int&,
      const int&, const int&,
#if (NDIM == 3)
      const int&, const int&,
#endif
      const int&,
      const double*,
      const double*,
#if (NDIM == 3)
      const double*,
#endif
      const int&, const int&,
      const int&, const int&,
#if (NDIM == 3)
      const int&, const int&,
#endif
      const double*,
      const int&, const int&,
      const int&, const int&,
#if (NDIM == 3)
      const int&, const int&,
#endif
      const char*,
      const double&
      );

   void FORT_QUATGRAD_MODULUS_FROM_SIDES_COMPACT(
      const int&, const int&,
      const int&, const int&,
#if (NDIM == 3)
      const int&, const int&,
#endif
      const int&,
      const double*,
      const double*,
#if (NDIM == 3)
      const double*,
#endif
      const int&, const int&,
      const int&, const int&,
#if (NDIM == 3)
      const int&, const int&,
#endif
      const double*,
      const int&, const int&,
      const int&, const int&,
#if (NDIM == 3)
      const int&, const int&,
#endif
      const char*,
      const double&
      );

   void FORT_QUATDIFFUSION(
      const int& , const int& ,
      const int& , const int& ,
#if (NDIM == 3)
      const int& , const int& ,
#endif
      const double&,
      const double*, const int&,
      const double*, const int&,
      const double*, const double*,
#if (NDIM == 3)
      const double*,
#endif
      const int&,
      const char*, const char*);

   void FORT_QUATDIFFUSIONDERIV(
      const int& , const int& ,
      const int& , const int& ,
#if (NDIM == 3)
      const int& , const int& ,
#endif
      const double&,
      const double*, const int&,
      const double*, const int&,
      const int&,
      const double*, const double*,
#if (NDIM == 3)
      const double*,
#endif
      const int&,
      const double*, const double*,
#if (NDIM == 3)
      const double*,
#endif
      const int&,
      const double&,
      const char*, const char*, const char* );

   void FORT_CORRECTRHSQUAT(
      const int& lo0, const int& hi0,
      const int& lo1, const int& hi1,
#if (NDIM == 3)
      const int& lo2, const int& hi2,
#endif
      const int& depth,
      const double* dx,
      const double* nonsymm_diff_x,
      const double* nonsymm_diff_y,
#if (NDIM == 3)
      const double* nonsymm_diff_z,
#endif
      const double* symm_diff_x,
      const double* symm_diff_y,
#if (NDIM == 3)
      const double* symm_diff_z,
#endif
      const int& ngdiff,
      const double* rhs, const int& ngrhs,
      const double* quat, const int& ngq,
      const double* diffusion_x,
      const double* diffusion_y,
#if (NDIM == 3)
      const double* diffusion_z,
#endif
      const int& ngdiffusion ,
      const double* mobility, const int& ngmob,
      const int* iqrot_x,
      const int* iqrot_y,
#if (NDIM == 3)
      const int* iqrot_z,
#endif
      const int& ngiq );

   void FORT_CORRECTRHSQUATWITHCONSTRAINT(
      const int&, const int&,
      const int&, const int&,
#if (NDIM == 3)
      const int& , const int& ,
#endif
      const int&,
      const double *, const double *,
      const double *, const double *,
#if (NDIM == 3)
      const double*, const double *,
#endif
      const int&,
      const double *,
      const int&, const int&,
      const int&, const int&,
#if (NDIM == 3)
      const int& , const int& ,
#endif
      const double *,
      const int&, const int&,
      const int&, const int&,
#if (NDIM == 3)
      const int& , const int& ,
#endif
      const double *, const double *,
#if (NDIM == 3)
      const double*,
#endif
      const double *,
      const double * );

   void FORT_TAG_GRAD(
      const int&, const int&, const int&, const int&,
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
      const double*,
      const int*,
      const int&, const int&, const int&, const int&,
#if (NDIM == 3)
      const int&, const int&,
#endif
      const int&,
      const double& );

   void FORT_TAG_QUAT_GRAD(
      const int&, const int&, const int&, const int&,
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
      const int&,
      const double*,
      const int*,
      const int&, const int&, const int&, const int&,
#if (NDIM == 3)
      const int&, const int&,
#endif
      const int&,
      const double& );

   void FORT_QUATMOBILITY(
      const int& ifirst0, const int& ilast0,
      const int& ifirst1, const int& ilast1,
#if (NDIM == 3)
      const int& ifirst2, const int& ilast2,
#endif
      const double* phase,
      const int& ngphase,
      const double* mobility,
      const int& ngmobility,
      const double& scale_mobility,
      const double& min_mobility,
      const char* func_type,
      const double& alt_scale_factor );

   void FORT_QUATMOBILITYDERIV(
      const int& ifirst0, const int& ilast0,
      const int& ifirst1, const int& ilast1,
#if (NDIM == 3)
      const int& ifirst2, const int& ilast2,
#endif
      const double* phase,
      const int& ngphase,
      const double* dmobility,
      const int& ngmobility,
      const double& scale_mobility,
      const double& min_mobility,
      const char* func_type,
      const double& alt_scale_factor );

   void FORT_QUATENERGY(
      const int& ifirst0, const int& ilast0,
      const int& ifirst1, const int& ilast1,
#if (NDIM == 3)
      const int& ifirst2, const int& ilast2,
#endif
      const int& depth,
      const double* const dx,
      const double* const gqx,
      const double* const gqy,
#if (NDIM == 3)
      const double* const gqz,
#endif
      const int& nggq,
      const double* const phi, const int& ngphi,
      const double* const eta, const int& ngeta,
      const double& epsilon_phi,
      const double& epsilon_eta,
      const double& epsilon_q,
      const double& misorientation_factor,
      const double* const temperature, const int& ngtemp,
      const double& phi_well_scale, 
      const double& eta_well_scale, 
      const double* const fl,
      const double* const fa,
      const double* const fb,
      const int& three_phase,
      const double* const weight,
      double& total_energy,
      double& total_phi_e,
      double& total_eta_e,
      double& total_orient_e,
      double& total_qint_e,
      double& total_well_e,
      double& total_free_e,
      const double* const energy,
      const int& eval_per_cell,
      const char* phi_interp_type,
      const char* eta_interp_type,
      const char* phi_well_type,
      const char* eta_well_type,
      const char* orient_interp_type,
      const char* avg_type,
      const char* floor_type,
      const double& floor2 );

   void FORT_QUATAVG(
      const int&, const int&,
      const int&, const int&,
#if (NDIM == 3)
      const int&, const int&,
#endif
      const int&,
      const double *,
      const int&, const int&,
      const int&, const int&,
#if (NDIM == 3)
      const int&, const int&,
#endif
      const double *,
      const int&, const int&,
      const int&, const int&
#if (NDIM == 3)
      , const int&, const int&
#endif
      );

   void FORT_VELOCITY(
      const int&, const int&,
      const int&, const int&,
#if (NDIM == 3)
      const int&, const int&,
#endif
      const double&,
      const double*,
      const double*,
#if (NDIM == 3)
      const double*,
#endif
      const double*,
      const double*
      );

#if NDIM==2

  // 2d prototypes

  void compute_face_coef2d_(const int&, const int&, const int&, const int&,
                            const int&,
                            const double &,
                            const double *, const int&, const int&, const int&, const int&,
                            const double *, const int&, const int&, const int&, const int&,
                            const double *, const int&, const int&, const int&, const int&,
                            const double *, const int&, const int&, const int&, const int&,
                            const double *, const int&, const int&, const int&, const int&,
                            const double *, const int&, const int&, const int&, const int&,
                            const double &, const char*);
  void compute_dquatdphi_face_coef2d_(const int&, const int&, const int&, const int&,
                                      const int&,
                                      const double *, const int&, const int&, const int&, const int&,
                                      const double *, const int&, const int&, const int&, const int&,
                                      const double *, const int&, const int&, const int&, const int&,
                                      const double *, const int&, const int&, const int&, const int&,
                                      const double *, const int&, const int&, const int&, const int&);
  void compute_flux2d_(const int&, const int&, const int&, const int&,
                       const int&,
                       const double *, const int&, const int&, const int&, const int&,
                       const double *, const int&, const int&, const int&, const int&,
                       const double *, const int&, const int&, const int&, const int&,
                       const double *,
                       const double *, const int&, const int&, const int&, const int&,
                       const double *, const int&, const int&, const int&, const int&);
  void compute_flux2d_from_gradq_(const int&, const int&, const int&, const int&,
                       const int&,
                       const double *, const int&, const int&, const int&, const int&,
                       const double *, const int&, const int&, const int&, const int&,
                       const double *, //const int&, const int&, const int&, const int&,
                       const double *, //const int&, const int&, const int&, const int&,
                       const double *, const int&, const int&, const int&, const int&,
                       const double *, const int&, const int&, const int&, const int&);
  void compute_sym_flux2d_(const int&, const int&, const int&, const int&,
                           const int&,
                           const double *, const int&, const int&, const int&, const int&,
                           const double *, const int&, const int&, const int&, const int&,
                           const double *, const int&, const int&, const int&, const int&,
                           const double *, const int&, const int&, const int&, const int&,
                           const double *,
                           const double *, const int&, const int&, const int&, const int&,
                           const double *, const int&, const int&, const int&, const int&);
  void compute_q_residual2d_(const int&, const int&, const int&, const int&,
                             const int&,
                             const double *, const int&, const int&, const int&, const int&,
                             const double *, const int&, const int&, const int&, const int&,
                             const double *, const int&, const int&, const int&, const int&,
                             const double *, const int&, const int&, const int&, const int&,
                             const double *, const double &,
                             const double *, const int&, const int&, const int&, const int&,
                             const double *, const int&, const int&, const int&, const int&);
  void compute_q_residual2d_symm_(const int&, const int&, const int&, const int&,
                             const int&,
                             const double *, const int&, const int&, const int&, const int&,
                             const double *, const int&, const int&, const int&, const int&,
                             const double *, const int&, const int&, const int&, const int&,
                             const double *, const int&, const int&, const int&, const int&,
                             const double *, const double &,
                             const double *, const int&, const int&, const int&, const int&,
                             const double *, const int&, const int&, const int&, const int&,
                             const int*, const int*, const int&);
  void add_quat_op2d_(const int&, const int&, const int&, const int&,
                      const int&,
                      const double *, const int&, const int&, const int&, const int&,
                      const double *, const int&, const int&, const int&, const int&,
                      const double *, const int&, const int&, const int&, const int&,
                      const double *,
                      const double *, const int&, const int&, const int&, const int&);
  void add_quat_proj_op2d_(const int&, const int&, const int&, const int&,
                           const int&,
                           const double *, const int&, const int&, const int&, const int&,
                           const double *, const int&, const int&, const int&, const int&,
                           const double *, const int&, const int&, const int&, const int&,
                           const double *, const int&, const int&, const int&, const int&,
                           const double *, const int&, const int&, const int&, const int&,
                           const double *,
                           const double *, const int&, const int&, const int&, const int&);
  void add_quat_proj_op2d_symm_(const int&, const int&, const int&, const int&,
                           const int&,
                           const double *, const int&, const int&, const int&, const int&,
                           const double *, const int&, const int&, const int&, const int&,
                           const double *, const int&, const int&, const int&, const int&,
                           const double *, const int&, const int&, const int&, const int&,
                           const double *, const int&, const int&, const int&, const int&,
                           const double *,
                           const double *, const int&, const int&, const int&, const int&,
                           const int*, const int*, const int&);
  void compute_lambda_flux2d_(const int&, const int&, const int&, const int&,
                              const int&,
                              const double *, const int&, const int&, const int&, const int&,
                              const double *, const int&, const int&, const int&, const int&,
                              const double *, const int&, const int&, const int&, const int&,
                              const double *,
                              const double *, const int&, const int&, const int&, const int&);
  void compute_lambda_flux2d_symm_(const int&, const int&, const int&, const int&,
                              const int&,
                              const double *, const int&, const int&, const int&, const int&,
                              const double *, const int&, const int&, const int&, const int&,
                              const double *, const int&, const int&, const int&, const int&,
                              const double *,
                              const double *, const int&, const int&, const int&, const int&,
                              const int*, const int*, const int&);
  void fixflux2d_(double *, double *, const int &, const int &, const int &,
                  const double *, const double *, const int &, const int &, const int &,
                  const double *, const int &, const int &, const int &,
                  const int &, const int &, const int &, const int &, const int &,
                  const int &,
                  const int *,
                  const int *, const int *,
                  const double *);
  void project2d_(const int&, const int&, const int&, const int&,
                  const int&,
                  const double *, const int&, const int&, const int&, const int&,
                  const double *, const int&, const int&, const int&, const int&,
                  const double *, const int&, const int&, const int&, const int&);
  void take_square_root2d_(const int&, const int&, const int&, const int&,
                           double *, const int&, const int&, const int&, const int&);
  void multicomponent_multiply2d_(const int&, const int&, const int&, const int&,
                                  const double *, const int&, const int&, const int&, const int&,
                                  double *, const int&, const int&, const int&, const int&, const int&);
  void multicomponent_divide2d_(const int&, const int&, const int&, const int&,
                                const double *, const int&, const int&, const int&, const int&,
                                double *, const int&, const int&, const int&, const int&, const int&);
#endif
#if NDIM==3

  // 3d prototypes

  void compute_face_coef3d_(const int&, const int&, const int&, const int&, const int&, const int&,
                            const int&,
                            const double &,
                            const double *, const int&, const int&, const int&, const int&, const int&, const int&,
                            const double *, const int&, const int&, const int&, const int&, const int&, const int&,
                            const double *, const int&, const int&, const int&, const int&, const int&, const int&,
                            const double *, const int&, const int&, const int&, const int&, const int&, const int&,
                            const double *, const int&, const int&, const int&, const int&, const int&, const int&,
                            const double *, const int&, const int&, const int&, const int&, const int&, const int&,
                            const double *, const int&, const int&, const int&, const int&, const int&, const int&,
                            const double *, const int&, const int&, const int&, const int&, const int&, const int&,
                            const double *, const int&, const int&, const int&, const int&, const int&, const int&,
                            const double &, const char*);
  void compute_dquatdphi_face_coef3d_(const int&, const int&, const int&, const int&, const int&, const int&,
                                      const int&,
                                      const double *, const int&, const int&, const int&, const int&, const int&, const int&,
                                      const double *, const int&, const int&, const int&, const int&, const int&, const int&,
                                      const double *, const int&, const int&, const int&, const int&, const int&, const int&,
                                      const double *, const int&, const int&, const int&, const int&, const int&, const int&,
                                      const double *, const int&, const int&, const int&, const int&, const int&, const int&,
                                      const double *, const int&, const int&, const int&, const int&, const int&, const int&,
                                      const double *, const int&, const int&, const int&, const int&, const int&, const int&);
  void compute_flux3d_(const int&, const int&, const int&, const int&, const int&, const int&,
                       const int&,
                       const double *, const int&, const int&, const int&, const int&, const int&, const int&,
                       const double *, const int&, const int&, const int&, const int&, const int&, const int&,
                       const double *, const int&, const int&, const int&, const int&, const int&, const int&,
                       const double *, const int&, const int&, const int&, const int&, const int&, const int&,
                       const double *,
                       const double *, const int&, const int&, const int&, const int&, const int&, const int&,
                       const double *, const int&, const int&, const int&, const int&, const int&, const int&,
                       const double *, const int&, const int&, const int&, const int&, const int&, const int&);
  void compute_flux3d_from_gradq_(const int&, const int&, const int&, const int&,const int&, const int&,
                       const int&,
                       const double *, const int&, const int&, const int&, const int&,const int&, const int&,
                       const double *, const int&, const int&, const int&, const int&,const int&, const int&,
                       const double *, const int&, const int&, const int&, const int&,const int&, const int&,
                       const double *, 
                       const double *, 
                       const double *, 
                       const double *, const int&, const int&, const int&, const int&, const int&, const int&,
                       const double *, const int&, const int&, const int&, const int&, const int&, const int&,
                       const double *, const int&, const int&, const int&, const int&, const int&, const int&);
  void compute_sym_flux3d_(const int&, const int&, const int&, const int&, const int&, const int&,
                           const int&,
                           const double *, const int&, const int&, const int&, const int&, const int&, const int&,
                           const double *, const int&, const int&, const int&, const int&, const int&, const int&,
                           const double *, const int&, const int&, const int&, const int&, const int&, const int&,
                           const double *, const int&, const int&, const int&, const int&, const int&, const int&,
                           const double *, const int&, const int&, const int&, const int&, const int&, const int&,
                           const double *,
                           const double *, const int&, const int&, const int&, const int&, const int&, const int&,
                           const double *, const int&, const int&, const int&, const int&, const int&, const int&,
                           const double *, const int&, const int&, const int&, const int&, const int&, const int&);
  void compute_q_residual3d_(const int&, const int&, const int&, const int&, const int&, const int&,
                             const int&,
                             const double *, const int&, const int&, const int&, const int&, const int&, const int&,
                             const double *, const int&, const int&, const int&, const int&, const int&, const int&,
                             const double *, const int&, const int&, const int&, const int&, const int&, const int&,
                             const double *, const int&, const int&, const int&, const int&, const int&, const int&,
                             const double *, const int&, const int&, const int&, const int&, const int&, const int&,
                             const double *, const double &,
                             const double *, const int&, const int&, const int&, const int&, const int&, const int&,
                             const double *, const int&, const int&, const int&, const int&, const int&, const int&);
  void compute_q_residual3d_symm_(const int&, const int&, const int&, const int&, const int&, const int&,
                             const int&,
                             const double *, const int&, const int&, const int&, const int&, const int&, const int&,
                             const double *, const int&, const int&, const int&, const int&, const int&, const int&,
                             const double *, const int&, const int&, const int&, const int&, const int&, const int&,
                             const double *, const int&, const int&, const int&, const int&, const int&, const int&,
                             const double *, const int&, const int&, const int&, const int&, const int&, const int&,
                             const double *, const double &,
                             const double *, const int&, const int&, const int&, const int&, const int&, const int&,
                             const double *, const int&, const int&, const int&, const int&, const int&, const int&,
                             const int*, const int*, const int*, const int&);
  void add_quat_op3d_(const int&, const int&, const int&, const int&, const int&, const int&,
                      const int&,
                      const double *, const int&, const int&, const int&, const int&, const int&, const int&,
                      const double *, const int&, const int&, const int&, const int&, const int&, const int&,
                      const double *, const int&, const int&, const int&, const int&, const int&, const int&,
                      const double *, const int&, const int&, const int&, const int&, const int&, const int&,
                      const double *,
                      const double *, const int&, const int&, const int&, const int&, const int&, const int&);
  void add_quat_proj_op3d_(const int&, const int&, const int&, const int&, const int&, const int&,
                           const int&,
                           const double *, const int&, const int&, const int&, const int&, const int&, const int&,
                           const double *, const int&, const int&, const int&, const int&, const int&, const int&,
                           const double *, const int&, const int&, const int&, const int&, const int&, const int&,
                           const double *, const int&, const int&, const int&, const int&, const int&, const int&,
                           const double *, const int&, const int&, const int&, const int&, const int&, const int&,
                           const double *, const int&, const int&, const int&, const int&, const int&, const int&,
                           const double *,
                           const double *, const int&, const int&, const int&, const int&, const int&, const int&);
  void compute_lambda_flux3d_(const int&, const int&, const int&, const int&, const int&, const int&,
                              const int&,
                              const double *, const int&, const int&, const int&, const int&, const int&, const int&,
                              const double *, const int&, const int&, const int&, const int&, const int&, const int&,
                              const double *, const int&, const int&, const int&, const int&, const int&, const int&,
                              const double *, const int&, const int&, const int&, const int&, const int&, const int&,
                              const double *,
                              const double *, const int&, const int&, const int&, const int&, const int&, const int&);
  void fixflux3d_(double *, double *, double *,
                  const int &, const int &, const int &, const int &,
                  const double *, const double *, const double *,
                  const int &, const int &, const int &, const int &,
                  const double *, const int &, const int &, const int &, const int &,
                  const int &, const int &, const int &, const int &,
                  const int &, const int &, const int &,
                  const int &,
                  const int *,
                  const int *, const int *,
                  const double *);
  void project3d_(const int&, const int&, const int&, const int&, const int&, const int&,
                  const int&,
                  const double *, const int&, const int&, const int&, const int&, const int&, const int&,
                  const double *, const int&, const int&, const int&, const int&, const int&, const int&,
                  const double *, const int&, const int&, const int&, const int&, const int&, const int&);
  void take_square_root3d_(const int&, const int&, const int&, const int&, const int&, const int&,
                           double *, const int&, const int&, const int&, const int&, const int&, const int&);
  void multicomponent_multiply3d_(const int&, const int&, const int&, const int&, const int&, const int&,
                                  const double *, const int&, const int&, const int&, const int&, const int&, const int&,
                                  double *, const int&, const int&, const int&, const int&, const int&, const int&, const int&);
  void multicomponent_divide3d_(const int&, const int&, const int&, const int&, const int&, const int&,
                                const double *, const int&, const int&, const int&, const int&, const int&, const int&,
                                double *, const int&, const int&, const int&, const int&, const int&, const int&, const int&);

#endif

}
