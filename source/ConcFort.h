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
#define FORT_CONCENTRATION_FLUX concentrationflux_
#define FORT_CONCENTRATION_FLUX_EBS concentrationflux_ebs_
#define FORT_ADD_CONCENTRATION_FLUX_EBS add_concentrationflux_ebs_
#define FORT_CONCENTRATION_FLUX_SPINODAL concentrationflux_spinodal_
#define FORT_CONCENTRATIONDIFFUSION0 concentrationdiffusion0_
#define FORT_CONCENTRATIONDIFFUSIONBECKERMANN concentrationdiffusion_beckermann_
#define FORT_COMPUTERHSCONCENTRATION computerhsconcentration_
#define FORT_PHASERHS_FENERGY      phaserhs_fenergy_
#define FORT_ETARHS_FENERGY      etarhs_fenergy_
#define FORT_CALPHAD_CONC_SOLV_THREE calphad_c_three_
#define FORT_CALPHAD_CONC_SOLV_TWO calphad_c_two_
#define FORT_CALPHAD_CONC_CEQ_TWO calphad_ceq_two_
#define FORT_ADD_CONCENTRATION_FLUX_FROM_GRADT addconcentrationfluxfromgradt_
#define FORT_ADD_CONCENTRATION_FLUX_FROM_AT addconcentrationfluxfromantitrapping_
#define FORT_INITGAUSSIAN initgaussian_
#define FORT_INITGAUSSIAN_SOURCE initgaussiansource_
#define FORT_LINEARMELTINGLINE linearmeltingline_
#define FORT_COMPUTE_CONCENTRATION_FROM_PHASE_CONCENTRATIONS compute_concentration_from_phase_concentrations_


// Function argument list interfaces
extern "C" {

   void FORT_CONCENTRATION_FLUX(
      const int& ifirst0, const int& ilast0,
      const int& ifirst1, const int& ilast1,
#if (NDIM == 3)
      const int& ifirst2, const int& ilast2,
#endif
      const double* dx,
      const double* conc, const int& ngconc,
      const double* phase, const int& ngphase,
      const double* eta, const int& ngeta,
      const double* diffconc0, const double* diffconc1,
#if (NDIM == 3)
      const double* diffconc2,
#endif
      const int& ngdiffconc,
      const double* dphicoupl0, const double* dphicoupl1,
#if (NDIM == 3)
      const double* dphicoupl2,
#endif
      const int& ngdphicoupl,
      const double* detacoupl0, const double* detacoupl1,
#if (NDIM == 3)
      const double* detacoupl2,
#endif
      const int& ngdetacoupl,
      const double* flux0, const double* flux1,
#if (NDIM == 3)
      const double* flux2,
#endif
      const int& ngflux,
      const int& three_phase
      );

   void FORT_ADD_CONCENTRATION_FLUX_FROM_GRADT(
      const int& ifirst0, const int& ilast0,
      const int& ifirst1, const int& ilast1,
#if (NDIM == 3)
      const int& ifirst2, const int& ilast2,
#endif
      const double* dx,
      const double* temperature, const int& ngtemperature,
      const double* mq0, const double* mq1,
#if (NDIM == 3)
      const double* mq2,
#endif
      const int& ngmq,
      const double* flux0, const double* flux1,
#if (NDIM == 3)
      const double* flux2,
#endif
      const int& ngflux,
      const char* const avg_type
      );

   void FORT_ADD_CONCENTRATION_FLUX_FROM_AT(
      const int& ifirst0, const int& ilast0,
      const int& ifirst1, const int& ilast1,
#if (NDIM == 3)
      const int& ifirst2, const int& ilast2,
#endif
      const double* dx,
      const double* phase, const int& ngphase,
      const double* cl, const double* ca, const int& ngc,
      const double* dphi, const int& ngdphi,
      const double& alpha,
      const double* flux0, const double* flux1,
#if (NDIM == 3)
      const double* flux2,
#endif
      const int& ngflux
      );

   void FORT_ADD_CONCENTRATION_FLUX_EBS(
      const int& ifirst0, const int& ilast0,
      const int& ifirst1, const int& ilast1,
#if (NDIM == 3)
      const int& ifirst2, const int& ilast2,
#endif
      const double* dx,
      const double* conc, const int& ngconc,
      const int& ncomp,
      const double* diffconc0, const double* diffconc1,
#if (NDIM == 3)
      const double* diffconc2,
#endif
      const int& ngdiffconc,
      const double* flux0, const double* flux1,
#if (NDIM == 3)
      const double* flux2,
#endif
      const int& ngflux
      );

   void FORT_CONCENTRATION_FLUX_SPINODAL(
      const int& ifirst0, const int& ilast0,
      const int& ifirst1, const int& ilast1,
#if (NDIM == 3)
      const int& ifirst2, const int& ilast2,
#endif
      const double* dx,
      const double* conc, const int& ngconc,
      const int& ncomp,
      const double* conca, const int& ngconca,
      const double* concb, const int& ngconcb,
      const double* diffconc0, const double* diffconc1,
#if (NDIM == 3)
      const double* diffconc2,
#endif
      const int& ngdiffconcl,
      const double* eta, const int& ngeta,
      const double& kappa,
      const double* flux0, const double* flux1,
#if (NDIM == 3)
      const double* flux2,
#endif
      const int& ngflux
      );

   void FORT_CONCENTRATIONDIFFUSION0(
      const int& ifirst0, const int& ilast0,
      const int& ifirst1, const int& ilast1,
#if (NDIM == 3)
      const int& ifirst2, const int& ilast2,
#endif
      const double* phi, const int& ngphi,
      const double* eta, const int& ngeta,
      const double* diff0, const double* diff1,
#if (NDIM == 3)
      const double* diff2,
#endif
      const int& ngdiff,
      const double* t, const int& ngt,
      const double& d_phase0, const double& q0_phase0,
      const double& d_phase1, const double& q0_phase1,
      const double& d_phase2, const double& q0_phase2,
      const double& gas_constant_R,
      const char* phi_interp_type,
      const char* avg_func_type,
      const int& three_phase );

   void FORT_CONCENTRATIONDIFFUSIONBECKERMANN(
      const int& ifirst0, const int& ilast0,
      const int& ifirst1, const int& ilast1,
#if (NDIM == 3)
      const int& ifirst2, const int& ilast2,
#endif
      const double* phi, const int& ngphi,
      const double* diff0, 
      const double* diff1,
#if (NDIM == 3)
      const double* diff2,
#endif
      const int& ngdiff,
      const double* partition_coeff, const int& ngk,
      const double& d_phase0, 
      const double& d_phase1, 
      const char* phi_interp_type,
      const char* avg_func_type );

   void FORT_COMPUTERHSCONCENTRATION(
      const int& ifirst0, const int& ilast0,
      const int& ifirst1, const int& ilast1,
#if (NDIM == 3)
      const int& ifirst2, const int& ilast2,
#endif
      const double* dx,
      const double* flux0, const double* flux1,
#if (NDIM == 3)
      const double* flux2,
#endif
      const int& ngflux,
      const double& conc_mobility,
      const double* rhs, const int& ngrhs);

   void FORT_PHASERHS_FENERGY(
      const int& ifirst0, const int& ilast0,
      const int& ifirst1, const int& ilast1,
#if (NDIM == 3)
      const int& ifirst2, const int& ilast2,
#endif
      const double* fl,
      const double* fa,
      const double* fb,
      const double* phi, const int& ngphi,
      const double* eta, const int& ngeta,
      const double* rhs, const int& ngrhs,
      const char* phi_interp_type,
      const char* eta_interp_type,
      const int& three_phase );

   void FORT_ETARHS_FENERGY(
      const int& ifirst0, const int& ilast0,
      const int& ifirst1, const int& ilast1,
#if (NDIM == 3)
      const int& ifirst2, const int& ilast2,
#endif
      const double* fl,
      const double* fa,
      const double* fb,
      const double* phi, const int& ngphi,
      const double* eta, const int& ngeta,
      const double* rhs, const int& ngrhs,
      const char* phi_interp_type,
      const char* eta_interp_type );

   void FORT_CALPHAD_CONC_SOLV_TWO(
      const double* x,
      const double& c,
      const double& hphi,
      const double& d_gas_constant_R,
      const double* L0,
      const double* L1,
      const double* L2,
      const double* fA,
      const double* fB );
   
   int FORT_CALPHAD_CONC_CEQ_TWO(
      const double* x,
      const double& cl,
      const double& cs,
      const double& d_gas_constant_R,
      const double* L0,
      const double* L1,
      const double* L2,
      const double* fA,
      const double* fB );
   
   void FORT_CALPHAD_CONC_SOLV_THREE(
      const double* x,
      const double& c,
      const double& hphi,
      const double& heta,
      const double& d_gas_constant_R,
      const double* L0,
      const double* L1,
      const double* L2,
      const double* fA,
      const double* fB );
   
   void FORT_INITGAUSSIAN(
      const double*,const double*,const double*,
      const int& ifirst0, const int& ilast0,
      const int& ifirst1, const int& ilast1,
#if (NDIM == 3)
      const int& ifirst2, const int& ilast2,
#endif
      const int&, const int&,
#if (NDIM == 3)
      const int&,
#endif
      const double*,
      const double*,const double*,
      const double&, const double&, const double&);

   void FORT_INITGAUSSIAN_SOURCE(
      const double*,const double*,const double*,
      const int& ifirst0, const int& ilast0,
      const int& ifirst1, const int& ilast1,
#if (NDIM == 3)
      const int&,const int&,
#endif
      const double*,const int&,
      const double*,const int&,
      const double*,const double*,
      const double&, const double&);

   void FORT_LINEARMELTINGLINE(
      const int& ifirst0, const int& ilast0,
      const int& ifirst1, const int& ilast1,
#if (NDIM == 3)
      const int&,const int&,
#endif
      const double*,const int&,
      const double*,const int&,
      const double&, const double&, const double&);


   void FORT_COMPUTE_CONCENTRATION_FROM_PHASE_CONCENTRATIONS(
      const int& ifirst0, const int& ilast0,
      const int& ifirst1, const int& ilast1,
#if (NDIM == 3)
      const int&,const int&,
#endif
      const double*,const int&,
      const double*,const double*,const int&,
      const double*,const int&,
      const char* const);
}

