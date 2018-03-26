/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   F77 external declarations for SAMRAI Euler gas dynamics ex.
 *
 ************************************************************************/

extern "C" {

void F77_FUNC(eulerinit, EULERINIT) (
   const int&, const double *, const double *, const double *,
   const int&, const int&,
   const int&, const int&,
#if (NDIM > 2)
   const int&, const int&,
#endif
   const int&,
   const int&,
#if (NDIM > 2)
   const int&,
#endif
   const double&,
   double *, double *, double *,
   const int&,
   const double *,
   const double *, const double *, const double *);

void F77_FUNC(eulerinitsphere, EULERINITSPHERE) (
   const int&, const double *, const double *, const double *,
   const int&, const int&,
   const int&, const int&,
#if (NDIM > 2)
   const int&, const int&,
#endif
   const int&,
   const int&,
#if (NDIM > 2)
   const int&,
#endif
   const double&,
   double *, double *, double *,
   const double&, const double *, const double&,
   const double&, const double *, const double&,
   const double *, const double&);

void F77_FUNC(stabledt, STABLEDT) (
   const double *,
   const int&, const int&,
   const int&, const int&,
#if (NDIM > 2)
   const int&, const int&,
#endif
   const int&,
   const int&,
#if (NDIM > 2)
   const int&,
#endif
   const double&,
   const double *, const double *, const double *, double&);

void F77_FUNC(inittraceflux, INITTRACEFLUX) (
   const int&, const int&,
   const int&, const int&,
#if (NDIM > 2)
   const int&, const int&,
#endif
   const double *, const double *, const double *,
   double *, double *, double *,
#if (NDIM > 2)
   double *, double *, double *,
#endif
   double *, double *, double *);

void F77_FUNC(computesound, COMPUTESOUND) (
   const int&, const int&,
   const int&, const int&,
#if (NDIM > 2)
   const int&, const int&,
#endif
   const double&,
   const double *, const double *, const double *,
   double *);

void F77_FUNC(chartracing0, CHARTRACING0) (
   const double&,
   const int&, const int&,
   const int&, const int&,
#if (NDIM > 2)
   const int&, const int&,
#endif
   const int&, const double&, const double&, const int&,
   const double *,
   double *, double *,
   double *, double *,
   double *,
   double *, double *);

void F77_FUNC(chartracing1, CHARTRACING1) (
   const double&,
   const int&, const int&,
   const int&, const int&,
#if (NDIM > 2)
   const int&, const int&,
#endif
   const int&, const double&, const double&, const int&,
   const double *,
   double *, double *,
   double *, double *,
   double *,
   double *, double *);

#if (NDIM == 3)
void F77_FUNC(chartracing2, CHARTRACING2) (
   const double&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const double&, const double&, const int&,
   const double *,
   double *, double *,
   double *, double *,
   double *,
   double *, double *);
#endif

void F77_FUNC(fluxcalculation, FLUXCALCULATION) (
   const double&, const int&,
#if (NDIM > 2)
   const int&,
#endif
   const int&,
   const double *,
   const int&, const int&,
   const int&, const int&,
#if (NDIM > 2)
   const int&, const int&,
#endif
   const double&,
   const int&,
   const double *, const double *, const double *,
#if (NDIM > 2)
   double *, double *, double *,
#endif
   double *, double *, double *,
   double *, double *, double *);

#if (NDIM == 3)
void F77_FUNC(fluxcorrec2d, FLUXCORREC2D) (
   const double&,
   const int&, const int&, const int&, const int&, const int&, const int&,
   const double *, const double&, const int&,
   const double *,
   const double *,
   const double *,
   const double *, const double *, const double *,
   const double *, const double *, const double *,
   const double *, const double *, const double *,
   double *, double *, double *,
   double *, double *, double *);

void F77_FUNC(fluxcorrec3d, FLUXCORREC3D) (
   const double&,
   const int&, const int&, const int&, const int&, const int&, const int&,
   const double *, const double&,
   const double *, const double *, const double *,
   const double *, const double *, const double *,
   const double *, const double *, const double *,
   double *, double *, double *,
   double *, double *, double *);
#endif

#if (NDIM == 2)
void F77_FUNC(fluxcorrec, FLUXCORREC) (
   const double&,
   const int&, const int&, const int&, const int&,
   const double *, const double&,
   const double *, const double *, const double *,
   double *, double *,
   double *, double *,
   double *, double *);
#endif

void F77_FUNC(consdiff, CONSDIFF) (
   const int&, const int&,
   const int&, const int&,
#if (NDIM > 2)
   const int&, const int&,
#endif
   const double *,
   const double *,
   const double *,
#if (NDIM > 2)
   const double *,
#endif
   const double&,
   double *, double *, double *);

#if (NDIM > 2)
void F77_FUNC(onethirdstate, ONETHIRDSTATE) (
   const double&, const double *, const int&,
   const int&, const int&, const int&, const int&, const int&, const int&,
   const double&,
   const double *, const double *, const double *,
   const double *, const double *, const double *,
   double *);

void F77_FUNC(fluxthird, FLUXTHIRD) (
   const double&, const double *, const int&,
   const int&, const int&, const int&, const int&, const int&, const int&,
   const double&,
   const int&,
   const double *, const double *, const double *, const double *,
   double *, double *, double *);

void F77_FUNC(fluxcorrecjt, FLUXCORRECJT) (
   const double&, const double *, const int&,
   const int&, const int&, const int&, const int&, const int&, const int&,
   const double&,
   const double *, const double *, const double *,
   const double *, const double *, const double *,
   double *, double *, double *,
   double *, double *, double *);
#endif

#if (NDIM == 2)
void F77_FUNC(conservlinint2d, CONSERVLININT2D) (
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int *, const double *, const double *, const double&,
   const double *, const double *,
   const double *, const double *,
   double *, double *,
   double *,
   double *, double *, double *, const int&,
   double *, double *, double *,
   double *, double *,
   double *, double *, double *, double *);

void F77_FUNC(conservavg2d, CONSERVAVG2D) (
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int *, const double *, const double *, const double&,
   const double *, const double *,
   const double *, const double *,
   double *, double *,
   double *);
#endif
#if (NDIM == 3)
void F77_FUNC(conservlinint3d, CONSERVLININT3D) (
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int *, const double *, const double *, const double&,
   const double *, const double *,
   const double *, const double *,
   double *, double *,
   double *,
   double *, double *, double *, const int&,
   double *, double *, double *,
   double *, double *, double *,
   double *, double *, double *, double *, double *, double *);

void F77_FUNC(conservavg3d, CONSERVAVG3D) (
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int *, const double *, const double *, const double&,
   const double *, const double *,
   const double *, const double *,
   double *, double *,
   double *);
#endif

void F77_FUNC(detectgrad, DETECTGRAD) (
   const int&, const int&,
   const int&, const int&,
#if (NDIM > 2)
   const int&, const int&,
#endif
   const int&, const int&, const int&,
   const int&, const int&, const int&,
#if (NDIM > 2)
   const int&, const int&, const int&,
#endif
   const double *,
   const double&,
   const int&, const int&,
   const double *,
   int *, int *);

void F77_FUNC(detectshock, DETECTSHOCK) (
   const int&, const int&,
   const int&, const int&,
#if (NDIM > 2)
   const int&, const int&,
#endif
   const int&, const int&, const int&,
   const int&, const int&, const int&,
#if (NDIM > 2)
   const int&, const int&, const int&,
#endif
   const double *,
   const double&, const double&,
   const int&, const int&,
   const double *,
   int *, int *);

void F77_FUNC(stufprobc, STUFPROBC) (
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&,
   const int&, const int&, const int&);
}
