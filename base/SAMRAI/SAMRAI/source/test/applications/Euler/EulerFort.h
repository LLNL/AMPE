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

void F77_FUNC(eulerinit2d, EULERINIT2D) (
   const int&, const double *, const double *, const double *,
   const int&, const int&,
   const int&, const int&,
   const int&,
   const int&,
   const double&,
   double *, double *, double *,
   const int&,
   const double *,
   const double *, const double *, const double *);

void F77_FUNC(eulerinit3d, EULERINIT3D) (
   const int&, const double *, const double *, const double *,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&,
   const int&,
   const int&,
   const double&,
   double *, double *, double *,
   const int&,
   const double *,
   const double *, const double *, const double *);

void F77_FUNC(eulerinitsphere2d, EULERINITSPHERE2D) (
   const int&, const double *, const double *, const double *,
   const int&, const int&,
   const int&, const int&,
   const int&,
   const int&,
   const double&,
   double *, double *, double *,
   const double&, const double *, const double&,
   const double&, const double *, const double&,
   const double *, const double&);

void F77_FUNC(eulerinitsphere3d, EULERINITSPHERE3D) (
   const int&, const double *, const double *, const double *,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&,
   const int&,
   const int&,
   const double&,
   double *, double *, double *,
   const double&, const double *, const double&,
   const double&, const double *, const double&,
   const double *, const double&);

void F77_FUNC(stabledt2d, STABLEDT2D) (
   const double *,
   const int&, const int&,
   const int&, const int&,
   const int&,
   const int&,
   const double&,
   const double *, const double *, const double *, double&);

void F77_FUNC(stabledt3d, STABLEDT3D) (
   const double *,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&,
   const int&,
   const int&,
   const double&,
   const double *, const double *, const double *, double&);

void F77_FUNC(inittraceflux2d, INITTRACEFLUX2D) (
   const int&, const int&,
   const int&, const int&,
   const double *, const double *, const double *,
   double *, double *, double *,
   double *, double *, double *);

void F77_FUNC(inittraceflux3d, INITTRACEFLUX3D) (
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const double *, const double *, const double *,
   double *, double *, double *,
   double *, double *, double *,
   double *, double *, double *);

void F77_FUNC(computesound2d, COMPUTESOUND2D) (
   const int&, const int&,
   const int&, const int&,
   const double&,
   const double *, const double *, const double *,
   double *);

void F77_FUNC(computesound3d, COMPUTESOUND3D) (
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const double&,
   const double *, const double *, const double *,
   double *);

void F77_FUNC(chartracing2d0, CHARTRACING2D0) (
   const double&,
   const int&, const int&,
   const int&, const int&,
   const int&, const double&, const double&, const int&,
   const double *,
   double *, double *,
   double *, double *,
   double *,
   double *, double *);

void F77_FUNC(chartracing3d0, CHARTRACING3D0) (
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

void F77_FUNC(chartracing2d1, CHARTRACING2D1) (
   const double&,
   const int&, const int&,
   const int&, const int&,
   const int&, const double&, const double&, const int&,
   const double *,
   double *, double *,
   double *, double *,
   double *,
   double *, double *);

void F77_FUNC(chartracing3d1, CHARTRACING3D1) (
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

void F77_FUNC(chartracing3d2, CHARTRACING3D2) (
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

void F77_FUNC(fluxcalculation2d, FLUXCALCULATION2D) (
   const double&, const int&,
   const int&,
   const double *,
   const int&, const int&,
   const int&, const int&,
   const double&,
   const int&,
   const double *, const double *, const double *,
   double *, double *, double *,
   double *, double *, double *);

void F77_FUNC(fluxcalculation3d, FLUXCALCULATION3D) (
   const double&, const int&,
   const int&,
   const int&,
   const double *,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const double&,
   const int&,
   const double *, const double *, const double *,
   double *, double *, double *,
   double *, double *, double *,
   double *, double *, double *);

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

void F77_FUNC(fluxcorrec, FLUXCORREC) (
   const double&,
   const int&, const int&, const int&, const int&,
   const double *, const double&,
   const double *, const double *, const double *,
   double *, double *,
   double *, double *,
   double *, double *);

void F77_FUNC(consdiff2d, CONSDIFF2D) (
   const int&, const int&,
   const int&, const int&,
   const double *,
   const double *,
   const double *,
   const double&,
   double *, double *, double *);

void F77_FUNC(consdiff3d, CONSDIFF3D) (
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const double *,
   const double *,
   const double *,
   const double *,
   const double&,
   double *, double *, double *);

void F77_FUNC(onethirdstate, onethirdstate) (
   const double&, const double *, const int&,
   const int&, const int&, const int&, const int&, const int&, const int&,
   const double&,
   const double *, const double *, const double *,
   const double *, const double *, const double *,
   double *);

void F77_FUNC(fluxthird, fluxthird) (
   const double&, const double *, const int&,
   const int&, const int&, const int&, const int&, const int&, const int&,
   const double&,
   const int&,
   const double *, const double *, const double *, const double *,
   double *, double *, double *);

void F77_FUNC(fluxcorrecjt, fluxcorrecjt) (
   const double&, const double *, const int&,
   const int&, const int&, const int&, const int&, const int&, const int&,
   const double&,
   const double *, const double *, const double *,
   const double *, const double *, const double *,
   double *, double *, double *,
   double *, double *, double *);

void F77_FUNC(conservlinint2d, conservlinint2d) (
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

void F77_FUNC(conservavg2d, conservavg2d) (
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int *, const double *, const double *, const double&,
   const double *, const double *,
   const double *, const double *,
   double *, double *,
   double *);
void F77_FUNC(conservlinint3d, conservlinint3d) (
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

void F77_FUNC(conservavg3d, conservavg3d) (
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

void F77_FUNC(detectgrad2d, detectgrad2d) (
   const int&, const int&,
   const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const double *,
   const double&,
   const int&, const int&,
   const double *,
   int *, int *);

void F77_FUNC(detectgrad3d, detectgrad3d) (
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const double *,
   const double&,
   const int&, const int&,
   const double *,
   int *, int *);

void F77_FUNC(detectshock2d, detectshock2d) (
   const int&, const int&,
   const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const double *,
   const double&, const double&,
   const int&, const int&,
   const double *,
   int *, int *);

void F77_FUNC(detectshock3d, detectshock3d) (
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const double *,
   const double&, const double&,
   const int&, const int&,
   const double *,
   int *, int *);

void F77_FUNC(stufprobc, stufprobc) (
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&,
   const int&, const int&, const int&);
}
