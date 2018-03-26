/*******************************************************************
 *                                                                 *
 * File          : nvector_cplx.h                                  *
 * Programmers   : Andrei Schaffer and Radu Serban @ LLNL          *
 * Version of    : 03 July 2003                                    *
 *-----------------------------------------------------------------*
 * Copyright (c) 2002, The Regents of the University of California *
 * Produced at the Lawrence Livermore National Laboratory          *
 * All rights reserved                                             *
 * For details, see sundials/shared/LICENSE                        *
 *-----------------------------------------------------------------*
 * This is the header file for a generic complex NVECTOR package.  *
 *******************************************************************/

#ifndef included_nvector_cplx_h
#define included_nvector_cplx_h

#include "nvector.h"
#include "complexify.h"

/****************************************************************
 * Functions exported by nvector_cplx                           *
 ****************************************************************/

cplx_vector N_VNew_Cplx(NV_Spec nvSpec);

void N_VFree_Cplx(N_Vector v_re, N_Vector v_im);

cplx_vector N_VMake_Cplx(realtype *v_re, realtype *v_im , NV_Spec nvSpec);

void N_VDispose_Cplx(N_Vector v_re, N_Vector v_im);

cplx_vector N_VGetData_Cplx(N_Vector v_re, N_Vector v_im);

void N_VSetData_Cplx(realtype *v_data_re, realtype *v_data_im, N_Vector_Cplx v);

/*--------------------------------------------------------------*
 * Function  : N_VLinearSum                                     *
 * Operation : z = a x + b y                                    *
 *--------------------------------------------------------------*/

void N_VLinearSum_Cplx(cplx a, N_Vector x_re, N_Vector x_im, 
                       cplx b, N_Vector y_re, N_Vector y_im,
                       N_Vector z_re, N_Vector z_im);
void N_VLinearSum_Cplx(realtype a, N_Vector x, 
                       realtype b, N_Vector y, 
                       N_Vector z_re, N_Vector z_im);
void N_VLinearSum_Cplx(cplx a, N_Vector x, 
                       realtype b, N_Vector y, 
                       N_Vector z_re, N_Vector z_im);
void N_VLinearSum_Cplx(realtype a, N_Vector x_re, N_Vector x_im, 
                       realtype b, N_Vector y, 
                       N_Vector z_re, N_Vector z_im);
void N_VLinearSum_Cplx(realtype a, N_Vector x, 
                       cplx b, N_Vector y, 
                       N_Vector z_re, N_Vector z_im);
void N_VLinearSum_Cplx(realtype a, N_Vector x, 
                       realtype b, N_Vector y_re, N_Vector y_im, 
                       N_Vector z_re, N_Vector z_im);
void N_VLinearSum_Cplx(cplx a, N_Vector x_re, N_Vector x_im, 
                       realtype b, N_Vector y, 
                       N_Vector z_re, N_Vector z_im);
void N_VLinearSum_Cplx(cplx a, N_Vector x, 
                       cplx b, N_Vector y, 
                       N_Vector z_re, N_Vector z_im);
void N_VLinearSum_Cplx(realtype a, N_Vector x_re, N_Vector x_im, 
                       cplx b, N_Vector y, 
                       N_Vector z_re, N_Vector z_im);
void N_VLinearSum_Cplx(cplx a, N_Vector x, 
                       realtype b, N_Vector y_re, N_Vector y_im, 
                       N_Vector z_re, N_Vector z_im);
void N_VLinearSum_Cplx(realtype a, N_Vector x_re, N_Vector x_im, 
                       realtype b, N_Vector y_re, N_Vector y_im, 
                       N_Vector z_re, N_Vector z_im);
void N_VLinearSum_Cplx(realtype a, N_Vector x, 
                       cplx b, N_Vector y_re, N_Vector y_im, 
                       N_Vector z_re, N_Vector z_im);
void N_VLinearSum_Cplx(cplx a, N_Vector x_re, N_Vector x_im, 
                       cplx b, N_Vector y, 
                       N_Vector z_re, N_Vector z_im);
void N_VLinearSum_Cplx(cplx a, N_Vector x_re, N_Vector x_im, 
                       realtype b, N_Vector y_re, N_Vector y_im, 
                       N_Vector z_re, N_Vector z_im);
void N_VLinearSum_Cplx(cplx a, N_Vector x, 
                       cplx b, N_Vector y_re, N_Vector y_im, 
                       N_Vector z_re, N_Vector z_im);
void N_VLinearSum_Cplx(realtype a, N_Vector x_re, N_Vector x_im, 
                       cplx b, N_Vector y_re, N_Vector y_im, 
                       N_Vector z_re, N_Vector z_im);

/*--------------------------------------------------------------*
 * Function  : N_VConst                                         *
 * Operation : z[i] = c for i=0, 1, ..., N-1                    *
 *--------------------------------------------------------------*/

void N_VConst_Cplx(cplx c, N_Vector z_re, N_Vector z_im);
void N_VConst_Cplx(realtype c, N_Vector z_re, N_Vector z_im);

/*--------------------------------------------------------------*
 * Function  : N_VProd                                          *
 * Operation : z[i] = x[i] * y[i] for i=0, 1, ..., N-1          *
 *--------------------------------------------------------------*/

void N_VProd_Cplx(N_Vector_Cplx x, N_Vector_Cplx y, N_Vector_Cplx z);

void N_VProd_Cplx(N_Vector_Cplx x, N_Vector y, N_Vector_Cplx z);

void N_VProd_Cplx(N_Vector x, N_Vector_Cplx y, N_Vector_Cplx z);

void N_VProd_Cplx(N_Vector x, N_Vector y, N_Vector_Cplx z);

/*--------------------------------------------------------------*
 * Function  : N_VDiv                                           *
 * Operation : z[i] = x[i] / y[i] for i=0, 1, ..., N-1          *
 *--------------------------------------------------------------*/

void N_VDiv_Cplx(N_Vector_Cplx x, N_Vector_Cplx y, N_Vector_Cplx z);

void N_VDiv_Cplx(N_Vector_Cplx x, N_Vector y, N_Vector_Cplx z);

void N_VDiv_Cplx(N_Vector x, N_Vector_Cplx y, N_Vector_Cplx z);

void N_VDiv_Cplx(N_Vector x, N_Vector y, N_Vector_Cplx z);

/*--------------------------------------------------------------*
 * Function  : N_VScale                                         *
 * Operation : z = c x                                          *
 *--------------------------------------------------------------*/

void N_VScale_Cplx(cplx c, N_Vector_Cplx x, N_Vector_Cplx z);

void N_VScale_Cplx(cplx c, N_Vector x, N_Vector_Cplx z);

void N_VScale_Cplx(realtype c, N_Vector_Cplx x, N_Vector_Cplx z);

void N_VScale_Cplx(realtype c, N_Vector x, N_Vector_Cplx z);

/*--------------------------------------------------------------*
 * Function  : N_VAbs                                           *
 * Operation : z[i] = |x[i]|,   for i=0, 1, ..., N-1            *
 *--------------------------------------------------------------*/

void N_VAbs_Cplx(N_Vector_Cplx x, N_Vector z);

/*--------------------------------------------------------------*
 * Function  : N_VInv                                           *
 * Operation : z[i] = 1.0 / x[i] for i = 0, 1, ..., N-1         *
 *--------------------------------------------------------------*

void N_VInv_Cplx(N_Vector_Cplx x, N_Vector_Cplx z);

void N_VInv_Cplx(N_Vector x, N_Vector_Cplx z);

/*--------------------------------------------------------------*
 * Function  : N_VAddConst                                      *
 * Operation : z[i] = x[i] + b   for i = 0, 1, ..., N-1         *
 *--------------------------------------------------------------*/

void N_VAddConst_Cplx(N_Vector_Cplx x, cplx b, N_Vector_Cplx z);

void N_VAddConst_Cplx(N_Vector x, cplx b, N_Vector_Cplx z);

void N_VAddConst_Cplx(N_Vector_Cplx x, realtype b, N_Vector_Cplx z);

void N_VAddConst_Cplx(N_Vector x, realtype b, N_Vector_Cplx z);

/*--------------------------------------------------------------*
 * Function : N_VDotProd                                        *
 * Usage    : dotprod = N_VDotProd_Cplx(x, y);                  *
 *--------------------------------------------------------------*/

cplx N_VDotProd_Cplx(N_Vector_Cplx x, N_Vector_Cplx y);

cplx N_VDotProd_Cplx(N_Vector x, N_Vector_Cplx y);

cplx N_VDotProd_Cplx(N_Vector_Cplx x, N_Vector y);


#endif

