/*******************************************************************
 * File          : nvector_cplx.c                                  *
 * Programmers   : Andrei Schaffer and Radu Serban @ LLNL          *
 * Version of    : 07 July 2003                                    *
 *-----------------------------------------------------------------*
 * Copyright (c) 2002, The Regents of the University of California *
 * Produced at the Lawrence Livermore National Laboratory          *
 * All rights reserved                                             *
 * For details, see sundials/shared/LICENSE                        *
 *-----------------------------------------------------------------*
 * This is the implementation file for a generic complex N_Vector  *
 * package. It contains the implementation of the kernels listed   *
 * kernels listed in nvector_cplx.hpp.                             *
 *                                                                 *
 *******************************************************************/
 
#include <stdlib.h>
#include "nvector_cplx.h"

cplx_vector N_VNew_Cplx(NV_Spec nvSpec)
{
  cplx_vector tmp;

  tmp.re = N_VNew(nvSpec);
  if(tmp.re == NULL) {
    return(NULL);
  }

  tmp.im = N_VNew(nvSpec);
  if(tmp.re == NULL) {
    N_VFree(tmp.re);
    return(NULL);
  }

  return(tmp);
}

void N_VFree_Cplx(N_Vector v_re, N_Vector v_im)
{
  N_VFree(v_re);
  N_VFree(v_im);
}


cplx_vector N_VMake_Cplx(realtype *v_re, realtype *v_im, NV_Spec nvSpec)
{
  cplx_vector tmp;

  tmp.re = N_VMake(v_re, nvSpec);
  if(tmp.re == NULL) {
    return(NULL);
  }

  tmp.im = N_VMake(v_im, nvSpec);
  if(tmp.re == NULL) {
    N_VDispose(tmp.re);
    return(NULL);
  }

  return(tmp);
}

void N_VDispose_Cplx(N_Vector v_re, N_Vector v_im)
{
  N_VDispose(v_re);
  N_VDispose(v_im);
}

cplx_vector N_VGetData_Cplx(N_Vector v_re, N_Vector v_im)
{
  cplx_vector tmp;
  tmp.re = N_VGetData(v_re);
  tmp.im = N_VGetData(v_im);
  return(tmp);
}

void N_VSetData_Cplx(realtype *v_data_re, realtype *v_data_im , 
                     N_Vector v_re, N_Vector v_im)
{
  N_VSetData(v_data_re,v_re);
  N_VSetData(v_data_im,v_im);
}


/*---------------------------------------------------------------------------------------*/

void N_VLinearSum_Cplx(cplx a, N_Vector x_re, N_Vector x_im, 
                       cplx b, N_Vector y_re, N_Vector y_im,
                       N_Vector z_re, N_Vector z_im)
{
  realtype a_re=a.real();
  realtype a_im=a.imag();
  
  realtype b_re=b.real();
  realtype b_im=b.imag();
  
  N_VLinearSum(a_re, x_re, -a_im, x_im, z_re);
  N_VLinearSum(1.0, z_re, b_re, y_re, z_re);
  N_VLinearSum(1.0, z_re, -b_im, y_im, z_re);
  
  N_VLinearSum(a_re, x_im, a_im, x_re, z_im);
  N_VLinearSum(1.0, z_im, b_re, y_im, z_im);
  N_VLinearSum(1.0, z_im, b_im, y_re, z_im);
}

void N_VLinearSum_Cplx(realtype a, N_Vector x, 
                       realtype b, N_Vector y, 
                       N_Vector z_re, N_Vector z_im)
{ 
  N_VLinearSum(a, x, b, y, z_re);
  N_VConst(0.0, z_im);
}

void N_VLinearSum_Cplx(cplx a, N_Vector x, 
                       realtype b, N_Vector y, 
                       N_Vector z_re, N_Vector z_im)
{
  realtype a_re=a.real();
  realtype a_im=a.imag();
  
  N_VLinearSum(a_re, x, b, y, z_re);
  N_VScale(a_im, x, z_im);
}

void N_VLinearSum_Cplx(realtype a, N_Vector x_re, N_Vector x_im, 
                       realtype b, N_Vector y, 
                       N_Vector z_re, N_Vector z_im)
{ 
  N_VLinearSum(a, x_re, b, y, z_re);
  N_VScale(a, x_im, z_im);
}

void N_VLinearSum_Cplx(realtype a, N_Vector x, 
                       cplx b, N_Vector y, 
                       N_Vector z_re, N_Vector z_im)
{
  realtype b_re=b.real();
  realtype b_im=b.imag();
  
  N_VLinearSum(a, x, b_re, y, z_re);
  N_VScale(b_im, y, z_im);
}

void N_VLinearSum_Cplx(realtype a, N_Vector x, 
                       realtype b, N_Vector y_re, N_Vector y_im, 
                       N_Vector z_re, N_Vector z_im)
{ 
  N_VLinearSum(a, x, b, y_re, z_re);
  N_VScale(b, y_im, z_im);
}

void N_VLinearSum_Cplx(cplx a, N_Vector x_re, N_Vector x_im, 
                       realtype b, N_Vector y, 
                       N_Vector z_re, N_Vector z_im)
{
  realtype a_re=a.real();
  realtype a_im=a.imag();
  
  N_VLinearSum(a_re, x_re, b, y, z_re);
  N_VLinearSum(-a_im, x_im, 1.0, z_re, z_re);
  N_VLinearSum(a_re, x_im, a_im, x_re, z_im);
}

void N_VLinearSum_Cplx(cplx a, N_Vector x, 
                       cplx b, N_Vector y, 
                       N_Vector z_re, N_Vector z_im)
{
  realtype a_re=a.real();
  realtype a_im=a.imag();
  realtype b_re=b.real();
  realtype b_im=b.imag();
  
  N_VLinearSum(a_re, x, b_re, y, z_re);
  N_VLinearSum(a_im, x, b_im, y, z_im);
}

void N_VLinearSum_Cplx(realtype a, N_Vector x_re, N_Vector x_im, 
                       cplx b, N_Vector y, 
                       N_Vector z_re, N_Vector z_im)
{ 
  realtype b_re=b.real();
  realtype b_im=b.imag();
  
  N_VLinearSum(a, x_re, b_re, y, z_re);
  N_VLinearSum(a, x_im, b_im, y, z_im);
}

void N_VLinearSum_Cplx(cplx a, N_Vector x, 
                       realtype b, N_Vector y_re, N_Vector y_im, 
                       N_Vector z_re, N_Vector z_im)
{
  realtype a_re=a.real();
  realtype a_im=a.imag();
  
  N_VLinearSum(a_re, x, b, y_re, z_re);
  N_VLinearSum(a_im, x, b, y_im, z_im);
}

void N_VLinearSum_Cplx(realtype a, N_Vector x_re, N_Vector x_im, 
                       realtype b, N_Vector y_re, N_Vector y_im, 
                       N_Vector z_re, N_Vector z_im)
{ 
  N_VLinearSum(a, x_re, b, y_re, z_re);
  N_VLinearSum(a, x_im, b, y_im, z_im);
}

void N_VLinearSum_Cplx(realtype a, N_Vector x, 
                       cplx b, N_Vector y_re, N_Vector y_im, 
                       N_Vector z_re, N_Vector z_im)
{
  realtype b_re=b.real();
  realtype b_im=b.imag();
  
  N_VLinearSum(a, x, b_re, y_re, z_re);
  N_VLinearSum(-b_im, y_im, 1.0, z_re, z_re);
  N_VLinearSum(b_re, y_im, b_im, y_re, z_im);
}

void N_VLinearSum_Cplx(cplx a, N_Vector x_re, N_Vector x_im, 
                       cplx b, N_Vector y, 
                       N_Vector z_re, N_Vector z_im)
{
  realtype a_re=a.real();
  realtype a_im=a.imag();
  
  realtype b_re=b.real();
  realtype b_im=b.imag();
  
  N_VLinearSum(a_re, x_re, b_re, y, z_re);
  N_VLinearSum(-a_im, x_im, 1.0, z_re, z_re);
  N_VLinearSum(a_re, x_im, a_im, x_re, z_im);
  N_VLinearSum(b_im, y, 1.0, z_im, z_im);
}

void N_VLinearSum_Cplx(cplx a, N_Vector x_re, N_Vector x_im, 
                       realtype b, N_Vector y_re, N_Vector y_im, 
                       N_Vector z_re, N_Vector z_im)
{
  realtype a_re=a.real();
  realtype a_im=a.imag();
  
  N_VLinearSum(a_re, x_re, b, y_re, z_re);
  N_VLinearSum(-a_im, x_im, 1.0, z_re, z_re);
  N_VLinearSum(a_re, x_im, a_im, x_re, z_im);
  N_VLinearSum(b, y_im, 1.0, z_im, z_im);
}

void N_VLinearSum_Cplx(cplx a, N_Vector x, 
                       cplx b, N_Vector y_re, N_Vector y_im, 
                       N_Vector z_re, N_Vector z_im)
{
  realtype a_re=a.real();
  realtype a_im=a.imag();
  
  realtype b_re=b.real();
  realtype b_im=b.imag();
  
  N_VLinearSum(b_re, y_re, a_re, x, z_re);
  N_VLinearSum(-b_im, y_im, 1.0, z_re, z_re);
  N_VLinearSum(b_re, y_im, b_im, y_re, z_im);
  N_VLinearSum(a_im, x, 1.0, z_im, z_im);
}

void N_VLinearSum_Cplx(realtype a, N_Vector x_re, N_Vector x_im, 
                       cplx b, N_Vector y_re, N_Vector y_im, 
                       N_Vector z_re, N_Vector z_im)
{  
  realtype b_re=b.real();
  realtype b_im=b.imag();
  
  N_VLinearSum(a, x_re, b_re, y_re, z_re);
  N_VLinearSum(-b_im, y_im, 1.0, z_re, z_re);
  N_VLinearSum(b_re, y_im, b_im, y_re, z_im);
  N_VLinearSum(a, x_im, 1.0, z_im, z_im);
}

/*-------------------------------------------------------------------------*/

void N_VConst_Cplx(cplx c, N_Vector z_re, N_Vector z_im)
{
  N_VConst(c.real(),z_re);
  N_VConst(c.imag(),z_im);
}

void N_VConst_Cplx(realtype c, N_Vector z_re, N_Vector z_im)
{
  N_VConst(c,z_re);
  N_VConst(0.0,z_im);
}

/*-------------------------------------------------------------------------*/


/* RADU---- GOT HERE!!! */

void N_VProd_Cplx(N_Vector_Cplx x, N_Vector_Cplx y, N_Vector_Cplx z)
{
  N_Vector tmp1=N_VNew(x->re->nvspec);
  
  N_VProd(x->re,y->re,tmp1);
  N_VProd(x->im,y->im,z->re);
  N_VLinearSum(1.0, tmp1, -1.0, z->re, z->re);
  
  N_VProd(x->re,y->im,tmp1);
  N_VProd(x->im,y->re,z->im);
  N_VLinearSum(1.0, tmp1, 1.0, z->im, z->im);
  
  N_VFree(tmp1);
}

void N_VProd_Cplx(N_Vector_Cplx x, N_Vector y, N_Vector_Cplx z)
{
  N_VProd(x->re,y,z->re);
  N_VProd(x->im,y,z->im);
}

void N_VProd_Cplx(N_Vector x, N_Vector_Cplx y, N_Vector_Cplx z)
{
  N_VProd(x,y->re,z->re);
  N_VProd(x,y->im,z->im);
}

void N_VProd_Cplx(N_Vector x, N_Vector y, N_Vector_Cplx z)
{
  N_VProd(x,y,z->re);
  N_VConst(0.0,z->im);
}



extern void map_f_to_N_Vector(realtype (*f)(realtype ), N_Vector x, N_Vector z);
extern double sqrt(double);

void N_VDiv_Cplx(N_Vector_Cplx x, N_Vector_Cplx y, N_Vector_Cplx z)
{
  N_VInv_Cplx(y,z);
  N_VProd_Cplx(x, z, z);
}

void N_VDiv_Cplx(N_Vector_Cplx x, N_Vector y, N_Vector_Cplx z)
{
  N_VDiv(x->re,y,z->re);
  N_VDiv(x->im,y,z->im);
}

void N_VDiv_Cplx(N_Vector x, N_Vector_Cplx y, N_Vector_Cplx z)
{
  N_VInv_Cplx(y,z);
  N_VProd(x, z->re, z->re);
  N_VProd(x, z->im, z->im);
}

void N_VDiv_Cplx(N_Vector x, N_Vector y, N_Vector_Cplx z)
{
  N_VDiv(x,y,z->re);
  N_VConst(0.0,z->im);
}

void N_VScale_Cplx(cplx c, N_Vector_Cplx x, N_Vector_Cplx z) 
{
  N_VLinearSum(c.real(), x->re, (-1.0)*(c.imag()), x->im, z->re);
  N_VLinearSum(c.real(), x->im, c.imag(), x->re, z->im);
}

void N_VScale_Cplx(cplx c, N_Vector x, N_Vector_Cplx z) 
{
  N_VScale(c.real(), x, z->re);
  N_VScale(c.imag(), x, z->im);
}

void N_VScale_Cplx(realtype c, N_Vector_Cplx x, N_Vector_Cplx z) 
{
  N_VScale(c, x->re, z->re);
  N_VScale(c, x->im, z->im);
}

void N_VScale_Cplx(realtype c, N_Vector x, N_Vector_Cplx z) 
{
  N_VScale(c, x, z->re);
  N_VConst(0.0,z->im);
}

void N_VAbs_Cplx(N_Vector_Cplx x, N_Vector z)
{
  N_Vector tmp1=N_VNew(z->nvspec);
  
  N_VProd(x->re,x->re,tmp1);
  N_VProd(x->im,x->im,z);
  
  N_VLinearSum(1.0, tmp1, 1.0, z, z);
  map_f_to_N_Vector(sqrt, z, z);
  
  N_VFree(tmp1);
}

void N_VInv_Cplx(N_Vector_Cplx x, N_Vector_Cplx z)
{  
  N_VAbs_Cplx(x, z->re);
  N_VScale(-1.0,x->im,z->im);
  N_VDiv(z->im,z->re,z->im);
  N_VDiv(x->re,z->re,z->re);
}

void N_VInv_Cplx(N_Vector x, N_Vector_Cplx z)
{  
  N_VInv(x,z->re);
  N_VConst(0.0,z->im);
}

void N_VAddConst_Cplx(N_Vector_Cplx x, cplx b, N_Vector_Cplx z)
{
  N_VAddConst(x->re,b.real(),z->re);
  N_VAddConst(x->im,b.imag(),z->im);
}

void N_VAddConst_Cplx(N_Vector x, cplx b, N_Vector_Cplx z)
{
  N_VAddConst(x,b.real(),z->re);
  N_VConst(b.imag(),z->im);
}

void N_VAddConst_Cplx(N_Vector_Cplx x, realtype b, N_Vector_Cplx z)
{
  N_VAddConst(x->re,b,z->re);
  N_VScale(1.0,x->im,z->im);
}

void N_VAddConst_Cplx(N_Vector x, realtype b, N_Vector_Cplx z)
{
  N_VAddConst(x,b,z->re);
  N_VConst(0.0,z->im);
}

cplx N_VDotProd_Cplx(N_Vector_Cplx x, N_Vector_Cplx y)
{
  realtype p1_re;
  realtype p2_re;
  
  realtype p1_im;
  realtype p2_im;
  
  p1_re=N_VDotProd(x->re, y->re);
  p2_re=N_VDotProd(x->im, y->im);
  
  p1_im=N_VDotProd(x->re, y->im);
  p2_im=N_VDotProd(x->im, y->re);
  
  cplx prod(p1_re-p2_re,p1_im+p2_im);
  
  return(prod);
}

cplx N_VDotProd_Cplx(N_Vector x, N_Vector_Cplx y)
{
  realtype p1_re;
  
  realtype p1_im;
  
  p1_re=N_VDotProd(x, y->re);  
  p1_im=N_VDotProd(x, y->im);
  
  cplx prod(p1_re,p1_im);
  
  return(prod);
}

cplx N_VDotProd_Cplx(N_Vector_Cplx x, N_Vector y)
{  
  return(N_VDotProd_Cplx(y,x));
}

void N_VPrint_Cplx(N_Vector_Cplx x)
{
  N_VPrint(x->re);
  N_VPrint(x->im);
}

