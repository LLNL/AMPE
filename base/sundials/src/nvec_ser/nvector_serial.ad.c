/************************** DISCLAIMER ********************************/
/*                                                                    */
/*   This file was generated on 06/22/05 15:09:57 by the version of   */
/*   ADIC compiled on  06/10/05 18:10:38                              */
/*                                                                    */
/*   ADIC was prepared as an account of work sponsored by an          */
/*   agency of the United States Government and the University of     */
/*   Chicago.  NEITHER THE AUTHOR(S), THE UNITED STATES GOVERNMENT    */
/*   NOR ANY AGENCY THEREOF, NOR THE UNIVERSITY OF CHICAGO, INCLUDING */
/*   ANY OF THEIR EMPLOYEES OR OFFICERS, MAKES ANY WARRANTY, EXPRESS  */
/*   OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY FOR */
/*   THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY INFORMATION OR  */
/*   PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE */
/*   PRIVATELY OWNED RIGHTS.                                          */
/*                                                                    */
/**********************************************************************/
#include "ad_deriv.h"
#include <stdio.h>
#include <stdlib.h>
#if !defined(AD_INCLUDE_nvector_46_h)
#define AD_INCLUDE_nvector_46_h
#include "nvector.ad.h"
#endif
struct ad__N_VectorContent_Serial {
long int  length;
int  own_data;
DERIV_TYPE  *data;
}
;
typedef struct ad__N_VectorContent_Serial *ad_N_VectorContent_Serial;
ad_N_Vector  ad_N_VNew_Serial(long int  vec_length);
ad_N_Vector  ad_N_VNewEmpty_Serial(long int  vec_length);
ad_N_Vector  ad_N_VMake_Serial(long int  vec_length,DERIV_TYPE  *v_data);
ad_N_Vector  *ad_N_VCloneVectorArray_Serial(int  count,ad_N_Vector  w);
ad_N_Vector  *ad_N_VCloneVectorArrayEmpty_Serial(int  count,ad_N_Vector  w);
void   ad_N_VDestroyVectorArray_Serial(ad_N_Vector  *vs,int  count);
void   ad_N_VPrint_Serial(ad_N_Vector  v);
ad_N_Vector  ad_N_VCloneEmpty_Serial(ad_N_Vector  w);
ad_N_Vector  ad_N_VClone_Serial(ad_N_Vector  w);
void   ad_N_VDestroy_Serial(ad_N_Vector  v);
void   ad_N_VSpace_Serial(ad_N_Vector  v,long int  *lrw,long int  *liw);
DERIV_TYPE  *ad_N_VGetArrayPointer_Serial(ad_N_Vector  v);
void   ad_N_VSetArrayPointer_Serial(DERIV_TYPE  *v_data,ad_N_Vector  v);
void   ad_N_VLinearSum_Serial(DERIV_TYPE  a,ad_N_Vector  x,DERIV_TYPE  b,ad_N_Vector  y,ad_N_Vector  z);
void   ad_N_VConst_Serial(DERIV_TYPE  c,ad_N_Vector  z);
void   ad_N_VProd_Serial(ad_N_Vector  x,ad_N_Vector  y,ad_N_Vector  z);
void   ad_N_VDiv_Serial(ad_N_Vector  x,ad_N_Vector  y,ad_N_Vector  z);
void   ad_N_VScale_Serial(DERIV_TYPE  c,ad_N_Vector  x,ad_N_Vector  z);
void   ad_N_VAbs_Serial(ad_N_Vector  x,ad_N_Vector  z);
void   ad_N_VInv_Serial(ad_N_Vector  x,ad_N_Vector  z);
void   ad_N_VAddConst_Serial(ad_N_Vector  x,DERIV_TYPE  b,ad_N_Vector  z);
void   ad_N_VDotProd_Serial(DERIV_TYPE  *ad_var_ret,ad_N_Vector  x,ad_N_Vector  y);
void   ad_N_VMaxNorm_Serial(DERIV_TYPE  *ad_var_ret,ad_N_Vector  x);
void   ad_N_VWrmsNorm_Serial(DERIV_TYPE  *ad_var_ret,ad_N_Vector  x,ad_N_Vector  w);
void   ad_N_VWrmsNormMask_Serial(DERIV_TYPE  *ad_var_ret,ad_N_Vector  x,ad_N_Vector  w,ad_N_Vector  id);
void   ad_N_VMin_Serial(DERIV_TYPE  *ad_var_ret,ad_N_Vector  x);
void   ad_N_VWL2Norm_Serial(DERIV_TYPE  *ad_var_ret,ad_N_Vector  x,ad_N_Vector  w);
void   ad_N_VL1Norm_Serial(DERIV_TYPE  *ad_var_ret,ad_N_Vector  x);
void   ad_N_VCompare_Serial(DERIV_TYPE  c,ad_N_Vector  x,ad_N_Vector  z);
int  ad_N_VInvTest_Serial(ad_N_Vector  x,ad_N_Vector  z);
int  ad_N_VConstrMask_Serial(ad_N_Vector  c,ad_N_Vector  x,ad_N_Vector  m);
void   ad_N_VMinQuotient_Serial(DERIV_TYPE  *ad_var_ret,ad_N_Vector  num,ad_N_Vector  denom);
void   ad_RPowerI(DERIV_TYPE  *ad_var_ret,DERIV_TYPE  base,int  exponent);
void   ad_RPowerR(DERIV_TYPE  *ad_var_ret,DERIV_TYPE  base,DERIV_TYPE  exponent);
void   ad_RSqrt(DERIV_TYPE  *ad_var_ret,DERIV_TYPE  x);
void   ad_RAbs(DERIV_TYPE  *ad_var_ret,DERIV_TYPE  x);
void   ad_RPower2(DERIV_TYPE  *ad_var_ret,DERIV_TYPE  x);
static void   ad_VCopy_Serial(ad_N_Vector  x,ad_N_Vector  z);
static void   ad_VSum_Serial(ad_N_Vector  x,ad_N_Vector  y,ad_N_Vector  z);
static void   ad_VDiff_Serial(ad_N_Vector  x,ad_N_Vector  y,ad_N_Vector  z);
static void   ad_VNeg_Serial(ad_N_Vector  x,ad_N_Vector  z);
static void   ad_VScaleSum_Serial(DERIV_TYPE  c,ad_N_Vector  x,ad_N_Vector  y,ad_N_Vector  z);
static void   ad_VScaleDiff_Serial(DERIV_TYPE  c,ad_N_Vector  x,ad_N_Vector  y,ad_N_Vector  z);
static void   ad_VLin1_Serial(DERIV_TYPE  a,ad_N_Vector  x,ad_N_Vector  y,ad_N_Vector  z);
static void   ad_VLin2_Serial(DERIV_TYPE  a,ad_N_Vector  x,ad_N_Vector  y,ad_N_Vector  z);
static void   ad_Vaxpy_Serial(DERIV_TYPE  a,ad_N_Vector  x,ad_N_Vector  y);
static void   ad_VScaleBy_Serial(DERIV_TYPE  a,ad_N_Vector  x);
ad_N_Vector  ad_N_VNewEmpty_Serial(long int  length) {
ad_N_Vector  v;
ad_N_Vector_Ops  ops;
ad_N_VectorContent_Serial  content;
void   *ad_var_0, *ad_var_1, *ad_var_2, *ad_var_3, *ad_var_5, *ad_var_6, *ad_var_7, *ad_var_9, *ad_var_10, *ad_var_11, *ad_var_13, *ad_var_14, *ad_var_15;
ad_N_Vector  ad_var_4;
ad_N_Vector_Ops  ad_var_8;
ad_N_VectorContent_Serial  ad_var_12;
    ad_var_0 = (void  * )0;
    v = ad_var_0;
    ad_var_1 = (void  * )0;
    ops = ad_var_1;
    ad_var_2 = (void  * )0;
    content = ad_var_2;
    ad_var_3 = malloc( sizeof *v);
    ad_var_4 = (ad_N_Vector )ad_var_3;
    v = ad_var_4;
    ad_var_5 = (void  * )0;
    if (v == ad_var_5)     {
        ad_var_6 = (void  * )0;
        return ad_var_6;
    }
    ad_var_7 = malloc( sizeof (struct ad__generic_N_Vector_Ops ));
    ad_var_8 = (ad_N_Vector_Ops )ad_var_7;
    ops = ad_var_8;
    ad_var_9 = (void  * )0;
    if (ops == ad_var_9)     {
        free(v);
        ad_var_10 = (void  * )0;
        return ad_var_10;
    }
    ops->nvclone = ad_N_VClone_Serial();
    ops->nvcloneempty = ad_N_VCloneEmpty_Serial();
    ops->nvdestroy = ad_N_VDestroy_Serial();
    ops->nvspace = ad_N_VSpace_Serial();
    ops->nvgetarraypointer = ad_N_VGetArrayPointer_Serial();
    ops->nvsetarraypointer = ad_N_VSetArrayPointer_Serial();
    ops->nvlinearsum = ad_N_VLinearSum_Serial();
    ops->nvconst = ad_N_VConst_Serial();
    ops->nvprod = ad_N_VProd_Serial();
    ops->nvdiv = ad_N_VDiv_Serial();
    ops->nvscale = ad_N_VScale_Serial();
    ops->nvabs = ad_N_VAbs_Serial();
    ops->nvinv = ad_N_VInv_Serial();
    ops->nvaddconst = ad_N_VAddConst_Serial();
    ops->nvdotprod = ad_N_VDotProd_Serial(DERIV_TYPE  *ad_var_ret);
    ops->nvmaxnorm = ad_N_VMaxNorm_Serial(DERIV_TYPE  *ad_var_ret);
    ops->nvwrmsnormmask = ad_N_VWrmsNormMask_Serial(DERIV_TYPE  *ad_var_ret);
    ops->nvwrmsnorm = ad_N_VWrmsNorm_Serial(DERIV_TYPE  *ad_var_ret);
    ops->nvmin = ad_N_VMin_Serial(DERIV_TYPE  *ad_var_ret);
    ops->nvwl2norm = ad_N_VWL2Norm_Serial(DERIV_TYPE  *ad_var_ret);
    ops->nvl1norm = ad_N_VL1Norm_Serial(DERIV_TYPE  *ad_var_ret);
    ops->nvcompare = ad_N_VCompare_Serial();
    ops->nvinvtest = ad_N_VInvTest_Serial();
    ops->nvconstrmask = ad_N_VConstrMask_Serial();
    ops->nvminquotient = ad_N_VMinQuotient_Serial(DERIV_TYPE  *ad_var_ret);
    ad_var_11 = malloc( sizeof (struct ad__N_VectorContent_Serial ));
    ad_var_12 = (ad_N_VectorContent_Serial )ad_var_11;
    content = ad_var_12;
    ad_var_13 = (void  * )0;
    if (content == ad_var_13)     {
        free(ops);
        free(v);
        ad_var_14 = (void  * )0;
        return ad_var_14;
    }
    content->length = length;
    content->own_data = 0;
    ad_var_15 = (void  * )0;
    content->data = ad_var_15;
    v->content = content;
    v->ops = ops;
    return v;
}
ad_N_Vector  ad_N_VNew_Serial(long int  length) {
ad_N_Vector  v;
DERIV_TYPE  *data;
void   *ad_var_0, *ad_var_1, *ad_var_3, *ad_var_4, *ad_var_5, *ad_var_7, *ad_var_8;
ad_N_Vector  ad_var_2;
DERIV_TYPE  *ad_var_6;
ad_N_VectorContent_Serial  ad_var_9, ad_var_10;
    ad_var_0 = (void  * )0;
    v = ad_var_0;
    ad_var_1 = (void  * )0;
    data = ad_var_1;
    ad_var_2 = ad_N_VNewEmpty_Serial(length);
    v = ad_var_2;
    ad_var_3 = (void  * )0;
    if (v == ad_var_3)     {
        ad_var_4 = (void  * )0;
        return ad_var_4;
    }
    if (length > 0)     {
        ad_var_5 = malloc(length *  sizeof (DERIV_TYPE ));
        ad_var_6 = (DERIV_TYPE * )ad_var_5;
        data = ad_var_6;
        ad_var_7 = (void  * )0;
        if (data == ad_var_7)         {
            ad_N_VDestroy_Serial(v);
            ad_var_8 = (void  * )0;
            return ad_var_8;
        }
        ad_var_9 = (ad_N_VectorContent_Serial )v->content;
        ad_var_9->own_data = 1;
        ad_var_10 = (ad_N_VectorContent_Serial )v->content;
        ad_var_10->data = data;
    }
    return v;
}
ad_N_Vector  ad_N_VMake_Serial(long int  length,DERIV_TYPE  *v_data) {
ad_N_Vector  v;
void   *ad_var_0, *ad_var_2, *ad_var_3;
ad_N_Vector  ad_var_1;
ad_N_VectorContent_Serial  ad_var_4, ad_var_5;
    ad_var_0 = (void  * )0;
    v = ad_var_0;
    ad_var_1 = ad_N_VNewEmpty_Serial(length);
    v = ad_var_1;
    ad_var_2 = (void  * )0;
    if (v == ad_var_2)     {
        ad_var_3 = (void  * )0;
        return ad_var_3;
    }
    if (length > 0)     {
        ad_var_4 = (ad_N_VectorContent_Serial )v->content;
        ad_var_4->own_data = 0;
        ad_var_5 = (ad_N_VectorContent_Serial )v->content;
        ad_var_5->data = v_data;
    }
    return v;
}
ad_N_Vector  *ad_N_VCloneVectorArray_Serial(int  count,ad_N_Vector  w) {
ad_N_Vector  *vs;
int  j;
void   *ad_var_0, *ad_var_1, *ad_var_2, *ad_var_4, *ad_var_5, *ad_var_8, *ad_var_9;
ad_N_Vector  *ad_var_3, ad_var_7;
int  ad_var_6;
    ad_var_0 = (void  * )0;
    vs = ad_var_0;
    if (count <= 0)     {
        ad_var_1 = (void  * )0;
        return ad_var_1;
    }
    ad_var_2 = malloc(count *  sizeof (ad_N_Vector ));
    ad_var_3 = (ad_N_Vector * )ad_var_2;
    vs = ad_var_3;
    ad_var_4 = (void  * )0;
    if (vs == ad_var_4)     {
        ad_var_5 = (void  * )0;
        return ad_var_5;
    }
    for (j = 0;     j < count;     )    {
        ad_var_7 = ad_N_VClone_Serial(w);
        vs[j] = ad_var_7;
        ad_var_8 = (void  * )0;
        if (vs[j] == ad_var_8)         {
            ad_N_VDestroyVectorArray_Serial(vs, j - 1);
            ad_var_9 = (void  * )0;
            return ad_var_9;
        }
        ad_var_6 = j++;
    }
    return vs;
}
ad_N_Vector  *ad_N_VCloneVectorArrayEmpty_Serial(int  count,ad_N_Vector  w) {
ad_N_Vector  *vs;
int  j;
void   *ad_var_0, *ad_var_1, *ad_var_2, *ad_var_4, *ad_var_5, *ad_var_8, *ad_var_9;
ad_N_Vector  *ad_var_3, ad_var_7;
int  ad_var_6;
    ad_var_0 = (void  * )0;
    vs = ad_var_0;
    if (count <= 0)     {
        ad_var_1 = (void  * )0;
        return ad_var_1;
    }
    ad_var_2 = malloc(count *  sizeof (ad_N_Vector ));
    ad_var_3 = (ad_N_Vector * )ad_var_2;
    vs = ad_var_3;
    ad_var_4 = (void  * )0;
    if (vs == ad_var_4)     {
        ad_var_5 = (void  * )0;
        return ad_var_5;
    }
    for (j = 0;     j < count;     )    {
        ad_var_7 = ad_N_VCloneEmpty_Serial(w);
        vs[j] = ad_var_7;
        ad_var_8 = (void  * )0;
        if (vs[j] == ad_var_8)         {
            ad_N_VDestroyVectorArray_Serial(vs, j - 1);
            ad_var_9 = (void  * )0;
            return ad_var_9;
        }
        ad_var_6 = j++;
    }
    return vs;
}
void   ad_N_VDestroyVectorArray_Serial(ad_N_Vector  *vs,int  count) {
int  j;
int  ad_var_0;
    for (j = 0;     j < count;     )    {
        ad_N_VDestroy_Serial(vs[j]);
        ad_var_0 = j++;
    }
    free(vs);
    return;
}
void   ad_N_VPrint_Serial(ad_N_Vector  x) {
long int  i, N;
DERIV_TYPE  *xd;
void   *ad_var_0;
ad_N_VectorContent_Serial  ad_var_1, ad_var_2;
long int  ad_var_3;
    ad_var_0 = (void  * )0;
    xd = ad_var_0;
    ad_var_1 = (ad_N_VectorContent_Serial )x->content;
    N = ad_var_1->length;
    ad_var_2 = (ad_N_VectorContent_Serial )x->content;
    xd = ad_var_2->data;
    for (i = 0;     i < N;     )    {
        printf("%11.8lg\n", xd[i]);
        ad_var_3 = i++;
    }
    printf("\n");
    return;
}
ad_N_Vector  ad_N_VCloneEmpty_Serial(ad_N_Vector  w) {
ad_N_Vector  v;
ad_N_Vector_Ops  ops;
ad_N_VectorContent_Serial  content;
void   *ad_var_0, *ad_var_1, *ad_var_2, *ad_var_3, *ad_var_4, *ad_var_5, *ad_var_7, *ad_var_8, *ad_var_9, *ad_var_11, *ad_var_12, *ad_var_13, *ad_var_15, *ad_var_16, *ad_var_18;
ad_N_Vector  ad_var_6;
ad_N_Vector_Ops  ad_var_10;
ad_N_VectorContent_Serial  ad_var_14, ad_var_17;
    ad_var_0 = (void  * )0;
    v = ad_var_0;
    ad_var_1 = (void  * )0;
    ops = ad_var_1;
    ad_var_2 = (void  * )0;
    content = ad_var_2;
    ad_var_3 = (void  * )0;
    if (w == ad_var_3)     {
        ad_var_4 = (void  * )0;
        return ad_var_4;
    }
    ad_var_5 = malloc( sizeof *v);
    ad_var_6 = (ad_N_Vector )ad_var_5;
    v = ad_var_6;
    ad_var_7 = (void  * )0;
    if (v == ad_var_7)     {
        ad_var_8 = (void  * )0;
        return ad_var_8;
    }
    ad_var_9 = malloc( sizeof (struct ad__generic_N_Vector_Ops ));
    ad_var_10 = (ad_N_Vector_Ops )ad_var_9;
    ops = ad_var_10;
    ad_var_11 = (void  * )0;
    if (ops == ad_var_11)     {
        free(v);
        ad_var_12 = (void  * )0;
        return ad_var_12;
    }
    ops->nvclone = (w->ops)->nvclone;
    ops->nvcloneempty = (w->ops)->nvcloneempty;
    ops->nvdestroy = (w->ops)->nvdestroy;
    ops->nvspace = (w->ops)->nvspace;
    ops->nvgetarraypointer = (w->ops)->nvgetarraypointer;
    ops->nvsetarraypointer = (w->ops)->nvsetarraypointer;
    ops->nvlinearsum = (w->ops)->nvlinearsum;
    ops->nvconst = (w->ops)->nvconst;
    ops->nvprod = (w->ops)->nvprod;
    ops->nvdiv = (w->ops)->nvdiv;
    ops->nvscale = (w->ops)->nvscale;
    ops->nvabs = (w->ops)->nvabs;
    ops->nvinv = (w->ops)->nvinv;
    ops->nvaddconst = (w->ops)->nvaddconst;
    ops->nvdotprod = (w->ops)->nvdotprod;
    ops->nvmaxnorm = (w->ops)->nvmaxnorm;
    ops->nvwrmsnormmask = (w->ops)->nvwrmsnormmask;
    ops->nvwrmsnorm = (w->ops)->nvwrmsnorm;
    ops->nvmin = (w->ops)->nvmin;
    ops->nvwl2norm = (w->ops)->nvwl2norm;
    ops->nvl1norm = (w->ops)->nvl1norm;
    ops->nvcompare = (w->ops)->nvcompare;
    ops->nvinvtest = (w->ops)->nvinvtest;
    ops->nvconstrmask = (w->ops)->nvconstrmask;
    ops->nvminquotient = (w->ops)->nvminquotient;
    ad_var_13 = malloc( sizeof (struct ad__N_VectorContent_Serial ));
    ad_var_14 = (ad_N_VectorContent_Serial )ad_var_13;
    content = ad_var_14;
    ad_var_15 = (void  * )0;
    if (content == ad_var_15)     {
        free(ops);
        free(v);
        ad_var_16 = (void  * )0;
        return ad_var_16;
    }
    ad_var_17 = (ad_N_VectorContent_Serial )w->content;
    content->length = ad_var_17->length;
    content->own_data = 0;
    ad_var_18 = (void  * )0;
    content->data = ad_var_18;
    v->content = content;
    v->ops = ops;
    return v;
}
ad_N_Vector  ad_N_VClone_Serial(ad_N_Vector  w) {
ad_N_Vector  v;
DERIV_TYPE  *data;
long int  length;
void   *ad_var_0, *ad_var_1, *ad_var_3, *ad_var_4, *ad_var_6, *ad_var_8, *ad_var_9;
ad_N_Vector  ad_var_2;
ad_N_VectorContent_Serial  ad_var_5, ad_var_10, ad_var_11;
DERIV_TYPE  *ad_var_7;
    ad_var_0 = (void  * )0;
    v = ad_var_0;
    ad_var_1 = (void  * )0;
    data = ad_var_1;
    ad_var_2 = ad_N_VCloneEmpty_Serial(w);
    v = ad_var_2;
    ad_var_3 = (void  * )0;
    if (v == ad_var_3)     {
        ad_var_4 = (void  * )0;
        return ad_var_4;
    }
    ad_var_5 = (ad_N_VectorContent_Serial )w->content;
    length = ad_var_5->length;
    if (length > 0)     {
        ad_var_6 = malloc(length *  sizeof (DERIV_TYPE ));
        ad_var_7 = (DERIV_TYPE * )ad_var_6;
        data = ad_var_7;
        ad_var_8 = (void  * )0;
        if (data == ad_var_8)         {
            ad_N_VDestroy_Serial(v);
            ad_var_9 = (void  * )0;
            return ad_var_9;
        }
        ad_var_10 = (ad_N_VectorContent_Serial )v->content;
        ad_var_10->own_data = 1;
        ad_var_11 = (ad_N_VectorContent_Serial )v->content;
        ad_var_11->data = data;
    }
    return v;
}
void   ad_N_VDestroy_Serial(ad_N_Vector  v) {
ad_N_VectorContent_Serial  ad_var_0, ad_var_1;
    ad_var_0 = (ad_N_VectorContent_Serial )v->content;
    if (ad_var_0->own_data == 1)     {
        ad_var_1 = (ad_N_VectorContent_Serial )v->content;
        free(ad_var_1->data);
    }
    free(v->content);
    free(v->ops);
    free(v);
    return;
}
void   ad_N_VSpace_Serial(ad_N_Vector  v,long int  *lrw,long int  *liw) {
ad_N_VectorContent_Serial  ad_var_0;
    ad_var_0 = (ad_N_VectorContent_Serial )v->content;
    *lrw = ad_var_0->length;
    *liw = 1;
    return;
}
DERIV_TYPE  *ad_N_VGetArrayPointer_Serial(ad_N_Vector  v) {
ad_N_VectorContent_Serial  ad_var_0;
DERIV_TYPE  *ad_var_1;
    ad_var_0 = (ad_N_VectorContent_Serial )v->content;
    ad_var_1 = (DERIV_TYPE * )ad_var_0->data;
    return ad_var_1;
}
void   ad_N_VSetArrayPointer_Serial(DERIV_TYPE  *v_data,ad_N_Vector  v) {
ad_N_VectorContent_Serial  ad_var_0, ad_var_1;
    ad_var_0 = (ad_N_VectorContent_Serial )v->content;
    if (ad_var_0->length > 0)     {
        ad_var_1 = (ad_N_VectorContent_Serial )v->content;
        ad_var_1->data = v_data;
    }
    return;
}
void   ad_N_VLinearSum_Serial(DERIV_TYPE  a,ad_N_Vector  x,DERIV_TYPE  b,ad_N_Vector  y,ad_N_Vector  z) {
long int  i, N;
DERIV_TYPE  c, *xd, *yd, *zd;
ad_N_Vector  v1, v2;
int  test;
void   *ad_var_0;
DERIV_TYPE  ad_var_1, ad_var_2;
ad_N_VectorContent_Serial  ad_var_3, ad_var_4, ad_var_5, ad_var_6;
long int  ad_var_7;
    double  ad_loc_0;
    double  ad_loc_1;
    double  ad_loc_2;
    ad_var_0 = (void  * )0;
    xd = yd = zd = ad_var_0;
    if (DERIV_val(b) == 1.0 && z == y)     {
        ad_Vaxpy_Serial(a, x, y);
        return;
    }
    if (DERIV_val(a) == 1.0 && z == x)     {
        ad_Vaxpy_Serial(b, y, x);
        return;
    }
    if (DERIV_val(a) == 1.0 && DERIV_val(b) == 1.0)     {
        ad_VSum_Serial(x, y, z);
        return;
    }
    if ((test = DERIV_val(a) == 1.0 && DERIV_val(b) ==  -1.0) || DERIV_val(a) ==  -1.0 && DERIV_val(b) == 1.0)     {
        v1 = test? y : x;
        v2 = test? x : y;
        ad_VDiff_Serial(v2, v1, z);
        return;
    }
    if ((test = DERIV_val(a) == 1.0) || DERIV_val(b) == 1.0)     {
        if (test)         {
            {
                ad_grad_axpy_copy(&(ad_var_1), &(b));
                DERIV_val(ad_var_1) = DERIV_val(b);
            }
        }
        else         {
            {
                ad_grad_axpy_copy(&(ad_var_1), &(a));
                DERIV_val(ad_var_1) = DERIV_val(a);
            }
        }
        {
            ad_grad_axpy_copy(&(c), &(ad_var_1));
            DERIV_val(c) = DERIV_val(ad_var_1);
        }
        v1 = test? y : x;
        v2 = test? x : y;
        ad_VLin1_Serial(c, v1, v2, z);
        return;
    }
    if ((test = DERIV_val(a) ==  -1.0) || DERIV_val(b) ==  -1.0)     {
        if (test)         {
            {
                ad_grad_axpy_copy(&(ad_var_2), &(b));
                DERIV_val(ad_var_2) = DERIV_val(b);
            }
        }
        else         {
            {
                ad_grad_axpy_copy(&(ad_var_2), &(a));
                DERIV_val(ad_var_2) = DERIV_val(a);
            }
        }
        {
            ad_grad_axpy_copy(&(c), &(ad_var_2));
            DERIV_val(c) = DERIV_val(ad_var_2);
        }
        v1 = test? y : x;
        v2 = test? x : y;
        ad_VLin2_Serial(c, v1, v2, z);
        return;
    }
    if (DERIV_val(a) == DERIV_val(b))     {
        ad_VScaleSum_Serial(a, x, y, z);
        return;
    }
    if (DERIV_val(a) ==  -DERIV_val(b))     {
        ad_VScaleDiff_Serial(a, x, y, z);
        return;
    }
    ad_var_3 = (ad_N_VectorContent_Serial )x->content;
    N = ad_var_3->length;
    ad_var_4 = (ad_N_VectorContent_Serial )x->content;
    xd = ad_var_4->data;
    ad_var_5 = (ad_N_VectorContent_Serial )y->content;
    yd = ad_var_5->data;
    ad_var_6 = (ad_N_VectorContent_Serial )z->content;
    zd = ad_var_6->data;
    for (i = 0;     i < N;     )    {
        {
            ad_loc_0 = DERIV_val(a) * DERIV_val(xd[i]);
            ad_loc_1 = DERIV_val(b) * DERIV_val(yd[i]);
            ad_loc_2 = ad_loc_0 + ad_loc_1;
            ad_grad_axpy_4(&(zd[i]), DERIV_val(xd[i]), &(a), DERIV_val(a), &(xd[i]), DERIV_val(yd[i]), &(b), DERIV_val(b), &(yd[i]));
            DERIV_val(zd[i]) = ad_loc_2;
        }
        ad_var_7 = i++;
    }
    return;
}
void   ad_N_VConst_Serial(DERIV_TYPE  c,ad_N_Vector  z) {
long int  i, N;
DERIV_TYPE  *zd;
void   *ad_var_0;
ad_N_VectorContent_Serial  ad_var_1, ad_var_2;
long int  ad_var_3;
    ad_var_0 = (void  * )0;
    zd = ad_var_0;
    ad_var_1 = (ad_N_VectorContent_Serial )z->content;
    N = ad_var_1->length;
    ad_var_2 = (ad_N_VectorContent_Serial )z->content;
    zd = ad_var_2->data;
    for (i = 0;     i < N;     )    {
        {
            ad_grad_axpy_copy(&(zd[i]), &(c));
            DERIV_val(zd[i]) = DERIV_val(c);
        }
        ad_var_3 = i++;
    }
    return;
}
void   ad_N_VProd_Serial(ad_N_Vector  x,ad_N_Vector  y,ad_N_Vector  z) {
long int  i, N;
DERIV_TYPE  *xd, *yd, *zd;
void   *ad_var_0;
ad_N_VectorContent_Serial  ad_var_1, ad_var_2, ad_var_3, ad_var_4;
long int  ad_var_5;
    double  ad_loc_0;
    ad_var_0 = (void  * )0;
    xd = yd = zd = ad_var_0;
    ad_var_1 = (ad_N_VectorContent_Serial )x->content;
    N = ad_var_1->length;
    ad_var_2 = (ad_N_VectorContent_Serial )x->content;
    xd = ad_var_2->data;
    ad_var_3 = (ad_N_VectorContent_Serial )y->content;
    yd = ad_var_3->data;
    ad_var_4 = (ad_N_VectorContent_Serial )z->content;
    zd = ad_var_4->data;
    for (i = 0;     i < N;     )    {
        {
            ad_loc_0 = DERIV_val(xd[i]) * DERIV_val(yd[i]);
            ad_grad_axpy_2(&(zd[i]), DERIV_val(yd[i]), &(xd[i]), DERIV_val(xd[i]), &(yd[i]));
            DERIV_val(zd[i]) = ad_loc_0;
        }
        ad_var_5 = i++;
    }
    return;
}
void   ad_N_VDiv_Serial(ad_N_Vector  x,ad_N_Vector  y,ad_N_Vector  z) {
long int  i, N;
DERIV_TYPE  *xd, *yd, *zd;
void   *ad_var_0;
ad_N_VectorContent_Serial  ad_var_1, ad_var_2, ad_var_3, ad_var_4;
long int  ad_var_5;
    double  ad_loc_0;
    double  ad_adj_0;
    double  ad_adj_1;
    ad_var_0 = (void  * )0;
    xd = yd = zd = ad_var_0;
    ad_var_1 = (ad_N_VectorContent_Serial )x->content;
    N = ad_var_1->length;
    ad_var_2 = (ad_N_VectorContent_Serial )x->content;
    xd = ad_var_2->data;
    ad_var_3 = (ad_N_VectorContent_Serial )y->content;
    yd = ad_var_3->data;
    ad_var_4 = (ad_N_VectorContent_Serial )z->content;
    zd = ad_var_4->data;
    for (i = 0;     i < N;     )    {
        {
            ad_loc_0 = DERIV_val(xd[i]) / DERIV_val(yd[i]);
            ad_adj_0 =  -ad_loc_0 / DERIV_val(yd[i]);
            ad_adj_1 = 1.000000000000000e+00 / DERIV_val(yd[i]);
            ad_grad_axpy_2(&(zd[i]), ad_adj_1, &(xd[i]), ad_adj_0, &(yd[i]));
            DERIV_val(zd[i]) = ad_loc_0;
        }
        ad_var_5 = i++;
    }
    return;
}
void   ad_N_VScale_Serial(DERIV_TYPE  c,ad_N_Vector  x,ad_N_Vector  z) {
long int  i, N;
DERIV_TYPE  *xd, *zd;
void   *ad_var_0;
ad_N_VectorContent_Serial  ad_var_1, ad_var_2, ad_var_3;
long int  ad_var_4;
    double  ad_loc_0;
    ad_var_0 = (void  * )0;
    xd = zd = ad_var_0;
    if (z == x)     {
        ad_VScaleBy_Serial(c, x);
        return;
    }
    if (DERIV_val(c) == 1.0)     {
        ad_VCopy_Serial(x, z);
    }
    else     {
        if (DERIV_val(c) ==  -1.0)         {
            ad_VNeg_Serial(x, z);
        }
        else         {
            ad_var_1 = (ad_N_VectorContent_Serial )x->content;
            N = ad_var_1->length;
            ad_var_2 = (ad_N_VectorContent_Serial )x->content;
            xd = ad_var_2->data;
            ad_var_3 = (ad_N_VectorContent_Serial )z->content;
            zd = ad_var_3->data;
            for (i = 0;             i < N;             )            {
                {
                    ad_loc_0 = DERIV_val(c) * DERIV_val(xd[i]);
                    ad_grad_axpy_2(&(zd[i]), DERIV_val(xd[i]), &(c), DERIV_val(c), &(xd[i]));
                    DERIV_val(zd[i]) = ad_loc_0;
                }
                ad_var_4 = i++;
            }
        }
    }
    return;
}
void   ad_N_VAbs_Serial(ad_N_Vector  x,ad_N_Vector  z) {
long int  i, N;
DERIV_TYPE  *xd, *zd;
void   *ad_var_0;
ad_N_VectorContent_Serial  ad_var_1, ad_var_2, ad_var_3;
long int  ad_var_4;
    ad_var_0 = (void  * )0;
    xd = zd = ad_var_0;
    ad_var_1 = (ad_N_VectorContent_Serial )x->content;
    N = ad_var_1->length;
    ad_var_2 = (ad_N_VectorContent_Serial )x->content;
    xd = ad_var_2->data;
    ad_var_3 = (ad_N_VectorContent_Serial )z->content;
    zd = ad_var_3->data;
    for (i = 0;     i < N;     )    {
        ad_RAbs( &zd[i], xd[i]);
        ad_var_4 = i++;
    }
    return;
}
void   ad_N_VInv_Serial(ad_N_Vector  x,ad_N_Vector  z) {
long int  i, N;
DERIV_TYPE  *xd, *zd;
void   *ad_var_0;
ad_N_VectorContent_Serial  ad_var_1, ad_var_2, ad_var_3;
long int  ad_var_4;
    double  ad_loc_0;
    double  ad_adj_0;
    ad_var_0 = (void  * )0;
    xd = zd = ad_var_0;
    ad_var_1 = (ad_N_VectorContent_Serial )x->content;
    N = ad_var_1->length;
    ad_var_2 = (ad_N_VectorContent_Serial )x->content;
    xd = ad_var_2->data;
    ad_var_3 = (ad_N_VectorContent_Serial )z->content;
    zd = ad_var_3->data;
    for (i = 0;     i < N;     )    {
        {
            ad_loc_0 = 1.0 / DERIV_val(xd[i]);
            ad_adj_0 =  -ad_loc_0 / DERIV_val(xd[i]);
            ad_grad_axpy_1(&(zd[i]), ad_adj_0, &(xd[i]));
            DERIV_val(zd[i]) = ad_loc_0;
        }
        ad_var_4 = i++;
    }
    return;
}
void   ad_N_VAddConst_Serial(ad_N_Vector  x,DERIV_TYPE  b,ad_N_Vector  z) {
long int  i, N;
DERIV_TYPE  *xd, *zd;
void   *ad_var_0;
ad_N_VectorContent_Serial  ad_var_1, ad_var_2, ad_var_3;
long int  ad_var_4;
    double  ad_loc_0;
    ad_var_0 = (void  * )0;
    xd = zd = ad_var_0;
    ad_var_1 = (ad_N_VectorContent_Serial )x->content;
    N = ad_var_1->length;
    ad_var_2 = (ad_N_VectorContent_Serial )x->content;
    xd = ad_var_2->data;
    ad_var_3 = (ad_N_VectorContent_Serial )z->content;
    zd = ad_var_3->data;
    for (i = 0;     i < N;     )    {
        {
            ad_loc_0 = DERIV_val(xd[i]) + DERIV_val(b);
            ad_grad_axpy_2(&(zd[i]), 1.000000000000000e+00, &(xd[i]), 1.000000000000000e+00, &(b));
            DERIV_val(zd[i]) = ad_loc_0;
        }
        ad_var_4 = i++;
    }
    return;
}
void   ad_N_VDotProd_Serial(DERIV_TYPE  *ad_var_ret,ad_N_Vector  x,ad_N_Vector  y) {
long int  i, N;
DERIV_TYPE  sum, *xd, *yd;
void   *ad_var_0;
ad_N_VectorContent_Serial  ad_var_1, ad_var_2, ad_var_3;
long int  ad_var_4;
    double  ad_loc_0;
    double  ad_loc_1;
    {
        ad_grad_axpy_0(&(sum));
        DERIV_val(sum) = 0.0;
    }
    ad_var_0 = (void  * )0;
    xd = yd = ad_var_0;
    ad_var_1 = (ad_N_VectorContent_Serial )x->content;
    N = ad_var_1->length;
    ad_var_2 = (ad_N_VectorContent_Serial )x->content;
    xd = ad_var_2->data;
    ad_var_3 = (ad_N_VectorContent_Serial )y->content;
    yd = ad_var_3->data;
    for (i = 0;     i < N;     )    {
        {
            ad_loc_0 = DERIV_val(xd[i]) * DERIV_val(yd[i]);
            ad_loc_1 = DERIV_val(sum) + ad_loc_0;
            ad_grad_axpy_3(&(sum), 1.000000000000000e+00, &(sum), DERIV_val(yd[i]), &(xd[i]), DERIV_val(xd[i]), &(yd[i]));
            DERIV_val(sum) = ad_loc_1;
        }
        ad_var_4 = i++;
    }
    {
        ad_grad_axpy_copy(&(*ad_var_ret), &(sum));
        DERIV_val(*ad_var_ret) = DERIV_val(sum);
    }
    return;
}
void   ad_N_VMaxNorm_Serial(DERIV_TYPE  *ad_var_ret,ad_N_Vector  x) {
long int  i, N;
DERIV_TYPE  max, *xd;
void   *ad_var_0;
ad_N_VectorContent_Serial  ad_var_1, ad_var_2;
long int  ad_var_3;
DERIV_TYPE  ad_var_4;
    {
        ad_grad_axpy_0(&(max));
        DERIV_val(max) = 0.0;
    }
    ad_var_0 = (void  * )0;
    xd = ad_var_0;
    ad_var_1 = (ad_N_VectorContent_Serial )x->content;
    N = ad_var_1->length;
    ad_var_2 = (ad_N_VectorContent_Serial )x->content;
    xd = ad_var_2->data;
    for (i = 0;     i < N;     )    {
        ad_RAbs( &ad_var_4, xd[i]);
        if (DERIV_val(ad_var_4) > DERIV_val(max))         {
            ad_RAbs( &max, xd[i]);
        }
        ad_var_3 = i++;
    }
    {
        ad_grad_axpy_copy(&(*ad_var_ret), &(max));
        DERIV_val(*ad_var_ret) = DERIV_val(max);
    }
    return;
}
void   ad_N_VWrmsNorm_Serial(DERIV_TYPE  *ad_var_ret,ad_N_Vector  x,ad_N_Vector  w) {
long int  i, N;
DERIV_TYPE  sum, prodi, *xd, *wd;
void   *ad_var_0;
ad_N_VectorContent_Serial  ad_var_1, ad_var_2, ad_var_3;
long int  ad_var_4;
DERIV_TYPE  ad_var_5, ad_var_6, ad_var_7;
    double  ad_loc_0;
    double  ad_adj_0;
    {
        ad_grad_axpy_0(&(sum));
        DERIV_val(sum) = 0.0;
    }
    ad_var_0 = (void  * )0;
    xd = wd = ad_var_0;
    ad_var_1 = (ad_N_VectorContent_Serial )x->content;
    N = ad_var_1->length;
    ad_var_2 = (ad_N_VectorContent_Serial )x->content;
    xd = ad_var_2->data;
    ad_var_3 = (ad_N_VectorContent_Serial )w->content;
    wd = ad_var_3->data;
    for (i = 0;     i < N;     )    {
        {
            ad_loc_0 = DERIV_val(xd[i]) * DERIV_val(wd[i]);
            ad_grad_axpy_2(&(prodi), DERIV_val(wd[i]), &(xd[i]), DERIV_val(xd[i]), &(wd[i]));
            DERIV_val(prodi) = ad_loc_0;
        }
        ad_RPower2( &ad_var_5, prodi);
        {
            ad_loc_0 = DERIV_val(sum) + DERIV_val(ad_var_5);
            ad_grad_axpy_2(&(sum), 1.000000000000000e+00, &(sum), 1.000000000000000e+00, &(ad_var_5));
            DERIV_val(sum) = ad_loc_0;
        }
        ad_var_4 = i++;
    }
    {
        ad_loc_0 = DERIV_val(sum) / N;
        ad_adj_0 = 1.000000000000000e+00 / N;
        ad_grad_axpy_1(&(ad_var_6), ad_adj_0, &(sum));
        DERIV_val(ad_var_6) = ad_loc_0;
    }
    ad_RSqrt( &ad_var_7, ad_var_6);
    {
        ad_grad_axpy_copy(&(*ad_var_ret), &(ad_var_7));
        DERIV_val(*ad_var_ret) = DERIV_val(ad_var_7);
    }
    return;
}
void   ad_N_VWrmsNormMask_Serial(DERIV_TYPE  *ad_var_ret,ad_N_Vector  x,ad_N_Vector  w,ad_N_Vector  id) {
long int  i, N;
DERIV_TYPE  sum, prodi, *xd, *wd, *idd;
void   *ad_var_0;
ad_N_VectorContent_Serial  ad_var_1, ad_var_2, ad_var_3, ad_var_4;
long int  ad_var_5;
DERIV_TYPE  ad_var_6, ad_var_7, ad_var_8;
    double  ad_loc_0;
    double  ad_adj_0;
    {
        ad_grad_axpy_0(&(sum));
        DERIV_val(sum) = 0.0;
    }
    ad_var_0 = (void  * )0;
    xd = wd = idd = ad_var_0;
    ad_var_1 = (ad_N_VectorContent_Serial )x->content;
    N = ad_var_1->length;
    ad_var_2 = (ad_N_VectorContent_Serial )x->content;
    xd = ad_var_2->data;
    ad_var_3 = (ad_N_VectorContent_Serial )w->content;
    wd = ad_var_3->data;
    ad_var_4 = (ad_N_VectorContent_Serial )id->content;
    idd = ad_var_4->data;
    for (i = 0;     i < N;     )    {
        if (DERIV_val(idd[i]) > 0.0)         {
            {
                {
                    ad_loc_0 = DERIV_val(xd[i]) * DERIV_val(wd[i]);
                    ad_grad_axpy_2(&(prodi), DERIV_val(wd[i]), &(xd[i]), DERIV_val(xd[i]), &(wd[i]));
                    DERIV_val(prodi) = ad_loc_0;
                }
                ad_RPower2( &ad_var_6, prodi);
                {
                    ad_loc_0 = DERIV_val(sum) + DERIV_val(ad_var_6);
                    ad_grad_axpy_2(&(sum), 1.000000000000000e+00, &(sum), 1.000000000000000e+00, &(ad_var_6));
                    DERIV_val(sum) = ad_loc_0;
                }
            }
        }
        ad_var_5 = i++;
    }
    {
        ad_loc_0 = DERIV_val(sum) / N;
        ad_adj_0 = 1.000000000000000e+00 / N;
        ad_grad_axpy_1(&(ad_var_7), ad_adj_0, &(sum));
        DERIV_val(ad_var_7) = ad_loc_0;
    }
    ad_RSqrt( &ad_var_8, ad_var_7);
    {
        ad_grad_axpy_copy(&(*ad_var_ret), &(ad_var_8));
        DERIV_val(*ad_var_ret) = DERIV_val(ad_var_8);
    }
    return;
}
void   ad_N_VMin_Serial(DERIV_TYPE  *ad_var_ret,ad_N_Vector  x) {
long int  i, N;
DERIV_TYPE  min, *xd;
void   *ad_var_0;
ad_N_VectorContent_Serial  ad_var_1, ad_var_2;
long int  ad_var_3;
    ad_var_0 = (void  * )0;
    xd = ad_var_0;
    ad_var_1 = (ad_N_VectorContent_Serial )x->content;
    N = ad_var_1->length;
    ad_var_2 = (ad_N_VectorContent_Serial )x->content;
    xd = ad_var_2->data;
    {
        ad_grad_axpy_copy(&(min), &(xd[0]));
        DERIV_val(min) = DERIV_val(xd[0]);
    }
    for (i = 1;     i < N;     )    {
        if (DERIV_val(xd[i]) < DERIV_val(min))         {
            {
                ad_grad_axpy_copy(&(min), &(xd[i]));
                DERIV_val(min) = DERIV_val(xd[i]);
            }
        }
        ad_var_3 = i++;
    }
    {
        ad_grad_axpy_copy(&(*ad_var_ret), &(min));
        DERIV_val(*ad_var_ret) = DERIV_val(min);
    }
    return;
}
void   ad_N_VWL2Norm_Serial(DERIV_TYPE  *ad_var_ret,ad_N_Vector  x,ad_N_Vector  w) {
long int  i, N;
DERIV_TYPE  sum, prodi, *xd, *wd;
void   *ad_var_0;
ad_N_VectorContent_Serial  ad_var_1, ad_var_2, ad_var_3;
long int  ad_var_4;
DERIV_TYPE  ad_var_5, ad_var_6;
    double  ad_loc_0;
    {
        ad_grad_axpy_0(&(sum));
        DERIV_val(sum) = 0.0;
    }
    ad_var_0 = (void  * )0;
    xd = wd = ad_var_0;
    ad_var_1 = (ad_N_VectorContent_Serial )x->content;
    N = ad_var_1->length;
    ad_var_2 = (ad_N_VectorContent_Serial )x->content;
    xd = ad_var_2->data;
    ad_var_3 = (ad_N_VectorContent_Serial )w->content;
    wd = ad_var_3->data;
    for (i = 0;     i < N;     )    {
        {
            ad_loc_0 = DERIV_val(xd[i]) * DERIV_val(wd[i]);
            ad_grad_axpy_2(&(prodi), DERIV_val(wd[i]), &(xd[i]), DERIV_val(xd[i]), &(wd[i]));
            DERIV_val(prodi) = ad_loc_0;
        }
        ad_RPower2( &ad_var_5, prodi);
        {
            ad_loc_0 = DERIV_val(sum) + DERIV_val(ad_var_5);
            ad_grad_axpy_2(&(sum), 1.000000000000000e+00, &(sum), 1.000000000000000e+00, &(ad_var_5));
            DERIV_val(sum) = ad_loc_0;
        }
        ad_var_4 = i++;
    }
    ad_RSqrt( &ad_var_6, sum);
    {
        ad_grad_axpy_copy(&(*ad_var_ret), &(ad_var_6));
        DERIV_val(*ad_var_ret) = DERIV_val(ad_var_6);
    }
    return;
}
void   ad_N_VL1Norm_Serial(DERIV_TYPE  *ad_var_ret,ad_N_Vector  x) {
long int  i, N;
DERIV_TYPE  sum, *xd;
void   *ad_var_0;
ad_N_VectorContent_Serial  ad_var_1, ad_var_2;
long int  ad_var_3;
DERIV_TYPE  ad_var_4;
    double  ad_loc_0;
    {
        ad_grad_axpy_0(&(sum));
        DERIV_val(sum) = 0.0;
    }
    ad_var_0 = (void  * )0;
    xd = ad_var_0;
    ad_var_1 = (ad_N_VectorContent_Serial )x->content;
    N = ad_var_1->length;
    ad_var_2 = (ad_N_VectorContent_Serial )x->content;
    xd = ad_var_2->data;
    for (i = 0;     i < N;     )    {
        ad_RAbs( &ad_var_4, xd[i]);
        {
            ad_loc_0 = DERIV_val(sum) + DERIV_val(ad_var_4);
            ad_grad_axpy_2(&(sum), 1.000000000000000e+00, &(sum), 1.000000000000000e+00, &(ad_var_4));
            DERIV_val(sum) = ad_loc_0;
        }
        ad_var_3 = i++;
    }
    {
        ad_grad_axpy_copy(&(*ad_var_ret), &(sum));
        DERIV_val(*ad_var_ret) = DERIV_val(sum);
    }
    return;
}
void   ad_N_VCompare_Serial(DERIV_TYPE  c,ad_N_Vector  x,ad_N_Vector  z) {
long int  i, N;
DERIV_TYPE  *xd, *zd;
void   *ad_var_0;
ad_N_VectorContent_Serial  ad_var_1, ad_var_2, ad_var_3;
long int  ad_var_4;
DERIV_TYPE  ad_var_5;
DERIV_TYPE  ad_var_6;
    ad_var_0 = (void  * )0;
    xd = zd = ad_var_0;
    ad_var_1 = (ad_N_VectorContent_Serial )x->content;
    N = ad_var_1->length;
    ad_var_2 = (ad_N_VectorContent_Serial )x->content;
    xd = ad_var_2->data;
    ad_var_3 = (ad_N_VectorContent_Serial )z->content;
    zd = ad_var_3->data;
    for (i = 0;     i < N;     )    {
        ad_RAbs( &ad_var_5, xd[i]);
        if (DERIV_val(ad_var_5) >= DERIV_val(c))         {
            {
                ad_grad_axpy_0(&(ad_var_6));
                DERIV_val(ad_var_6) = 1.0;
            }
        }
        else         {
            {
                ad_grad_axpy_0(&(ad_var_6));
                DERIV_val(ad_var_6) = 0.0;
            }
        }
        {
            ad_grad_axpy_copy(&(zd[i]), &(ad_var_6));
            DERIV_val(zd[i]) = DERIV_val(ad_var_6);
        }
        ad_var_4 = i++;
    }
    return;
}
int  ad_N_VInvTest_Serial(ad_N_Vector  x,ad_N_Vector  z) {
long int  i, N;
DERIV_TYPE  *xd, *zd;
void   *ad_var_0;
ad_N_VectorContent_Serial  ad_var_1, ad_var_2, ad_var_3;
long int  ad_var_4;
    double  ad_loc_0;
    double  ad_adj_0;
    ad_var_0 = (void  * )0;
    xd = zd = ad_var_0;
    ad_var_1 = (ad_N_VectorContent_Serial )x->content;
    N = ad_var_1->length;
    ad_var_2 = (ad_N_VectorContent_Serial )x->content;
    xd = ad_var_2->data;
    ad_var_3 = (ad_N_VectorContent_Serial )z->content;
    zd = ad_var_3->data;
    for (i = 0;     i < N;     )    {
        if (DERIV_val(xd[i]) == 0.0)         {
            return 0;
        }
        {
            ad_loc_0 = 1.0 / DERIV_val(xd[i]);
            ad_adj_0 =  -ad_loc_0 / DERIV_val(xd[i]);
            ad_grad_axpy_1(&(zd[i]), ad_adj_0, &(xd[i]));
            DERIV_val(zd[i]) = ad_loc_0;
        }
        ad_var_4 = i++;
    }
    return 1;
}
int  ad_N_VConstrMask_Serial(ad_N_Vector  c,ad_N_Vector  x,ad_N_Vector  m) {
long int  i, N;
int  test;
DERIV_TYPE  *cd, *xd, *md;
void   *ad_var_0;
ad_N_VectorContent_Serial  ad_var_1, ad_var_2, ad_var_3, ad_var_4;
long int  ad_var_5;
    ad_var_0 = (void  * )0;
    cd = xd = md = ad_var_0;
    ad_var_1 = (ad_N_VectorContent_Serial )x->content;
    N = ad_var_1->length;
    ad_var_2 = (ad_N_VectorContent_Serial )x->content;
    xd = ad_var_2->data;
    ad_var_3 = (ad_N_VectorContent_Serial )c->content;
    cd = ad_var_3->data;
    ad_var_4 = (ad_N_VectorContent_Serial )m->content;
    md = ad_var_4->data;
    test = 1;
    for (i = 0;     i < N;     )    {
        {
            ad_grad_axpy_0(&(md[i]));
            DERIV_val(md[i]) = 0.0;
        }
        if (DERIV_val(cd[i]) == 0.0)         {
            continue;
        }
        if (DERIV_val(cd[i]) > 1.5 || DERIV_val(cd[i]) <  -1.5)         {
            if (DERIV_val(xd[i]) * DERIV_val(cd[i]) <= 0.0)             {
                {
                    test = 0;
                    {
                        ad_grad_axpy_0(&(md[i]));
                        DERIV_val(md[i]) = 1.0;
                    }
                }
            }
            continue;
        }
        if (DERIV_val(cd[i]) > 0.5 || DERIV_val(cd[i]) <  -0.5)         {
            if (DERIV_val(xd[i]) * DERIV_val(cd[i]) < 0.0)             {
                {
                    test = 0;
                    {
                        ad_grad_axpy_0(&(md[i]));
                        DERIV_val(md[i]) = 1.0;
                    }
                }
            }
        }
        ad_var_5 = i++;
    }
    return test;
}
void   ad_N_VMinQuotient_Serial(DERIV_TYPE  *ad_var_ret,ad_N_Vector  num,ad_N_Vector  denom) {
int  notEvenOnce;
long int  i, N;
DERIV_TYPE  *nd, *dd, min;
void   *ad_var_0;
ad_N_VectorContent_Serial  ad_var_1, ad_var_2, ad_var_3;
long int  ad_var_4;
DERIV_TYPE  ad_var_5;
    double  ad_loc_0;
    double  ad_adj_0;
    double  ad_adj_1;
    ad_var_0 = (void  * )0;
    nd = dd = ad_var_0;
    ad_var_1 = (ad_N_VectorContent_Serial )num->content;
    N = ad_var_1->length;
    ad_var_2 = (ad_N_VectorContent_Serial )num->content;
    nd = ad_var_2->data;
    ad_var_3 = (ad_N_VectorContent_Serial )denom->content;
    dd = ad_var_3->data;
    notEvenOnce = 1;
    {
        ad_grad_axpy_0(&(min));
        DERIV_val(min) = 1.7976931348623157e+308;
    }
    for (i = 0;     i < N;     )    {
        if (DERIV_val(dd[i]) == 0.0)         {
            continue;
        }
        else         {
            if ( !notEvenOnce)             {
                if (DERIV_val(min) < DERIV_val(nd[i]) / DERIV_val(dd[i]))                 {
                    {
                        ad_grad_axpy_copy(&(ad_var_5), &(min));
                        DERIV_val(ad_var_5) = DERIV_val(min);
                    }
                }
                else                 {
                    {
                        ad_loc_0 = DERIV_val(nd[i]) / DERIV_val(dd[i]);
                        ad_adj_0 =  -ad_loc_0 / DERIV_val(dd[i]);
                        ad_adj_1 = 1.000000000000000e+00 / DERIV_val(dd[i]);
                        ad_grad_axpy_2(&(ad_var_5), ad_adj_1, &(nd[i]), ad_adj_0, &(dd[i]));
                        DERIV_val(ad_var_5) = ad_loc_0;
                    }
                }
                {
                    ad_grad_axpy_copy(&(min), &(ad_var_5));
                    DERIV_val(min) = DERIV_val(ad_var_5);
                }
            }
            else             {
                {
                    ad_loc_0 = DERIV_val(nd[i]) / DERIV_val(dd[i]);
                    ad_adj_0 =  -ad_loc_0 / DERIV_val(dd[i]);
                    ad_adj_1 = 1.000000000000000e+00 / DERIV_val(dd[i]);
                    ad_grad_axpy_2(&(min), ad_adj_1, &(nd[i]), ad_adj_0, &(dd[i]));
                    DERIV_val(min) = ad_loc_0;
                }
                notEvenOnce = 0;
            }
        }
        ad_var_4 = i++;
    }
    {
        ad_grad_axpy_copy(&(*ad_var_ret), &(min));
        DERIV_val(*ad_var_ret) = DERIV_val(min);
    }
    return;
}
static void   ad_VCopy_Serial(ad_N_Vector  x,ad_N_Vector  z) {
long int  i, N;
DERIV_TYPE  *xd, *zd;
void   *ad_var_0;
ad_N_VectorContent_Serial  ad_var_1, ad_var_2, ad_var_3;
long int  ad_var_4;
    ad_var_0 = (void  * )0;
    xd = zd = ad_var_0;
    ad_var_1 = (ad_N_VectorContent_Serial )x->content;
    N = ad_var_1->length;
    ad_var_2 = (ad_N_VectorContent_Serial )x->content;
    xd = ad_var_2->data;
    ad_var_3 = (ad_N_VectorContent_Serial )z->content;
    zd = ad_var_3->data;
    for (i = 0;     i < N;     )    {
        {
            ad_grad_axpy_copy(&(zd[i]), &(xd[i]));
            DERIV_val(zd[i]) = DERIV_val(xd[i]);
        }
        ad_var_4 = i++;
    }
    return;
}
static void   ad_VSum_Serial(ad_N_Vector  x,ad_N_Vector  y,ad_N_Vector  z) {
long int  i, N;
DERIV_TYPE  *xd, *yd, *zd;
void   *ad_var_0;
ad_N_VectorContent_Serial  ad_var_1, ad_var_2, ad_var_3, ad_var_4;
long int  ad_var_5;
    double  ad_loc_0;
    ad_var_0 = (void  * )0;
    xd = yd = zd = ad_var_0;
    ad_var_1 = (ad_N_VectorContent_Serial )x->content;
    N = ad_var_1->length;
    ad_var_2 = (ad_N_VectorContent_Serial )x->content;
    xd = ad_var_2->data;
    ad_var_3 = (ad_N_VectorContent_Serial )y->content;
    yd = ad_var_3->data;
    ad_var_4 = (ad_N_VectorContent_Serial )z->content;
    zd = ad_var_4->data;
    for (i = 0;     i < N;     )    {
        {
            ad_loc_0 = DERIV_val(xd[i]) + DERIV_val(yd[i]);
            ad_grad_axpy_2(&(zd[i]), 1.000000000000000e+00, &(xd[i]), 1.000000000000000e+00, &(yd[i]));
            DERIV_val(zd[i]) = ad_loc_0;
        }
        ad_var_5 = i++;
    }
    return;
}
static void   ad_VDiff_Serial(ad_N_Vector  x,ad_N_Vector  y,ad_N_Vector  z) {
long int  i, N;
DERIV_TYPE  *xd, *yd, *zd;
void   *ad_var_0;
ad_N_VectorContent_Serial  ad_var_1, ad_var_2, ad_var_3, ad_var_4;
long int  ad_var_5;
    double  ad_loc_0;
    ad_var_0 = (void  * )0;
    xd = yd = zd = ad_var_0;
    ad_var_1 = (ad_N_VectorContent_Serial )x->content;
    N = ad_var_1->length;
    ad_var_2 = (ad_N_VectorContent_Serial )x->content;
    xd = ad_var_2->data;
    ad_var_3 = (ad_N_VectorContent_Serial )y->content;
    yd = ad_var_3->data;
    ad_var_4 = (ad_N_VectorContent_Serial )z->content;
    zd = ad_var_4->data;
    for (i = 0;     i < N;     )    {
        {
            ad_loc_0 = DERIV_val(xd[i]) - DERIV_val(yd[i]);
            ad_grad_axpy_2(&(zd[i]), 1.000000000000000e+00, &(xd[i]), -1.000000000000000e+00, &(yd[i]));
            DERIV_val(zd[i]) = ad_loc_0;
        }
        ad_var_5 = i++;
    }
    return;
}
static void   ad_VNeg_Serial(ad_N_Vector  x,ad_N_Vector  z) {
long int  i, N;
DERIV_TYPE  *xd, *zd;
void   *ad_var_0;
ad_N_VectorContent_Serial  ad_var_1, ad_var_2, ad_var_3;
long int  ad_var_4;
    double  ad_loc_0;
    ad_var_0 = (void  * )0;
    xd = zd = ad_var_0;
    ad_var_1 = (ad_N_VectorContent_Serial )x->content;
    N = ad_var_1->length;
    ad_var_2 = (ad_N_VectorContent_Serial )x->content;
    xd = ad_var_2->data;
    ad_var_3 = (ad_N_VectorContent_Serial )z->content;
    zd = ad_var_3->data;
    for (i = 0;     i < N;     )    {
        {
            ad_loc_0 =  -DERIV_val(xd[i]);
            ad_grad_axpy_1(&(zd[i]), -1.000000000000000e+00, &(xd[i]));
            DERIV_val(zd[i]) = ad_loc_0;
        }
        ad_var_4 = i++;
    }
    return;
}
static void   ad_VScaleSum_Serial(DERIV_TYPE  c,ad_N_Vector  x,ad_N_Vector  y,ad_N_Vector  z) {
long int  i, N;
DERIV_TYPE  *xd, *yd, *zd;
void   *ad_var_0;
ad_N_VectorContent_Serial  ad_var_1, ad_var_2, ad_var_3, ad_var_4;
long int  ad_var_5;
    double  ad_loc_0;
    double  ad_loc_1;
    ad_var_0 = (void  * )0;
    xd = yd = zd = ad_var_0;
    ad_var_1 = (ad_N_VectorContent_Serial )x->content;
    N = ad_var_1->length;
    ad_var_2 = (ad_N_VectorContent_Serial )x->content;
    xd = ad_var_2->data;
    ad_var_3 = (ad_N_VectorContent_Serial )y->content;
    yd = ad_var_3->data;
    ad_var_4 = (ad_N_VectorContent_Serial )z->content;
    zd = ad_var_4->data;
    for (i = 0;     i < N;     )    {
        {
            ad_loc_0 = DERIV_val(xd[i]) + DERIV_val(yd[i]);
            ad_loc_1 = DERIV_val(c) * ad_loc_0;
            ad_grad_axpy_3(&(zd[i]), ad_loc_0, &(c), DERIV_val(c), &(xd[i]), DERIV_val(c), &(yd[i]));
            DERIV_val(zd[i]) = ad_loc_1;
        }
        ad_var_5 = i++;
    }
    return;
}
static void   ad_VScaleDiff_Serial(DERIV_TYPE  c,ad_N_Vector  x,ad_N_Vector  y,ad_N_Vector  z) {
long int  i, N;
DERIV_TYPE  *xd, *yd, *zd;
void   *ad_var_0;
ad_N_VectorContent_Serial  ad_var_1, ad_var_2, ad_var_3, ad_var_4;
long int  ad_var_5;
    double  ad_loc_0;
    double  ad_loc_1;
    double  ad_adj_0;
    ad_var_0 = (void  * )0;
    xd = yd = zd = ad_var_0;
    ad_var_1 = (ad_N_VectorContent_Serial )x->content;
    N = ad_var_1->length;
    ad_var_2 = (ad_N_VectorContent_Serial )x->content;
    xd = ad_var_2->data;
    ad_var_3 = (ad_N_VectorContent_Serial )y->content;
    yd = ad_var_3->data;
    ad_var_4 = (ad_N_VectorContent_Serial )z->content;
    zd = ad_var_4->data;
    for (i = 0;     i < N;     )    {
        {
            ad_loc_0 = DERIV_val(xd[i]) - DERIV_val(yd[i]);
            ad_loc_1 = DERIV_val(c) * ad_loc_0;
            ad_adj_0 =  -DERIV_val(c);
            ad_grad_axpy_3(&(zd[i]), ad_loc_0, &(c), DERIV_val(c), &(xd[i]), ad_adj_0, &(yd[i]));
            DERIV_val(zd[i]) = ad_loc_1;
        }
        ad_var_5 = i++;
    }
    return;
}
static void   ad_VLin1_Serial(DERIV_TYPE  a,ad_N_Vector  x,ad_N_Vector  y,ad_N_Vector  z) {
long int  i, N;
DERIV_TYPE  *xd, *yd, *zd;
void   *ad_var_0;
ad_N_VectorContent_Serial  ad_var_1, ad_var_2, ad_var_3, ad_var_4;
long int  ad_var_5;
    double  ad_loc_0;
    double  ad_loc_1;
    ad_var_0 = (void  * )0;
    xd = yd = zd = ad_var_0;
    ad_var_1 = (ad_N_VectorContent_Serial )x->content;
    N = ad_var_1->length;
    ad_var_2 = (ad_N_VectorContent_Serial )x->content;
    xd = ad_var_2->data;
    ad_var_3 = (ad_N_VectorContent_Serial )y->content;
    yd = ad_var_3->data;
    ad_var_4 = (ad_N_VectorContent_Serial )z->content;
    zd = ad_var_4->data;
    for (i = 0;     i < N;     )    {
        {
            ad_loc_0 = DERIV_val(a) * DERIV_val(xd[i]);
            ad_loc_1 = ad_loc_0 + DERIV_val(yd[i]);
            ad_grad_axpy_3(&(zd[i]), DERIV_val(xd[i]), &(a), DERIV_val(a), &(xd[i]), 1.000000000000000e+00, &(yd[i]));
            DERIV_val(zd[i]) = ad_loc_1;
        }
        ad_var_5 = i++;
    }
    return;
}
static void   ad_VLin2_Serial(DERIV_TYPE  a,ad_N_Vector  x,ad_N_Vector  y,ad_N_Vector  z) {
long int  i, N;
DERIV_TYPE  *xd, *yd, *zd;
void   *ad_var_0;
ad_N_VectorContent_Serial  ad_var_1, ad_var_2, ad_var_3, ad_var_4;
long int  ad_var_5;
    double  ad_loc_0;
    double  ad_loc_1;
    ad_var_0 = (void  * )0;
    xd = yd = zd = ad_var_0;
    ad_var_1 = (ad_N_VectorContent_Serial )x->content;
    N = ad_var_1->length;
    ad_var_2 = (ad_N_VectorContent_Serial )x->content;
    xd = ad_var_2->data;
    ad_var_3 = (ad_N_VectorContent_Serial )y->content;
    yd = ad_var_3->data;
    ad_var_4 = (ad_N_VectorContent_Serial )z->content;
    zd = ad_var_4->data;
    for (i = 0;     i < N;     )    {
        {
            ad_loc_0 = DERIV_val(a) * DERIV_val(xd[i]);
            ad_loc_1 = ad_loc_0 - DERIV_val(yd[i]);
            ad_grad_axpy_3(&(zd[i]), DERIV_val(xd[i]), &(a), DERIV_val(a), &(xd[i]), -1.000000000000000e+00, &(yd[i]));
            DERIV_val(zd[i]) = ad_loc_1;
        }
        ad_var_5 = i++;
    }
    return;
}
static void   ad_Vaxpy_Serial(DERIV_TYPE  a,ad_N_Vector  x,ad_N_Vector  y) {
long int  i, N;
DERIV_TYPE  *xd, *yd;
void   *ad_var_0;
ad_N_VectorContent_Serial  ad_var_1, ad_var_2, ad_var_3;
long int  ad_var_4, ad_var_5, ad_var_6;
    double  ad_loc_0;
    double  ad_loc_1;
    ad_var_0 = (void  * )0;
    xd = yd = ad_var_0;
    ad_var_1 = (ad_N_VectorContent_Serial )x->content;
    N = ad_var_1->length;
    ad_var_2 = (ad_N_VectorContent_Serial )x->content;
    xd = ad_var_2->data;
    ad_var_3 = (ad_N_VectorContent_Serial )y->content;
    yd = ad_var_3->data;
    if (DERIV_val(a) == 1.0)     {
        for (i = 0;         i < N;         )        {
            {
                ad_loc_0 = DERIV_val(yd[i]) + DERIV_val(xd[i]);
                ad_grad_axpy_2(&(yd[i]), 1.000000000000000e+00, &(yd[i]), 1.000000000000000e+00, &(xd[i]));
                DERIV_val(yd[i]) = ad_loc_0;
            }
            ad_var_4 = i++;
        }
        return;
    }
    if (DERIV_val(a) ==  -1.0)     {
        for (i = 0;         i < N;         )        {
            {
                ad_loc_0 = DERIV_val(yd[i]) - DERIV_val(xd[i]);
                ad_grad_axpy_2(&(yd[i]), 1.000000000000000e+00, &(yd[i]), -1.000000000000000e+00, &(xd[i]));
                DERIV_val(yd[i]) = ad_loc_0;
            }
            ad_var_5 = i++;
        }
        return;
    }
    for (i = 0;     i < N;     )    {
        {
            ad_loc_0 = DERIV_val(a) * DERIV_val(xd[i]);
            ad_loc_1 = DERIV_val(yd[i]) + ad_loc_0;
            ad_grad_axpy_3(&(yd[i]), 1.000000000000000e+00, &(yd[i]), DERIV_val(xd[i]), &(a), DERIV_val(a), &(xd[i]));
            DERIV_val(yd[i]) = ad_loc_1;
        }
        ad_var_6 = i++;
    }
    return;
}
static void   ad_VScaleBy_Serial(DERIV_TYPE  a,ad_N_Vector  x) {
long int  i, N;
DERIV_TYPE  *xd;
void   *ad_var_0;
ad_N_VectorContent_Serial  ad_var_1, ad_var_2;
long int  ad_var_3;
    double  ad_loc_0;
    ad_var_0 = (void  * )0;
    xd = ad_var_0;
    ad_var_1 = (ad_N_VectorContent_Serial )x->content;
    N = ad_var_1->length;
    ad_var_2 = (ad_N_VectorContent_Serial )x->content;
    xd = ad_var_2->data;
    for (i = 0;     i < N;     )    {
        {
            ad_loc_0 = DERIV_val(xd[i]) * DERIV_val(a);
            ad_grad_axpy_2(&(xd[i]), DERIV_val(a), &(xd[i]), DERIV_val(xd[i]), &(a));
            DERIV_val(xd[i]) = ad_loc_0;
        }
        ad_var_3 = i++;
    }
    return;
}
void   ad_AD_Init(int  arg0) {
    ad_AD_GradInit(arg0);

}
void   ad_AD_Final() {
    ad_AD_GradFinal();

}
