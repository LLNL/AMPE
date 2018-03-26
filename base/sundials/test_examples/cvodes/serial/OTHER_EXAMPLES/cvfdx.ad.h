/************************** DISCLAIMER ********************************/
/*                                                                    */
/*   This file was generated on 06/22/05 14:41:00 by the version of   */
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
#include <stdlib.h>
#include <float.h>
typedef DERIV_TYPE  ad_realtype;
typedef struct ad__generic_N_Vector_Ops  *ad_N_Vector_Ops;
typedef struct ad__generic_N_Vector  *ad_N_Vector;
typedef ad_N_Vector  *ad_N_Vector_S;
struct ad__generic_N_Vector_Ops {
ad_N_Vector  (*nvclone)(ad_N_Vector  );
ad_N_Vector  (*nvcloneempty)(ad_N_Vector  );
void   (*nvdestroy)(ad_N_Vector  );
void   (*nvspace)(ad_N_Vector  ,long int  *,long int  *);
DERIV_TYPE  *(*nvgetarraypointer)(ad_N_Vector  );
void   (*nvsetarraypointer)(DERIV_TYPE  *,ad_N_Vector  );
void   (*nvlinearsum)(DERIV_TYPE  ,ad_N_Vector  ,DERIV_TYPE  ,ad_N_Vector  ,ad_N_Vector  );
void   (*nvconst)(DERIV_TYPE  ,ad_N_Vector  );
void   (*nvprod)(ad_N_Vector  ,ad_N_Vector  ,ad_N_Vector  );
void   (*nvdiv)(ad_N_Vector  ,ad_N_Vector  ,ad_N_Vector  );
void   (*nvscale)(DERIV_TYPE  ,ad_N_Vector  ,ad_N_Vector  );
void   (*nvabs)(ad_N_Vector  ,ad_N_Vector  );
void   (*nvinv)(ad_N_Vector  ,ad_N_Vector  );
void   (*nvaddconst)(ad_N_Vector  ,DERIV_TYPE  ,ad_N_Vector  );
void   (*nvdotprod)(DERIV_TYPE  *ad_var_ret,ad_N_Vector  ,ad_N_Vector  );
void   (*nvmaxnorm)(DERIV_TYPE  *ad_var_ret,ad_N_Vector  );
void   (*nvwrmsnorm)(DERIV_TYPE  *ad_var_ret,ad_N_Vector  ,ad_N_Vector  );
void   (*nvwrmsnormmask)(DERIV_TYPE  *ad_var_ret,ad_N_Vector  ,ad_N_Vector  ,ad_N_Vector  );
void   (*nvmin)(DERIV_TYPE  *ad_var_ret,ad_N_Vector  );
void   (*nvwl2norm)(DERIV_TYPE  *ad_var_ret,ad_N_Vector  ,ad_N_Vector  );
void   (*nvl1norm)(DERIV_TYPE  *ad_var_ret,ad_N_Vector  );
void   (*nvcompare)(DERIV_TYPE  ,ad_N_Vector  ,ad_N_Vector  );
int  (*nvinvtest)(ad_N_Vector  ,ad_N_Vector  );
int  (*nvconstrmask)(ad_N_Vector  ,ad_N_Vector  ,ad_N_Vector  );
void   (*nvminquotient)(DERIV_TYPE  *ad_var_ret,ad_N_Vector  ,ad_N_Vector  );
}
;
struct ad__generic_N_Vector {
void   *content;
    struct ad__generic_N_Vector_Ops  *ops;
}
;
ad_N_Vector  ad_N_VClone(ad_N_Vector  w);
ad_N_Vector  ad_N_VCloneEmpty(ad_N_Vector  w);
void   ad_N_VDestroy(ad_N_Vector  v);
void   ad_N_VSpace(ad_N_Vector  v,long int  *lrw,long int  *liw);
DERIV_TYPE  *ad_N_VGetArrayPointer(ad_N_Vector  v);
void   ad_N_VSetArrayPointer(DERIV_TYPE  *v_data,ad_N_Vector  v);
void   ad_N_VLinearSum(DERIV_TYPE  a,ad_N_Vector  x,DERIV_TYPE  b,ad_N_Vector  y,ad_N_Vector  z);
void   ad_N_VConst(DERIV_TYPE  c,ad_N_Vector  z);
void   ad_N_VProd(ad_N_Vector  x,ad_N_Vector  y,ad_N_Vector  z);
void   ad_N_VDiv(ad_N_Vector  x,ad_N_Vector  y,ad_N_Vector  z);
void   ad_N_VScale(DERIV_TYPE  c,ad_N_Vector  x,ad_N_Vector  z);
void   ad_N_VAbs(ad_N_Vector  x,ad_N_Vector  z);
void   ad_N_VInv(ad_N_Vector  x,ad_N_Vector  z);
void   ad_N_VAddConst(ad_N_Vector  x,DERIV_TYPE  b,ad_N_Vector  z);
void   ad_N_VDotProd(DERIV_TYPE  *ad_var_ret,ad_N_Vector  x,ad_N_Vector  y);
void   ad_N_VMaxNorm(DERIV_TYPE  *ad_var_ret,ad_N_Vector  x);
void   ad_N_VWrmsNorm(DERIV_TYPE  *ad_var_ret,ad_N_Vector  x,ad_N_Vector  w);
void   ad_N_VWrmsNormMask(DERIV_TYPE  *ad_var_ret,ad_N_Vector  x,ad_N_Vector  w,ad_N_Vector  id);
void   ad_N_VMin(DERIV_TYPE  *ad_var_ret,ad_N_Vector  x);
void   ad_N_VWL2Norm(DERIV_TYPE  *ad_var_ret,ad_N_Vector  x,ad_N_Vector  w);
void   ad_N_VL1Norm(DERIV_TYPE  *ad_var_ret,ad_N_Vector  x);
void   ad_N_VCompare(DERIV_TYPE  c,ad_N_Vector  x,ad_N_Vector  z);
int  ad_N_VInvTest(ad_N_Vector  x,ad_N_Vector  z);
int  ad_N_VConstrMask(ad_N_Vector  c,ad_N_Vector  x,ad_N_Vector  m);
void   ad_N_VMinQuotient(DERIV_TYPE  *ad_var_ret,ad_N_Vector  num,ad_N_Vector  denom);
ad_N_Vector  *ad_N_VCloneEmptyVectorArray(int  count,ad_N_Vector  w);
ad_N_Vector  *ad_N_VCloneVectorArray(int  count,ad_N_Vector  w);
void   ad_N_VDestroyVectorArray(ad_N_Vector  *vs,int  count);
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
typedef struct  {
DERIV_TYPE  p[3];
}
*ad_UserData;
void   ad_f(DERIV_TYPE  t,ad_N_Vector  y,ad_N_Vector  ydot,void   *f_data);
