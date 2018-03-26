#include <stdlib.h>
#include <nvector_serial.h>

#define Ith(v,i) NV_Ith_S(v,i-1)

typedef struct {
  realtype p[3];
} *UserData;

void f(realtype t, N_Vector y, N_Vector ydot, void *f_data);
