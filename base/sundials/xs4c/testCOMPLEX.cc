#include <stdio.h>
#include <math.h>
#include "complexify.h"

int main() {
     cplx a, b, c, x1, x2;
     a = 1;
     b = -4;
     c=13;

     x1 = (-b + sqrt(b*b-4*a*c))/(2*a);
     x2 = (-b - sqrt(b*b-4*a*c))/(2*a);

     printf("x1 = complex(%f,%f)\n", real(x1), imag(x1));
     printf("x2 = complex(%f,%f)\n", real(x2), imag(x2));
  }

