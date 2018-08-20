#include "CALPHADSpeciesPhaseGibbsEnergyExpansion.h"

#include <math.h>

CALPHADSpeciesPhaseGibbsEnergyExpansion::
   CALPHADSpeciesPhaseGibbsEnergyExpansion(
      const double a,
      const double b,
      const double c,
      const double d2,
      const double d3,
      const double d4,
      const double d7,
      const double dm1,
      const double dm9):
         d_a(a),
         d_b(b),
         d_c(c),
         d_d2(d2),
         d_d3(d3),
         d_d4(d4),
         d_d7(d7),
         d_dm1(dm1),
         d_dm9(dm9)
{
}

double CALPHADSpeciesPhaseGibbsEnergyExpansion::value(
   const double temperature)const
{
   const double t2=temperature*temperature;
   const double t4=t2*t2;

   return d_a
         +d_b*temperature
         +d_c*temperature*log(temperature)
         +d_d2*t2
         +d_d3*t2*temperature
         +d_d4*t4
         +d_d7*t4*t2*temperature
         +d_dm1/temperature
         +d_dm9/(t4*t4*temperature);
}

