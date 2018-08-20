#ifndef included_CALPHADSpeciesPhaseGibbsEnergyExpansion
#define included_CALPHADSpeciesPhaseGibbsEnergyExpansion

#include <string>

class CALPHADSpeciesPhaseGibbsEnergyExpansion
{
public:
   CALPHADSpeciesPhaseGibbsEnergyExpansion(
      const double a,
      const double b,
      const double c,
      const double d2,
      const double d3,
      const double d4,
      const double d7,
      const double dm1,
      const double dm9);

   double value(const double temperature)const;
      
private:
   /*
    * Expansion coefficient for energy of species as a function of temperature
    */
   const double d_a;
   const double d_b;
   const double d_c;
   const double d_d2;
   const double d_d3;
   const double d_d4;
   const double d_d7;
   const double d_dm1;
   const double d_dm9;
};

#endif

