#ifndef TbasedCompositionDiffusionStrategy_H
#define TbasedCompositionDiffusionStrategy_H

#include "CompositionDiffusionStrategy.h"

#include <cstring>

using namespace std;

class TbasedCompositionDiffusionStrategy :
   public CompositionDiffusionStrategy
{
public:

   TbasedCompositionDiffusionStrategy(
      const int diffusion_l_id,
      const int diffusion_a_id,
      const double D_liquid, const double Q0_liquid,
      const double D_solid_A, const double Q0_solid_A,
      const string& phase_interp_func_type,
      const string& avg_func_type
);


/*
 * compute actual diffusion by weighting diffusion in each pahse
 * using phase variable
 */
   virtual void setDiffCoeff(
      const boost::shared_ptr< hier::PatchHierarchy > hierarchy,
      const int temperature_id,
      const int phase_id,
      const int eta_id);

private:
   /*!
    * holds data for diffusion coefficients in composition equation
    * weighted by phase fraction
    */
   int d_diffusion_l_id;
   int d_diffusion_a_id;

   double d_D_liquid;
   double d_D_solid_A;
   double d_Q0_liquid;
   double d_Q0_solid_A;

   std::string d_phase_interp_func_type;
   std::string d_avg_func_type;

   bool d_with_third_phase;
};

#endif
