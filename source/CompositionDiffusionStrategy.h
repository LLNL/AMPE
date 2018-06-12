#ifndef CompositionDiffusionStrategy_H
#define CompositionDiffusionStrategy_H

#include "SAMRAI/hier/PatchHierarchy.h"

#include <boost/make_shared.hpp>

using namespace SAMRAI;

enum class DiffusionInterpolationType
{
   LINEAR,
   PBG,
   UNDEFINED
};

class CompositionDiffusionStrategy
{
public:
   CompositionDiffusionStrategy(DiffusionInterpolationType interp_func_type)
      : d_interp_func_type(interp_func_type)
   {};

/*
 * compute actual diffusion by weighting diffusion in each phase
 * using phase variable
 */
   virtual void setDiffusion(
      const boost::shared_ptr< hier::PatchHierarchy > hierarchy,
      const int temperature_id,
      const int phase_id,
      const int eta_id)=0;

protected:
   char interpChar()const
   {
      return ((d_interp_func_type
       ==DiffusionInterpolationType::LINEAR) ? 'l' : 'p');
   }

private:
   DiffusionInterpolationType d_interp_func_type;
};

#endif
