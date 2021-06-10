#ifndef included_PhaseFluxStrategyFactory
#define included_PhaseFluxStrategyFactory

#include "PhaseFluxStrategyAnisotropy.h"
#include "PhaseFluxStrategyIsotropic.h"
#include "PhaseFluxStrategySimple.h"

class PhaseFluxStrategyFactory
{
 public:
   static std::shared_ptr<PhaseFluxStrategy> create(QuatModelParameters& model_parameters)
   {
   std::shared_ptr<PhaseFluxStrategy> phase_flux_strategy;

   double epsilon_anisotropy = model_parameters.epsilon_anisotropy();

   if (epsilon_anisotropy >= 0.)
      phase_flux_strategy.reset(
          new PhaseFluxStrategyAnisotropy(model_parameters.epsilon_phase(),
                                          epsilon_anisotropy, 4));
   else if (model_parameters.useIsotropicStencil()) {
      phase_flux_strategy.reset(
          new PhaseFluxStrategyIsotropic(model_parameters.epsilon_phase()));
   } else {
      phase_flux_strategy.reset(
          new PhaseFluxStrategySimple(model_parameters.epsilon_phase()));
   }

   return phase_flux_strategy;
   }
};
#endif
