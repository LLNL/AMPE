#ifndef included_PhaseFluxStrategyFactory
#define included_PhaseFluxStrategyFactory

#include "PhaseFluxStrategyAnisotropy.h"
#include "PhaseFluxStrategyIsotropic.h"
#include "PhaseFluxStrategySimple.h"
#include "PhaseFluxStrategyQuaternion.h"

class PhaseFluxStrategyFactory
{
 public:
   static std::shared_ptr<PhaseFluxStrategy> create(
       QuatModelParameters& model_parameters)
   {
      std::shared_ptr<PhaseFluxStrategy> phase_flux_strategy;

      double epsilon_anisotropy = model_parameters.epsilon_anisotropy();
      double epsilon = model_parameters.epsilon_phase();
      // if multi order parameters, divide d_epsilon_phase by sqrt(2.)
      // to compensate for double counting
      if (model_parameters.norderp() > 1 &&
          !model_parameters.with_three_phases())
         epsilon /= sqrt(2.);

      if (epsilon_anisotropy >= 0.) {
         if (model_parameters.with_orientation()) {
            phase_flux_strategy.reset(
                new PhaseFluxStrategyAnisotropy(epsilon, epsilon_anisotropy,
                                                4));
         } else {
            std::vector<std::array<double, 4>> quat =
                model_parameters.orderp_quat();
            phase_flux_strategy.reset(
                new PhaseFluxStrategyQuaternion(epsilon, epsilon_anisotropy, 4,
                                                quat));
         }
      } else if (model_parameters.useIsotropicStencil()) {
         phase_flux_strategy.reset(new PhaseFluxStrategyIsotropic(epsilon));
      } else {
         phase_flux_strategy.reset(new PhaseFluxStrategySimple(epsilon));
      }

      return phase_flux_strategy;
   }
};
#endif
