#ifndef included_PhaseFluxStrategyFactory
#define included_PhaseFluxStrategyFactory

#include "PhaseFluxStrategyAnisotropy.h"
#include "PhaseFluxStrategyIsotropic.h"
#include "PhaseFluxStrategySimple.h"
#include "PhaseFluxStrategyAnisotropyMultiOrder.h"

class PhaseFluxStrategyFactory
{
 public:
   static std::shared_ptr<PhaseFluxStrategy> create(
       QuatModelParameters& model_parameters)
   {
      std::shared_ptr<PhaseFluxStrategy> phase_flux_strategy;

      double epsilon_anisotropy = model_parameters.epsilon_anisotropy();
      double epsilon = model_parameters.epsilon_phase();

      bool use_epsilon = false;
      if (model_parameters.with_concentration())
         if (model_parameters.isConcentrationModelWangSintering()) {
            // kappa = epsilon^2/2
            // so just use epsilon*epsilon = 6*sigma*delta = 2*kappa
            // in flux computation
            use_epsilon = true;
         }
      if (!use_epsilon) {
         // if multi order parameters, divide d_epsilon_phase by sqrt(2.)
         // to compensate for double counting
         if (model_parameters.norderp() > 1 &&
             !model_parameters.with_three_phases())
            epsilon /= sqrt(2.);
      }

      if (epsilon_anisotropy >= 0.) {
         if (!model_parameters.with_orientation())
            TBOX_ERROR("Phase anisotropy requires quaternion orientation");
         if (model_parameters.norientations() > 0) {
            phase_flux_strategy.reset(new PhaseFluxStrategyAnisotropyMultiOrder(
                epsilon, epsilon_anisotropy, 4,
                model_parameters.orientations()));
         } else {
            phase_flux_strategy.reset(
                new PhaseFluxStrategyAnisotropy(epsilon, epsilon_anisotropy,
                                                4));
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
