#ifndef included_PhaseFluxStrategyIsotropic
#define included_PhaseFluxStrategyIsotropic

#include "PhaseFluxStrategy.h"

class PhaseFluxStrategyIsotropic:
   public PhaseFluxStrategy
{
public:
   PhaseFluxStrategyIsotropic(const double epsilon_phase):
      d_epsilon_phase(epsilon_phase)
   {
      tbox::plog<<"Uses PhaseFluxStrategyIsotropic class..."<<std::endl;
   }
   
   void computeFluxes(const boost::shared_ptr<hier::PatchLevel> level,
                      const int phase_id,
                      const int quat_id,
                      const int flux_id);

private:
   const double d_epsilon_phase;
};

#endif
