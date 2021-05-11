#ifndef included_ConcFreeEnergyStrategy
#define included_ConcFreeEnergyStrategy

#include "FreeEnergyStrategy.h"

class ConcFreeEnergyStrategy : public FreeEnergyStrategy
{
 public:
   virtual bool computeCeqT(const double temperature, const PhaseIndex pi0,
                            const PhaseIndex pi1, double* ceq) = 0;

   virtual void energyVsPhiAndC(const double temperature,
                                const double* const ceq, const bool found_ceq,
                                const double phi_well_scale,
                                const std::string& phi_well_type,
                                const int npts_phi, const int npts_c) = 0;
};

#endif
