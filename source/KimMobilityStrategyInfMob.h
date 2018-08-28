#ifndef included_KimMobilityStrategyInfMob
#define included_KimMobilityStrategyInfMob

#include "KimMobilityStrategy.h"

class KimMobilityStrategyInfMob:
   public KimMobilityStrategy
{
public:

   KimMobilityStrategyInfMob(
      QuatModel* quat_model,
      const int conc_l_id,
      const int conc_s_id,
      const int temp_id,
      const double epsilon,
      const double phase_well_scale,
      const std::string& energy_interp_func_type,
      const std::string& conc_interp_func_type,
      boost::shared_ptr<tbox::Database> calphad_db,
      boost::shared_ptr<tbox::Database> newton_db,
      const unsigned ncompositions,
      const double DL,
      const double Q0,
      const double mv);

private:

   double evaluateMobility(const double temp,
      const std::vector<double>&  phaseconc);

   std::vector<double> d_d2fdc2;

   double d_factor; 
};

#endif


