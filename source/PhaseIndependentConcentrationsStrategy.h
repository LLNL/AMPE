#ifndef included_PhaseIndependentConcentrationsStrategy
#define included_PhaseIndependentConcentrationsStrategy

#include "PhaseConcentrationsStrategy.h"

/*
 * Simply set c_l and c_s equal to c
 */
class PhaseIndependentConcentrationsStrategy:public PhaseConcentrationsStrategy
{
public:

   PhaseIndependentConcentrationsStrategy(
      const int conc_l_id,
      const int conc_a_id,
      const int conc_b_id):
         PhaseConcentrationsStrategy(
            conc_l_id,
            conc_a_id,
            conc_b_id,
            false)
   {
   };
   
   ~PhaseIndependentConcentrationsStrategy(){};

   virtual void computePhaseConcentrationsOnPatch(
      boost::shared_ptr< pdat::CellData<double> > cd_temperature,
      boost::shared_ptr< pdat::CellData<double> > cd_phi,
      boost::shared_ptr< pdat::CellData<double> > cd_eta,
      boost::shared_ptr< pdat::CellData<double> > cd_concentration,
      boost::shared_ptr< pdat::CellData<double> > cd_c_l,
      boost::shared_ptr< pdat::CellData<double> > cd_c_a,
      boost::shared_ptr< pdat::CellData<double> > cd_c_b,
      boost::shared_ptr<hier::Patch > patch );

private:

};

#endif
