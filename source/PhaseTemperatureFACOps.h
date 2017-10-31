#ifndef included_PhaseTemperatureFACOps
#define included_PhaseTemperatureFACOps

#include "EllipticFACOps.h"

#include <string>

using namespace SAMRAI;

class PhaseTemperatureFACOps
   : public EllipticFACOps
{

public:

   PhaseTemperatureFACOps(
      const std::string &object_name ,
      const boost::shared_ptr<tbox::Database> database=
         boost::shared_ptr<tbox::Database>() );

   void setOperatorCoefficients(
      const int phase_id,
      const int mobility_id,
      const double epsilon_phase, 
      const double latent_heat,
      const std::string phase_interp_func_type,
      const double phase_well_scale,
      const std::string phase_well_func_type);

   void multiplyDTDPhiBlock(const int,const int);

private:

   void setC(
      const int phi_id,
      const double latent_heat,
      const std::string phi_interp_func_type,
      const double phi_well_scale,
      const std::string phi_well_func_type);

   void setCOnPatchForPreconditionODE(
      boost::shared_ptr< pdat::CellData<double> > cd_phi,
      boost::shared_ptr< pdat::CellData<double> > cd_m,
      boost::shared_ptr< pdat::CellData<double> > cd_c,
      const double latent_heat,
      const char* phi_interp_func_type,
      const double phi_well_scale,
      const char* phi_well_func_type,
      const hier::Box& pbox );

};

#endif // included_PhaseTemperatureFACOps
