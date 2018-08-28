#include "KimMobilityStrategyInfMob.h"

KimMobilityStrategyInfMob::KimMobilityStrategyInfMob(
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
   const double mv)
   :  KimMobilityStrategy(quat_model,conc_l_id,conc_s_id,temp_id,
         energy_interp_func_type,conc_interp_func_type,
         calphad_db, newton_db,
         ncompositions,
         DL, Q0)
{
   assert( epsilon>0. );
   assert( phase_well_scale>=0. );
   assert( mv>0. );

   const double a2 = 47./60.;
   const double xi = epsilon/sqrt(32.*phase_well_scale);

   d_factor = 3.*(2.*xi*xi)*a2;
   d_factor *= (1.e-6/mv); // convert from J/mol to pJ/um^3

   d_d2fdc2.resize(d_ncompositions*d_ncompositions);
}

double KimMobilityStrategyInfMob::evaluateMobility(
   const double temp,
   const std::vector<double>&  phaseconc)
{
   const PHASE_INDEX pi0=phaseL;

   d_calphad_fenergy->computeSecondDerivativeFreeEnergy(
      temp,&phaseconc[0],pi0,d_d2fdc2);

   const double* const cl=&phaseconc[0];
   const double* const cs=&phaseconc[d_ncompositions];

   double zeta=0.;
   for(unsigned i=0;i<d_ncompositions;i++)
   for(unsigned j=0;j<d_ncompositions;j++)
      zeta+=(cl[i]-cs[i])*d_d2fdc2[2*i+j]*(cl[j]-cs[j]);
   const double DL=d_DL*exp(-d_Q0/(gas_constant_R_JpKpmol*temp));

   return DL/(d_factor*zeta);
}
