// Ref: Eiken, Boettger, Steinbach, PRE 73, 066122 (2006)
#ifndef included_EBSCompositionRHSStrategy
#define included_EBSCompositionRHSStrategy 

#include "CompositionRHSStrategy.h"
#include "CALPHADMobility.h"
#include "CALPHADFreeEnergyStrategyBinary.h"
#include "CompositionDiffusionStrategy.h"

#include <vector>
#include <string>

class CompositionDiffusionStrategy;
class CompositionStrategyMobilities;

class EBSCompositionRHSStrategy:
   public CompositionRHSStrategy
{
public:
   EBSCompositionRHSStrategy(
      const int phase_scratch_id,
      const int eta_scratch_id,
      const unsigned short ncompositions,
      const int conc_l_scratch_id,
      const int conc_a_scratch_id,
      const int conc_b_scratch_id,
      const int temperature_scratch_id,
      const int diffusion_l_id,
      const int diffusion_a_id,
      const int diffusion_b_id,
      const int Mq_id,
      const std::vector<double>& Q_heat_transport,
      const std::vector<int> diffusion_precond_id,
      const std::string& avg_func_type,
      FreeEnergyStrategy* free_energy_strategy,
      CompositionStrategyMobilities* mobilities_strategy,
      CompositionDiffusionStrategy* diffusion_for_conc_in_phase
   );
   
   ~EBSCompositionRHSStrategy(){};

   void addFluxFromAntitrappingonPatch(
      hier::Patch& patch,
      const int phase_scratch_id,
      const int dphidt_id,
      const double alpha,
      const int flux_id);

   void computeFluxOnPatch(
      hier::Patch& patch,
      const int flux_id);
   
   void setDiffusionCoeff(
      const boost::shared_ptr< hier::PatchHierarchy > hierarchy,
      const double                                               time);
   
   //void printDiagnostics(const boost::shared_ptr<hier::PatchHierarchy > hierarchy);
   
private:
   unsigned short d_ncompositions;

   std::vector<double> d_Q_heat_transport;
   
   int d_phase_scratch_id;
   int d_eta_scratch_id;

   int d_conc_l_scratch_id;
   int d_conc_a_scratch_id;
   int d_conc_b_scratch_id;
   
   /*!
    * holds data for diffusion coefficients in composition equation
    * weighted by phase fraction
    */
   int d_diffusion_l_id;
   int d_diffusion_a_id;
   int d_diffusion_b_id;
   
   int d_Mq_id;
   
   int d_temperature_scratch_id;
   
   bool d_with_third_phase;
   
   bool d_with_gradT;
   
   bool d_with_diffusion_for_preconditioner;
   std::vector<int> d_diffusion_precond_id;
   
   // free energy needed to compute diffusion in each phase
   FreeEnergyStrategy* d_free_energy_strategy;

   CompositionDiffusionStrategy* d_diffusion_for_conc_in_phase;

   CompositionStrategyMobilities* d_mobilities_strategy;

   void addFluxFromGradTonPatch(
      hier::Patch& patch,
      const int temperature_id,
      const int flux_id);

   /*
    * Take sum of diffusion coefficients in each phase
    * to get diffusion used in preconditioner
    * (same diffusion coefficient as used in KKS equations)
    * Assumes coefficients in each pahse includes a phase fraction weight
    */
   void setDiffusionCoeffForPreconditioner(
      const boost::shared_ptr< hier::PatchHierarchy > hierarchy);
   void setDiffusionCoeffForPreconditionerOnPatch(
      boost::shared_ptr< pdat::SideData<double> > sd_d_l,
      boost::shared_ptr< pdat::SideData<double> > sd_d_a,
      boost::shared_ptr< pdat::SideData<double> > sd_d_b,
      const int depth_in_Dmatrix,
      boost::shared_ptr< pdat::SideData<double> > sd_d_coeff,
      const hier::Box& pbox,
      const int depth );
   void setDiffusionCoeffForTOnPatch(
      boost::shared_ptr< pdat::CellData<double> > cd_c,
      boost::shared_ptr< pdat::CellData<double> > cd_temp,
      boost::shared_ptr< pdat::SideData<double> > sd_d_coeff,
      const hier::Box& pbox );
    void setDiffusionCoeffForTOnPatch(
      boost::shared_ptr< pdat::CellData<double> > cd_c_l,
      boost::shared_ptr< pdat::CellData<double> > cd_c_a,
      boost::shared_ptr< pdat::CellData<double> > cd_c_b,
      boost::shared_ptr< pdat::CellData<double> > cd_temp,
      boost::shared_ptr< pdat::CellData<double> > cd_phi,
      boost::shared_ptr< pdat::CellData<double> > cd_eta,
      boost::shared_ptr< pdat::SideData<double> > mq,
      const hier::Box& pbox );
   void setDiffusionCoeffForT(
      const boost::shared_ptr< hier::PatchHierarchy > hierarchy,
      const int concentration_l_id,
      const int concentration_a_id,
      const int concentration_b_id,
      const int temperature_id);
};

#endif

