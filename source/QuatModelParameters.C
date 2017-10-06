// Copyright (c) 2018, Lawrence Livermore National Security, LLC.
// Produced at the Lawrence Livermore National Laboratory
// Written by M.R. Dorr, J.-L. Fattebert and M.E. Wickett
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE. 
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
// - Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the disclaimer below.
// - Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the disclaimer (as noted below) in the
//   documentation and/or other materials provided with the distribution.
// - Neither the name of the LLNS/LLNL nor the names of its contributors may be
//   used to endorse or promote products derived from this software without
//   specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
// LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
// OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
// 
#include "QuatModelParameters.h"

using namespace std;

static double def_val = tbox::IEEE::getSignalingNaN();

static
void readSpeciesCP(boost::shared_ptr<tbox::Database> cp_db,
                   map<short,double>& cp)
{
   double tmp=cp_db->getDouble("a");
   cp.insert( std::pair<short,double>(0,tmp) );
   if ( cp_db->keyExists( "b" ) )
   {
      tmp = cp_db->getDouble("b");
      cp.insert( std::pair<short,double>(1,tmp) );
   }
   if ( cp_db->keyExists( "dm2" ) )
   {
      tmp = cp_db->getDouble("dm2");
      cp.insert( std::pair<short,double>(-2,tmp) );
   }
}
      
QuatModelParameters::QuatModelParameters()
{

   d_H_parameter = def_val;
   d_epsilon_phase = def_val;
   d_epsilon_eta = def_val;
   d_epsilon_q = def_val;
   d_phase_mobility = def_val;
   d_eta_mobility = def_val;
   d_quat_mobility = def_val;
   d_min_eta_mobility = def_val;
   d_min_quat_mobility = def_val;
   d_max_quat_mobility = def_val;
   d_exp_scale_quat_mobility = def_val;
   d_q0_phase_mobility = def_val;
   d_q0_eta_mobility = def_val;
   d_phase_well_scale = def_val;
   d_eta_well_scale = def_val;
   d_free_energy_liquid = def_val;
   d_free_energy_solid_A = def_val;
   d_free_energy_solid_B = def_val;
   d_molar_volume_liquid = def_val;
   d_molar_volume_solid_A = def_val;
   d_molar_volume_solid_B = def_val;
   d_D_liquid = def_val;
   d_D_solid_A = def_val;
   d_D_solid_B = def_val;
   d_Q0_liquid = def_val;
   d_Q0_solid_A = def_val;
   d_Q0_solid_B = def_val;
   d_conc_mobility = def_val;
   d_thermal_diffusivity = def_val;
   d_latent_heat = def_val;
   d_meltingT = def_val;

   d_well_bias_alpha = 0.;
   d_well_bias_gamma = def_val;
   d_liquidus_slope = def_val;
   d_average_concentration = def_val;
   
   d_avg_func_type = "";
   d_diffq_avg_func_type = "";
   d_phase_well_func_type = "";

   d_orient_interp_func_type = "";
   d_phase_interp_func_type = "";
   d_diffusion_interp_func_type = "";
   d_eta_well_func_type = "";
   d_eta_interp_func_type = "";
   d_eta_well_func_type = "";
   d_quat_mobility_func_type = "";

   d_conc_avg_func_type = "";
   
   d_conc_model = UNDEFINED;
   d_conc_rhs_strategy = UNKNOWN;

   d_temperature_type = SCALAR;

   d_with_concentration = false;
   d_ncompositions=-1;

   d_with_phase = true;
   d_with_orientation = false;
   d_with_third_phase = false;
   d_with_heat_equation = false;
   d_with_steady_temperature = false;
   d_with_gradT = false;
   d_with_antitrapping = false;
   d_with_bias_well = false;
   d_with_rescaled_temperature = false;
   
   d_partition_coeff = "";
   d_phase_concentration_model = "";
   d_vd = -1.;
   d_with_velocity = false;
   d_use_diffs_to_compute_flux = false;
}

//=======================================================================

void QuatModelParameters::readMolarVolumes(boost::shared_ptr<tbox::Database> db)
{
   bool data_read=false;
   
   if ( db->keyExists( "molar_volume" ) ) {
      d_molar_volume_solid_A = db->getDouble( "molar_volume" );
      d_molar_volume_solid_B = d_molar_volume_solid_A;
      d_molar_volume_liquid  = d_molar_volume_solid_A;
      data_read=true;
   }
   else if ( db->keyExists( "molar_volume_solid_A" ) ) {
      d_molar_volume_liquid  = db->getDouble( "molar_volume_liquid" );
      d_molar_volume_solid_A = db->getDouble( "molar_volume_solid_A" );
      if ( db->keyExists( "molar_volume_solid_B" ) )
         d_molar_volume_solid_B = db->getDouble( "molar_volume_solid_B" );
      else
         d_molar_volume_solid_B = d_molar_volume_solid_A;
      data_read=true;
   }else if ( db->keyExists( "ConcentrationModel" ) ){
      readMolarVolumes( db->getDatabase( "ConcentrationModel" ) );
   }

   if( data_read )
      tbox::plog<<"Molar volume liquid [m^3/mol]: "<<d_molar_volume_liquid<<endl;
}

//=======================================================================

void QuatModelParameters::readNumberSpecies(boost::shared_ptr<tbox::Database> conc_db)
{
   int nspecies = conc_db->getIntegerWithDefault( "nspecies", 2 );
   d_ncompositions=nspecies-1;
}

//=======================================================================

void QuatModelParameters::readConcDB(boost::shared_ptr<tbox::Database> conc_db)
{
   d_with_concentration = true;
   
   //if concentration is ON, it means we have at least two species
   assert( d_ncompositions>0 );
   
   string conc_model =
      conc_db->getStringWithDefault( "model", "undefined" );
   if ( conc_model[0] == 'c' ) {
      d_conc_model = CALPHAD;
   }
   else if ( conc_model[0] == 'h' ) {
      d_conc_model = HBSM;
   }
   else if ( conc_model[0] == 'l' ) {
      d_conc_model = LINEAR;
   }
   else if ( conc_model[0] == 'i' ) {
      d_conc_model = INDEPENDENT; //energy independent of composition
   }
   else {
      TBOX_ERROR( "Error: unknown concentration model in QuatModelParameters" );
   }
   string conc_rhs_strategy =
      conc_db->getStringWithDefault( "rhs_form", "kks" );
   if ( conc_rhs_strategy[0] == 'k' ) {
      d_conc_rhs_strategy = KKS;
   }
   else if ( conc_rhs_strategy[0] == 'e' ) {
      d_conc_rhs_strategy = EBS;
   }
   else if ( conc_rhs_strategy[0] == 's' ) {
      d_conc_rhs_strategy = SPINODAL;
   }
   else if ( conc_rhs_strategy[0] == 'u' 
          || conc_rhs_strategy[0] == 'B'
          || conc_rhs_strategy[0] == 'b' ) {
      tbox::plog<<"Using Beckermann's model"<<endl;
      d_conc_rhs_strategy = Beckermann;
   }
   else {
      TBOX_ERROR( "Error: unknown concentration r.h.s. strategy" );
   }
   
   if ( d_conc_rhs_strategy == Beckermann ){
      d_D_liquid = conc_db->getDouble( "D_liquid" );
      d_D_solid_A = conc_db->getDouble( "D_solid_A" );
   }
   if ( d_conc_rhs_strategy == KKS ){
      d_D_liquid = conc_db->getDouble( "D_liquid" );
      d_Q0_liquid = conc_db->getDouble( "Q0_liquid" );

      if ( ! d_with_third_phase ) {
         if( conc_db->keyExists( "D_solid_A" ) )
            d_D_solid_A = conc_db->getDouble( "D_solid_A" );
         else
            d_D_solid_A = conc_db->getDouble( "D_solid" );
         if( conc_db->keyExists( "Q0_solid_A" ) )
            d_Q0_solid_A = conc_db->getDouble( "Q0_solid_A" );
         else
            d_Q0_solid_A = conc_db->getDouble( "Q0_solid" );
         d_D_solid_B = 0.;
         d_Q0_solid_B = 0.;
      }
      else {
         d_D_solid_A = conc_db->getDouble( "D_solid_A" );
         d_D_solid_B = conc_db->getDouble( "D_solid_B" );
         d_Q0_solid_A = conc_db->getDouble( "Q0_solid_A" );
         d_Q0_solid_B = conc_db->getDouble( "Q0_solid_B" );
      }
   }
   d_conc_mobility = conc_db->getDoubleWithDefault("mobility", 1.);
   
   if( d_conc_rhs_strategy == SPINODAL ){
      d_kappa = conc_db->getDouble("kappa");
   }

   d_conc_avg_func_type =
      conc_db->getStringWithDefault( "avg_func_type", d_avg_func_type );
   if ( d_conc_avg_func_type[0] != 'a' &&
        d_conc_avg_func_type[0] != 'h' ) {
      TBOX_ERROR( "Error: invalid value for avg_func_type" );
   }
   
   if ( conc_db->keyExists("gradT_Q0") ){
      d_with_gradT = true;
      //read heat of transport
      d_Q_heat_transport.resize(2);
      d_Q_heat_transport[0] = conc_db->getDouble( "gradT_Q0" );
      d_Q_heat_transport[1] = conc_db->getDouble( "gradT_Q1" );
   }else{
      d_with_gradT = false;
   }
   
   d_with_antitrapping = conc_db->getBoolWithDefault( "antitrapping", false );

   d_grand_potential = conc_db->getBoolWithDefault( "gc", false );
   
   d_partition_coeff = conc_db->getStringWithDefault( "partition_coeff", "none" );
   tbox::plog<<"Partition coefficient type: "<<d_partition_coeff<<endl;

   if( d_partition_coeff.compare("Aziz")==0 ){
      d_vd = conc_db->getDouble( "vd" );
      
      d_with_velocity = true;
   
      // a value of -1 for k_eq means it needs to be computed on the fly
      d_keq = conc_db->getDoubleWithDefault( "keq", -1. );
      tbox::plog<<"Aziz partition coefficient with Keq: "<<d_keq<<endl;
   }
   if( d_partition_coeff.compare("uniform")==0 ){
      d_keq = conc_db->getDouble( "keq" );
      tbox::plog<<"Uniform Keq: "<<d_keq<<endl;      
   }
   
   assert( d_partition_coeff.compare("Aziz")==0 
        || d_partition_coeff.compare("uniform")==0 
        || d_partition_coeff.compare("none")==0 );
   
   string default_model="none";
   if( d_conc_model==CALPHAD || d_conc_model==HBSM )default_model="kks";
   d_phase_concentration_model = 
      conc_db->getStringWithDefault( "phase_concentration_model", default_model );
   assert( d_phase_concentration_model.compare("none")==0 
        || d_phase_concentration_model.compare("kks")==0 
        || d_phase_concentration_model.compare("partition")==0 );
   tbox::plog<<"phase_concentration_model: "<<d_phase_concentration_model<<endl;

   if( d_phase_concentration_model.compare("partition")==0 )
      assert( d_partition_coeff.compare("none")!=0 );

   if( d_conc_model==LINEAR ){
      assert( d_meltingT == d_meltingT );

      d_liquidus_slope = conc_db->getDoubleWithDefault( "liquidus_slope", 0. );
      if( fabs(d_liquidus_slope)>0. )
         d_average_concentration    = conc_db->getDouble( "average_concentration" );
      d_well_bias_alpha = conc_db->getDouble( "alpha" );
      if( d_well_bias_alpha>0. ){
         d_with_bias_well = true;
      }
      d_well_bias_gamma = conc_db->getDoubleWithDefault( "gamma", -1. );
      if( d_well_bias_gamma<0. ){
         d_bias_well_beckermann = true;
      }else{
         d_bias_well_beckermann = false;
      }
      if( d_with_rescaled_temperature ){
         tbox::plog<<"Rescale liquidus_slope and gamma..."<<endl;
         d_liquidus_slope /= d_meltingT;
         d_well_bias_gamma *= d_meltingT;
      }
   }
}

//=======================================================================

void QuatModelParameters::readVisitOptions(
   boost::shared_ptr<tbox::Database> visit_db)
{
   d_extra_visit_output =
      visit_db->getBoolWithDefault( "extra_output", false );
   d_visit_energy_output =
      visit_db->getBoolWithDefault( "energy_output", false );
   d_visit_grain_output =
      visit_db->getBoolWithDefault( "grain_output", false );
   d_with_velocity =
      visit_db->getBoolWithDefault( "velocity_output", d_with_velocity);
   d_rhs_visit_output =
      visit_db->getBoolWithDefault( "rhs_output", false );
}

//=======================================================================

void QuatModelParameters::readTemperatureModel(
   boost::shared_ptr<tbox::Database> model_db)
{
   boost::shared_ptr<tbox::Database> temperature_db;
   string temperature_type = "";
   if ( model_db->keyExists( "Temperature" ) ) {
      temperature_db = model_db->getDatabase( "Temperature" );
      temperature_type =
         temperature_db->getStringWithDefault( "type", "scalar" );
   }else{
      temperature_db = model_db;
      temperature_type =
         temperature_db->getStringWithDefault( "temperature_type", "scalar" );
   }
   
   if ( temperature_type[0] != 's' &&  // scalar
        temperature_type[0] != 'S' &&
        temperature_type[0] != 'g' &&  // Gaussian
        temperature_type[0] != 'G' &&
        temperature_type[0] != 'c' &&  // constant
        temperature_type[0] != 'C'&&
        temperature_type[0] != 'h' &&  // heat equation
        temperature_type[0] != 'H' ) {
      TBOX_ERROR( "Error: invalid value for temperature_type" );
   }

   d_meltingT = temperature_db->getDoubleWithDefault( "meltingT", -1. ); // in [K]

   if ( temperature_type[0] == 's' ||
        temperature_type[0] == 'S' ) {
      
      d_temperature_type   = SCALAR;
      d_with_heat_equation = false;
      
   }
   else if(temperature_type[0] == 'g' ||
           temperature_type[0] == 'G'){
      
      d_temperature_type   = GAUSSIAN;
      d_with_heat_equation = false;
      
   }
   else if(temperature_type[0] == 'h' ||
           temperature_type[0] == 'H'){
      
      assert( d_molar_volume_liquid==d_molar_volume_liquid );
      assert( d_molar_volume_liquid>0. );
      assert( d_molar_volume_liquid<1.e15 );
      
      d_temperature_type   = CONSTANT;
      d_with_heat_equation = true;

      tbox::plog<<"Read heat equation parameters..."<<endl;

      string method = temperature_db->getStringWithDefault( "equation_type", "steady" );
      
      // heat capacity value
      boost::shared_ptr<tbox::Database> cp_db = temperature_db->getDatabase( "cp" );
      map<short,double> empty_map;
     
      d_with_rescaled_temperature = ( ( method != "steady" ) && (d_meltingT>0.) );
      if( d_with_rescaled_temperature )
         tbox::plog<<"Solve temperature equation with rescaled T"<<endl; 
      d_cp.push_back(empty_map);
      readSpeciesCP(cp_db->getDatabase( "SpeciesA" ), d_cp[0]);
      if( d_ncompositions>0 )
      {
         assert( cp_db->keyExists( "SpeciesB" ) );
         d_cp.push_back(empty_map);
         readSpeciesCP(cp_db->getDatabase("SpeciesB" ), d_cp[1]);
      }
      if( d_ncompositions>1 )
      {
         assert( cp_db->keyExists( "SpeciesC" ) );
         d_cp.push_back(empty_map);
         readSpeciesCP(cp_db->getDatabase("SpeciesC" ), d_cp[2]);
      }
      
      tbox::plog<<"Cp for each species: "<<endl;
      for(vector< map<short,double> >::iterator it=d_cp.begin(); it!=d_cp.end(); ++it)
      {
         for(map<short,double>::iterator itm= it->begin(); itm!=it->end(); ++itm)
         {
            itm->second *= ( 1.e-6 / d_molar_volume_liquid ); // conversion from [J/mol*K] to [pJ/(mu m)^3*K]
            tbox::plog<<"Cp [pJ/(mu m)^3*K]: "<<itm->second<<endl;
         }
      }
      if( d_with_rescaled_temperature ){
         for(vector< map<short,double> >::iterator it=d_cp.begin(); it!=d_cp.end(); ++it)
         {
            for(map<short,double>::iterator itm= it->begin(); itm!=it->end(); ++itm)
            {
               itm->second *= d_meltingT;
               tbox::plog<<"rescaled Cp: "<<itm->second<<endl;
            }
         }
      }
      
      d_thermal_diffusivity = temperature_db->getDouble( "thermal_diffusivity" );
      tbox::plog<<"Thermal diffusivity [cm^2/s]:        "<<d_thermal_diffusivity<<endl;
      d_thermal_diffusivity *= 1.e8; // cm^2/s -> um^2/s
      tbox::plog<<"Thermal diffusivity [um^2/s]:        "<<d_thermal_diffusivity<<endl;

      if ( method == "steady" ) {
         d_with_steady_temperature = true;
         
         if ( temperature_db->keyExists( "heat_source_type" ) ) {
            d_heat_source_type = temperature_db->getString( "heat_source_type" );
            if( d_heat_source_type=="composition")
            {
               size_t nterms=temperature_db->getArraySize("source");
               d_T_source.resize(nterms);
               temperature_db->getDoubleArray("source",&d_T_source[0],nterms);
               for(size_t i=0;i<nterms;++i)
                  d_T_source[i] *= ( 1.e-6 / d_molar_volume_liquid ); // conversion from [J/mol] to [pJ/(mu m)^3]
            }
            else
            {
               assert( d_heat_source_type=="gaussian");
            }
         }
      }else{
         d_with_steady_temperature = false;
      }

      if( d_with_rescaled_temperature ){   //rescale units
         tbox::plog<<"Rescale thermal diffusivity by factor 1./d_meltingT = "<<1./d_meltingT<<endl;
         d_thermal_diffusivity/=d_meltingT;
         tbox::plog<<"Rescaled thermal diffusivity: "<<d_thermal_diffusivity<<endl;
         d_H_parameter *= d_meltingT;
      }
      
   }
   else {
      d_temperature_type = CONSTANT;
   }

   if( temperature_db->keyExists("latent_heat" ) ){
      d_latent_heat = temperature_db->getDouble( "latent_heat" ); // in [J/mol]
      d_latent_heat *= ( 1.e-6 / d_molar_volume_liquid ); // conversion from [J/mol] to [pJ/(mu m)^3]
      tbox::plog<<"Latent heat [pJ/(mu m)^3]:           "<<d_latent_heat<<endl;
      assert( d_latent_heat>0. );
      assert( d_latent_heat<1.e32 );
   }else{
      d_latent_heat = def_val;
   }

}

//=======================================================================

void QuatModelParameters::initializeOrientation(
   boost::shared_ptr<tbox::Database> model_db)
{
   
   if ( model_db->keyExists( "orient_mobility" ) ) {
      d_quat_mobility = model_db->getDouble( "orient_mobility" );
   }
   else if ( model_db->keyExists( "quat_mobility" ) ) {
      d_quat_mobility = model_db->getDouble( "quat_mobility" );
      printDeprecated( "quat_mobility", "orient_mobility" );
   }
   else if ( model_db->keyExists( "tau_quat" ) ) {
      double tau = model_db->getDouble( "tau_quat" );
      d_quat_mobility = 1. / tau;
      printDeprecated( "tau_quat", "orient_mobility" );
   }
   else {
      TBOX_ERROR( "Error: quaternion mobility not specified" );
   }      

   d_min_quat_mobility = 1.e-6;
   if ( model_db->keyExists( "min_orient_mobility" ) ) {
      d_min_quat_mobility = model_db->getDouble( "min_orient_mobility" );
   }
   else if ( model_db->keyExists( "min_quat_mobility" ) ) {
      d_min_quat_mobility = model_db->getDouble( "min_quat_mobility" );
      printDeprecated( "min_quat_mobility", "min_orient_mobility" );
   }
      
   if ( model_db->keyExists( "epsilon_orient" ) ) {
      d_epsilon_q = model_db->getDouble( "epsilon_orient" );
   }
   else if ( model_db->keyExists( "epsilon_q" ) ) {
      d_epsilon_q = model_db->getDouble( "epsilon_q" );
      printDeprecated( "epsilon_q", "epsilon_orient" );
   }
   else if ( model_db->keyExists( "epsilon_quat" ) ) {
      d_epsilon_q = model_db->getDouble( "epsilon_quat" );
      printDeprecated( "epsilon_quat", "epsilon_orient" );
   }
   else {
      TBOX_ERROR( "Error: epsilon_quat not specified" );
   }

   d_quat_grad_floor = model_db->getDoubleWithDefault(
      "orient_grad_floor", 1.e-2 );
   if ( !model_db->keyExists( "orient_grad_floor" ) )
   if ( model_db->keyExists( "quat_grad_floor" ) ) {
      d_quat_grad_floor = model_db->getDouble( "quat_grad_floor" );
      printDeprecated( "quat_grad_floor", "orient_grad_floor" );
   }

   // options for "smooth floor" are: 1./|grad(q)| =...
   // "max": max( 1./|grad(q)|, 1./d_quat_grad_floor )
   // "tanh": tanh(grad_norm/d_quat_grad_floor)/ grad_norm
   // "sqrt": 1./sqrt( |grad(q)|**2+d_quat_grad_floor**2 )
   d_quat_grad_floor_type = model_db->getStringWithDefault( "orient_grad_floor_type", "max" );
   if ( model_db->keyExists( "quat_smooth_floor" ) ) {
      printDeprecated( "quat_smooth_floor", "orient_grad_floor_type" );
   }
   if ( model_db->keyExists( "orient_smooth_floor" ) ) {
      printDeprecated( "quat_smooth_floor", "orient_grad_floor_type" );
   }

   d_quat_grad_modulus_type = model_db->getStringWithDefault( "quat_grad_modulus_type", "cells" );

   if ( model_db->keyExists( "diff_interp_func_type" ) ) {
      d_orient_interp_func_type =
         model_db->getString( "diff_interp_func_type" );
      printDeprecated( "diff_interp_func_type", "orient_interp_func_type" );
   }else{
      d_orient_interp_func_type =
         model_db->getStringWithDefault( "orient_interp_func_type", "quadratic" );
   }
   if ( d_orient_interp_func_type[0] != 'q' &&
        d_orient_interp_func_type[0] != 'Q' &&
        d_orient_interp_func_type[0] != 'w' &&
        d_orient_interp_func_type[0] != 'W' &&
        d_orient_interp_func_type[0] != 'p' &&
        d_orient_interp_func_type[0] != 'P' &&
        d_orient_interp_func_type[0] != 'l' &&
        d_orient_interp_func_type[0] != 'L' &&
        d_orient_interp_func_type[0] != 't' &&
        d_orient_interp_func_type[0] != 'T' &&
        d_orient_interp_func_type[0] != 's' &&
        d_orient_interp_func_type[0] != 'S' &&
        d_orient_interp_func_type[0] != '3' &&
        d_orient_interp_func_type[0] != 'c' &&
        d_orient_interp_func_type[0] != 'C' ) {
      TBOX_ERROR( "Error: invalid value for orient_interp_func_type" );
   }

   if ( model_db->keyExists( "quat_mobility_func_type" ) ) {
      d_quat_mobility_func_type =
         model_db->getString( "quat_mobility_func_type" );
      printDeprecated(
         "quat_mobility_func_type",
         "orient_mobility_func_type" );
   }else{
      d_quat_mobility_func_type =
         model_db->getStringWithDefault( "orient_mobility_func_type", "pbg");
   }
   if ( d_quat_mobility_func_type[0] != 'p' &&
        d_quat_mobility_func_type[0] != 'P' &&
        d_quat_mobility_func_type[0] != 'i' &&
        d_quat_mobility_func_type[0] != 'I' &&
        d_quat_mobility_func_type[0] != 'e' &&
        d_quat_mobility_func_type[0] != 'E' ) {
      TBOX_ERROR( "Error: invalid value for orient_mobility_func_type" );
   }

   if ( d_quat_mobility_func_type[0] == 'i' ||
        d_quat_mobility_func_type[0] == 'I' ) {
      d_max_quat_mobility = 1.e6;
      if ( model_db->keyExists( "max_orient_mobility" ) ) {
         d_max_quat_mobility = model_db->getDouble( "max_orient_mobility" );
      }
      else if ( model_db->keyExists( "max_quat_mobility" ) ) {
         d_max_quat_mobility = model_db->getDouble( "max_quat_mobility" );
         printDeprecated( "max_quat_mobility", "max_orient_mobility" );
      }
   }
      
   if ( d_quat_mobility_func_type[0] == 'e' ||
        d_quat_mobility_func_type[0] == 'E' ) {
      d_exp_scale_quat_mobility = 1.e6;
      if ( model_db->keyExists( "exp_scale_orient_mobility" ) ) {
         d_exp_scale_quat_mobility = model_db->getDouble( "exp_scale_orient_mobility" );
      }
      else if ( model_db->keyExists( "exp_scale_quat_mobility" ) ) {
         d_exp_scale_quat_mobility = model_db->getDouble( "exp_scale_quat_mobility" );
         printDeprecated( "exp_scale_quat_mobility", "exp_scale_orient_mobility" );
      }
   }
}

//=======================================================================

void QuatModelParameters::initializeEta(boost::shared_ptr<tbox::Database> model_db)
{
   d_epsilon_eta = model_db->getDouble( "epsilon_eta" );

   d_eta_mobility = model_db->getDouble( "eta_mobility" );
   d_min_eta_mobility = model_db->getDoubleWithDefault(
      "min_eta_mobility", d_eta_mobility );

   d_q0_eta_mobility = model_db->getDoubleWithDefault(
      "q0_eta_mobility", 0.0 );

   d_eta_well_scale = model_db->getDoubleWithDefault( "eta_well_scale", 1.0 );

   d_eta_well_func_type =
      model_db->getStringWithDefault( "eta_well_func_type", "double" );
   if ( d_eta_well_func_type[0] != 's' &&
        d_eta_well_func_type[0] != 'S' &&
        d_eta_well_func_type[0] != 'd' &&
        d_eta_well_func_type[0] != 'D' ) {
      TBOX_ERROR( "Error: invalid value for eta_well_func_type" );
   }

   d_eta_interp_func_type =
      model_db->getStringWithDefault( "eta_interp_func_type", "pbg" );
   if ( d_eta_interp_func_type[0] != 'q' &&
        d_eta_interp_func_type[0] != 'Q' &&
        d_eta_interp_func_type[0] != 'h' &&
        d_eta_interp_func_type[0] != 'H' &&
        d_eta_interp_func_type[0] != 'p' &&
        d_eta_interp_func_type[0] != 'P' ) {
      TBOX_ERROR( "Error: invalid value for eta_interp_func_type" );
   }
}

//=======================================================================

void QuatModelParameters::readModelParameters(boost::shared_ptr<tbox::Database> model_db)
{
   // Set d_H_parameter to negative value, to turn off orientation terms
   d_H_parameter = model_db->getDoubleWithDefault( "H_parameter", -1. );
   if ( d_H_parameter >= 0. ) {
      d_with_orientation = true;
   }
   else {
      d_with_orientation = false;
   }

   // Interface width
   if ( model_db->keyExists( "epsilon_phi" ) ) {
      d_epsilon_phase = model_db->getDouble( "epsilon_phi" );
   }
   else if ( model_db->keyExists( "epsilon_phase" ) ) {
      d_epsilon_phase = model_db->getDouble( "epsilon_phase" );
      printDeprecated( "epsilon_phase", "epsilon_phi" );
   }
   else if ( model_db->keyExists( "epsilon_parameter" ) ) {
      d_epsilon_phase = model_db->getDouble( "epsilon_parameter" );
      printDeprecated( "epsilon_parameter", "epsilon_phi" );
   }
   else {
      d_with_phase = false;
      tbox::pout<<"No epsilon specified -> run without phase..."<<endl;
   }

   // Mobility
   if( d_with_phase ){
   if ( model_db->keyExists( "phi_mobility" ) ) {
      d_phase_mobility = model_db->getDouble( "phi_mobility" );
   }
   else if ( model_db->keyExists( "phase_mobility" ) ) {
      d_phase_mobility = model_db->getDouble( "phase_mobility" );
      printDeprecated( "phase_mobility", "phi_mobility" );
   }
   else if ( model_db->keyExists( "tau_phase" ) ) {
      double tau = model_db->getDouble( "tau_phase" );
      d_phase_mobility = 1. / tau;
      printDeprecated( "tau_phase", "phi_mobility" );
   }
   else {
      TBOX_ERROR( "Error: phi_mobility not specified" );
   }      
   }

   d_q0_phase_mobility = model_db->getDoubleWithDefault(
      "q0_phi_mobility", 0.0 );

   // Well energy
   if( d_with_phase ){
   d_phase_well_scale = 1.0;
   if ( model_db->keyExists( "phi_well_scale" ) ) {
      d_phase_well_scale = model_db->getDouble( "phi_well_scale" );
   }
   else if ( model_db->keyExists( "scale_energy_well" ) ) {
      d_phase_well_scale = model_db->getDouble( "scale_energy_well" );
      printDeprecated( "scale_energy_well", "phi_well_scale" );
   }

   d_phase_well_func_type = "double";
   if ( model_db->keyExists( "phi_well_func_type" ) ) {
      d_phase_well_func_type =
         model_db->getString( "phi_well_func_type" );
   }
   else if ( model_db->keyExists( "energy_well_func_type" ) ) {
      d_phase_well_func_type =
         model_db->getString( "energy_well_func_type" );
      printDeprecated( "energy_well_func_type", "phi_well_func_type" );
   }
   if ( d_phase_well_func_type[0] != 's' &&
        d_phase_well_func_type[0] != 'S' &&
        d_phase_well_func_type[0] != 'd' &&
        d_phase_well_func_type[0] != 'D' ) {
      TBOX_ERROR( "Error: invalid value for phi_well_func_type" );
   }
   }
   
   if ( !model_db->keyExists( "ConcentrationModel" ) ){
      if ( model_db->keyExists( "bias_well_alpha" ) ) {
         d_well_bias_alpha = model_db->getDouble( "bias_well_alpha" );
         if( d_well_bias_alpha>0. ){
            d_with_bias_well = true;
         }
      }
      if ( model_db->keyExists( "bias_well_gamma" ) ) {
         d_well_bias_gamma = model_db->getDouble( "bias_well_gamma" );
      }
   }
   

   // Molar volumes
   readMolarVolumes(model_db);

   // Interpolation

   d_phase_interp_func_type = "pbg";
   if ( model_db->keyExists( "phi_interp_func_type" ) ) {
      d_phase_interp_func_type =
         model_db->getString( "phi_interp_func_type" );
   }
   else if ( model_db->keyExists( "energy_interp_func_type" ) ) {
      d_phase_interp_func_type =
         model_db->getString( "energy_interp_func_type" );
      printDeprecated( "energy_interp_func_type", "phi_interp_func_type" );
   }
   if ( d_phase_interp_func_type[0] != 'q' &&
        d_phase_interp_func_type[0] != 'Q' &&
        d_phase_interp_func_type[0] != 'h' &&
        d_phase_interp_func_type[0] != 'H' &&
        d_phase_interp_func_type[0] != 'l' &&
        d_phase_interp_func_type[0] != 'L' &&
        d_phase_interp_func_type[0] != 'c' &&
        d_phase_interp_func_type[0] != 'C' &&
        d_phase_interp_func_type[0] != 'p' &&
        d_phase_interp_func_type[0] != 'P' ) {
      TBOX_ERROR( "Error: invalid value for phase_interp_func_type" );
   }
   
   d_diffusion_interp_func_type =
      model_db->getStringWithDefault( "diffusion_interp_func_type", d_phase_interp_func_type );

   // Currently "arithmetic" or "harmonic"
   // arithmetic: (x1+x2)/2
   // harmonic:   2/(1/x1+1/x2)
   d_avg_func_type =
      model_db->getStringWithDefault( "avg_func_type", "harmonic" );
   if ( d_avg_func_type[0] != 'a' &&
        d_avg_func_type[0] != 'h' ) {
      TBOX_ERROR( "Error: invalid value for avg_func_type" );
   }

   d_diffq_avg_func_type =
      model_db->getStringWithDefault( "diffq_avg_func_type", d_avg_func_type );
   if ( d_diffq_avg_func_type[0] != 'a' &&
        d_diffq_avg_func_type[0] != 'h' ) {
      TBOX_ERROR( "Error: invalid value for avg_func_type" );
   }
   
   d_use_diffs_to_compute_flux =
      model_db->getBoolWithDefault( "use_diffs_to_compute_flux", false );      
   d_stencil_type =
      model_db->getStringWithDefault( "stencil_type", "normal");

   d_with_third_phase =
      model_db->getBoolWithDefault( "three_phase", false );
   
   
   if ( with_third_phase() )
      initializeEta(model_db);

   if ( with_orientation() )
      initializeOrientation(model_db);

   //we need to know how many species we have before reading temperature model
   //which may include species dependent Cp 
   if ( model_db->keyExists( "ConcentrationModel" ) ) {
      readNumberSpecies( model_db->getDatabase( "ConcentrationModel" ) );
   }else{
      d_ncompositions=0;
   } 

   readTemperatureModel(model_db);
 
   if ( model_db->keyExists( "ConcentrationModel" ) ) {      
      boost::shared_ptr<tbox::Database> conc_db(model_db->getDatabase( "ConcentrationModel" ));
      readConcDB(conc_db);
   }
}

void QuatModelParameters::readFreeEnergies(boost::shared_ptr<tbox::Database> model_db)
{
   boost::shared_ptr<tbox::Database> db( 
      model_db->keyExists( "FreeEnergyModel" ) ?
      model_db->getDatabase( "FreeEnergyModel" ) :
      model_db );

   d_free_energy_type = db->getStringWithDefault("type","none");

   if( d_free_energy_type[0]=='s'){
      d_free_energy_liquid =
         db->getDouble( "free_energy_liquid" );

      if ( ! with_third_phase() ) {
         d_free_energy_solid_A =
            db->getDouble( "free_energy_solid" );
      }
      else {
         d_free_energy_solid_A =
            db->getDouble( "free_energy_solid_A" );
         d_free_energy_solid_B =
            db->getDouble( "free_energy_solid_B" );
      }
   }else{
      if( d_free_energy_type[0]=='l'){
         d_free_energy_solid_A = 0.;
      }
   }
}

double QuatModelParameters::quatMobilityScaleFactor()const
{
   if ( d_quat_mobility_func_type[0] == 'e' ||
        d_quat_mobility_func_type[0] == 'E' ) {
      return d_exp_scale_quat_mobility;
   }
   if ( d_quat_mobility_func_type[0] == 'i' ||
        d_quat_mobility_func_type[0] == 'I' ) {
      return d_max_quat_mobility;
   }
   
   return tbox::IEEE::getSignalingNaN();
}
