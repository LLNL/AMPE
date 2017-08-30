#ifndef included_CALPHADFreeEnergyFunctions
#define included_CALPHADFreeEnergyFunctions

#include "FreeEnergyFunctions.h"

class CALPHADFreeEnergyFunctions:
   public FreeEnergyFunctions
{
public:

   CALPHADFreeEnergyFunctions(){};

   virtual ~CALPHADFreeEnergyFunctions(){};

   virtual int computePhaseConcentrations(
      const double temperature, const double* conc, const double phi, const double eta,
      double* x)=0;

};

#endif

