/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   RefineSchedule's implementation of PatchHierarchy
 *
 ************************************************************************/

#ifndef included_xfer_RefineScheduleConnectorWidthRequestor
#define included_xfer_RefineScheduleConnectorWidthRequestor

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/PatchHierarchy.h"

namespace SAMRAI {
namespace xfer {

/*!
 * @brief Implementation of the strategy class
 * hier::PatchHierarchy::ConnectorWidthRequestorStrategy to tell the
 * hier::PatchHierarchy how wide RefineSchedule needs Connectors
 * between hierarchy levels to be.
 */
class RefineScheduleConnectorWidthRequestor:
   public hier::PatchHierarchy::ConnectorWidthRequestorStrategy
{

public:
   /*!
    * @brief Constructor.
    */
   RefineScheduleConnectorWidthRequestor();

   /*!
    * @brief Compute Connector widths that this class requires in
    * order to work properly on a given hierarchy.
    *
    * Implements the pure virtual method
    * hier::PatchHierarchy::ConnectorWidthRequestorStrategy::computeRequiredConnectorWidths()
    *
    * @param[out] self_connector_widths Array of widths for Connectors
    * from a level to itself.
    *
    * @param[out] fine_connector_widths Array of widths for Connectors
    * from a level to the next finer level.
    *
    * @param[in]  patch_hierarchy
    */
   void
   computeRequiredConnectorWidths(
      std::vector<hier::IntVector>& self_connector_widths,
      std::vector<hier::IntVector>& fine_connector_widths,
      const hier::PatchHierarchy& patch_hierarchy) const;

   /*!
    * @brief Set the factor by which the ghost data widths in the
    * PatchHierarchy are multiplied to increase their effective values.
    * for the purpose of handling more ghost data than registered.
    *
    * We support values of ghost data width that are larger than those
    * registered by some factor.  An example of the need is when we
    * fill a PatchLevel that is a coarsened version of a current level
    * in the hierarchy.  Coarsening the level effectively increases
    * the size of the ghost data region because the coarse cells are
    * bigger.
    *
    * Note that the Connector widths are not linear functions of this
    * factor.  Increasing the factor is NOT the same as multiplying
    * the Connector widths by the same factor.
    *
    * @param[in] factor By default, @c factor=1.
    */
   void
   setGhostCellWidthFactor(
      int factor);

private:
   /*!
    * @brief Set up things for the entire class.
    *
    * Only called by StartupShutdownManager.
    */
   static void
   initializeCallback()
   {
      hier::PatchHierarchy::registerAutoConnectorWidthRequestorStrategy(
         s_auto_registered_connector_width_requestor);
   }

   /*!
    * Free static timers.
    *
    * Only called by StartupShutdownManager.
    */
   static void
   finalizeCallback()
   {
   }

   /*!
    * @brief Static object that is auto-registered in PatchHierarchy
    * by default.
    */
   static RefineScheduleConnectorWidthRequestor
      s_auto_registered_connector_width_requestor;

   static tbox::StartupShutdownManager::Handler
      s_initialize_finalize_handler;

   /*
    * @brief The factor by which the ghost data widths in the
    * PatchHierarchy are multiplied to increase their effective
    * values.  for the purpose of handling more ghost data than
    * registered.
    */
   int d_gcw_factor;

};

}
}

#endif
