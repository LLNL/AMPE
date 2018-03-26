/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Registry of PersistentOverlapConnectorss incident from a common BoxLevel.
 *
 ************************************************************************/
#ifndef included_hier_PersistentOverlapConnectors_C
#define included_hier_PersistentOverlapConnectors_C

#include "SAMRAI/hier/PersistentOverlapConnectors.h"
#include "SAMRAI/hier/Connector.h"
#include "SAMRAI/hier/BoxLevel.h"
#include "SAMRAI/hier/OverlapConnectorAlgorithm.h"
#include "SAMRAI/tbox/InputManager.h"

namespace SAMRAI {
namespace hier {

char PersistentOverlapConnectors::s_check_created_connectors('\0');
char PersistentOverlapConnectors::s_check_accessed_connectors('\0');
bool PersistentOverlapConnectors::s_always_create_missing_connector(true);

/*
 ************************************************************************
 * This private constructor can only be used by the friend
 * class BoxLevel.
 ************************************************************************
 */
PersistentOverlapConnectors::PersistentOverlapConnectors(
   const BoxLevel& my_mapped_box_level):
   d_my_mapped_box_level(my_mapped_box_level)
{
   if (s_check_created_connectors == '\0') {
      boost::shared_ptr<tbox::Database> idb(
         tbox::InputManager::getInputDatabase());
      if (idb && idb->isDatabase("PersistentOverlapConnectors")) {
         boost::shared_ptr<tbox::Database> rsdb(
            idb->getDatabase("PersistentOverlapConnectors"));

         const bool check_created_connectors(
            rsdb->getBoolWithDefault("check_created_connectors", false));
         s_check_created_connectors = check_created_connectors ? 'y' : 'n';

         const bool check_accessed_connectors(
            rsdb->getBoolWithDefault("check_accessed_connectors", false));
         s_check_accessed_connectors = check_accessed_connectors ? 'y' : 'n';

         s_always_create_missing_connector =
            rsdb->getBoolWithDefault("always_create_missing_connector",
               s_always_create_missing_connector);
      }

   }
}

/*
 ************************************************************************
 *
 ************************************************************************
 */
PersistentOverlapConnectors::~PersistentOverlapConnectors()
{
   clear();
}

/*
 ************************************************************************
 * Create Connector using global search for edges.
 ************************************************************************
 */
const Connector&
PersistentOverlapConnectors::createConnector(
   const BoxLevel& head,
   const IntVector& connector_width)
{
   TBOX_ASSERT(d_my_mapped_box_level.isInitialized());
   TBOX_ASSERT(head.isInitialized());

   for (int i = 0; i < d_cons_from_me.size(); ++i) {
      if (&d_cons_from_me[i]->getHead() == &head &&
          d_cons_from_me[i]->getConnectorWidth() == connector_width) {
         TBOX_ERROR(
            "PersistentOverlapConnectors::createConnector:\n"
            << "Cannot create duplicate Connectors.");
      }
   }

   Connector* new_connector = new Connector(
         d_my_mapped_box_level,
         head,
         connector_width,
         BoxLevel::DISTRIBUTED);
   OverlapConnectorAlgorithm oca;
   oca.findOverlaps(*new_connector, head.getGlobalizedVersion());

   d_cons_from_me.push_back(new_connector);
   head.getPersistentOverlapConnectors().d_cons_to_me.push_back(new_connector);

   return *d_cons_from_me.back();
}

/*
 ************************************************************************
 * Create Connector with user-provided for edges.
 ************************************************************************
 */
const Connector&
PersistentOverlapConnectors::createConnector(
   const BoxLevel& head,
   const IntVector& connector_width,
   const Connector& relationships)
{
   TBOX_ASSERT(d_my_mapped_box_level.isInitialized());
   TBOX_ASSERT(head.isInitialized());

   for (int i = 0; i < d_cons_from_me.size(); ++i) {
      TBOX_ASSERT(d_cons_from_me[i]->isFinalized());
      TBOX_ASSERT(d_cons_from_me[i]->getBase().isInitialized());
      TBOX_ASSERT(d_cons_from_me[i]->getHead().isInitialized());
      TBOX_ASSERT(d_cons_from_me[i]->getBase().getBoxLevelHandle());
      TBOX_ASSERT(d_cons_from_me[i]->getHead().getBoxLevelHandle());
      TBOX_ASSERT(&d_cons_from_me[i]->getBase() ==
         &d_cons_from_me[i]->getBase().getBoxLevelHandle()->
         getBoxLevel());
      TBOX_ASSERT(&d_cons_from_me[i]->getHead() ==
         &d_cons_from_me[i]->getHead().getBoxLevelHandle()->
         getBoxLevel());
      if (&(d_cons_from_me[i]->getHead()) == &head &&
          d_cons_from_me[i]->getConnectorWidth() == connector_width) {
         TBOX_ERROR(
            "PersistentOverlapConnectors::createConnector:\n"
            << "Cannot create duplicate Connectors.");
      }
   }

   Connector* new_connector = new Connector(relationships);
   new_connector->setBase(d_my_mapped_box_level);
   new_connector->setHead(head);
   new_connector->setWidth(connector_width, true);
   if (s_check_created_connectors == 'y') {
      OverlapConnectorAlgorithm oca;
      TBOX_ASSERT(oca.checkOverlapCorrectness(*new_connector) == 0);
   }

   d_cons_from_me.push_back(new_connector);
   head.getPersistentOverlapConnectors().d_cons_to_me.push_back(new_connector);

   new_connector = NULL; // Help Insure++ avoid false positive dangling pointer.

   return *d_cons_from_me.back();
}

/*
 ************************************************************************
 * Cache the user-provided Connector.
 ************************************************************************
 */
void
PersistentOverlapConnectors::cacheConnector(
   const BoxLevel& head,
   Connector* connector)
{
   TBOX_ASSERT(d_my_mapped_box_level.isInitialized());

   for (int i = 0; i < d_cons_from_me.size(); ++i) {
      TBOX_ASSERT(d_cons_from_me[i]->isFinalized());
      TBOX_ASSERT(d_cons_from_me[i]->getBase().isInitialized());
      TBOX_ASSERT(d_cons_from_me[i]->getHead().isInitialized());
      TBOX_ASSERT(d_cons_from_me[i]->getBase().getBoxLevelHandle());
      TBOX_ASSERT(d_cons_from_me[i]->getHead().getBoxLevelHandle());
      TBOX_ASSERT(&d_cons_from_me[i]->getBase() ==
         &d_cons_from_me[i]->getBase().getBoxLevelHandle()->
         getBoxLevel());
      TBOX_ASSERT(&d_cons_from_me[i]->getHead() ==
         &d_cons_from_me[i]->getHead().getBoxLevelHandle()->
         getBoxLevel());
      if (&(d_cons_from_me[i]->getHead()) == &head &&
          d_cons_from_me[i]->getConnectorWidth() == connector->getConnectorWidth()) {
         TBOX_ERROR(
            "PersistentOverlapConnectors::createConnector:\n"
            << "Cannot create duplicate Connectors.");
      }
   }

   connector->setBase(d_my_mapped_box_level);
   connector->setHead(head, true);

   if (s_check_created_connectors == 'y') {
      OverlapConnectorAlgorithm oca;
      if (oca.checkOverlapCorrectness(*connector) != 0) {
         TBOX_ERROR("PersistentOverlapConnectors::cacheConnector errror:\n"
                    <<"Bad overlap Connector found.");
      }
   }

   d_cons_from_me.push_back(connector);
   head.getPersistentOverlapConnectors().d_cons_to_me.push_back(connector);

   return;
}

/*
 ************************************************************************
 *
 ************************************************************************
 */
const Connector&
PersistentOverlapConnectors::findConnector(
   const BoxLevel& head,
   const IntVector& min_connector_width,
   bool exact_width_only)
{
   if (s_always_create_missing_connector) {
      return findOrCreateConnector(head, min_connector_width, exact_width_only);
   }

   TBOX_ASSERT(d_my_mapped_box_level.isInitialized());
   TBOX_ASSERT(head.isInitialized());

   const Connector* found = NULL;
   for (int i = 0; i < d_cons_from_me.size(); ++i) {
      TBOX_ASSERT(d_cons_from_me[i]->isFinalized());
      TBOX_ASSERT(d_cons_from_me[i]->getBase().isInitialized());
      TBOX_ASSERT(d_cons_from_me[i]->getHead().isInitialized());
      TBOX_ASSERT(d_cons_from_me[i]->getBase().getBoxLevelHandle());
      TBOX_ASSERT(d_cons_from_me[i]->getHead().getBoxLevelHandle());
      TBOX_ASSERT(&d_cons_from_me[i]->getBase() ==
         &d_cons_from_me[i]->getBase().getBoxLevelHandle()->
         getBoxLevel());
      TBOX_ASSERT(&d_cons_from_me[i]->getHead() ==
         &d_cons_from_me[i]->getHead().getBoxLevelHandle()->
         getBoxLevel());

      if (&(d_cons_from_me[i]->getHead()) == &head) {
         if (d_cons_from_me[i]->getConnectorWidth() >= min_connector_width) {
            if (found == NULL) {
               found = d_cons_from_me[i];
            } else {
               IntVector vdiff =
                  d_cons_from_me[i]->getConnectorWidth()
                  - found->getConnectorWidth();

               TBOX_ASSERT(vdiff != IntVector::getZero(vdiff.getDim()));

               int diff = 0;
               for (int j = 0; j < vdiff.getDim().getValue(); ++j) {
                  diff += vdiff(j);
               }
               if (diff < 0) {
                  found = d_cons_from_me[i];
               }
            }
            if (found->getConnectorWidth() == min_connector_width) {
               break;
            }
         }
      }
   }

   OverlapConnectorAlgorithm oca;

   if (found == NULL) {

      TBOX_ERROR(
         "PersistentOverlapConnectors::findConnector: Failed to find Connector\n"
         << &d_my_mapped_box_level << "--->" << &head
         << " with " << (exact_width_only ? "exact" : "min")
         << " width of " << min_connector_width << ".\n"
         << "base:\n" << d_my_mapped_box_level.format("B: ")
         << "head:\n" << head.format("H: ")
         << "To automatically create the missing\n"
         << "connector, use findOrCreateConnector.");

   } else if (exact_width_only &&
              found->getConnectorWidth() != min_connector_width) {

      /*
       * Found a sufficient Connector, but it is too wide.  Extract
       * relevant neighbors from it to make a Connector with the exact
       * width.  This is scalable!
       */

      Connector* new_connector = new Connector(
         d_my_mapped_box_level,
         head,
         min_connector_width);
      oca.extractNeighbors(*new_connector, *found, min_connector_width);
      /*
       * Remove empty neighborhood sets.  They are not essential to an
       * overlap Connector.
       */
      new_connector->eraseEmptyNeighborSets();

      d_cons_from_me.push_back(new_connector);
      head.getPersistentOverlapConnectors().d_cons_to_me.push_back(
         new_connector);

      found = new_connector;

   }

   if (s_check_accessed_connectors == 'y') {
      if (oca.checkOverlapCorrectness(*found) != 0) {
         TBOX_ERROR("PersistentOverlapConnectors::findConnector errror:\n"
                    <<"Bad overlap Connector found.");
      }
   }

   return *found;
}

/*
 ************************************************************************
 *
 ************************************************************************
 */
const Connector&
PersistentOverlapConnectors::findOrCreateConnector(
   const BoxLevel& head,
   const IntVector& min_connector_width,
   bool exact_width_only)
{
   TBOX_ASSERT(d_my_mapped_box_level.isInitialized());
   TBOX_ASSERT(head.isInitialized());

   const Connector* found = NULL;
   for (int i = 0; i < d_cons_from_me.size(); ++i) {
      TBOX_ASSERT(d_cons_from_me[i]->isFinalized());
      TBOX_ASSERT(d_cons_from_me[i]->getBase().isInitialized());
      TBOX_ASSERT(d_cons_from_me[i]->getHead().isInitialized());
      TBOX_ASSERT(d_cons_from_me[i]->getBase().getBoxLevelHandle());
      TBOX_ASSERT(d_cons_from_me[i]->getHead().getBoxLevelHandle());
      TBOX_ASSERT(&d_cons_from_me[i]->getBase() ==
         &d_cons_from_me[i]->getBase().getBoxLevelHandle()->
         getBoxLevel());
      TBOX_ASSERT(&d_cons_from_me[i]->getHead() ==
         &d_cons_from_me[i]->getHead().getBoxLevelHandle()->
         getBoxLevel());

      if (&(d_cons_from_me[i]->getHead()) == &head) {
         if (d_cons_from_me[i]->getConnectorWidth() >= min_connector_width) {
            if (found == NULL) {
               found = d_cons_from_me[i];
            } else {
               IntVector vdiff =
                  d_cons_from_me[i]->getConnectorWidth()
                  - found->getConnectorWidth();

               TBOX_ASSERT(vdiff != IntVector::getZero(vdiff.getDim()));

               int diff = 0;
               for (int j = 0; j < vdiff.getDim().getValue(); ++j) {
                  diff += vdiff(j);
               }
               if (diff < 0) {
                  found = d_cons_from_me[i];
               }
            }
            if (found->getConnectorWidth() == min_connector_width) {
               break;
            }
         }
      }
   }

   OverlapConnectorAlgorithm oca;

   if (found == NULL) {

      Connector* new_connector = new Connector(
            d_my_mapped_box_level,
            head,
            min_connector_width,
            BoxLevel::DISTRIBUTED);
      oca.findOverlaps(*new_connector, head.getGlobalizedVersion());
      found = new_connector;

      d_cons_from_me.push_back(new_connector);
      head.getPersistentOverlapConnectors().d_cons_to_me.push_back(
         new_connector);

   } else if (exact_width_only &&
              found->getConnectorWidth() != min_connector_width) {

      /*
       * Found a sufficient Connector, but it is too wide.  Extract
       * relevant neighbors from it to make a Connector with the exact
       * width.  This is scalable!
       */

      Connector* new_connector = new Connector(
         d_my_mapped_box_level,
         head,
         min_connector_width);
      oca.extractNeighbors(*new_connector, *found, min_connector_width);
      /*
       * Remove empty neighborhood sets.  They are not essential to an
       * overlap Connector.
       */
      new_connector->eraseEmptyNeighborSets();

      d_cons_from_me.push_back(new_connector);
      head.getPersistentOverlapConnectors().d_cons_to_me.push_back(
         new_connector);

      found = new_connector;

   }

   if (s_check_accessed_connectors == 'y') {
      if (oca.checkOverlapCorrectness(*found) != 0) {
         TBOX_ERROR("PersistentOverlapConnectors::findOrCreateConnector errror:\n"
                    <<"Bad overlap Connector found.");
      }
   }

   return *found;
}

/*
 ************************************************************************
 *
 ************************************************************************
 */
bool
PersistentOverlapConnectors::hasConnector(
   const BoxLevel& head,
   const IntVector& min_connector_width,
   bool exact_width_only) const
{
   if (exact_width_only) {
      for (int i = 0; i < d_cons_from_me.size(); ++i) {
         if (&d_cons_from_me[i]->getHead() == &head &&
             d_cons_from_me[i]->getConnectorWidth() == min_connector_width) {
            return true;
         }
      }
   } else {
      for (int i = 0; i < d_cons_from_me.size(); ++i) {
         if (&d_cons_from_me[i]->getHead() == &head &&
             d_cons_from_me[i]->getConnectorWidth() >= min_connector_width) {
            return true;
         }
      }
   }
   return false;
}

/*
 ************************************************************************
 *
 ************************************************************************
 */
void
PersistentOverlapConnectors::clear()
{
   if (d_cons_from_me.empty() && d_cons_to_me.empty()) {
      return;
   }

   /*
    * Delete Connectors from me.
    */
   for (int i = 0; i < d_cons_from_me.size(); ++i) {

      const Connector* delete_me = d_cons_from_me[i];

      ConVect& cons_at_head =
         delete_me->getHead().getPersistentOverlapConnectors().d_cons_to_me;

      for (int j = 0; j < cons_at_head.size(); ++j) {
         if (cons_at_head[j] == delete_me) {
            cons_at_head.erase(j);
            break;
         }
      }

#ifdef DEBUG_CHECK_ASSERTIONS

      for (int j = 0; j < cons_at_head.size(); ++j) {
         TBOX_ASSERT(cons_at_head[j] != delete_me);
      }
#endif

      delete d_cons_from_me[i];
      d_cons_from_me[i] = NULL;
   }
   d_cons_from_me.clear();

   /*
    * Delete Connectors to me.
    */
   for (int i = 0; i < d_cons_to_me.size(); ++i) {

      const Connector* delete_me = d_cons_to_me[i];

      // Remove reference held by other end of Connector.
      ConVect& cons_at_base =
         delete_me->getBase().getPersistentOverlapConnectors().d_cons_from_me;

      for (int j = 0; j < cons_at_base.size(); ++j) {
         if (cons_at_base[j] == delete_me) {
            cons_at_base.erase(j);
            break;
         }
      }

#ifdef DEBUG_CHECK_ASSERTIONS

      for (int j = 0; j < cons_at_base.size(); ++j) {
         TBOX_ASSERT(cons_at_base[j] != delete_me);
      }

#endif

      delete d_cons_to_me[i];
      d_cons_to_me[i] = NULL;

   }
   d_cons_to_me.clear();
}

}
}
#endif
