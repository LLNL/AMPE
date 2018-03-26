/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Manager of Connectors incident from a common BoxLevel.
 *
 ************************************************************************/
#ifndef included_hier_PersistentOverlapConnectors
#define included_hier_PersistentOverlapConnectors

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/hier/IntVector.h"

#include <vector>

namespace SAMRAI {
namespace hier {

class Connector;
class BoxLevel;

/*!
 * @brief A managager of overlap Connectors incident from a
 * BoxLevel, used to store and, if needed, generate overlap
 * Connectors in BoxLevel.
 *
 * PersistantOverlapConnectors provides a mechanism for objects using
 * the BoxLevel to look up overlap Connectors for the
 * BoxLevel.  Connectors created or returned by this class are
 * complete overlap Connectors that always contain the correct overlap
 * neighbor data.
 *
 * For improved scalability, Connectors can be constructed externally
 * and copied into the collection.  Connectors can also be
 * automatically computed using a non-scalable global search.
 *
 * If the input database contains a PersistentOverlapConnectors database,
 * the following inputs are recognized:
 *
 * <b>bool check_created_connectors:</b> When true, checks Connectors when
 * they are created.  The check is an non-scalable operation and is
 * meant for debugging.
 *
 * <b>bool check_accessed_connectors:</b> When true, check Connectors when
 * they are accessed.  The check is an non-scalable operation and is
 * meant for debugging.
 *
 * <b>bool always_create_missing_connector:</b> When true, override
 * findConnector() to behave like findOrCreateConnector().  This
 * essentially ensures that any Connectors sought are always found.
 *
 * @note
 * Creating Connectors this way (setting always_create_missing_connector
 * to true) is non-scalable. Nevertheless, the default is true, so that
 * application writers need not to worry about creating Connectors in
 * a scalable way.  For performance, this should be set to false.  To
 * selectively enable automatic Connector generation, set this to false and
 * use findOrCreateConnector() instead of findConnector() where one is
 * unsure if the Connector has been created.
 *
 * @see findConnector()
 * @see findOrCreateConnector()
 * @see hier::Connector
 */
class PersistentOverlapConnectors
{

public:
   /*!
    * @brief Deletes all Connectors to and from this object
    */
   ~PersistentOverlapConnectors();

   /*!
    * @brief Create an overlap Connector, computing relationships by
    * globalizing data.
    *
    * The base will be the BoxLevel that owns this object.
    * Find Connector relationships using a (non-scalable) global search.
    *
    * @see hier::Connector
    * @see hier::Connector::initialize()
    *
    * @param[in] head
    * @param[in] connector_width
    *
    * @return A const reference to the newly created overlap Connector.
    */
   const Connector&
   createConnector(
      const BoxLevel& head,
      const IntVector& connector_width);

   /*!
    * @brief Create an overlap Connector using externally
    * computed relationships.
    *
    * Create the Connector initialized with the arguments.
    * The base will be the BoxLevel that owns this object.
    *
    * @see hier::Connector
    * @see hier::Connector::initialize()
    *
    * @param[in] head
    * @param[in] connector_width
    * @param[in] relationships
    */
   const Connector&
   createConnector(
      const BoxLevel& head,
      const IntVector& connector_width,
      const Connector& relationships);

   /*!
    * @brief Cache the supplied overlap Connector.
    *
    * The head will be the supplied head and the base will be the
    * BoxLevel that owns this object.
    *
    * @param[in] head
    * @param[in] connector
    */
   void
   cacheConnector(
      const BoxLevel& head,
      Connector* connector);

   /*!
    * @brief Find an overlap Connector with the given head and minimum
    * Connector width.
    *
    * If multiple Connectors fit the criteria, the one with the
    * smallest ghost cell width (based on the algebraic sum of the
    * components) is selected.
    *
    * TODO: The criterion for selecting a single Connector is
    * arbitrary and should be re-examined.
    *
    * @par Assertions
    * If no Connector fits the criteria and @c
    * always_create_missing_connector is false, an assertion is
    * thrown.  To automatically create the Connector instead, use
    * findOrCreateConnector() or set @c
    * always_create_missing_connector to true.
    *
    * @param[in] head Find the overlap Connector with this specified head.
    * @param[in] min_connector_width Find the overlap Connector satisfying
    *      this minimum Connector width.
    * @param[in] exact_width_only If true, reject Connectors that do not
    *      match the requested width exactly.
    *
    * @return The Connector which matches the search criterion.
    */
   const Connector&
   findConnector(
      const BoxLevel& head,
      const IntVector& min_connector_width,
      bool exact_width_only = false);

   /*!
    * @brief Find or create an overlap Connectors with the
    * given head and minimum Connector width.
    *
    * If multiple Connectors fit the criteria, the one with the
    * smallest ghost cell width (based on the algebraic sum of the
    * components) is selected.
    *
    * TODO: The criterion for selecting a
    * single Connector is arbitrary and should be re-examined.
    *
    * If no Connector fits the criteria, a new one is created using
    * global search for edges.
    *
    * @param[in] head Find the overlap Connector with this specified head.
    * @param[in] min_connector_width Find the overlap Connector satisfying
    *      this minimum ghost cell width.
    * @param[in] exact_width_only If true, reject Connectors that do not
    *      match the requested width exactly.
    */
   const Connector&
   findOrCreateConnector(
      const BoxLevel& head,
      const IntVector& min_connector_width,
      bool exact_width_only = false);

   /*!
    * @brief Returns whether the object has overlap
    * Connectors with the given head and minimum Connector
    * width.
    *
    * TODO:  does the following comment mean that this must be called
    * before the call to findConnector?
    *
    * If this returns true, the Connector fitting the specification
    * exists and findConnector() will not throw an assertion.
    *
    * @param[in] head Find the overlap Connector with this specified head.
    * @param[in] min_connector_width Find the overlap Connector satisfying
    *      this minimum ghost cell width.
    * @param[in] exact_width_only If true, reject Connectors that do not
    *      match the requested width exactly.
    *
    * @return True if a Connector is found, otherwise false.
    */
   bool
   hasConnector(
      const BoxLevel& head,
      const IntVector& min_connector_width,
      bool exact_width_only = false) const;

   /*!
    * @brief Delete stored Connectors.
    */
   void
   clear();

private:
   //@{ @name Methods meant only for BoxLevel to use.

   /*!
    * @brief Constructor, to be called from the BoxLevel
    * allocating the object.
    *
    * @param my_mapped_box_level The BoxLevel served by this
    * object.
    */
   explicit PersistentOverlapConnectors(
      const BoxLevel& my_mapped_box_level);

   //@}

   //@{
   /*!
    * @brief Only BoxLevels are allowed to construct a
    * PersistentOverlapConnectors.
    */
   friend class BoxLevel;
   //@}

   typedef tbox::Array<const Connector *> ConVect;

   /*!
    * @brief Persistent overlap Connectors incident from me.
    */
   ConVect d_cons_from_me;

   /*!
    * @brief Persistent overlap Connectors incident to me.
    */
   ConVect d_cons_to_me;

   /*!
    * @brief Reference to the BoxLevel served by this object.
    */
   const BoxLevel& d_my_mapped_box_level;

   /*!
    * @brief Whether to check overlap Connectors when they are created.
    */
   static char s_check_created_connectors;

   /*!
    * @brief Whether to check overlap Connectors when they are accessed.
    */
   static char s_check_accessed_connectors;

   /*!
    * @brief Whether to force Connector finding functions to create
    * connectors that are missing.
    */
   static bool s_always_create_missing_connector;

};

}
}

#endif // included_hier_PersistentOverlapConnectors
