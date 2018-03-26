/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Interface for factory objects that create transactions for
 *                oarsen schedules.
 *
 ************************************************************************/

#ifndef included_xfer_CoarsenTransactionFactory
#define included_xfer_CoarsenTransactionFactory

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/Transaction.h"
#include "SAMRAI/hier/BaseGridGeometry.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/xfer/CoarsenClasses.h"

namespace SAMRAI {
namespace xfer {

/*!
 * @brief Abstract base class defining the interface for all concrete transaction
 * factory objects that generate data transaction objects used with a CoarsenSchedule
 * object.  A concrete subclass will allocate new transaction objects.  This class
 * is an example of the ``Abstract Factory'' method described in the Design Patterns
 * book by Gamma, et al.
 *
 * To add a new type of Transaction object MyCoarsenTransaction:
 *
 * -# Implement a concrete CoarsenTransactionFactory object as a subclass
 *       that is derived from this CoarsenTransactionFactory base class.
 *       Implement the abstract virtual functions as appropriate for the
 *       concrete subclass; in particular, the allocate() function must return
 *       a new instance of the desired transaction object.
 * -# The type of the transaction allocated by the concrete factory is a
 *       Transaction<DIM>.  Thus, the new transaction object must be derived
 *       from the Transaction<DIM> base class and implement the abstract
 *       virtual functions declared by the base class.
 *
 * @see tbox::Transaction
 */

class CoarsenTransactionFactory
{
public:
   /*!
    * @brief Default constructor.
    */
   CoarsenTransactionFactory();

   /*!
    * @brief Virtual destructor.
    */
   virtual ~CoarsenTransactionFactory();

   /*!
    * @brief Pure virtual function to set the array of CoarsenClass::Data items
    * associated with the coarsen schedule.  Typical concrete transactions used by
    * the schedule use this information to communicate data.  This operation
    * is called by the coarsen schedule during the execution of the
    * CoarsenSchedule::fillData() routine before data communication
    * operations begin.
    */
   virtual void
   setCoarsenItems(
      const CoarsenClasses::Data ** coarsen_items,
      int num_coarsen_items) = 0;

   /*!
    * @brief Pure virtual function to clear the array of CoarsenClass::Data items
    * associated with the coarsen schedule.  This operation is called by the
    * coarsen schedule after data communication operations are complete.
    */
   virtual void
   unsetCoarsenItems() = 0;

   /*!
    * @brief Pure virtual function to allocate a concrete coarsen transaction object.
    * This routine is called by the coarsen schedule during construction of the
    * schedule.
    *
    * @param dst_level      boost::shared_ptr to destination patch level.
    * @param src_level      boost::shared_ptr to source patch level.
    * @param overlap        boost::shared_ptr to overlap region between patches.
    * @param dst_mapped_box Destination Box in destination patch level.
    * @param src_mapped_box Source Box in source patch level.
    * @param citem_id       Integer index of CoarsenClass::Data item associated
    *                       with transaction.
    */
   virtual boost::shared_ptr<tbox::Transaction>
   allocate(
      const boost::shared_ptr<hier::PatchLevel>& dst_level,
      const boost::shared_ptr<hier::PatchLevel>& src_level,
      const boost::shared_ptr<hier::BoxOverlap>& overlap,
      const hier::Box& dst_mapped_box,
      const hier::Box& src_mapped_box,
      int citem_id) const = 0;

private:
   CoarsenTransactionFactory(
      const CoarsenTransactionFactory&);                        // not implemented
   void
   operator = (
      const CoarsenTransactionFactory&);                  // not implemented

};

}
}
#endif
