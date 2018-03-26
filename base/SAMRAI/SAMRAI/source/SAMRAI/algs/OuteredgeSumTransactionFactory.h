/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Factory for creating outeredge sum transaction objects
 *
 ************************************************************************/

#ifndef included_algs_OuteredgeSumTransactionFactory
#define included_algs_OuteredgeSumTransactionFactory

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/hier/ComponentSelector.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/xfer/RefineClasses.h"
#include "SAMRAI/xfer/RefineTransactionFactory.h"

#include <boost/shared_ptr.hpp>

namespace SAMRAI {
namespace algs {

/*!
 * @brief Concrete subclass of the xfer::RefineTransactionFactory base class that
 * allocates outeredge sum transaction objects for a xfer::RefineSchedule object.
 *
 * @see xfer::RefineTransactionFactory
 * @see xfer::OuteredgeSumTransaction
 */

class OuteredgeSumTransactionFactory:public xfer::RefineTransactionFactory
{
public:
   /*!
    * @brief Default constructor.
    */
   OuteredgeSumTransactionFactory();

   /*!
    * @brief Virtual destructor for base class.
    */
   virtual ~OuteredgeSumTransactionFactory();

   /*!
    * @brief Set the array of xfer::RefineClass<DIM>::Data items used by the transactions.
    */
   void
   setRefineItems(
      const xfer::RefineClasses::Data ** refine_items,
      int num_refine_items);

   /*!
    * @brief Clear the array of xfer::RefineClass<DIM>::Data items used by the transactions.
    */
   void
   unsetRefineItems();

   /*!
    * @brief Allocate an OuteredgeSumTransaction object.
    *
    * @param dst_level      boost::shared_ptr to destination patch level.
    * @param src_level      boost::shared_ptr to source patch level.
    * @param overlap        boost::shared_ptr to overlap region between patches.
    * @param dst_node       Destination Box in destination patch level.
    * @param src_node       Source Box in source patch level.
    * @param ritem_id       Integer index of xfer::RefineClass<DIM>::Data item
    *                       associated with transaction.
    * @param box            Optional const reference to box defining region of
    *                       refine transaction.  Use next method if not required.
    * @param use_time_interpolation  Optional boolean flag indicating whether the
    *                       refine transaction involves time interpolation.
    *                       Default is false.
    */
   boost::shared_ptr<tbox::Transaction>
   allocate(
      const boost::shared_ptr<hier::PatchLevel>& dst_level,
      const boost::shared_ptr<hier::PatchLevel>& src_level,
      const boost::shared_ptr<hier::BoxOverlap>& overlap,
      const hier::Box& dst_node,
      const hier::Box& src_node,
      int ritem_id,
      const hier::Box& box,
      bool use_time_interpolation = false) const;

   /*!
    * @brief Allocate an OuteredgeSumTransaction object.
    *
    * Same as previous allocate routine but with default empty box and no
    * timer interpolation.
    */
   boost::shared_ptr<tbox::Transaction>
   allocate(
      const boost::shared_ptr<hier::PatchLevel>& dst_level,
      const boost::shared_ptr<hier::PatchLevel>& src_level,
      const boost::shared_ptr<hier::BoxOverlap>& overlap,
      const hier::Box& dst_node,
      const hier::Box& src_node,
      int ritem_id) const;

   /*!
    * @brief Function to initialize scratch space data for the sum transactions
    * (patch data components indicated by the component selector) to zero.
    *
    * @param level        boost::shared_ptr to patch level holding scratch data.
    * @param fill_time    Double value of simulation time at which preprocess
    *                     operation is called.
    * @param preprocess_vector Const reference to hier::ComponentSelector indicating
    *                     patch data array indices of scratch patch data objects
    *                     to preprocess.
    */
   void
   preprocessScratchSpace(
      const boost::shared_ptr<hier::PatchLevel>& level,
      double fill_time,
      const hier::ComponentSelector& preprocess_vector) const;

private:
   // The following two functions are not implemented
   OuteredgeSumTransactionFactory(
      const OuteredgeSumTransactionFactory&);
   void
   operator = (
      const OuteredgeSumTransactionFactory&);

   const xfer::RefineClasses::Data** d_refine_items;
   int d_number_refine_items;

};

}
}
#endif
