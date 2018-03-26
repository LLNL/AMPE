/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Singleton registry for all tranfer operators.
 *
 ************************************************************************/

#ifndef included_hier_TransferOperatorRegistry_C
#define included_hier_TransferOperatorRegistry_C

#include "SAMRAI/hier/TransferOperatorRegistry.h"
#include "SAMRAI/hier/BaseGridGeometry.h"

#include <typeinfo>

namespace SAMRAI {
namespace hier {

/*
 *************************************************************************
 *
 * Constructor and destructor for TransferOperatorRegistry objects.
 *
 *************************************************************************
 */

TransferOperatorRegistry::TransferOperatorRegistry(
   const tbox::Dimension& dim):
   d_min_stencil_width(dim, 0),
   d_dim(dim),
   d_max_op_stencil_width_req(false)
{
}

TransferOperatorRegistry::~TransferOperatorRegistry()
{
}

/*
 *************************************************************************
 *
 * Add operator to appropriate lookup hash map.
 *
 *************************************************************************
 */

void
TransferOperatorRegistry::addCoarsenOperator(
   const char* var_type_name,
   const boost::shared_ptr<CoarsenOperator>& coarsen_op)
{
   if (d_max_op_stencil_width_req &&
       (coarsen_op->getStencilWidth() > getMaxTransferOpStencilWidth())) {
      TBOX_WARNING(
         "Adding coarsen operator " << coarsen_op->getOperatorName()
                                    << "\nwith stencil width greater than current maximum\n"
                                    << "after call to getMaxTransferOpStencilWidth.\n");
   }
   boost::unordered_map<std::string, boost::unordered_map<std::string,
      boost::shared_ptr<CoarsenOperator> > >::iterator coarsen_ops =
         d_coarsen_operators.find(coarsen_op->getOperatorName());
   if (coarsen_ops == d_coarsen_operators.end()) {
      coarsen_ops = d_coarsen_operators.insert(
         std::make_pair(coarsen_op->getOperatorName(),
            boost::unordered_map<std::string,
               boost::shared_ptr<CoarsenOperator> >(0))).first;
   }
   coarsen_ops->second.insert(std::make_pair(var_type_name, coarsen_op));
}

void
TransferOperatorRegistry::addRefineOperator(
   const char* var_type_name,
   const boost::shared_ptr<RefineOperator>& refine_op)
{
   if (d_max_op_stencil_width_req &&
       (refine_op->getStencilWidth() > getMaxTransferOpStencilWidth())) {
      TBOX_WARNING(
         "Adding refine operator " << refine_op->getOperatorName()
                                   << "\nwith stencil width greater than current maximum\n"
                                   << "after call to getMaxTransferOpStencilWidth.\n");
   }
   boost::unordered_map<std::string, boost::unordered_map<std::string,
      boost::shared_ptr<RefineOperator> > >::iterator refine_ops =
         d_refine_operators.find(refine_op->getOperatorName());
   if (refine_ops == d_refine_operators.end()) {
      refine_ops = d_refine_operators.insert(
         std::make_pair(refine_op->getOperatorName(),
            boost::unordered_map<std::string,
               boost::shared_ptr<RefineOperator> >(0))).first;
   }
   refine_ops->second.insert(std::make_pair(var_type_name, refine_op));
}

void
TransferOperatorRegistry::addTimeInterpolateOperator(
   const char* var_type_name,
   const boost::shared_ptr<TimeInterpolateOperator>& time_op)
{
   boost::unordered_map<std::string, boost::unordered_map<std::string,
      boost::shared_ptr<TimeInterpolateOperator> > >::iterator time_ops =
         d_time_operators.find(time_op->getOperatorName());
   if (time_ops == d_time_operators.end()) {
      time_ops = d_time_operators.insert(
         std::make_pair(time_op->getOperatorName(),
            boost::unordered_map<std::string,
               boost::shared_ptr<TimeInterpolateOperator> >(0))).first;
   }
   time_ops->second.insert(std::make_pair(var_type_name, time_op));
}

/*
 *************************************************************************
 *
 * Search operator hash maps for operator matching request.
 *
 *************************************************************************
 */

boost::shared_ptr<CoarsenOperator>
TransferOperatorRegistry::lookupCoarsenOperator(
   BaseGridGeometry& grid_geometry,
   const boost::shared_ptr<Variable>& var,
   const std::string& op_name)
{
   TBOX_ASSERT(var);
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, *var);

   boost::shared_ptr<CoarsenOperator> coarsen_op;

   if ((op_name == "NO_COARSEN") ||
       (op_name == "USER_DEFINED_COARSEN") ||
       (op_name.empty())) {
   }
   else {
      if (d_coarsen_operators.empty()) {
         grid_geometry.buildOperators();
      }

      boost::unordered_map<std::string, boost::unordered_map<std::string,
         boost::shared_ptr<CoarsenOperator> > >::iterator coarsen_ops =
            d_coarsen_operators.find(op_name);
      if (coarsen_ops == d_coarsen_operators.end()) {
         TBOX_ERROR(
            "TransferOperatorRegistry::lookupCoarsenOperator" <<
            " could not find any operators with name " << op_name <<
            std::endl);
      }
      boost::unordered_map<std::string,
         boost::shared_ptr<CoarsenOperator> >::iterator the_op =
            coarsen_ops->second.find(typeid(*var).name());
      if (the_op == coarsen_ops->second.end()) {
         TBOX_ERROR(
            "TransferOperatorRegistry::lookupCoarsenOperator" <<
            " could not find operator with name " << op_name <<
            " for variable named " << typeid(*var).name() << std::endl);
      }
      coarsen_op = the_op->second;
   }

   return coarsen_op;
}

boost::shared_ptr<RefineOperator>
TransferOperatorRegistry::lookupRefineOperator(
   BaseGridGeometry& grid_geometry,
   const boost::shared_ptr<Variable>& var,
   const std::string& op_name)
{
   TBOX_ASSERT(var);
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, *var);

   boost::shared_ptr<RefineOperator> refine_op;

   if ((op_name == "NO_REFINE") ||
       (op_name == "USER_DEFINED_REFINE") ||
       (op_name.empty())) {
   }
   else {
      if (d_refine_operators.empty()) {
         grid_geometry.buildOperators();
      }

      boost::unordered_map<std::string, boost::unordered_map<std::string,
         boost::shared_ptr<RefineOperator> > >::iterator refine_ops =
            d_refine_operators.find(op_name);
      if (refine_ops == d_refine_operators.end()) {
         TBOX_ERROR(
            "TransferOperatorRegistry::lookupRefineOperator" <<
            " could not find any operators with name " << op_name <<
            std::endl);
      }
      boost::unordered_map<std::string,
         boost::shared_ptr<RefineOperator> >::iterator the_op =
            refine_ops->second.find(typeid(*var).name());
      if (the_op == refine_ops->second.end()) {
         TBOX_ERROR(
            "TransferOperatorRegistry::lookupRefineOperator" <<
            " could not find operator with name " << op_name <<
            " for variable named " << typeid(*var).name() << std::endl);
      }
      refine_op = the_op->second;
   }

   return refine_op;
}

boost::shared_ptr<TimeInterpolateOperator>
TransferOperatorRegistry::lookupTimeInterpolateOperator(
   BaseGridGeometry& grid_geometry,
   const boost::shared_ptr<Variable>& var,
   const std::string& op_name)
{
   TBOX_ASSERT(var);
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, *var);

   boost::shared_ptr<TimeInterpolateOperator> time_op;

   if ((op_name == "NO_TIME_INTERPOLATE") ||
       (op_name.empty())) {
   }
   else {
      if (d_time_operators.empty()) {
         grid_geometry.buildOperators();
      }

      boost::unordered_map<std::string, boost::unordered_map<std::string,
         boost::shared_ptr<TimeInterpolateOperator> > >::iterator time_ops =
            d_time_operators.find(op_name);
      if (time_ops == d_time_operators.end()) {
         TBOX_ERROR(
            "TransferOperatorRegistry::lookupTimeInterpolateOperator" <<
            " could not find any operators with name " << op_name <<
            std::endl);
      }
      boost::unordered_map<std::string,
         boost::shared_ptr<TimeInterpolateOperator> >::iterator the_op =
            time_ops->second.find(typeid(*var).name());
      if (the_op == time_ops->second.end()) {
         TBOX_ERROR(
            "TransferOperatorRegistry::lookupTimeInterpolateOperator" <<
            " could not find operator with name " << op_name <<
            " for variable named " << typeid(*var).name() << std::endl);
      }
      time_op = the_op->second;
   }

   return time_op;
}

/*
 *************************************************************************
 * Compute the max operator stencil width from all constructed
 * refine and coarsen operators and the user-specified minimum value.
 *************************************************************************
 */
IntVector
TransferOperatorRegistry::getMaxTransferOpStencilWidth()
{
   IntVector max_width(d_min_stencil_width);
   max_width.max(RefineOperator::getMaxRefineOpStencilWidth(getDim()));
   max_width.max(CoarsenOperator::getMaxCoarsenOpStencilWidth(getDim()));
   d_max_op_stencil_width_req = true;
   return max_width;
}

/*
 *************************************************************************
 *
 * Print TransferOperatorRegistry class data.
 *
 *************************************************************************
 */

void
TransferOperatorRegistry::printClassData(
   std::ostream& os) const
{
   os << "printing TransferOperatorRegistry data..." << std::endl;
   os << "TransferOperatorRegistry: this = "
      << (TransferOperatorRegistry *)this << std::endl;

   os << "Coarsen operators: " << std::endl;
   boost::unordered_map<std::string, boost::unordered_map<std::string,
      boost::shared_ptr<CoarsenOperator> > >::const_iterator cop =
         d_coarsen_operators.begin();
   while (cop != d_coarsen_operators.end()) {
      os << cop->first << std::endl;
      boost::unordered_map<std::string,
         boost::shared_ptr<CoarsenOperator> >::const_iterator ccop =
            cop->second.begin();
      while (ccop != cop->second.end()) {
         os << ccop->second.get() << std::endl;
         ++ccop;
      }
      ++cop;
   }

   os << "Refine operators: " << std::endl;
   boost::unordered_map<std::string, boost::unordered_map<std::string,
      boost::shared_ptr<RefineOperator> > >::const_iterator rop =
         d_refine_operators.begin();
   while (rop != d_refine_operators.end()) {
      os << rop->first << std::endl;
      boost::unordered_map<std::string,
         boost::shared_ptr<RefineOperator> >::const_iterator rrop =
            rop->second.begin();
      while (rrop != rop->second.end()) {
         os << rrop->second.get() << std::endl;
         ++rrop;
      }
      ++rop;
   }

   os << "Time interpolate operators: " << std::endl;
   boost::unordered_map<std::string, boost::unordered_map<std::string,
      boost::shared_ptr<TimeInterpolateOperator> > >::const_iterator top =
         d_time_operators.begin();
   while (top != d_time_operators.end()) {
      os << top->first << std::endl;
      boost::unordered_map<std::string,
         boost::shared_ptr<TimeInterpolateOperator> >::const_iterator ttop =
            top->second.begin();
      while (ttop != top->second.end()) {
         os << ttop->second.get() << std::endl;
         ++ttop;
      }
      ++top;
   }
}

}
}
#endif
