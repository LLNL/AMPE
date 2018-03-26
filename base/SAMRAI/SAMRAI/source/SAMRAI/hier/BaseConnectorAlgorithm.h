/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Algorithms to work with maping Connectors.
 *
 ************************************************************************/
#ifndef included_hier_BaseConnectorAlgorithm
#define included_hier_BaseConnectorAlgorithm

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxLevel.h"
#include "SAMRAI/hier/Connector.h"
#include "SAMRAI/tbox/AsyncCommPeer.h"

namespace SAMRAI {
namespace hier {

class BaseConnectorAlgorithm
{
protected:
   /*!
    * @brief Constructor
    */
   BaseConnectorAlgorithm();

   /*!
    * @brief Destructor.
    */
   virtual ~BaseConnectorAlgorithm();

   /*!
    * @brief Set up communication objects for use in privateBridge/Modify.
    */
   void
   setupCommunication(
      tbox::AsyncCommPeer<int> *& all_comms,
      tbox::AsyncCommStage& comm_stage,
      const tbox::SAMRAI_MPI& mpi,
      const std::set<int>& incoming_ranks,
      const std::set<int>& outgoing_ranks,
      const boost::shared_ptr<tbox::Timer>& mpi_wait_timer,
      int& operation_mpi_tag) const;

   //! @brief Send discovery to one processor during privateBridge/Modify.
   void
   sendDiscoveryToOneProcess(
      std::vector<int>& send_mesg,
      const int idx_offset_to_ref,
      BoxContainer& referenced_new_head_nabrs,
      BoxContainer& referenced_new_base_nabrs,
      tbox::AsyncCommPeer<int>& outgoing_comm,
      const tbox::Dimension& dim) const;

   /*!
    * @brief Receive messages and unpack info sent from other processes.
    */
   void
   receiveAndUnpack(
      Connector& new_base_to_new_head,
      Connector* new_head_to_new_base,
      std::set<int>& incoming_ranks,
      tbox::AsyncCommPeer<int> all_comms[],
      tbox::AsyncCommStage& comm_stage,
      const boost::shared_ptr<tbox::Timer>& receive_and_unpack_timer) const;

   //! @brief Unpack message sent by sendDiscoverytoOneProcess().
   void
   unpackDiscoveryMessage(
      const tbox::AsyncCommPeer<int>* incoming_comm,
      Connector& west_to_east,
      Connector* east_to_west) const;

   // Extra checks independent of optimization/debug.
   static char s_print_steps;

private:
   /*
    * Data length limit on first message of a communication.
    */
   static const int BASE_CONNECTOR_ALGORITHM_FIRST_DATA_LENGTH;
};

}
}

#endif // included_hier_BaseConnectorAlgorithm
