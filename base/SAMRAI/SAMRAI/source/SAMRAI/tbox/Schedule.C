/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Schedule of communication transactions between processors
 *
 ************************************************************************/
#include "SAMRAI/tbox/Schedule.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/TimerManager.h"

#include <cstring>

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace tbox {

typedef std::list<boost::shared_ptr<Transaction> >::iterator Iterator;
typedef std::list<boost::shared_ptr<Transaction> >::const_iterator ConstIterator;

const int Schedule::s_default_first_tag = 0;
const int Schedule::s_default_second_tag = 1;
/*
 * TODO: Set the default first message length to the maximum value
 * possible without incurring any additional cost associated with the
 * MPI communication.  This parameter should be dependent on the MPI
 * implementation.
 */
const size_t Schedule::s_default_first_message_length = 1000;

const std::string Schedule::s_default_timer_prefix("tbox::Schedule");
std::map<std::string, Schedule::TimerStruct> Schedule::s_static_timers;

StartupShutdownManager::Handler
Schedule::s_initialize_finalize_handler(
   Schedule::initializeCallback,
   0,
   0,
   Schedule::finalizeCallback,
   StartupShutdownManager::priorityTimers);

/*
 *************************************************************************
 *************************************************************************
 */

Schedule::Schedule():
   d_coms(NULL),
   d_com_stage(),
   d_mpi(SAMRAI_MPI::getSAMRAIWorld()),
   d_first_tag(s_default_first_tag),
   d_second_tag(s_default_second_tag),
   d_first_message_length(s_default_first_message_length)
{
   setTimerPrefix("tbox::Schedule");
}

/*
 *************************************************************************
 * Note that the destructor should not be called during a communication
 * phase.
 *************************************************************************
 */
Schedule::~Schedule()
{
   if (d_coms != NULL) {
      TBOX_ERROR("Destructing a schedule while communication is pending\n"
         << "leads to lost messages.  Aborting.");
   }
}

/*
 *************************************************************************
 * Add a transaction to the head of a list of data transactions in
 * this schedule. The assignment of the transaction to a list depends
 * on the source and destination processors of the transaction.
 *************************************************************************
 */
void
Schedule::addTransaction(
   const boost::shared_ptr<Transaction>& transaction)
{
   const int src_id = transaction->getSourceProcessor();
   const int dst_id = transaction->getDestinationProcessor();

   if ((d_mpi.getRank() == src_id) && (d_mpi.getRank() == dst_id)) {
      d_local_set.push_front(transaction);
   } else {
      if (d_mpi.getRank() == dst_id) {
         d_recv_sets[src_id].push_front(transaction);
      } else if (d_mpi.getRank() == src_id) {
         d_send_sets[dst_id].push_front(transaction);
      }
   }
}

/*
 *************************************************************************
 * Append a transaction to the tail of a list of data transactions in
 * this schedule.  The assignment of the transaction to a list depends
 * on the source and destination processors of the transaction.
 *************************************************************************
 */
void
Schedule::appendTransaction(
   const boost::shared_ptr<Transaction>& transaction)
{
   const int src_id = transaction->getSourceProcessor();
   const int dst_id = transaction->getDestinationProcessor();

   if ((d_mpi.getRank() == src_id) && (d_mpi.getRank() == dst_id)) {
      d_local_set.push_back(transaction);
   } else {
      if (d_mpi.getRank() == dst_id) {
         d_recv_sets[src_id].push_back(transaction);
      } else if (d_mpi.getRank() == src_id) {
         d_send_sets[dst_id].push_back(transaction);
      }
   }
}

/*
 *************************************************************************
 * Access number of send transactions.
 *************************************************************************
 */
int
Schedule::getNumSendTransactions(
   const int rank) const
{
   int size = 0;
   TransactionSets::const_iterator mi = d_send_sets.find(rank);
   if (mi != d_send_sets.end()) {
      size = static_cast<int>(mi->second.size());
   }
   return size;
}

/*
 *************************************************************************
 * Access number of receive transactions.
 *************************************************************************
 */
int
Schedule::getNumRecvTransactions(
   const int rank) const
{
   int size = 0;
   TransactionSets::const_iterator mi = d_recv_sets.find(rank);
   if (mi != d_recv_sets.end()) {
      size = static_cast<int>(mi->second.size());
   }
   return size;
}

/*
 *************************************************************************
 * Perform the communication described by the schedule.
 *************************************************************************
 */
void
Schedule::communicate()
{
   d_object_timers->t_communicate->start();
   beginCommunication();
   finalizeCommunication();
   d_object_timers->t_communicate->stop();
}

/*
 *************************************************************************
 * Begin communication but do not wait for it to finish.  This routine
 * posts receives, and sends outgoing messages.  Since we do not wait
 * for message completion, use finalizeCommunication() to ensure that
 * communication has finished.
 *************************************************************************
 */
void
Schedule::beginCommunication()
{
   d_object_timers->t_begin_communication->start();
   allocateCommunicationObjects();
   postReceives();
   postSends();
   d_object_timers->t_begin_communication->stop();
}

/*
 *************************************************************************
 * Perform the local data copies, complete receive operations and
 * unpack received data into their destinations.
 *************************************************************************
 */
void
Schedule::finalizeCommunication()
{
   d_object_timers->t_finalize_communication->start();
   performLocalCopies();
   processCompletedCommunications();
   deallocateCommunicationObjects();
   d_object_timers->t_finalize_communication->stop();
}

/*
 *************************************************************************
 * Post receives.
 *
 * Where message lengths can be locally computed, use the correct
 * message lengths to avoid overheads due to unknown lengths.
 *************************************************************************
 */
void
Schedule::postReceives()
{
   if (d_recv_sets.empty()) {
      /*
       * Short cut because some looping logic in this method assumes
       * non-empty d_recv_sets.
       */
      return;
   }

   int rank = d_mpi.getRank();

   /*
    * We loop through d_recv_sets starting with the highest rank that
    * is lower than local process.  We loop backwards, continuing at
    * the opposite end when we run out of sets.  This ordering is in
    * the reverse direction of the message send ordering so that a
    * send posted earlier is paired with a receive that is also posted
    * earlier.
    */
   AsyncCommPeer<char>* recv_coms = d_coms;

   // Initialize iterators to where we want to start looping.
   size_t icom = 0; // Index into recv_coms.
   while (icom < d_recv_sets.size() &&
          recv_coms[icom].getPeerRank() < rank) {
      ++icom;
   }
   icom = icom > 0 ? icom - 1 : d_recv_sets.size() - 1;

   // Map iterator mi corresponds to recv_coms[icom].
   TransactionSets::const_iterator mi =
      d_recv_sets.find(recv_coms[icom].getPeerRank());

   for (size_t counter = 0;
        counter < d_recv_sets.size();
        ++counter, --mi, --icom) {

      TBOX_ASSERT(mi->first == recv_coms[icom].getPeerRank());

      // Compute incoming message size, if possible.
      const std::list<boost::shared_ptr<Transaction> >& transactions =
         mi->second;
      unsigned int byte_count = 0;
      bool can_estimate_incoming_message_size = true;
      for (ConstIterator r = transactions.begin();
           r != transactions.end(); r++) {
         if (!(*r)->canEstimateIncomingMessageSize()) {
            can_estimate_incoming_message_size = false;
            break;
         }
         byte_count +=
            static_cast<unsigned int>((*r)->computeIncomingMessageSize());
      }

      // Set AsyncCommPeer to receive known message length.
      if (can_estimate_incoming_message_size) {
         recv_coms[icom].limitFirstDataLength(byte_count);
      }

      // Begin non-blocking receive operation.
      d_object_timers->t_post_receives->start();
      recv_coms[icom].beginRecv();
      if (recv_coms[icom].isDone()) {
         recv_coms[icom].pushToCompletionQueue();
      }
      d_object_timers->t_post_receives->stop();

      if (mi == d_recv_sets.begin()) {
         // Continue loop at the opposite end.
         mi = d_recv_sets.end();
         icom = d_recv_sets.size();
      }
   }
}

/*
 *************************************************************************
 * Allocate the send buffer, pack the data, and initiate the message
 * sends.
 *************************************************************************
 */
void
Schedule::postSends()
{
   d_object_timers->t_post_sends->start();
   /*
    * We loop through d_send_sets starting with the first set with
    * rank higher than the local process, continuing at the opposite
    * end when we run out of sets.  This ordering tends to spread out
    * the communication traffic over the entire network to reduce the
    * potential network contention.
    */

   int rank = d_mpi.getRank();

   AsyncCommPeer<char>* send_coms = d_coms + d_recv_sets.size();

   // Initialize iterators to where we want to start looping.
   TransactionSets::const_iterator mi = d_send_sets.upper_bound(rank);
   size_t icom = 0; // send_coms[icom] corresponds to mi.
   while (icom < d_send_sets.size() &&
          send_coms[icom].getPeerRank() < rank) {
      ++icom;
   }

   for (size_t counter = 0;
        counter < d_send_sets.size();
        ++counter, ++mi, ++icom) {

      if (mi == d_send_sets.end()) {
         // Continue loop at the opposite end.
         mi = d_send_sets.begin();
         icom = 0;
      }
      TBOX_ASSERT(mi->first == send_coms[icom].getPeerRank());

      // Compute message size and whether receiver can estimate it.
      const std::list<boost::shared_ptr<Transaction> >& transactions =
         mi->second;
      size_t byte_count = 0;
      bool can_estimate_incoming_message_size = true;
      for (ConstIterator pack = transactions.begin();
           pack != transactions.end(); pack++) {
         if (!(*pack)->canEstimateIncomingMessageSize()) {
            can_estimate_incoming_message_size = false;
         }
         byte_count += (*pack)->computeOutgoingMessageSize();
      }

      // Pack outgoing data into a message.
      MessageStream outgoing_stream(byte_count, MessageStream::Write);
      d_object_timers->t_pack_stream->start();
      for (ConstIterator pack = transactions.begin();
           pack != transactions.end(); pack++) {
         (*pack)->packStream(outgoing_stream);
      }
      d_object_timers->t_pack_stream->stop();

      if (can_estimate_incoming_message_size) {
         // Receiver knows message size, so set it exactly.
         send_coms[icom].limitFirstDataLength(byte_count);
      }

      // Begin non-blocking send operation.
      send_coms[icom].beginSend(
         (const char *)outgoing_stream.getBufferStart(),
         static_cast<int>(outgoing_stream.getCurrentSize()));
      if (send_coms[icom].isDone()) {
         send_coms[icom].pushToCompletionQueue();
      }
   }

   d_object_timers->t_post_sends->stop();
}

/*
 *************************************************************************
 * Perform all of the local memory-to-memory copies for this processor.
 *************************************************************************
 */
void
Schedule::performLocalCopies()
{
   d_object_timers->t_local_copies->start();
   for (Iterator local = d_local_set.begin();
        local != d_local_set.end(); local++) {
      (*local)->copyLocalData();
   }
   d_object_timers->t_local_copies->stop();
}

/*
 *************************************************************************
 * Process completed operations as they come in.  Initially, completed
 * operations are placed in d_completed_comm.  Process these first,
 * then check for next set of completed operations.  Repeat until all
 * operations are completed.
 *************************************************************************
 */
void
Schedule::processCompletedCommunications()
{
   d_object_timers->t_process_incoming_messages->start();

   while ( d_com_stage.numberOfCompletedMembers() > 0 ||
           d_com_stage.advanceSome() ) {

      AsyncCommPeer<char>* completed_comm =
         dynamic_cast<AsyncCommPeer<char> *>(d_com_stage.popCompletionQueue());

      TBOX_ASSERT(completed_comm != NULL);
      TBOX_ASSERT(completed_comm->isDone());
      if (static_cast<size_t>(completed_comm - d_coms) < d_recv_sets.size()) {

         const int sender = completed_comm->getPeerRank();

         // Copy message into stream.
         MessageStream incoming_stream(
            completed_comm->getRecvSize() * sizeof(char),
            MessageStream::Read,
            completed_comm->getRecvData());
         completed_comm->clearRecvData();

         d_object_timers->t_unpack_stream->start();
         for (Iterator recv = d_recv_sets[sender].begin();
              recv != d_recv_sets[sender].end(); recv++) {
            (*recv)->unpackStream(incoming_stream);
         }
         d_object_timers->t_unpack_stream->stop();
      } else {
         // No further action required for completed send.
      }
   }

   d_object_timers->t_process_incoming_messages->stop();
}

/*
 *************************************************************************
 * Allocate communication objects, set them up on the stage and get
 * them ready to send/receive.
 *************************************************************************
 */
void
Schedule::allocateCommunicationObjects()
{
   const size_t length = d_recv_sets.size() + d_send_sets.size();
   d_coms = new AsyncCommPeer<char>[length];

   size_t counter = 0;
   for (TransactionSets::iterator ti = d_recv_sets.begin();
        ti != d_recv_sets.end();
        ++ti) {
      d_coms[counter].initialize(&d_com_stage);
      d_coms[counter].setPeerRank(ti->first);
      d_coms[counter].setMPITag(d_first_tag, d_second_tag);
      d_coms[counter].setMPI(d_mpi);
      d_coms[counter].limitFirstDataLength(d_first_message_length);
      ++counter;
   }
   for (TransactionSets::iterator ti = d_send_sets.begin();
        ti != d_send_sets.end();
        ++ti) {
      d_coms[counter].initialize(&d_com_stage);
      d_coms[counter].setPeerRank(ti->first);
      d_coms[counter].setMPITag(d_first_tag, d_second_tag);
      d_coms[counter].setMPI(d_mpi);
      d_coms[counter].limitFirstDataLength(d_first_message_length);
      ++counter;
   }
}

/*
 *************************************************************************
 * Print class data to the specified output stream.
 *************************************************************************
 */
void
Schedule::printClassData(
   std::ostream& stream) const
{
   stream << "Schedule::printClassData()" << std::endl;
   stream << "-------------------------------" << std::endl;

   stream << "Number of sends: " << d_send_sets.size() << std::endl;
   stream << "Number of recvs: " << d_recv_sets.size() << std::endl;

   for (TransactionSets::const_iterator ss = d_send_sets.begin();
        ss != d_send_sets.end(); ++ss) {
      const std::list<boost::shared_ptr<Transaction> >& send_set = ss->second;
      stream << "Send Set: " << ss->first << std::endl;
      for (ConstIterator send = send_set.begin();
           send != send_set.end(); send++) {
         (*send)->printClassData(stream);
      }
   }

   for (TransactionSets::const_iterator rs = d_recv_sets.begin();
        rs != d_recv_sets.end(); ++rs) {
      const std::list<boost::shared_ptr<Transaction> >& recv_set = rs->second;
      stream << "Recv Set: " << rs->first << std::endl;
      for (ConstIterator recv = recv_set.begin();
           recv != recv_set.end(); recv++) {
         (*recv)->printClassData(stream);
      }
   }

   stream << "Local Set" << std::endl;
   for (ConstIterator local = d_local_set.begin();
        local != d_local_set.end(); local++) {
      (*local)->printClassData(stream);
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
Schedule::setTimerPrefix(
   const std::string& timer_prefix)
{
   std::map<std::string, TimerStruct>::iterator ti(
      s_static_timers.find(timer_prefix));
   if (ti == s_static_timers.end()) {
      d_object_timers = &s_static_timers[timer_prefix];
      getAllTimers(timer_prefix, *d_object_timers);
   } else {
      d_object_timers = &(ti->second);
   }
   d_com_stage.setCommunicationWaitTimer(d_object_timers->t_MPI_wait);
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
Schedule::getAllTimers(
   const std::string& timer_prefix,
   TimerStruct& timers)
{
   timers.t_communicate = TimerManager::getManager()->
      getTimer(timer_prefix + "::communicate()");
   timers.t_begin_communication = TimerManager::getManager()->
      getTimer(timer_prefix + "::beginCommunication()");
   timers.t_finalize_communication = TimerManager::getManager()->
      getTimer(timer_prefix + "::finalizeCommunication()");
   timers.t_post_receives = TimerManager::getManager()->
      getTimer(timer_prefix + "::postReceives()");
   timers.t_post_sends = TimerManager::getManager()->
      getTimer(timer_prefix + "::postSends()");
   timers.t_process_incoming_messages = TimerManager::getManager()->
      getTimer(timer_prefix + "::processIncomingMessages()");
   timers.t_MPI_wait = TimerManager::getManager()->
      getTimer(timer_prefix + "::MPI_wait");
   timers.t_pack_stream = TimerManager::getManager()->
      getTimer(timer_prefix + "::pack_stream");
   timers.t_unpack_stream = TimerManager::getManager()->
      getTimer(timer_prefix + "::unpack_stream");
   timers.t_local_copies = TimerManager::getManager()->
      getTimer(timer_prefix + "::performLocalCopies()");
}

}
}

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Unsuppress XLC warnings
 */
#pragma report(enable, CPPC5334)
#pragma report(enable, CPPC5328)
#endif
