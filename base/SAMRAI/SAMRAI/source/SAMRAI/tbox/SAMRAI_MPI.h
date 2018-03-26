/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Simple utility class for interfacing with MPI
 *
 ************************************************************************/

#ifndef included_tbox_SAMRAI_MPI
#define included_tbox_SAMRAI_MPI

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/Complex.h"
#include "SAMRAI/tbox/Utilities.h"

#ifdef HAVE_MPI
#include "mpi.h"
#else

/*!
 * @brief Enumeration to define MPI constants when compiling without MPI.
 *
 * These are defined in the global namespace, and it does not matter
 * what values these take because they are not used.  (They are just place
 * holders to let code compile without MPI without adding excessive
 * preprocessor guards to the code.)
 *
 * This is not a complete set.  Developers should add as needed to extend
 * SAMRAI_MPI's functionality.
 */
enum {
   MPI_COMM_WORLD,
   MPI_COMM_NULL,
   // Special values:
   MPI_SUCCESS = 0,
   MPI_CONGRUENT,
   MPI_REQUEST_NULL,
   MPI_ERR_IN_STATUS,
   MPI_UNDEFINED,
   MPI_ANY_SOURCE,
   MPI_ANY_TAG,
   // Data types:
   MPI_BYTE,
   MPI_CHAR,
   MPI_DOUBLE,
   MPI_FLOAT,
   MPI_INT,
   MPI_DOUBLE_COMPLEX,
   MPI_2INT,
   MPI_DOUBLE_INT,
   MPI_FLOAT_INT,
   // Operations:
   MPI_MIN,
   MPI_MINLOC,
   MPI_MAX,
   MPI_MAXLOC,
   MPI_SUM,
   // Attributes:
   MPI_TAG_UB
};
#endif

namespace SAMRAI {
namespace tbox {

/*!
 * @brief Provides C++ wrapper around MPI routines.
 *
 * SAMRAI_MPI provides a single point of access in SAMRAI for using MPI
 * function calls and for the run-time decision whether to use MPI.
 * The purpose of having a single point is to facilitate writing code
 * that is not cluttered by the repetitive logic absorbed into this class.
 * Codes accessing MPI through this class should work whether SAMRAI is
 * configured with or without MPI and whether MPI is enabled at run time.
 *
 * This class provides two sets of interfaces:
 *
 * -# Methods closely matching the MPI interfaces associated with a
 * communicator.  The communicator is set during object construction
 * and is removed from the argument list of the MPI-like interfaces.
 *
 * -# Static methods exactly matching the MPI interfaces not associated
 * with any specific communicator.
 *
 * If SAMRAI has MPI enabled, the MPI wrappers in this class delegate to
 * the real MPI methods.  Otherwise, they are no-ops and will throw an
 * assertion if used, except for those specifically identified as able to
 * work without MPI enabled.  The purpose of these interfaces are strictly
 * to check whether MPI was configured into SAMRAI and SAMRAI is using
 * MPI at run time.  They are not "single-process implementations" of the
 * MPI standard.
 *
 * This class also manages a central SAMRAI_MPI object meant for
 * SAMRAI's general use.  See getSAMRAIWorld() and init().
 */

class SAMRAI_MPI
{

public:
   /*!
    * @brief Aliases for MPI data type names.
    *
    * Define aliases for MPI types that can be used whether SAMRAI is
    * configured with or without MPI.  Without MPI, these types are "dummies",
    * used to allow the code to be compiled without changes.
    */
#ifdef HAVE_MPI
   typedef MPI_Comm Comm;
   typedef MPI_Datatype Datatype;
   typedef MPI_Group Group;
   typedef MPI_Op Op;
   typedef MPI_Request Request;
   typedef MPI_Status Status;
#else
   typedef int Comm;
   typedef int Datatype;
   typedef int Group;
   typedef int Op;
   typedef int Request;

   /*!
    * @brief Dummy definition of Status to match the MPI standard.
    *
    * Codes are allowed to access the members of this struct, so they must
    * exist for compilation.  Without MPI, we won't actually use them.
    */
   struct Status {
      Status();
      int MPI_SOURCE;
      int MPI_TAG;
      int MPI_ERROR;
   };
#endif

   /*!
    * @brief MPI communicator constants
    */
   static Comm commWorld;
   static Comm commNull;

   /*!
    * @brief Get the primary SAMRAI_MPI object owned by SAMRAI.
    *
    * This is SAMRAI's primary communication object set up when
    * SAMRAI_MPI is initialized.  It is used for basic communications
    * that are not associated with another communication object.
    * Various parts of the SAMRAI library may create and use (and
    * destroy) communicators that are derived from this object.
    *
    * The use of this object outside of the SAMRAI library should be
    * carefully limited to avoid mixing messages.
    *
    * @see init()
    */
   static const SAMRAI_MPI&
   getSAMRAIWorld()
   {
      return s_samrai_world;
   }

   /*!
    * @brief Get a static invalid rank number.
    *
    * This value is intended to be used by other classes as an invalid rank
    * number rather than using a hard-coded "magic" negative integer value.
    */
   static int
   getInvalidRank()
   {
      return s_invalid_rank;
   }

   /*!
    * @brief Constructor.
    *
    * The given MPI communicator will be used for all communications invoked
    * by this object.
    *
    * Note that the object will NOT automatically free the communicator.
    * To manually free the communicator, use freeCommunicator().
    *
    * If the specified communicator is MPI_COMM_NULL, the rank (see getRank())
    * and size (see getSize()) are set to invalid values.  Otherwise, the
    * rank and size are set using the given communicator.  If MPI is enabled
    * but MPI has not been initialized, the rank and size are set to invalid
    * values.  If MPI is not enabled, the rank is set to 0 and size to 1.
    *
    * @param[in] comm
    */
   explicit SAMRAI_MPI(
      const Comm& comm);

   /*!
    * @brief Get the local process rank from the last time the
    * internal communicator was set.
    */
   int
   getRank() const
   {
      return d_rank;
   }

   /*!
    * @brief Get the size (number of processes) of the internal
    * communicator the last time it was set.
    */
   int
   getSize() const
   {
      return d_size;
   }

   /*!
    * @brief Get the internal communicator.
    */
   const Comm&
   getCommunicator() const
   {
      return d_comm;
   }

   /*!
    * @brief Set the internal communicator.
    *
    * Note that this call does does automatically free the
    * existing communicator.  To manually free communicators, see
    * freeCommunicator().
    *
    * @param[in] comm
    */
   void
   setCommunicator(
      const Comm& comm);

   /*!
    * @brief Duplicate and internally use the communicator of a given
    * SAMRAI_MPI object.
    *
    * Note that this call does not automatically free the
    * existing communicator.  The duplicate communicator will also NOT
    * be automatically freed.  To manually free communicators, see
    * freeCommunicator().
    *
    * @param[in] other  Contains the communicator to be duplicated.
    *
    */
   void
   dupCommunicator(
      const SAMRAI_MPI& other);

   /*!
    * @brief Free the internal communicator and set it to MPI_COMM_NULL.
    *
    * If the internal communicator is already MPI_COMM_NULL, do nothing.
    */
   void
   freeCommunicator();

   /*!
    * @brief Assignment operator.
    *
    * @param[in] rhs
    */
   const SAMRAI_MPI&
   operator = (
      const SAMRAI_MPI& rhs)
   {
      d_comm = rhs.d_comm;
      d_rank = rhs.d_rank;
      d_size = rhs.d_size;
      return *this;
   }

   /*!
    * @brief Equality comparison operator (compares MPI communicator).
    *
    * @param[in] rhs
    */
   bool
   operator == (
      const SAMRAI_MPI& rhs) const
   {
      return d_comm == rhs.d_comm;
   }

   /*!
    * @brief Inequality comparison operator (compares MPI communicator).
    *
    * @param[in] rhs
    */
   bool
   operator != (
      const SAMRAI_MPI& rhs) const
   {
      return d_comm != rhs.d_comm;
   }

   //@{
   //! @name Static MPI wrappers matching MPI interfaces.

   /*!
    * @brief MPI wrappers that don't use the internal communicator.
    *
    * The purpose of these wrappers is to provide a single place for
    * compile- and run-time toggling of MPI code.  The signatures and
    * return values of these methods are identical to the MPI C bindings.
    * Thes methods will throw an assertion if they are called when MPI is
    * not enabled.
    */

   static int
   Comm_compare(
      Comm comm1,
      Comm comm2,
      int* result);

   static int
   Comm_free(
      Comm* comm);

   static int
   Finalized(
      int* flag);

   static int
   Get_count(
      Status* status,
      Datatype datatype,
      int* count);

   static int
   Test(
      Request* request,
      int* flag,
      Status* status);

   static int
   Test_cancelled(
      Status* status,
      int* flag);

   static int
   Wait(
      Request* request,
      Status* status);

   static int
   Waitall(
      int count,
      Request* reqs,
      Status* stats);

   static int
   Waitany(
      int count,
      Request* array_of_requests,
      int* index,
      Status* status);

   static int
   Waitsome(
      int incount,
      Request* array_of_requests,
      int* outcount,
      int* array_of_indices,
      Status* array_of_statuses);

   /*!
    * @brief MPI Wtime (if MPI is enabled) or an alternate time
    * calculation.
    *
    * @return If MPI is running, use Wtime.  If not, but POSIX time is
    * available, use that.  Else return 0.0.
    */
   static double
   Wtime();

   //@}

   //@{
   //! @name MPI wrappers bound to SAMRAI_MPI objects.

   /*!
    * @brief MPI wrappers for methods associated with an MPI communicator.
    *
    * The purpose of these wrappers is to provide a single place for compile-
    * and run-time toggling of MPI code.  The signatures and return values
    * of these methods are identical to the MPI C bindings, except that the
    * communicators are omitted.  The communicator used is that associated
    * with the SAMRAI_MPI object (typically passed to the constructor).
    * These methods throw an assertion if called while MPI is not enabled,
    * except where noted.
    */

   int
   Allgather(
      void* sendbuf,
      int sendcount,
      Datatype sendtype,
      void* recvbuf,
      int recvcount,
      Datatype recvtype) const;

   int
   Allgatherv(
      void* sendbuf,
      int sendcvount,
      Datatype sendtype,
      void* recvbuf,
      int* recvcounts,
      int* displs,
      Datatype recvtype) const;

   int
   Allreduce(
      void* sendbuf,
      void* recvbuf,
      int count,
      Datatype datatype,
      Op op) const;

   int
   Attr_get(
      int keyval,
      void* attribute_val,
      int* flag) const;

   /*!
    * @brief MPI Barrier (does nothing when MPI is disabled).
    */
   int
   Barrier() const;

   /*!
    * @brief MPI Bcast (does nothing when MPI is disabled).
    */
   int
   Bcast(
      void* buffer,
      int count,
      Datatype datatype,
      int root) const;

   int
   Comm_dup(
      Comm* newcomm) const;

   /*!
    * @brief MPI Comm_rank (Set rank to 0 when MPI is disabled).
    */
   int
   Comm_rank(
      int* rank) const;

   /*!
    * @brief MPI Comm_size (Set size to 1 when MPI is disabled).
    */
   int
   Comm_size(
      int* size) const;

   int
   Gather(
      void* sendbuf,
      int sendcount,
      Datatype sendtype,
      void* recvbuf,
      int recvcount,
      Datatype recvtype,
      int root) const;

   int
   Gatherv(
      void* sendbuf,
      int sendcount,
      Datatype sendtype,
      void* recvbuf,
      int* recvcounts,
      int* displs,
      Datatype recvtype,
      int root) const;

   int
   Iprobe(
      int source,
      int tag,
      int* flag,
      Status* status) const;

   int
   Isend(
      void* buf,
      int count,
      Datatype datatype,
      int dest,
      int tag,
      Request* req) const;

   int
   Irecv(
      void* buf,
      int count,
      Datatype datatype,
      int source,
      int tag,
      Request* request) const;

   int
   Probe(
      int source,
      int tag,
      Status* status) const;

   int
   Recv(
      void* buf,
      int count,
      Datatype datatype,
      int source,
      int tag,
      Status* status) const;

   int
   Reduce(
      void* sendbuf,
      void* recvbuf,
      int count,
      Datatype datatype,
      Op op,
      int root) const;

   int
   Send(
      void* buf,
      int count,
      Datatype datatype,
      int dest,
      int tag) const;

   //@}

   //@{

   //! @name Specialized reductions for certain common types.

   /*!
    * @brief Specialized Allreduce for integers.
    *
    * If MPI_MINLOC or MPI_MAXLOC is given as the operator, the
    * ranks_of_extrema argument must be provided with space allocated for
    * @c count integers.
    *
    * @return MPI error code
    *
    * @param[in,out] x   Array of integers to reduce.
    * @param[in]  count  Number of items in x.
    * @param[in]  op     A valid MPI reduce operation.
    * @param[out] ranks_of_extrema  Ranks associated with min or max of x
    *                               (if op indicates min or max operation).
    */
   int
   AllReduce(
      int* x,
      int count,
      Op op,
      int* ranks_of_extrema = NULL) const;

   /*!
    * @brief Specialized Allreduce for doubles.
    *
    * If MPI_MINLOC or MPI_MAXLOC is given as the operator, the
    * ranks_of_extrema argument must be provided with space allocated for
    * @c count integers.
    *
    * @param[in,out] x   Array of doubles to reduce.
    * @param[in]  count  Number of items in x.
    * @param[in]  op     A valid MPI reduce operation.
    * @param[out] ranks_of_extrema  Ranks associated with min or max of x
    *                               (if op indicates min or max operation).
    */
   int
   AllReduce(
      double* x,
      int count,
      Op op,
      int* ranks_of_extrema = NULL) const;

   /*!
    * @brief Specialized Allreduce for floats.
    *
    * If MPI_MINLOC or MPI_MAXLOC is given as the operator, the
    * ranks_of_extrema argument must be provided with space allocated for
    * @c count integers.
    *
    * @param[in,out] x   Array of floats to reduce.
    * @param[in]  count  Number of items in x.
    * @param[in]  op     A valid MPI reduce operation.
    * @param[out] ranks_of_extrema  Ranks associated with min or max of x
    *                               (if op indicates min or max operation).
    */
   int
   AllReduce(
      float* x,
      int count,
      Op op,
      int* ranks_of_extrema = NULL) const;

   //@}

   //@{

   //! @name Generic high-level operations not in MPI interfaces.

   /*!
    * @brief Parallel prefix sum for integers.
    *
    * Given an input x, the output is the sum of x from process 0 up
    * to and including the local process's x.  This implementation
    * allows an array of x, each of which is summed independently of
    * the other.
    *
    * @param[in,out] x   Array of integers to sum.
    *
    * @param[in] count Number of items in x.  Must be the same on all
    * processes.
    *
    * @param[in] tag MPI tag for communication.
    *
    * @return MPI_SUCCESS or an MPI error code.
    */
   int
   parallelPrefixSum(
      int* x,
      int count,
      int tag) const;

   // @}

   /*!
    * @brief Set flag indicating whether exit or MPI_Abort is called
    * when running with one processor.
    *
    * Calling this function influences the behavior of calls to
    * SAMRAI_MPI::abort().  If the value set is true, it means that
    * system abort() will be called.  Passing false means exit(-1) will be
    * called.
    */
   static void
   setCallAbortInSerialInsteadOfExit(
      bool flag = true)
   {
      s_call_abort_in_serial_instead_of_exit = flag;
   }

   /*!
    * @brief Set flag indicating whether MPI_Abort or abort is called
    * when running with more than one processor.
    *
    * Calling this function influences the behavior of calls to
    * SAMRAI_MPI::abort().  If the value set is true, it means that
    * system abort() will be called.  Passing false means MPI_Abort will be
    * called.
    */
   static void
   setCallAbortInParallelInsteadOfMPIAbort(
      bool flag = true)
   {
      s_call_abort_in_parallel_instead_of_mpiabort = flag;
   }

   /*!
    * @brief Call MPI_Abort or exit depending on whether running with one
    * or more processes and value set by function above, if called.
    *
    * The default is to call exit(-1) if running with one processor and to
    * call MPI_Abort() otherwise.  This function avoids having to guard abort
    * calls in application code.
    */
   static void
   abort();

   /*!
    * @brief Disable MPI usage and run in sequential mode.
    *
    * The default setting is to enable MPI at run-time if SAMRAI was
    * configured with MPI and to disable if SAMRAI was not configured
    * with MPI.  This method disables MPI in the former case.
    *
    * If MPI is disabled for SAMRAI, the communication interfaces in
    * this class are no-ops and will throw an assertion if used,
    * except for those interfaces specifically described as working in
    * sequential mode.
    *
    * Disabling MPI must be done before calling SAMRAI_MPI::init().
    */
   static void
   disableMPI();

   /*!
    * @brief Whether SAMRAI is using MPI.
    *
    * @see disableMPI().
    */
   static bool
   usingMPI()
   {
      return s_mpi_is_initialized;
   }

   /*!
    * @brief Initialize MPI and SAMRAI_MPI.
    *
    * This initialization sets up getSAMRAIWorld() to return an
    * SAMRAI_MPI object with a communicator duplicated from
    * MPI_COMM_WORLD.
    *
    * Use only one of the three initialization methods.
    *
    * @param[in]  argc  Pointer to parameter from main()
    * @param[in]  argv  Pointer to parameter from main()
    */
   static void init(
      int* argc,
      char** argv[]);

   /*!
    * @brief Initialize SAMRAI_MPI when MPI is already initialized.
    *
    * This initialization sets up getSAMRAIWorld() to return an
    * SAMRAI_MPI object with a communicator duplicated from
    * the one passed in.
    *
    * Use only one of the three initialization methods.
    *
    * @param comm  An MPI Communicator to be duplicated for SAMRAI use.
    */
   static void
   init(
      Comm comm);

   /*!
    * @brief Initialize SAMRAI_MPI to work without MPI.
    *
    * This initialization sets up SAMRAI to run in sequential mode.
    *
    * Use only one of the three initialization methods.
    */
   static void
   initMPIDisabled();

   /*!
    * @brief Shut down SAMRAI_MPI and, if appropriate, shut down MPI.
    *
    * This method frees the MPI communicator duplicated for SAMRAI
    * use.  If this class has started MPI with the first init() method,
    * finalize() will also shut down MPI.
    */
   static void
   finalize();

private:
   //@{

   /*!
    * @brief Flags to control program abort in SAMRAI_MPI::abort().
    */
   static bool s_call_abort_in_serial_instead_of_exit;
   static bool s_call_abort_in_parallel_instead_of_mpiabort;
   //@}

   /*!
    * @brief Invalid (negative) rank number for getInvalidRank().
    */
   static int s_invalid_rank;

   //@{
   //@name Structs for passing arguments to MPI
   struct DoubleIntStruct { double d;
                            int i;
   };
   struct FloatIntStruct { float f;
                           int i;
   };
   struct IntIntStruct { int j;
                         int i;
   };
   //@}

   /*!
    * @brief Internal communicator
    */
   SAMRAI_MPI::Comm d_comm;

   /*!
    * @brief The local rank.
    */
   int d_rank;

   /*!
    * @brief The size of the communicator.
    */
   int d_size;

   /*!
    * @brief Whether MPI is initialized.
    */
   static bool s_mpi_is_initialized;

   /*!
    * @brief Whether this class started up MPI.
    */
   static bool s_we_started_mpi;

   /*!
    * @brief Primary SAMRAI_MPI object.
    */
   static SAMRAI_MPI s_samrai_world;

};

}
}

#endif
