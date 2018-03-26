/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Singleton manager class for statistic objects.
 *
 ************************************************************************/

#include "SAMRAI/tbox/Statistician.h"

#include "SAMRAI/tbox/IOStream.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/RestartManager.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"
#include "SAMRAI/tbox/Schedule.h"
#include "SAMRAI/tbox/StatTransaction.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/MathUtilities.h"

#include <boost/make_shared.hpp>

namespace SAMRAI {
namespace tbox {

const int Statistician::DEFAULT_NUMBER_OF_TIMERS_INCREMENT = 128;
const int StatisticRestartDatabase::TBOX_STATISTICRESTARTDATABASE_VERSION = 1;

Statistician * Statistician::s_statistician_instance = (Statistician *)NULL;

StartupShutdownManager::Handler
Statistician::s_finalize_handler(
   Statistician::initializeCallback,
   0,
   Statistician::shutdownCallback,
   Statistician::finalizeCallback,
   StartupShutdownManager::priorityStatistician);

/*
 *************************************************************************
 *
 * Static statistician member functions.
 *
 *************************************************************************
 */

Statistician*
Statistician::createStatistician(
   bool register_for_restart,
   bool read_from_restart)
{
   makeStatisticianInstance(register_for_restart, read_from_restart);
   return s_statistician_instance;
}

Statistician*
Statistician::getStatistician()
{
   /* Should have instance constructed in initializeCallback */
   TBOX_ASSERT(s_statistician_instance);

   return s_statistician_instance;
}

void
Statistician::initializeCallback()
{
   if (!s_statistician_instance) {
      makeStatisticianInstance();
   }
}

void
Statistician::shutdownCallback()
{
   if (s_statistician_instance) {
      if (s_statistician_instance->d_restart_database_instance) {
         delete s_statistician_instance->d_restart_database_instance;
      }
      s_statistician_instance->d_restart_database_instance = 0;
   }
}

void
Statistician::finalizeCallback()
{
   if (s_statistician_instance) {
      delete s_statistician_instance;
      s_statistician_instance = 0;
   }
}

void
Statistician::registerSingletonSubclassInstance(
   Statistician* subclass_instance)
{
   if (!s_statistician_instance) {
      s_statistician_instance = subclass_instance;
   } else {
      TBOX_ERROR("Statistician internal error...\n"
         << "Attemptng to set Singleton instance to subclass instance,"
         << "\n but Singleton instance already set." << std::endl);
   }
}

/*
 *************************************************************************
 *
 * Protected statistician constructor and destructor.
 *
 *************************************************************************
 */

Statistician::Statistician():
   d_has_gathered_stats(false)
{
   d_restart_database_instance = (StatisticRestartDatabase *)NULL;

   d_must_call_finalize = true;

   setMaximumNumberOfStatistics(DEFAULT_NUMBER_OF_TIMERS_INCREMENT);

   d_num_proc_stats = 0;
   d_num_patch_stats = 0;

}

Statistician::~Statistician()
{
   if (d_restart_database_instance) delete d_restart_database_instance;

   d_proc_statistics.resizeArray(0);
   d_patch_statistics.resizeArray(0);

   d_num_proc_stats = 0;
   d_num_patch_stats = 0;
}

/*
 *************************************************************************
 *
 * Private members for creating and managing the singleton instance.
 *
 *************************************************************************
 */

void
Statistician::makeStatisticianInstance(
   bool register_for_restart,
   bool read_from_restart)
{
   /* If reading from restart then force new instance
    * to be created from existing state in the
    * restart file */
   if (read_from_restart) {
      delete s_statistician_instance;
      s_statistician_instance = 0;
   }

   if (!s_statistician_instance) {
      s_statistician_instance = new Statistician();
      s_statistician_instance->initRestartDatabase(register_for_restart,
         read_from_restart);
   }
}

void
Statistician::initRestartDatabase(
   bool register_for_restart,
   bool read_from_restart)
{
   d_restart_database_instance =
      new StatisticRestartDatabase("StatisticRestartDatabase",
         register_for_restart,
         read_from_restart);
}

/*
 *************************************************************************
 *
 * Utility functions for getting statistics, adding them to the
 * database, checking whether a particular statistic exists, resetting.
 *
 *************************************************************************
 */

boost::shared_ptr<Statistic>
Statistician::getStatistic(
   const std::string& name,
   const std::string& stat_type)
{
   TBOX_ASSERT(!name.empty());
   TBOX_ASSERT(!stat_type.empty());

   boost::shared_ptr<Statistic> stat;

   bool found = false;

   if (stat_type == "PROC_STAT") {

      for (int i = 0; i < d_num_proc_stats; i++) {
         if (d_proc_statistics[i]->getName() == name) {
            stat = d_proc_statistics[i];
            found = true;
            break;
         }
      }

      if (!found) {
         if (d_num_proc_stats ==
             Statistician::getMaximumNumberOfStatistics()) {
            setMaximumNumberOfStatistics(
               Statistician::getMaximumNumberOfStatistics()
               + DEFAULT_NUMBER_OF_TIMERS_INCREMENT);
         }
         stat.reset(new Statistic(name, stat_type, d_num_proc_stats));
         d_proc_statistics[d_num_proc_stats] = stat;
         d_num_proc_stats++;
         d_must_call_finalize = true;
      }

   } else if (stat_type == "PATCH_STAT") {

      for (int i = 0; i < d_num_patch_stats; i++) {
         if (d_patch_statistics[i]->getName() == name) {
            stat = d_patch_statistics[i];
            found = true;
            break;
         }
      }

      if (!found) {
         if (d_num_patch_stats ==
             Statistician::getMaximumNumberOfStatistics()) {
            setMaximumNumberOfStatistics(
               Statistician::getMaximumNumberOfStatistics()
               + DEFAULT_NUMBER_OF_TIMERS_INCREMENT);
         }
         stat.reset(new Statistic(name, stat_type, d_num_patch_stats));
         d_patch_statistics[d_num_patch_stats] = stat;
         d_num_patch_stats++;
         d_must_call_finalize = true;
      }

   } else {
      TBOX_ERROR("Statistician::getStatistic error ..."
         << "\n   Unrecognized stat type string " << stat_type
         << "passed to routine." << std::endl);
   }

   return stat;

}

bool
Statistician::checkStatisticExists(
   boost::shared_ptr<Statistic>& stat,
   const std::string& name) const
{
   TBOX_ASSERT(!name.empty());

   if (checkProcStatExists(stat, name)) {
      return true;
   }

   if (checkPatchStatExists(stat, name)) {
      return true;
   }

   return false;

}

bool
Statistician::checkProcStatExists(
   boost::shared_ptr<Statistic>& stat,
   const std::string& name) const
{
   stat.reset();

   bool stat_found = false;
   for (int i = 0; i < d_num_proc_stats; i++) {
      if (d_proc_statistics[i]->getName() == name) {
         stat_found = true;
         stat = d_proc_statistics[i];
         break;
      }
   }

   return stat_found;
}

bool
Statistician::checkPatchStatExists(
   boost::shared_ptr<Statistic>& stat,
   const std::string& name) const
{
   stat.reset();

   bool stat_found = false;
   for (int i = 0; i < d_num_patch_stats; i++) {
      if (d_patch_statistics[i]->getName() == name) {
         stat_found = true;
         stat = d_patch_statistics[i];
         break;
      }
   }

   return stat_found;
}

void
Statistician::resetProcessorStatistics()
{
   for (int i = 0; i < d_num_proc_stats; i++) {
      d_proc_statistics[i]->reset();
   }
   d_must_call_finalize = true;
}

void
Statistician::resetPatchStatistics()
{
   for (int i = 0; i < d_num_patch_stats; i++) {
      d_patch_statistics[i]->reset();
   }
   d_must_call_finalize = true;
}

void
Statistician::resetStatistics()
{
   resetProcessorStatistics();
   resetPatchStatistics();
}

/*
 *************************************************************************
 *
 * Utility functions to retrieve global statistic data from
 * statistician database.
 *
 *************************************************************************
 */

int
Statistician::getProcStatId(
   const std::string& name) const
{
   int ret_val = -1;

   if (!name.empty()) {
      boost::shared_ptr<Statistic> stat;
      if (checkProcStatExists(stat, name)) {
         ret_val = stat->getInstanceId();
      }
   }

   return ret_val;

}

int
Statistician::getGlobalProcStatSequenceLength(
   int proc_stat_id)
{
   const SAMRAI_MPI& mpi(SAMRAI_MPI::getSAMRAIWorld());
   int seq_len = -1;

   if (mpi.getRank() == 0) {

      if (d_must_call_finalize) {
         TBOX_ERROR("Statistician::getGlobalProcStatSeqLength ..."
            << "\n   The finalize() method to construct global data "
            "must be called BEFORE this method." << std::endl);

      }
      if (!d_has_gathered_stats) {
         TBOX_ERROR("Statistician::getGlobalProcStatSeqLength ..."
            << "\n   The finalize() method to construct global data "
            << "must be called with the argument to gather global "
            << "stats data BEFORE this metho." << std::endl);
      }

      TBOX_ASSERT(proc_stat_id >= 0 &&
         proc_stat_id < d_global_proc_stat_data.getSize());

      seq_len = d_global_proc_stat_data[proc_stat_id].getSize();
   }

   return seq_len;

}

double
Statistician::getGlobalProcStatValue(
   int proc_stat_id,
   int seq_num,
   int proc_num)
{
   const SAMRAI_MPI& mpi(SAMRAI_MPI::getSAMRAIWorld());
   double val = 0.;

   if (mpi.getRank() == 0) {

      if (d_must_call_finalize) {
         TBOX_ERROR("Statistician::getGlobalProcStatValue ..."
            << "\n   The finalize() method to construct global data "
            "must be called BEFORE this method." << std::endl);

      }
      if (!d_has_gathered_stats) {
         TBOX_ERROR("Statistician::getGlobalProcStatSeqLength ..."
            << "\n   The finalize() method to construct global data "
            << "must be called with the argument to gather global "
            << "stats data BEFORE this metho." << std::endl);
      }

      TBOX_ASSERT(proc_stat_id >= 0 &&
         proc_stat_id < d_global_proc_stat_data.getSize());
      TBOX_ASSERT(seq_num >= 0 &&
         seq_num < d_global_proc_stat_data[proc_stat_id].getSize());
      TBOX_ASSERT(proc_num < mpi.getSize());

      val = d_global_proc_stat_data[proc_stat_id][seq_num][proc_num];
   }

   return val;
}

double
Statistician::getGlobalProcStatSum(
   int proc_stat_id,
   int seq_num)
{
   const SAMRAI_MPI& mpi(SAMRAI_MPI::getSAMRAIWorld());
   double sum = 0.;

   if (mpi.getRank() == 0) {

      if (d_must_call_finalize) {
         TBOX_ERROR("Statistician::getGlobalProcStatSum ..."
            << "\n   The finalize() method to construct global data "
            "must be called BEFORE this method." << std::endl);

      }

      TBOX_ASSERT(proc_stat_id >= 0 &&
         proc_stat_id < d_global_proc_stat_sum.getSize());
      TBOX_ASSERT(seq_num >= 0 &&
         seq_num < d_global_proc_stat_sum[proc_stat_id].getSize());

      sum = d_global_proc_stat_sum[proc_stat_id][seq_num];
   }

   return sum;
}

double
Statistician::getGlobalProcStatMax(
   int proc_stat_id,
   int seq_num)
{
   const SAMRAI_MPI& mpi(SAMRAI_MPI::getSAMRAIWorld());
   double pmax = -(MathUtilities<double>::getMax());

   if (mpi.getRank() == 0) {

      if (d_must_call_finalize) {
         TBOX_ERROR("Statistician::getGlobalProcStatMax ..."
            << "\n   The finalize() method to construct global data "
            "must be called BEFORE this method." << std::endl);

      }

      TBOX_ASSERT(proc_stat_id >= 0 &&
         proc_stat_id < d_global_proc_stat_max.getSize());
      TBOX_ASSERT(seq_num >= 0 &&
         seq_num < d_global_proc_stat_max[proc_stat_id].getSize());

      pmax = d_global_proc_stat_max[proc_stat_id][seq_num];
   }

   return pmax;
}

int
Statistician::getGlobalProcStatMaxProcessorId(
   int proc_stat_id,
   int seq_num)
{
   const SAMRAI_MPI& mpi(SAMRAI_MPI::getSAMRAIWorld());
   int id = -1;

   if (mpi.getRank() == 0) {

      if (d_must_call_finalize) {
         TBOX_ERROR("Statistician::getGlobalProcStatMaxProcId ..."
            << "\n   The finalize() method to construct global data "
            "must be called BEFORE this method." << std::endl);

      }

      TBOX_ASSERT(proc_stat_id >= 0 &&
         proc_stat_id < d_global_proc_stat_data.getSize());
      TBOX_ASSERT(seq_num >= 0 &&
         seq_num < d_global_proc_stat_imax[proc_stat_id].getSize());

      id = d_global_proc_stat_imax[proc_stat_id][seq_num];

   }

   return id;
}

double
Statistician::getGlobalProcStatMin(
   int proc_stat_id,
   int seq_num)
{
   const SAMRAI_MPI& mpi(SAMRAI_MPI::getSAMRAIWorld());
   double pmin = MathUtilities<double>::getMax();

   if (mpi.getRank() == 0) {

      if (d_must_call_finalize) {
         TBOX_ERROR("Statistician::getGlobalProcStatMin ..."
            << "\n   The finalize() method to construct global data "
            "must be called BEFORE this method." << std::endl);

      }

      TBOX_ASSERT(proc_stat_id >= 0 &&
         proc_stat_id < d_global_proc_stat_min.getSize());
      TBOX_ASSERT(seq_num >= 0 &&
         seq_num < d_global_proc_stat_min[proc_stat_id].getSize());

      pmin = d_global_proc_stat_min[proc_stat_id][seq_num];
   }

   return pmin;
}

int
Statistician::getGlobalProcStatMinProcessorId(
   int proc_stat_id,
   int seq_num)
{
   const SAMRAI_MPI& mpi(SAMRAI_MPI::getSAMRAIWorld());
   int id = -1;

   if (mpi.getRank() == 0) {

      if (d_must_call_finalize) {
         TBOX_ERROR("Statistician::getGlobalProcStatMinProcId ..."
            << "\n   The finalize() method to construct global data "
            "must be called BEFORE this method." << std::endl);

      }

      TBOX_ASSERT(proc_stat_id >= 0 &&
         proc_stat_id < d_global_proc_stat_imin.getSize());
      TBOX_ASSERT(seq_num >= 0 &&
         seq_num < d_global_proc_stat_imin[proc_stat_id].getSize());

      id = d_global_proc_stat_imin[proc_stat_id][seq_num];
   }

   return id;
}

void
Statistician::printGlobalProcStatData(
   int proc_stat_id,
   std::ostream& os,
   int precision)
{
   const SAMRAI_MPI& mpi(SAMRAI_MPI::getSAMRAIWorld());

   if (mpi.getRank() == 0) {

      if (d_must_call_finalize) {
         TBOX_ERROR("Statistician::printGlobalProcStatData ..."
            << "\n   The finalize() method to construct global data "
            "must be called BEFORE this method." << std::endl);

      }
      if (!d_has_gathered_stats) {
         TBOX_ERROR("Statistician::getGlobalProcStatSeqLength ..."
            << "\n   The finalize() method to construct global data "
            << "must be called with the argument to gather global "
            << "stats data BEFORE this metho." << std::endl);
      }

      TBOX_ASSERT(proc_stat_id >= 0);
      TBOX_ASSERT(precision > 0);

      os.precision(precision);

      os << "\n   " << proc_stat_id << ":    "
         << d_proc_statistics[proc_stat_id]->getName() << std::endl;

      int nnodes = mpi.getSize();
      const Array<Array<double> >& sdata =
         d_global_proc_stat_data[proc_stat_id];

      for (int ipsl = 0; ipsl < sdata.getSize(); ipsl++) {
         os << "      Seq # " << ipsl << std::endl;
         os << "         proc : value" << std::endl;
         for (int ip = 0; ip < nnodes; ip++) {
            /*
             * Write out data only if data entry is NOT an "empty"
             * entry, defined by the Statistic::s_empty_seq_tag_entry
             * value.
             */
            if (!(MathUtilities<double>::equalEps(sdata[ipsl][ip],
                     Statistic::s_empty_seq_tag_entry))) {
               os << "         " << ip << "    : " << sdata[ipsl][ip]
                  << std::endl;
            }
         }
      }
   }

}

void
Statistician::printGlobalProcStatDataFormatted(
   int proc_stat_id,
   std::ostream& os,
   int precision)
{
   const SAMRAI_MPI& mpi(SAMRAI_MPI::getSAMRAIWorld());

   if (mpi.getRank() == 0) {

      if (d_must_call_finalize) {
         TBOX_ERROR("Statistician::printGlobalPatchStatDataFormatted"
            << "\n   The finalize() method to construct global data "
            "must be called BEFORE this method." << std::endl);

      }
      if (!d_has_gathered_stats) {
         TBOX_ERROR("Statistician::getGlobalProcStatSeqLength ..."
            << "\n   The finalize() method to construct global data "
            << "must be called with the argument to gather global "
            << "stats data BEFORE this metho." << std::endl);
      }

      TBOX_ASSERT(proc_stat_id >= 0);
      TBOX_ASSERT(precision > 0);

      os.precision(precision);

      /*
       * Output in C++ is by default right justified with the setw()
       * option e.g. cout << "[" << setw(5) << 1 << "]" will output
       * [   1].  Using setf(ios::left) makes it left justified, which
       * is more convenient to output columns of tables.
       */
      os.setf(std::ios::left);
      int s, n;
      int nnodes = mpi.getSize();

      // heading - line 1
      os << "Seq#\t";
      for (n = 0; n < nnodes; n++) {
         os << "Proc\t";
      }
      os << std::endl;

      // heading - line 2
      os << "    \t";
      for (n = 0; n < nnodes; n++) {
         os << n << "\t";
      }
      os << std::endl;

      /*
       * Now print values.
       */
      const Array<Array<double> >& sdata =
         d_global_proc_stat_data[proc_stat_id];
      for (s = 0; s < sdata.getSize(); s++) {
         os << s << "\t";
         for (n = 0; n < nnodes; n++) {
            /*
             * Write out data only if data entry is NOT an "empty"
             * entry, defined by the Statistic::s_empty_seq_tag_entry
             * value.
             */
            if (MathUtilities<double>::equalEps(sdata[s][n],
                   Statistic::s_empty_seq_tag_entry)) {
               os << "  " << "\t";
            } else {
               os << sdata[s][n] << "\t";
            }

         }
         os << std::endl;
      }

   }

}

void
Statistician::printGlobalProcStatDataFormatted(
   int proc_stat_id,
   int proc_id,
   std::ostream& os,
   int precision)
{
   const SAMRAI_MPI& mpi(SAMRAI_MPI::getSAMRAIWorld());

   if (mpi.getRank() == 0) {

      if (d_must_call_finalize) {
         TBOX_ERROR("Statistician::printGlobalPatchStatDataFormatted"
            << "\n   The finalize() method to construct global data "
            "must be called BEFORE this method." << std::endl);

      }
      if (!d_has_gathered_stats) {
         TBOX_ERROR("Statistician::getGlobalProcStatSeqLength ..."
            << "\n   The finalize() method to construct global data "
            << "must be called with the argument to gather global "
            << "stats data BEFORE this metho." << std::endl);
      }

      TBOX_ASSERT(proc_stat_id >= 0);
      TBOX_ASSERT(proc_id >= 0);
      TBOX_ASSERT(precision > 0);

      os.precision(precision);

      /*
       * Output in C++ is by default right justified with the setw()
       * option e.g. cout << "[" << setw(5) << 1 << "]" will output
       * [   1].  Using setf(ios::left) makes it left justified, which
       * is more convenient to output columns of tables.
       */
      os.setf(std::ios::left);
      int s;

      // heading - line 1
      os << "Seq#\t" << "Proc\t" << std::endl;

      // heading - line 2
      os << "    \t" << proc_id << "\t" << std::endl;

      /*
       * Now print values.
       */
      const Array<Array<double> >& sdata =
         d_global_proc_stat_data[proc_stat_id];
      for (s = 0; s < sdata.getSize(); s++) {
         os << s << "\t";
         /*
          * Write out data only if data entry is NOT an "empty"
          * entry, defined by the Statistic::s_empty_seq_tag_entry
          * value.
          */
         if (MathUtilities<double>::equalEps(sdata[s][proc_id],
                Statistic::s_empty_seq_tag_entry)) {
            os << "  " << "\t";
         } else {
            os << sdata[s][proc_id] << "\t";
         }
         os << std::endl;
      }

   }

}

int
Statistician::getPatchStatId(
   const std::string& name) const
{
   int ret_val = -1;

   if (!name.empty()) {
      boost::shared_ptr<Statistic> stat;
      if (checkPatchStatExists(stat, name)) {
         ret_val = stat->getInstanceId();
      }
   }

   return ret_val;

}

int
Statistician::getGlobalPatchStatSequenceLength(
   int patch_stat_id)
{
   const SAMRAI_MPI& mpi(SAMRAI_MPI::getSAMRAIWorld());
   int seq_len = -1;

   if (mpi.getRank() == 0) {

      if (d_must_call_finalize) {
         TBOX_ERROR("Statistician::getGlobalProcStatSeqLen ..."
            << "\n   The finalize() method to construct global data "
            "must be called BEFORE this method." << std::endl);

      }
      if (!d_has_gathered_stats) {
         TBOX_ERROR("Statistician::getGlobalProcStatSeqLength ..."
            << "\n   The finalize() method to construct global data "
            << "must be called with the argument to gather global "
            << "stats data BEFORE this metho." << std::endl);
      }

      TBOX_ASSERT(patch_stat_id >= 0 &&
         patch_stat_id < d_global_patch_stat_data.getSize());

      seq_len = d_global_patch_stat_data[patch_stat_id].getSize();
   }

   return seq_len;

}

int
Statistician::getGlobalPatchStatNumberPatches(
   int patch_stat_id,
   int seq_num)
{
   const SAMRAI_MPI& mpi(SAMRAI_MPI::getSAMRAIWorld());

   int num_patches = -1;

   if (mpi.getRank() == 0) {

      if (d_must_call_finalize) {
         TBOX_ERROR("Statistician::getGlobalProcStatNumPatches ..."
            << "\n   The finalize() method to construct global data "
            "must be called BEFORE this method." << std::endl);

      }
      if (!d_has_gathered_stats) {
         TBOX_ERROR("Statistician::getGlobalProcStatSeqLength ..."
            << "\n   The finalize() method to construct global data "
            << "must be called with the argument to gather global "
            << "stats data BEFORE this metho." << std::endl);
      }

      TBOX_ASSERT(patch_stat_id >= 0 &&
         patch_stat_id < d_global_patch_stat_data.getSize());
      TBOX_ASSERT(seq_num >= 0 &&
         seq_num < d_global_patch_stat_data[patch_stat_id].getSize());

      num_patches =
         d_global_patch_stat_data[patch_stat_id][seq_num].getSize();
   }

   return num_patches;

}

int
Statistician::getGlobalPatchStatPatchMapping(
   int patch_stat_id,
   int seq_num,
   int patch_num)
{
   const SAMRAI_MPI& mpi(SAMRAI_MPI::getSAMRAIWorld());
   int mapping = -1;

   if (mpi.getRank() == 0) {

      if (d_must_call_finalize) {
         TBOX_ERROR("Statistician::getGlobalPatchStatPatchMapping ..."
            << "\n   The finalize() method to construct global data "
            "must be called BEFORE this method." << std::endl);

      }
      if (!d_has_gathered_stats) {
         TBOX_ERROR("Statistician::getGlobalProcStatSeqLength ..."
            << "\n   The finalize() method to construct global data "
            << "must be called with the argument to gather global "
            << "stats data BEFORE this metho." << std::endl);
      }

      TBOX_ASSERT(patch_stat_id >= 0 &&
         patch_stat_id < d_global_patch_stat_mapping.getSize());
      TBOX_ASSERT(seq_num >= 0 &&
         seq_num < d_global_patch_stat_mapping[patch_stat_id].getSize());
      TBOX_ASSERT(patch_num >= 0 &&
         patch_num < d_global_patch_stat_mapping[patch_stat_id][seq_num].
         getSize());

      mapping =
         d_global_patch_stat_mapping[patch_stat_id][seq_num][patch_num];
   }

   return mapping;
}

double
Statistician::getGlobalPatchStatValue(
   int patch_stat_id,
   int seq_num,
   int patch_num)
{
   const SAMRAI_MPI& mpi(SAMRAI_MPI::getSAMRAIWorld());
   double val = 0.;

   if (mpi.getRank() == 0) {

      if (d_must_call_finalize) {
         TBOX_ERROR("Statistician::getGlobalPatchStatValue ..."
            << "\n   The finalize() method to construct global data "
            "must be called BEFORE this method." << std::endl);

      }
      if (!d_has_gathered_stats) {
         TBOX_ERROR("Statistician::getGlobalProcStatSeqLength ..."
            << "\n   The finalize() method to construct global data "
            << "must be called with the argument to gather global "
            << "stats data BEFORE this metho." << std::endl);
      }

      TBOX_ASSERT(patch_stat_id >= 0 &&
         patch_stat_id < d_global_patch_stat_data.getSize());
      TBOX_ASSERT(seq_num >= 0 &&
         seq_num < d_global_patch_stat_data[patch_stat_id].getSize());
      TBOX_ASSERT(patch_num >= 0 &&
         patch_num < d_global_patch_stat_data[patch_stat_id][seq_num].
         getSize());
      val = d_global_patch_stat_data[patch_stat_id][seq_num][patch_num];
   }

   return val;
}

double
Statistician::getGlobalPatchStatSum(
   int patch_stat_id,
   int seq_num)
{
   const SAMRAI_MPI& mpi(SAMRAI_MPI::getSAMRAIWorld());

   double sum = 0.;

   if (mpi.getRank() == 0) {

      if (d_must_call_finalize) {
         TBOX_ERROR("Statistician::gettGlobalPatchStatSum ..."
            << "\n   The finalize() method to construct global data "
            "must be called BEFORE this method." << std::endl);

      }
      if (!d_has_gathered_stats) {
         TBOX_ERROR("Statistician::getGlobalProcStatSeqLength ..."
            << "\n   The finalize() method to construct global data "
            << "must be called with the argument to gather global "
            << "stats data BEFORE this metho." << std::endl);
      }

      TBOX_ASSERT(patch_stat_id >= 0 &&
         patch_stat_id < d_global_patch_stat_data.getSize());
      TBOX_ASSERT(seq_num >= 0 &&
         seq_num < d_global_patch_stat_data[patch_stat_id].getSize());

      int num_patches = d_global_patch_stat_data[patch_stat_id][seq_num].
         getSize();

      for (int np = 0; np < num_patches; np++) {
         sum += d_global_patch_stat_data[patch_stat_id][seq_num][np];
      }
   }

   return sum;
}

double
Statistician::getGlobalPatchStatMax(
   int patch_stat_id,
   int seq_num)
{
   const SAMRAI_MPI& mpi(SAMRAI_MPI::getSAMRAIWorld());
   double pmax = -(MathUtilities<double>::getMax());

   if (mpi.getRank() == 0) {

      if (d_must_call_finalize) {
         TBOX_ERROR("Statistician::getGlobalPatchStatMax ..."
            << "\n   The finalize() method to construct global data "
            "must be called BEFORE this method." << std::endl);

      }
      if (!d_has_gathered_stats) {
         TBOX_ERROR("Statistician::getGlobalProcStatSeqLength ..."
            << "\n   The finalize() method to construct global data "
            << "must be called with the argument to gather global "
            << "stats data BEFORE this metho." << std::endl);
      }

      TBOX_ASSERT(patch_stat_id >= 0 &&
         patch_stat_id < d_global_patch_stat_data.getSize());
      TBOX_ASSERT(seq_num >= 0 &&
         seq_num < d_global_patch_stat_data[patch_stat_id].getSize());

      int num_patches = d_global_patch_stat_data[patch_stat_id][seq_num].
         getSize();
      double val = pmax;
      for (int np = 0; np < num_patches; np++) {
         val = d_global_patch_stat_data[patch_stat_id][seq_num][np];
         pmax = (val > pmax ? val : pmax);
      }
   }

   return pmax;
}

int
Statistician::getGlobalPatchStatMaxPatchId(
   int patch_stat_id,
   int seq_num)
{
   const SAMRAI_MPI& mpi(SAMRAI_MPI::getSAMRAIWorld());
   double pmax = -(MathUtilities<double>::getMax());
   int id = -1;

   if (mpi.getRank() == 0) {

      if (d_must_call_finalize) {
         TBOX_ERROR("Statistician::getGlobalPatchStatMaxPatchId ..."
            << "\n   The finalize() method to construct global data "
            "must be called BEFORE this method." << std::endl);

      }
      if (!d_has_gathered_stats) {
         TBOX_ERROR("Statistician::getGlobalProcStatSeqLength ..."
            << "\n   The finalize() method to construct global data "
            << "must be called with the argument to gather global "
            << "stats data BEFORE this metho." << std::endl);
      }

      TBOX_ASSERT(patch_stat_id >= 0 &&
         patch_stat_id < d_global_patch_stat_data.getSize());
      TBOX_ASSERT(seq_num >= 0 &&
         seq_num < d_global_patch_stat_data[patch_stat_id].getSize());

      int num_patches = d_global_patch_stat_data[patch_stat_id][seq_num].
         getSize();
      double val = pmax;
      for (int np = 0; np < num_patches; np++) {
         val = d_global_patch_stat_data[patch_stat_id][seq_num][np];
         if (val > pmax) {
            id = np;
            pmax = val;
         }
      }

   }
   return id;
}

double
Statistician::getGlobalPatchStatMin(
   int patch_stat_id,
   int seq_num)
{
   const SAMRAI_MPI& mpi(SAMRAI_MPI::getSAMRAIWorld());
   double pmin = MathUtilities<double>::getMax();

   if (mpi.getRank() == 0) {

      if (d_must_call_finalize) {
         TBOX_ERROR("Statistician::getGlobalPatchStatMin ..."
            << "\n   The finalize() method to construct global data "
            "must be called BEFORE this method." << std::endl);

      }
      if (!d_has_gathered_stats) {
         TBOX_ERROR("Statistician::getGlobalProcStatSeqLength ..."
            << "\n   The finalize() method to construct global data "
            << "must be called with the argument to gather global "
            << "stats data BEFORE this metho." << std::endl);
      }

      TBOX_ASSERT(patch_stat_id >= 0 &&
         patch_stat_id < d_global_patch_stat_data.getSize());
      TBOX_ASSERT(seq_num >= 0 &&
         seq_num < d_global_patch_stat_data[patch_stat_id].getSize());

      int num_patches = d_global_patch_stat_data[patch_stat_id][seq_num].
         getSize();
      double val = pmin;
      for (int np = 0; np < num_patches; np++) {
         val = d_global_patch_stat_data[patch_stat_id][seq_num][np];
         pmin = (val < pmin ? val : pmin);
      }
   }

   return pmin;
}

int
Statistician::getGlobalPatchStatMinPatchId(
   int patch_stat_id,
   int seq_num)
{
   const SAMRAI_MPI& mpi(SAMRAI_MPI::getSAMRAIWorld());
   double pmin = MathUtilities<double>::getMax();
   int id = -1;

   if (mpi.getRank() == 0) {

      if (d_must_call_finalize) {
         TBOX_ERROR("Statistician::getGlobalPatchStatMinPatchId ..."
            << "\n   The finalize() method to construct global data "
            "must be called BEFORE this method." << std::endl);

      }
      if (!d_has_gathered_stats) {
         TBOX_ERROR("Statistician::getGlobalProcStatSeqLength ..."
            << "\n   The finalize() method to construct global data "
            << "must be called with the argument to gather global "
            << "stats data BEFORE this metho." << std::endl);
      }

      TBOX_ASSERT(patch_stat_id >= 0 &&
         patch_stat_id < d_global_patch_stat_data.getSize());
      TBOX_ASSERT(seq_num >= 0 &&
         seq_num < d_global_patch_stat_data[patch_stat_id].getSize());

      int num_patches = d_global_patch_stat_data[patch_stat_id][seq_num].
         getSize();
      double val = pmin;
      for (int np = 0; np < num_patches; np++) {
         val = d_global_patch_stat_data[patch_stat_id][seq_num][np];
         if (val < pmin) {
            id = np;
            pmin = val;
         }
      }
   }

   return id;

}

double
Statistician::getGlobalPatchStatProcessorSum(
   int patch_stat_id,
   int processor_id,
   int seq_num)
{
   const SAMRAI_MPI& mpi(SAMRAI_MPI::getSAMRAIWorld());
   double sum = -1.;

   if (mpi.getRank() == 0) {

      if (d_must_call_finalize) {
         TBOX_ERROR("Statistician::getGlobalPatchStatProcSum ..."
            << "\n   The finalize() method to construct global data "
            "must be called BEFORE this method." << std::endl);

      }
      if (!d_has_gathered_stats) {
         TBOX_ERROR("Statistician::getGlobalProcStatSeqLength ..."
            << "\n   The finalize() method to construct global data "
            << "must be called with the argument to gather global "
            << "stats data BEFORE this metho." << std::endl);
      }

      TBOX_ASSERT(patch_stat_id >= 0 &&
         patch_stat_id < d_global_patch_stat_proc_data.getSize());
      TBOX_ASSERT(processor_id >= 0 &&
         processor_id < mpi.getSize());
      TBOX_ASSERT(seq_num >= 0 &&
         seq_num < d_global_patch_stat_proc_data[patch_stat_id].getSize());
      sum =
         d_global_patch_stat_proc_data[patch_stat_id][seq_num][processor_id];
   }

   return sum;
}

double
Statistician::getGlobalPatchStatProcessorSumMax(
   int patch_stat_id,
   int seq_num)
{
   const SAMRAI_MPI& mpi(SAMRAI_MPI::getSAMRAIWorld());
   double pmax = -(MathUtilities<double>::getMax());

   if (mpi.getRank() == 0) {

      if (d_must_call_finalize) {
         TBOX_ERROR("Statistician::gettGlobalPatchStatProcSumMax ..."
            << "\n   The finalize() method to construct global data "
            "must be called BEFORE this method." << std::endl);

      }
      if (!d_has_gathered_stats) {
         TBOX_ERROR("Statistician::getGlobalProcStatSeqLength ..."
            << "\n   The finalize() method to construct global data "
            << "must be called with the argument to gather global "
            << "stats data BEFORE this metho." << std::endl);
      }

      TBOX_ASSERT(patch_stat_id >= 0 &&
         patch_stat_id < d_global_patch_stat_proc_data.getSize());
      TBOX_ASSERT(seq_num >= 0 &&
         seq_num < d_global_patch_stat_proc_data[patch_stat_id].getSize());

      double val = pmax;
      for (int np = 0; np < mpi.getSize(); np++) {
         val = d_global_patch_stat_proc_data[patch_stat_id][seq_num][np];
         pmax = (val > pmax ? val : pmax);
      }
   }

   return pmax;
}

int
Statistician::getGlobalPatchStatProcessorSumMaxId(
   int patch_stat_id,
   int seq_num)
{
   const SAMRAI_MPI& mpi(SAMRAI_MPI::getSAMRAIWorld());
   double pmax = -(MathUtilities<double>::getMax());
   int id = -1;

   if (mpi.getRank() == 0) {

      if (d_must_call_finalize) {
         TBOX_ERROR("Statistician::getGlobalPatchStatMaxProcSumId ..."
            << "\n   The finalize() method to construct global data "
            "must be called BEFORE this method." << std::endl);

      }
      if (!d_has_gathered_stats) {
         TBOX_ERROR("Statistician::getGlobalProcStatSeqLength ..."
            << "\n   The finalize() method to construct global data "
            << "must be called with the argument to gather global "
            << "stats data BEFORE this metho." << std::endl);
      }

      TBOX_ASSERT(patch_stat_id >= 0 &&
         patch_stat_id < d_global_patch_stat_proc_data.getSize());
      TBOX_ASSERT(seq_num >= 0 &&
         seq_num < d_global_patch_stat_proc_data[patch_stat_id].getSize());

      double val = pmax;
      for (int np = 0; np < mpi.getSize(); np++) {
         val = d_global_patch_stat_proc_data[patch_stat_id][seq_num][np];
         if (val > pmax) {
            id = np;
            pmax = val;
         }
      }
   }
   return id;
}

double
Statistician::getGlobalPatchStatProcessorSumMin(
   int patch_stat_id,
   int seq_num)
{
   const SAMRAI_MPI& mpi(SAMRAI_MPI::getSAMRAIWorld());
   double pmin = MathUtilities<double>::getMax();

   if (mpi.getRank() == 0) {

      if (d_must_call_finalize) {
         TBOX_ERROR("Statistician::getGlobalPatchStatProcSumMin ..."
            << "\n   The finalize() method to construct global data "
            "must be called BEFORE this method." << std::endl);

      }
      if (!d_has_gathered_stats) {
         TBOX_ERROR("Statistician::getGlobalProcStatSeqLength ..."
            << "\n   The finalize() method to construct global data "
            << "must be called with the argument to gather global "
            << "stats data BEFORE this metho." << std::endl);
      }

      TBOX_ASSERT(patch_stat_id >= 0 &&
         patch_stat_id < d_global_patch_stat_proc_data.getSize());
      TBOX_ASSERT(seq_num >= 0 &&
         seq_num < d_global_patch_stat_proc_data[patch_stat_id].getSize());

      double val = pmin;

      for (int np = 0; np < mpi.getSize(); np++) {
         val = d_global_patch_stat_proc_data[patch_stat_id][seq_num][np];
         pmin = (val < pmin ? val : pmin);
      }
   }

   return pmin;
}

int
Statistician::getGlobalPatchStatProcessorSumMinId(
   int patch_stat_id,
   int seq_num)
{
   const SAMRAI_MPI& mpi(SAMRAI_MPI::getSAMRAIWorld());
   double pmin = MathUtilities<double>::getMax();
   int id = -1;

   if (mpi.getRank() == 0) {

      if (d_must_call_finalize) {
         TBOX_ERROR("Statistician::getGlobalPatchStatProcSumMinId ..."
            << "\n   The finalize() method to construct global data "
            "must be called BEFORE this method." << std::endl);

      }
      if (!d_has_gathered_stats) {
         TBOX_ERROR("Statistician::getGlobalProcStatSeqLength ..."
            << "\n   The finalize() method to construct global data "
            << "must be called with the argument to gather global "
            << "stats data BEFORE this metho." << std::endl);
      }

      TBOX_ASSERT(patch_stat_id >= 0 &&
         patch_stat_id < d_global_patch_stat_proc_data.getSize());
      TBOX_ASSERT(seq_num >= 0 &&
         seq_num < d_global_patch_stat_proc_data[patch_stat_id].getSize());

      double val = pmin;
      for (int np = 0; np < mpi.getSize(); np++) {
         val = d_global_patch_stat_proc_data[patch_stat_id][seq_num][np];
         if (val < pmin) {
            id = np;
            pmin = val;
         }
      }
   }
   return id;
}

int
Statistician::getGlobalPatchStatNumberPatchesOnProc(
   int patch_stat_id,
   int seq_num,
   int proc_id)
{
   const SAMRAI_MPI& mpi(SAMRAI_MPI::getSAMRAIWorld());
   TBOX_ASSERT(patch_stat_id >= 0 &&
      patch_stat_id < d_global_patch_stat_data.getSize());
   TBOX_ASSERT(seq_num >= 0 &&
      seq_num < d_global_patch_stat_data[patch_stat_id].getSize());
   TBOX_ASSERT(proc_id >= 0 &&
      proc_id < mpi.getRank());
   int num_patches = -1;

   if (mpi.getRank() == 0) {

      if (d_must_call_finalize) {
         TBOX_ERROR("Statistician::getGlobalPatchStatNumberPatchesOnProc"
            << "\n   The finalize() method to construct global data "
            << "must be called BEFORE this method." << std::endl);
      }
      if (!d_has_gathered_stats) {
         TBOX_ERROR("Statistician::getGlobalProcStatSeqLength ..."
            << "\n   The finalize() method to construct global data "
            << "must be called with the argument to gather global "
            << "stats data BEFORE this metho." << std::endl);
      }

      num_patches =
         d_global_patch_stat_mapping[patch_stat_id][seq_num][proc_id];
   }

   return num_patches;

}

int
Statistician::getGlobalPatchStatMaxPatchesPerProc(
   int patch_stat_id,
   int seq_num)
{
   const SAMRAI_MPI& mpi(SAMRAI_MPI::getSAMRAIWorld());
   TBOX_ASSERT(patch_stat_id >= 0 &&
      patch_stat_id < d_global_patch_stat_proc_data.getSize());
   TBOX_ASSERT(seq_num >= 0 &&
      seq_num < d_global_patch_stat_proc_data[patch_stat_id].getSize());

   int pmax = -999999;

   if (mpi.getRank() == 0) {

      if (d_must_call_finalize) {
         TBOX_ERROR("Statistician::getGlobalPatchStatMaxPatchesPerProc"
            << "\n   The finalize() method to construct global data "
            "must be called BEFORE this method." << std::endl);

      }
      if (!d_has_gathered_stats) {
         TBOX_ERROR("Statistician::getGlobalProcStatSeqLength ..."
            << "\n   The finalize() method to construct global data "
            << "must be called with the argument to gather global "
            << "stats data BEFORE this metho." << std::endl);
      }

      int val = pmax;
      Array<int> patches_per_proc;
      patches_per_proc.resizeArray(mpi.getSize());

      int num_patches = d_global_patch_stat_data[patch_stat_id][seq_num].
         getSize();

      int np, p;
      for (np = 0; np < mpi.getSize(); np++) {
         patches_per_proc[np] = 0;
      }

      for (p = 0; p < num_patches; p++) {
         np = d_global_patch_stat_mapping[patch_stat_id][seq_num][np];
         patches_per_proc[np]++;
      }

      for (np = 0; np < mpi.getSize(); np++) {
         val = patches_per_proc[np];
         pmax = (val > pmax ? val : pmax);
      }
   }

   return pmax;
}

int
Statistician::getGlobalPatchStatMaxPatchesPerProcId(
   int patch_stat_id,
   int seq_num)
{
   const SAMRAI_MPI& mpi(SAMRAI_MPI::getSAMRAIWorld());
   TBOX_ASSERT(patch_stat_id >= 0 &&
      patch_stat_id < d_global_patch_stat_proc_data.getSize());
   TBOX_ASSERT(seq_num >= 0 &&
      seq_num < d_global_patch_stat_proc_data[patch_stat_id].getSize());

   int pmax = -999999;
   int id = -1;

   if (mpi.getRank() == 0) {

      if (d_must_call_finalize) {
         TBOX_ERROR("Statistician::getGlobalPatchStatMaxPatchesPerProcId"
            << "\n   The finalize() method to construct global data "
            "must be called BEFORE this method." << std::endl);

      }
      if (!d_has_gathered_stats) {
         TBOX_ERROR("Statistician::getGlobalProcStatSeqLength ..."
            << "\n   The finalize() method to construct global data "
            << "must be called with the argument to gather global "
            << "stats data BEFORE this metho." << std::endl);
      }

      int val = pmax;
      Array<int> patches_per_proc;
      patches_per_proc.resizeArray(mpi.getSize());

      int num_patches = d_global_patch_stat_data[patch_stat_id][seq_num].
         getSize();

      int np, p;
      for (np = 0; np < mpi.getSize(); np++) {
         patches_per_proc[np] = 0;
      }

      for (p = 0; p < num_patches; p++) {
         np = d_global_patch_stat_mapping[patch_stat_id][seq_num][np];
         patches_per_proc[np]++;
      }

      for (np = 0; np < mpi.getSize(); np++) {
         val = patches_per_proc[np];
         if (val > pmax) {
            id = np;
            pmax = val;
         }

      }
   }

   return id;
}

int
Statistician::getGlobalPatchStatMinPatchesPerProc(
   int patch_stat_id,
   int seq_num)
{
   const SAMRAI_MPI& mpi(SAMRAI_MPI::getSAMRAIWorld());
   TBOX_ASSERT(patch_stat_id >= 0 &&
      patch_stat_id < d_global_patch_stat_proc_data.getSize());
   TBOX_ASSERT(seq_num >= 0 &&
      seq_num < d_global_patch_stat_proc_data[patch_stat_id].getSize());

   int pmin = 999999;

   if (mpi.getRank() == 0) {

      if (d_must_call_finalize) {
         TBOX_ERROR("Statistician::getGlobalPatchStatMinPatchesPerProc"
            << "\n   The finalize() method to construct global data "
            "must be called BEFORE this method." << std::endl);

      }
      if (!d_has_gathered_stats) {
         TBOX_ERROR("Statistician::getGlobalProcStatSeqLength ..."
            << "\n   The finalize() method to construct global data "
            << "must be called with the argument to gather global "
            << "stats data BEFORE this metho." << std::endl);
      }

      TBOX_ASSERT(patch_stat_id >= 0 &&
         patch_stat_id < d_global_patch_stat_data.getSize());
      TBOX_ASSERT(seq_num >= 0 &&
         seq_num < d_global_patch_stat_data[patch_stat_id].getSize());

      int val = pmin;
      Array<int> patches_per_proc;
      patches_per_proc.resizeArray(mpi.getSize());

      int num_patches = d_global_patch_stat_data[patch_stat_id][seq_num].
         getSize();

      int np, p;
      for (np = 0; np < mpi.getSize(); np++) {
         patches_per_proc[np] = 0;
      }

      for (p = 0; p < num_patches; p++) {
         np = d_global_patch_stat_mapping[patch_stat_id][seq_num][np];
         patches_per_proc[np]++;
      }

      for (np = 0; np < mpi.getSize(); np++) {
         val = patches_per_proc[np];
         pmin = (val < pmin ? val : pmin);
      }

   }

   return pmin;
}

int
Statistician::getGlobalPatchStatMinPatchesPerProcId(
   int patch_stat_id,
   int seq_num)
{
   const SAMRAI_MPI& mpi(SAMRAI_MPI::getSAMRAIWorld());
   TBOX_ASSERT(patch_stat_id >= 0 &&
      patch_stat_id < d_global_patch_stat_proc_data.getSize());
   TBOX_ASSERT(seq_num >= 0 &&
      seq_num < d_global_patch_stat_proc_data[patch_stat_id].getSize());

   int pmin = 9999999;
   int id = -1;

   if (mpi.getRank() == 0) {

      if (d_must_call_finalize) {
         TBOX_ERROR("Statistician::getGlobalPatchStatMinPatchesPerProcId"
            << "\n   The finalize() method to construct global data "
            "must be called BEFORE this method." << std::endl);

      }
      if (!d_has_gathered_stats) {
         TBOX_ERROR("Statistician::getGlobalProcStatSeqLength ..."
            << "\n   The finalize() method to construct global data "
            << "must be called with the argument to gather global "
            << "stats data BEFORE this metho." << std::endl);
      }

      int val = pmin;
      Array<int> patches_per_proc;
      patches_per_proc.resizeArray(mpi.getSize());

      int num_patches = d_global_patch_stat_data[patch_stat_id][seq_num].
         getSize();

      int np, p;
      for (np = 0; np < mpi.getSize(); np++) {
         patches_per_proc[np] = 0;
      }

      for (p = 0; p < num_patches; p++) {
         np = d_global_patch_stat_mapping[patch_stat_id][seq_num][np];
         patches_per_proc[np]++;
      }

      for (np = 0; np < mpi.getSize(); np++) {
         val = patches_per_proc[np];
         if (val < pmin) {
            id = np;
            pmin = val;
         }

      }
   }

   return id;
}

void
Statistician::printGlobalPatchStatData(
   int patch_stat_id,
   std::ostream& os,
   int precision)
{
   const SAMRAI_MPI& mpi(SAMRAI_MPI::getSAMRAIWorld());

   if (mpi.getRank() == 0) {

      if (d_must_call_finalize) {
         TBOX_ERROR("Statistician::printGlobalPatchStatData ..."
            << "\n   The finalize() method to construct global data "
            "must be called BEFORE this method." << std::endl);
      }
      if (!d_has_gathered_stats) {
         TBOX_ERROR("Statistician::getGlobalProcStatSeqLength ..."
            << "\n   The finalize() method to construct global data "
            << "must be called with the argument to gather global "
            << "stats data BEFORE this metho." << std::endl);
      }

      TBOX_ASSERT(patch_stat_id >= 0);
      TBOX_ASSERT(precision > 0);

      os.precision(precision);

      os << "\n   " << patch_stat_id << ":    "
         << d_patch_statistics[patch_stat_id]->getName() << std::endl;

      const Array<Array<double> >& sdata =
         d_global_patch_stat_data[patch_stat_id];
      const Array<Array<int> >& spmap =
         d_global_patch_stat_mapping[patch_stat_id];
      for (int ipsl = 0; ipsl < sdata.getSize(); ipsl++) {
         os << "      Seq # " << ipsl << std::endl;
         os << "         patch[proc]: value" << std::endl;
         for (int ip = 0; ip < sdata[ipsl].getSize(); ip++) {
            /*
             * Write out data only if data entry is NOT an "empty"
             * entry, defined by the Statistic::s_empty_seq_tag_entry
             * value.
             */
            if (!(MathUtilities<double>::equalEps(sdata[ipsl][ip],
                     Statistic::s_empty_seq_tag_entry))) {
               os << "         " << ip << "    [" << spmap[ipsl][ip]
                  << "]:    " << sdata[ipsl][ip] << std::endl;
            }
         }
      }
   }
}

void
Statistician::printGlobalPatchStatDataFormatted(
   int patch_stat_id,
   std::ostream& os,
   int precision)
{
   const SAMRAI_MPI& mpi(SAMRAI_MPI::getSAMRAIWorld());
   if (mpi.getRank() == 0) {

      if (d_must_call_finalize) {
         TBOX_ERROR("Statistician::printGlobalPatchStatFormatted ..."
            << "\n   The finalize() method to construct global data "
            "must be called BEFORE this method." << std::endl);

      }
      if (!d_has_gathered_stats) {
         TBOX_ERROR("Statistician::getGlobalProcStatSeqLength ..."
            << "\n   The finalize() method to construct global data "
            << "must be called with the argument to gather global "
            << "stats data BEFORE this metho." << std::endl);
      }

      TBOX_ASSERT(patch_stat_id >= 0);
      TBOX_ASSERT(precision > 0);

      os.precision(precision);

      /*
       * Output in C++ is by default right justified with the setw()
       * option e.g. cout << "[" << setw(5) << 1 << "]" will output
       * [   1].  Using setf(ios::left) makes it left justified, which
       * is more convenient to output columns of tables.
       */
      os.setf(std::ios::left);
      int s, n;
      int npatches = 0;

      /*
       * Print tab-separated table header
       */
      const Array<Array<double> >& sdata =
         d_global_patch_stat_data[patch_stat_id];
      for (s = 0; s < sdata.getSize(); s++) {
         int np_at_s = sdata[s].getSize();
         npatches = (np_at_s > npatches ? np_at_s : npatches);
      }

      // heading - line 1
      os << "Seq#\t";
      for (n = 0; n < npatches; n++) {
         os << "Patch\t";
      }
      os << std::endl;

      // heading - line 2
      os << "    \t";
      for (n = 0; n < npatches; n++) {
         os << n << "\t";
      }
      os << std::endl;

      /*
       * Now print values.
       */
      for (s = 0; s < sdata.getSize(); s++) {
         os << s << "\t";
         for (n = 0; n < sdata[s].getSize(); n++) {
            /*
             * Write out data only if data entry is NOT an "empty"
             * entry, defined by the Statistic::s_empty_seq_tag_entry
             * value.
             */
            if (MathUtilities<double>::equalEps(sdata[s][n],
                   Statistic::s_empty_seq_tag_entry)) {
               os << "  " << "\t";
            } else {
               os << sdata[s][n] << "\t";
            }
         }
         os << std::endl;
      }
   }

}

/*
 *************************************************************************
 *
 * Gather all statistic data to processor zero for analysis.
 *
 *************************************************************************
 */

void
Statistician::finalize(
   bool gather_individual_stats_on_proc_0)
{
   const SAMRAI_MPI& mpi(SAMRAI_MPI::getSAMRAIWorld());
   int nnodes = mpi.getSize();
   int my_rank = mpi.getRank();

   d_global_proc_stat_data.resizeArray(0);
   d_global_patch_stat_data.resizeArray(0);
   d_global_patch_stat_proc_data.resizeArray(0);
   d_global_patch_stat_mapping.resizeArray(0);

   Array<int> total_patches;
   checkStatsForConsistency(total_patches);

   reduceGlobalStatistics();

   if (gather_individual_stats_on_proc_0) {

      /*
       * On processor zero, we create array structures to assemble
       * global statistic information.   Note that the array "total_patches"
       * holds the total number of patches for each sequence entry for
       * each patch statistic.
       */

      int is, ip, ipsl, seq_len;

      /*
       * Create and initialize global_proc and global_patch stats arrays.
       * These are only needed on processor 0 but we construct them on
       * all processors to avoid compiler warnings.
       */
      Array<boost::shared_ptr<Statistic> >* global_proc_stats =
         d_num_proc_stats > 0
         ? new Array<boost::shared_ptr<Statistic> >[d_num_proc_stats]
         : (Array<boost::shared_ptr<Statistic> > *)NULL;
      Array<boost::shared_ptr<Statistic> >* global_patch_stats =
         d_num_patch_stats > 0
         ? new Array<boost::shared_ptr<Statistic> >[d_num_patch_stats]
         : (Array<boost::shared_ptr<Statistic> > *)NULL;

      if (my_rank == 0) {

         d_global_proc_stat_data.resizeArray(d_num_proc_stats);
         for (is = 0; is < d_num_proc_stats; is++) {

            seq_len = d_proc_statistics[is]->getStatSequenceLength();
            d_global_proc_stat_data[is].resizeArray(seq_len);

            for (ip = 0; ip < seq_len; ip++) {
               d_global_proc_stat_data[is][ip].resizeArray(nnodes);
            }

         }

         d_global_patch_stat_data.resizeArray(d_num_patch_stats);
         d_global_patch_stat_mapping.resizeArray(d_num_patch_stats);
         d_global_patch_stat_proc_data.resizeArray(d_num_patch_stats);

         ipsl = 0;
         for (is = 0; is < d_num_patch_stats; is++) {

            seq_len = d_patch_statistics[is]->getStatSequenceLength();
            d_global_patch_stat_data[is].resizeArray(seq_len);
            d_global_patch_stat_mapping[is].resizeArray(seq_len);
            d_global_patch_stat_proc_data[is].resizeArray(seq_len);

            for (ip = 0; ip < seq_len; ip++) {
               int npatches = total_patches[ipsl];
               d_global_patch_stat_data[is][ip].
               resizeArray(npatches);
               d_global_patch_stat_mapping[is][ip].
               resizeArray(npatches);
               d_global_patch_stat_proc_data[is][ip].
               resizeArray(nnodes);
               ipsl++;
            }

         }

      }

      /*
       * Now we need to gather statistic information from each processor
       * to fill arrays.  If my_rank == 0, then I am a destination processor;
       * otherwise, I am a source processor.
       */

      if (nnodes > 1) {

         Schedule stat_schedule;
         stat_schedule.setTimerPrefix("tbox::Statistician");

         if (d_num_proc_stats > 0) {

            if (my_rank == 0) {

               /*
                * Create storage for remote statistic information from each
                * processor and add transaction to schedule to receive stat.
                */

               for (is = 0; is < d_num_proc_stats; is++) {
                  const std::string& sname = d_proc_statistics[is]->getName();
                  const std::string& stype = d_proc_statistics[is]->getType();
                  global_proc_stats[is].resizeArray(nnodes);
                  for (ip = 1; ip < nnodes; ip++) {
                     global_proc_stats[is][ip].reset(
                        new Statistic(sname, stype, is));
                     stat_schedule.addTransaction(
                        boost::make_shared<StatTransaction>(
                           global_proc_stats[is][ip], ip, 0));
                  }
               }

            } else {

               /*
                * Add transaction to schedule to send statistic information
                * to processor zero.
                */

               for (is = 0; is < d_num_proc_stats; is++) {
                  stat_schedule.addTransaction(
                     boost::make_shared<StatTransaction>(
                        d_proc_statistics[is], my_rank, 0));
               }

            }

         }

         if (d_num_patch_stats > 0) {

            if (my_rank == 0) {

               /*
                * Create storage for remote statistic information from each
                * processor and add transaction to schedule to receive stat.
                */

               for (is = 0; is < d_num_patch_stats; is++) {
                  const std::string& sname = d_patch_statistics[is]->getName();
                  const std::string& stype = d_patch_statistics[is]->getType();
                  global_patch_stats[is].resizeArray(nnodes);
                  for (ip = 1; ip < nnodes; ip++) {
                     global_patch_stats[is][ip].reset(
                        new Statistic(sname, stype, is));
                     stat_schedule.addTransaction(
                        boost::make_shared<StatTransaction>(
                           global_patch_stats[is][ip], ip, 0));
                  }
               }

            } else {

               /*
                * Add transaction to schedule to send statistic information
                * to processor zero.
                */

               for (is = 0; is < d_num_patch_stats; is++) {
                  stat_schedule.addTransaction(
                     boost::make_shared<StatTransaction>(
                        d_patch_statistics[is], my_rank, 0));
               }

            }

         }

         /*
          * Communicate all statistic data to processor zero.
          */

         stat_schedule.communicate();

      }

      /*
       * Fill arrays with global statistic data.
       */

      if (my_rank == 0) {

         for (is = 0; is < d_num_proc_stats; is++) {

            Array<Array<double> >& sdata =
               d_global_proc_stat_data[is];

            for (ip = 0; ip < nnodes; ip++) {

               boost::shared_ptr<Statistic> stat(((ip == 0) ?
                                          d_proc_statistics[is] :
                                          global_proc_stats[is][ip]));

               for (ipsl = 0; ipsl < stat->getStatSequenceLength(); ipsl++) {
                  sdata[ipsl][ip] = stat->getProcStatSeqArray()[ipsl].value;
               }

            } // iterate over processors

         } // iterate over proc stats

         for (is = 0; is < d_num_patch_stats; is++) {

            Array<Array<double> >& sdata =
               d_global_patch_stat_data[is];
            Array<Array<int> >& spmap =
               d_global_patch_stat_mapping[is];

            for (ip = 0; ip < nnodes; ip++) {

               boost::shared_ptr<Statistic> stat(((ip == 0) ?
                                          d_patch_statistics[is] :
                                          global_patch_stats[is][ip]));

               for (ipsl = 0; ipsl < stat->getStatSequenceLength(); ipsl++) {
                  const std::list<Statistic::PatchStatRecord>& psrl =
                     stat->getPatchStatSeqArray()[ipsl].patch_records;
                  std::list<Statistic::PatchStatRecord>::const_iterator ilr =
                     psrl.begin();
                  for ( ; ilr != psrl.end(); ilr++) {
                     int patch_id = ilr->patch_id;
                     sdata[ipsl][patch_id] = ilr->value;
                     spmap[ipsl][patch_id] = ip;
                  }

               }

            } // iterate over processors

         } // iterate over patch stats

         /*
          * Construct patch stat procesor array, which holds the sum of
          * values on each processor for the patch stats.  In effect, this
          * just constructs a processor stat from the patch stat data.
          */

         for (is = 0; is < d_num_patch_stats; is++) {

            const Array<Array<double> >& pdata =
               d_global_patch_stat_data[is];

            for (int s = 0; s < pdata.getSize(); s++) {

               for (ip = 0; ip < nnodes; ip++) {
                  d_global_patch_stat_proc_data[is][s][ip] = 0.;
               }

               for (int p = 0; p < pdata[s].getSize(); p++) {
                  ip = d_global_patch_stat_mapping[is][s][p];
                  d_global_patch_stat_proc_data[is][s][ip] +=
                     d_global_patch_stat_data[is][s][p];
               }

            } // iterate over sequence of patch stat is

         } // iterate over patch stats

      } // if I am processor zero

      delete[] global_proc_stats;
      delete[] global_patch_stats;

      d_has_gathered_stats = true;

   } else {
      d_has_gathered_stats = false;
   }

   d_must_call_finalize = false;
}

/*
 *************************************************************************
 * Globally reduce all statistic data.
 *************************************************************************
 */

void
Statistician::reduceGlobalStatistics()
{
   std::vector<double> proc_stat_values;
   for (int istat = 0; istat < d_num_proc_stats; istat++) {
      const Statistic& stat(*d_proc_statistics[istat]);
      for (int iseq = 0; iseq < stat.getStatSequenceLength(); iseq++) {
         proc_stat_values.push_back(stat.getProcStatSeqArray()[iseq].value);
      }
   }

   const SAMRAI_MPI& mpi(SAMRAI_MPI::getSAMRAIWorld());

   std::vector<double> sum_proc_stat_values(proc_stat_values);
   if (mpi.getSize() > 1) {
      mpi.Allreduce(&proc_stat_values[0],
         &sum_proc_stat_values[0],
         static_cast<int>(proc_stat_values.size()),
         MPI_DOUBLE,
         MPI_SUM);
   }

   std::vector<double> max_proc_stat_values(proc_stat_values);
   std::vector<int> imax_proc_stat_values(proc_stat_values.size(), mpi.getRank());
   if (mpi.getSize() > 1) {
      mpi.AllReduce(&max_proc_stat_values[0],
         static_cast<int>(max_proc_stat_values.size()),
         MPI_MAXLOC,
         &imax_proc_stat_values[0]);
   }

   std::vector<double> min_proc_stat_values(proc_stat_values);
   std::vector<int> imin_proc_stat_values(proc_stat_values.size(), mpi.getRank());
   if (mpi.getSize() > 1) {
      mpi.AllReduce(&min_proc_stat_values[0],
         static_cast<int>(min_proc_stat_values.size()),
         MPI_MINLOC,
         &imin_proc_stat_values[0]);
   }

   d_global_proc_stat_sum.setNull();
   d_global_proc_stat_sum.resizeArray(d_num_proc_stats);
   d_global_proc_stat_max.setNull();
   d_global_proc_stat_max.resizeArray(d_num_proc_stats);
   d_global_proc_stat_min.setNull();
   d_global_proc_stat_min.resizeArray(d_num_proc_stats);
   d_global_proc_stat_imax.setNull();
   d_global_proc_stat_imax.resizeArray(d_num_proc_stats);
   d_global_proc_stat_imin.setNull();
   d_global_proc_stat_imin.resizeArray(d_num_proc_stats);

   size_t isum(0);
   size_t imax(0);
   size_t imin(0);
   for (int istat = 0; istat < d_num_proc_stats; istat++) {
      const Statistic& stat(*d_proc_statistics[istat]);
      if (stat.getStatSequenceLength() != 0) {

         d_global_proc_stat_sum[istat].resizeArray(stat.getStatSequenceLength());
         memcpy(d_global_proc_stat_sum[istat].getPointer(),
            &sum_proc_stat_values[isum],
            sizeof(double) * stat.getStatSequenceLength());
         isum += stat.getStatSequenceLength();

         d_global_proc_stat_max[istat].resizeArray(stat.getStatSequenceLength());
         d_global_proc_stat_imax[istat].resizeArray(stat.getStatSequenceLength());
         memcpy(d_global_proc_stat_max[istat].getPointer(),
            &max_proc_stat_values[imax],
            sizeof(double) * stat.getStatSequenceLength());
         memcpy(d_global_proc_stat_imax[istat].getPointer(),
            &imax_proc_stat_values[imax],
            sizeof(int) * stat.getStatSequenceLength());
         imax += stat.getStatSequenceLength();

         d_global_proc_stat_min[istat].resizeArray(stat.getStatSequenceLength());
         d_global_proc_stat_imin[istat].resizeArray(stat.getStatSequenceLength());
         memcpy(d_global_proc_stat_min[istat].getPointer(),
            &min_proc_stat_values[imin],
            sizeof(double) * stat.getStatSequenceLength());
         memcpy(d_global_proc_stat_imin[istat].getPointer(),
            &imin_proc_stat_values[imin],
            sizeof(int) * stat.getStatSequenceLength());
         imin += stat.getStatSequenceLength();

      }

   }
   TBOX_ASSERT(isum == sum_proc_stat_values.size());
   TBOX_ASSERT(imax == max_proc_stat_values.size());
   TBOX_ASSERT(imin == min_proc_stat_values.size());
}

/*
 *************************************************************************
 *
 * Check statistic information on all processors for consistency.
 * If the number of patch statistics is non-zero, the array will be
 * resized to the sum of sequence lengths over all patch statistics.
 * Then, each entry will be the total number of global patches for
 * that particular patch stat sequence.
 *
 *************************************************************************
 */

void
Statistician::checkStatsForConsistency(
   Array<int>& total_patches)
{
   TBOX_ASSERT(total_patches.getSize() == 0);

   const SAMRAI_MPI& mpi(SAMRAI_MPI::getSAMRAIWorld());

   int ip, is, consistent, n_stats;

   /*
    * First consistency check for statistics.
    *
    * We continue if number of processor stats and number of patch stats
    * are same on every processor.  Otherwise, we abort with an error message.
    */

   n_stats = d_num_proc_stats;
   if (mpi.getSize() > 1) {
      mpi.Bcast(&n_stats, 1, MPI_INT, 0);
   }
   consistent = ((n_stats == d_num_proc_stats) ? 1 : 0);
   if (mpi.getSize() > 1) {
      int consistent1(consistent);
      mpi.Allreduce(&consistent1, &consistent, 1, MPI_INT, MPI_MIN);
   }

   if (!consistent) {
      TBOX_ERROR("Statistician::finalize error ..."
         << "\n   Number of processor stats inconsistent"
         << " across processors.  Cannot gather information." << std::endl);
   }

   n_stats = d_num_patch_stats;
   if (mpi.getSize() > 1) {
      mpi.Bcast(&n_stats, 1, MPI_INT, 0);
   }
   consistent = ((n_stats == d_num_patch_stats) ? 1 : 0);
   if (mpi.getSize() > 1) {
      int consistent1(consistent);
      mpi.Allreduce(&consistent1, &consistent, 1, MPI_INT, MPI_MIN);
   }

   if (!consistent) {
      TBOX_ERROR("Statistician::finalize error ..."
         << "\n   Number of patch stats inconsistent"
         << " across processors.  Cannot gather information." << std::endl);
   }

   /*
    * Second consistency check for statistics.
    *
    * We continue if time sequence length for each processor stat and
    * each patch stat is same on every processor.  Otherwise, we abort
    * with an error message.  Note that we implicitly assume that the
    * statistics are created in the same order on each processor.
    */

   n_stats = d_num_proc_stats + d_num_patch_stats;

   int n_patch_stat_seq_items = 0;

   if (n_stats > 0) {

      std::vector<int> my_seq_lengths(n_stats);
      std::vector<int> max_seq_lengths(n_stats);
      for (is = 0; is < d_num_proc_stats; is++) {
         my_seq_lengths[is] = d_proc_statistics[is]->getStatSequenceLength();
         max_seq_lengths[is] = my_seq_lengths[is];
      }
      for (is = 0; is < d_num_patch_stats; is++) {
         int mark = d_num_proc_stats + is;
         my_seq_lengths[mark] =
            d_patch_statistics[is]->getStatSequenceLength();
         max_seq_lengths[mark] = my_seq_lengths[mark];
         n_patch_stat_seq_items += my_seq_lengths[mark];
      }

      if (mpi.getSize() > 1) {
         std::vector<int> max_seq_lengths1(max_seq_lengths);
         mpi.Allreduce(&max_seq_lengths1[0], &max_seq_lengths[0],
            static_cast<int>(max_seq_lengths.size()), MPI_INT, MPI_MAX);
      }

      consistent = 1;
      for (is = 0; is < n_stats; is++) {
         if (max_seq_lengths[is] != my_seq_lengths[is]) {
            consistent = 0;
            break;
         }
      }
      if (mpi.getSize() > 1) {
         int consistent1(consistent);
         mpi.Allreduce(&consistent1, &consistent, 1, MPI_INT, MPI_MIN);
      }

      if (!consistent) {
         TBOX_ERROR("Statistician::finalize error ..."
            << "\n   Sequence length for some statistic inconsistent"
            << " across processors.  Cannot gather information." << std::endl);
      }

   }

   /*
    * Third, we gather the total number of patches for each sequence item
    * for each patch statitic.
    */

   total_patches.resizeArray(n_patch_stat_seq_items);
   if (d_num_patch_stats > 0) {
      int ipsl = 0;
      for (is = 0; is < d_num_patch_stats; is++) {
         for (ip = 0; ip < d_patch_statistics[is]->getStatSequenceLength();
              ip++) {
            total_patches[ipsl] = static_cast<int>(
               d_patch_statistics[is]->getPatchStatSeqArray()[ip].patch_records.size());
            ipsl++;
         }
      }

      if (mpi.getSize() > 1) {
         std::vector<int> total_patches1(total_patches.getPointer(),
                                         total_patches.getPointer() + total_patches.size());
         mpi.Allreduce(&total_patches1[0], &total_patches[0],
            total_patches.size(), MPI_INT, MPI_SUM);
      }

   }
}

/*
 *************************************************************************
 *
 * Utility functions to print local statistic information and global
 * statistic information (after calling finalize()).
 *
 *************************************************************************
 */

void
Statistician::printLocalStatData(
   std::ostream& os,
   int precision) const
{
   int is;

   os << "\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++"
      << std::endl;
   os << "Printing local statistic information...";
   os << "\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++"
      << std::endl;

   for (is = 0; is < d_num_proc_stats; is++) {
      os << std::endl;
      d_proc_statistics[is]->printClassData(os, precision);
   }

   for (is = 0; is < d_num_patch_stats; is++) {
      os << std::endl;
      d_patch_statistics[is]->printClassData(os, precision);
   }
}

void
Statistician::printAllGlobalStatData(
   std::ostream& os,
   int precision)
{
   const SAMRAI_MPI& mpi(SAMRAI_MPI::getSAMRAIWorld());

   if (mpi.getRank() == 0) {

      if (d_must_call_finalize) {
         TBOX_ERROR("Statistician::printAllGlobalStatData ..."
            << "\n   The finalize() method to construct global data "
            "must be called BEFORE this method." << std::endl);
      }
      if (!d_has_gathered_stats) {
         TBOX_ERROR("Statistician::getGlobalProcStatSeqLength ..."
            << "\n   The finalize() method to construct global data "
            << "must be called with the argument to gather global "
            << "stats data BEFORE this metho." << std::endl);
      }

      int is;
      os << "\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++"
         << "\nPrinting global statistic information..."
         << "\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++"
         << std::endl;

      os << "\n  ---------------------" << std::endl;
      os << "  Processor Statistics:" << std::endl;
      os << "  ---------------------" << std::endl;
      os << "   Stat #: name" << std::endl;
      for (is = 0; is < d_num_proc_stats; is++) {
         printGlobalProcStatData(is, os, precision);
      }

      os << "\n  ---------------------" << std::endl;
      os << "  Patch Statistics:" << std::endl;
      os << "  ---------------------" << std::endl;
      os << "   Stat #: name" << std::endl;
      for (is = 0; is < d_num_patch_stats; is++) {
         printGlobalPatchStatData(is, os, precision);
      }

   }

}

/*
 *************************************************************************
 *
 * Print global summed output for all statistics to the specified file.
 *
 *************************************************************************
 */
void
Statistician::printAllSummedGlobalStatData(
   const std::string& filename,
   int precision)
{
   const SAMRAI_MPI& mpi(SAMRAI_MPI::getSAMRAIWorld());
   if (mpi.getRank() == 0) {
      std::ofstream file(filename.c_str());
      printAllSummedGlobalStatData(file, precision);
      file.close();

   }
}

void
Statistician::printAllSummedGlobalStatData(
   std::ostream& os,
   int precision)
{
   const SAMRAI_MPI& mpi(SAMRAI_MPI::getSAMRAIWorld());

   if (mpi.getRank() == 0) {
      os.precision(precision);

      int is, id, num_sequences, n;
      double sum;
      for (is = 0; is < d_num_proc_stats; is++) {
         std::string procstat_name = d_proc_statistics[is]->getName();
         os << "PROCESSOR STAT: " << procstat_name << std::endl;
         id = d_proc_statistics[is]->getInstanceId();
         num_sequences = getGlobalProcStatSequenceLength(id);
         for (n = 0; n < num_sequences; n++) {
            sum = getGlobalProcStatSum(id, n);
            os << "\t" << n << "\t" << sum << std::endl;
         }
         os << "\n" << std::endl;
      }

      for (is = 0; is < d_num_patch_stats; is++) {
         std::string patchstat_name = d_patch_statistics[is]->getName();
         os << "PATCH STAT: " << patchstat_name << std::endl;
         id = d_patch_statistics[is]->getInstanceId();
         num_sequences = getGlobalPatchStatSequenceLength(id);
         for (n = 0; n < num_sequences; n++) {
            sum = getGlobalPatchStatSum(id, n);
            os << "\t" << n << "\t" << sum << std::endl;
         }
         os << "\n" << std::endl;
      }
   }

}

/*
 *************************************************************************
 *
 * Print formatted Statistician information to formatted output that
 * may be read by Excel, or other spreadsheet type programs. Each
 * processor statistic is printed to the file "<statname>-proc.txt",
 * where <statname> is the name of the statistic.  Likewise, each patch
 * statistic is written to the file "<statname>-patch.txt".  If a
 * directory name is supplied, the method will create the directory and
 * the put the files in it.  If no directory name is supplied, the files
 * will be written to the directory where the application is run.
 *
 *************************************************************************
 */

void
Statistician::printSpreadSheetOutput(
   const std::string& dirname,
   int precision)
{
   const SAMRAI_MPI& mpi(SAMRAI_MPI::getSAMRAIWorld());

   /*
    * Generate directory to write output
    */
   bool write_to_dir = false;
   if (dirname.size() > 0) {
      write_to_dir = true;
      Utilities::recursiveMkdir(dirname);
   }

   /*
    * Processor 0 writes the output files.
    */
   if (mpi.getRank() == 0) {

      if (d_must_call_finalize) {
         TBOX_ERROR("Statistician::printSpreadSheetOutput() ..."
            << "\n   The finalize() method to construct global data "
            "must be called BEFORE this method." << std::endl);
      }
      if (!d_has_gathered_stats) {
         TBOX_ERROR("Statistician::getGlobalProcStatSeqLength ..."
            << "\n   The finalize() method to construct global data "
            << "must be called with the argument to gather global "
            << "stats data BEFORE this metho." << std::endl);
      }

      /*
       * Create file in the directory for each proc and patch
       * statistic.  Files will be named <procstat-name>.procstat
       * and <patchstat-name>.patchstat.
       */
      int is;
      for (is = 0; is < d_num_proc_stats; is++) {
         std::string filename = d_proc_statistics[is]->getName();
         filename = filename + "-proc.txt";
         if (write_to_dir) filename = dirname + "/" + filename;
         std::ofstream file(filename.c_str());
         printGlobalProcStatDataFormatted(is, file, precision);
         file.close();
      }

      for (is = 0; is < d_num_patch_stats; is++) {
         std::string filename = d_patch_statistics[is]->getName();
         filename = filename + "-patch.txt";
         if (write_to_dir) filename = dirname + "/" + filename;
         std::ofstream file(filename.c_str());
         printGlobalPatchStatDataFormatted(is, file, precision);
         file.close();
      }
   }

}

/*
 *************************************************************************
 *
 * Print formatted Statistician information to formatted output that
 * may be read by Excel, or other spreadsheet type programs. This method
 * prints statistic information from only a single specified processor
 * which may be useful for information that is the same across all
 * processors.
 *
 *************************************************************************
 */

void
Statistician::printSpreadSheetOutputForProcessor(
   const int proc_id,
   const std::string& dirname,
   int precision)
{
   const SAMRAI_MPI& mpi(SAMRAI_MPI::getSAMRAIWorld());

   /*
    * Generate directory to write output
    */
   bool write_to_dir = false;
   if (dirname.size() > 0) {
      write_to_dir = true;
      Utilities::recursiveMkdir(dirname);
   }

   /*
    * Processor 0 writes the output files.  Even though we are writing
    * data from only a single processor, proc 0 will still do the
    * writing.
    */
   if (mpi.getRank() == 0) {

      if (d_must_call_finalize) {
         TBOX_ERROR("Statistician::printSpreadSheetOutput() ..."
            << "\n   The finalize() method to construct global data "
            "must be called BEFORE this method." << std::endl);
      }
      if (!d_has_gathered_stats) {
         TBOX_ERROR("Statistician::getGlobalProcStatSeqLength ..."
            << "\n   The finalize() method to construct global data "
            << "must be called with the argument to gather global "
            << "stats data BEFORE this metho." << std::endl);
      }

      /*
       * Create file in the directory for each proc and patch
       * statistic.  Files will be named <name>-<type>-proc_id.txt
       * where <name> is the name of the stat and <type> is either
       * "proc" or "patch" depending on the type of statistic.
       */
      int is;
      for (is = 0; is < d_num_proc_stats; is++) {
         std::string name = d_proc_statistics[is]->getName();
         name = name + "-proc-";
         std::string filename = name + Utilities::processorToString(
               proc_id) + ".txt";
         if (write_to_dir) filename = dirname + "/" + filename;
         std::ofstream file(filename.c_str());
         printGlobalProcStatDataFormatted(is, proc_id, file, precision);
         file.close();
      }

   }

}

/*
 *************************************************************************
 *
 * Implementation of StatisticRestartDatabase class.
 *
 *************************************************************************
 */

StatisticRestartDatabase::StatisticRestartDatabase(
   const std::string& object_name,
   bool register_for_restart,
   bool read_from_restart)
{
   TBOX_ASSERT(!object_name.empty());

   d_object_name = object_name;
   d_registered_for_restart = register_for_restart;

   if (d_registered_for_restart) {
      RestartManager::getManager()->registerRestartItem(d_object_name, this);
   }

   bool is_from_restart = RestartManager::getManager()->isFromRestart();
   if (is_from_restart && read_from_restart) {
      getFromRestart();
   }
}

StatisticRestartDatabase::~StatisticRestartDatabase()
{
   if (d_registered_for_restart) {
      RestartManager::getManager()->unregisterRestartItem(d_object_name);
   }
}

void StatisticRestartDatabase::putToDatabase(
   const boost::shared_ptr<Database>& db) const
{
   TBOX_ASSERT(db);

   db->putInteger("TBOX_STATISTICRESTARTDATABASE_VERSION",
      TBOX_STATISTICRESTARTDATABASE_VERSION);

   /*
    * Get pointer to Statistician object.
    */
   Statistician* statistician = Statistician::getStatistician();

   /*
    * Write the number of statistics
    */
   int number_of_procstats = statistician->getNumberProcessorStats();
   db->putInteger("number_of_procstats", number_of_procstats);

   int number_of_patchstats = statistician->getNumberPatchStats();
   db->putInteger("number_of_patchstats", number_of_patchstats);

   /*
    * Iterate through the list of statistics and write out a
    * sub-database for each stat, containing the data values
    * Store the name of each stat in the string arrays
    * "proc_stat_names" and "patch_stat_names".
    */
   Array<std::string> proc_stat_names(number_of_procstats);
   Array<std::string> patch_stat_names(number_of_patchstats);

   /*
    * Write procstat and patchstats to database.
    */
   boost::shared_ptr<Statistic> stat;
   boost::shared_ptr<Database> stat_database;
   int n;
   for (n = 0; n < number_of_procstats; n++) {
      stat = statistician->d_proc_statistics[n];
      proc_stat_names[n] = stat->getName();
      stat_database = db->putDatabase(proc_stat_names[n]);
      stat->putUnregisteredToDatabase(stat_database);
   }

   for (n = 0; n < number_of_patchstats; n++) {
      stat = statistician->d_patch_statistics[n];
      patch_stat_names[n] = stat->getName();
      stat_database = db->putDatabase(patch_stat_names[n]);
      stat->putUnregisteredToDatabase(stat_database);
   }

   /*
    * Write out string array containing names of stats.  This will be
    * used upon restart to specify the name of the sub-database from
    * which to read the stat info.
    */
   if (number_of_procstats > 0) {
      db->putStringArray("proc_stat_names", proc_stat_names);
   }

   if (number_of_patchstats > 0) {
      db->putStringArray("patch_stat_names", patch_stat_names);
   }
}

void StatisticRestartDatabase::getFromRestart()
{
   boost::shared_ptr<Database> root_db(
      RestartManager::getManager()->getRootDatabase());

   if (!root_db->isDatabase(d_object_name)) {
      TBOX_ERROR("Restart database corresponding to "
         << d_object_name << " not found in restart file" << std::endl);
   }
   boost::shared_ptr<Database> db(root_db->getDatabase(d_object_name));

   int ver = db->getInteger("TBOX_STATISTICRESTARTDATABASE_VERSION");
   if (ver != TBOX_STATISTICRESTARTDATABASE_VERSION) {
      TBOX_WARNING(
         d_object_name << ":  "
         "Restart file version different than class version. \n"
                       << "Cannot read statistic information from restart file so"
                       << "all statistics will be reset.");
   }

   int number_of_procstats = db->getInteger("number_of_procstats");
   int number_of_patchstats = db->getInteger("number_of_patchstats");

   /*
    * Read in the list of sub-database names.
    */
   Array<std::string> proc_stat_names;
   if (number_of_procstats > 0) {
      proc_stat_names = db->getStringArray("proc_stat_names");
   }
   Array<std::string> patch_stat_names;
   if (number_of_patchstats > 0) {
      patch_stat_names = db->getStringArray("patch_stat_names");
   }

   /*
    * Read in each stat from restart database.
    */
   Statistician* statistician = Statistician::getStatistician();
   boost::shared_ptr<Database> sub_database;
   std::string sub_database_name;
   boost::shared_ptr<Statistic> stat;
   int i;
   for (i = 0; i < number_of_procstats; i++) {
      sub_database = db->getDatabase(proc_stat_names[i]);
      stat = statistician->getStatistic(proc_stat_names[i], "PROC_STAT");
      stat->getFromRestart(sub_database);
   }

   for (i = 0; i < number_of_patchstats; i++) {
      sub_database = db->getDatabase(patch_stat_names[i]);
      stat = statistician->getStatistic(patch_stat_names[i], "PATCH_STAT");
      stat->getFromRestart(sub_database);
   }

}

}
}
