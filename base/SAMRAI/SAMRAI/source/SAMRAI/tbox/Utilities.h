/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Utility functions for error reporting, file manipulation, etc.
 *
 ************************************************************************/

#ifndef included_tbox_Utilities
#define included_tbox_Utilities

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/IOStream.h"
#include "SAMRAI/tbox/Logger.h"

#include <string>

#include <sys/types.h>
#include <sys/stat.h>

namespace SAMRAI {
namespace tbox {

#ifdef _MSC_VER
#include <direct.h>
typedef int mode_t;
#define  S_ISDIR(m) (((m) & S_IFMT) == S_IFDIR)
#define S_IRUSR 0
#define S_IWUSR 0
#define S_IXUSR 0
#endif

/*!
 * A statement that does nothing, for insure++ make it something
 * more complex than a simple C null statement to avoid a warning.
 */
#ifdef __INSURE__
#define NULL_STATEMENT if (0) int nullstatement = 0
#else
#define NULL_STATEMENT
#endif

/*!
 * A null use of a variable, use to avoid GNU compiler
 * warnings about unused variables.
 */
#define NULL_USE(variable)                               \
   do {                                                  \
      if (0) { char* temp = (char *)&variable; temp++; } \
   } while (0)

/*!
 * Throw an error assertion from within any C++ source code.  The
 * macro argument may be any standard ostream expression.  The file and
 * line number of the abort are also printed.
 */
#define TBOX_ERROR(X)                                           \
   do {                                                         \
      std::ostringstream tboxos;                                \
      tboxos << X << std::ends;                                 \
      tbox::Utilities::abort(tboxos.str(), __FILE__, __LINE__); \
   } while (0)

/*!
 * Print a warning without exit.  Print file and line number of the warning.
 */
#define TBOX_WARNING(X)                        \
   do {                                        \
      std::ostringstream tboxos;               \
      tboxos << X << std::ends;                \
      tbox::Logger::getInstance()->logWarning( \
         tboxos.str(), __FILE__, __LINE__);    \
   } while (0)

/*!
 * Print a debug without exit.  Print file and line number of the debug.
 */
#define TBOX_DEBUG(X)                        \
   do {                                      \
      std::ostringstream tboxos;             \
      tboxos << X << std::ends;              \
      tbox::Logger::getInstance()->logDebug( \
         tboxos.str(), __FILE__, __LINE__);  \
   } while (0)

/*!
 * Throw an error assertion from within any C++ source code if the
 * given expression is not true.  This is a parallel-friendly version
 * of assert.
 * The file and line number of the abort are also printed.
 */
#ifdef DEBUG_CHECK_ASSERTIONS

#define TBOX_ASSERT(EXP)                                           \
   do {                                                            \
      if (!(EXP)) {                                                \
         std::ostringstream tboxos;                                \
         tboxos << "Failed assertion: " << # EXP << std::ends;     \
         tbox::Utilities::abort(tboxos.str(), __FILE__, __LINE__); \
      }                                                            \
   } while (0)
#else

/*
 * No assertion checking
 */
#define TBOX_ASSERT(EXP)

#endif

/*!
 * Throw an error assertion from within any C++ source code if the
 * given expression is not true.  This is a parallel-friendly version
 * of assert.
 * The file and line number of the abort are also printed along with the
 * supplied message.
 */
#ifdef DEBUG_CHECK_ASSERTIONS

#define TBOX_ASSERT_MSG(EXP, MSG)                                  \
   do {                                                            \
      if (!(EXP)) {                                                \
         std::ostringstream tboxos;                                \
         tboxos << "Failed assertion: " << # EXP << std::endl << # MSG << std::ends; \
         tbox::Utilities::abort(tboxos.str(), __FILE__, __LINE__); \
      }                                                            \
   } while (0)
#else

/*
 * No assertion checking
 */
#define TBOX_ASSERT_MSG(EXP, MSG)

#endif

/*!
 * Throw an error assertion from within any C++ source code if the
 * given expression is not true.  This version is used for assertions
 * that are useful checking internal SAMRAI state for developers working
 * on SAMRAI.  User level assertions should use TBOX_ASSERT.
 *
 * This is a parallel-friendly version of assert.  The file and line
 * number of the abort are also printed.
 */
#ifdef DEBUG_CHECK_DEV_ASSERTIONS

#define TBOX_DEV_ASSERT(EXP)                                       \
   do {                                                            \
      if (!(EXP)) {                                                \
         std::ostringstream tboxos;                                \
         tboxos << "Failed assertion: " << # EXP << std::ends;     \
         tbox::Utilities::abort(tboxos.str(), __FILE__, __LINE__); \
      }                                                            \
   } while (0)
#else

/*
 * No SAMRAI internal developer assertion checking
 */
#define TBOX_DEV_ASSERT(EXP)

#endif

#define TBOX_DIM_ASSERT_CHECK_DIM(dim)  \
   TBOX_DIM_ASSERT(                     \
   (dim).isValid()                      \
   )

#define TBOX_DIM_ASSERT_CHECK_DIM_ALLOW_UNINITIALIZED(dim)  \
   TBOX_DIM_ASSERT(                                         \
   !(dim).isInitialized || (dim).isValid()                  \
   )

#define TBOX_DIM_ASSERT_CHECK_ARGS1(arg1)  \
   TBOX_DIM_ASSERT(                        \
   (arg1).getDim().isValid()               \
   )

#define TBOX_DIM_ASSERT_CHECK_ARGS2(arg1, arg2)  \
   TBOX_DIM_ASSERT(                              \
   (arg1).getDim().isValid() &&                  \
   ((arg1).getDim() == (arg2).getDim())          \
   )

#define TBOX_DIM_ASSERT_CHECK_ARGS3(arg1, arg2, arg3)   \
   TBOX_DIM_ASSERT(                                     \
   (arg1).getDim().isValid() &&                         \
   ((arg1).getDim() == (arg2).getDim()) &&              \
   ((arg1).getDim() == (arg3).getDim())                 \
   )

#define TBOX_DIM_ASSERT_CHECK_ARGS4(arg1, arg2, arg3, arg4)  \
   TBOX_DIM_ASSERT(                                          \
   (arg1).getDim().isValid() &&                              \
   ((arg1).getDim() == (arg2).getDim()) &&                   \
   ((arg1).getDim() == (arg3).getDim()) &&                   \
   ((arg1).getDim() == (arg4).getDim())                      \
   )

#define TBOX_DIM_ASSERT_CHECK_ARGS5(arg1, arg2, arg3, arg4, arg5)  \
   TBOX_DIM_ASSERT(                                                \
   (arg1).getDim().isValid() &&                                    \
   ((arg1).getDim() == (arg2).getDim()) &&                         \
   ((arg1).getDim() == (arg3).getDim()) &&                         \
   ((arg1).getDim() == (arg4).getDim()) &&                         \
   ((arg1).getDim() == (arg5).getDim())                            \
   )

#define TBOX_DIM_ASSERT_CHECK_ARGS6(arg1, arg2, arg3, arg4, arg5, arg6)  \
   TBOX_DIM_ASSERT(                                                      \
   (arg1).getDim().isValid() &&                                          \
   ((arg1).getDim() == (arg2).getDim()) &&                               \
   ((arg1).getDim() == (arg3).getDim()) &&                               \
   ((arg1).getDim() == (arg4).getDim()) &&                               \
   ((arg1).getDim() == (arg5).getDim()) &&                               \
   ((arg1).getDim() == (arg6).getDim())                                  \
   )

#define TBOX_DIM_ASSERT_CHECK_ARGS7(arg1, arg2, arg3, arg4, arg5, arg6, arg7) \
   TBOX_DIM_ASSERT(                                                           \
   (arg1).getDim().isValid() &&                                               \
   ((arg1).getDim() == (arg2).getDim()) &&                                    \
   ((arg1).getDim() == (arg3).getDim()) &&                                    \
   ((arg1).getDim() == (arg4).getDim()) &&                                    \
   ((arg1).getDim() == (arg5).getDim()) &&                                    \
   ((arg1).getDim() == (arg6).getDim()) &&                                    \
   ((arg1).getDim() == (arg7).getDim())                                       \
   )

#define TBOX_DIM_ASSERT_CHECK_DIM_ARGS1(dim, arg1)  \
   TBOX_DIM_ASSERT(                                 \
   (dim).isValid() &&                               \
   ((dim) == (arg1).getDim())                       \
   )                                                \

#define TBOX_DIM_ASSERT_CHECK_DIM_ARGS2(dim, arg1, arg2)  \
   TBOX_DIM_ASSERT(                                       \
   (dim).isValid() &&                                     \
   ((dim) == (arg1).getDim()) &&                          \
   ((dim) == (arg2).getDim())                             \
   )

#define TBOX_DIM_ASSERT_CHECK_DIM_ARGS3(dim, arg1, arg2, arg3)  \
   TBOX_DIM_ASSERT(                                             \
   (dim).isValid() &&                                           \
   ((dim) == (arg1).getDim()) &&                                \
   ((dim) == (arg2).getDim()) &&                                \
   ((dim) == (arg3).getDim())                                   \
   )

#define TBOX_DIM_ASSERT_CHECK_DIM_ARGS4(dim, arg1, arg2, arg3, arg4)  \
   TBOX_DIM_ASSERT(                                                   \
   (dim).isValid() &&                                                 \
   ((dim) == (arg1).getDim()) &&                                      \
   ((dim) == (arg2).getDim()) &&                                      \
   ((dim) == (arg3).getDim()) &&                                      \
   ((dim) == (arg4).getDim())                                         \
   )

#define TBOX_DIM_ASSERT_CHECK_DIM_ARGS5(dim, arg1, arg2, arg3, arg4, arg5)  \
   TBOX_DIM_ASSERT(                                                         \
   (dim).isValid() &&                                                       \
   ((dim) == (arg1).getDim()) &&                                            \
   ((dim) == (arg2).getDim()) &&                                            \
   ((dim) == (arg3).getDim()) &&                                            \
   ((dim) == (arg4).getDim()) &&                                            \
   ((dim) == (arg5).getDim())                                               \
   )

#define TBOX_DIM_ASSERT_CHECK_DIM_ARGS6(dim, arg1, arg2, arg3, arg4, arg5, arg6) \
   TBOX_DIM_ASSERT(                                                  \
   (dim).isValid() &&                                                \
   ((dim) == (arg1).getDim()) &&                                     \
   ((dim) == (arg2).getDim()) &&                                     \
   ((dim) == (arg3).getDim()) &&                                     \
   ((dim) == (arg4).getDim()) &&                                     \
   ((dim) == (arg5).getDim()) &&                                     \
   ((dim) == (arg6).getDim())                                        \
   )

#define TBOX_DIM_ASSERT_CHECK_DIM_ARGS7(dim,  \
                                        arg1, \
                                        arg2, \
                                        arg3, \
                                        arg4, \
                                        arg5, \
                                        arg6, \
                                        arg7) \
   TBOX_DIM_ASSERT(                           \
   (dim).isValid() &&                         \
   ((dim) == (arg1).getDim()) &&              \
   ((dim) == (arg2).getDim()) &&              \
   ((dim) == (arg3).getDim()) &&              \
   ((dim) == (arg4).getDim()) &&              \
   ((dim) == (arg5).getDim()) &&              \
   ((dim) == (arg6).getDim()) &&              \
   ((dim) == (arg7).getDim())                 \
   )

#define TBOX_DIM_ASSERT_CHECK_DIM_ARGS8(dim,  \
                                        arg1, \
                                        arg2, \
                                        arg3, \
                                        arg4, \
                                        arg5, \
                                        arg6, \
                                        arg7, \
                                        arg8) \
   TBOX_DIM_ASSERT(                           \
   (dim).isValid() &&                         \
   ((dim) == (arg1).getDim()) &&              \
   ((dim) == (arg2).getDim()) &&              \
   ((dim) == (arg3).getDim()) &&              \
   ((dim) == (arg4).getDim()) &&              \
   ((dim) == (arg5).getDim()) &&              \
   ((dim) == (arg6).getDim()) &&              \
   ((dim) == (arg7).getDim()) &&              \
   ((dim) == (arg8).getDim())                 \
   )

#ifdef DEBUG_CHECK_DIM_ASSERTIONS

#define TBOX_DIM_ASSERT(EXP)                                             \
   do {                                                                  \
      if (!(EXP)) {                                                      \
         std::ostringstream tboxos;                                      \
         tboxos << "Failed dimension assertion: " << # EXP << std::ends; \
         tbox::Utilities::abort(tboxos.str(), __FILE__, __LINE__);       \
      }                                                                  \
   } while (0)

#else

/*
 * No dimensional assertion checking
 */
#define TBOX_DIM_ASSERT(EXP)

#endif

/**
 * Throw an error assertion from within any C++ source code.  This is
 * is similar to TBOX_ERROR(), but is designed to be invoked after a
 * call to a PETSc library function.  In other words, it acts similarly
 * to the PETSc CHKERRQ(ierr) macro.
 */
#ifdef HAVE_PETSC

/*
 * In the following, "CHKERRCONTINUE(ierr);" will cause PETSc to print out
 * a stack trace that led to the error; this may be useful for debugging.
 */

#define PETSC_SAMRAI_ERROR(ierr)                                   \
   do {                                                            \
      if (ierr) {                                                  \
         std::ostringstream tboxos;                                \
         tbox::Utilities::abort(tboxos.str(), __FILE__, __LINE__); \
      }                                                            \
   } while (0)
#endif

/*!
 * Utilities is a Singleton class containing basic routines for error
 * reporting, file manipulations, etc.
 */

struct Utilities {
   /*!
    * Create the directory specified by the path string.  Permissions are set
    * by default to rwx by user.  The intermediate directories in the
    * path are created if they do not already exist.  When
    * only_node_zero_creates is true, only node zero creates the
    * directories.  Otherwise, all nodes create the directories.
    */
   static void
   recursiveMkdir(
      const std::string& path,
      mode_t mode = (S_IRUSR | S_IWUSR | S_IXUSR),
      bool only_node_zero_creates = true);

   /*!
    * Rename a file from old file name to new file name.
    */
   static void
   renameFile(
      const std::string& old_filename,
      const std::string& new_filename)
   {
      TBOX_ASSERT(!old_filename.empty());
      TBOX_ASSERT(!new_filename.empty());
      rename(old_filename.c_str(), new_filename.c_str());
   }

   /*!
    * Convert an integer to a string.
    *
    * The returned string is padded with zeros as needed so that it
    * contains at least the number of characters indicated by the
    * minimum width argument.  When the number is positive, the
    * string is padded on the left. When the number is negative,
    * the '-' sign appears first, followed by the integer value
    * padded on the left with zeros.  For example, the statement
    * intToString(12, 5) returns "00012" and the statement
    * intToString(-12, 5) returns "-0012".
    */
   static std::string
   intToString(
      int num,
      int min_width = 1);

   /*!
    * Convert common integer values to strings.
    *
    * These are simply wrappers around intToString that ensure the
    * same width is uniformally used when converting to string
    * representations.
    */
   static std::string
   nodeToString(
      int num)
   {
      return intToString(num, s_node_width);
   }
   static std::string
   processorToString(
      int num)
   {
      return intToString(num, s_processor_width);
   }
   static std::string
   patchToString(
      int num)
   {
      return intToString(num, s_patch_width);
   }
   static std::string
   levelToString(
      int num)
   {
      return intToString(num, s_level_width);
   }
   static std::string
   blockToString(
      int num)
   {
      return intToString(num, s_block_width);
   }

   /*!
    * Aborts the run after printing an error message with file and
    * linenumber information.
    */
   static void
   abort(
      const std::string& message,
      const std::string& filename,
      const int line);

private:
   /*
    * Sizes for converting integers to fixed width strings
    * for things like filenames etc.
    */
   static const int s_node_width = 7;
   static const int s_processor_width = 7;
   static const int s_patch_width = 7;
   static const int s_level_width = 4;
   static const int s_block_width = 7;

};

}
}

#endif
