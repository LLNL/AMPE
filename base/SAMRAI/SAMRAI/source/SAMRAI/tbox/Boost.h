/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2016 Lawrence Livermore National Security, LLC
 * Description:   SAMRAI support utilities for Boost
 *
 ************************************************************************/

#ifndef included_tbox_boost
#define included_tbox_boost

/*!
 * Suppress warning messages from boost.
 *
 * The boost implementation may generate a large number of warnings
 * with the default SAMRAI compiler options.   These macros may
 * be used to suppress these warnings.
 *
 * These are by nature compiler specific as there is no standard
 * warning suppression mechanism.
 */

#if __GNUC__ >= 5 || (__GNUC__ == 4 && (__GNUC_MINOR__ > 6))

#define BEGIN_BOOST_WARNING_SUPPRESSION \
   _Pragma("GCC diagnostic push") \
   _Pragma("GCC diagnostic ignored \"-Wconversion\"")

#define END_BOOST_WARNING_SUPPRESSION \
   _Pragma("GCC diagnostic pop")

#else

/*
 * Default to suppress nothing.
 */
#define BEGIN_BOOST_WARNING_SUPPRESSION
#define END_BOOST_WARNING_SUPPRESSION

#endif

#endif
