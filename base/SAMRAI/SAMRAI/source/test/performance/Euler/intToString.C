/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2016 Lawrence Livermore National Security, LLC
 * Description:   $Description
 *
 ************************************************************************/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "intToString.h"

#include <string>

#include <sstream>
#include <iomanip>

using namespace std;

string intToString(
   int i,
   int min_length)
{
#if 1
   ostringstream co;
   co << setw(min_length) << setfill('0') << i;
   return string(co.str());

#else
   char f[50];
   sprintf(f, "%%%dd", min_length);
   char s[50];
   sprintf(s, f, i);
   return string(s);

#endif
}
