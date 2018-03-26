/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Globally unique identifier that can be locally determined.
 *
 ************************************************************************/
#ifndef included_hier_GlobalId_C
#define included_hier_GlobalId_C

#include "SAMRAI/hier/GlobalId.h"

#include <iostream>

namespace SAMRAI {
namespace hier {

/*
 ************************************************************************
 ************************************************************************
 */
GlobalId::GlobalId():
   d_owner_rank(tbox::SAMRAI_MPI::getInvalidRank()),
   d_local_id(LocalId::getInvalidId())
{
}

/*
 ************************************************************************
 ************************************************************************
 */
GlobalId::GlobalId(
   const LocalId& local_id,
   const int owner):
   d_owner_rank(owner),
   d_local_id(local_id)
{
}

/*
 ************************************************************************
 ************************************************************************
 */
GlobalId::GlobalId(
   const GlobalId& other):
   d_owner_rank(other.d_owner_rank),
   d_local_id(other.d_local_id)
{
}

/*
 ************************************************************************
 ************************************************************************
 */
GlobalId::~GlobalId()
{
}

std::ostream&
operator << (
   std::ostream& co,
   const GlobalId& r)
{
   co << r.d_owner_rank << '#' << r.d_local_id;
   return co;
}

}
}
#endif
