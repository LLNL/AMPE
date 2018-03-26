/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Fixed-size message buffer used in interprocessor communication
 *
 ************************************************************************/

#ifndef included_tbox_MessageStream_C
#define included_tbox_MessageStream_C

#include "SAMRAI/tbox/MessageStream.h"
#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace tbox {

/*
 *************************************************************************
 *
 * The constructor and destructor for MessageStream.
 *
 *************************************************************************
 */

MessageStream::MessageStream(
   const size_t num_bytes,
   const StreamMode mode,
   const void *data_to_read,
   bool deep_copy):
   d_mode(mode),
   d_buffer(),
   d_buffer_access(NULL),
   d_buffer_size(0),
   d_buffer_index(0),
   d_grow_as_needed(false)
{
   TBOX_ASSERT(num_bytes >= 1);
   d_buffer.reserve(num_bytes);

   if ( mode == Read ) {
      if ( num_bytes > 0 && data_to_read == NULL ) {
         TBOX_ERROR("MessageStream::MessageStream: error:\n"
                    <<"No data_to_read was given to a Read-mode MessageStream.\n");
      }
      if ( deep_copy ) {
         d_buffer.insert( d_buffer.end(),
                          static_cast<const char*>(data_to_read),
                          static_cast<const char*>(data_to_read)+num_bytes );
         d_buffer_access = &d_buffer[0];
      }
      else {
         d_buffer_access = static_cast<const char*>(data_to_read);
      }
      d_buffer_size = num_bytes;
   }
   return;
}

MessageStream::MessageStream()
   : d_mode(Write),
     d_buffer(),
     d_buffer_access(NULL),
     d_buffer_size(0),
     d_buffer_index(0),
     d_grow_as_needed(true)
{
   d_buffer.reserve(10);
   return;
}

MessageStream::~MessageStream()
{
   d_buffer_access = NULL;
}

/*
 *************************************************************************
 *
 * Print out class data if an assertion is thrown.
 *
 *************************************************************************
 */

void
MessageStream::printClassData(
   std::ostream& os) const
{
   os << "Maximum buffer size = " << d_buffer_size << std::endl;
   os << "Current buffer index = " << d_buffer_index << std::endl;
   os << "Pointer to buffer data = " << static_cast<const void *>(d_buffer_access) << std::endl;
}

}
}

#endif
