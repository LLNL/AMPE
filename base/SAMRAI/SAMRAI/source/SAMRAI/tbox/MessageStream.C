/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2016 Lawrence Livermore National Security, LLC
 * Description:   Fixed-size message buffer used in interprocessor communication
 *
 ************************************************************************/
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
   const void* data_to_read,
   bool deep_copy):
   d_mode(mode),
   d_write_buffer(),
   d_read_buffer(0),
   d_buffer_size(0),
   d_buffer_index(0),
   d_grow_as_needed(false),
   d_deep_copy_read(deep_copy)
{
   TBOX_ASSERT(num_bytes >= 1);

   if (mode == Read) {
      if (data_to_read == 0) {
         TBOX_ERROR("MessageStream::MessageStream: error:\n"
            << "No data_to_read was given to a Read-mode MessageStream.\n");
      }
      if (deep_copy) {
         d_read_buffer = new char[num_bytes];
         memcpy(const_cast<char *>(d_read_buffer), data_to_read, num_bytes);
      } else {
         d_read_buffer = static_cast<const char *>(data_to_read);
      }
      d_buffer_size = num_bytes;
   } else {
      d_write_buffer.reserve(num_bytes);
   }
}

MessageStream::MessageStream():
   d_mode(Write),
   d_write_buffer(),
   d_read_buffer(0),
   d_buffer_size(0),
   d_buffer_index(0),
   d_grow_as_needed(true),
   d_deep_copy_read(false)
{
   d_write_buffer.reserve(10);
}

MessageStream::~MessageStream()
{
   if (d_mode == Read && d_deep_copy_read) {
      delete[] d_read_buffer;
   }
   d_read_buffer = 0;
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
   os << "Pointer to buffer data = " << static_cast<const void *>(d_read_buffer) << std::endl;
}

}
}
