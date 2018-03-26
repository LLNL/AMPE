/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Fixed-size message buffer used in interprocessor communication
 *
 ************************************************************************/

#ifndef included_tbox_MessageStream
#define included_tbox_MessageStream

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/Complex.h"
#include "SAMRAI/tbox/Utilities.h"

#include <cstring>
#include <iostream>
#include <vector>

namespace SAMRAI {
namespace tbox {

/*!
 * @brief Class to provide buffers for communication of data.
 *
 * MessageStream provides a message buffer that can hold data of any
 * type.  It is used by communication routines in the Schedule class.
 *
 * TODO: Because this class supports both read and write modes, it has
 * extra data and methods that don't make sense, depending on the
 * mode.  It should be rewritten as two classes, like std::cin and
 * std::cout are.  BTNG.
 *
 * @see tbox::Schedule
 */

class MessageStream
{
public:
   /*!
    * @brief Enumeration to identify if a buffer is being used to read or
    * write data.
    */
   enum StreamMode { Read, Write };

   /*!
    * @brief Create a message stream of the specified size and mode
    *
    * @param[in] bytes   Number of bytes in the stream.
    *
    * @param[in] mode    MessageStream::Read or MessageStream::Write.
    *
    * @param[in] data_to_read    Data for unpacking, should be num_bytes bytes long.
    *   This is used when mode == MessageStream::Read, ignored in write mode.
    *
    * @param[in] deep_copy Whether to make deep copy of data_to_read.
    * The default is to make a deep copy, which is safer but slower
    * than a shallow (pointer) copy.  This is used when mode ==
    * MessageStream::Read, ignored in write mode.  In shallow copy mode,
    * you cannot call growBufferAsNeeded().
    */
   MessageStream(
      const size_t bytes,
      const StreamMode mode,
      const void *data_to_read = NULL,
      bool deep_copy = true);

   /*!
    * @brief Default constructor creates a message stream with a
    * buffer that automatically grows as needed, for writing.
    */
   MessageStream();

   /*!
    * Destructor for a message stream.
    */
   ~MessageStream();

   /*!
    * @brief Static method to get amount of message stream space needed to
    * communicate data type indicated by template parameter.
    *
    * IMPORTANT:  All size information given to the message stream should
    * be based on values returned by this method.
    *
    * TODO:  Implementation should be moved out of header?  If we do this,
    * Then we need to create another implementation file to include in this
    * header.  I don't think it's worth it. RDH
    *
    * @return The number of bytes for num_items of type DATA_TYPE.
    *
    * @param[in] num_items
    */
   template<typename DATA_TYPE>
   static unsigned int getSizeof(
      unsigned int num_items = 1)
   {
      return num_items * static_cast<unsigned int>(sizeof(DATA_TYPE));
   }

   /*!
    * @brief Return a pointer to the start of the message buffer.
    */
   const void *
   getBufferStart() const
   {
      TBOX_ASSERT( d_buffer_access != NULL );
      return static_cast<const void *>(d_buffer_access);
   }

   /*!
    * @brief Return the current size of the buffer in bytes.
    */
   size_t
   getCurrentSize() const
   {
      return d_buffer_index;
   }

   /*!
    * @brief Tell a Write-mode stream to allocate more buffer
    * as needed for data.
    *
    * It is an error to use this method for a Read-mode stream.
    */
   void
   growBufferAsNeeded()
   {
      TBOX_ASSERT( d_mode == Write );
      d_grow_as_needed = true;
      return;
   }

   /*!
    * @brief Whether a Read-mode MessageStream has reached the end of
    * its data.
    */
   bool endOfData() const
   {
      TBOX_ASSERT( d_mode == Read );
      return d_buffer_index >= d_buffer_size;
   }

   /*!
    * @brief Pack a single data item into message stream.
    *
    * @param[in] data  Single item of type DATA_TYPE to be copied
    * into the stream.
    */
   template<typename DATA_TYPE>
   MessageStream&
   operator << (
      const DATA_TYPE& data)
   {
      TBOX_ASSERT(d_mode == MessageStream::Write);
      static const unsigned int nbytes =
         MessageStream::getSizeof<DATA_TYPE>(1);
      copyDataIn(static_cast<const void *>(&data), nbytes);
      return *this;
   }

   /*!
    * @brief Pack an array of data items into message stream.
    *
    * @param[in] data  Pointer to an array of data of type DATA_TYPE
    *                  to be copied into the stream.
    * @param[in] size  Number of items to pack.
    */
   template<typename DATA_TYPE>
   void
   pack(
      const DATA_TYPE* data,
      unsigned int size = 1)
   {
      TBOX_ASSERT(d_mode == MessageStream::Write);
      if (data && (size > 0)) {
         const unsigned int nbytes = MessageStream::getSizeof<DATA_TYPE>(size);
         copyDataIn(static_cast<const void *>(data), nbytes);
      }
   }

   /*!
    * @brief Unpack a single data item from message stream.
    *
    * @param[out] data  Single item of type DATA_TYPE that will be
    *                   copied from the stream.
    */
   template<typename DATA_TYPE>
   MessageStream&
   operator >> (
      DATA_TYPE& data)
   {
      TBOX_ASSERT(d_mode == MessageStream::Read);
      static const unsigned int nbytes =
         MessageStream::getSizeof<DATA_TYPE>(1);
      copyDataOut(static_cast<void *>(&data), nbytes);
      return *this;
   }

   /*!
    * @brief Unpack an array of data items from message stream.
    *
    * @param[out] data  Pointer to an array of data of type DATA_TYPE
    *                   that will receive data copied from
    *                   the stream.
    * @param[out] size  Number of items that will be copied.
    */
   template<typename DATA_TYPE>
   void
   unpack(
      DATA_TYPE * data,
      unsigned int size = 1)
   {
      TBOX_ASSERT(d_mode == MessageStream::Read);
      if (data && (size > 0)) {
         const unsigned int nbytes = MessageStream::getSizeof<DATA_TYPE>(size);
         copyDataOut(static_cast<void *>(data), nbytes);
      }
   }

   /*!
    * @brief Print out internal object data.
    *
    * @param[out] os  Output stream.
    */
   void
   printClassData(
      std::ostream& os) const;

private:

   /*!
    * @brief Copy data into the stream, advancing the stream pointer.
    *
    * @param[in]  data
    * @param[in]  num_bytes
    */
   void copyDataIn(
      const void *input_data,
      const size_t num_bytes)
      {
         if ( !d_grow_as_needed ) {
            TBOX_ASSERT(d_buffer_index + num_bytes <= d_buffer.capacity());
         }
         if ( num_bytes > 0 ) {
            d_buffer.insert( d_buffer.end(),
                             static_cast<const char*>(input_data),
                             static_cast<const char*>(input_data) + num_bytes );
            d_buffer_size = d_buffer.size();
            d_buffer_index += num_bytes;
            d_buffer_access = &d_buffer[0];
         }
         return;
      }

   /*!
    * @brief Copy data out of the stream, advancing the stream pointer.
    *
    * @param[in]  output_data
    * @param[in]  num_bytes
    */
   void copyDataOut(
      void *output_data,
      const size_t num_bytes)
      {
         TBOX_ASSERT( d_buffer_index + num_bytes <= d_buffer_size );
         memcpy(output_data, &d_buffer_access[d_buffer_index], num_bytes);
         d_buffer_index += num_bytes;
         return;
      }

   MessageStream(
      const MessageStream&);            // not implemented
   void
   operator = (
      const MessageStream&);            // not implemented

   /*!
    * @brief  Read/write mode of the stream.
    */
   const StreamMode d_mode;

   /*!
    * The buffer for the streamed data.
    */
   std::vector<char> d_buffer;

   /*!
    * @brief Pointer to either d_buffer space or, in shallow-copy Read
    * mode, external memory.
    */
   const char *d_buffer_access;

   /*!
    * @brief Number of bytes in the buffer.
    *
    * Equal to d_buffer.size() if using internal buffer.  Otherwixe,
    * equal to external buffer size.
    */
   size_t d_buffer_size;

   /*!
    * Current index into the buffer used when traversing.
    */
   size_t d_buffer_index;

   /*!
    * @brief Whether to grow buffer as needed in a Write-mode stream.
    */
   bool d_grow_as_needed;

};

}
}

#endif
