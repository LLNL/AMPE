/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Boundary cell struct for embedded boundary implementations
 *
 ************************************************************************/

#ifndef included_SampleIndexDataXD
#define included_SampleIndexDataXD

#include "SAMRAI/SAMRAI_config.h"

//#include "SAMRAI/pdat/CellIndex.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/tbox/MessageStream.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/IOStream.h"

#include <boost/shared_ptr.hpp>

/**
 * The SampleClass struct holds some dummy data and methods.  It's intent
 * is to indicate how a user could construct their own index data type.
 */

using namespace SAMRAI;

class SampleIndexData
{
public:
   SampleIndexData(
      const tbox::Dimension& d);
   /**
    * Copy constructor.
    */
   SampleIndexData(
      const SampleIndexData& data);

   /**
    * Constructor supplying cell index where data is defined.
    */
   SampleIndexData(
      const hier::Index& ic);
//      const pdat::CellIndex& ic);

   /**
    * The assignment operator copies the data of the argument cell.
    */
   SampleIndexData&
   operator = (
      const SampleIndexData& cell);

   /**
    * The destructor for SampleIndexData.
    */
   ~SampleIndexData();

   /**
    * Sets a dummy integer in this class.
    */
   void
   setInt(
      const int dummy);

   void
   setIndex(
      const hier::Index& index);

   /**
    * Returns a dummy integer in this class.
    */
   int
   getInt() const;

   /**
    * Returns the cell index where the index data is stored.
    */
   //pdat::CellIndex
   const hier::Index&
   getIndex() const;

   /**
    * Print class data representation when an unrecoverable run-time
    * assertion is thrown. Or, when desired.
    */
   void
   printClassData(
      std::ostream& os) const;

   /**
    * The copySourceItem() method allows SampleIndexData to be a templated
    * data type for IndexData - i.e. IndexData<SampleIndexData>.  In
    * addition to this method, the other methods that must be defined are
    * getDataStreamSize(), packStream(), unpackStream() for communication,
    * putToDatabase(), getFromDatabase for restart.  These are described
    * below.
    */
   void
   copySourceItem(
      hier::Index& index,
      SampleIndexData& src_item);

   /**
    * The following functions enable parallel communication with SampleIndexDatas.
    * They are used in SAMRAI communication infrastructure to
    * specify the number of bytes of data stored in each SampleIndexData object,
    * and to pack and unpack the data to the specified stream.
    */
   size_t
   getDataStreamSize();
   void
   packStream(
      tbox::MessageStream& stream);
   void
   unpackStream(
      tbox::MessageStream& stream,
      const hier::IntVector& offset);

   /**
    * These functions are used to read/write SampleIndexData data to/from
    * restart.
    */
   void
   getFromDatabase(
      const boost::shared_ptr<tbox::Database>& database);
   void
   putUnregisteredToDatabase(
      const boost::shared_ptr<tbox::Database>& database) const;

private:
   /*
    * Cell index where SampleIndexData is defined.
    */
   //pdat::CellIndex d_index;
   hier::Index d_index;

   /*
    * Dummy int data
    */
   int d_dummy_int;

   // ADD ANY OTHER DATA HERE
};
#endif
