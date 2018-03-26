/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   SampleIndexData example demonstrating IndexData type.
 *
 ************************************************************************/

#include "SampleIndexData.h"

#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/pdat/CellIndex.h"
#include "SAMRAI/tbox/Dimension.h"

#include <iostream>

/*
 *************************************************************************
 *
 * Constructor providing cell index.
 *
 *************************************************************************
 */

using namespace SAMRAI;

SampleIndexData::SampleIndexData(
   const tbox::Dimension& dim):
   d_index(dim),
   d_dummy_int(0)
{
}

SampleIndexData::SampleIndexData(
   const hier::Index& ic):
   d_index(ic),
   d_dummy_int(0)
{
}

/*
 *************************************************************************
 *
 * Copy Constructor
 *
 *************************************************************************
 */

SampleIndexData::SampleIndexData(
   const SampleIndexData& data):
   d_index(data.d_index),
   d_dummy_int(data.d_dummy_int)
{

}

/*
 *************************************************************************
 *
 * Assignment operator
 *
 *************************************************************************
 */

SampleIndexData& SampleIndexData::operator = (
   const SampleIndexData& data)
{
   d_index = data.d_index;
   d_dummy_int = data.d_dummy_int;
   return *this;
}

/*
 *************************************************************************
 *
 * Destructor
 *
 *************************************************************************
 */

SampleIndexData::~SampleIndexData()
{
}

/*
 *************************************************************************
 *
 * Set dummy int data
 *
 *************************************************************************
 */
void SampleIndexData::setInt(
   const int dummy)
{
   d_dummy_int = dummy;
}

void SampleIndexData::setIndex(
   const hier::Index& index)
{
   d_index = index;
}

/*
 *************************************************************************
 *
 *  Return dummy int data
 *
 *************************************************************************
 */
int SampleIndexData::getInt() const
{
   return d_dummy_int;
}

/*
 *************************************************************************
 *
 *  Return index
 *
 *************************************************************************
 */
const hier::Index& SampleIndexData::getIndex() const
{
   return d_index;
}

/*
 *************************************************************************
 *
 *  Print contents
 *
 *************************************************************************
 */
void SampleIndexData::printClassData(
   std::ostream& os) const
{
   os << "printing data";
}

/*
 *************************************************************************
 *
 * The copySourceItem() method allows SampleIndexData to be a templated
 * data type for IndexData - i.e. IndexData<SampleIndexData>.
 *
 *************************************************************************
 */
void SampleIndexData::copySourceItem(
   hier::Index& index,
   SampleIndexData& src_item)
{
   d_index = index;
   d_dummy_int = src_item.d_dummy_int;
}

/*
 *************************************************************************
 *
 * The getDataStreamSize(), packStream(), and unpackStream() methods
 * are required to template SampleIndexData as IndexData type - i.e.
 * IndexData<SampleIndexData>.  They are used to communicate SampleIndexData,
 * specifying how many bytes will be packed during the "packStream()"
 * method.
 *
 *************************************************************************
 */

size_t SampleIndexData::getDataStreamSize()
{
   /*
    * #bytes =
    *   d_index           (int[NDIM]) +
    *   d_dummy_int       (int)
    */

   int dim = d_index.getDim().getValue();
   size_t bytes = (dim + 1) * tbox::MessageStream::getSizeof<int>();

   return bytes;
}

void SampleIndexData::packStream(
   tbox::MessageStream& stream)
{
   int counter = 0;
   int dim = d_index.getDim().getValue();
   int ibuffer[dim + 1];
   for (int i = 0; i < dim; i++) {
      ibuffer[i] = d_index(i);
      counter++;
   }
   ibuffer[counter] = d_dummy_int;
   stream.pack(ibuffer, dim + 1);

}

void SampleIndexData::unpackStream(
   tbox::MessageStream& stream,
   const hier::IntVector& offset)
{
   int counter = 0;
   int dim = d_index.getDim().getValue();
   int ibuffer[dim + 1];
   stream.unpack(ibuffer, dim);
   //pdat::CellIndex index;
   hier::Index* index = new pdat::CellIndex(d_index.getDim());
   for (int i = 0; i < dim; i++) {
      (*index)[i] = ibuffer[i];
      counter++;
   }
   *index += offset;
   d_index = *index;

   d_dummy_int = ibuffer[counter];

}

/*
 *************************************************************************
 *
 * The putToDatabase() and getFromDatabase() methods
 * are required to template SampleIndexData as IndexData type - i.e.
 * IndexData<SampleIndexData>.  They are used to write/read SampleIndexData,
 * data to/from the restart database.
 *
 *************************************************************************
 */

void SampleIndexData::putUnregisteredToDatabase(
   const boost::shared_ptr<tbox::Database>& database) const
{

   int counter = 0;
   int dim = d_index.getDim().getValue();
   int ibuffer[dim + 1];
   for (int i = 0; i < dim; i++) {
      ibuffer[i] = d_index(i);
      counter++;
   }
   ibuffer[counter] = d_dummy_int;

   database->putIntegerArray("ibuffer", ibuffer, dim + 1);

}

void SampleIndexData::getFromDatabase(
   const boost::shared_ptr<tbox::Database>& database)
{
   int dim = d_index.getDim().getValue();
   int ibuffer[dim + 1];
   database->getIntegerArray("ibuffer", ibuffer, dim + 1);
   hier::Index* index = new pdat::CellIndex(d_index.getDim());
   for (int i = 0; i < dim; i++) {
      (*index)(i) = ibuffer[i];
   }
   d_index = *index;
   d_dummy_int = ibuffer[dim];
}

/*
 *****************************************************************
 *
 *  Templates used for SampleIndexData
 *
 *****************************************************************
 */

//#include "SampleIndexData.h"
//#include "SAMRAI/tbox/Array.C"
//#include "SAMRAI/pdat/IndexData.C"
//#include "SAMRAI/pdat/IndexDataFactory.C"
//#include "SAMRAI/pdat/IndexVariable.C"
//#include "SAMRAI/pdat/CellGeometry.h"
//
//namespace SAMRAI {
//
//template class pdat::SparseData<SampleIndexData, pdat::CellGeometry>;
//template class pdat::SparseDataFactory<SampleIndexData, pdat::CellGeometry>;
//template class pdat::IndexData<SampleIndexData, pdat::CellGeometry>;
//template class pdat::IndexDataFactory<SampleIndexData, pdat::CellGeometry>;
//template class pdat::IndexDataNode<NDIM, SampleIndexData, pdat::CellGeometry>;
//template class pdat::IndexIterator<NDIM, SampleIndexData, pdat::CellGeometry>;
//template class pdat::IndexVariable<SampleIndexData, pdat::CellGeometry>;
//template class tbox::Array<SampleIndexData>;
//template class tbox::Array<pdat::IndexDataNode<NDIM, SampleIndexData,
//                                               pdat::CellGeometry> >;
//template class boost::shared_ptr<pdat::IndexData<NDIM, SampleIndexData,
//                                                 pdat::CellGeometry> >;
//template class boost::shared_ptr<pdat::IndexVariable<SampleIndexData,
//                                                     pdat::CellGeometry> >;
//template class boost::shared_ptr<pdat::IndexDataFactory<SampleIndexData,
//                                                        pdat::CellGeometry> >;
//
//}
