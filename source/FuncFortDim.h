#ifndef FuncFortDim_H
#define FuncFortDim_H

// Link between C/C++ and Fortran files
//       name in             name in
//      C/C++ code            Fortran code
//      ----------            ------------
#if (NDIM == 2)
#define FORT_SINGLE_INDEX_STIFFNESS single_index_stiffness2d_
#define FORT_STORAGE_INDEX_STIFFNESS storage_index_stiffness2d_
#endif
#if (NDIM == 3)
#define FORT_SINGLE_INDEX_STIFFNESS single_index_stiffness3d_
#define FORT_STORAGE_INDEX_STIFFNESS storage_index_stiffness3d_
#endif

// Function argument list interfaces
extern "C" {
   int    FORT_SINGLE_INDEX_STIFFNESS( const int&, const int& );
   int    FORT_STORAGE_INDEX_STIFFNESS( const int&, const int& );
}

#endif

