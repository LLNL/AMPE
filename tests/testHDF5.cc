// simply test HDF5 is working properly by writing an array of data into
// a file and then reading it.
#include "hdf5.h"
#include <iostream>

#define FILE "dset.h5"

#define NX 4
#define NY 6

int main()
{
   /* Initialize a dataset. */
   int dset_in[NX][NY];
   for (int i = 0; i < NX; i++)
      for (int j = 0; j < NY; j++)
         dset_in[i][j] = i * NY + j + 1;

   /* Open a new file */
   hid_t file_id = H5Fcreate(FILE, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
   if (file_id < 0) {
      std::cerr << "H5Fopen failed." << std::endl;
      return 1;
   }

   hsize_t dims[2] = {NX, NY};
   hid_t dataspace_id = H5Screate_simple(2, dims, NULL);

   /* Open a new dataset */
   hid_t dataset_id = H5Dcreate2(file_id, "/dset", H5T_STD_I32BE, dataspace_id,
                                 H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
   if (dataset_id < 0) {
      std::cerr << "H5Dopen2 failed." << std::endl;
      return 1;
   }

   /* Write dataset */
   herr_t status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
                            H5P_DEFAULT, dset_in);
   if (status < 0) {
      std::cerr << "H5Dwrite failed." << std::endl;
      return 1;
   }

   int dset_out[NX][NY];
   status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                    dset_out);
   if (status < 0) {
      std::cerr << "H5Dread failed." << std::endl;
      return 1;
   }

   /* analyse data read */
   int count = 0;
   for (int i = 0; i < NX; i++)
      for (int j = 0; j < NY; j++) {
         std::cout << "in: " << dset_in[i][j] << ", out: " << dset_out[i][j]
                   << std::endl;
         if (std::abs(dset_out[i][j] - dset_in[i][j]) > 0) {
            std::cout << "Difference larger than tol!" << std::endl;
            count++;
         }
      }
   if (count > 0) return 1;

   /* Close dataset */
   status = H5Dclose(dataset_id);
   if (status < 0) {
      std::cerr << "H5Dclose failed." << std::endl;
      return 1;
   }

   /* Close file */
   status = H5Fclose(file_id);
   if (status < 0) {
      std::cerr << "H5Fclose failed." << std::endl;
      return 1;
   }

   std::cout << "TEST successful!" << std::endl;

   return 0;
}
