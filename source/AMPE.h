#ifndef included_AMPE
#define included_AMPE

#include "SAMRAI/tbox/TimerManager.h"

#include <cstring>

class PFModel;

class AMPE
{
public:
   AMPE(MPI_Comm comm);

   ~AMPE();

   std::string gitCommitID()const
   {
      #define xstr2(x) #x
      #define xstr(x) xstr2(x)

      return xstr(GITVERSION);
   }

   void initialize(const std::string input_filename,
                   const std::string restart_read_dirname="",
                   const int restore_num=-1);

   void run();

private:

   PFModel* d_pfm;

   SAMRAI::tbox::TimerManager* d_time_man;
};

#endif


