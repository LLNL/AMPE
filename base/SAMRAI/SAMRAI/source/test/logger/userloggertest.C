/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Test program to demonstrate/test a user defined logger appender
 *
 ************************************************************************/

#include "SAMRAI/SAMRAI_config.h"

#include <string>
#include <iostream>
using namespace std;

#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/Logger.h"

using namespace SAMRAI;

/*
 * Simple appender that sends log messages to a file
 */
class StreamAppender:public tbox::Logger::Appender
{

public:
   StreamAppender(
      ostream* stream) {
      d_stream = stream;
   }

   void logMessage(
      const std::string& message,
      const std::string& filename,
      const int line)
   {
      (*d_stream) << "At :" << filename << " line :" << line
                  << " message: " << message << std::endl;
   }

private:
   ostream* d_stream;
};

int main(
   int argc,
   char* argv[])
{
   int fail_count = 0;

   tbox::SAMRAI_MPI::init(&argc, &argv);
   tbox::SAMRAIManager::initialize();
   tbox::SAMRAIManager::startup();

   /*
    * Create block to force pointer deallocation.  If this is not done
    * then there will be memory leaks reported.
    */
   {

      fstream file("user.log", fstream::out);

      boost::shared_ptr<tbox::Logger::Appender> appender(
         new StreamAppender(&file));

      tbox::Logger::getInstance()->setWarningAppender(appender);
      tbox::Logger::getInstance()->setAbortAppender(appender);
      tbox::Logger::getInstance()->setDebugAppender(appender);

      /*
       * Write a test warning message.
       */
      TBOX_WARNING("Test warning");

      /*
       * Write a test debug message. Shouldn't see this since
       * Debug messages are off by default.
       */
      TBOX_DEBUG("Test debug1 : should not show up");

      tbox::Logger::getInstance()->setDebug(true);

      /*
       * Write a test debug message. Should see this
       * one since we have turned on debug messages.
       */
      TBOX_DEBUG("Test debug2 : should show up");
   }

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAIManager::finalize();
   tbox::SAMRAI_MPI::finalize();
   return fail_count;
}
