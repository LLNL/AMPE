/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Registry of shutdown routines to be called at program exit
 *
 ************************************************************************/

#include "SAMRAI/tbox/StartupShutdownManager.h"

#include "SAMRAI/tbox/Utilities.h"

#include <cstdlib>

namespace SAMRAI {
namespace tbox {

bool StartupShutdownManager::s_singleton_initialized = false;

bool StartupShutdownManager::s_in_initialize = false;
bool StartupShutdownManager::s_in_startup = false;
bool StartupShutdownManager::s_in_shutdown = false;
bool StartupShutdownManager::s_in_finalize = false;

bool StartupShutdownManager::s_initialized = false;
bool StartupShutdownManager::s_startuped = false;
bool StartupShutdownManager::s_shutdowned = false;
bool StartupShutdownManager::s_finalized = false;

StartupShutdownManager::ListElement *
StartupShutdownManager::s_manager_list[s_number_of_priorities];
StartupShutdownManager::ListElement *
StartupShutdownManager::s_manager_list_last[s_number_of_priorities];
int StartupShutdownManager::s_num_manager_items[s_number_of_priorities];

void
StartupShutdownManager::registerHandler(
   AbstractHandler* handler)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(handler);

   // Don't allow registering handlers when we are looping and the
   // handler needs to be called in that loop.  This would create the
   // possibility that a handler is registered that needs to get
   // called.
   //
   // SGS Ideally this would not be needed and maybe with some
   // additional work this could be made more clean.
   assert(!(s_in_initialize && handler->hasInitialize()));
   assert(!(s_in_startup && handler->hasStartup()));
   assert(!(s_in_shutdown && handler->hasShutdown()));
   assert(!s_in_finalize);
#endif

   if (!s_singleton_initialized) {
      setupSingleton();
   }

   ListElement* item = new ListElement;
   item->handler = handler;

   unsigned char priority = handler->getPriority();

   item->next = (ListElement *)NULL;
   if (s_num_manager_items[priority] == 0) {
      s_manager_list[priority] = item;
   } else {
      s_manager_list_last[priority]->next = item;
   }
   s_manager_list_last[priority] = item;
   s_num_manager_items[priority]++;
}

void
StartupShutdownManager::initialize()
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!s_initialized);
#endif

   s_initialized = true;
   // only shutdown if something was registered
   if (s_singleton_initialized) {

      s_in_initialize = true;

      for (int priority = 0;
           priority < s_number_of_priorities;
           ++priority) {
         ListElement* item = s_manager_list[priority];
         while (item) {
            if (item->handler) {
               item->handler->initialize();
            }
            item = item->next;
         }
      }

      s_in_initialize = false;
   }
}

void
StartupShutdownManager::startup()
{
#ifdef DEBUG_CHECK_ASSERTIONS
   // If this is thrown you need to make sure SAMRAIManger::initialize
   // is called before startup.
   assert(s_initialized);
   assert(!s_startuped);
#endif

   s_startuped = true;

   // only shutdown if something was registered
   if (s_singleton_initialized) {
      s_in_startup = true;

      for (int priority = 0;
           priority < s_number_of_priorities;
           ++priority) {
         ListElement* item = s_manager_list[priority];
         while (item) {
            if (item->handler) {
               item->handler->startup();
            }
            item = item->next;
         }
      }

      s_in_startup = false;
   }

   s_shutdowned = false;
}

void
StartupShutdownManager::shutdown()
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(s_initialized);
   assert(s_startuped);
   assert(!s_shutdowned);
#endif

   s_shutdowned = true;

   // only shutdown if something was registered
   if (s_singleton_initialized) {
      s_in_shutdown = true;

      for (int priority = s_number_of_priorities - 1;
           priority > -1;
           priority--) {
         ListElement* item = s_manager_list[priority];
         while (item) {
            if (item->handler) {
               item->handler->shutdown();
            }
            item = item->next;
         }
      }
      s_in_shutdown = false;
   }

   s_startuped = false;

}

void
StartupShutdownManager::setupSingleton()
{
   for (int priority = s_number_of_priorities - 1; priority > -1; priority--) {
      s_manager_list[priority] = (ListElement *)NULL;
      s_manager_list_last[priority] = (ListElement *)NULL;
      s_num_manager_items[priority] = 0;
   }

   s_singleton_initialized = true;
}

void
StartupShutdownManager::finalize()
{

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(s_initialized);
   assert(s_shutdowned);
   assert(!s_finalized);
#endif

   s_finalized = true;

   // only finalize if something was registered
   if (s_singleton_initialized) {
      s_in_finalize = true;

      for (int priority = s_number_of_priorities - 1;
           priority > -1;
           priority--) {
         ListElement* item = s_manager_list[priority];
         while (item) {
            if (item->handler) {
               item->handler->finalize();
            }
            item = item->next;
         }
      }

      for (int priority = 0;
           priority < s_number_of_priorities;
           ++priority) {
         ListElement* item = s_manager_list[priority];
         while (item) {
            ListElement* to_delete = item;
            item = item->next;
            delete to_delete;
         }
      }

      s_in_finalize = false;
   }

   s_initialized = false;
}

StartupShutdownManager::AbstractHandler::~AbstractHandler()
{
}

StartupShutdownManager::Handler::Handler()
{
}

StartupShutdownManager::Handler::Handler(
   void(*initialize)(),
   void(*startup)(),
   void(*shutdown)(),
   void(*finalize)(),
   unsigned char priority):
   d_initialize(initialize),
   d_startup(startup),
   d_shutdown(shutdown),
   d_finalize(finalize),
   d_priority(priority)
{
   StartupShutdownManager::registerHandler(this);
}

StartupShutdownManager::Handler::~Handler()
{
}

void
StartupShutdownManager::Handler::initialize()
{
   if (d_initialize) {
      (*d_initialize)();
   }
}

void
StartupShutdownManager::Handler::startup()
{
   if (d_startup) {
      (*d_startup)();
   }
}

void
StartupShutdownManager::Handler::shutdown()
{
   if (d_shutdown) {
      (*d_shutdown)();
   }
}

void
StartupShutdownManager::Handler::finalize()
{
   if (d_finalize) {
      (*d_finalize)();
   }
}

unsigned char
StartupShutdownManager::Handler::getPriority()
{
   return d_priority;
}

bool
StartupShutdownManager::Handler::hasInitialize()
{
   return d_initialize != 0;
}

bool
StartupShutdownManager::Handler::hasStartup()
{
   return d_startup != 0;
}

bool
StartupShutdownManager::Handler::hasShutdown()
{
   return d_shutdown != 0;
}

bool
StartupShutdownManager::Handler::hasFinalize()
{
   return d_finalize != 0;
}

StartupShutdownManager::ListElement::ListElement():
   handler(0),
   next(0)
{
}

}
}
