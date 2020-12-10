#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/MemoryDatabase.h"

#include <boost/property_tree/ptree.hpp>

void copyDatabase(std::shared_ptr<SAMRAI::tbox::Database> database,
                  boost::property_tree::ptree& new_database);

void copyDatabase(std::shared_ptr<SAMRAI::tbox::Database> database,
                  std::shared_ptr<SAMRAI::tbox::Database> new_database);
