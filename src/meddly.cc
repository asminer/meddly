
// $Id$

/*
    Implementation of the "global" functions given in
    meddly.h and meddly_expert.h
*/

#include "defines.h"
#include "revision.h"


namespace MEDDLY {
  // "global" variables

  settings meddlySettings;

  statistics meddlyStats;

  expert_compute_manager* ECM = 0;

};

// helpers

void initStats(MEDDLY::statistics &s)
{
}

//----------------------------------------------------------------------
// front end

void MEDDLY::initialize(settings s)
{
  if (ECM) throw error(error::ALREADY_INITIALIZED);
  meddlySettings = s;
  initStats(meddlyStats);

  ECM = new expert_compute_manager(meddlySettings);
}

void MEDDLY::cleanup()
{
  if (0==ECM) throw error(error::UNINITIALIZED);
  delete ECM;
  ECM = 0;
}

const MEDDLY::settings& MEDDLY::getLibrarySettings()
{
  return meddlySettings;
}

const MEDDLY::statistics& MEDDLY::getLibraryStats()
{
  return meddlyStats;
}

const char* MEDDLY::getLibraryInfo(int what)
{
  static char* title = 0;
  switch (what) {
    case 0:
      if (!title) {
        title = new char[80];
        if (REVISION_NUMBER) {
          snprintf(title, 80, 
            "%s version %s.%d (32-bit and 64-bit compatible)", 
            PACKAGE_NAME, VERSION, REVISION_NUMBER
          );
        } else {
          snprintf(title, 80, 
            "%s version %s (32-bit and 64-bit compatible)", 
            PACKAGE_NAME, VERSION
          );
        }
      }
      return title;

    case 1:
      return "Copyright (C) 2009, Iowa State University Research Foundation, Inc.";

    case 2:
      return "Released under the GNU Lesser General Public License, version 3";
 
    case 3:
      return "http://meddly.sourceforge.net/";

    case 4:
      return "Data Structures and operations available:\n\
(1) MDDs: Union, Intersection, Difference.\n\
(2) Matrix Diagrams (MXDs): Union, Intersection, Difference.\n\
(3) Multi-Terminal MDDs (MTMDDs) with integer or real terminals:\n\
    Arithmetic: Plus, Minus, Multiply, Divide, Min, Max.\n\
    Logical: <, <=, >, >=, ==, !=.\n\
    Conversion to and from MDDs.\n\
(4) Multi-Terminal MXDs (MTMXDs) with integer or real terminals:\n\
    Arithmetic: Plus, Minus, Multiply, Divide, Min, Max.\n\
    Logical: <, <=, >, >=, ==, !=.\n\
    Conversion to and from MXDs.\n\
";
  }
  return 0;
}

MEDDLY::domain* MEDDLY::createDomain()
{
  if (0==ECM) throw error(error::UNINITIALIZED);
  return new expert_domain();
}



MEDDLY::compute_manager* MEDDLY::getComputeManager()
{
  if (0==ECM) throw error(error::UNINITIALIZED);
  return ECM;
}



