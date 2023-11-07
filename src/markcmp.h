
#ifndef MEDDLY_MARKCMP_H
#define MEDDLY_MARKCMP_H
#include "defines.h"
#include "forest.h"

namespace MEDDLY {
    class markcmp;
    class forest;
    class expert_forest;
};

class MEDDLY::markcmp {
public:
  virtual int compare(int i,int j,int k);
  virtual void mprint();
  virtual bool hasOmega(int i,int k);

};
#endif
