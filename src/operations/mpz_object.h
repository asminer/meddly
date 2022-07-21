
/*
    Meddly: Multi-terminal and Edge-valued Decision Diagram LibrarY.
    Copyright (C) 2009, Iowa State University Research Foundation, Inc.

    This library is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published 
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this library.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef MPZ_OBJECT
#define MPZ_OBJECT

#ifdef HAVE_LIBGMP
#include <gmp.h>

namespace MEDDLY {

class mpz_object : public ct_object {
  mpz_t value;
public:
    double rdvalue; //used for double reminder of a division
    int index;
    bool isIn;
    bool neverIn;

  mpz_object();
  mpz_object(const mpz_t &v);
  mpz_object(const mpz_object &v);

 mpz_object& operator=(const mpz_object& that) {
     // printf("Copy Operator\n" );

     if (this == &that)
        return *this;
        mpz_set(value, that.value);
        index=that.index;
        rdvalue=that.rdvalue;
        isIn=that.isIn;
        neverIn= that.neverIn;
      return *this;
 }

  virtual ~mpz_object();
  virtual opnd_type getType();

  inline void copyInto(mpz_t &x) const {
    mpz_set(x, value);
  }
  inline void copyInto(mpz_object &x) const {
    mpz_set(x.value, value);
  }
  inline void setIndex(int i) {
    index=i;
  }
  inline void setReminder(double i) {
    rdvalue=i;
  }
  inline void copyIntowithReminder(mpz_object &x)  {
    mpz_set(x.value, value);
    index=x.index;
    isIn=x.isIn;
    neverIn=x.neverIn;
    // setReminder(x.rdvalue);
    rdvalue=x.rdvalue;
  }

  inline void setValue(long i) {
    mpz_set_si(value, i);
  }

  void show(output &strm) const;
  void showwithreminder(output &strm) const;

  inline void multiply(long i) {
    mpz_mul_si(value, value, i);
  }
  inline void division(unsigned long i){
      int rivalue=mpz_cdiv_q_ui(value, value, i);
     rdvalue=rivalue/(double)i;
      // mpf_div(value,value,)
  }
  inline int compare(mpz_object x, mpz_object y){
      int qcomparision =mpz_cmp(x.value,y.value);
     if(qcomparision !=0) return qcomparision ;
     if(x.rdvalue> y.rdvalue) return 1;
     if(x.rdvalue< y.rdvalue) return -1;
     else return 0;
  }
  inline void add(const mpz_object &x) {
    mpz_add(value, value, x.value);
  }

  static void initBuffer();
  static void clearBuffer();

private:
  
  static char* buffer;  // for show().
  static int bufsize;   // for show().
  static void enlargeBuffer(int digits);
};

}

#endif
#endif

