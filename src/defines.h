  "06/15/2026"
  "06/15/2026"
  "06/15/2026"
  "06/15/2026"
  "06/15/2026"
  "06/15/2026"
  "06/15/2026"
  "06/15/2026"
  "06/15/2026"
  "06/15/2026"
  "06/15/2026"
  "06/15/2026"
  "06/15/2026"
  "06/15/2026"
  "06/15/2026"
  "06/15/2026"
  "06/15/2026"
  "06/15/2026"
  "06/15/2026"
  "06/15/2026"
  "06/15/2026"
  "06/15/2026"
  "06/15/2026"
  "06/15/2026"
  "06/15/2026"
  "06/15/2026"
  "06/15/2026"
  "06/15/2026"
  "06/15/2026"
  "06/15/2026"
  "06/15/2026"
  "06/15/2026"
  "06/15/2026"
  "06/15/2026"
  "06/15/2026"
// ^ This line is important for the release script, don't remove it

// Reference counts on for now
#define REFCOUNTS_ON

/*

    There are now two modes of code generation:
    "DEVELOPMENT_CODE" and "RELEASE_CODE".

    If "DEVELOPMENT_CODE" is defined (usually done in the makefile)
    then debugging macros and assertions will be turned on.  Otherwise
    we assume that we have "RELEASE_CODE" and they are turned off.

*/

#ifdef DEVELOPMENT_CODE
// ***********************************************************************
// *                                                                     *
// *                      Development code settings                      *
// *                                                                     *
// ***********************************************************************

#include <cassert>
#include <iostream>
#include <cassert>

#define smart_cast  dynamic_cast
#define MEDDLY_DCASSERT(X)  assert(X)
#define MEDDLY_CHECK_RANGE(L, X, U)  \
    do { \
        if ((X) < (L) || (X) >= (U) ) { \
            std::cerr << "Check range at " << __FILE__ \
                      << " line " << __LINE__ \
                      << " failed:\n    min: " << (L) \
                      << "\n    val: " << (X) \
                      << "\n    max: " << (U) << '\n'; \
            exit(1); \
        } \
    } while (0)


// ***********************************************************************
#else
// ***********************************************************************
// *                                                                     *
// *                        Release code settings                        *
// *                                                                     *
// ***********************************************************************

#define smart_cast  static_cast
#define MEDDLY_DCASSERT(X)
#define MEDDLY_CHECK_RANGE(L, X, U)

// ***********************************************************************
#endif


#include <limits>

namespace MEDDLY {
    // ***********************************************************************
    // *                               typedefs                              *
    // ***********************************************************************

    /** Handles for nodes.
        This should be either int or long, and effectively limits
        the number of possible nodes per forest.
        As an int, we get 2^32-1 possible nodes per forest,
        which should be enough for most applications.
        As a long on a 64-bit machine, we get 2^64-1 possible nodes
        per forest, at the expense of nearly doubling the memory used.
        This also specifies the incoming count range for each node.
    */
    typedef int  node_handle;
    // typedef long node_handle;

    /** Handles for relation nodes.
        TBD: can we just use node_handle everywhere?
     */
    typedef node_handle rel_node_handle;

    /** Node addresses.
        This is used for internal storage of a node,
        and should probably not be changed.
        The typedef is given simply to clarify the code
        (hopefully :^)
    */
    typedef unsigned long node_address;


    // ***********************************************************************
    // *                             Handy macros                            *
    // ***********************************************************************

    template<typename T> inline T Inf() {
        return std::numeric_limits<T>::max();
    }

    /// Standard MAX "macro".
    template <class T> inline T MAX(T X,T Y) {
        return ((X>Y)?X:Y);
    }

    /// Standard MIN "macro".
    template <class T> inline T MIN(T X,T Y) {
        return ((X<Y)?X:Y);
    }

    /// Standard ABS "macro".
    template <class T> inline T ABS(T X) {
        return ((X<0)?(-X):(X));
    }

    /// SWAP "macro".
    template <class T> inline void SWAP(T &x, T &y) {
        T tmp=x;
        x=y;
        y=tmp;
    }

    /// POSITIVE "macro".
    template <class T> inline bool POSITIVE(T X) {
        return (X>0);
    }

    /// Update a maximum
    template <class T> inline void UPDATEMAX(T &X, T Y) {
        if (Y > X) X = Y;
    }

    // Number of digits
    template <class T>
    inline unsigned digits(T a) {
        unsigned d;
        for (d=1; a; d++) { a /= 10; }
        return d;
    }

    /// Determine if level k1 is above k2.  Works for primed & unprimed.
    inline bool isLevelAbove(int k1, int k2) {
        if (ABS(k1) > ABS(k2)) return true;
        if (ABS(k2) > ABS(k1)) return false;
        return k1 > k2;
    }

    // ***********************************************************************
    // *                         Sanity check  macros                        *
    // ***********************************************************************
    inline void FAIL(const char* fn, unsigned ln, const char* why = nullptr)
    {
#ifdef DEVELOPMENT_CODE
        if (why) {
            std::cerr   << "\nFATAL: " << why << " at " << fn << " line "
                        << ln << "\n";
        } else {
            std::cerr   << "\nFATAL error at " << fn << " line "
                        << ln << "\n";
        }
        exit(1);
#endif
    }

} // namespace MEDDLY


#endif // #include guard

