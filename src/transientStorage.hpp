#ifndef CIBART_TRANSIENT_STORAGE_HPP
#define CIBART_TRANSIENT_STORAGE_HPP

// Memory for the duration of one .Call, from R's transient storage (R_alloc):
// R reclaims it when the call returns or is jumped out of. The dbarts entries
// raise R errors - the run an interrupt, or a warning a handler turns into
// one - as do the glm fit and the argument checks, and a raise longjmps past
// every C++ frame without running a destructor or a delete, so what the
// sensitivity analysis holds across those calls comes from here.

#include <cstddef> // size_t
#include <new> // placement new
#if __cplusplus >= 201103L
#  include <type_traits> // is_trivially_destructible
#endif

#include <external/R.h> // R_alloc

namespace cibart {
  template <typename T>
  T* allocateTransient(std::size_t length)
  {
    return reinterpret_cast<T*>(R_alloc(length, sizeof(T)));
  }

  // a transient object is never destroyed, so T must not need to be
  template <typename T, typename A>
  T* createTransient(const A& argument)
  {
#if __cplusplus >= 201103L
    static_assert(std::is_trivially_destructible<T>::value, "a transient object is never destroyed");
#endif
    return new (static_cast<void*>(allocateTransient<T>(1))) T(argument);
  }

  template <typename T, typename A1, typename A2>
  T* createTransient(const A1& argument1, const A2& argument2)
  {
#if __cplusplus >= 201103L
    static_assert(std::is_trivially_destructible<T>::value, "a transient object is never destroyed");
#endif
    return new (static_cast<void*>(allocateTransient<T>(1))) T(argument1, argument2);
  }
}

#endif // CIBART_TRANSIENT_STORAGE_HPP
