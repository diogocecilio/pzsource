// Emacs will be in -*- Mode: c++ -*-
//
// ************ DO NOT REMOVE THIS BANNER ****************
//
//  Nicolas Di Cesare <Nicolas.Dicesare@ann.jussieu.fr>
//  http://www.ann.jussieu.fr/~dicesare
//
//            CEMRACS 98 : C++ courses, 
//         templates : new C++ techniques 
//            for scientific computing 
// 
//********************************************************
//
//  A short implementation ( not all operators and 
//  functions are overloaded ) of 1st order Automatic
//  Differentiation in forward mode (FAD) using
//  EXPRESSION TEMPLATES.
//
//********************************************************
#ifndef _fadlog_h_
#define _fadlog_h_

#include "pzreal.h"


#define FAD_LOG_MACRO(OP)                                                      \
  template <class T>                                                           \
  inline bool operator OP(const Fad<T> &a, const Fad<T> &b) {                  \
    return (a.val() OP b.val());                                               \
  }                                                                            \
                                                                               \
  template <                                                                   \
      class A, class B,                                                        \
      typename std::enable_if<((std::is_integral<B>::value ||                  \
                                is_arithmetic_pz<B>::value) &&     \
                               (std::is_integral<A>::value ||                  \
                                is_arithmetic_pz<A>::value)),      \
                              int>::type * = nullptr>                          \
  inline bool operator OP(const Fad<A> &a, const B &b) {                       \
    return (a.val() OP b);                                                     \
  }                                                                            \
                                                                               \
  template <                                                                   \
      class A, class B,                                                        \
      typename std::enable_if<((std::is_integral<B>::value ||                  \
                                is_arithmetic_pz<B>::value) &&     \
                               (std::is_integral<A>::value ||                  \
                                is_arithmetic_pz<A>::value)),      \
                              int>::type * = nullptr>                          \
  inline bool operator OP(const A &a, const Fad<B> &b) {                       \
    return (a OP b.val());                                                     \
  }                                                                            \
                                                                               \
  template <class T>                                                           \
  inline bool operator OP(const FadExpr<T> &a, const FadExpr<T> &b) {          \
    return (a.val() OP b.val());                                               \
  }                                                                            \
                                                                               \
  template <class T>                                                           \
  inline bool operator OP(const T &a, const FadExpr<T> &b) {                   \
    return (a OP b.val());                                                     \
  }                                                                            \
                                                                               \
  template <class T>                                                           \
  inline bool operator OP(const FadExpr<T> &a, const T &b) {                   \
    return (a.val() OP b);                                                     \
  }
                                                 

FAD_LOG_MACRO(==)
FAD_LOG_MACRO(!=)
FAD_LOG_MACRO(<)
FAD_LOG_MACRO(>)
FAD_LOG_MACRO(<=)
FAD_LOG_MACRO(>=)
FAD_LOG_MACRO(<<=)
FAD_LOG_MACRO(>>=)
FAD_LOG_MACRO(&)

#undef FAD_LOG_MACRO


template <class T> inline bool operator !(const Fad<T> &a) {
    return ( !a.val() );
}

template <class T> inline bool operator !(const FadExpr<T> &a) {
    return ( !a.val() );
}

#endif
