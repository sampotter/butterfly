#pragma once

#include <assert.h>

/* `BF_DEFINE_VTABLE_STRUCT` is an X-macro which defines the virtual
 * table struct for an interface. To use it, first #define a macro of
 * the form:
 *
 *   <MACRONAME>(Type, Subtype, _) \
 *     _(Type, Subtype, ReturnType1, Args1...) \
 *     _(Type, Subtype, ReturnType2, Args2...) \
 *     ...
 *
 * where <MACRONAME> can be whatever you please. Then write:
 *
 *   #define INTERFACE MACRONAME
 *   BF_DEFINE_VTABLE_STRUCT(InterfaceName)
 *   #undef INTERFACE
 *
 * to have BF_DEFINE_VTABLE_STRUCT stamp out the vtable struct. */
#define BF_DEFINE_VTABLE_STRUCT_impl(Type, _2, X)   \
  typedef struct Bf##Type##Vtable {                 \
    INTERFACE(Type, _1, X)                          \
  } Bf##Type##Vtable;
#define BF_EMIT_VTABLE_STRUCT_ENTRY(_1, _2, ReturnType, FuncName, ...)  \
  ReturnType (*FuncName)(__VA_ARGS__);
#define BF_DEFINE_VTABLE_STRUCT(Type)                               \
  BF_DEFINE_VTABLE_STRUCT_impl(Type, , BF_EMIT_VTABLE_STRUCT_ENTRY)

#define BF_DEFINE_VTABLE_impl(Type, Subtype, _)    \
  static Bf##Type##Vtable Type##Vtbl = {           \
    INTERFACE(Type, Subtype, _)                    \
  };
#define BF_EMIT_VTABLE_ENTRY(Type, Subtype, _3, FuncName, ...)          \
  .FuncName = (__typeof__(&bf##Type##FuncName))bf##Subtype##FuncName,
#define BF_DEFINE_VTABLE(Type, Subtype) \
  BF_DEFINE_VTABLE_impl(Type, Subtype, BF_EMIT_VTABLE_ENTRY)

#define BF_DECLARE_INTERFACE_impl(Subtype, X)   \
  INTERFACE(_1, Subtype, X)
#define BF_EMIT_INTERFACE_ENTRY(_1, Subtype, ReturnType, FuncName, ...) \
  ReturnType bf##Subtype##FuncName(__VA_ARGS__);
#define BF_DECLARE_INTERFACE(Subtype) \
  BF_DECLARE_INTERFACE_impl(Subtype, BF_EMIT_INTERFACE_ENTRY)

#define BF_STUB(ReturnType, FuncName, ...)                  \
  _Pragma("GCC diagnostic push")                            \
  _Pragma("GCC diagnostic ignored \"-Wpedantic\"")  \
  ReturnType bf##FuncName(__VA_ARGS__) {                    \
    assert(false);                                          \
  }                                                         \
  _Pragma("GCC diagnostic pop")
