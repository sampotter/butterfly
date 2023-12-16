Implementation
==============

*This is a living design document---not all of the code reflects these guidelines yet.*

Style guidelines
----------------

A few simple naming conventions should be adopted throughout ``butterfly``'s source:

- Apart from macro, enum, and global constant names, symbols should be in camelCase.
- Macros, enums, and global constants should be written in ALL_CAPS with words separated using underscores.
- All *publicly exported* symbols should be prefixed with ``bf``, ``Bf`` or ``BF_`` (*namespaced*), depending on whether:

  1. The symbol is a struct, enum, union, or typedef (``Bf``).
  2. It's a variable or a function name (``bf``).
  3. It's a macro, enum, or global constant name (``BF_``).

- Static symbols should not be namespaced, but instead should be in camelCase, PascalCase, or ALL_CAPS, depending on whether the symbol belongs to case 1, 2, or 3 above.

Since `butterfly <https://butterfly.github.io>`_ uses its own ad-hoc object system, it's important to follow consistent naming conventions and to order function arguments predictably. If we want to define a ``Frobnicate`` "method" on the type ``BfWidget``, we would declare::

  void bfWidgetFrobnicate(BfWidget *widget);

It's also helpful to adhere to the following naming convention for functions used to set up and tear down objects:

- ``Init*`` is used for functions which initialize an object.
- ``New*`` is used for functions which allocate memory for an object and possibly initialize it.

For example::

  BfWidget *bfWidgetNewWithDefaultParams();
  void bfWidgetInitWithDefaultParams(BfWidget *widget);

For some types, both types of functions may be provided, but it isn't necessary. For example, if ``BfWidget`` is opaque, then an ``Init`` function doesn't make sense, since we can't use ``sizeof(BfWidget)`` to determine how much memory to allocate. On the other hand, ``bfWidgetNew*`` is fine, since ``sizeof(BfWidget)`` will be available statically to the module defining ``BfWidget``. This lets us hide the implementation of ``BfWidget``, which is one tool that can be used to manage complexity.

On the other hand, ``BfWidget`` may not be that complex, in which case making it opaque would be its own source of complexity. In this case, for performance reasons, it might be preferable to allocate ``BfWidget`` on the stack and initialize it using an ``Init`` function.

If both functions are included, then the convention is that the ``New`` function first allocates and then initializes the type.

In summary:

- If a type is opaque, then ``New`` functions provide the only way of instantiating them. Otherwise, ``New`` functions are convenience functions which both allocate and initialize a type.
- If a type if opaque, it won't have any ``Init`` functions. If it isn't opaque, then ``Init`` functions can be useful for initializing stack-allocated types.

Having both ``New`` and ``Init`` functions leaves the door open to other useful allocation strategies---for instance, for pool allocation, only ``Init`` will be useful.

There are three styles of teardown functions: ``Deinit``, ``Dealloc``, ``DeinitAndDealloc``, and ``Delete``. The ``Deinit`` and ``Dealloc`` should just be the inverse of ``Init`` and ``New``. However, while there can be many different ``Init`` and ``New`` functions per type, there will only ever be one ``Deinit`` and ``Dealloc`` function. The ``DeinitAndDealloc`` functions combine ``Deinit`` and ``Dealloc``. In addition, to cope with object inheritance, we have ``Delete`` functions, which recursively deinitialize the parent types (i.e., a "virtual destructor") but which otherwise take the place of ``DeinitAndDealloc``.

Since memory allocation/deallocation is separated from object initialize/deinitialization, objects can be allocate and deallocated once, and then initialized and deinitialized in a loop repeatedly.

Position of the ``const`` keyword
`````````````````````````````````

Where ``const`` should be placed in C is contentious. The most consistent thing to do is to place it *after* the type, pointers, or qualifier it modifies. For all qualifiers and pointers, it *must* be placed to the right, but it can be placed either before or after the leading type. Many C programmers prefer to place the ``const`` keyword before the leading type since it follows the rules for adjective placement in English: "a constant integer" (``const int``) rather than "an integer constant" (``int const``). Obviously not all programmers are native English speakers, so the choice of whether to place ``const`` before or after the leading type is a subjective one.

We choose to place the ``const`` *before* the leading type (*or at least we will once we update the code...*) for one practical reason: Cython does not allow ``const`` to be placed after the leading type, making it annoying to copy and paste C headers into Cython `.pxd` files. Every ``BfType const`` must be converted to ``const BfType`` by hand.

We also never use ``const`` to modify pointers or qualfiers, so consistency is unimportant.

Coding style
------------


Ownership
---------

*copy vs view vs steal*

Error handling
--------------

Object system
-------------

Although butterfly is written in the C programming language, it makes use of object-oriented programming (OOP). The primary motivation for this is the inclusion of a hierarchy of matrix types. Since butterfly factorizations can be built up recursively in terms of operations on block matrices, having matrix types which can be nested hierarchically, concatenated, or blocked together in arbitrary ways at runtime is essential.

The object system supports single inheritance which can be multiple levels deep. Objects methods are implemented using virtual tables. Since C doesn't directly support OOP, we adopt a system of ad hoc conventions to implement the desired features.

A class consists of a struct defining the type, a virtual table (or "vtable"), a set of methods defined on it, and its subclasses. For a simple example, consider a base class ``Animal``::

  struct Animal {
    struct AnimalVtable vtable;
    ...
  };

A subclass is a struct whose first member is named ``super`` and which has the same type as the base class. Continuing the example, consider two subclasses ``Dog`` and ``Cat``::

  struct Dog {
    struct Animal super;
    ...
  };

  struct Cat {
    struct Animal super;
    ...
  };

Let's assume ``Animal`` has one method, ``MakeNoise``. To implement this method, we first define the vtable::

  struct AnimalVtable {
    void (*MakeNoise)(struct Animal *);
  };

We can implement this method as follows::

  void animalMakeNoise(struct Animal *animal) {
    animal->vtable->MakeNoise(animal);
  }

  void dogMakeNoise(struct Dog *dog) {
    printf("woof!\n");
  }

  void catMakeNoise(struct Cat *cat) {
    printf("meow!\n");
  }

When we create an "instance" of ``Dog`` or ``Cat``, we first define a variable of type ``AnimalVtable`` for each subclass, and then create the instance. For example::

  struct AnimalVtable DogVtable = {
    .MakeNoise = __typeof__(&animalMakeNoise)dogMakeNoise
  };

  struct Dog *newDog() {
    struct Dog *dog = malloc(sizeof(struct Dog));
    dog->vtable = &DogVtable;
    ...
    return dog;
  }

To illustrate runtime polymorphic behavior, we simple create a new instance of ``Dog``, "cast" it to ``Animal``, and call ``animalMakeNoise``. The call will be dispatched to ``dogMakeNoise`` so that ``woof!`` is print to standard out::

  struct Dog *dog = newDog();
  struct Animal *animal = &dog->super;
  animalMakeNoise(animal); // woof!

Bear in mind that all of this is completely ad hoc and relies on the programmer to establish and follow their own set of conventions. None of this is enforced syntactically by the compiler. We will describe the conventions adopted in butterfly in the following sections so they are less opaque.
