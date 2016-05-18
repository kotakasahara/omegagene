==================================
C++ Language Overview and Guidance
==================================

.. contents::

|Celeste| uses **C11** for C files and **C++14** for C++ files.  C++ has a lot
of features, but to keep the source code maintainable and easy to read, we will
avoid using some of them in the |Celeste| codebase.

In general, most programmers who claim to program in C++ do not actually
understand the language semantics of C++, and so are prone to programming C++
as if it were C.  This limits code maintainability and prevents the compiler
from making good optimizations wherever possible.  Despite its length, this page
is **only introductory**, as C++ is very complex, and the intent is to provide a
*high level* overview of just **some** of the core principles and important
ideas/features of the language.  Hopefully, this will serve tp bootstrap
onboarding |Celeste| developers into the terse programming mindset required for
effective C++ programming.


--------------------------------------------------------------
Select Core Language Ideas/Features For Better C++ Programming
--------------------------------------------------------------

^^^^^^^^^^^^
C++ and RAII
^^^^^^^^^^^^

**What is RAII?**

RAII is short for **Resource Acquisition Is Initialization**, and is also
synonymous with **Scope-based Resource Management**.  This means that holding a
resource is a class invariant.  Resource allocation (acquisition) is done
during object creation (specifically initialization), by the constructor, while
resource deallocation (release) is done during object destruction (specifically
finalization), by the destructor.

To start off with an example:

.. code-block:: c++

    struct Foo { /*...*/ }

    int function() {
        Foo foo1;               // foo1 is allocated in RAII-style
        Foo *foo2 = new Foo();  // foo2 is allocated in new/malloc()-style
        // do something with foo1 and foo2
        return 10;
    }

When ``function()`` is called, ``foo1`` and ``foo2`` are both allocated in
memory.  However, as soon as the code **leaves** the scope of the function,
``foo1`` is freed while ``foo2`` is not, unless a **delete** statement is
explicitly called before the return statement.

This is not the only issue.  One language feature of C++ is Exceptions, which
allow code to easily short-circuit and exit a scope when an invariant is broken.
In the above example, even if the programmer adds a ``delete`` statement
before the return statment, it may be possible for exceptions to be thrown
somewhere inside `function()` such that the ``delete`` is never called, leading
to a memory leak.  One can circumvent the problem by encapsulating each and
every code scope in large ``try``-``catch`` blocks, but that is not an optimal
solution.


**Is this the same as garbage collection or reference counting?**

**No**.  There is no garbage collection at runtime, nor is there reference
counting at compile-time.  All the compiler does is that it automatically
inserts cleanup code to free up resources at the exit point of a scope, provided
that the resources were declared/initialized in RAII style.


**Then why is there also a `delete` operator in the language?**

This is for compatibility with C, and this is for programmers who, for edge-case
reasons, don't want to use RAII.


**Why is this powerful?**

Nothing about RAII mentions the word "memory"; in fact, a "resource" can be
**anything**: a file handle, a thread, a mutex lock, a network socket, etc.
This is what makes RAII much more powerful than the garbage collection concept
in other programming languages - garbage collection is limited to memory
management, while RAII can be defined over anything that is potentially a
resource.

Consider the following example:

.. code-block:: c++

    void write_to_file (const std::string & message) {
        // mutex to protect file access (shared across threads)
        static std::mutex mutex;

        // lock mutex before accessing file
        std::lock_guard<std::mutex> lock(mutex);

        // try to open file
        std::ofstream file("example.txt");
        if (!file.is_open())
            throw std::runtime_error("unable to open file");

        // write message to file
        file << message << std::endl;

        // file will be automatically closed when leaving scope
        // (regardless of exceptions being thrown)
        // mutex will be unlocked (from lock destructor) when leaving
        // scope (regardless of exceptions being thrown)
    }

In a non-RAII programming language, the programmer needs to worry about writing
if-statements at every line in this function to check for errors and exceptions,
and manually insert method calls to unlock the mutex and close the file handle,
in order to prevent deadlocks and dangling file handles.

Very few programming languages support RAII (Ada, D, Rust), but this simple
idea is extremely powerful for keeping resource maintanence scalable.


**What is the important message, and how can I get started using RAII?**

Using RAII greatly simplifies resource management (allows for "exception-safe"
code) and **reduces overall code size**.  Consider the following class:

WIthout RAII:

.. code-block:: c++

    struct MiniCell {
        int *n_cells;
        int *idx_atom_cell;
        int *idx_cell_n_atoms;
        int *flg_dummy;
        int *idx_atom_cell_xy;
        int *idx_xy_head_cell;
        int *idx_xy_head_atom;

        // define default constructor
        MiniCell() {
            // allocate n_cells and set all elements to 0;
            n_cells = new int[100];
            for (int i=0; i < 100; ++i) {
                n_cells[i] = 0;
            }

            // allocate idx_atom_cell
            // ...
            // ...
        }

        // define default destructor
        ~MiniCell() {
            delete n_cells;
            // deallocate idx_atom_cell
            // ...
            // ...
        }

        // define default copy assignment operator
        MiniCell& operator= (const MiniCell &other) {
            // ...
            // ...
        }
    }

With RAII:

.. code-block:: c++

    struct MiniCell {
        std::vector<int> n_cells(100, 0),
                        idx_atom_cell,
                        idx_cell_n_atoms,
                        flg_dummy,
                        idx_atom_cell_xy,
                        idx_xy_head_cell,
                        idx_xy_head_atom;
    }

That's it!  This is provably resource-safe code; furthermore, manually writing
default constructors, destructors, and copy operators also becomes
**unneccessary** (more on this in the Rule of Five section).

Much of the C++ Standard Library uses the RAII idiom, and so it is easy to start
RAII-ifying the current |Celeste| codebase by identifying code that is manually
allocating raw memory, and replacing them with STL container counterparts
(mainly ``std::vector``).


**What custom/cool things can I do with RAII?**

Consider the following example:

.. code-block:: c++

    void dangerous_function() {
        turn_off_logging();

        try {
            if (not doSomething()) {
                turn_on_logging();
                return;
            }

            if (doNextThing() != 0) {
                turn_on_logging();
                return;
            }

            doLastThing();
            turn_on_logging();
            return;

        } catch {
            turn_on_logging();
        }
    }

To summarize, what we want is to turn off logging prior to running the body of
``dangerous_function()``, but we would like to turn the logging back on whenever
we exit this scope, regardless of *how* we exited the scope.

Using RAII:

.. code-block:: c++

    void dangerous_function() {
        struct Guard {
            ~guard() {
                turn_on_logging();
            }
        };
        Guard g;

        turn_off_logging();
        if (not doSomething()) return;
        if (doNextThing() != 0) return;
        doLastThing();
        return;
    }

This "always do this when exiting the scope" trick known as the
**Scope Guard Pattern**.


^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Rvalue References and Move Semantics
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**What are rvalues?**

This question is very difficult to answer with just a brief summary.  Consider
the following expression:

.. code-block:: c++

    Foo f = generate_foo_object();

In this example, ``f`` is considered an *lvalue* because it appears on the
*left* side of the assignment, while the value of calling
``generate_foo_object()`` is an *rvalue*, since it appears on the *right* side
of the assignment.  rvalues usually refer to the **unnamed temporary object/value**
that is created **before** it is assigned a binding to an lvalue.  Accessing
lvalues for use is very straightforward - the programmer just refers to the
variable's name.  However, accessing rvalues is possible only with modern C++.

This explanation is a gross oversimplification, because modern C++ categorizes
expressions into **5** types, some of them being subcategories of others.  They
are as follows:

* *lvalues*
* *glvalues*
* *xvalues*
* *prvalues*
* *rvalues*

A Venn diagram of the categorization will look like this:

.. code-block:: text

    ______ ______
   /      X      \
  /      / \      \
 |   l  | x |  pr  |
  \      \ /      /
   \______X______/
       gl    r

Excellent in-depth explanations can be found in the following links:

* http://stackoverflow.com/questions/3601602/what-are-rvalues-lvalues-xvalues-glvalues-and-prvalues
* https://www.justsoftwaresolutions.co.uk/cplusplus/core-c++-lvalues-and-rvalues.html


**What are rvalue references?**

Rvalue references enable code to distinguish an lvalue from an rvalue.  Consider
the following example:

.. code-block:: c++

    void print_string(const std::string &s) {
        std::cout << s;
    }

    void print_string(std::string &&s) {
        std::cout << s;
    }

    int main() {
        std::string s1("hello");
        print_string(s1);
        print_string("world"s);
    }

In C++98, only the first ``print_string`` definition is legal syntax, and that
function accepts both ``s1`` (const lvalue) and ``s"world"`` (rvalue) as
arguments.  However, since C++11, the second ``print_string`` definition is
allowed, and calling ``print_string(s"world");`` will invoke the second overload
of the function.

Excellent in-depth explanation of rvalue references can be found here:

* https://msdn.microsoft.com/en-us/library/dd293668.aspx
* http://www.cprogramming.com/c++11/rvalue-references-and-move-semantics-in-c++11.html

Together with the complex categorization of values, rvalue references allow C++
to define **move semantics**.


**What are move semantics?**

Consider the following code (simplified for illustrative purposes):

.. code-block:: c++

    vector<int> prepare_vector() {
        vector<int> v;
        // perform complex operations on v, such as adding and removing elements
        return v;
    }

    int main() {
        vector<int> v2 = prepare_vector();
    }

In C++98, the code will run the following steps to construct ``f2``:

1. Inside the body of ``prepare_vector()``, construct ``v``.
2. Inside the body of ``main()``, allocate memory for a ``Foo`` object.
3. **Copy** the memory referenced by ``v`` (the underlying heap buffer) inside the ``prepare_vector()`` scope into the memory referenced by ``v2`` inside the ``main()`` scope.

Having an implicit copy under the hood is extremely inefficient.  C++11 solves
this issue by introducing **move** semantics.  In the **very same** above code
example under C++11, the assignment operation in ``main()`` will directly "move"
the pointer from the underlying data buffer referenced by ``v``into the that of
``v2``, such that the extra work of copying the heap contents of ``v`` to ``v2``
is unnecessary.  Specifically, this compiler optimization enabled by C++11 is
called **copy elision**.  Optimizations such as this is made possible only with
C++11's categorization of lvalues and rvalues.


**When are rvalue references and move semantics most useful?**

One common use pattern of rvalue references is in defining move constructors and
move assignment operators for classes. A move constructor, like a copy
constructor, takes an instance of an object as its argument and creates a new
instance based on the original object. However, the move constructor can avoid
memory reallocation because the compiler knows it has been provideda temporary
object, so rather than copy the fields of the object, it can simply move them.

All the STL classes have move constructors and move assignment operators
defined; however, move constructors and move assignment operators may need to be
defined for user-defined classes.  Please see the Rule of Five section below for
further elaboration.


**How much should I care about all this / what do I need to do in my code?**

Formerly, in C++98, code such as the ``prepare_vector()`` example above was
discouraged on the grounds that extra copy operations were not desirable for
performance reasons. Thus, programmers wrote code similar to the following to
avoid the extra-copy problem:

.. code-block:: c++

    vector<int> prepare_vector(vector<int> &v) {
        // perform complex operations on v, such as adding and removing elements
        return;
    }

    int main() {
        vector<int> v2;
        prepare_vector(&v2);
    }

There are a few issues with this code:

1. The code employs value semantics as opposed to reference semantics, which make it difficult for the programmer to reasona about and maintain.
2. It prevents the compiler from making analysis and useful optimizations.

Code such as above example should be avoided now that C++11 solves the
extra-copy problem with move semantics.


^^^^^^^^^^^^^^^^
The Rule of Five
^^^^^^^^^^^^^^^^

**Automatically generated definitions for user-defined types in C++**

Consider the following code:

.. code-block:: c++

    struct Foo {
        std::string str;
        std::vector<int>
    }

When this code is compiled as is, the compiler will **automatically** generate
the following implicit definitions for the user defined type:

* default constructor (the constructor without arguments)
* destructor
* copy constructor
* move constructor
* copy assignment operator
* move assignment operator

This allows programmers to immediately be able write code such as the following:

.. code-block:: c++

    Foo f;                      // default constructor called
    f = Foo();                  // default move assignment called
    Foo g = f;                  // default copy assignment called
    Foo h(Foo());               // default move constructor called
    Foo j(f);                   // default copy constructor called

    vector<Foo> foos;
    foos.emplace_back(Foo());   // default move constructor called


**What is the Rule of Five?**

The automatic definitions generation is only guaranteed **under certain
conditions**.  For example, if a programmer explicitly defines a custom default
copy constructor, then the compiler will **not** automatically generate a
default constuctor nor a default move constructor, nor a default move assignment
operator.  For a complete table on what definitions are guaranteed to be
generated under which circumstances, please refer to the diagram in:

* http://stackoverflow.com/a/24512883

The Rule of Five (or the `Rule of Three <https://en.wikipedia.org/wiki/Rule_of_three_(C%2B%2B_programming)>`_
as it was known in C++98, before the introduction of move semantics) simply
states that if a programmer explicitly defines, for a user-defined type, a
method that is one of the five in the above list (excluding the default
destructor), then the programmer should also explicitly define the other four.
This prevents confusion and undefined behavior (resource leaks) when the type
is used in contexts that require any of these methods.  For example, if the
programmer defined a custom default copy constructor for a `Foo` class without
also explicitly defining a custom default move constructor, then calling
``emplace_back()`` on ``vector<Foo>`` will actually fail to compile.


**Why is this important to understand?**

Typically, the need to explicitly define any one of the six methods in the above
list is due to the fact that the attributes in the class were declared as
pointers and not in RAII style.  Synonymous with the Rule of Five is the Rule of
Zero, which states that if a class has all member resources whose types are RAII
and support the appropriate copy/move semantics, then the class itself does not
need to have any of the above 6 methods explicitly defined.  This reduces code
bloat *and* avoids undefined/unexpected behavior.  Please see the ``MiniCell``
code example from the RAII section of this guide.

To summarize: use RAII; don't use ``new``-allocated objects. A lot of code bloat
will go away, and and resource safety guarantees (i.e. deep copies on objects)
will come for free.


**What about cases where a custom constructor that takes in arguments is declared and defined?**

Consider the following code:

.. code-block:: c++

    struct Foo {
        int a=1, b=2;

        Foo(int tmp_a) {
            a = tmp_a;
        }
    }

    int main() {
        Foo f;      // ERROR!  Foo() not defined!
    }

If a constructor that takes in arguments is declared, then the default
constructor will not be implicitly *declared*, and so will not be implicitly
*defined*.  To fix this issue, simply add the following declaration:

.. code-block:: c++

    struct Foo {
        int a=1, b=2;

        Foo(int tmp_a) {
            a = tmp_a;
        }

        Foo() = default;    // declare this and the compiler will generate the definition appropriately
    }


^^^^^^^^^^^^^^
Smart Pointers
^^^^^^^^^^^^^^

**What are smart pointers?**

They are wrapper classes that enable RAII-style resource management for
manually-allocated resources (i.e. raw pointers).  The most common types of
smart pointers are ``std::shared_ptr`` and ``std::unique_ptr``.

Consider the following example:

.. code-block:: c++

    struct Foo { /*...*/ }

    void function() {
        Foo *foo1 = new Foo(42, "bar");                 // foo1 is allocated in new/malloc()-style; UNSAFE CODE
        std::unique_ptr<Foo> foo2(new Foo(42, "bar"));  // foo2 owns a Foo object; SAFE CODE
        auto foo3 = std::make_unique<Foo>(42, "bar");   // alternate syntax
        foo3->bar();                                    // Call Foo object's methods as if using a raw pointer
    }

``foo1`` can be a resource leak if the body of ``function()`` throws exceptions
later on, while ``foo2`` and ``foo3`` are automatically resource-managed.

**Should I use `std::unique_ptr` or `std::shared_ptr`?**

``std::shared_ptr`` should be used if there are two entities that are sharing
**ownership** responsibilities over a resource.  In most use cases, though, only
one entity ever owns/manages a resource, though that entity may pass the pointer
around for **use** in other code.  Therefore, ``std::unique_ptr`` will
suffice for most use cases.

**But I am using a C library that comes with its own special methods for
managing its defined datatypes. How can I use RAII here?**

Both ``std::shared_ptr`` and ``std::unique_ptr`` allow users to specify a custom
deallocator, which is particularly helpful when dealing with C-library
functions, which usually provide custom allocator and deallocator methods for
their data structures.

Consider the following real-world example from ``libgit``
(taken from https://dun.gs/posts/2014-11-06-c-cpp-memory.html):

Using raw pointers:

.. code-block:: c++

    void foo() {
        git_repository *repo = NULL;
        git_repository_open(&repo, "./testrepo");
        /* Do something with repo. */
        git_repository_free(repo);
    }

Using RAII with ``std::unique_ptr``:

.. code-block:: c++

    void foo() {
        git_repository *rawRepo = nullptr;
        git_repository_open(&rawRepo, "./testrepo");
        std::unique_ptr<git_repository, decltype(&git_repository_free)> repo { std::move(rawRepo), git_repository_free };
        /* Do something with repo. */
    }

**I can use `std::unique_ptr`, but how do I transfer resource ownership responsibility (i.e. outside of scope)?**

This is possible using ``std::move`` (more on this below).

**What mistakes can I make with smart pointers?**

Some common mistakes with smart pointer usage are listed in
http://www.acodersjourney.com/2016/05/top-10-dumb-mistakes-avoid-c-11-smart-pointers/


^^^^^^^^^^
Exceptions
^^^^^^^^^^

(This is a large topic to be filled in at a later date)

http://arne-mertz.de/2016/01/modern-c-features-keyword-noexcept/
http://arne-mertz.de/2016/05/raii-vs-exceptions/


^^^^^^^^^^^^^^^^^^^^^
``const`` Correctness
^^^^^^^^^^^^^^^^^^^^^

(To be filled in at a later date)


^^^^^^^^
``auto``
^^^^^^^^

(To be filled in at a later date)


^^^^^^^^^^^^^^^^^^^^
C++ Lambda Functions
^^^^^^^^^^^^^^^^^^^^

(To be filled in at a later date)


-------------------------------------------------------------
Other Language Bits For Better Understanding and Usage of C++
-------------------------------------------------------------

^^^^^^^^^^^^^^^^^^^^^
User-Defined Literals
^^^^^^^^^^^^^^^^^^^^^

(To be filled in a later date)

^^^^^^^^^^^^^^^^^^^^^^^^^^
Variable templates (C++14)
^^^^^^^^^^^^^^^^^^^^^^^^^^

(To be filled in a later date)
https://isocpp.org/wiki/faq/cpp14-language

^^^^^^^^^^^^
Initializers
^^^^^^^^^^^^

(To be filled in a later date)


^^^^^^^^^^^^^
``constexpr``
^^^^^^^^^^^^^

**What are constant expressions?**

Constant expressions are expressions that can be evaluated at compile-time.
Consider the following code:

.. code-block:: c++

    int factorial(int n) {
        return n <= 1? 1 : (n * factorial(n - 1));
    }

    int main() {
        cout << "4 factorial is: " << factorial(4);
    }

There is nothing preventing ``factorial(4)`` from being pre-computed, since the
argument is known at compile-time. Previously in C++98, there was no way to tell
the compiler to do this; however, C++11 provides the ``constexpr`` keyword,
which tells the compiler to evaluate the function if possible.

.. code-block:: c++

    constexpr int factorial(int n) {
        return n <= 1? 1 : (n * factorial(n - 1));
    }

Thus, in the above example, ``factorial(4)`` will be compiled away, and the
generated code will only reference the value ``24``.

**How are they different from C++ templates?**

C++ templates are a compile-time-only language feature can be thought of as
something similar to Lisp macros.  On the other hand, constexpr functions can
be used as regular functions at runtime, though if enough information is
provided, then the compiler is able to fully evaluate the code away at
compile-time.

However, C++ templates provide type-genericity in the function definitions,
while constexpr does not.  C++ templates enable functions over types as opposed
to values, while constexpr is limited to being defined over the space of
functions over values.

**When are constexpr functions useful?**

Having more code be evaluated at compile-time means will reduce the overall
run-time of the compiled program.  However, because ``constexpr`` optimization
only works when the values can be known at compile-time, its use can be limited.
One good use of constexpr is in defining functions over constants (i.e.
computing ``PI/4``, which can certainly be computed at compile-time).


^^^^^^^^^^^^^^^^^^^^^^^
``struct`` vs ``class``
^^^^^^^^^^^^^^^^^^^^^^^

**What's the difference?**

In C++, ``struct`` and ``class`` are the same and generate the same assembly
code, with the exception that all fields in a ``struct`` are **public** by default,
and all fields in a ``class`` are **private** by default.

**When should I be using which?**

If the intent of the data structure is to carry a lot of *data* around for
functions to openly read from and use, then consider using ``struct``.  If the
intent of the data structure is to carry *state* around, then consider using
``class`` instead.  They are always interchangeable with with the appropriate
``public`` and ``private`` keywords in place anyway.


^^^^^^^^^^^^^^^
``std::vector``
^^^^^^^^^^^^^^^

**Why should I use `std::vector`?**

``std::vector`` is a superior low-overhead replacement for manually-managed
memory (arrays).  It is an RAII data structure, comes with many convenience
methods, and allows for code-efficient copying and dynamic resizing.

.. code-block:: c++

    std::vector<Foo> foos = {{ foo1, foo2, foo3, ... }};
    auto copy_of_foos = foos;   // That's it!



**Why should I almost always use `std::vector` (over other STL data structures)?**

Of all the STL data structures, ``std::vector`` is the most well-defined and
plays well with existing STL algorithms (i.e. ``std::sort`` and ``std::find``).
In fact, other STL data structures (i.e. ``std::list`` or ``std::deque`` are
really implemented as arrays/vectors underneath for performance reasons.

**Should I use `push_back()` or `emplace_back()`?**

Beginning with C++11, ``std::vector`` comes with an ``emplace_back()`` method,
which behaves similarly to ``push_back()``.  However the subtle differences can
be seen in this example:

.. code-block:: c++

    std::vector<Foo> foos;
    Foo foo;

    foos.push_back(foo);                // Copies foo to the back of the vector
    foos.emplace_back(foo);             // Copies foo to the back of the vector
    foos.emplace_back(std::move(foo));  // Moves foo to the back of the vector; foo will be in an undefined state afterwards

    foos.push_back(Foo());              // Builds a Foo object somwhere else, then **copies** it over to the back of the vector
    foos.emplace_back(Foo());           // Builds a Foo object **directly** on the back of the vector; saves one copy operation

To conclude, ``emplace_back()`` is highly encouraged, but may not be useful in all use cases.

**What mistakes can I make with `std::vector`?**

The most common mistake with using ``std::vector``s is the
**dangling references mistake**.  Consider the following example:

.. code-block:: c++

    std::vector<int> vec = {{ 1, 2, 3, 4 }};  // vec has size of 4 and **capacity** of 4
    int *arr = &vec[0];                       // Use the underlying raw memory (array)
    vec.push_back(5)                          // vec is out of capacity, and so allocates itself a **new underlying buffer** of size 8,
                                              // copies the old elements over, and adds 5 to the end of the vector, and
                                              // **deletes the old underlying buffer**
    arr[0] = 100;                             // ERROR!  arr is now pointing to an invalid memory buffer since it has been freed!

Compilation for this code will pass, but the program will exhibit segfaults at
runtime. Running compiler diagnostics may be able to look for these kinds of
issues, but the best method for avoiding this error by making sure the size of
the ``std::vector`` is fixed (Using raw underlying pointers will be unavoidable
in the |Celeste| codebase anyway, since underlying libraries require raw
pointers to be passed in).


^^^^^^^^^^^^^^
``std::tuple``
^^^^^^^^^^^^^^

**What are tuples?**

Tuples are fixed-size ordered collections of **heterogeneous-typed** values, and
are a generalization of pairs (``std::pair`` in C++).  Tuples were first
directly implemented as product types in functional programming languages (such
as ML, OCaml, Haskell) as a non-object-oriented programming way to express a
composite set of values.

.. code-block:: c++

    auto x = make_tuple(42, Foo(), "hello")     // generates the <int,Foo,std::string> tuple { 42, Foo(), "hello" }

**Why should I use tuples?**

Tuples are useful for packaging values together into a data structure without
having to explicitly define a new struct/class.  In particular, if the component
values have little in relation to each other, and if there are no methods that
are bound to using them together in conjunction (i.e. an instance method), then
tuples can be used in place of structs/classes.

**When should I use tuples?**

If you find yourself writing functions that need to return multiple values, but
those values are not strongly related enough to justify packaging them together
into a new struct/class, then *consider* using ``std::tuple``s as opposed
to passing pointer handles.  If there is a large number of values to return, and
these values tend to be used together in many functions, then consider using a
struct/class to present the data instead.

Without tuples:

.. code-block:: c++

    int return_two_vals(float *float_return_value, std::string &string_return_value) {
        *float_return_value = 3.14;
        string_return_value = "hello";
        return 42;
    }

    int main() {
        float f;
        std::string s;
        int i = return_two_vals(&f, s);
    }

Using tuples:

.. code-block:: c++

    auto return_two_vals() {
        ruturn make_tuple(42, 3.14, "hello"s);
    }

    int main() {
        int i1;
        float f1;
        string s1;
        std::tie(i1, f1, s1) = return_two_vals();

        auto { i2, f2, s2 } = return_two_vals();    // C++17 structured-bindings (future) syntax
    }


-------------------------------------------------
General Rules to Follow for |Celeste| Development
-------------------------------------------------

Most of the rules below have been shamelessly copied from the `Allowed C++
Features for Gromacs Development
<http://www.gromacs.org/index.php?title=Developer_Zone/Programming_Guide/Allowed_C%2b%2b_Features>`_
but have been extended.  Most of these rules are not strict rules, but
developers should have a very good reason for deviating from them.

The `Google C++ Style Guide <https://google.github.io/styleguide/cppguide.html>`_
is a good source for additional points to think about. We might not
want to follow it point-to-point, but they have a thorough list of guidelines,
with justification for why they've made the choices.


^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Forbidden Language Features ("Do Not" Rules)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Don't use preprocessor defines for other than things directly related to
  preprocessing. Use templates and constexpr functions to generate code, and
  enums or const variables for constants.
* Don't use non-const references as function parameters. They make it impossible
  to tell whether a variable passed as a parameter may change as a result of a
  function call without looking up the prototype.
* Don't use C-style casts; use const_cast, static_cast or reinterpret_cast as
  appropriate.  Avoid dynamic_cast.
* Don't use ``std::auto_ptr``.  This is dangerous and deprecated since C++11.
* Don't use malloc, and limit the use of ``new``. Use container classes when
  appropriate instead of managing the memory everywhere manually.
* Don't use multiple inheritance. Inheriting from multiple pure interfaces is
  OK, as long as at most one base class (which should be the first base class)
  has any code.
* Don't write excessively deep inheritance graphs. Try to not inherit
  implementation just to save a bit of coding; follow the principle
  "inherit to be reused, not to reuse".
* Avoid using RTTI, even though the cost of it is not very high.
* Don't include unnecessary headers; his slows down compilation times.


^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Conditionally Allowed Language Features ("Proceed With Caution and Full Understanding" Rules)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Don't use Boost, except parts that all developers have agreed to be essential.
  Boost is a nice library, but but excessive template use slows down compilation
  significantly and may not work on all compilers.  Wherever possible, use the
  STL instead.
* Avoid overloading functions unless all variants really do the same thing,
  just with different types. Consider making the function names more descriptive
  instead.
* Avoid ``inline`` functions.  Modern C++ compilers are smart enough to optimize
  code away and largely ignore the ``inline`` keyword anyway
* Don't overload operators before thorough consideration whether it really is
  the best thing to do. Never overload &&, ||, or the comma operator, because
  it's impossible to keep their original behavior with respect to evaluation
  order.
* Avoid using default function arguments.
* Use exceptions for error handling.
* Avoid using complex templates, complex template specialization or
  techniques like SFINAE, unless you can show that it is actually useful.
  If nothing else, they can make the code more difficult to understand.


^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Encouraged Language Features ("Please Do" Rules)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Use namespaces. Everything in libceleste should be in a ``cst`` namespace.
  Don't use using in headers except possibly for aliasing some commonly-used
  names, and avoid file-level blanket ``using namespace cst`` and similar.
* Use the STL whereever possible.  They are more designed for optimized than
  what most programmers think.
* Write exception-safe code. All new code has to offer at least the basic or
  nothrow guarantee to make the above feasible.
* Use RAII for managing resources (memory, mutexes, file handles, ...).
* Use smart pointers (std::shared_ptr or std::unique_ptr) for manual memory
  management.
* Write const-correct code (no const_casts unless absolutely necessary).
* Use proper enums for variables that can only contain one of a limited set of
  values. C++ is much better than C in catching errors in such code.
* Follow the Rule of Five unless there is reason not to, such as making a class
  uncopyable.  If it is unnecessary to make the class copyable, then explicitly
  mark the copy constructor and assignment operator with the `delete` keyword.
* Declare all constructors with one parameter as explicit unless you really
  know what you are doing. Otherwise, they can be used for implicit type
  conversions, which can make the code difficult to understand, or even hide
  bugs that would be otherwise reported by the compiler. For the same reason,
  don't declare operators for converting your classes to other types without
  thorough consideration.


-----------------------
Compiler-Specific Notes
-----------------------

* When writing C++14 code, one must be aware that not all compilers fully
  support the C++14 specification (or even C++11 specification), and so
  developers must be vigilant in checking that new code be compileable under a
  relatively recent version for each of the major compilers on the platforms we
  support.
* MSVC supports only a subset of C99 and work-arounds are required in those
  cases.  One example is the required ``<ciso646>`` header to support the
  keyword versions of logical operators (i.e. ``and``, ``or`).

---------
Resources
---------

Books, etc