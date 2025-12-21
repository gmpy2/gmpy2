Release Notes
=============

.. currentmodule:: gmpy2

Changes in gmpy2 2.3.0
----------------------

* Improved support for mixed `mpfr` and `mpc` arithmetic.  (skirpichev)
* Preliminary support for free-threaded builds.  (skirpichev)
* Fix behavior and memory leaks for contexts.  (skirpichev)
* Fix memory leaks for `mpfr` caching.  (skirpichev)
* Fix rounding error in float(mpz).  (skirpichev)
* Add missing methods to `mpz` and `mpq`.  (skirpichev)
* Fix round(mpz, ndigits) if ndigits is negative.  (skirpichev)

Changes in gmpy2 2.2.2
----------------------

* Fix memory leak for mpfr/mpc cache, see issue 511.
* Fix memory leak in GMPy_CTXT_Exit(), see issue 515.
* Restore 'm' in library list, see issue 538.
* Adjust to_bytes() behavior to match CPython.
* Provide wheels for CPython v3.14.
* Provide wheels for Linux and Windows ARM64.

Changes in gmpy2 2.2.1
----------------------

* Fix internal use of char when int should be used. (jamesjer)
* Add :meth:`xmpz.bit_count`. (skirpichev)

Changes in gmpy2 2.2.0
----------------------

* Remove support for versions of Python < 3.7.  (skirpichev)
* Support more modern build tools.  (skirpichev)
* Use contextvars to manage gmpy2 contexts.  (casevh)
* _mpmath functions now use vectorcall protocol.  (casevh)
* Many documentation updates.  (skirpichev)
* Add :meth:`mpz.as_integer_ratio()` / :meth:`mpz.to_bytes()` and
  :meth:`mpz.from_bytes()`.  (skirpichev)
* Add :func:`is_probab_prime()` to directly expose the GMP
  behavior.  (skirpichev)
* :func:`gcd()`/:func:`lcm()` now uses vectorcall protocol.  (skirpichev)
* Expose :class:`context` type.  (skirpichev)
* Correct error in :func:`is_strong_bpsw_prp()`.  (casevh)
* Added :func:`prev_prime()` when GMP >= 6.3.  (sethtroisi)
* Change argument order of :func:`jn()` and :func:`yn()` to match
  MPFR.  (casevh)
* Fix documentation and code for
  :func:`is_extra_strong_lucas_prp()`.  (casevh)

Changes in gmpy2 2.1.5
----------------------

* Version bump to fix wheel issues.  No code changes.

Changes in gmpy2 2.1.4
----------------------

* Version bump to fix wheel issues.  No code changes.

Changes in gmpy2 2.1.3
----------------------

* Fix mpz(-3).is_prime().
* Add powmod_sec().
* Fix mpfr('inf') and mpfr('nan') if subnormalization is enabled.
* powmod() and powmod_sec() release the GIL.
* Fix error messages for iroot(x,n) for large n.
* Add powmod_base_list() and powmod_exp_list() (experimental).
* Fix gmpy2.mpq(mpq, int).
* Fix issues with INF, NAN, and mpfr("-0") when subnormalization is True

Changes in gmpy2 2.1.2
----------------------

* Code cleanup.
* Support Apple Silicon binary wheels.
* is_prime(-2) now returns False.  Issue #312.

Changes in gmpy2 2.1.1
----------------------

* Code cleanup.
* Properly return NOTIMPLEMENTED for unsupported arguments in
  ``**``.  Issue #319.

Changes in gmpy2 2.1.0
----------------------

* Improvements to setup.py.
* Add thread-safe contexts.
* MPFR and MPC are now required.
* Invalid Operation exception now raised for addition, etc.
* inverse() now raises exception if inverse does not exist.
* Add context methods.
* Major code refactoring required to properly support thread-safe contexts.
* `` __str__`` and ``__repr__`` no longer append "L" on Python 2.
* mpq(mpfr) now returns the exact result.
* Fix repr(mpc) for precision >325 bits.
* Intermediate conversions of Integer to mpfr are now done with the
  full precision of the Integer.
* Remove support for interaction with Decimal type.
* No longer attempt to override the memory allocation functions.
* Register gmpy2 types into the numeric tower.
* mpz(x) call int(x) if mpz() does not know how to convert x directly.
* Convert `mpz` to a type using ``__new__`` instead of a factory function.
* Bug fix for ``<<small mpfr>> ** <<small Python integer>>``.
* Compile with Python 3.11.

Changes in gmpy2 2.1.0rc2
-------------------------

* Documentation updates.
* Improvements to build environment.

Changes in gmpy2 2.1.0rc1
-------------------------

* Added support for embedded underscore characters in string literals.
* Allow GIL release for `mpz`/`xmpz`/`mpq` types only.

Changes in gmpy2 2.1.0b6
------------------------

* Improve argument type processing by saving type information to
  decrease the number of type check calls. Especially helpful
  for `mpfr` and `mpc` types. (Not complete but common operations
  are done.)
* Resolve bug in `mpfr` to `mpq` conversion; issue #287.
* Added limited support for releasing the GIL; disabled by default;
  see `context.allow_release_gil`.
* Refactored handling of inplace operations for `mpz` and `xmpz` types;
  inplace operations on `xmpz` will only return an `xmpz` result.
* Refactored handling of conversion to C integer types. Some
  exception types changes to reflect Python types.
* `gcd()` and `lcm()` now support more than two arguments to align with
  the corresponding functions in the math module.

Changes in gmpy2 2.1.0b5
------------------------

* Avoid MPFR bug in mfr_fac_ui (`factorial()`) on platforms where
  long is 32-bits and argument is >= 44787929.
* Fixed testing bugs with Python 2.7.
* Fixed ``mpz(0)`` to C long or long long.
* Fixed incorrect results in `f2q()`.
* Adjust test suite to reflect changes in output in MPFR 4.1.0.

Changes in gmpy2 2.1.0b4
------------------------

* Fix comparisons with `mpq` and custom rational objects.
* Fixes for some uncommon integer conversions scenarios.

Changes in gmpy2 2.1.0b3
------------------------

* Version bump only.

Changes in gmpy2 2.1.0b2
------------------------

* Many bug fixes.

Changes in gmpy2 2.1.0b1
------------------------

* Added `cmp()` and `cmp_abs()`.
* Improved compatibility with the :mod:`numbers` module protocol.
* Many bug fixes.

Changes in gmpy2 2.1.a05
------------------------

* Fix `qdiv()` not returning `mpz` when it should.
* Added `root_of_unity()`.

Changes in gmpy2 2.1.0a4
------------------------

* Fix issue 204; missing file for Cython.
* Additional support for MPFR 4

   - Add `fmma()` and `fmms()`.

Changes in gmpy2 2.1.0a3
------------------------

* Updates to ``setup.py``.
* Initial support for MPFR4

   - Add `mpfr_nrandom()`
   - `mpfr_grandom()` now calls nrandom twice; may return different values
     versus MPFR3.
   - Add `rootn()`; same as `root()` except different sign when taking even
     root of -0.0.

Changes in gmpy2 2.1.0a2
------------------------

* Revised build process.
* Removal of unused code/macros.
* Cleanup of Cython interface.

Changes in gmpy2 2.1.0a1
------------------------

* Thread-safe contexts are now supported. Properly integrating thread-safe
  contexts required an extensive rewrite of almost all internal functions.
* MPFR and MPC are now required. It is no longer possible to build a version
  of gmpy2 that only supports the GMP library.
* The function ``inverse()`` now raises an exception if the inverse does not
  exist.
* Context methods have been added for MPFR/MPC related functions.
* A new context option (`~context.rational_division`) has been added that
  changes the behavior of integer division involving `mpz` instances to return
  a rational result instead of a floating-point result.
* gmpy2 types are now registered in the numeric tower of the
  :mod:`numbers` module.
* In previous versions of gmpy2, ``mpz()`` was a factory function that
  returned an  `mpz` instance.  It is now an actual type. The same
  is true for the other gmpy2 types.
* If a Python object has an ``__mpz__`` method, it will be called bye ``mpz()``
  to allow an unrecognized type to be converted to an mpz instance. The same is
  true for the other gmpy2 types.
* A new C-API and Cython interface has been added.

Changes in gmpy2 2.0.4
----------------------

* Fix `bit_scan0()` for negative values.
* Changes to ``setup.py`` to allow static linking.
* Fix performance regression with mpmath and Python 3.

Changes in gmpy2 2.0.3
----------------------

* Fix `lucas2()` and `atanh()`; they were returning incorrect values.

Changes in gmpy2 2.0.2
----------------------

* Rebuild Windows binary installers due to MPIR 2.6.0 bug in `next_prime()`.
* Another fix for `is_extra_strong_lucas_prp()`.

Changes in gmpy2 2.0.1
----------------------

* Updated ``setup.py`` to work in more situations.
* Corrected exception handling in basic operations with `mpfr` type.
* Correct ``InvalidOperation`` exception not raised in certain circumstances.
* `invert()` now raises an exception if the modular inverse does not exist.
* Fixed internal exception in `is_bpsw_prp()` and `is_strong_bpsw_prp()`.
* Updated `is_extra_strong_lucas_prp()` to latest version.

Changes in gmpy2 2.0.0
----------------------

* Fix segmentation fault in ``_mpmath_normalize()`` (an undocumented helper
  function for mpmath).  (casevh)
* Fix issues when compiled without support for MPFR.  (casevh)
* Conversion of too large an `mpz` to `float` now raises `OverflowError`
  instead of returning ``inf``.  (casevh)
* Renamed ``min2()/max2()`` to `minnum()`/`maxnum()`.  (casevh)
* The build and install process (i.e. ``setup.py``) has been completely
  rewritten.  See the Installation section for more information.  (casevh)
* `get_context()` no longer accepts keyword arguments.  (casevh)

Known issues in gmpy2 2.0.0
-----------------------------

* The test suite is still incomplete.

Changes in gmpy2 2.0.0b4
------------------------

* Added ``__ceil__()``, ``__floor__()``, ``__trunc__()``, and ``__round__()``
  methods to `mpz` and `mpq` types.  (casevh)
* Added ``__complex__()`` to `mpc` type.  (casevh)
* ``round(mpfr)`` now correctly returns an `mpz` type.  (casevh)
* Add mpz.denominator and mpz.numerator.  (casevh)
* If no arguments are given to `mpz`, `mpq`, `mpfr`, `mpc`, and `xmpz`,
  return 0 of the appropriate type.  (casevh)
* Fix broken comparison between `mpz` and `mpq` when `mpz` is on
  the left.  (casevh)
* Added ``__sizeof__()`` to all types. Note: :func:`sys.getsizeof` calls
  ``__sizeof__()`` to get the memory size of a gmpy2 object. The returned
  value reflects the size of the allocated memory which may be larger than
  the actual minimum memory required by the object.  (casevh)

Known issues in gmpy2 2.0.0b4
-----------------------------

* The new test suite (``test/runtest.py``) is incomplete and some tests fail on
  Python 2.x due to formatting issues.

Changes in gmpy2 2.0.0b3
------------------------

* `mp_version()`, `mpc_version()`, and `mpfr_version()` now return normal
  strings on Python 2.x instead of Unicode strings.  (casevh)
* Fix warnings when shifting 32-bit integer by 32 bits.  (casevh)
* Faster conversion of the standard library `~fractions.Fraction` type
  to `mpq`.  (casevh)
* Improved conversion of the `~decimal.Decimal` type to `mpfr`.  (casevh)
* Consistently return `OverflowError` when converting ``inf``.  (casevh)
* Fix `mpz.__format__()` when the format code includes "#".  (casevh)
* Add `is_infinite()` and deprecate ``is_inf()``.  (casevh)
* Add `is_finite()` and deprecate ``is_number()``.  (casevh)
* Fixed the various ``is_XXX()`` tests when used with `mpc`.  (casevh)
* Fixed error handling with mpc(); mpc(1,"nan") is properly handled.  (casevh)
* Added caching for `mpc` objects.  (casevh)
* Faster code path for basic operation is both operands are `mpfr`
  or `mpc`.  (casevh)
* Fix `mpfr` + `float` segmentation fault.  (casevh)

Changes in gmpy2 2.0.0b2
------------------------

* Allow `xmpz` slice assignment to increase length of `xmpz` instance by
  specifying a value for stop.  (casevh)
* Fixed reference counting bug in several ``is_xxx_prp()`` tests.  (casevh)
* Added `~xmpz.iter_bits()`, `~xmpz.iter_clear()`, `~xmpz.iter_set()` methods
  to `xmpz`.  (casevh)
* Added `powmod()` for easy access to three argument :func:`pow()`.  (casevh)
* Removed ``addmul()`` and ``submul()`` which were added in 2.0.0b1 since they
  are slower than just using Python code.  (casevh)
* Bug fix in gcd_ext when both arguments are not `mpz`.  (casevh)
* Added `ieee()` to create contexts for 32, 64, or 128 bit `float`'s.  (casevh)
* Bug fix in `context()` not setting `~context.emax`/`~context.emin` correctly
  if they had been changed earlier.  (casevh)
* Contexts can be directly used in with statement without requiring
  `set_context()`/`local_context()` sequence.  (casevh)
* `local_context()` now accepts an optional context.  (casevh)

Changes in gmpy2 2.0.0b1
------------------------

* Rename to gmpy2 to allow backwards incompatible changes (casevh)
* Renamed 'mpf' to 'mpfr' to reflect use of MPFR (casevh)
* Renamed functions that manipulate individual bits to ``bit_XXX()`` to align
  with :meth:`~int.bit_length()`.
* Added caching for `mpq`.  (casevh)
* Added ``rootrem()``, `fib2()`, `lucas()`, `lucas2()`.  (casevh)
* Support changed hash function in Python 3.2.  (casevh)
* Added `is_even()`, `is_odd()`.  (casevh)
* Add caching of the calculated hash value.  (casevh)
* Add `xmpz` (mutable `mpz`) type.  (casevh)
* Fix `mpq` formatting issue.  (casevh)
* Add read/write bit access using slices to `xmpz`.  (casevh)
* Add read-only bit access using slices to `mpz`.  (casevh)
* Add `pack()`/`unpack()` methods to split/join an integer into n-bit
  chunks.  (casevh)
* Add support for MPFR (casevh)
* Removed fcoform float conversion modifier.  (casevh)
* Add support for MPC.  (casevh)
* Added context manager.  (casevh)
* Allow building with just GMP/MPIR if MPFR not available.  (casevh)
* Allow building with GMP/MPIR and MPFR if MPC not available.  (casevh)
* Removed most instance methods in favor of gmpy2.function. The general
  guideline is that *properties* of an instance can be done via instance
  methods but *functions* that return a new result are done using
  gmpy2.function.  (casevh)
* Added ``__ceil__()``, ``__floor__()``, and ``__trunc__()`` methods since they
  are called by :func:`math.ceil()`, :func:`math.floor()`, and
  :func:`math.trunc()`.  (casevh)
* Removed ``gmpy2.pow()`` to avoid conflicts.  (casevh)
* Removed ``gmpy2._copy()`` and added `xmpz.copy()`.  (casevh)
* Added support for ``__format__()``.  (casevh)
* Added ``as_integer_ratio()``, ``as_mantissa_exp()``,
  ``as_simple_fraction()``.  (casevh)
* Updated rich_compare.  (casevh)
* Require MPFR 3.1.0+ to get divby0 support.  (casevh)
* Added `fsum()`, `degrees()`, `radians()`.  (casevh)
* Updated random number generation support.  (casevh)
* Changed license to LGPL 3+.  (casevh)
* Added `lucasu()`, `lucasu_mod()`, `lucasv()`, and `lucasv_mod()`.  (casevh)
  *Based on code contributed by David Cleaver.*
* Added probable-prime tests.  (casevh)
  *Based on code contributed by David Cleaver.*
* Added `to_binary()`/`from_binary()`.  (casevh)
* Renamed ``numdigits()`` to `~mpz.num_digits()`.  (casevh)
* Added keyword precision to constants.  (casevh)
* Added ``addmul()`` and ``submul()``.  (casevh)
* Added ``__round__()``, `round2()`, `round_away()` for `mpfr`.  (casevh)
* ``round()`` is no longer a module level function.  (casevh)
* Renamed module functions ``min()/max()`` to ``min2()/max2()``.  (casevh)
  No longer conflicts with builtin :func:`min()` and :func:`max()`
* Removed ``set_debug()`` and related functionality.  (casevh)
* Removed mpf.setprec(), use mpf.round() (casevh)
* Fix test compatibility with Python 3.1.2 and 3.2 (casevh)
* Remove old random number functions, to be replaced later (casevh)
* Remove tagoff option (casevh)
* Debug messages only available if compiled with -DDEBUG (casevh)
* Renamed context() -> local_context(), new_context() -> context() (casevh)
* Added get_context() (casevh)

Changes in gmpy 1.11
--------------------

* Recognize True/False (bug in 1.10) (casevh)
* Optimize argument handling (casevh)
* Added caching for mpz (casevh)

Changes in gmpy 1.10
--------------------

* Remove dependancy on pymemcompat.h (casevh)
* Remove callback (casevh)
* Added support for -DMPIR to include MPIR instead of GMP (casevh)
* Major code revisions to add support for Python 3.x (casevh)
* Fixed bug in binary() and qbinary() (casevh)
* Fixed bug in rich comparisons (casevh)
* Added % and divmod support to mpq and mpf (casevh)
* Changed memory allocation functions to use PyMem (casevh)
* Removed small number interning (casevh)
* Added tdivmod, cdivmod, and fdivmod (casevh)
* Added more helper functions for mpmath (casevh)
* Faster mpz<>PyLong conversion (casevh)
* Faster hash(mpz) (casevh)

Changes in gmpy 1.04
--------------------

* Avoid GMP/mingw32 bug when converting very small floats to mpz. (casevh)
* Significant performance improvement for long->mpz and mpz->long. (casevh)
* Added "rich comparisons" to mpz, mpq and mpf types (aleaxit)
* Added additional tests (casevh, aleaxit)
* Fixed bug when converting very large mpz to str (casevh)
* Faster conversion from mpz->binary and binary->mpz (casevh)
* Added support for pickling (casevh)
* Added divexact (casevh)
* Fixed mpf comparisons by rounding mpf results when GMP returns
  a longer result. Added fround() (casevh)
* Added bit_length (Thanks Mario Pernici)
* Added helper functions for mpmath (casevh)
* Faster conversion from mpq->binary and binary->mpq (casevh)
* Recognize MPIR, mpir_version() (casevh)

Changes in gmpy 1.03
--------------------

* Fixed the bug that caused crashes on gmpy.mpf(float('inf')) and
  other such conversions, implicit and explicit
* Fixed a bug in get_zconst's prototype affecting 64-bit machines,
  thanks to Gary Bunting
* Fixed a bug in hashing on 64-bit systems. hash(long) now equals
  hash(mpz) for large values. (casevh)
* Changed int() to return a long value instead of OverFlowError.
  Complies with PEP 237. (casevh)
* Added support in setup.py for darwinports/macports build of GMP
  on MacOSX. (aleaxit)

Changes in gmpy 1.02
--------------------

* fix warning in comparison of mpq's
* added support of mpq('12.34') [[string w/o a slash, but with a dot]]
* fixes for 64-bit build (thanks to a patch by dmcooke)
* added experimental support for decimal.Decimal (and user-coded types)
  via wider use of special conversion methods (if present) and their
  sly insertion on-the-fly into the decimal.Decimal class (!)
* two bugfixes, thanks to Simon Burton
* Brought back into C89 compliance (thanks to Chip Turner), had
  drifted to C99 (declarations in the middle of the code).
* Python 2.5 support (Py_ssize_t, __index__) thanks to Chip Turner
* Pushed coverage to 93.3% (missing only "sanity check" level error
  tests [mostly for out-of-memory conditions], output to stderr
  conditioned by global.debug, & a couple of very obscure cases)

Changes in gmpy 1.01
--------------------

* cleanups, ensure support for Python 2.4.1 on MacOSX 10.4/XCode 2.1
  as well as Python 2.2 and 2.3 (on MacOSX and Linux)
* fixed memory leak on divm (thanks to mensanator@aol.com)
* fixed bug on mpq('123') [[str2mpq on string w/o a slash]]
* added floordiv and truediv operators, and tests for them
* NOT tested on GMP 3 (have none left around...), ONLY on GMP 4.*

Changes in gmpy 1.0
-------------------

* minor cleanups, ensure support for Python 2.3
* fixed misdiagnosis of some argument counts in macro
* SELF_ONE_ARG_CONVERTED (tx to Paul Rubin!)

Changes in gmpy 0.9
-------------------

* change ValueError to OverflowError for 'too-large' errors
* fix bug in mpq_pow (negative base, exp. with odd denominator)
  (fix now corrected -- _even_ denominator is the error!)
* fixed gcc warnings reported by K. Briggs
* support GMP 4 (but added no GMP4-only functionality yet)
* updated tests to 0.9, better coverage

Changes in gmpy 0.8
-------------------

(again, requests & suggestions by great Pearu!)

* raise test coverage 72.5% -> 90.0%
* introduced callbacks (not documented/tested for now;
  Pearu will test/support/document in PySymbolic)
* some errors went undiagnosed, caused crash: now fixed
* workaround for GMP bug(?s?) in mpz_fits\_... (?)
* added exposure of mpf\_ sqrt and pow_ui

Changes in gmpy 0.7
-------------------

Good feedback from Keith Briggs, some advice from Tim Peters and Fred Lundh ---
thanks all!

* fixed bug of '"%d" where "%ld" was meant' in many places
  and other sundry minor warnings given by gcc
* fixed hash (delegating to Python) so mp[nqz](x) will
  produce the same value as hash(x) for any Python number x
* workaround for GMP 3.1.1 bug, mpz_root wrongly returning
  'exact' for non-exact root if dest==source, which stopped
  needed value-error for inexact mpq**mpq operations
* determined correct 'actual precision' of floats
* explicitly stored precision with binary-form mpf's
* extended explicit-bits request to all ->mpf operations
  (good in itself, plus, preparing for future MPFR)
* removed the limitation of no binary-form for <0 mpz
* introduced macros to parse args, for conciseness

Changes in gmpy 0.6
-------------------

(lots of good ideas from Pearu once more!-)

* fixed silly bugs in kronecker and mpq_abs
* gmpy-level workaround for scan0/scan1 bugs (?) in gmp 3.1.1
* added qdiv; anynum->mpq substituted for all such conversions
  (also anynum->mpz and anynum->mpf by analogy, with care!)
* added global.fcoform for optional use of intermediate string in
  float2mpf (used for any float->mpf conversion)
* added set_fcoform function for global.fcoform access
* general cleanup of sources; added alloca for msvc++;
  - many sundry minor bugfixes & uniformization;
  - a little useful refactoring (more would be good...)
* added caching of mpq objects
* power for mpq
* stern-brocot algorithm for mpf->mpq (also exposed as f2q)
  - also used for float->mpq
  - with stricter tracking of mpf's requested-precision
  - added getrprec method to mpf, getrprec module-function
* exposed ceil, floor and trunc methods/functions for mpf's
* changed a couple exceptions from value to zerodivision
* added 'qual' and 'floa' options to gmpy.rand

Changes in gmpy 0.5
-------------------

* added jacobi, legendre, kronecker
* added random-number generation, seed set/save, shuffling
* added mpq (at last!-)

Changes in gmpy 0.4
-------------------

* split gmpy.c/gmpy.h introducing C-API interface (Pearu's suggestion)
* cleanup some casts using Pearu's new macros
* further cache-tweaks at Pearu's suggestion (macros introduced)
* added sign (Pearu's request), getbit, setbit
* added docstrings
* renamed copy functions to start with _ ('internal, private')
* added .comb as a synonym of .bincoef

Changes in gmpy 0.3
-------------------

* performance tweaks via mpz-caching & fixed-constants
* added get/set functions for zcache, zco min/max
* added get-only function for versions (of gmp, and of gmpy)
* removed all 'traces' of mutability (to be re-done... much later!)
* cleaned up all of the mpz_cmp_ui(X,0) to mpz_sgn(X)
* cleaned up Py_BuildValue usage (N vs O, explicit-() for tuples)
* added numdigits, lowbits, root, next_prime, invert, popcount,
* hamdist, scan0, scan1
* renamed bin to bincoef

Changes in gmpy 0.2
-------------------

15 Nov 2000

* pre-alpha: bugfixes re formatting (tx, Peanu!)
* no tags on oct() and hex() of mpz's
* insert 'tagoff' in options (gmpy.mpz() vs mpz() in repr) (for Peanu!)
* speedups for _nonzero & _cmp (tx, Peanu!)
* slight speedup (7/8%?) for excess reallocs 4<->8 bytes (Peanu's help!)
* added copy/fcopy; bin; fib; remove

Changes in gmpy 0.1
-------------------

6 Nov 2000

* pre-alpha --- first placed on sourceforge
