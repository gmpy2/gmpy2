gmpy2 is an optimized, C-coded Python extension module that supports fast
multiple-precision arithmetic.  gmpy2 is based on the original gmpy module.
gmpy2 adds support for correctly rounded multiple-precision real arithmetic
(using the MPFR library) and complex arithmetic (using the MPC library).

Version 2.3 (development)
-------------------------

* Preliminary support for free-threaded builds.
* Fix behavior and memory leaks for contexts.
* Fix memory leaks for `mpfr` caching.
* Fix rounding error in float(mpz).
* Add missing methods to `mpz` and `mpq`.
* Fix round(mpz, ndigits) if ndigits is negative.
* Improved support for mixed `mpfr` and `mpc` arithmetic.

gmpy2 2.2.2
-----------

* Fix memory leak for mpfr/mpc cache, see issue 511.
* Fix memory leak in GMPy_CTXT_Exit(), see issue 515.
* Restore 'm' in library list, see issue 538.
* Adjust to_bytes() behavior to match CPython.
* Provide wheels for CPython v3.14.
* Provide wheels for Linux and Windows ARM64.

gmpy2 2.2.1
-----------

* Bug fix: use C int instead of C char for some internal code.
* Bug fix: add xmpz.bit_count method.

gmpy2 2.2.0
-----------

gmpy2 2.2.0 is now available with support for Python 3.7 to 3.13.

* Support for thread-safe contexts and context methods has been improved.
* Interoperability with Cython extensions has been updated.
* Extensive improvements have been made to the build and testing processes.
* Many bug fixes.
* Extensive documentation cleanup.

Availability
------------

gmpy2 is available at https://pypi.python.org/pypi/gmpy2/

Documentation is available at https://gmpy2.readthedocs.io/en/latest/
