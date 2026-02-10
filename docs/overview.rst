Overview
========

.. currentmodule:: gmpy2

The `mpz` and `mpq` types support arbitrary precision integers and rationals
via the GMP library.  These types should be drop-in replacements for Python's
`int`'s and `~fractions.Fraction`'s, but are significantly faster for large
values.   The cutover point for performance varies, but can be as low as 20 to
40 digits.  All the special integer functions in the GMP are supported.

.. warning::

    gmpy2 can crash the Python interpreter in case of memory allocation
    failure.  To mitigate this feature of memory management in the GMP library,
    you should estimate the size of all results and prevent calculations that
    can exaust available memory.  In case you are interested only in
    the GMP bindings, you can use the
    `python-gmp <https://pypi.org/project/python-gmp/>`_ package.

The `mpfr` and `mpc` types provide support for correctly rounded multiple
precision real and complex arithmetic via the MPFR and MPC libraries.  The
`context` type is used to control precision, rounding modes, and exceptional
conditions.  For example, division by zero can either return an Infinity or
raise an exception.  It is possible to specify different precision and rounding
modes for both the real and imaginary components of an `mpc`.  The default
precision is 53 bits --- just same as for Python's `float` and `complex` types.

.. warning::

   The gmpy2 has support for :term:`free threading` builds of the CPython.  But
   this requires using thread-local storage for the MPFR data, such as as
   flags, the exponent range, the default precision, etc.  That is, the MPFR
   should be build with the ``--enable-thread-safe`` configure option.

Operator overloading is fully supported.  Conversion from native Python types is
optimized for performance.
