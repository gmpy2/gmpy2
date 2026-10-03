Integers (Advanced topics)
==========================

.. currentmodule:: gmpy2

gmpy2 provides access to an experimental integer type called `xmpz`. The
`xmpz` type is a mutable integer type. In-place operations (+=, //=,
etc.) modify the original object and do not create a new object. Instances of
`xmpz` cannot be used as dictionary keys.

.. doctest::

    >>> from gmpy2 import xmpz
    >>> a = xmpz(123)
    >>> b = a
    >>> a += 1
    >>> a
    xmpz(124)
    >>> b
    xmpz(124)

The ability to change an `xmpz` object in-place allows for efficient and
rapid bit manipulation.

Individual bits can be set or cleared:

.. doctest::

    >>> a[10]=1
    >>> a
    xmpz(1148)

Slice notation is supported. The bits referenced by a slice can be either 'read
from' or 'written to'. To clear a slice of bits, use a source value of 0. In
2s-complement format, 0 is represented by an arbitrary number of 0-bits. To set
a slice of bits, use a source value of ~0. The *tilde* operator inverts, or
complements the bits in an integer. (~0 is -1 so you can also use -1.) In
2s-complement format, -1 is represented by an arbitrary number of 1-bits.

For both `mpz` and `xmpz`, reading a slice of a negative integer with
nonnegative *start* and *stop* positions includes these sign bits beyond the
integer's bit-length. The selected bits are packed into an `mpz`, starting
with its least significant bit. For a positive *step*, *start* may be omitted
and defaults to zero. For a negative *step*, *start* must be explicit to
read beyond the bit-length. Slices with an omitted *stop*, a negative bound,
or an omitted *start* and negative *step* keep their finite-sequence behavior,
with bounds relative to the bit-length. Slicing nonnegative integers is
unchanged.

.. doctest::

    >>> xmpz(-1)[4:5]
    mpz(1)
    >>> xmpz(-2)[:8]
    mpz(254)
    >>> xmpz(-2)[4:0:-1]
    mpz(15)
    >>> xmpz(-2)[0:]
    mpz(2)

If a value for *stop* is specified in a slice assignment and the actual
bit-length of the `xmpz` is less than *stop*, then the destination
`xmpz` is logically padded with 0-bits to length *stop*.

.. doctest::

    >>> a=xmpz(0)
    >>> a[8:16] = ~0
    >>> bin(a)
    '0b1111111100000000'
    >>> a[4:12] = ~a[4:12]
    >>> bin(a)
    '0b1111000011110000'

Bits can be reversed:

.. doctest::

    >>> a = xmpz(1148)
    >>> bin(a)
    '0b10001111100'
    >>> a[::] = a[::-1]
    >>> bin(a)
    '0b111110001'

The `~xmpz.iter_bits()` method returns a generator that returns True or
False for each bit position. The methods `~xmpz.iter_clear()`, and
`~xmpz.iter_set()` return generators that return the bit positions that are
1 or 0. The methods support arguments *start* and *stop* that define the
beginning and ending bit positions that are used. To mimic the behavior of
slices. the bit positions checked include *start* but the last position checked
is *stop* - 1.

.. doctest::

    >>> a=xmpz(117)
    >>> bin(a)
    '0b1110101'
    >>> list(a.iter_bits())
    [True, False, True, False, True, True, True]
    >>> list(a.iter_clear())
    [1, 3]
    >>> list(a.iter_set())
    [0, 2, 4, 5, 6]
    >>> list(a.iter_bits(stop=12))
    [True, False, True, False, True, True, True, False, False, False, False, False]

The following program uses the Sieve of Eratosthenes to generate a list of
prime numbers.

.. code-block:: python

    import time
    import gmpy2

    def sieve(limit=1000000):
        '''Returns a generator that yields the prime numbers up to limit.'''

        # Increment by 1 to account for the fact that slices  do not include
        # the last index value but we do want to include the last value for
        # calculating a list of primes.
        sieve_limit = gmpy2.isqrt(limit) + 1
        limit += 1

        # Mark bit positions 0 and 1 as not prime.
        bitmap = gmpy2.xmpz(3)

        # Process 2 separately. This allows us to use p+p for the step size
        # when sieving the remaining primes.
        bitmap[4 : limit : 2] = -1

        # Sieve the remaining primes.
        for p in bitmap.iter_clear(3, sieve_limit):
            bitmap[p*p : limit : p+p] = -1

        return bitmap.iter_clear(2, limit)

    if __name__ == "__main__":
        start = time.time()
        result = list(sieve())
        print(time.time() - start)
        print(len(result))


The xmpz type
-------------

.. autoclass:: xmpz
   :special-members: __format__


Advanced Number Theory Functions
--------------------------------

The following functions are based on mpz_lucas.c and mpz_prp.c by David
Cleaver.

A good reference for probable prime testing is
http://www.pseudoprime.com/pseudo.html

.. autofunction:: is_bpsw_prp
.. autofunction:: is_euler_prp
.. autofunction:: is_extra_strong_lucas_prp
.. autofunction:: is_fermat_prp
.. autofunction:: is_fibonacci_prp
.. autofunction:: is_lucas_prp
.. autofunction:: is_selfridge_prp
.. autofunction:: is_strong_bpsw_prp
.. autofunction:: is_strong_lucas_prp
.. autofunction:: is_strong_prp
.. autofunction:: is_strong_selfridge_prp
.. autofunction:: lucasu
.. autofunction:: lucasu_mod
.. autofunction:: lucasv
.. autofunction:: lucasv_mod
