import platform
import sys

import pytest

import gmpy2


def test_misc():
    assert gmpy2.version().startswith('2.3')
    assert gmpy2.mp_limbsize() in (32,64)
    assert '5.0.0' <= gmpy2.mp_version()
    assert gmpy2.mpfr_version() and gmpy2.mpfr_version().startswith('MPFR')
    assert gmpy2.mpfr_version() and '3.1.0' <= gmpy2.mpfr_version().split()[1]
    assert gmpy2.mpc_version() and gmpy2.mpc_version().startswith('MPC')
    assert gmpy2.mpc_version() and '1.0' <= gmpy2.mpc_version().split()[1]
    assert gmpy2.license() == ('The GMPY2 source code is licensed under LGPL '
                               '3 or later. The supported versions of the GMP, '
                               'MPFR, and MPC libraries are also licensed '
                               'under LGPL 3 or later.')


@pytest.mark.skipif(platform.python_implementation() == "PyPy",
                    reason="sys.getsizeof raises TypeError")
def test_sizeof():
    assert sys.getsizeof(gmpy2.mpz(10)) > 0
    assert sys.getsizeof(gmpy2.mpfr('1.0')) > 0
