import sys

import pytest

from gmpy2 import mpz, xmpz


def bit_slice(n, positions):
    return sum(((n >> k) & 1) << j for j, k in enumerate(positions))


@pytest.mark.parametrize('integer_type', [mpz, xmpz])
@pytest.mark.parametrize('n', [-1, -2, -10, -(1 << 65), -((1 << 130) + 17)])
@pytest.mark.parametrize('start, stop, step', [
    (0, 1, 1), (4, 5, 1), (0, 8, 1), (4, 10, 1),
    (60, 75, 1), (126, 138, 1), (0, 10, 2), (4, 16, 3),
    (10, 2, -1), (10, 2, -3), (75, 60, -2), (138, 126, -3),
    (5, 5, 1), (6, 3, 2), (3, 6, -2),
])
def test_negative_bounded_bit_slice(integer_type, n, start, stop, step):
    value = integer_type(n)
    result = value[start:stop:step]
    expected = bit_slice(n, range(start, stop, step))
    assert result == expected
    assert type(result) is mpz
    assert int(value) == n
    if step == 1 and stop >= start:
        assert result == (n >> start) & ((1 << (stop - start)) - 1)


@pytest.mark.parametrize('integer_type', [mpz, xmpz])
@pytest.mark.parametrize('n', [-1, -2, -10])
@pytest.mark.parametrize('stop, step', [(0, 1), (8, 1), (11, 3)])
def test_negative_bit_slice_default_start(integer_type, n, stop, step):
    assert integer_type(n)[:stop:step] == bit_slice(n, range(0, stop, step))


@pytest.mark.parametrize('integer_type', [mpz, xmpz])
@pytest.mark.parametrize('n', [-1, -2, -10, -((1 << 130) + 17), 0, 1, 10,
                              (1 << 130) + 17])
@pytest.mark.parametrize('selection', [
    slice(None), slice(2, None), slice(None, None, -1),
    slice(140, None, -2), slice(-3, 140), slice(140, -1, -2),
    slice(None, 0, -1), slice(-1, None, -1),
])
def test_bit_slice_legacy_bounds(integer_type, n, selection):
    length = max(1, n.bit_length())
    positions = range(*selection.indices(length))
    assert integer_type(n)[selection] == bit_slice(n, positions)


@pytest.mark.parametrize('integer_type', [mpz, xmpz])
@pytest.mark.parametrize('n', [0, 1, 10, (1 << 130) + 17])
@pytest.mark.parametrize('selection', [
    slice(0, 140), slice(4, 140, 3), slice(140, 2, -3),
    slice(None, 140, 2),
])
def test_nonnegative_bit_slice_legacy_bounds(integer_type, n, selection):
    length = max(1, n.bit_length())
    positions = range(*selection.indices(length))
    assert integer_type(n)[selection] == bit_slice(n, positions)


@pytest.mark.parametrize('integer_type', [mpz, xmpz])
def test_bit_slice_index_protocol(integer_type):
    class Index:
        def __init__(self, value):
            self.value = value
            self.calls = 0

        def __index__(self):
            self.calls += 1
            return self.value

    start, stop, step = Index(4), Index(9), Index(2)
    assert integer_type(-1)[start:stop:step] == 7
    assert (start.calls, stop.calls, step.calls) == (1, 1, 1)


@pytest.mark.parametrize('integer_type', [mpz, xmpz])
@pytest.mark.parametrize('selection, error', [
    (slice(1, 5, 0), ValueError), (slice('x', 5), TypeError),
    (slice(1, 'x'), TypeError), (slice(1, 5, 'x'), TypeError),
])
def test_bit_slice_errors(integer_type, selection, error):
    with pytest.raises(error):
        integer_type(-1)[selection]


@pytest.mark.parametrize('integer_type', [mpz, xmpz])
def test_bit_slice_large_step(integer_type):
    # These slices select at most one bit, even when the step is clipped.
    assert integer_type(-1)[4:5:sys.maxsize * 2] == 1
    assert integer_type(-1)[4:3:-sys.maxsize * 2] == 1


@pytest.mark.parametrize('integer_type', [mpz, xmpz])
def test_bit_index_unchanged(integer_type):
    value = integer_type(-10)
    assert value[8] == 1
    assert value[-1] == 0
    with pytest.raises(IndexError):
        value[sys.maxsize + 1]
    with pytest.raises(TypeError):
        value['x']


def test_xmpz_slice_assignment_unchanged():
    value = xmpz(0)
    value[8:16] = -1
    assert value == 0xff00
    value[4:12:2] = 0
    assert value == 0xfa00
    value = xmpz(-2)
    value[:] = 0
    assert value == -4
