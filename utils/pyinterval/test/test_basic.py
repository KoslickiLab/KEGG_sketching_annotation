# Copyright (c) 2008-2016, Stefano Taschini <taschini@ieee.org>
# All rights reserved.
# See LICENSE for details.

import unittest
from interval import interval, fpu


class FpuTestCase(unittest.TestCase):

    def test_third(self):
        "Nearest rounding of 1/3 is downwards."
        from operator import truediv
        assert  1 / 3.0 == fpu.down(lambda: truediv(1.0,  3.0))
        assert  1 / 3.0 <  fpu.up  (lambda: truediv(1.0,  3.0))
        assert -1 / 3.0 == fpu.up  (lambda: truediv(1.0, -3.0))
        assert -1 / 3.0 >  fpu.down(lambda: truediv(1.0, -3.0))

    def test_fourth(self):
        " 1/4 is exact."
        from operator import truediv
        assert  1 / 4.0 == fpu.down(lambda: truediv(1.0,  4.0))
        assert  1 / 4.0 == fpu.up  (lambda: truediv(1.0,  4.0))
        assert -1 / 4.0 == fpu.up  (lambda: truediv(1.0, -4.0))
        assert -1 / 4.0 == fpu.down(lambda: truediv(1.0, -4.0))

    def test_fifth(self):
        "Nearest rounding of 1/5 is upwards."
        from operator import truediv
        assert  1 / 5.0 == fpu.up  (lambda: truediv(1.0,  5.0))
        assert  1 / 5.0 >  fpu.down(lambda: truediv(1.0,  5.0))
        assert -1 / 5.0 == fpu.down(lambda: truediv(1.0, -5.0))
        assert -1 / 5.0 <  fpu.up  (lambda: truediv(1.0, -5.0))

    def test_ieee754(self):
        "fpu.float respect ieee754 semantics."
        assert fpu.infinity + fpu.infinity == fpu.infinity
        assert fpu.isnan(fpu.nan)
        assert fpu.isnan(0.0 * fpu.infinity)
        assert fpu.isnan(fpu.infinity - fpu.infinity)

    def test_float_coercion(self):
        "Only real-number scalars should be able to coerce as fpu.float"
        self.assertRaises(Exception, lambda: float(1, 2))
        self.assertRaises(Exception, lambda: float((1, 2)))
        self.assertRaises(Exception, lambda: float([1, 2]))
        self.assertRaises(Exception, lambda: float('a'))
        self.assertRaises(Exception, lambda: float(1 + 1j))

    def test_min(self):
        "Verify corner cases with nan, -inf, +inf"
        assert fpu.min((1.0, 2.0))           == 1.0
        assert fpu.min((1.0, fpu.infinity))  == 1.0
        assert fpu.min((1.0, -fpu.infinity)) == -fpu.infinity
        assert fpu.isnan(fpu.min((1.0, -fpu.nan)))

    def test_max(self):
        "Verify corner cases with nan, -inf, +inf"
        assert fpu.max((1.0, 2.0))           == 2.0
        assert fpu.max((1.0, fpu.infinity))  == fpu.infinity
        assert fpu.max((1.0, -fpu.infinity)) == 1.0
        assert fpu.isnan(fpu.max((1.0, fpu.nan)))

    def test_power(self):
        x = 1 / 3.0
        # The cube of one third should depend on the rounding mode
        assert fpu.down(lambda: x * x * x) < fpu.up(lambda: x * x * x)
        # But using the built-in power operator, it doesn't necessarily do it
        # fpu.down(lambda: x ** 3) < fpu.up(lambda: x ** 3))
        # So we define an integer power methods that does
        assert fpu.power_rd( x, 3) < fpu.power_ru( x, 3)
        assert fpu.power_rd(-x, 3) < fpu.power_ru(-x, 3)
        assert fpu.power_rd( x, 4) < fpu.power_ru( x, 4)
        assert fpu.power_rd(-x, 4) < fpu.power_ru(-x, 4)

        assert (fpu.down(lambda: x * x * x), fpu.up(lambda: x * x * x)) == (fpu.power_rd(x, 3), fpu.power_ru(x, 3))


class ModuleTestCase(unittest.TestCase):

    def test_namespace(self):
        import interval
        assert [x for x in dir(interval) if not x.startswith('__')] == ['fpu', 'imath', 'inf', 'interval']


class IntervalTestCase(unittest.TestCase):

    def test_trivial_constructor(self):
        assert interval[1]              == ((1, 1),)
        assert interval(1)              == ((1, 1),)
        assert interval[1, 2]           == ((1, 2),)
        assert interval(1, 2)           == ((1, 1), (2, 2))
        assert interval([1, 2], [3, 4]) == ((1, 2), (3, 4))
        assert interval([1, 2])         == interval(interval([1, 2]))

    def test_nan_constructor(self):
        assert interval[2, fpu.nan]    == ((-fpu.infinity, fpu.infinity),)
        assert interval[2, fpu.nan]    == ((-fpu.infinity, fpu.infinity),)
        assert interval(2, fpu.nan, 9) == ((-fpu.infinity, fpu.infinity),)

    def test_failing_constructor(self):
        self.assertRaises(interval.ComponentError, lambda: interval[1, [2, 3]])
        self.assertRaises(interval.ComponentError, lambda: interval[1, 2, 3])
        self.assertRaises(interval.ComponentError, lambda: interval(0, [1, 2, 3]))
        self.assertRaises(interval.ComponentError, lambda: interval(0, [1, [2, 3]]))
        self.assertRaises(interval.ComponentError, lambda: interval['a', 1])

    def test_canonical_constructor(self):
        assert interval([1, 3], [4, 6], [2, 5], 9) == ((1, 6), (9, 9))
        assert interval[ 2 ** (52 + 1) - 1]        == interval[9007199254740991.0]
        assert interval[ 2 ** (52 + 1) + 1]        == interval[ 4503599627370496 * 2.0,  4503599627370497 * 2.0]
        assert interval[-2 ** (52 + 1) + 1]        == interval[-9007199254740991.0]
        assert interval[-2 ** (52 + 1) - 1]        == interval[-4503599627370497 * 2.0, -4503599627370496 * 2.0]
        assert interval[ 2 ** (52 + 2) + 1]        == interval[ 4503599627370496 * 4.0,  4503599627370497 * 4.0]
        assert interval[ 2 ** (52 + 2) + 2]        == interval[ 4503599627370496 * 4.0,  4503599627370497 * 4.0]
        assert interval[ 2 ** (52 + 2) + 3]        == interval[ 4503599627370496 * 4.0,  4503599627370497 * 4.0]
        assert interval[-2 ** (52 + 2) - 1]        == interval[-4503599627370497 * 4.0, -4503599627370496 * 4.0]
        assert interval[-2 ** (52 + 2) - 2]        == interval[-4503599627370497 * 4.0, -4503599627370496 * 4.0]
        assert interval[-2 ** (52 + 2) - 3]        == interval[-4503599627370497 * 4.0, -4503599627370496 * 4.0]

    def test_unary(self):
        assert interval[1, 2]   == +interval[1, 2]
        assert interval[-2, -1] == -interval[1, 2]

    def test_sum(self):
        assert interval[-fpu.infinity, +fpu.infinity]                            == interval[-fpu.infinity] + interval[fpu.infinity]
        assert interval[4, 6]                                                    == interval[1, 2] + interval[3, 4]
        assert interval[3, fpu.infinity]                                         == interval[1, fpu.infinity] + interval[2]
        assert interval[-fpu.infinity, +fpu.infinity]                            == interval[-fpu.infinity, -1] + interval[2, +fpu.infinity]
        assert interval[-fpu.infinity, +fpu.infinity]                            == interval[-fpu.infinity] + interval[8, +fpu.infinity]
        assert interval([1, 2], [10, fpu.infinity]) + interval([1, 9], [-2, -1]) == interval([-1, 1], [2, fpu.infinity])
        assert interval[1, 9] + interval([1, 2], [10, fpu.infinity])             == interval[2, fpu.infinity]

    def test_sum_coercion(self):
        assert interval[1, 2] + 2 == interval[3, 4]
        self.assertRaises(TypeError, lambda: interval[1, 2] + 1j)
        assert 1 + interval[4, 5] == interval[5, 6]
        self.assertRaises(TypeError, lambda: (1, 2) + interval[1, 2])
        assert fpu.infinity + interval[4, 5] == interval[fpu.infinity]

    def test_sub(self):
        assert interval[1, 2] - interval[3, 4] == interval[-3.0, -1.0]
        assert interval[1, 2] - 0.5            == interval[0.5, 1.5]
        assert 1.5 - interval[1, 2]            == interval[-0.5, 0.5]

    def test_mul(self):
        assert interval[-fpu.infinity, +fpu.infinity]      == fpu.infinity * interval[0]
        assert interval[+fpu.infinity]                     == interval[+fpu.infinity] * interval[3]
        assert interval[-8, +10]                           == interval[1, 2] * interval[-4, 5]
        assert interval[3, 8]                              == interval[1, 2] * interval[3, 4]
        assert interval[-fpu.infinity, +fpu.infinity]      == interval[0, 1] * interval[2, +fpu.infinity]
        assert interval[2, fpu.infinity]                   == interval[-fpu.infinity, -2] * interval[-fpu.infinity, -1]
        assert interval([1, 2], [3, 4]) * interval[0.5, 2] == interval[0.5, 8]
        assert interval[1, 2] * 2                          == interval[2, 4]

    def test_inverse(self):
        assert interval[0.5, 1]                                     == interval[1, 2].inverse()
        assert interval[-1, -0.5]                                   == (-interval[1, 2]).inverse()
        assert interval([-fpu.infinity, -1], [0.5, +fpu.infinity])  == interval[-1, 2].inverse()
        assert interval(-fpu.infinity, [1, +fpu.infinity])          == interval[0, 1].inverse()
        assert interval([-fpu.infinity, -2.0], [0.0, fpu.infinity]) == interval([-0.5, 0.5], [0.2, fpu.infinity]).inverse()

    def test_division(self):
        assert interval[-fpu.infinity, fpu.infinity] == interval[0, 1] / interval[0, 1]
        assert interval[0.5]                         == interval[1] / 2
        assert interval[0.5]                         == 1 / interval[2]

    def test_power(self):
        self.assertRaises(TypeError, lambda: interval[1, 2] ** (1.3))
        assert (-interval[1, 2]).inverse()           == (-interval[1, 2]) ** -1
        assert interval[0, 4]                        == interval[-1, 2] ** 2
        assert interval[-27, 8]                      == interval[-3, 2] ** 3
        assert interval[-1, 2]                       == (interval[-1, 2] ** -1) ** -1
        assert interval([-0.38712442133802405]) ** 3 == interval([-0.058016524353106828, -0.058016524353106808])

        from operator import truediv
        assert (
            interval[
                fpu.down(lambda: truediv(1, 3.0) * truediv(1, 3.0)),
                fpu.up  (lambda: truediv(1, 3.0) * truediv(1, 3.0))] ==
            (interval[1] / 3.0) ** 2)

        assert (
            interval[
                fpu.down(lambda: truediv(1, 3.0) * truediv(1, 3.0) * truediv(1, 3.0)),
                fpu.up(lambda: truediv(1, 3.0) * truediv(1, 3.0) * truediv(1, 3.0))] ==
            (interval[1] / 3.0) ** 3)

    def test_format(self):
        for x in interval[1], interval[1, 2], interval([1, 2], [3, 4]), interval[2, 2.0000000000000004]:
            assert x == eval(repr(x))
        x = interval([1, 2])
        assert str(x) == repr(x)

    def test_intersection(self):
        assert interval[1, 2] & interval[0, 3]             == interval[1, 2]
        assert interval[1.1, 1.9] & interval[1.3, 2.5]     == interval[1.3, 1.9]
        assert interval[1.1, 1.9] & interval[0.3, 0.7]     == interval()
        assert interval([1, 3], [4, 5]) & interval[2]      == interval[2]
        assert interval([1, 3], [4, 5]) & interval(2, 4.5) == interval(2, 4.5)
        assert interval[1, 2] & 1.2                        == interval(1.2)
        assert 2.1 & interval[1, 2]                        == interval()

    def test_union(self):
        assert interval([1, 6], 9)  == interval([1, 3], [4, 6]) | interval([2, 5], 9)
        assert interval[1, 2] | 2.1 == interval([1, 2], 2.1)
        assert 2.1 | interval[1, 2] == interval([1, 2], 2.1)
        self.assertRaises(TypeError, lambda: interval[1, 2] | 1j)

    def test_abs(self):
        assert interval([0, 3])     == abs(interval[-3, 2])
        assert interval([1, 6], 9)  == abs(interval([-9], [-5, -2], [1, 3], [4, 6]))
        assert interval([1, 6], 9)  == abs(interval([9], [2, 5], [-3, -1], [-6, -4]))

    def test_hull(self):
        assert interval([1, 9]) == interval.hull((interval([1, 3], [4, 6]), interval([2, 5], 9)))

    def test_inclusion(self):
        def verify_in(x, y):
            assert x in y
            assert x & y == interval(x)

        verify_in(1.5, interval[1, 2])
        verify_in(1,   interval[1, 2])
        verify_in(2,   interval[1, 2])
        verify_in(interval[1, 2],   interval[1, 2])
        verify_in(interval[1.1, 2], interval[1, 2])
        verify_in(interval[1, 1.8], interval[1, 2])
        verify_in(interval([1.1, 2.2], [3.3, 4.4]), interval(-1, [0, 2.5], [3, 5], [7, 9]))

        def verify_out(x, y):
            self.assertFalse(x in y)
            self.assertNotEqual(x & y, x)

        verify_out(0, interval[1, 2])
        verify_out(4, interval[1, 2])
        verify_out(interval[1, 3], interval[2, 4])
        verify_out(interval(1, 3), interval(2, 4))

    def test_extrema(self):
        assert interval(1, [2, 3], 4).extrema == interval(1, 2, 3, 4)

    def test_midpoint(self):
        assert interval[0, 4].midpoint == interval[2]
        assert interval(-1, 1, 4)      == interval(-1, [0, 2], [3, 5]).midpoint

    def test_pickle_copy(self):
        # https://github.com/taschini/pyinterval/issues/5
        import copy
        import pickle
        a = interval([-3, -2], [0, 1])
        assert a == pickle.loads(pickle.dumps(a, -1))
        assert a == copy.copy(a)
        assert a == copy.deepcopy(a)


class NewtonTestCase(unittest.TestCase):

    def test_opts(self):
        self.assertRaises(TypeError, lambda: interval(0, 1).newton(None, None, nonexisting=True))

    def test_cubic(self):
        assert interval[-2, 2].newton(lambda x: x ** 3 - x, lambda x: 3 * x ** 2 - 1)      == interval(-1, 0, 1)
        assert interval[-5, 5].newton(lambda x: x ** 3 + x - 10, lambda x: 3 * x ** 2 + 1) == interval[2]
        assert interval[-5, 5].newton(lambda x: x ** 3 + x - 15, lambda x: 3 * x ** 2 + 1) == interval[5249383869325653 * 2.0 ** -51, 5249383869325655 * 2.0 ** -51]
        # The sharpest result would be with 5249383869325654 * 2.0 ** -51 as sup.

    def test_sqrt2(self):
        import math
        f, p = lambda x: x**2 - 2, lambda x: 2 * x
        u, v = 6369051672525772 * 2.0 ** -52, 6369051672525773 * 2.0 ** -52
        s = interval[u, v]
        assert v          == math.sqrt(2)
        assert s          == interval[0.1, 5].newton(f, p)
        assert s          == interval[0, 2].newton(f, p)
        assert s          == interval[-1, 10].newton(f, p)
        assert interval() == interval[2, 5].newton(f, p)
        assert -s         == interval[-5, 0].newton(f, p)
        assert -s | s     == interval[-5, +5].newton(f, p)

        # Failure to converge in only three iterations:
        messages = []
        assert interval() == interval[0, 2].newton(
            f, p, maxiter=3,
            tracer_cb=lambda tag, interval: messages.append((tag, interval)))
        assert messages == [
            ('branch' , interval[0.0, 2.0]),
            ('step'   , interval[1.25, 2.0]),
            ('step'   , interval[1.36875, 1.46484375]),
            ('step'   , interval[1.4141253188320488, 1.4143005729166669]),
            ('abandon', interval[1.4141253188320488, 1.4143005729166669])
        ]
