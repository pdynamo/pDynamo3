"""Example 0."""

import math

# . Define a squaring function.
def Square ( x ):
    return x**2

# . Create a list of integers.
values = range ( 10 )

# . Loop over the integers.
for i in values:
    x = float ( i )
    print ( "{:5d}{:10.5f}{:10.5f}{:10.5f}".format ( i, x, math.sqrt ( x ), Square ( x ) ) )
