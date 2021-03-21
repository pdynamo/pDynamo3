"Orders of magnitude."

import math

from enum import Enum

#===================================================================================================================================
# . Definitions.
#===================================================================================================================================
# . Magnitudes.
class Magnitude ( Enum ):
    "An order of magnitude."
    Septillionth  = ( "quadrillionth", "septillionth" , "yocto", "y" , -24 )
    Sextillionth  = ( "trilliardth"  , "sextillionth" , "zepto", "z" , -21 )
    Quintillionth = ( "trillionth"   , "quintillionth", "atto" , "a" , -18 )
    Quadrillionth = ( "billiardth"   , "quadrillionth", "femto", "f" , -15 )
    Trillionth    = ( "billionth"    , "trillionth"   , "pico" , "p" , -12 )
    Billionth     = ( "milliardth"   , "billionth"    , "nano" , "n" ,  -9 )
    Millionth     = ( "millionth"    , "millionth"    , "micro", "mu",  -6 )
    Thousandth    = ( "thousandth"   , "thousandth"   , "milli", "m" ,  -3 )
    Hundredth     = ( "hundredth"    , "hundredth"    , "centi", "c" ,  -2 )
    Tenth         = ( "tenth"        , "tenth"        , "deci" , "d" ,  -1 )
    One           = ( "one"          , "one"          , ""     , ""  ,   0 )
    Ten           = ( "ten"          , "ten"          , "deca" , "da",   1 )
    Hundred       = ( "hundred"      , "hundred"      , "hecto", "h" ,   2 )
    Thousand      = ( "thousand"     , "thousand"     , "kilo" , "k" ,   3 )
    Million       = ( "million"      , "million"      , "mega" , "M" ,   6 )
    Billion       = ( "milliard"     , "billion"      , "giga" , "G" ,   9 )
    Trillion      = ( "billion"      , "trillion"     , "tera" , "T" ,  12 )
    Quadrillion   = ( "billiard"     , "quadrillion"  , "peta" , "P" ,  15 )
    Quintillion   = ( "trillion"     , "quintillion"  , "exa"  , "E" ,  18 )
    Sextillion    = ( "trilliard"    , "sextillion"   , "zetta", "Z" ,  21 )
    Septillion    = ( "quadrillion"  , "septillion"   , "yotta", "Y" ,  24 )

    def __init__ ( self, longName, shortName, prefix, symbol, power ):
        "Constructor."
        self.longName  = longName
        self.shortName = shortName
        self.prefix    = prefix
        self.symbol    = symbol
        self.power     = power

#===================================================================================================================================
# . Functions.
#===================================================================================================================================
def Magnitude_Adjust ( value, magnitude = Magnitude.One, powersOfThree = True ):
    """Change the magnitude of a value."""
    # . Power of value adjusting for existing magnitude.
    oldP  = magnitude.power
    power = math.floor ( math.log10 ( math.fabs ( value ) ) ) + oldP
    # . Find new magnitude.
    magnitude = None
    for new in Magnitude:
        if powersOfThree and ( new.power % 3 != 0 ): continue
        if new.power > power: break
        magnitude = new
    if magnitude is None: magnitude = Magnitude.Septillionth
    # . Adjust value.
    value *= math.pow ( 10.0, ( oldP - magnitude.power ) )
    return ( value, magnitude )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
