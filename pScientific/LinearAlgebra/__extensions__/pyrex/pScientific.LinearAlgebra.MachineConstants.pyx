from pCore.CPrimitiveTypes cimport CReal
from pCore                  import logFile       , \
                                   LogFileActive

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "f2clapack.h":

    # . Functions.
    cdef CReal dlamch_ ( char *cMach ) 

#===================================================================================================================================
# . Functions.
#===================================================================================================================================
def MachineConstants ( log = logFile ):
    """Return machine constants."""
    # . dlamch returns a float even though some of these are integers.
    cdef CReal value
    constants = []
    for ( c, label ) in ( ( "E" , "Relative Machine Precision" ) ,
                          ( "S" , "Safe Overflow Minimum"      ) ,
                          ( "B" , "Machine Base"               ) ,
                          ( "P" , "Precision"                  ) ,
                          ( "N" , "Mantissa Digits"            ) ,
                          ( "R" , "Rounding Value"             ) ,
                          ( "M" , "Minimum Underflow Exponent" ) ,
                          ( "U" , "Underflow Threshold"        ) ,
                          ( "L" , "Largest Overflow Exponent"  ) ,
                          ( "O" , "Overflow Threshold"         ) ):
        constants.append ( ( label, "{:g}".format ( dlamch_ ( c.encode ( u"ascii" ) ) ) ) )
    if LogFileActive ( log ): log.SummaryOfItems ( constants, title = "Machine Constants" )
    return constants

#===================================================================================================================================
