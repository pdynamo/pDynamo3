"""Scripts for MNDO parameter handling."""

import glob, math, os

from  pCore          import logFile                    , \
                            LogFileActive              , \
                            YAMLMappingFile_FromObject , \
                            YAMLMappingFile_ToObject   , \
                            YAMLPickleFileExtension
from  pScientific    import PeriodicTable
from .MNDOParameters import MNDOParameters

#===================================================================================================================================
# . Parameters for conversion between formats.
#===================================================================================================================================
# . Defined Hamiltonians.
_HAMILTONIANS = ( "am1", "mndo", "pddgmndo", "pddgpm3", "pmo1", "pm3", "pm6", "rm1" )

# . Hamiltonian specific data.
_HAMILTONIANMINHPP    = { "am1"       : 0.01  ,
                          "mndo"      : 0.01  ,
                          "pddgmndo"  : 0.01  ,
                          "pddgpm3"   : 0.008 ,
                          "pmo1"      : 0.1   ,
                          "pm3"       : 0.1   ,
                          "pm6"       : 0.1   ,
                          "rm1"       : 0.01  }
_HAMILTONIANTOEXCLUDE = { "am1"       : None                ,
                          "mndo"      : None                ,
                          "pddgmndo"  : set ( [ "eisol" ] ) ,
                          "pddgpm3"   : set ( [ "eisol" ] ) ,
                          "pmo1"      : None                ,
                          "pm3"       : None                ,
                          "pm6"       : None                ,
                          "rm1"       : None                }

# . File extensions.
_TEXTEXTENSION = ".txt"
_YAMLEXTENSION = YAMLPickleFileExtension
_YAMLTag       = "!MNDOParameters"

# . Directories.
_SCRATCH     = os.getenv ( "PDYNAMO3_SCRATCH" )
_TESTPATH    = os.path.join ( _SCRATCH, "mndoParameters" )
_TEXTPATH    = "text"
_YAMLPATHIN  = os.path.join ( os.getenv ( "PDYNAMO3_PARAMETERS" ), "mndoParameters" )
_YAMLPATHOUT = "yaml"

#===================================================================================================================================
# . Parameters indicating parameter names.
#===================================================================================================================================
# . Attributes.
_INTEGERATTRIBUTES = ( "atomicNumber", "orbitals", "iii", "iiid", "ir016", "ir066", "ir244", "ir266", "ir466", "qnd", "qnp", "qns" )
_REALATTRIBUTES    = ( "ad"   ,"alp"  ,"am"   ,"aq"   ,"betad","betap","betas","dd"   ,"eheat","eisol","f0sd" ,"gphot","gpp"  ,"gp2"  ,"gsp"  ,
                       "gss"  ,"g2sd" ,"hsp"  ,"pcore","qq"   ,"udd"  ,"upp"  ,"uss"  ,"zcore","zetad","zetap","zetas","zdn"  ,"zpn"  ,"zsn"  )
# . Array item names.
_AM1PM3GNAMES = ( "fn1", "fn2", "fn3" )
_PDDGNAMES    = ( "pddgc", "pddge"    )

# . Maximum deviation.
_MAXIMUMDEVIATION = 1.0e-5

# . Atomic units to eV.
_MOPAC_HARTREES_TO_ELECTRONVOLTS = 27.2113834

# . Units.
_UNITS = { "ad"        : "atomic"        ,
           "alp"       : "A^-1"          ,
           "am"        : "atomic"        ,
           "aq"        : "atomic"        ,
           "betad"     : "eV"            ,
           "betap"     : "eV"            ,
           "betas"     : "eV"            ,
           "dd"        : "atomic"        ,
           "diatomicA" : "A^-1"          ,
           "diatomicX" : "dimensionless" ,
           "eheat"     : "kcal/mole"     ,
           "eisol"     : "eV"            ,
           "f0sd"      : "eV"            ,
           "gphot"     : "dimensionless" ,
           "gpp"       : "eV"            ,
           "gp2"       : "eV"            ,
           "gsp"       : "eV"            ,
           "gss"       : "eV"            ,
           "g2sd"      : "eV"            ,
           "hsp"       : "eV"            ,
           "pcore"     : "atomic"        ,
           "qq"        : "atomic"        ,
           "udd"       : "eV"            ,
           "upp"       : "eV"            ,
           "uss"       : "eV"            ,
           "zcore"     : "atomic"        ,
           "zetad"     : "atomic"        ,
           "zetap"     : "atomic"        ,
           "zetas"     : "atomic"        ,
           "zdn"       : "atomic"        ,
           "zpn"       : "atomic"        ,
           "zsn"       : "atomic"        ,
           "fn1"       : "A_eV"          ,
           "fn2"       : "A^-2"          ,
           "fn3"       : "A"             ,
           "pddgc"     : "eV"            ,
           "pddge"     : "A"             }

# . Understood names.
_UNDERSTOODNAMES = set ( _INTEGERATTRIBUTES ) | set ( _AM1PM3GNAMES ) | set ( _PDDGNAMES ) | set ( _UNITS.keys ( ) )
for ( tags, number ) in ( ( _AM1PM3GNAMES, 4 ), ( _PDDGNAMES, 2 ) ):
    for tag in tags:
        for i in range ( number ): _UNDERSTOODNAMES.add ( tag + repr ( i ) )

#===================================================================================================================================
# . PM6 specific parameters.
#===================================================================================================================================
# . Extra values for rare gas atoms.
# . EISOL.
_EXTRAEISOL = { 10 : -272.229723, \
                18 : -264.225777, \
                36 : -259.001706  }

# . AD/DD.
_EXTRAADDD = {  2: ( 0.5 / 0.86864043, 0.24758191 ), \
               10: ( 0.5 / 0.47206084, 0.25922907 ), \
               18: ( 0.5 / 0.37835832, 0.33819095 ), \
               36: ( 0.5 / 0.70448068, 0.17977932 )  }

# . AQ/QQ.
_EXTRAAQQQ = { 2 : ( 0.5 / 0.51184367, 0.29953658 / math.sqrt ( 2.0 ) ) }

#===================================================================================================================================
# . Elemental parameter data that is more or less common to all MNDO methods.
#===================================================================================================================================
# . Class.
class UniversalParameters:
    """Class for universal parameters."""

    _attributable = { "atomicNumber" :   0, \
                      "eheat"        : 0.0, \
                      "iod"          :   0, \
                      "iop"          :   0, \
                      "ios"          :   0, \
                      "mainGroup"    : True }

    def __init__ ( self, **options ):
        """Constructor."""
        for ( key, value ) in self.__class__._attributable.items ( ): setattr ( self, key, value )
        for ( key, value ) in options.items ( )                     : setattr ( self, key, value )

#-----------------------------------------------------------------------------------------------------------------------------------
# . The parameter values.
#-----------------------------------------------------------------------------------------------------------------------------------
parameters      = [ None for i in range ( 108 ) ]
parameters[0]   = UniversalParameters ( atomicNumber = 0 )
parameters[1]   = UniversalParameters ( atomicNumber =   1, eheat =   52.102, iod = 0 , iop = 0, ios = 1, mainGroup =  True )
parameters[2]   = UniversalParameters ( atomicNumber =   2, eheat =    0.000, iod = 0 , iop = 0, ios = 2, mainGroup =  True )
parameters[3]   = UniversalParameters ( atomicNumber =   3, eheat =   38.410, iod = 0 , iop = 0, ios = 1, mainGroup =  True )
parameters[4]   = UniversalParameters ( atomicNumber =   4, eheat =   76.960, iod = 0 , iop = 0, ios = 2, mainGroup =  True )
parameters[5]   = UniversalParameters ( atomicNumber =   5, eheat =  135.700, iod = 0 , iop = 1, ios = 2, mainGroup =  True )
parameters[6]   = UniversalParameters ( atomicNumber =   6, eheat =  170.890, iod = 0 , iop = 2, ios = 2, mainGroup =  True )
parameters[7]   = UniversalParameters ( atomicNumber =   7, eheat =  113.000, iod = 0 , iop = 3, ios = 2, mainGroup =  True )
parameters[8]   = UniversalParameters ( atomicNumber =   8, eheat =   59.559, iod = 0 , iop = 4, ios = 2, mainGroup =  True )
parameters[9]   = UniversalParameters ( atomicNumber =   9, eheat =   18.890, iod = 0 , iop = 5, ios = 2, mainGroup =  True )
parameters[10]  = UniversalParameters ( atomicNumber =  10, eheat =    0.000, iod = 0 , iop = 6, ios = 0, mainGroup =  True )
parameters[11]  = UniversalParameters ( atomicNumber =  11, eheat =   25.650, iod = 0 , iop = 0, ios = 1, mainGroup =  True )
parameters[12]  = UniversalParameters ( atomicNumber =  12, eheat =   35.000, iod = 0 , iop = 0, ios = 2, mainGroup =  True )
parameters[13]  = UniversalParameters ( atomicNumber =  13, eheat =   79.490, iod = 0 , iop = 1, ios = 2, mainGroup =  True )
parameters[14]  = UniversalParameters ( atomicNumber =  14, eheat =  108.390, iod = 0 , iop = 2, ios = 2, mainGroup =  True )
parameters[15]  = UniversalParameters ( atomicNumber =  15, eheat =   75.570, iod = 0 , iop = 3, ios = 2, mainGroup =  True )
parameters[16]  = UniversalParameters ( atomicNumber =  16, eheat =   66.400, iod = 0 , iop = 4, ios = 2, mainGroup =  True )
parameters[17]  = UniversalParameters ( atomicNumber =  17, eheat =   28.990, iod = 0 , iop = 5, ios = 2, mainGroup =  True )
parameters[18]  = UniversalParameters ( atomicNumber =  18, eheat =    0.000, iod = 0 , iop = 6, ios = 0, mainGroup =  True )
parameters[19]  = UniversalParameters ( atomicNumber =  19, eheat =   21.420, iod = 0 , iop = 0, ios = 1, mainGroup =  True )
parameters[20]  = UniversalParameters ( atomicNumber =  20, eheat =   42.600, iod = 0 , iop = 0, ios = 2, mainGroup =  True )
parameters[21]  = UniversalParameters ( atomicNumber =  21, eheat =   90.300, iod = 1 , iop = 0, ios = 2, mainGroup = False ) # . Sc - Cu.
parameters[22]  = UniversalParameters ( atomicNumber =  22, eheat =  112.300, iod = 2 , iop = 0, ios = 2, mainGroup = False )
parameters[23]  = UniversalParameters ( atomicNumber =  23, eheat =  122.900, iod = 3 , iop = 0, ios = 2, mainGroup = False )
parameters[24]  = UniversalParameters ( atomicNumber =  24, eheat =   95.000, iod = 5 , iop = 0, ios = 1, mainGroup = False )
parameters[25]  = UniversalParameters ( atomicNumber =  25, eheat =   67.700, iod = 5 , iop = 0, ios = 2, mainGroup = False )
parameters[26]  = UniversalParameters ( atomicNumber =  26, eheat =   99.300, iod = 6 , iop = 0, ios = 2, mainGroup = False )
parameters[27]  = UniversalParameters ( atomicNumber =  27, eheat =  102.400, iod = 7 , iop = 0, ios = 2, mainGroup = False )
parameters[28]  = UniversalParameters ( atomicNumber =  28, eheat =  102.800, iod = 8 , iop = 0, ios = 2, mainGroup = False )
parameters[29]  = UniversalParameters ( atomicNumber =  29, eheat =   80.700, iod = 10, iop = 0, ios = 1, mainGroup = False )
parameters[30]  = UniversalParameters ( atomicNumber =  30, eheat =   31.170, iod = 0 , iop = 0, ios = 2, mainGroup =  True )
parameters[31]  = UniversalParameters ( atomicNumber =  31, eheat =   65.400, iod = 0 , iop = 1, ios = 2, mainGroup =  True )
parameters[32]  = UniversalParameters ( atomicNumber =  32, eheat =   89.500, iod = 0 , iop = 2, ios = 2, mainGroup =  True )
parameters[33]  = UniversalParameters ( atomicNumber =  33, eheat =   72.300, iod = 0 , iop = 3, ios = 2, mainGroup =  True )
parameters[34]  = UniversalParameters ( atomicNumber =  34, eheat =   54.300, iod = 0 , iop = 4, ios = 2, mainGroup =  True )
parameters[35]  = UniversalParameters ( atomicNumber =  35, eheat =   26.740, iod = 0 , iop = 5, ios = 2, mainGroup =  True )
parameters[36]  = UniversalParameters ( atomicNumber =  36, eheat =    0.000, iod = 0 , iop = 6, ios = 0, mainGroup =  True )
parameters[37]  = UniversalParameters ( atomicNumber =  37, eheat =   19.600, iod = 0 , iop = 0, ios = 1, mainGroup =  True )
parameters[38]  = UniversalParameters ( atomicNumber =  38, eheat =   39.100, iod = 0 , iop = 0, ios = 2, mainGroup =  True )
parameters[39]  = UniversalParameters ( atomicNumber =  39, eheat =  101.500, iod = 1 , iop = 0, ios = 2, mainGroup = False ) # . Y - Ag.
parameters[40]  = UniversalParameters ( atomicNumber =  40, eheat =  145.500, iod = 2 , iop = 0, ios = 2, mainGroup = False )
parameters[41]  = UniversalParameters ( atomicNumber =  41, eheat =  172.400, iod = 4 , iop = 0, ios = 1, mainGroup = False )
parameters[42]  = UniversalParameters ( atomicNumber =  42, eheat =  157.300, iod = 5 , iop = 0, ios = 1, mainGroup = False )
parameters[43]  = UniversalParameters ( atomicNumber =  43, eheat =  162.000, iod = 5 , iop = 0, ios = 2, mainGroup = False )
parameters[44]  = UniversalParameters ( atomicNumber =  44, eheat =  155.500, iod = 7 , iop = 0, ios = 1, mainGroup = False )
parameters[45]  = UniversalParameters ( atomicNumber =  45, eheat =  133.000, iod = 8 , iop = 0, ios = 1, mainGroup = False )
parameters[46]  = UniversalParameters ( atomicNumber =  46, eheat =   90.000, iod = 10, iop = 0, ios = 0, mainGroup = False )
parameters[47]  = UniversalParameters ( atomicNumber =  47, eheat =   68.100, iod = 10, iop = 0, ios = 1, mainGroup = False )
parameters[48]  = UniversalParameters ( atomicNumber =  48, eheat =   26.720, iod = 0 , iop = 0, ios = 2, mainGroup =  True )
parameters[49]  = UniversalParameters ( atomicNumber =  49, eheat =   58.000, iod = 0 , iop = 1, ios = 2, mainGroup =  True )
parameters[50]  = UniversalParameters ( atomicNumber =  50, eheat =   72.200, iod = 0 , iop = 2, ios = 2, mainGroup =  True )
parameters[51]  = UniversalParameters ( atomicNumber =  51, eheat =   63.200, iod = 0 , iop = 3, ios = 2, mainGroup =  True )
parameters[52]  = UniversalParameters ( atomicNumber =  52, eheat =   47.000, iod = 0 , iop = 4, ios = 2, mainGroup =  True )
parameters[53]  = UniversalParameters ( atomicNumber =  53, eheat =   25.517, iod = 0 , iop = 5, ios = 2, mainGroup =  True )
parameters[54]  = UniversalParameters ( atomicNumber =  54, eheat =    0.000, iod = 0 , iop = 6, ios = 0, mainGroup =  True )
parameters[55]  = UniversalParameters ( atomicNumber =  55, eheat =   18.700, iod = 0 , iop = 0, ios = 1, mainGroup =  True )
parameters[56]  = UniversalParameters ( atomicNumber =  56, eheat =   42.500, iod = 0 , iop = 0, ios = 2, mainGroup =  True )
parameters[57]  = UniversalParameters ( atomicNumber =  57, eheat =  928.900, iod = 1 , iop = 0, ios = 2, mainGroup = False ) # . Lanthanide eheat MNDOd values - for +3 cations (57-71).
parameters[58]  = UniversalParameters ( atomicNumber =  58, eheat =  944.700, iod = 0 , iop = 0, ios = 0, mainGroup = False ) # . Not mainGroup La - Au.
parameters[59]  = UniversalParameters ( atomicNumber =  59, eheat =  952.900, iod = 0 , iop = 0, ios = 0, mainGroup = False )
parameters[60]  = UniversalParameters ( atomicNumber =  60, eheat =  962.800, iod = 0 , iop = 0, ios = 0, mainGroup = False )
parameters[61]  = UniversalParameters ( atomicNumber =  61, eheat =  976.900, iod = 0 , iop = 0, ios = 0, mainGroup = False )
parameters[62]  = UniversalParameters ( atomicNumber =  62, eheat =  974.400, iod = 0 , iop = 0, ios = 0, mainGroup = False )
parameters[63]  = UniversalParameters ( atomicNumber =  63, eheat = 1006.600, iod = 0 , iop = 0, ios = 2, mainGroup = False )
parameters[64]  = UniversalParameters ( atomicNumber =  64, eheat =  991.370, iod = 0 , iop = 0, ios = 2, mainGroup = False )
parameters[65]  = UniversalParameters ( atomicNumber =  65, eheat =  999.000, iod = 0 , iop = 0, ios = 2, mainGroup = False )
parameters[66]  = UniversalParameters ( atomicNumber =  66, eheat = 1001.300, iod = 0 , iop = 0, ios = 2, mainGroup = False )
parameters[67]  = UniversalParameters ( atomicNumber =  67, eheat = 1009.600, iod = 0 , iop = 0, ios = 2, mainGroup = False )
parameters[68]  = UniversalParameters ( atomicNumber =  68, eheat = 1016.150, iod = 0 , iop = 0, ios = 2, mainGroup = False )
parameters[69]  = UniversalParameters ( atomicNumber =  69, eheat = 1022.060, iod = 0 , iop = 0, ios = 2, mainGroup = False )
parameters[70]  = UniversalParameters ( atomicNumber =  70, eheat = 1039.030, iod = 0 , iop = 0, ios = 2, mainGroup = False )
parameters[71]  = UniversalParameters ( atomicNumber =  71, eheat = 1031.200, iod = 1 , iop = 0, ios = 2, mainGroup = False )
parameters[72]  = UniversalParameters ( atomicNumber =  72, eheat =  148.000, iod = 2 , iop = 0, ios = 2, mainGroup = False )
parameters[73]  = UniversalParameters ( atomicNumber =  73, eheat =  186.900, iod = 3 , iop = 0, ios = 2, mainGroup = False )
parameters[74]  = UniversalParameters ( atomicNumber =  74, eheat =  203.100, iod = 5 , iop = 0, ios = 1, mainGroup = False )
parameters[75]  = UniversalParameters ( atomicNumber =  75, eheat =  185.000, iod = 5 , iop = 0, ios = 2, mainGroup = False )
parameters[76]  = UniversalParameters ( atomicNumber =  76, eheat =  188.000, iod = 6 , iop = 0, ios = 2, mainGroup = False )
parameters[77]  = UniversalParameters ( atomicNumber =  77, eheat =  160.000, iod = 7 , iop = 0, ios = 2, mainGroup = False )
parameters[78]  = UniversalParameters ( atomicNumber =  78, eheat =  135.200, iod = 9 , iop = 0, ios = 1, mainGroup = False )
parameters[79]  = UniversalParameters ( atomicNumber =  79, eheat =   88.000, iod = 10, iop = 0, ios = 1, mainGroup = False )
parameters[80]  = UniversalParameters ( atomicNumber =  80, eheat =   14.690, iod = 0 , iop = 0, ios = 2, mainGroup =  True )
parameters[81]  = UniversalParameters ( atomicNumber =  81, eheat =   43.550, iod = 0 , iop = 1, ios = 2, mainGroup =  True )
parameters[82]  = UniversalParameters ( atomicNumber =  82, eheat =   46.620, iod = 0 , iop = 2, ios = 2, mainGroup =  True )
parameters[83]  = UniversalParameters ( atomicNumber =  83, eheat =   50.100, iod = 0 , iop = 3, ios = 2, mainGroup =  True )
parameters[84]  = UniversalParameters ( atomicNumber =  84, eheat =    0.000, iod = 0 , iop = 4, ios = 2, mainGroup =  True )
parameters[85]  = UniversalParameters ( atomicNumber =  85, eheat =    0.000, iod = 0 , iop = 5, ios = 2, mainGroup =  True )
parameters[86]  = UniversalParameters ( atomicNumber =  86, eheat =    0.000, iod = 0 , iop = 6, ios = 0, mainGroup =  True )
parameters[87]  = UniversalParameters ( atomicNumber =  87, eheat =    0.000, iod = 0 , iop = 0, ios = 1, mainGroup =  True )
parameters[88]  = UniversalParameters ( atomicNumber =  88, eheat =    0.000, iod = 0 , iop = 0, ios = 1, mainGroup =  True )
parameters[89]  = UniversalParameters ( atomicNumber =  89, eheat =    0.000, iod = 1 , iop = 0, ios = 2, mainGroup =  True )
parameters[90]  = UniversalParameters ( atomicNumber =  90, eheat = 1674.640, iod = 0 , iop = 0, ios = 4, mainGroup =  True )
parameters[91]  = UniversalParameters ( atomicNumber =  91, eheat =    0.000, iod = 0 , iop = 0, ios = 2, mainGroup =  True )
parameters[92]  = UniversalParameters ( atomicNumber =  92, eheat =    0.000, iod = 0 , iop = 0, ios = 2, mainGroup =  True )
parameters[93]  = UniversalParameters ( atomicNumber =  93, eheat =    0.000, iod = 0 , iop = 0, ios = 2, mainGroup =  True )
parameters[94]  = UniversalParameters ( atomicNumber =  94, eheat =    0.000, iod = 0 , iop = 0, ios = 2, mainGroup =  True )
parameters[95]  = UniversalParameters ( atomicNumber =  95, eheat =    0.000, iod = 0 , iop = 0, ios = 2, mainGroup =  True )
parameters[96]  = UniversalParameters ( atomicNumber =  96, eheat =    0.000, iod = 0 , iop = 0, ios = 2, mainGroup =  True )
parameters[97]  = UniversalParameters ( atomicNumber =  97, eheat =    0.000, iod = 0 , iop = 0, ios = 2, mainGroup =  True )
parameters[98]  = UniversalParameters ( atomicNumber =  98, eheat =    0.000, iod = 0 , iop = 0, ios = 2, mainGroup =  True )
parameters[99]  = UniversalParameters ( atomicNumber =  99, eheat =    0.000, iod = 0 , iop = 0, ios = 0, mainGroup =  True )
parameters[100] = UniversalParameters ( atomicNumber = 100, eheat =    0.000, iod = 0 , iop = 0, ios = 0, mainGroup =  True )
parameters[101] = UniversalParameters ( atomicNumber = 101, eheat =    0.000, iod = 0 , iop = 0, ios = 0, mainGroup =  True )
parameters[102] = UniversalParameters ( atomicNumber = 102, eheat =  207.000, iod = 0 , iop = 0, ios = 0, mainGroup =  True )
parameters[103] = UniversalParameters ( atomicNumber = 103, eheat =    0.000, iod = 0 , iop = 0, ios = 0, mainGroup =  True )
parameters[104] = UniversalParameters ( atomicNumber = 104, eheat =    0.000, iod = 0 , iop = 0, ios = 0, mainGroup =  True )
parameters[105] = UniversalParameters ( atomicNumber = 105, eheat =    0.000, iod = 0 , iop = 0, ios = 0, mainGroup =  True )
parameters[106] = UniversalParameters ( atomicNumber = 106, eheat =    0.000, iod = 0 , iop = 0, ios = 0, mainGroup =  True )
parameters[107] = UniversalParameters ( atomicNumber = 107, eheat =    0.000, iod = 0 , iop = 0, ios = 0, mainGroup =  True )

#===================================================================================================================================
# . Functions for making derived parameters.
#===================================================================================================================================
def Factorial ( n ):
    """Factorial."""
    result = 1
    if n > 1:
        for i in range ( 2, n + 1 ): result *= i 
    return result 

def MakeExtraParameters ( n, p, minhpp, toExclude, isPMO ):
    """Make extra parameters for non-hydrogen elements."""
    data = {}
    if toExclude is None: toExclude = set ( )
    if ( n > 1 ) or isPMO:

        # . Preliminary checks.
# . Three lines with arbitrary parameters.
        if ( p["zetap"] > 0.0001 ) or ( p["zetas"] > 0.0001 ):
            if p["hsp"  ] < 1.0e-7: p["hsp"  ] = 1.0e-7
            if p["zetap"] < 0.3   : p["zetap"] = 0.3

            # . Set some intermediate variables.
            gpp = p["gpp"]
            gp2 = p["gp2"]
            gdd = p.get ( "gdd", 0.0 )
            gsd = p.get ( "gsd", 0.0 )
            gsp = p["gsp"]
            gss = p["gss"]
            hsp = p["hsp"]
            zp  = p["zetap"]
            zs  = p["zetas"]
            upp = p["upp"]
            uss = p["uss"]
            udd = p.get ( "udd", 0.0 )

            qns = p["qns"]
            qnp = p["qnp"]

            # . Start calculation.
            hpp = 0.5 * ( gpp - gp2 )
# . Arbitrary minhpp value - particularly bad for pm3!
            if math.fabs ( hpp ) < minhpp: print ( "\nLimiting hpp for element {:d}: {:15.5f} {:15.5f}.".format ( n, hpp, minhpp ) )
            hpp = max ( minhpp, hpp )

            # . Eisol - old.
            # . Mismatch for PM6 rare gases (OK) and for Tc/W/Pt/Au (Tc/Pt/Au - GSSC/gssc, W - USSC/ios, UDDC/iod, GSSC/gssc).
#            eisol = uss * USSC[n] + upp * UPPC[n] + udd * UDDC[n] + gss * GSSC[n] + gpp * GPPC[n] + gsp * GSPC[n] +  gp2 * GP2C[n] + hsp * HSPC[n] + gsd * GSDC[n] +  gdd * GDDC[n]

            # . Eisol - new.
            ios = parameters[n].ios
            iop = parameters[n].iop
            iod = parameters[n].iod
            # . GSSC is the number of two-electron terms of type <SS|SS>.
            gssc = max ( ios - 1, 0 ) 
            k    = iop 
            # . GSPC is the number of two-electron terms of type <SS|PP>.
            gspc = ios * iop 
            l    = min ( k, 6 - k ) 
            # . GP2C is the number of two-electron terms of type <PP|PP> plus 0.5 of the number of HPP integrals.
            gp2c = float ( (k*(k - 1))//2 ) + 0.5 * float ( (l*(l - 1))//2 ) 
            # . GPPC is minus 0.5 times the number of HPP integrals.
            gppc = float ( -0.5 * (l*(l - 1))//2 ) 
            # . HSPC is the number of two-electron terms of type <SP|SP>. (S and P must have the same spin.  In all cases, if P is non-zero, there are two S electrons).
            hspc  = -k 
            eisol = uss * ios + upp * iop + udd * iod + gss * gssc + gpp * gppc + gsp * gspc + gp2 * gp2c + hsp * hspc

            # . Print eisols.
            if n == 1:
                print ( "Eisols for {:d}: full - {:20.5f}; s contribution - {:20.5f}.".format ( n, eisol, ( uss * float ( ios ) ) ) )

            # . Charge separations.
            # . Dipole.
# . Old formula where qns = qnp.
#            dd    = ( 2.0 * qn + 1.0 ) * ( 4.0 * zs * zp )**( qn + 0.5 ) / ( zs + zp )**( 2.0 * qn + 2.0 ) / math.sqrt ( 3.0 )
# . New formula where qns != qnp.
            dd    = ( 2.0 * zs )**( qns + 0.5 ) * ( 2.0 * zp )**( qnp + 0.5 ) * float ( Factorial ( qns + qnp + 1 ) ) / ( math.sqrt ( 3.0 ) * ( zs + zp )**( qns + qnp + 2.0 ) * math.sqrt ( float ( Factorial ( 2 * qns ) * Factorial ( 2 * qnp ) ) ) )

            # . Quadrupole.
            qq    = math.sqrt ( ( 4.0 * qnp * qnp + 6.0 * qnp + 2.0 ) / 20.0 ) / zp

            # . Calculate additive terms in atomic units - iterative procedure.
            # . Parameters.
            jmax = 20
            P    = 2.0
            P2   = P*P
            P4   = P2*P2

            # . Iterate.
            gdd1 = ( P2 * hsp / ( _MOPAC_HARTREES_TO_ELECTRONVOLTS *  4.0 * dd**2 ) )**(1.0/3.0)
            gqq  = ( P4 * hpp / ( _MOPAC_HARTREES_TO_ELECTRONVOLTS * 48.0 * qq**4 ) )**0.2
            d1   = gdd1
            d2   = gdd1 + 0.04
            q1   = gqq
            q2   = gqq  + 0.04
            for j in range ( jmax ):
                df   = d2 - d1
                if math.fabs ( df ) <= 1.0e-10:
                    print ( "Breaking {:d} d loop after {:d} iterations with value {:.10f}.".format ( n, j, d2 ) )
                    break
                hsp1 = 2.0 * d1 - 2.0/math.sqrt(4.0 * dd**2 + 1.0/d1**2)
                hsp2 = 2.0 * d2 - 2.0/math.sqrt(4.0 * dd**2 + 1.0/d2**2)
                hsp1 /= P2
                hsp2 /= P2
                d3   = d1  +  df * ( hsp / _MOPAC_HARTREES_TO_ELECTRONVOLTS - hsp1 ) / ( hsp2 - hsp1 )
                d1   = d2
                d2   = d3
            # . HPP used here.
            for j in range ( jmax ):
                qf   = q2 - q1
                if math.fabs ( qf ) <= 1.0e-10:
                    print ( "Breaking {:d} q loop after {:d} iterations with value {:.10f}.".format ( n, j, q2 ) )
                    break
                hpp1 = 4.0 * q1 - 8.0/math.sqrt(4.0 * qq**2 + 1.0/q1**2) +  4.0/math.sqrt(8.0 * qq**2 + 1.0/q1**2)
                hpp2 = 4.0 * q2 - 8.0/math.sqrt(4.0 * qq**2 + 1.0/q2**2) +  4.0/math.sqrt(8.0 * qq**2 + 1.0/q2**2)
                hpp1 /= P4
                hpp2 /= P4
                q3   = q1 + qf * ( hpp / _MOPAC_HARTREES_TO_ELECTRONVOLTS - hpp1 ) / ( hpp2 - hpp1 )
                q1   = q2
                q2   = q3

    # . Hydrogen.
    elif n == 1:
        eisol = p["uss"]
        gss   = p["gss"]
        d2    = gss / _MOPAC_HARTREES_TO_ELECTRONVOLTS
        q2    = d2
        dd    = 0.0
        qq    = 0.0

        # . Zero values.
        for key in ( "betap", "dd", "gpp", "gp2", "gsp", "hsp", "hsp", "upp", "zetap" ): data[key] = 0.0

    # . Calculate norbitals and zcore.
    ns = p.get ( "qns", -1 )
    np = p.get ( "qnp", -1 )
    nd = p.get ( "qnd", -1 )
    norbitals = 0
    if ns >= 0: norbitals += 1
    if np >= 1: norbitals += 3
    if nd >= 2: norbitals += 5
    zcore = 0.0
    if ns >= 0: zcore += float ( parameters[n].ios )
    if np >= 0: zcore += float ( parameters[n].iop )
    if nd >= 0: zcore += float ( parameters[n].iod )

    # . Store the results.
    data["am"]        = gss / _MOPAC_HARTREES_TO_ELECTRONVOLTS
    extraaddd         = _EXTRAADDD.get ( n, None )
    if extraaddd is not None: ( d2, dd ) = extraaddd
    extraaqqq         = _EXTRAAQQQ.get ( n, None )
    if extraaqqq is not None: ( q2, qq ) = extraaqqq
    data["ad"]        = d2
    data["aq"]        = q2 # . Requires hpp.
    data["dd"]        = dd
    data["eheat"]     = parameters[n].eheat
    if "eisol" not in toExclude:
        extraeisol = _EXTRAEISOL.get ( n, None )
        if extraeisol is None: data["eisol"] = eisol
        else:                  data["eisol"] = extraeisol
    data["orbitals"]  = norbitals
    data["qq"]        = qq
    data["zcore"]     = zcore
    return data

#===================================================================================================================================
# . One-center two-electron integral data for PM6.
#===================================================================================================================================
class OCTEI:
    def __init__ ( self ):
        self.iii   = 0
        self.iiid  = 0
        self.ir016 = 0
        self.ir066 = 0
        self.ir244 = 0
        self.ir266 = 0
        self.ir466 = 0

def GetPM61CTEIData ( doPrint = False ):
    """Get the data."""

    octeidata = {}
    for i in range ( 107 ): octeidata[i+1] = OCTEI ( )

    # . iii and iiid.
    # data iii/ 2*1, 8*2, 8*3, 18*4, 18*5, 32*6, 21*0/  
    # data iiid/ 30*3, 18*4, 32*5, 6*6, 21*0/  
    iii  = 2 * [ 1 ] + 8 * [ 2 ] + 8 * [ 3 ] + 18 * [ 4 ] + 18 * [ 5 ] + 32 * [ 6 ] + 21 * [ 0 ]
    iiid = 30 * [ 3 ] + 18 * [ 4 ] + 32 * [ 5 ] + 6 * [ 6 ] + 21 * [ 0 ]
    for ( i, ( a, b ) ) in enumerate ( zip ( iii, iiid ) ):
        data = octeidata[i+1]
        data.iii  = a
        data.iiid = b

    #                                Sc  Ti   V  Cr  Mn  Fe  Co  Ni  Cu
    # Atomic orbital population   4s  2   2   2   1   2   2   2   2   1
    # of gaseous atom             3d  1   2   3   5   5   6   7   8  10
    #
    # State term:                    2D  3F  4F  7S  6S  5D  4F  3F  2S

    #      data (ir016(i),i=21,29)/ 2, 4, 6, 5, 10, 12, 14, 16, 10/  
    #      data (ir066(i),i=21,29)/ 0, 1, 3, 10, 10, 15, 21, 28, 45/  
    #      data (ir244(i),i=21,29)/ 1, 2, 3, 5, 5, 6, 7, 8, 5/  
    #      data (ir266(i),i=21,29)/ 0, 8, 15, 35, 35, 35, 43, 50, 70/  
    #      data (ir466(i),i=21,29)/ 0, 1, 8, 35, 35, 35, 36, 43, 70/  

    for ( i, value ) in enumerate ( [ 2, 4,  6,  5, 10, 12, 14, 16, 10 ] ): octeidata[i+21].ir016 = value
    for ( i, value ) in enumerate ( [ 0, 1,  3, 10, 10, 15, 21, 28, 45 ] ): octeidata[i+21].ir066 = value
    for ( i, value ) in enumerate ( [ 1, 2,  3,  5,  5,  6,  7,  8,  5 ] ): octeidata[i+21].ir244 = value
    for ( i, value ) in enumerate ( [ 0, 8, 15, 35, 35, 35, 43, 50, 70 ] ): octeidata[i+21].ir266 = value
    for ( i, value ) in enumerate ( [ 0, 1,  8, 35, 35, 35, 36, 43, 70 ] ): octeidata[i+21].ir466 = value

    #                                 Y  Zr  Nb  Mo  Tc  Ru  Rh  Pd  Ag
    # Atomic orbital population   5s  2   2   1   1   2   1   1   0   1
    # of gaseous atom             4d  1   2   4   5   5   7   8  10  10
    #
    # State term:                    2D  3F  6D  7S  6D  5F  4F  1S  2S

    #      data (ir016(i),i=39,47)/ 2, 4, 4, 5, 10, 7, 8, 0, 10/  
    #      data (ir066(i),i=39,47)/ 0, 1, 6, 10, 10, 21, 28, 45, 45/  
    #      data (ir244(i),i=39,47)/ 1, 2, 4, 5, 5, 5, 5, 0, 5/  
    #      data (ir266(i),i=39,47)/ 0, 8, 21, 35, 35, 43, 50, 70, 70/  
    #      data (ir466(i),i=39,47)/ 0, 1, 21, 35, 35, 36, 43, 70, 70/  

    for ( i, value ) in enumerate ( [ 2, 4,  4,  5, 10,  7,  8,  0, 10 ] ): octeidata[i+39].ir016 = value
    for ( i, value ) in enumerate ( [ 0, 1,  6, 10, 10, 21, 28, 45, 45 ] ): octeidata[i+39].ir066 = value
    for ( i, value ) in enumerate ( [ 1, 2,  4,  5,  5,  5,  5,  0,  5 ] ): octeidata[i+39].ir244 = value
    for ( i, value ) in enumerate ( [ 0, 8, 21, 35, 35, 43, 50, 70, 70 ] ): octeidata[i+39].ir266 = value
    for ( i, value ) in enumerate ( [ 0, 1, 21, 35, 35, 36, 43, 70, 70 ] ): octeidata[i+39].ir466 = value

    #                                    Hf  Ta   W  Re  Os  Ir  Pt  Au Hg
    # Atomic orbital population   6s      2   2   2   2   2   2   1   1  2
    # of gaseous atom             5d      2   3   4   5   6   7   9  10  0
    #
    # State term:                        3F  4F  5D  6S  5D  4F  3D  2S 1S

    #      data (ir016(i),i=72,80)/ 4, 6, 8, 10, 12, 14, 9, 10, 0/  
    #      data (ir066(i),i=72,80)/ 1, 3, 6, 10, 15, 21, 36, 45, 0/  
    #      data (ir244(i),i=72,80)/ 2, 3, 4, 5, 6, 7, 5, 5, 0/  
    #      data (ir266(i),i=72,80)/ 8, 15, 21, 35, 35, 43, 56, 70, 0/  
    #      data (ir466(i),i=72,80)/ 1, 8, 21, 35, 35, 36, 56, 70, 0/  

    for ( i, value ) in enumerate ( [ 4,  6,  8, 10, 12, 14,  9, 10, 0 ] ): octeidata[i+72].ir016 = value
    for ( i, value ) in enumerate ( [ 1,  3,  6, 10, 15, 21, 36, 45, 0 ] ): octeidata[i+72].ir066 = value
    for ( i, value ) in enumerate ( [ 2,  3,  4,  5,  6,  7,  5,  5, 0 ] ): octeidata[i+72].ir244 = value
    for ( i, value ) in enumerate ( [ 8, 15, 21, 35, 35, 43, 56, 70, 0 ] ): octeidata[i+72].ir266 = value
    for ( i, value ) in enumerate ( [ 1,  8, 21, 35, 35, 36, 56, 70, 0 ] ): octeidata[i+72].ir466 = value

    # . Print.
    if doPrint:
        print ( "\nData:\n" )
        keys = list ( octeidata.keys ( ) )
        keys.sort ( )
        for key in keys:
            data = octeidata[key]
            print ( "{:5d} {:<5s} {:5d} {:5d} {:5d} {:5d} {:5d} {:5d} {:5d}".format ( key, PeriodicTable.Symbol ( key ), data.iii, data.iiid, data.ir016, data.ir066, data.ir244, data.ir266, data.ir466 ) )

    return octeidata

#===================================================================================================================================
# . Convert YAML to text.
#===================================================================================================================================
def MNDOParameters_YAMLToText ( inPath, outPath ):
    """Convert YAML to text."""

    # . Read the file.
    parameters     = YAMLMappingFile_ToObject ( inPath, MNDOParameters )
    ( head, tail ) = os.path.split ( inPath )
    atomicNumber   = PeriodicTable.AtomicNumber ( tail[0:-4] )
    symbol         = PeriodicTable.Symbol ( atomicNumber )

    # . Write out the data.
    txtfile = open ( outPath, "w" )
    txtfile.write ( "# . Text parameter file for the element " + symbol + ".\n" )

    # . Get the data.
    state = parameters.__getstate__ ( )

    # . Scalar parameters.
    parameters = state.pop ( "Scalar Parameters", [] )
    scalars    = {}
    for ( label, value, units ) in parameters:
        scalars[label] = value

    # . Treat norbitals.
    norbitals = scalars.pop ( "orbitals", 0 )
    if ( norbitals > 0 ) and ( "qns" not in scalars ):
        eN = PeriodicTable[atomicNumber]
        scalars["qns"] = eN.sQuantumNumber
        if norbitals > 1: scalars["qnp"] = eN.pQuantumNumber
        if norbitals > 4: scalars["qnd"] = eN.dQuantumNumber

    # . Start output.
    txtfile.write ( "\n" )

    # . Integer scalars.
    for key in _INTEGERATTRIBUTES:
        value = scalars.pop ( key, None )
        if value is not None:
            txtfile.write ( "{:<5s}{:<15s}{:25d}\n".format ( symbol, key, value ) )

    # . Get ordered keys for remaining scalars.
    keys = list ( scalars.keys ( ) )
    keys.sort ( )

    # . Float scalars.
    for key in keys:
        if isinstance ( scalars[key], float ):
            value = scalars.pop ( key )
            txtfile.write ( "{:<5s}{:<15s}{:25.15f}  {:<s}\n".format ( symbol, key, value, _UNITS[key] ) )

    # . Array state.
    _OutputArray ( state.pop ( "AM1/PM3 Gaussian Parameters", [] ), txtfile, symbol, _AM1PM3GNAMES )
    _OutputArray ( state.pop ( "PDDG Gaussian Parameters"   , [] ), txtfile, symbol, _PDDGNAMES    )

    # . Diatomic interactions.
    diatomics = state.get ( "Diatomic Parameters", [] )
    if len ( diatomics ) > 0:
        for ( n, x0, a0 ) in diatomics:
            txtfile.write ( "{:<5s}{:<15s}{:<5s}{:20.15f}  {:<s}\n".format ( symbol, "diatomicA", PeriodicTable.Symbol ( n ), a0, _UNITS["diatomicA"] ) )
            txtfile.write ( "{:<5s}{:<15s}{:<5s}{:20.15f}  {:<s}\n".format ( symbol, "diatomicX", PeriodicTable.Symbol ( n ), x0, _UNITS["diatomicX"] ) )

# . Private functions.
def _OutputArray ( array, txtfile, symbol, names ):
    """Output an array."""
    for ( i, data ) in enumerate ( array ):
        for ( datum, name ) in zip ( data, names ):
            txtfile.write ( "{:<5s}{:<15s}{:25.15f}  {:<s}\n".format ( symbol, name + repr ( i ), datum, _UNITS[name] ) )

#===================================================================================================================================
# . Convert text to YAML.
#===================================================================================================================================
def MNDOParameters_TextToYAML ( inPath, outPath, isPMO = False, minhpp = 0.1, toExclude = None, QWRITE = True ):
    """Convert text to YAML."""

    # . Read the file.
    data = ReadTextFile ( inPath )

    # . Treat the data.
    if ( data is not None ) and ( len ( data ) > 0 ):

        # . Get OCTEI data.
        octeidata = GetPM61CTEIData ( )

        # . Loop over atomic numbers found on the file.
        atomicNumbers = list ( data.keys ( ) )
        atomicNumbers.sort ( )
        for atomicNumber in atomicNumbers:
            existing   = data[atomicNumber]
            extras     = MakeExtraParameters ( atomicNumber, existing, minhpp, toExclude, isPMO )
            keys       = list ( extras.keys ( ) )
            keys.sort ( )
            deviations = {}
            for key in keys:
                if key in existing:
                    deviation = extras[key] - existing[key]
                    if not isinstance ( deviation, float ): deviation = float ( deviation )
                    deviation = math.fabs ( deviation )
                    if deviation > _MAXIMUMDEVIATION: deviations[key] = ( extras[key], existing[key], deviation )
            if len ( deviations ) > 0:
                keys = list ( deviations.keys ( ) )
                keys.sort ( )
                print ( "\nExtra parameter deviations for element {:d} from {:s}:".format ( atomicNumber, inPath ) )
                for key in keys:
                    ( extra, old, deviation ) = deviations[key]
                    if isinstance ( extra, float ): print ( "{:<20s} {:25.15f} {:25.15f} {:25.15f}".format ( key, extra, old, deviation ) )
                    else:                           print ( "{:<20s} {:25d} {:25d} {:25.15f}"      .format ( key, extra, old, deviation ) )
            # . Update the parameters.
            existing.update ( extras )

            # . PM6 data.
            data = octeidata[atomicNumber]
            for attribute in ( "iii", "iiid", "ir016", "ir066", "ir244", "ir266", "ir466" ): existing[attribute] = getattr ( data, attribute )

            # . Construct scalar attributes.
            state      = {}
            parameters = []
            for label in _INTEGERATTRIBUTES:
                value = existing.pop ( label, None )
                if value is not None:
                    parameters.append ( [ label, value, "none" ] )
            for label in _REALATTRIBUTES:
                value = existing.pop ( label, None )
                if value is not None:
                    parameters.append ( [ label, value, _UNITS[label] ] )
            state["Scalar Parameters"] = parameters


            # . AM1/PM3 parameters.
            nam1pm3g = 0
            for i in range ( 4 ):
                if ( "fn1" + repr ( i ) ) in existing: nam1pm3g += 1
            if nam1pm3g > 0:
                parameters = []
                for i in range ( nam1pm3g ):
                    parameters.append ( [ existing["fn1" + repr ( i )] , 
                                          existing["fn2" + repr ( i )] ,
                                          existing["fn3" + repr ( i )] ] )
                state["AM1/PM3 Gaussian Parameters"] = parameters

            # . PDDG parameters.
            npddg = 0
            for i in range ( 2 ):
                if ( "pddgc" + repr ( i ) ) in existing: npddg += 1
            if npddg > 0:
                parameters = []
                for i in range ( npddg ):
                    parameters.append ( [ existing["pddgc" + repr ( i )] ,
                                          existing["pddge" + repr ( i )] ] )
                state["PDDG Gaussian Parameters"] = parameters

            # . Treat diatomic data - note ndiatomic is the highest atomic number for which a term is defined + 1.
            diatomicA = existing.get ( "diatomicA", {} )
            diatomicX = existing.get ( "diatomicX", {} )
            nA = len ( diatomicA )
            nX = len ( diatomicX )
            if ( nA > 0 ) and ( nA == nX ):
                partners = list ( set ( diatomicA.keys ( ) ) | set ( diatomicX.keys ( ) ) )
                partners.sort ( )
                parameters = []
                for partner in partners:
                    parameters.append ( [ partner, diatomicX[partner], diatomicA[partner] ] )
                state["Diatomic Parameters"] = parameters

            # . Create the MNDO parameters object.
            mndoparameters = MNDOParameters.Raw ( )
            mndoparameters.__setstate__ ( state )

            # . Create the yaml file.
            if QWRITE: YAMLMappingFile_FromObject ( os.path.join ( outPath, PeriodicTable.Symbol( atomicNumber ) + _YAMLEXTENSION ), _YAMLTag, mndoparameters )

#===================================================================================================================================
# . Read text file.
#===================================================================================================================================
def ReadTextFile ( inPath ):
    """Read a text file."""
    data = None
    if os.path.exists ( inPath ):
        data  = {}
        units = {}
        keys  = set ( )
        dfile = open ( inPath, "r" )
        for line in dfile:
            words = line.split ( )
            if ( len ( words ) > 0 ) and ( not words[0].startswith ( "#" ) ):
                if ( len ( words ) >= 3 ):
                    atomicNumber = PeriodicTable.AtomicNumber ( words[0] )
                    if atomicNumber not in data: data[atomicNumber] = {}
                    parameter    = words[1]
                    if len ( words ) == 3:
                        value = int   ( words[2] )
                        unit  = None
                    elif len ( words ) == 4:
                        value = float ( words[2] )
                        unit  = " ".join ( words[3:] )
                    elif len ( words ) == 5:
                        partner        = PeriodicTable.AtomicNumber ( words[2] )
                        value          = data[atomicNumber].get ( parameter, {} )
                        value[partner] = float ( words[3] )
                        unit           = " ".join ( words[4:] )
                    data[atomicNumber][parameter] = value
                    # . Check units.
                    if parameter in units:
                        oldunit = units[parameter]
                        if ( ( oldunit is None ) and ( unit is None ) ) or ( oldunit == unit ): pass
                        else: print ( "Unit mismatch for parameter: " + parameter + " - " + oldunit + " and " + unit + "." )
                    else:
                        units[parameter] = unit
                    keys.add ( parameter )
                else:
                    raise ValueError ( "Invalid text file line: " + line.strip ( ) + "." )
        # . Finish up.
        known   = _UNDERSTOODNAMES
        unknown = keys.difference ( known )
        if False:
            print ( "\nUnique parameters:" )
            keys = list ( keys )
            keys.sort ( )
            for key in keys: print ( key )
        if len ( unknown ) > 0:
            print ( "\nUnknown parameters:" )
            unknown = list ( unknown )
            unknown.sort ( )
            for key in unknown: print ( key )
    return data

#===================================================================================================================================
# . Convert text to YAML.
#===================================================================================================================================
def MNDOParametersTextToYAML ( ):
    """Convert MNDO parameters in text format to YAML."""

    # . Directories.
    textPathIn = os.path.join ( _TESTPATH, _TEXTPATH   )
    yamlPathOut = os.path.join ( _TESTPATH, _YAMLPATHOUT )
    if not os.path.exists ( yamlPathOut ): os.mkdir ( yamlPathOut )

    # . Loop over the files for each method.
    for hamiltonian in _HAMILTONIANS:

        # . Output directory.
        outPath = os.path.join ( yamlPathOut, hamiltonian )
        if not os.path.exists  ( outPath ): os.mkdir ( outPath )

        # . Get the text files.
        textPaths = glob.glob ( os.path.join ( textPathIn, hamiltonian, "*" + _TEXTEXTENSION ) )
        for textPath in textPaths:
            MNDOParameters_TextToYAML ( textPath ,
                                        outPath  ,
                                        isPMO     = hamiltonian.startswith ( "pmo" )                ,
                                        minhpp    = _HAMILTONIANMINHPP.get    ( hamiltonian, 0.1  ) ,
                                        toExclude = _HAMILTONIANTOEXCLUDE.get ( hamiltonian, None ) )

#===================================================================================================================================
# . Convert YAML to text.
#===================================================================================================================================
def MNDOParametersYAMLToText ( ):
    """Convert MNDO parameters in YAML format to text."""

    # . Directories.
    textPathOut = os.path.join ( _TESTPATH, _TEXTPATH )
    if not os.path.exists ( _TESTPATH   ): os.mkdir ( _TESTPATH   )
    if not os.path.exists ( textPathOut ): os.mkdir ( textPathOut )

    # . Loop over the files for each method.
    for hamiltonian in _HAMILTONIANS:

        # . Output directory.
        outPath = os.path.join ( textPathOut, hamiltonian )
        if not os.path.exists ( outPath ): os.mkdir ( outPath )

        # . Get the YAML files.
        yamlPaths = glob.glob ( os.path.join ( _YAMLPATHIN, hamiltonian, "*" + _YAMLEXTENSION ) )
        for yamlPath in yamlPaths:
            ( head, tail ) = os.path.split ( yamlPath )
            MNDOParameters_YAMLToText ( yamlPath, os.path.join ( outPath, tail[0:-5] + _TEXTEXTENSION ) )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :

    # . Testing.
    MNDOParametersYAMLToText ( )
    MNDOParametersTextToYAML ( )
