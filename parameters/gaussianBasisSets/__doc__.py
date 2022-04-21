"""Documentation information."""

# . Aliases.
Aliases = { "3-21g"             : "3-21G"                      ,
            "6-31g_st"          : "6-31G*"                     ,
            "def2-sv(p)"        : "Turbomole SV(P)"            ,
            "def2-svp"          : "Turbomole SVP"              ,
            "def2-sv(p)-rifit"  : "Turbomole SV(P) RI-fitting" ,
            "def2-svp-rifit"    : "Turbomole SVP RI-fitting"   ,
            "def2-tzvp"         : "Turbomole TZVP"             ,
            "def2-tzvp-rifit"   : "Turbomole TZVP RI-fitting"  ,
            "dgauss-a1-dftjfit" : "Dgauss    A1 DFT J-fitting" ,
            "dgauss-dzvp"       : "Dgauss    DZVP"             ,
            "sto-3g"            : "STO-3G"                     }

# . Descriptions.
Descriptions = { }

# . Notes.
Notes = """Bases are from the basis set exchange - https://www.basissetexchange.org.
           The same names are used. Full references for each of them can be found there.

           Only 3-21g and 6-31g_st employ Cartesian functions, the others use (real)
           spherical harmonics."""

# . Types.
Types = { "3-21g"             : "GaussianBasis" ,
          "6-31g_st"          : "GaussianBasis" ,
          "def2-sv(p)"        : "GaussianBasis" ,
          "def2-svp"          : "GaussianBasis" ,
          "def2-sv(p)-rifit"  : "GaussianBasis" ,
          "def2-svp-rifit"    : "GaussianBasis" ,
          "def2-tzvp"         : "GaussianBasis" ,
          "def2-tzvp-rifit"   : "GaussianBasis" ,
          "dgauss-a1-dftjfit" : "GaussianBasis" ,
          "dgauss-dzvp"       : "GaussianBasis" ,
          "sto-3g"            : "GaussianBasis" }
