"""Constants."""

import  os

YAMLPATHIN       = os.path.join ( os.getenv ( "PDYNAMO3_PARAMETERS" ), "pcetk" )

# . Maximum number of states for analytic treatment (set arbitraily to 2^26)
ANALYTIC_SITES   = 26
ANALYTIC_STATES  = 2**ANALYTIC_SITES

PREV_RESIDUE     = ("C", "O")
NEXT_RESIDUE     = ("N", "H",  "CA", "HA")
NEXT_RESIDUE_PRO = ("N", "CA", "HA", "CD",  "HD1", "HD2")
NEXT_RESIDUE_GLY = ("N", "H",  "CA", "HA1", "HA2", "HN")
TERM_REMOVE      = ("C", "O", "N", "H")

PROTEIN_RESIDUES = (
    "ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU",
    "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR",
    "HSP", "HSE", "HSD", "HIE", "HID"
)

REMOVE_RESIDUES  = (
    "WAT", "HOH", "TIP", "TIP3", "TP3M", "SOD",
)
