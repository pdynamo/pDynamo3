![pDynamo *logo*](logo.png)

# pDynamo3 (Fork for ORCA 6.0.0 Compatibility)

This repository is a **fork** of the original [pDynamo3](https://github.com/dynamo.modeling/pdynamo3) project, modified to work with **ORCA version 6.0.0**.  
The changes specifically adjust the coupling to ORCA to ensure compatibility with ORCA 6.0.0 syntax and behavior.

---

pDynamo is an open source program library designed for the simulation of molecular systems using quantum chemical (QC), molecular mechanical (MM), and hybrid QC/MM potential energy functions.  
pDynamo is written in Python with the computationally intensive parts of the code implemented in C and Cython. The current version of pDynamo, pDynamo3, uses Python 3.

Principal author: Martin Field  
Development team: <dynamo.modeling@gmail.com>  

Released under the [GNU General Public License](gpl-3.0.txt)

---

## Features of pDynamo3
- Density functional theory and Hartree-Fock QC methods employing Gaussian basis sets
- Semi-empirical QC methods of the MNDO type, including AM1, MNDO, PDDG, PM3, RM1 and PM6
- Support for some standard MM force fields, including AMBER, CHARMM and OPLS-AA
- Hybrid QC/MM methods using any combination of the QC and MM potentials implemented in the library
- Coupling to third-party programs
- Energy calculations
- Geometry optimizations
- Transition state searches
- Reaction path calculations
- Normal mode analyses
- Property calculations, such as charges and dipoles
- Molecular dynamics simulations
- Monte Carlo simulations
- Various geometrical restraints
- The ability to handle various common molecular file formats
- Miscellaneous analysis tools

---

## Citation
If you use pDynamo or this modified version in your work, please cite:

- M. J. Field, ["The *pDynamo* Library for Molecular Simulations using Hybrid Quantum Mechanical and Molecular Mechanical Potentials"](https://pubs.acs.org/doi/10.1021/ct800092p), *J. Chem. Theo. Comp.* **2008**, *4*, 1151-1161.
- M. J. Field, ["A Practical Introduction to the Simulation of Molecular Systems"](https://www.cambridge.org/core/books/practical-introduction-to-the-simulation-of-molecular-systems/E91B9A8E90237C3D63F6A589105FF38B), *Cambridge University Press.* **2007**.

---

## Installation instructions
Required software and libraries:

- Python 3, version 3.5 or higher (including header files; python3-dev package in Debian)
- [Cython](https://cython.org/) (cython3 package in Debian)
- PyYAML (python3-yaml package in Debian)
- C compiler (gcc or other)

Download this modified version of pDynamo3 from GitHub either as a zip file or by cloning the repository.  
Unpack the zip file or place the repository where pDynamo3 is to be installed.

Several modules within pDynamo3's packages and subpackages are written in C and Cython and have to be compiled before use.  
To do this, go to the installation directory:

```bash
cd installation
python3 Install.py -f
```

Installation should take a few minutes.

Some environment variables need to be set. Example bash and csh scripts are generated at the end of installation in the `shellScripts` subdirectory.  
Minimal environment variables for operation:

```csh
# Example for cshell users:
setenv PDYNAMO3_HOME          $HOME/pDynamo3
setenv PDYNAMO3_SCRATCH       $PDYNAMO3_HOME/scratch
setenv PYTHONPATH             .:$PDYNAMO3_HOME
setenv PDYNAMO3_PARAMETERS    $PDYNAMO3_HOME/parameters
setenv PDYNAMO3_PYTHONCOMMAND python3
setenv PDYNAMO3_STYLE         $PDYNAMO3_PARAMETERS/ccsStyleSheets/defaultStyle.css
```

You can test the installation by running:

```bash
python3 RunExamples.py book
```

A full list of available tests can be listed with:

```bash
python3 RunExamples.py -l
```

Output logs and generated files will be placed in the `examples` subdirectory of `PDYNAMO3_SCRATCH`.

---

## Additional packages and programs

pDynamo3 includes third-party extension packages in the `addOns` directory, including:

- [pcetk](https://github.com/mfx9/pcetk): proton binding energetics in proteins.
- [pyCPR](http://www.bisb.uni-bayreuth.de/index.php?page=data/PyCPR/PyCPR): transition state search tool.

Optional modules in pDynamo3 require external programs:

- [DFTB+](https://dftbplus.org/)
- [ORCA](https://orcaforum.kofo.mpg.de/app.php/portal) (this fork supports ORCA 6.0.0)
- [PySCF](https://github.com/pyscf/pyscf)
- [extended-MEAD](http://www.bisb.uni-bayreuth.de/index.php?page=downloads)

You must set additional environment variables for these programs:

```csh
# Example for cshell users:
setenv PDYNAMO3_DFTBCOMMAND /path/to/dftb+/dftb+
setenv PDYNAMO3_ORCACOMMAND /path/to/orca
setenv PDYNAMO3_MEADPATH /path/to/mead/bin
setenv PDYNAMO3_PYSCFPATH /path/to/pyscf
setenv PYTHONPATH .:$PDYNAMO3_HOME/:$PDYNAMO3_PYSCFPATH/
```

---

## Further information
For more details and examples, visit the original [pDynamo website](https://sites.google.com/site/pdynamomodeling/home).

---
