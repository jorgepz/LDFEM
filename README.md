# LDFEM (currently Unmaintained)

LDFEM is a simple code for Large Deformations hyperelastic solids analysis using the Finite Element Method and Fortran90.


![example deformation](https://github.com/jorgepz/LDFEM/blob/main/def.png)

## About

The code recieves as input a [gmsh](https://gmsh.info/#Download) legacy file and a specific .inp file and, after solving using the Newton-Raphson method, generates results and [VTK](https://vtk.org/) files.

The code was developed between 2014 and 2015 by Jorge Pérez Zerpa ([@jorgepz](https://github.com/jorgepz)), with contributions (on stress tensors output and compressible case tangent matrix) made by Pablo Castrillo ([@pablocgp](https://github.com/pablocgp)) in December 2014, at the Instituto de Estructuras y Transporte de la [Facultad de Ingeniería](https://www.fing.edu.uy/) de la Universidad de la República, Uruguay.

## Using LDFEM

### compilation

To compile the code you need to have installed gcc's fortran compiler and also the BLAS libraries in a linux environment. You can execute the script `build.sh` to generate the `ldfem.lnx` binary.

### input
LDFEM uses a specific .inp format file and

### output

Example terminal output:

```
====================================
-------   Welcome to LDFEM   -------
 Los factores de carga a resolver son:   0.20000000000000001       0.26666666666666666       0.26666666666666666       0.26666666666666666     
 barra.msh           
Totales:
  nodos:    14
  tetra:    24
  trian:    24
-------------------------------
Control node number          1/         1
Closest node to coordinates entered is:          2 distance:    0.000E+00
Coordinates:    0.100E+01   0.000E+00   0.000E+00
 boundary condtitions: nneumdofs:          22  Ndiridofs:          20
 Nnonzero diridofs:           5
    0.00    0.00    0.00    0.00    0.00    0.00    0.00
    0.00    1.00    0.00    0.00    0.00    0.00    0.00
    0.00    0.00    1.00    0.00    0.00    0.00    0.00
    0.00    0.00    0.00    1.00    0.00    0.00    0.00
    0.00    1.00    0.00    0.00    0.20    0.00    0.00

-------  New time step: -------
 indtime:    1 factor:   0.20E+00 (delta:   0.27E+00)
--- Starts Newton-Raphson -------- 
#it |    gradient norm     |  residual load norm  |
---------------------------------------------------
  1 |   0.134159445656E-01 |   0.612370188774E-17 | 
  2 |   0.371405878766E-01 |   0.207104676539E-16 | 
  3 |   0.668707674562E-03 |   0.207542086232E-18 | 
  4 |   0.120980990213E-05 |   0.299251950680E-21 | 
  5 |   0.656169760904E-11 |   0.361196727392E-26 | 
 linear system solver:            1    time:   5.00400085E-03 seconds
---   Loop end   --------

-------  New time step: -------
 indtime:    2 factor:   0.47E+00 (delta:   0.27E+00)
--- Starts Newton-Raphson -------- 
#it |    gradient norm     |  residual load norm  |
---------------------------------------------------
  1 |   0.193503906761E-01 |   0.729300710064E-17 | 
  2 |   0.560448879535E-01 |   0.195827276197E-16 | 
  3 |   0.123638868363E-02 |   0.556945901834E-18 | 
  4 |   0.387333163243E-05 |   0.116767435777E-20 | 
  5 |   0.651562152966E-10 |   0.163256883516E-25 | 
 linear system solver:            1    time:   4.48699854E-03 seconds
---   Loop end   --------

-------  New time step: -------
 indtime:    3 factor:   0.73E+00 (delta:   0.27E+00)
--- Starts Newton-Raphson -------- 
#it |    gradient norm     |  residual load norm  |
---------------------------------------------------
  1 |   0.213966727814E-01 |   0.812549298955E-17 | 
  2 |   0.627522074871E-01 |   0.332441647205E-16 | 
  3 |   0.128510063193E-02 |   0.372833076419E-18 | 
  4 |   0.411804510961E-05 |   0.113233744743E-20 | 
  5 |   0.720892904071E-10 |   0.169831873709E-25 | 
 linear system solver:            1    time:   4.38899919E-03 seconds
---   Loop end   --------

-------  New time step: -------
 indtime:    4 factor:   0.10E+01 (delta:   0.27E+00)
--- Starts Newton-Raphson -------- 
#it |    gradient norm     |  residual load norm  |
---------------------------------------------------
  1 |   0.235520435462E-01 |   0.624624186222E-17 | 
  2 |   0.697852729115E-01 |   0.223404664773E-16 | 
  3 |   0.133280761750E-02 |   0.524033484032E-18 | 
  4 |   0.436388479103E-05 |   0.217358638487E-20 | 
  5 |   0.792851238956E-10 |   0.222906112046E-25 | 
 linear system solver:            1    time:   7.38699734E-03 seconds
---   Loop end   --------
 All files written, ldfem finished. Go and open the vtks with Paraview!
```



