# semiGPU - Semi-empirical methods for quantum chemistry on GPUs
This project aims to implement semi-empirical methods for quantum chemistry MNDO on GPUs. The main goal is to create a starting point for other semi-empirical methods. The project is based on the MNDO method, but the code is designed to be easily extended to other semi-empirical methods.

## Installation
The project is written in C++ and uses OpenACC for GPU acceleration. The code is tested with the HPC SDK 23.3 from NVIDIA. The project uses Make for building the code. 

The following commands build the project without OpenACC support:
```bash
git clone https://github.com/LuisAlfredoNu/semiGPU.git
cd semiGPU
make
```

The following commands build the project with OpenACC support for multicore:
```bash
make OPENACC=1 DEVICE=multicore
```

The following commands build the project with OpenACC support for GPU:
```bash
make OPENACC=1 DEVICE=tesla
```

## Usage
The executable file is generated in the `App` directory. The following command runs the project:
```bash
./App/semiGPU geometry.xyz
```
The input file is the following:
```bash
XYZ file : Geometries4Test/04_pyridine.xyz
Read XYZ file successful
Geometry :
  N   -1.4026   -0.0003         0
  C    1.3927    0.0002         0
  C    0.6945    1.2027    0.0001
  C     0.695   -1.2025         0
  C     -0.69    1.1487   -0.0001
  C   -0.6896   -1.1489         0
  H    2.4786    0.0004   -0.0001
  H    1.2175    2.1521    0.0001
  H    1.2184   -2.1517   -0.0001
  H   -1.2813    2.0585   -0.0001
  H   -1.2804    -2.059         0
********************************************************************************
Correct Alloc: all2CenterIntegrals
Correct Alloc: hcoreMatrix
Correct Alloc: SCFData
********************************************************************************
SCF step : 0  Energy = -3091.23
SCF step : 1  Energy = -3359.43  LapackTime = 49  PmatTime = 0  FmatTime = 6  EnerTime = 0  StepTime = 55
SCF step : 2  Energy = -3360.78  LapackTime = 12  PmatTime = 0  FmatTime = 4  EnerTime = 0  StepTime = 17
SCF step : 3  Energy = -3360.87  LapackTime = 9  PmatTime = 0  FmatTime = 6  EnerTime = 0  StepTime = 16
SCF step : 4  Energy = -3360.88  LapackTime = 0  PmatTime = 0  FmatTime = 7  EnerTime = 0  StepTime = 7
SCF step : 5  Energy = -3360.88  LapackTime = 1  PmatTime = 0  FmatTime = 5  EnerTime = 0  StepTime = 6
********************************************************************************
  Electronic Energy =   -3360.88 eV
********************************************************************************
Core-Core Repulsion =    2444.36 eV
********************************************************************************
       Total Energy =   -916.525 eV
********************************************************************************
Dealloc array: SCFData
Dealloc array: hcoreMatrix
Dealloc array: all2CenterIntegral
********************************************************************************
Time taken for calculation: 0.116 s
```

## Tests
The project has a general test for the MNDO method which compares the results with the MNDO results from the MOPAC package. The test are the following:

- `test_twoCenterIntegral.x` : Test for the two-center integrals
- `test_hcore.x` : Test for the core Hamiltonian matrix
- `test_fockmatrix.x` : Test for the Fock matrix
- `test_corecorerepulsion.x` : Test for the core-core repulsion
- `test_scfcalculation.x` : Test for the SCF calculation

The following command runs the tests:
```bash
cd test
bash fullTester.bsh ../Geometries4Test/<geometry>.xyz
```

The output of the tests prints all the values of the matrices and the energies. Also, highlights the differences with the MOPAC. At the end, the test prints the following message:
```bash
...
.
.
.
...
********************************************************************************
Compute Electronic Energy = -873.482
    Ref Electronic Energy = -873.483
                    Error < 0.001 %
********************************************************************************
Test pass
===================================================
 test_twoCenterIntegral.x : Pass
             test_hcore.x : Pass
        test_fockmatrix.x : Fail
 test_corecorerepulsion.x : Pass
                    Error : 0.0001
    test_scfcalculation.x : Pass
                    Error : 0.001
===================================================
```

## TODO
- [] Implement a geometry optimizer
- [] Implement a PM6 method
- [] Implement a PM7 method

## About This Project

This project was part of my Master's degree project in Theoretical Chemistry at the *Universidad Autónoma Metropolitana Unidad Iztapalapa*. 

<p align="right">
    <img src="./img/uam.jpg" alt="UAMI" width="auto" height="100"> <!-- adjust width and height as needed -->
</p>

