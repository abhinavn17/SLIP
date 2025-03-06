# SLIP

Spectral Line Imaging Pipeline for radiop interferomtric data. Uses CASA's tclean to image and SOFIA to make moment maps. Tclean is parallelised using casampi.

## Prerequisites

SOFIA should be installed and running from the command line. To utiltise casampi, an mpirun implementation should be up and running properly. Install casampi using:

```bash
pip install mpi4py==3.1.4
pip install casampi==0.5.0
```

## How to use

1. To install, make a python environment and do:

```bash
pip install .
```

2. To run you need to make a config.ini file. It can be made using:

```bash 
slip_config -c <config.ini name> <msfile>
```

3. Edit the parameters in the config file to your use case. Additionally, a SOFIA .par is also copied to the folder. Parameters in that file can also be edited to your use case.

4. To run the selfcal loop, do:

```bash
slip <config.ini>
```

