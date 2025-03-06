# SLIP

Spectral Line Imaging Pipeline for radio interferometric data. Accepts cleaned and calibrated MS files as input. Uses CASA's tclean to image and SOFIA to make moment maps. Tclean is parallelized using casampi.

## Prerequisites

SOFIA should be installed and running from the command line. To utilize casampi, an mpirun implementation should be up and running properly. Install casampi using:

```bash
pip install mpi4py==3.1.4
pip install casampi==0.5.0
```

## How to use

1. To install, create a Python environment and run:

```bash
pip install .
```

2. To run, you need to create a config.ini file. It can be created using:

```bash
slip_config -c <config.ini name> <msfile>
```

3. Edit the parameters in the config file to your use case. Additionally, a SOFIA .par file is also copied to the folder. Parameters in that file can also be edited to your use case.

4. To run the selfcal loop, execute:

```bash
slip <config.ini>
```
