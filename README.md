# SLIP

Spectral Line Imaging Pipeline

## How to use

1. To install, make a python environment and do:

```bash
pip install .
```

2. To run you need to make a config.ini file. It can be made using:

```bash 
slip_config -c <config.ini name> <msfile>
```

3. Edit the parameters in the config file to your use case.

4. To run the selfcal loop, do:

```bash
slip <config.ini>
```

