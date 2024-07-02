# AlterAlterBBN
A modified (altered) version of AlterBBN (https://alterbbn.hepforge.org/) that can be used for a broader class of BSM scenarios.

## How to compile and run
Compile via
```
./build.sh
```
and afterwards run with
```
./bin/alteralterbbn <io_directory>
```
Here, the command-line argument ``<io_directory>`` is optional. In case it is not provided, ``io/sm`` is used as default.

An example output would be
```
INFO: Using 'io/sm' to read and write data
INFO: Using eta = 6.13700e-10 from 'io/sm/param_file.dat'
INFO: Using cosmological data from 'io/sm/cosmo_file.dat'
INFO: Running nucleosynthesis...Done!
INFO: The final abundances are:

mean         high         low
0.000000e+00 0.000000e+00 0.000000e+00
7.532623e-01 7.533308e-01 7.531942e-01
1.882976e-05 1.842939e-05 1.924010e-05
6.028483e-08 5.846802e-08 6.218854e-08
7.718020e-06 7.778272e-06 7.663938e-06
6.166916e-02 6.165221e-02 6.168603e-02
8.273446e-15 2.691357e-14 1.282042e-15
2.201507e-11 2.229536e-11 2.161805e-11
3.765606e-10 4.021328e-10 3.532401e-10

INFO: The final abundances have been written to 'io/sm/abundance_file.dat'
```

## The param-file
``AlterAlterBBN`` expects a param-file with name ``param_file.dat`` in ``<io_directory>`` with at least one line that reads ``eta=<eta>``, e.g. ``eta=6.137e-10``

## The cosmo-file
``AlterAlterBBN`` expects a cosmo-file with name ``cosmo_file.dat`` in ``<io_directory>`` with the same column structure and units as ``ACROPOLIS``, i.e.
* time in s
* temperature in MeV
* dTdt in MeV²
* neutrino temperature in MeV
* hubble rate in MeV
* nb/eta final ratio in MeV³

Note that this is different from how the original version of ``AlterAlterBBN`` handled the cosmo-file!
