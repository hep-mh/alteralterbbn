# AlterAlterBBN

![Language: C](https://img.shields.io/badge/Language-C-%23D0D0D0.svg?style=flat-square)

A modified (altered) version of ``AlterBBN`` (https://alterbbn.hepforge.org/) that can be used to calculate BBN constraints for a broader class of BSM scenarios.

| Input          | Output             |
| -------------- | ------------------ |
| param_file.dat | abundance_file.dat |
| cosmo_file.dat |                    |

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

An example output of running ``./bin/alteralterbbn`` would be
```
AlterAlterBBN v2.0

INFO   : Using the directory 'io/sm' to read and write data.
INFO   : Using eta = 6.13700e-10 from 'io/sm/param_file.dat'.
INFO   : Using cosmological data from 'io/sm/cosmo_file.dat'.
INFO   : Running nucleosynthesis...Done!
INFO   : The final abundances are:

     |     mean     |     high     |     low
-------------------------------------------------
   n | 0.000000e+00 | 0.000000e+00 | 0.000000e+00
   p | 7.532623e-01 | 7.533308e-01 | 7.531942e-01
  H2 | 1.882976e-05 | 1.842938e-05 | 1.924009e-05
  H3 | 6.028469e-08 | 5.846805e-08 | 6.218953e-08
 He3 | 7.718019e-06 | 7.778271e-06 | 7.663937e-06
 He4 | 6.166916e-02 | 6.165221e-02 | 6.168603e-02
 Li6 | 8.273444e-15 | 2.691357e-14 | 1.282041e-15
 Li7 | 2.201511e-11 | 2.229534e-11 | 2.161780e-11
 Be7 | 3.765608e-10 | 4.021330e-10 | 3.532403e-10

INFO   : The final abundances have been written to 'io/sm/abundance_file.dat'.
```

## The param-file
``AlterAlterBBN`` expects a param-file with name ``param_file.dat`` in ``<io_directory>`` with at least one line that reads ``eta=<eta>``, e.g. ``eta=6.137e-10``

## The cosmo-file
``AlterAlterBBN`` expects a cosmo-file with name ``cosmo_file.dat`` in ``<io_directory>`` with the column structure
* ``t`` in s
* ``T`` in MeV
* ``dTdt`` in MeV²
* ``Tnu`` in MeV
* ``H`` in MeV
* ``nb_etaf`` in MeV³

Note that this is different from how the original version of ``AlterAlterBBN`` handled the cosmo-file!
