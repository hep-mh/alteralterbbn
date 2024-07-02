# AlterAlterBBN
A modified (altered) version of AlterBBN (https://alterbbn.hepforge.org/) that can be used for a broader class of BSM scenarios.

## How to compile and run
Compile via
```
./build.sh
```
and afterwards run with
```
./bin/alteralterbbn <eta_1e10> <io_directory>
```
Here, all command-line arguments are optional and if any is not provided a default value is used.

## The cosmo-file
``AlterAlterBBN`` expects a cosmo-file with name ``cosmo_file.dat`` in ``<io_directory>`` with the same column structure and units as ``ACROPOLIS``, i.e.
* time in s
* temperature in MeV
* dTdt in MeV²
* neutrino temperature in MeV
* hubble rate in MeV
* nb/eta final ratio in MeV³

Note that this is different from how the original version of ``AlterAlterBBN`` handled the cosmo-file!
