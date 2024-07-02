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
Here, the command-line argument ``<io_directory>`` is optional. In case it is not provided, ``io/sm``is used as default.

## The paran-file
``AlterAlterBBN`` expects a param-file with name ``param_file.dat`` in ``<io_directory>`` with at least one line that reads ``eta=<eta>``

## The cosmo-file
``AlterAlterBBN`` expects a cosmo-file with name ``cosmo_file.dat`` in ``<io_directory>`` with the same column structure and units as ``ACROPOLIS``, i.e.
* time in s
* temperature in MeV
* dTdt in MeV²
* neutrino temperature in MeV
* hubble rate in MeV
* nb/eta final ratio in MeV³

Note that this is different from how the original version of ``AlterAlterBBN`` handled the cosmo-file!
