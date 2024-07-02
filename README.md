# AlterAlterBBN
A modified (altered) version of AlterBBN (/https://alterbbn.hepforge.org/) that can be used for a broader class of BSM scenarios.

## How to compile and run
First create a ``bin`` direcory and then compile via
```
./build.sh
```
Afterwards run the code with
```
./bin/alteralterbbn <eta_1e10> <cosmo_file> <cosmo_file_lines> <abundance_file>
```
Here, all command-line arguments are optional and if any is not provided a default value is used.

## The cosmo-file
``AlterAlterBBN`` expects a cosmo-file with the same column structure and units as ``ACROPOLIS``, i.e.
* time in s
* temperature in MeV
* dTdt in MeV²
* neutrino temperature in MeV
* hubble rate in MeV
* nb/eta final ratio in MeV³

Note that this is different from how the original version of ``AlterAlterBBN`` handled the cosmo-file!
