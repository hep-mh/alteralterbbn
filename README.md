# alteralterbbn
A modified (altered) version of AlterBBN (/https://alterbbn.hepforge.org/) that can be used for a broader class of BSM scenarios.

First create a ``bin`` direcory and then compile via
```
clang src/*.c -lm -O3 -o bin/alteralterbbn
```
or
```
./build.sh
```
Afterwards run the code with
```
./bin/alteralterbbn <eta_1e10> <cosmo_file> <cosmo_file_lines> <abundance_file>
```
Here, all command-line arguments are optional
