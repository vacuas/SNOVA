SNOVA
=======
This repository contains the latest official Reference, Optimized, and AVX2 implementations in C of the SNOVA signature scheme.

The `src` directory contains the updated version of SNOVA. The `./round2_src` subdirectory contains the version as submitted to NIST round 2. 

The `dist` directory contains a makefile that will create subdirectories containing all the recommended SNOVA instances. The makefile in the top level directory will also create subdirectories in the `dist` directory.

## Build instructions

Building SNOVA requires a C compiler, `make` and the OpenSSL library (for AES). The SNOVA parameters are set in `src/snova_params.h`. The SNOVA parameters can also be changed by the command line parameters of the `make` command.

To create the recommended SNOVA instances enter `make` in this folder or in the `dist`folder. This will create the instances in subdirectories of `dist`.

To create the KAT files and their digests
```
cd dist/ref
make kat
make digest
```

To obtain benchmarks
```
cd dist/avx2
OPT=AVX2 make speed
OPT=AVX2 make speed
```
