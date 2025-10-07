SNOVA
=======
This repository contains the latest official Reference, Optimized, and AVX2 implementations in C of the SNOVA signature scheme.

The `src` directory contains the updated version of SNOVA. The `dist` directory contains a makefile that will create subdirectories containing all the recommended SNOVA instances.

## Build instructions

Building SNOVA requires a C compiler, `make` and the OpenSSL library (for AES). The SNOVA parameters are set in `src/snova_params.h`. The SNOVA parameters can also be changed by the command line parameters of the `make` command.

To create the recommended SNOVA instances enter `make` in the `dist` folder.
```
cd dist
make
```

To create the KAT files and their digests, starting in the `dist` folder
```
cd ref
make kat
make digest
```

To obtain benchmarks
```
cd avx2
OPT=AVX2 make speed
OPT=AVX2 make speed
```

Use
```
OPT=OPT make speed
```
when AVX2 instructions are not available.
