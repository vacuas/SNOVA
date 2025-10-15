SNOVA
=======
This directory contains the official constant-time implementation of the SNOVA signature scheme.


Building
-------

Building SNOVA requires a C compiler and `make`. The SNOVA parameters are set in `snova_params.h`. The SNOVA parameters can also be changed by the command line parameters of the `make` command, e.g.
```
make clean all P="-D SNOVA_v=37 -D SNOVA_o=8 -D SNOVA_q=19 -D SNOVA_l=4"
```
An example command line build for $q=16$, using AES-CTR for the public key expansion, is
```
make clean all P="-D SNOVA_v=24 -D SNOVA_o=5 -D SNOVA_q=16 -D SNOVA_l=4 -D AESCTR"
```

Available optimization options are:
1. Use `make OPT=REF` to build the reference implementation.
2. Use `make OPT=OPT` for the plain-C optimized version.
3. Use `make OPT=AVX2` for a version using AVX2 vectorization instructions.


Recommended parameters
-------

Our recommended parameter sets for odd prime $q$ all have a matrix rank $l=4$. The following are the recommended parameters:

| SL |        Name      |  v |  o |   q |  l |  sk size |  pk size |  sign size |
|----|------------------|----|----|-----|----|----------|----------|------------|
|  1 |  SNOVA_24_5_23_4 | 24 |  5 |  23 |  4 |      48  |      616 |        282 |
|  1 |  SNOVA_24_5_16_4 | 24 |  5 |  16 |  4 |      48  |     1016 |        248 |
|  1 | SNOVA_43_17_16_2 | 43 | 17 |  16 |  2 |      48  |     9842 |        136 |
|  3 |  SNOVA_37_8_19_4 | 37 |  8 |  19 |  4 |      48  |     2269 |        400 |
|  3 |  SNOVA_37_8_16_4 | 37 |  8 |  16 |  4 |      48  |     4112 |        376 |
|  3 | SNOVA_69_25_16_2 | 69 | 25 |  16 |  2 |      48  |    31266 |        204 |
|  5 | SNOVA_60_10_23_4 | 60 | 10 |  23 |  4 |      48  |     4702 |        656 |
|  5 | SNOVA_60_10_16_4 | 60 | 10 |  16 |  4 |      48  |     8016 |        576 |
|  5 | SNOVA_99_25_16_2 | 99 | 25 |  16 |  2 |      48  |    71890 |        280 |


Performance
-------

While the optimized versions use only C statements, the compiler will actually vectorize the code. We found that the level of vectorization that the compiler produces depends significantly on the compiler used, and also the version of the compiler used. The best performance was obtained using gcc version 15.2.1 20250813 on Arch Linux.

We have observed the following cycle counts for the optimized `OPT=OPT` version for SNOVA_24_5_23_4:

| Compiler | version |   Genkey  |     Sign   |  Verify  |
|----------|---------|-----------|------------|----------|
|   gcc    |  15.2.1 |   456,966 |    758,046 |  356,460 |
|   gcc    |  13.3.0 |   459,138 |    977,592 |  364,516 |
|   gcc    |  11.4.1 |   504,406 |  1,449,538 |  380,404 |
|  clang   |  20.1.8 |   953,932 |  1,960,312 |  425,523 |

These cycle counts were collected on an Intel(R) Core(TM) Ultra 7 155H (Meteor Lake) laptop running Arch Linux and using gcc version 15.2.1 20250813. We report the median over 2048 benchmark runs.

The performance of the sign on other compilers can be made to be close to that of gcc 15.1.1 by using explicit AVX2 vectorization instructions. Benchmark results for the vectorized `OPT=AVX2` version of SNOVA_24_5_23_4 are:

| Compiler | version |   Genkey  |     Sign   |  Verify  |
|----------|---------|-----------|------------|----------|
|   gcc    |  15.2.1 |   456,966 |    758,046 |  356,460 |
|   gcc    |  13.3.0 |   459,376 |    949,297 |  366,787 |
|   gcc    |  11.4.1 |   503,318 |  1,068,568 |  379,610 |
|  clang   |  20.1.8 |   970,835 |  1,022,163 |  422,670 |
