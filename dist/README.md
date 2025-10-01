SNOVA
=======
This directory contains a tool to create SNOVA instances.

Use `make` to create `ref`, `opt`, and `avx2` directories.

In those directories use one of
```
OPT=REF make clean kat speed
make clean kat speed
OPT=AVX2 make clean kat speed
```

The default is `OPT=OPT`
