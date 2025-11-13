SNOVA in SageMath
=============

This directory contains a sage script `snova.sage` implementing SNOVA. It generates correct KAT files. The scripts supports the Rectangular SNOVA variant described in the appendix of `SNOVA_21.pdf`.

In addition an `emat_snova.sage` script is provided that uses an alternative representation of the SNOVA public map. This version will create valid SNOVA KAT files as well as recreate the MAYO-1 KAT files if so configured. This script demonstrates that both SNOVA and MAYO are instances of Rectangular SNOVA.

The scripts `smat_q.sage` and `smat16.sage` have been used to find and verify $S$ matrices that are irreducible for $l \in \{ 2, 3, 4, 5 \}$.
