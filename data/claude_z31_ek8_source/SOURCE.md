# Elkies-Kumar disc-8 ancillary data

Source: arXiv:1209.3527 (Elkies & Kumar, "K3 surfaces and equations for Hilbert
modular surfaces", Algebra & Number Theory 8 (2014) 2297-2411).
The arXiv e-print source tarball https://arxiv.org/e-print/1209.3527 IS the
ancillary data: one directory per fundamental discriminant (5, 8, 12, ..., 97),
each with the derivation transcript (<D>.txt) and the Igusa-Clebsch map
(igusa<D>.txt) in PARI/GP-readable form. Downloaded 2026-07-30.

Files here: 8.txt, igusa8.txt (directory "8/" of the tarball, gzip dated
2015-01-27, i.e. arXiv v3).

Disc-8 chart: rational plane (r,s); Y_-(8) is z^2 = 2*(16rs^2+32r^2 s-40rs-s
+16r^3+24r^2+12r+2) but the Igusa-Clebsch point depends on (r,s) only:
  [I2:I4:I6:I10] = [-24*B1/A1, -12*A, 96*(A/A1)*B1 - 36*B, -4*A1*B2]
  A1 = 2rs^2, A = -(9rs+4r^2+4r+1)/3, B1 = rs^2(3s+8r-2)/3,
  B = -(54r^2 s+81rs-16r^3-24r^2-12r-2)/27, B2 = r^2.
badlocus = r*s*(16rs^2+32r^2s-40rs-s+16r^3+24r^2+12r+2)
             *(27rs^2+36r^2s+18rs-72s+16r^3+48r^2+48r+16).

Note (public snapshot): the two data files are omitted from this
notebook precisely because they are citable published data; fetch them
from the arXiv tarball as described above.
