#!/bin/env python

import numpy as np
MuonEtaLut = (
(
( -2 , -2 , -2 , -2 , -2 , -2 , -2 , -2 ),
( -5 , -5 , -5 , -5 , -5 , -5 , -5 , -5 ),
( -8 , -8 , -8 , -8 , -8 , -8 , -8 , -8 ),
( -12 , -12 , -11 , -11 , -11 , -11 , -12 , -12 ),
( -15 , -15 , -15 , -15 , -15 , -15 , -15 , -15 ),
( -18 , -18 , -18 , -18 , -18 , -18 , -18 , -18 ),
( -22 , -22 , -22 , -22 , -22 , -22 , -22 , -22 ) ),
(
( -2 , -2 , -2 , -2 , -2 , -2 , -2 , -2 ),
( -5 , -5 , -5 , -5 , -5 , -5 , -5 , -5 ),
( -8 , -8 , -8 , -8 , -8 , -8 , -8 , -8 ),
( -12 , -12 , -11 , -11 , -11 , -11 , -12 , -12 ),
( -15 , -15 , -15 , -15 , -15 , -15 , -15 , -15 ),
( -18 , -18 , -18 , -18 , -18 , -18 , -18 , -18 ),
( -22 , -22 , -22 , -22 , -22 , -22 , -22 , -22 ) ),
(
( -2 , -2 , -2 , -2 , -2 , -2 , -2 , -2 ),
( -5 , -5 , -5 , -5 , -5 , -5 , -5 , -5 ),
( -8 , -8 , -8 , -8 , -8 , -8 , -8 , -8 ),
( -12 , -12 , -11 , -11 , -11 , -11 , -12 , -12 ),
( -15 , -15 , -15 , -15 , -15 , -15 , -15 , -15 ),
( -18 , -18 , -18 , -18 , -18 , -18 , -18 , -18 ),
( -22 , -22 , -22 , -22 , -22 , -22 , -22 , -22 ) ),
(
( -2 , -2 , -2 , -2 , -2 , -2 , -2 , -2 ),
( -5 , -5 , -5 , -5 , -5 , -5 , -5 , -5 ),
( -8 , -8 , -8 , -8 , -8 , -8 , -8 , -8 ),
( -12 , -12 , -11 , -11 , -11 , -11 , -12 , -12 ),
( -15 , -15 , -15 , -15 , -15 , -15 , -15 , -15 ),
( -18 , -18 , -18 , -18 , -18 , -18 , -18 , -18 ),
( -22 , -22 , -22 , -22 , -22 , -22 , -22 , -22 ) ),
(
( -2 , -2 , -2 , -2 , -2 , -2 , -2 , -2 ),
( -5 , -5 , -5 , -5 , -5 , -5 , -5 , -5 ),
( -8 , -8 , -8 , -8 , -8 , -8 , -8 , -8 ),
( -12 , -12 , -11 , -11 , -11 , -11 , -12 , -12 ),
( -15 , -15 , -15 , -15 , -15 , -15 , -15 , -15 ),
( -18 , -18 , -18 , -18 , -18 , -18 , -18 , -18 ),
( -22 , -22 , -22 , -22 , -22 , -22 , -22 , -22 ) ),
(
( -2 , -2 , -2 , -2 , -2 , -2 , -2 , -2 ),
( -5 , -5 , -5 , -5 , -5 , -5 , -6 , -6 ),
( -8 , -8 , -8 , -8 , -8 , -8 , -8 , -8 ),
( -12 , -12 , -11 , -11 , -11 , -11 , -12 , -12 ),
( -15 , -15 , -15 , -15 , -15 , -15 , -15 , -15 ),
( -18 , -18 , -18 , -18 , -18 , -18 , -18 , -18 ),
( -22 , -22 , -22 , -22 , -22 , -22 , -22 , -22 ) ),
(
( -2 , -2 , -2 , -2 , -2 , -2 , -2 , -2 ),
( -6 , -6 , -5 , -5 , -5 , -5 , -6 , -6 ),
( -8 , -8 , -9 , -9 , -9 , -9 , -8 , -8 ),
( -12 , -12 , -11 , -11 , -11 , -11 , -12 , -12 ),
( -15 , -15 , -15 , -15 , -15 , -15 , -15 , -15 ),
( -18 , -18 , -18 , -18 , -18 , -18 , -18 , -18 ),
( -22 , -22 , -22 , -22 , -22 , -22 , -22 , -22 ) ),
(
( -2 , -2 , -2 , -2 , -2 , -2 , -2 , -2 ),
( -6 , -6 , -5 , -5 , -5 , -5 , -5 , -5 ),
( -8 , -8 , -8 , -8 , -8 , -8 , -8 , -8 ),
( -12 , -12 , -11 , -11 , -11 , -11 , -12 , -12 ),
( -15 , -15 , -15 , -15 , -15 , -15 , -15 , -15 ),
( -18 , -18 , -18 , -18 , -18 , -18 , -18 , -18 ),
( -22 , -22 , -22 , -22 , -22 , -22 , -22 , -22 ) ),
(
( 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 ),
( 5 , 5 , 5 , 5 , 5 , 5 , 5 , 5 ),
( 8 , 8 , 8 , 8 , 8 , 8 , 8 , 8 ),
( 12 , 12 , 11 , 11 , 11 , 11 , 12 , 12 ),
( 15 , 15 , 15 , 15 , 15 , 15 , 15 , 15 ),
( 18 , 18 , 18 , 18 , 18 , 18 , 18 , 18 ),
( 22 , 22 , 22 , 22 , 22 , 22 , 22 , 22 ) ),
(
( 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 ),
( 5 , 5 , 5 , 5 , 5 , 5 , 5 , 5 ),
( 8 , 8 , 8 , 8 , 8 , 8 , 8 , 8 ),
( 12 , 12 , 11 , 11 , 11 , 11 , 12 , 12 ),
( 15 , 15 , 15 , 15 , 15 , 15 , 15 , 15 ),
( 18 , 18 , 18 , 18 , 18 , 18 , 18 , 18 ),
( 22 , 22 , 22 , 22 , 22 , 22 , 22 , 22 ) ),
(
( 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 ),
( 5 , 5 , 5 , 5 , 5 , 5 , 5 , 5 ),
( 8 , 8 , 8 , 8 , 8 , 8 , 8 , 8 ),
( 12 , 12 , 11 , 11 , 11 , 11 , 12 , 12 ),
( 15 , 15 , 15 , 15 , 15 , 15 , 15 , 15 ),
( 18 , 18 , 18 , 18 , 18 , 18 , 18 , 18 ),
( 22 , 22 , 22 , 22 , 22 , 22 , 22 , 22 ) ),
(
( 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 ),
( 5 , 5 , 5 , 5 , 5 , 5 , 5 , 5 ),
( 8 , 8 , 8 , 8 , 8 , 8 , 8 , 8 ),
( 12 , 12 , 11 , 11 , 11 , 11 , 12 , 12 ),
( 15 , 15 , 15 , 15 , 15 , 15 , 15 , 15 ),
( 18 , 18 , 18 , 18 , 18 , 18 , 18 , 18 ),
( 22 , 22 , 22 , 22 , 22 , 22 , 22 , 22 ) ),
(
( 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 ),
( 5 , 5 , 5 , 5 , 5 , 5 , 5 , 5 ),
( 8 , 8 , 8 , 8 , 8 , 8 , 8 , 8 ),
( 12 , 12 , 11 , 11 , 11 , 11 , 12 , 12 ),
( 15 , 15 , 15 , 15 , 15 , 15 , 15 , 15 ),
( 18 , 18 , 18 , 18 , 18 , 18 , 18 , 18 ),
( 22 , 22 , 22 , 22 , 22 , 22 , 22 , 22 ) ),
(
( 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 ),
( 5 , 5 , 5 , 5 , 5 , 5 , 6 , 6 ),
( 8 , 8 , 8 , 8 , 8 , 8 , 8 , 8 ),
( 12 , 12 , 11 , 11 , 11 , 11 , 12 , 12 ),
( 15 , 15 , 15 , 15 , 15 , 15 , 15 , 15 ),
( 18 , 18 , 18 , 18 , 18 , 18 , 18 , 18 ),
( 22 , 22 , 22 , 22 , 22 , 22 , 22 , 22 ) ),
(
( 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 ),
( 6 , 6 , 5 , 5 , 5 , 5 , 6 , 6 ),
( 8 , 8 , 9 , 9 , 9 , 9 , 8 , 8 ),
( 12 , 12 , 11 , 11 , 11 , 11 , 12 , 12 ),
( 15 , 15 , 15 , 15 , 15 , 15 , 15 , 15 ),
( 18 , 18 , 18 , 18 , 18 , 18 , 18 , 18 ),
( 22 , 22 , 22 , 22 , 22 , 22 , 22 , 22 ) ),
(
( 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 ),
( 6 , 6 , 5 , 5 , 5 , 5 , 5 , 5 ),
( 8 , 8 , 8 , 8 , 8 , 8 , 8 , 8 ),
( 12 , 12 , 11 , 11 , 11 , 11 , 12 , 12 ),
( 15 , 15 , 15 , 15 , 15 , 15 , 15 , 15 ),
( 18 , 18 , 18 , 18 , 18 , 18 , 18 , 18 ),
( 22 , 22 , 22 , 22 , 22 , 22 , 22 , 22 ) )
)

MuonPhiLut = (
(
( 59 , 60 , 61 , 62 , 1 , 2 , 3 , 3 ),
( 59 , 60 , 61 , 62 , 1 , 2 , 3 , 3 ),
( 59 , 60 , 61 , 62 , 1 , 2 , 3 , 3 ),
( 59 , 60 , 61 , 62 , 1 , 2 , 2 , 3 ),
( 59 , 60 , 61 , 62 , 0 , 1 , 2 , 3 ),
( 59 , 60 , 61 , 62 , 0 , 1 , 2 , 3 ),
( 61 , 62 , 31 , 1 , 2 , 3 , 4 , 5 ) ),
(
( 4 , 5 , 6 , 7 , 8 , 10 , 10 , 11 ),
( 4 , 5 , 6 , 7 , 8 , 10 , 10 , 11 ),
( 4 , 5 , 6 , 7 , 8 , 10 , 10 , 11 ),
( 4 , 5 , 6 , 7 , 8 , 9 , 10 , 11 ),
( 4 , 5 , 6 , 7 , 8 , 9 , 10 , 11 ),
( 4 , 5 , 6 , 7 , 8 , 9 , 10 , 11 ),
( 6 , 7 , 8 , 9 , 10 , 11 , 12 , 13 ) ),
(
( 12 , 13 , 14 , 15 , 16 , 17 , 18 , 19 ),
( 12 , 13 , 14 , 15 , 16 , 17 , 18 , 19 ),
( 12 , 13 , 14 , 15 , 16 , 17 , 18 , 19 ),
( 12 , 13 , 14 , 15 , 16 , 17 , 18 , 19 ),
( 12 , 13 , 14 , 15 , 16 , 17 , 18 , 19 ),
( 12 , 13 , 14 , 15 , 16 , 17 , 18 , 19 ),
( 14 , 15 , 16 , 17 , 18 , 19 , 20 , 21 ) ),
(
( 20 , 21 , 22 , 23 , 24 , 25 , 26 , 27 ),
( 20 , 21 , 22 , 23 , 24 , 25 , 26 , 27 ),
( 20 , 21 , 22 , 23 , 24 , 25 , 26 , 27 ),
( 20 , 21 , 22 , 23 , 24 , 25 , 26 , 27 ),
( 20 , 21 , 22 , 23 , 24 , 25 , 26 , 27 ),
( 20 , 21 , 22 , 23 , 24 , 25 , 26 , 27 ),
( 22 , 23 , 24 , 25 , 26 , 27 , 27 , 28 ) ),
(
( 28 , 29 , 30 , 31 , 32 , 33 , 34 , 35 ),
( 28 , 29 , 30 , 31 , 32 , 33 , 34 , 35 ),
( 28 , 29 , 30 , 31 , 32 , 33 , 34 , 35 ),
( 28 , 29 , 30 , 31 , 32 , 33 , 34 , 35 ),
( 28 , 29 , 30 , 31 , 32 , 33 , 34 , 35 ),
( 28 , 29 , 30 , 31 , 32 , 33 , 34 , 35 ),
( 29 , 30 , 31 , 32 , 33 , 34 , 35 , 36 ) ),
(
( 36 , 37 , 38 , 39 , 40 , 41 , 42 , 43 ),
( 36 , 37 , 38 , 39 , 40 , 41 , 42 , 43 ),
( 36 , 37 , 38 , 39 , 40 , 41 , 42 , 43 ),
( 36 , 37 , 38 , 39 , 40 , 41 , 42 , 43 ),
( 36 , 37 , 38 , 39 , 40 , 41 , 42 , 43 ),
( 36 , 37 , 38 , 39 , 40 , 41 , 42 , 43 ),
( 37 , 38 , 39 , 40 , 41 , 42 , 43 , 44 ) ),
(
( 44 , 45 , 45 , 47 , 48 , 49 , 50 , 51 ),
( 44 , 45 , 45 , 47 , 48 , 49 , 50 , 51 ),
( 44 , 45 , 45 , 47 , 48 , 49 , 50 , 51 ),
( 44 , 45 , 46 , 47 , 48 , 49 , 50 , 51 ),
( 44 , 45 , 46 , 47 , 48 , 49 , 50 , 51 ),
( 44 , 45 , 46 , 47 , 48 , 49 , 50 , 51 ),
( 45 , 46 , 47 , 48 , 49 , 50 , 51 , 52 ) ),
(
( 51 , 52 , 53 , 54 , 56 , 57 , 58 , 58 ),
( 51 , 52 , 53 , 54 , 56 , 57 , 58 , 58 ),
( 51 , 52 , 53 , 54 , 56 , 57 , 58 , 58 ),
( 52 , 53 , 53 , 54 , 56 , 57 , 57 , 58 ),
( 52 , 53 , 54 , 54 , 55 , 56 , 57 , 58 ),
( 52 , 53 , 54 , 54 , 55 , 56 , 57 , 58 ),
( 53 , 54 , 55 , 56 , 57 , 58 , 59 , 60 ) ),
(
( 59 , 60 , 61 , 62 , 1 , 2 , 3 , 3 ),
( 59 , 60 , 61 , 62 , 1 , 2 , 3 , 3 ),
( 59 , 60 , 61 , 62 , 1 , 2 , 3 , 3 ),
( 59 , 60 , 61 , 62 , 1 , 2 , 2 , 3 ),
( 59 , 60 , 61 , 62 , 0 , 1 , 2 , 3 ),
( 59 , 60 , 61 , 62 , 0 , 1 , 2 , 3 ),
( 61 , 62 , 31 , 1 , 2 , 3 , 4 , 5 ) ),
(
( 4 , 5 , 6 , 7 , 8 , 10 , 10 , 11 ),
( 4 , 5 , 6 , 7 , 8 , 10 , 10 , 11 ),
( 4 , 5 , 6 , 7 , 8 , 10 , 10 , 11 ),
( 4 , 5 , 6 , 7 , 8 , 9 , 10 , 11 ),
( 4 , 5 , 6 , 7 , 8 , 9 , 10 , 11 ),
( 4 , 5 , 6 , 7 , 8 , 9 , 10 , 11 ),
( 6 , 7 , 8 , 9 , 10 , 11 , 12 , 13 ) ),
(
( 12 , 13 , 14 , 15 , 16 , 17 , 18 , 19 ),
( 12 , 13 , 14 , 15 , 16 , 17 , 18 , 19 ),
( 12 , 13 , 14 , 15 , 16 , 17 , 18 , 19 ),
( 12 , 13 , 14 , 15 , 16 , 17 , 18 , 19 ),
( 12 , 13 , 14 , 15 , 16 , 17 , 18 , 19 ),
( 12 , 13 , 14 , 15 , 16 , 17 , 18 , 19 ),
( 14 , 15 , 16 , 17 , 18 , 19 , 20 , 21 ) ),
(
( 20 , 21 , 22 , 23 , 24 , 25 , 26 , 27 ),
( 20 , 21 , 22 , 23 , 24 , 25 , 26 , 27 ),
( 20 , 21 , 22 , 23 , 24 , 25 , 26 , 27 ),
( 20 , 21 , 22 , 23 , 24 , 25 , 26 , 27 ),
( 20 , 21 , 22 , 23 , 24 , 25 , 26 , 27 ),
( 20 , 21 , 22 , 23 , 24 , 25 , 26 , 27 ),
( 22 , 23 , 24 , 25 , 26 , 27 , 27 , 28 ) ),
(
( 28 , 29 , 30 , 31 , 32 , 33 , 34 , 35 ),
( 28 , 29 , 30 , 31 , 32 , 33 , 34 , 35 ),
( 28 , 29 , 30 , 31 , 32 , 33 , 34 , 35 ),
( 28 , 29 , 30 , 31 , 32 , 33 , 34 , 35 ),
( 28 , 29 , 30 , 31 , 32 , 33 , 34 , 35 ),
( 28 , 29 , 30 , 31 , 32 , 33 , 34 , 35 ),
( 29 , 30 , 31 , 32 , 33 , 34 , 35 , 36 ) ),
(
( 36 , 37 , 38 , 39 , 40 , 41 , 42 , 43 ),
( 36 , 37 , 38 , 39 , 40 , 41 , 42 , 43 ),
( 36 , 37 , 38 , 39 , 40 , 41 , 42 , 43 ),
( 36 , 37 , 38 , 39 , 40 , 41 , 42 , 43 ),
( 36 , 37 , 38 , 39 , 40 , 41 , 42 , 43 ),
( 36 , 37 , 38 , 39 , 40 , 41 , 42 , 43 ),
( 37 , 38 , 39 , 40 , 41 , 42 , 43 , 44 ) ),
(
( 44 , 45 , 45 , 47 , 48 , 49 , 50 , 51 ),
( 44 , 45 , 45 , 47 , 48 , 49 , 50 , 51 ),
( 44 , 45 , 45 , 47 , 48 , 49 , 50 , 51 ),
( 44 , 45 , 46 , 47 , 48 , 49 , 50 , 51 ),
( 44 , 45 , 46 , 47 , 48 , 49 , 50 , 51 ),
( 44 , 45 , 46 , 47 , 48 , 49 , 50 , 51 ),
( 45 , 46 , 47 , 48 , 49 , 50 , 51 , 52 ) ),
(
( 51 , 52 , 53 , 54 , 56 , 57 , 58 , 58 ),
( 51 , 52 , 53 , 54 , 56 , 57 , 58 , 58 ),
( 51 , 52 , 53 , 54 , 56 , 57 , 58 , 58 ),
( 52 , 53 , 53 , 54 , 56 , 57 , 57 , 58 ),
( 52 , 53 , 54 , 54 , 55 , 56 , 57 , 58 ),
( 52 , 53 , 54 , 54 , 55 , 56 , 57 , 58 ),
( 53 , 54 , 55 , 56 , 57 , 58 , 59 , 60 ) )
)


MuonPtLut =( 4 , 6 , 10 , 11 , 15 , 20)



MuonEtaLut = np.array(MuonEtaLut)
MuonPhiLut = np.array(MuonPhiLut)
MuonPtLut = np.array(MuonPtLut)

print "eta: ",MuonEtaLut.shape
print "phi: ",MuonPhiLut.shape
print "pt: ",MuonPtLut.shape

print("looping with 'octant' in 0-7")
for side in range(2):
    for part in range(8):
        print("{0:d}".format(side*8+part))

print("looping with 'octant'/ in 0-15")
for side in range(2):
    for part in range(16):
        print("{0:d}".format(side*8+part/2))
