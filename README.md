--------------------------------------------------------------------------------
Overview
--------------------------------------------------------------------------------
bls48 is an Optimal Ate Pairing(OAP) on BLS Curve implementation in C++.
This implementation is based on the drafts below;

- Y. Kiyomura, et al. "Secure and Eifficient Apiaring at 256-Bit Secure Level" (2017)
- C. Costello, K. Lauter, M. Naehrig, "Attractive Subfamilies of BLS Curves for Implementing High-Security Pairings." (2011)
- draft-by-Kato(not submitted)

Furthermore, this implementation uses libsnark library.
- https://github.com/scipr-lab/libsnark/tree/deprecated-master  (branch : deprecated-master)

We are motivated by the security consideration for pairing on any curve give by T. Kim.
- T. Kim et al. "The extended tower number field sieve: A new complexity for the medium prime case. In Advances in Cryptology", (2016).

Kim's attack made us to update the security parameters or to use alternative curves in OAP.
BLS-48 cuerve is one of the pairing friendly curves, and recommended parameters for realizing 256-bit security is written in [Kiyomura].

--------------------------------------------------------------------------------
Build instructions
--------------------------------------------------------------------------------
This implementation relies on the following:

- C++ build environment
- GMP for certain bit-integer arithmetic
- libprocps for reporting memory usage
- GTest for some of the unit tests

We have tested these only Linux so far(Ubuntu14.04).
For example, on a fresh install of Ubuntu 14.04, install the following packages:

    $ sudo apt-get install build-essential git libgmp3-dev libprocps3-dev libgtest-dev python-markdown libboost-all-dev libssl-dev

Then, to compile, run:

    $ make
