--------------------------------------------------------------------------------
Overview
--------------------------------------------------------------------------------
bls48 is an Optimal Ate Pairing(OAP) on BLS Curve implementation in C++.
This implementation is based on the drafts below;

- Y. Kiyomura, et al. "Secure and Eifficient Pairing at 256-Bit Secure Level" (2017)
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

In the case of Ubuntu 16.04, install the following packages:

    $ sudo apt-get install build-essential git libgmp3-dev libprocps4-dev libgtest-dev python-markdown libboost-all-dev libssl-dev

Then, to compile, run:

    $ make


--------------------------------------------------------------------------------
Timing Test
--------------------------------------------------------------------------------

* Ubuntu14.04　　
Intel(R) Core(TM) i5-4210M CPU @ 2.60GHz, 4 cores RAM 8GB

<table>
    <tr>
        <td>Security Level</td>
        <td>100</td>
        <td>256</td>
    </tr>
    <tr>
        <td>Curve Parameter</td>
        <td>BN256</td>
        <td>BLS48</td>
    </tr>
    <tr>
        <td>Scalar Multiplication on G1</td>
        <td>0.08[ms]</td>
        <td>2.32[ms]</td>
    </tr>
    <tr>
        <td>Scalar Multiplication on G2</td>
        <td>1.11[ms]</td>
        <td>160.09[ms]</td>
    </tr>
    <tr>
        <td>Pairing</td>
        <td>2.81[ms]</td>
        <td>509[ms]</td>
    </tr>
</table>



* Ubuntu16.04　　
Intel(R) Core(TM) i7 CPU 950 @ 3.07GHz, 8 cores RAM 24GB


<table>
    <tr>
        <td>Security Level</td>
        <td>100</td>
        <td>256</td>
    </tr>
    <tr>
        <td>Curve Parameter</td>
        <td>BN256</td>
        <td>BLS48</td>
    </tr>
    <tr>
        <td>Scalar Multiplication on G1</td>
        <td>0.06[ms]</td>
        <td>3.05[ms]</td>
    </tr>
    <tr>
        <td>Scalar Multiplication on G2</td>
        <td>1.24[ms]</td>
        <td>230.61[ms]</td>
    </tr>
    <tr>
        <td>Pairing</td>
        <td>4.05[ms]</td>
        <td>740.69[ms]</td>
    </tr>
</table>



* Raspberry Pi 3 model B
Linux devuan 4.6.3-gb48d47a(64bit OS) CPU Cortex-A53 ARMv8 64bit @ 1.2GHz, 4 cores RAM 1GB


<table>
    <tr>
        <td>Security Level</td>
        <td>100</td>
        <td>256</td>
    </tr>
    <tr>
        <td>Curve Parameter</td>
        <td>BN256</td>
        <td>BLS48</td>
    </tr>
    <tr>
        <td>Scalar Multiplication on G1</td>
        <td>6.68[ms]</td>
        <td>41.15[ms]</td>
    </tr>
    <tr>
        <td>Scalar Multiplication on G2</td>
        <td>24.23[ms]</td>
        <td>2643.23[ms]</td>
    </tr>
    <tr>
        <td>Pairing</td>
        <td>58[ms]</td>
        <td>8549.02[ms]</td>
    </tr>
</table>


* Raspberry Pi model B
Ubuntu16.04.2 LTS xenial(32bit OS) CPU Cortex-A53 ARMv8 64bit @ 1.2GHz, 4 cores RAM 1GB


<table>
    <tr>
        <td>Security Level</td>
        <td>100</td>
        <td>256</td>
    </tr>
    <tr>
        <td>Curve Parameter</td>
        <td>BN256</td>
        <td>BLS48</td>
    </tr>
    <tr>
        <td>Scalar Multiplication on G1</td>
        <td>5.39[ms]</td>
        <td>38.42[ms]</td>
    </tr>
    <tr>
        <td>Scalar Multiplication on G2</td>
        <td>21.81[ms]</td>
        <td>2660.55[ms]</td>
    </tr>
    <tr>
        <td>Pairing</td>
        <td>53.03[ms]</td>
        <td>8585.53[ms]</td>
    </tr>
</table>


* Mac
Mac OS Sierra Intel Core i7 @ 3.5 GHz, 2 cores RAM 16 GB

<table>
    <tr>
        <td>Security Level</td>
        <td>100</td>
        <td>256</td>
    </tr>
    <tr>
        <td>Curve Parameter</td>
        <td>BN256</td>
        <td>BLS48</td>
    </tr>
    <tr>
        <td>Scalar Multiplication on G1</td>
        <td>0[ms]</td>
        <td>1.07[ms]</td>
    </tr>
    <tr>
        <td>Scalar Multiplication on G2</td>
        <td>0.14[ms]</td>
        <td>112.86[ms]</td>
    </tr>
    <tr>
        <td>Pairing</td>
        <td>2.11[ms]</td>
        <td>373.02[ms]</td>
    </tr>
</table>


* Windows (Using MSYS2)
Windows10 Intel(R) Core(TM) i7-7500U CPU @ 2.70GHz, 4 cores RAM 16 GB

<table>
    <tr>
        <td>Security Level</td>
        <td>100</td>
        <td>256</td>
    </tr>
    <tr>
        <td>Curve Parameter</td>
        <td>BN256</td>
        <td>BLS48</td>
    </tr>
    <tr>
        <td>Scalar Multiplication on G1</td>
        <td>0[ms]</td>
        <td>1.5[ms]</td>
    </tr>
    <tr>
        <td>Scalar Multiplication on G2</td>
        <td>0.9[ms]</td>
        <td>127.74[ms]</td>
    </tr>
    <tr>
        <td>Pairing</td>
        <td>2.25[ms]</td>
        <td>409.8[ms]</td>
    </tr>
</table>
