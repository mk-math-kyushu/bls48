/** @file
 *****************************************************************************
 Declaration of arithmetic in the finite field F[(p^2)^3]
 *****************************************************************************
 * @author     This file is part of libsnark, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef FP24_2OVER2OVER2OVER3_HPP_
#define FP24_2OVER2OVER2OVER3_HPP_
#include "algebra/fields/fp.hpp"
#include "algebra/fields/fp2.hpp"
#include "algebra/fields/fp4.hpp"
#include "algebra/fields/fp8.hpp"
#include <vector>

namespace libsnark {

template<mp_size_t n, const bigint<n>& modulus>
class Fp24_2over2over2over3_model;

template<mp_size_t n, const bigint<n>& modulus>
std::ostream& operator<<(std::ostream &, const Fp24_2over2over2over3_model<n, modulus> &);

template<mp_size_t n, const bigint<n>& modulus>
std::istream& operator>>(std::istream &, Fp24_2over2over2over3_model<n, modulus> &);

/**
 * Arithmetic in the finite field F[(((p^2)^2)^2)^3].
 *
 * Let p := modulus. This interface provides arithmetic for the extension field
 *  Fp24 = Fp8[V]/(Z^3-non_residue) where non_residue is in Fp4.
 *
 * ASSUMPTION: p = 1 (mod 6)
 */
template<mp_size_t n, const bigint<n>& modulus>
class Fp24_2over2over2over3_model {
public:
    typedef Fp_model<n, modulus> my_Fp;
    typedef Fp2_model<n, modulus> my_Fp2;
    typedef Fp4_model<n, modulus> my_Fp4;
    typedef Fp8_model<n, modulus> my_Fp8;

    static my_Fp8 non_residue;
    static my_Fp8 Frobenius_coeffs_c1[24]; // non_residue^((modulus^i-1)/3)   for i=0,1,2,..,23
    static my_Fp8 Frobenius_coeffs_c2[24]; // non_residue^((2*modulus^i-2)/3) for i=0,1,2,..,23

    my_Fp8 c0, c1, c2;
    Fp24_2over2over2over3_model() {};
    Fp24_2over2over2over3_model(const my_Fp8& c0, const my_Fp8& c1, const my_Fp8& c2) : c0(c0), c1(c1), c2(c2) {};

    void clear() { c0.clear(); c1.clear(); c2.clear(); }
    void print() const { printf("c0/c1/c2:\n"); c0.print(); c1.print(); c2.print(); }

    static Fp24_2over2over2over3_model<n, modulus> zero();
    static Fp24_2over2over2over3_model<n, modulus> one();
    static Fp24_2over2over2over3_model<n, modulus> random_element();

    bool is_zero() const { return c0.is_zero() && c1.is_zero() && c2.is_zero(); }
    bool operator==(const Fp24_2over2over2over3_model &other) const;
    bool operator!=(const Fp24_2over2over2over3_model &other) const;

    Fp24_2over2over2over3_model operator+(const Fp24_2over2over2over3_model &other) const;
    Fp24_2over2over2over3_model operator-(const Fp24_2over2over2over3_model &other) const;
    Fp24_2over2over2over3_model operator*(const Fp24_2over2over2over3_model &other) const;
    Fp24_2over2over2over3_model operator-() const;
    Fp24_2over2over2over3_model squared() const;
    Fp24_2over2over2over3_model inverse() const;
    Fp24_2over2over2over3_model Frobenius_map(unsigned long power) const;

    static my_Fp8 mul_by_non_residue(const my_Fp8 &elt);

    template<mp_size_t m>
    Fp24_2over2over2over3_model operator^(const bigint<m> &other) const;

    static bigint<n> base_field_char() { return modulus; }
    static size_t extension_degree() { return 24; }

    friend std::ostream& operator<< <n, modulus>(std::ostream &out, const Fp24_2over2over2over3_model<n, modulus> &el);
    friend std::istream& operator>> <n, modulus>(std::istream &in, Fp24_2over2over2over3_model<n, modulus> &el);
};

template<mp_size_t n, const bigint<n>& modulus>
std::ostream& operator<<(std::ostream& out, const std::vector<Fp24_2over2over2over3_model<n, modulus> > &v);

template<mp_size_t n, const bigint<n>& modulus>
std::istream& operator>>(std::istream& in, std::vector<Fp24_2over2over2over3_model<n, modulus> > &v);

template<mp_size_t n, const bigint<n>& modulus>
Fp24_2over2over2over3_model<n, modulus> operator*(const Fp_model<n, modulus> &lhs, const Fp24_2over2over2over3_model<n, modulus> &rhs);

template<mp_size_t n, const bigint<n>& modulus>
Fp24_2over2over2over3_model<n, modulus> operator*(const Fp2_model<n, modulus> &lhs, const Fp24_2over2over2over3_model<n, modulus> &rhs);

template<mp_size_t n, const bigint<n>& modulus>
Fp24_2over2over2over3_model<n, modulus> operator*(const Fp4_model<n, modulus> &lhs, const Fp24_2over2over2over3_model<n, modulus> &rhs);

template<mp_size_t n, const bigint<n>& modulus>
Fp24_2over2over2over3_model<n, modulus> operator*(const Fp8_model<n, modulus> &lhs, const Fp24_2over2over2over3_model<n, modulus> &rhs);

template<mp_size_t n, const bigint<n>& modulus>
Fp8_model<n, modulus> Fp24_2over2over2over3_model<n, modulus>::non_residue;

template<mp_size_t n, const bigint<n>& modulus>
Fp8_model<n, modulus> Fp24_2over2over2over3_model<n, modulus>::Frobenius_coeffs_c1[24];

template<mp_size_t n, const bigint<n>& modulus>
Fp8_model<n, modulus> Fp24_2over2over2over3_model<n, modulus>::Frobenius_coeffs_c2[24];

} // libsnark
#include "algebra/fields/fp24_2over2over2over3.tcc"

#endif // FP24_2OVER2OVER2OVER3_HPP_
