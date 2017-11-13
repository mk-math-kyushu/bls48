/** @file
 *****************************************************************************
 Declaration of arithmetic in the finite field F[((p^2)^3)^2].
 *****************************************************************************
 * @author     This file is part of libsnark, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef FP48_2OVER2OVER2OVER3OVER2_HPP_
#define FP48_2OVER2OVER2OVER3OVER2_HPP_
#include "algebra/fields/fp.hpp"
#include "algebra/fields/fp2.hpp"
#include "algebra/fields/fp4.hpp"
#include "algebra/fields/fp8.hpp"
#include "algebra/fields/fp24_2over2over2over3.hpp"
#include <vector>

namespace libsnark {

template<mp_size_t n, const bigint<n>& modulus>
class Fp48_2over2over2over3over2_model;

template<mp_size_t n, const bigint<n>& modulus>
std::ostream& operator<<(std::ostream &, const Fp48_2over2over2over3over2_model<n, modulus> &);

template<mp_size_t n, const bigint<n>& modulus>
std::istream& operator>>(std::istream &, Fp48_2over2over2over3over2_model<n, modulus> &);

/**
 * Arithmetic in the finite field F[((((p^2)^2)^2)^3)^2].
 *
 * Let p := modulus. This interface provides arithmetic for the extension field
 * Fp48 = Fp24[S]/(S^2+Z)
 *
 * ASSUMPTION: p = 1 (mod 6)
 */
template<mp_size_t n, const bigint<n>& modulus>
class Fp48_2over2over2over3over2_model {
public:
    typedef Fp_model<n, modulus> my_Fp;
    typedef Fp2_model<n, modulus> my_Fp2;
    typedef Fp4_model<n, modulus> my_Fp4;
    typedef Fp8_model<n, modulus> my_Fp8;
    typedef Fp24_2over2over2over3_model<n, modulus> my_Fp24_2over2over2over3;

    static Fp24_2over2over2over3_model<n, modulus> non_residue;
    static Fp24_2over2over2over3_model<n, modulus> Frobenius_coeffs_c1[48]; // non_residue^((modulus^i-1)/48) for i=0,...,47

    my_Fp24_2over2over2over3 c0, c1;
    Fp48_2over2over2over3over2_model() {};
    Fp48_2over2over2over3over2_model(const my_Fp24_2over2over2over3& c0, const my_Fp24_2over2over2over3& c1) : c0(c0), c1(c1) {};

    void clear() { c0.clear(); c1.clear(); }
    void print() const { printf("c0/c1:\n"); c0.print(); c1.print(); }

    static Fp48_2over2over2over3over2_model<n, modulus> zero();
    static Fp48_2over2over2over3over2_model<n, modulus> one();
    static Fp48_2over2over2over3over2_model<n, modulus> random_element();

    bool is_zero() const { return c0.is_zero() && c1.is_zero(); }
    bool operator==(const Fp48_2over2over2over3over2_model &other) const;
    bool operator!=(const Fp48_2over2over2over3over2_model &other) const;

    Fp48_2over2over2over3over2_model operator+(const Fp48_2over2over2over3over2_model &other) const;
    Fp48_2over2over2over3over2_model operator-(const Fp48_2over2over2over3over2_model &other) const;
    Fp48_2over2over2over3over2_model operator*(const Fp48_2over2over2over3over2_model &other) const;
    Fp48_2over2over2over3over2_model operator-() const;
    Fp48_2over2over2over3over2_model squared() const; // default is squared_complex
    Fp48_2over2over2over3over2_model squared_karatsuba() const;
    Fp48_2over2over2over3over2_model squared_complex() const;
    Fp48_2over2over2over3over2_model inverse() const;
    Fp48_2over2over2over3over2_model Frobenius_map(unsigned long power) const;
    Fp48_2over2over2over3over2_model unitary_inverse() const;
    //Fp48_2over2over2over3over2_model cyclotomic_squared() const;

    static my_Fp24_2over2over2over3 mul_by_non_residue(const my_Fp24_2over2over2over3 &elt);

    //template<mp_size_t m>
    //Fp48_2over2over2over3over2_model cyclotomic_exp(const bigint<m> &exponent) const;

    static bigint<n> base_field_char() { return modulus; }
    static size_t extension_degree() { return 48; }

    friend std::ostream& operator<< <n, modulus>(std::ostream &out, const Fp48_2over2over2over3over2_model<n, modulus> &el);
    friend std::istream& operator>> <n, modulus>(std::istream &in, Fp48_2over2over2over3over2_model<n, modulus> &el);
};

template<mp_size_t n, const bigint<n>& modulus>
std::ostream& operator<<(std::ostream& out, const std::vector<Fp48_2over2over2over3over2_model<n, modulus> > &v);

template<mp_size_t n, const bigint<n>& modulus>
std::istream& operator>>(std::istream& in, std::vector<Fp48_2over2over2over3over2_model<n, modulus> > &v);

template<mp_size_t n, const bigint<n>& modulus>
Fp48_2over2over2over3over2_model<n, modulus> operator*(const Fp_model<n, modulus> &lhs, const Fp48_2over2over2over3over2_model<n, modulus> &rhs);

template<mp_size_t n, const bigint<n>& modulus>
Fp48_2over2over2over3over2_model<n, modulus> operator*(const Fp2_model<n, modulus> &lhs, const Fp48_2over2over2over3over2_model<n, modulus> &rhs);

template<mp_size_t n, const bigint<n>& modulus>
Fp48_2over2over2over3over2_model<n, modulus> operator*(const Fp4_model<n, modulus> &lhs, const Fp48_2over2over2over3over2_model<n, modulus> &rhs);

template<mp_size_t n, const bigint<n>& modulus>
Fp48_2over2over2over3over2_model<n, modulus> operator*(const Fp8_model<n, modulus> &lhs, const Fp48_2over2over2over3over2_model<n, modulus> &rhs);

template<mp_size_t n, const bigint<n>& modulus>
Fp48_2over2over2over3over2_model<n, modulus> operator*(const Fp24_2over2over2over3_model<n, modulus> &lhs, const Fp48_2over2over2over3over2_model<n, modulus> &rhs);

template<mp_size_t n, const bigint<n>& modulus, mp_size_t m>
Fp48_2over2over2over3over2_model<n, modulus> operator^(const Fp48_2over2over2over3over2_model<n, modulus> &self, const bigint<m> &exponent);

template<mp_size_t n, const bigint<n>& modulus, mp_size_t m, const bigint<m>& exp_modulus>
Fp48_2over2over2over3over2_model<n, modulus> operator^(const Fp48_2over2over2over3over2_model<n, modulus> &self, const Fp_model<m, exp_modulus> &exponent);

template<mp_size_t n, const bigint<n>& modulus>
Fp24_2over2over2over3_model<n, modulus> Fp48_2over2over2over3over2_model<n, modulus>::non_residue;

template<mp_size_t n, const bigint<n>& modulus>
Fp24_2over2over2over3_model<n, modulus> Fp48_2over2over2over3over2_model<n, modulus>::Frobenius_coeffs_c1[48];

} // libsnark
#include "algebra/fields/fp48_2over2over2over3over2.tcc"
#endif // FP48_2OVER2OVER2OVER3OVER2_HPP_
