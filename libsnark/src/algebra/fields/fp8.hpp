/** @file
 *****************************************************************************

////////////
 Declaration of interfaces for the (extension) field Fp8.

 The field Fp8 equals Fp4[W/(W^2-V) where Fp4 = Fp2[V]/(V^2-U).

 ASSUMPTION: the modulus p is 1 mod 6.
////////////
///
////////////

 *****************************************************************************
 * @author     This file is part of libsnark, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/


#ifndef FP8_HPP_
#define FP8_HPP_
////////////
///
////////////

#include "algebra/fields/fp.hpp"
#include "algebra/fields/fp2.hpp"
#include "algebra/fields/fp4.hpp"//


namespace libsnark {

template<mp_size_t n, const bigint<n>& modulus>
class Fp8_model;//

template<mp_size_t n, const bigint<n>& modulus>
std::ostream& operator<<(std::ostream &, const Fp8_model<n, modulus> &);//

template<mp_size_t n, const bigint<n>& modulus>
std::istream& operator>>(std::istream &, Fp8_model<n, modulus> &);

template<mp_size_t n, const bigint<n>& modulus>
class Fp8_model {//
public:
    typedef Fp_model<n, modulus> my_Fp;
    typedef Fp2_model<n, modulus> my_Fp2;
    typedef Fp4_model<n, modulus> my_Fp4;//add
    typedef my_Fp4 my_Fpe;

    static bigint<8*n> euler; // (modulus^2-1)/2
    static size_t s;       // modulus^2 = 2^s * t + 1
    static bigint<8*n> t;  // with t odd
    static bigint<8*n> t_minus_1_over_2; // (t-1)/2
    static Fp8_model<n, modulus> nqr; // a quadratic nonresidue in Fp2
    static Fp8_model<n, modulus> nqr_to_t; // nqr^t

    static my_Fp4 non_residue;//
    static my_Fp4 Frobenius_coeffs_c1[8]; //non_residue^((modulus^i-1)/4) for i=0,1,2,3

    my_Fp4 c0, c1;//
    Fp8_model() {};
    Fp8_model(const my_Fp4& c0, const my_Fp4& c1) : c0(c0), c1(c1) {};

    void print() const { printf("c0/c1:\n"); c0.print(); c1.print(); }
    void clear() { c0.clear(); c1.clear(); }

    static Fp8_model<n, modulus> zero();//
    static Fp8_model<n, modulus> one();//
    static Fp8_model<n, modulus> random_element();//

    bool is_zero() const { return c0.is_zero() && c1.is_zero(); }
    bool operator==(const Fp8_model &other) const;//
    bool operator!=(const Fp8_model &other) const;//

    Fp8_model operator+(const Fp8_model &other) const;//
    Fp8_model operator-(const Fp8_model &other) const;//
    Fp8_model operator*(const Fp8_model &other) const;//
    Fp8_model mul_by_023(const Fp8_model &other) const;
    Fp8_model operator-() const;//
    Fp8_model squared() const;//

    Fp8_model inverse() const;//
    Fp8_model Frobenius_map(unsigned long power) const;
    Fp8_model unitary_inverse() const;
    Fp8_model sqrt() const; // HAS TO BE A SQUARE (else does not terminate)
    //Fp4_model cyclotomic_squared() const;

    static my_Fp4 mul_by_non_residue(const my_Fp4 &elt);//

    //template<mp_size_t m>
    //Fp4_model cyclotomic_exp(const bigint<m> &exponent) const;

    static size_t size_in_bits() { return 8*my_Fp::size_in_bits(); }
    static bigint<n> base_field_char() { return modulus; }
    static constexpr size_t extension_degree() { return 8; }

    friend std::ostream& operator<< <n, modulus>(std::ostream &out, const Fp8_model<n, modulus> &el);
    friend std::istream& operator>> <n, modulus>(std::istream &in, Fp8_model<n, modulus> &el);
};

template<mp_size_t n, const bigint<n>& modulus>
Fp8_model<n, modulus> operator*(const Fp_model<n, modulus> &lhs, const Fp8_model<n, modulus> &rhs);

template<mp_size_t n, const bigint<n>& modulus>
Fp8_model<n, modulus> operator*(const Fp2_model<n, modulus> &lhs, const Fp8_model<n, modulus> &rhs);

////////
//add the case which lhs is in Fp4
////////
template<mp_size_t n, const bigint<n>& modulus>
Fp8_model<n, modulus> operator*(const Fp4_model<n, modulus> &lhs, const Fp8_model<n, modulus> &rhs);

template<mp_size_t n, const bigint<n>& modulus, mp_size_t m>
Fp8_model<n, modulus> operator^(const Fp8_model<n, modulus> &self, const bigint<m> &exponent);

template<mp_size_t n, const bigint<n>& modulus, mp_size_t m, const bigint<m>& modulus_p>
Fp8_model<n, modulus> operator^(const Fp8_model<n, modulus> &self, const Fp_model<m, modulus_p> &exponent);

template<mp_size_t n, const bigint<n>& modulus>
Fp4_model<n, modulus> Fp8_model<n, modulus>::non_residue;

template<mp_size_t n, const bigint<n>& modulus>
Fp4_model<n, modulus> Fp8_model<n, modulus>::Frobenius_coeffs_c1[8];




template<mp_size_t n, const bigint<n>& modulus>
bigint<8*n> Fp8_model<n, modulus>::euler;

template<mp_size_t n, const bigint<n>& modulus>
size_t Fp8_model<n, modulus>::s;

template<mp_size_t n, const bigint<n>& modulus>
bigint<8*n> Fp8_model<n, modulus>::t;

template<mp_size_t n, const bigint<n>& modulus>
bigint<8*n> Fp8_model<n, modulus>::t_minus_1_over_2;

template<mp_size_t n, const bigint<n>& modulus>
Fp8_model<n, modulus> Fp8_model<n, modulus>::nqr;

template<mp_size_t n, const bigint<n>& modulus>
Fp8_model<n, modulus> Fp8_model<n, modulus>::nqr_to_t;


////

} // libsnark

#include "algebra/fields/fp8.tcc"

#endif // FP8_HPP_
