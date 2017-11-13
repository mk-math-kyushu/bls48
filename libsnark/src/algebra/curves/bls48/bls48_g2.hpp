/** @file
 *****************************************************************************
 * @author     This file is part of libsnark, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef BLS48_G2_HPP_
#define BLS48_G2_HPP_
#include <vector>
#include "algebra/curves/bls48/bls48_init.hpp"
#include "algebra/curves/curve_utils.hpp"

namespace libsnark {

class bls48_G2;
std::ostream& operator<<(std::ostream &, const bls48_G2&);
std::istream& operator>>(std::istream &, bls48_G2&);

class bls48_G2 {
public:
#ifdef PROFILE_OP_COUNTS
    static long long add_cnt;
    static long long dbl_cnt;
#endif
    static std::vector<size_t> wnaf_window_table;
    static std::vector<size_t> fixed_base_exp_window_table;
    static bls48_G2 G2_zero;
    static bls48_G2 G2_one;

    typedef bls48_Fq base_field;
    typedef bls48_Fq8 twist_field;
    typedef bls48_Fr scalar_field;

    bls48_Fq8 X, Y, Z;

    // using Jacobian coordinates
    bls48_G2();
    bls48_G2(const bls48_Fq8& X, const bls48_Fq8& Y, const bls48_Fq8& Z) : X(X), Y(Y), Z(Z) {};

    //static bls48_Fq8 mul_by_b(const bls48_Fq8 &elt);

    void print() const;
    void print_coordinates() const;

    void to_affine_coordinates();
    void to_special();
    bool is_special() const;

    bool is_zero() const;

    bool operator==(const bls48_G2 &other) const;
    bool operator!=(const bls48_G2 &other) const;

    bls48_G2 operator+(const bls48_G2 &other) const;
    bls48_G2 operator-() const;
    bls48_G2 operator-(const bls48_G2 &other) const;

    bls48_G2 add(const bls48_G2 &other) const;
    bls48_G2 mixed_add(const bls48_G2 &other) const;
    bls48_G2 dbl() const;
    //bls48_G2 mul_by_q() const;

    bool is_well_formed() const;

    static bls48_G2 zero();
    static bls48_G2 one();
    static bls48_G2 random_element();

    static size_t size_in_bits() { return twist_field::size_in_bits() + 1; }
    static bigint<base_field::num_limbs> base_field_char() { return base_field::field_char(); }
    static bigint<scalar_field::num_limbs> order() { return scalar_field::field_char(); }

    friend std::ostream& operator<<(std::ostream &out, const bls48_G2 &g);
    friend std::istream& operator>>(std::istream &in, bls48_G2 &g);
};

template<mp_size_t m>
bls48_G2 operator*(const bigint<m> &lhs, const bls48_G2 &rhs)
{
    return scalar_mul<bls48_G2, m>(rhs, lhs);
}

template<mp_size_t m, const bigint<m>& modulus_p>
bls48_G2 operator*(const Fp_model<m,modulus_p> &lhs, const bls48_G2 &rhs)
{
    return scalar_mul<bls48_G2, m>(rhs, lhs.as_bigint());
}

template<typename T>
void batch_to_special_all_non_zeros(std::vector<T> &vec);
template<>
void batch_to_special_all_non_zeros<bls48_G2>(std::vector<bls48_G2> &vec);

} // libsnark
#endif // BLS48_G2_HPP_
