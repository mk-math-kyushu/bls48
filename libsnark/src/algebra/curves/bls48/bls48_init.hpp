/** @file
 *****************************************************************************
 * @author     This file is part of libsnark, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef BLS48_INIT_HPP_
#define BLS48_INIT_HPP_
#include "algebra/curves/public_params.hpp"
#include "algebra/fields/fp.hpp"
#include "algebra/fields/fp2.hpp"
#include "algebra/fields/fp4.hpp"
#include "algebra/fields/fp8.hpp"
#include "algebra/fields/fp24_2over2over2over3.hpp"
#include "algebra/fields/fp48_2over2over2over3over2.hpp"


namespace libsnark {

const mp_size_t bls48_r_bitcount = 518;//518
const mp_size_t bls48_q_bitcount = 581;//581

const mp_size_t bls48_r_limbs = (bls48_r_bitcount+GMP_NUMB_BITS-1)/GMP_NUMB_BITS;
const mp_size_t bls48_q_limbs = (bls48_q_bitcount+GMP_NUMB_BITS-1)/GMP_NUMB_BITS;

extern bigint<bls48_r_limbs> bls48_modulus_r;
extern bigint<bls48_q_limbs> bls48_modulus_q;

typedef Fp_model<bls48_r_limbs, bls48_modulus_r> bls48_Fr;
typedef Fp_model<bls48_q_limbs, bls48_modulus_q> bls48_Fq;
typedef Fp2_model<bls48_q_limbs, bls48_modulus_q> bls48_Fq2;
typedef Fp4_model<bls48_q_limbs, bls48_modulus_q> bls48_Fq4;
typedef Fp8_model<bls48_q_limbs, bls48_modulus_q> bls48_Fq8;
typedef Fp24_2over2over2over3_model<bls48_q_limbs, bls48_modulus_q> bls48_Fq24;
typedef Fp48_2over2over2over3over2_model<bls48_q_limbs, bls48_modulus_q> bls48_Fq48;
typedef bls48_Fq48 bls48_GT;



// parameters for BLS48 E/Fq : y^2 = x^3 + b
extern bls48_Fq bls48_coeff_b;
// parameters for twisted curve for BLS48 E'/Fq8 : y^2 = x^3 + b/omega
extern bls48_Fq8 bls48_twist;
extern bls48_Fq8 bls48_twist_coeff_b;
//extern bls48_Fq4 bls48_twist_mul_by_b_c0;
//extern bls48_Fq4 bls48_twist_mul_by_b_c1;
//extern bls48_Fq8 bls48_twist_mul_by_q_X;
//extern bls48_Fq8 bls48_twist_mul_by_q_Y;



// parameters for pairing
extern bigint<bls48_q_limbs> bls48_ate_loop_count;
extern bool bls48_ate_is_loop_count_neg;
extern bigint<48*bls48_q_limbs> bls48_final_exponent;
extern bigint<bls48_q_limbs> bls48_final_exponent_z;
extern bool bls48_final_exponent_is_z_neg;



void init_bls48_params();

class bls48_G1;
class bls48_G2;

} // libsnark
#endif // BLS48_INIT_HPP_
