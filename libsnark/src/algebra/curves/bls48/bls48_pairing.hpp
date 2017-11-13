/** @file
 *****************************************************************************
 * @author     This file is part of libsnark, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef BLS48_PAIRING_HPP_
#define BLS48_PAIRING_HPP_
#include <vector>
#include "algebra/curves/bls48/bls48_init.hpp"

namespace libsnark {

/* final exponentiation */

bls48_GT bls48_final_exponentiation(const bls48_Fq48 &elt);

/* ate pairing */



struct bls48_ate_G1_precomp {
    bls48_Fq PX;
    bls48_Fq PY;

    bool operator==(const bls48_ate_G1_precomp &other) const;
    friend std::ostream& operator<<(std::ostream &out, const bls48_ate_G1_precomp &prec_P);
    friend std::istream& operator>>(std::istream &in, bls48_ate_G1_precomp &prec_P);
};

/*
struct bls48_ate_ell_coeffs {
    bls48_Fq8 minus_ramda;
    bls48_Fq8 c;

    bool operator==(const bls48_ate_ell_coeffs &other) const;
    friend std::ostream& operator<<(std::ostream &out, const bls48_ate_ell_coeffs &dc);
    friend std::istream& operator>>(std::istream &in, bls48_ate_ell_coeffs &dc);
};
*/
struct bls48_ate_G2_precomp {
    bls48_Fq8 QX;
    bls48_Fq8 QY;
    //std::vector<bls48_ate_ell_coeffs> coeffs;

    bool operator==(const bls48_ate_G2_precomp &other) const;
    friend std::ostream& operator<<(std::ostream &out, const bls48_ate_G2_precomp &prec_Q);
    friend std::istream& operator>>(std::istream &in, bls48_ate_G2_precomp &prec_Q);
};

bls48_ate_G1_precomp bls48_ate_precompute_G1(const bls48_G1& P);
bls48_ate_G2_precomp bls48_ate_precompute_G2(const bls48_G2& Q);

bls48_Fq48 bls48_ate_miller_loop(const bls48_ate_G1_precomp &prec_P,
                              const bls48_ate_G2_precomp &prec_Q);
bls48_Fq48 bls48_ate_double_miller_loop(const bls48_ate_G1_precomp &prec_P1,
                                     const bls48_ate_G2_precomp &prec_Q1,
                                     const bls48_ate_G1_precomp &prec_P2,
                                     const bls48_ate_G2_precomp &prec_Q2);

bls48_Fq48 bls48_ate_pairing(const bls48_G1& P,
                          const bls48_G2 &Q);
bls48_GT bls48_ate_reduced_pairing(const bls48_G1 &P,
                                 const bls48_G2 &Q);

/* choice of pairing */

typedef bls48_ate_G1_precomp bls48_G1_precomp;
typedef bls48_ate_G2_precomp bls48_G2_precomp;

bls48_G1_precomp bls48_precompute_G1(const bls48_G1& P);

bls48_G2_precomp bls48_precompute_G2(const bls48_G2& Q);

bls48_Fq48 bls48_miller_loop(const bls48_G1_precomp &prec_P,
                          const bls48_G2_precomp &prec_Q);
/*
bls48_Fq48 bls48_double_miller_loop(const bls48_G1_precomp &prec_P1,
                                 const bls48_G2_precomp &prec_Q1,
                                 const bls48_G1_precomp &prec_P2,
                                 const bls48_G2_precomp &prec_Q2);
*/
bls48_Fq48 bls48_pairing(const bls48_G1& P,
                      const bls48_G2 &Q);

bls48_GT bls48_reduced_pairing(const bls48_G1 &P,
                             const bls48_G2 &Q);

//bls48_GT bls48_affine_reduced_pairing(const bls48_G1 &P,
//                                    const bls48_G2 &Q);

} // libsnark
#endif // BLS48_PAIRING_HPP_
