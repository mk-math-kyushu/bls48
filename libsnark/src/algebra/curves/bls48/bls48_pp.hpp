/** @file
 *****************************************************************************
 * @author     This file is part of libsnark, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef BLS48_PP_HPP_
#define BLS48_PP_HPP_
#include "algebra/curves/public_params.hpp"
#include "algebra/curves/bls48/bls48_init.hpp"
#include "algebra/curves/bls48/bls48_g1.hpp"
#include "algebra/curves/bls48/bls48_g2.hpp"
#include "algebra/curves/bls48/bls48_pairing.hpp"


namespace libsnark {

class bls48_pp {
public:
    typedef bls48_Fr Fp_type;
    typedef bls48_G1 G1_type;
    typedef bls48_G2 G2_type;
    typedef bls48_G1_precomp G1_precomp_type;
    typedef bls48_G2_precomp G2_precomp_type;
    typedef bls48_Fq Fq_type;
    typedef bls48_Fq8 Fqe_type;
    typedef bls48_Fq48 Fqk_type;
    typedef bls48_Fq24 Fqt_type;
    typedef bls48_GT GT_type;

    static const bool has_affine_pairing = false;


    static void init_public_params();
    static bls48_GT final_exponentiation(const bls48_Fq48 &elt);
    static bls48_G1_precomp precompute_G1(const bls48_G1 &P);
    static bls48_G2_precomp precompute_G2(const bls48_G2 &Q);
    static bls48_Fq48 miller_loop(const bls48_G1_precomp &prec_P,
                                      const bls48_G2_precomp &prec_Q);
    static bls48_Fq48 pairing(const bls48_G1 &P,
                                  const bls48_G2 &Q);
    static bls48_Fq48 reduced_pairing(const bls48_G1 &P,
                                          const bls48_G2 &Q);

};

} // libsnark
#endif // BLS48_PP_HPP_
