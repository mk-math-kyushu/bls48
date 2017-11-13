/** @file
 *****************************************************************************
 * @author     This file is part of libsnark, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#include "algebra/curves/bls48/bls48_pp.hpp"
#include "common/profiling.hpp"

namespace libsnark {

void bls48_pp::init_public_params()
{
    init_bls48_params();
}

bls48_GT bls48_pp::final_exponentiation(const bls48_Fq48 &elt)
{
    return bls48_final_exponentiation(elt);
}

bls48_G1_precomp bls48_pp::precompute_G1(const bls48_G1 &P)
{
    return bls48_precompute_G1(P);
}

bls48_G2_precomp bls48_pp::precompute_G2(const bls48_G2 &Q)
{
    return bls48_precompute_G2(Q);
}

bls48_Fq48 bls48_pp::miller_loop(const bls48_G1_precomp &prec_P,
                                         const bls48_G2_precomp &prec_Q)
{
    return bls48_miller_loop(prec_P, prec_Q);
}


bls48_Fq48 bls48_pp::pairing(const bls48_G1 &P,
                                     const bls48_G2 &Q)
{
    return bls48_pairing(P, Q);

}

bls48_Fq48 bls48_pp::reduced_pairing(const bls48_G1 &P,
                                             const bls48_G2 &Q)
{
    return bls48_reduced_pairing(P, Q);
}

} // libsnark
