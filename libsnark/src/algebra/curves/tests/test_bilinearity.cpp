/**
 *****************************************************************************
 * @author     This file is part of libsnark, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/
#include "common/profiling.hpp"

#include "algebra/curves/alt_bn128/alt_bn128_pp.hpp"
#include "algebra/curves/bls48/bls48_pp.hpp"
#include "iostream"

using namespace libsnark;

template<typename ppT>
void pairing_test()
{

    GT<ppT> GT_one = GT<ppT>::one();

    printf("Running bilinearity tests:\n");
    G1<ppT> P = (Fr<ppT>::random_element()) * G1<ppT>::one();
    //G1<ppT> P = Fr<ppT>("2") * G1<ppT>::one();
    G2<ppT> Q = (Fr<ppT>::random_element()) * G2<ppT>::one();//G2<ppT> Q = G2<ppT>::one();//
    //std::cout << "Q is in G2"<< Q.is_well_formed() << std::endl;
    //assert(Q.Y.squared() - Q.X.squared()*Q.X == (Q.Z.squared()*Q.Z).squared());

    //G2<ppT> Q = Fr<ppT>("3") * G2<ppT>::one();

    //printf("GT_one:\n");
    //GT_one.print();


    printf("P:\n");
    P.print();
    P.print_coordinates();
    printf("Q:\n");
    Q.print();
    Q.print_coordinates();
    printf("\n\n");




    Fr<ppT> s = Fr<ppT>::random_element();
    //Fr<ppT> s = Fr<ppT>("2");
    G1<ppT> sP = s * P;//G1<ppT>::G1_zero;//
    G2<ppT> sQ = s * Q;//G2<ppT>::G2_one;//
    //assert(Q + Q == sQ);
    //std::cout << "sQ is in G2  = "<< sQ.is_well_formed() << std::endl;


    printf("Pairing bilinearity tests (three must match):\n");
    //printf("e(sP, Q):\n");
    //GT<ppT> pre_ans1 = ppT::pairing(sP, Q);
    //pre_ans1.print();

    GT<ppT> ans1 = ppT::reduced_pairing(sP, Q);


    //printf("e(P, sQ):\n");
    //GT<ppT> pre_ans2 = ppT::pairing(P, sQ);
    //pre_ans2.print();
    GT<ppT> ans2 = ppT::reduced_pairing(P, sQ);
    //printf("e(P, Q)^s:\n");
    GT<ppT> ans3 = ppT::reduced_pairing(P, Q)^s;

    //std::cout << "Fr<ppT>::field_char()" << Fr<ppT>::field_char() << std::endl;

    assert(ans1 != GT_one);
    assert((ans1^Fr<ppT>::field_char()) == GT_one);


    ans1.print();
    ans2.print();
    ans3.print();
    assert(ans1 == ans3);
    assert(ans2 == ans3);




    printf("\n\n");


}


template<typename ppT>
void double_miller_loop_test()
{
    const G1<ppT> P1 = (Fr<ppT>::random_element()) * G1<ppT>::one();
    const G1<ppT> P2 = (Fr<ppT>::random_element()) * G1<ppT>::one();
    const G2<ppT> Q1 = (Fr<ppT>::random_element()) * G2<ppT>::one();
    const G2<ppT> Q2 = (Fr<ppT>::random_element()) * G2<ppT>::one();

    const G1_precomp<ppT> prec_P1 = ppT::precompute_G1(P1);
    const G1_precomp<ppT> prec_P2 = ppT::precompute_G1(P2);
    const G2_precomp<ppT> prec_Q1 = ppT::precompute_G2(Q1);
    const G2_precomp<ppT> prec_Q2 = ppT::precompute_G2(Q2);

    const Fqk<ppT> ans_1 = ppT::miller_loop(prec_P1, prec_Q1);
    const Fqk<ppT> ans_2 = ppT::miller_loop(prec_P2, prec_Q2);
    const Fqk<ppT> ans_12 = ppT::double_miller_loop(prec_P1, prec_Q1, prec_P2, prec_Q2);
    assert(ans_1 * ans_2 == ans_12);
}


template<typename ppT>
void affine_pairing_test()
{
    GT<ppT> GT_one = GT<ppT>::one();

    printf("Running bilinearity tests:\n");
    G1<ppT> P = (Fr<ppT>::random_element()) * G1<ppT>::one();
    G2<ppT> Q = (Fr<ppT>::random_element()) * G2<ppT>::one();

/*
    printf("P:\n");
    P.print();
    printf("Q:\n");
    Q.print();
    printf("\n\n");
*/
    Fr<ppT> s = Fr<ppT>::random_element();
    G1<ppT> sP = s * P;
    G2<ppT> sQ = s * Q;

    printf("Pairing bilinearity tests (three must match):\n");
    GT<ppT> ans1 = ppT::affine_reduced_pairing(sP, Q);
    GT<ppT> ans2 = ppT::affine_reduced_pairing(P, sQ);
    GT<ppT> ans3 = ppT::affine_reduced_pairing(P, Q)^s;
    ans1.print();
    ans2.print();
    //ans3.print();
    assert(ans1 == ans2);
    //assert(ans2 == ans3);

    //assert(ans1 != GT_one);
    //assert((ans1^Fr<ppT>::field_char()) == GT_one);
    printf("\n\n");
}

int main(void)
{
    start_profiling();
/*
    edwards_pp::init_public_params();
    pairing_test<edwards_pp>();
    double_miller_loop_test<edwards_pp>();

    mnt6_pp::init_public_params();
    pairing_test<mnt6_pp>();
    double_miller_loop_test<mnt6_pp>();
    affine_pairing_test<mnt6_pp>();

    mnt4_pp::init_public_params();
    pairing_test<mnt4_pp>();
    double_miller_loop_test<mnt4_pp>();
    affine_pairing_test<mnt4_pp>();





    alt_bn128_pp::init_public_params();
    pairing_test<alt_bn128_pp>();
    std::cout << "alt" << std::endl;
    //double_miller_loop_test<alt_bn128_pp>();
*/

    std::cout << "bls" << std::endl;
    bls48_pp::init_public_params();
    pairing_test<bls48_pp>();


/*
#ifdef CURVE_BN128       // BN128 has fancy dependencies so it may be disabled
    bn128_pp::init_public_params();
    pairing_test<bn128_pp>();
    double_miller_loop_test<bn128_pp>();
#endif
*/
}
