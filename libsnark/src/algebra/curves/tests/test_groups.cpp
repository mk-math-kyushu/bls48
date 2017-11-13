/**
 *****************************************************************************
 * @author     This file is part of libsnark, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/
#include "common/profiling.hpp"
#include "algebra/curves/alt_bn128/alt_bn128_pp.hpp"
#include "algebra/curves/bls48/bls48_pp.hpp"
#include <sstream>
#include <iostream>

using namespace libsnark;

template<typename GroupT>
void test_mixed_add()
{
    GroupT base, el, result;

    base = GroupT::zero();
    el = GroupT::zero();
    el.to_special();
    result = base.mixed_add(el);
    assert(result == base + el);

    //std::cout << "y^2 - x^3 = " << (GroupT::one().Y.squared() - GroupT::one().X.squared()*GroupT::one().X) << std::endl;

    //std::cout << "Fq*::one = " << GroupT::one().Z << std::endl;

    base = GroupT::zero();
    el = GroupT::random_element();
    el.to_special();
    result = base.mixed_add(el);
    assert(result == base + el);

    base = GroupT::random_element();
    el = GroupT::zero();
    el.to_special();
    result = base.mixed_add(el);
    assert(result == base + el);

    base = GroupT::random_element();
    el = GroupT::random_element();
    el.to_special();
    result = base.mixed_add(el);
    assert(result == base + el);

    base = GroupT::random_element();
    el = base;
    el.to_special();
    result = base.mixed_add(el);
    assert(result == base.dbl());
}

template<typename GroupT>
void test_group()
{
    bigint<1> rand1 = bigint<1>("76749407");
    bigint<1> rand2 = bigint<1>("44410867");
    bigint<1> randsum = bigint<1>("121160274");

    GroupT zero = GroupT::zero();
    assert(zero == zero);
    //std::cout << "GroupT::zero " << zero << std::endl;
    //zero.print_coordinates();
    GroupT one = GroupT::one();
    assert(one == one);
    //std::cout << "GroupT::one " << one << std::endl;
    //one.print_coordinates();
    GroupT two = bigint<1>(2l) * GroupT::one();
    assert(two == two);
    //std::cout << "GroupT::two " << two << std::endl;
    GroupT five = bigint<1>(5l) * GroupT::one();

    GroupT three = bigint<1>(3l) * GroupT::one();
    GroupT four = bigint<1>(4l) * GroupT::one();

    assert(two+five == three+four);

    GroupT a = GroupT::random_element();
    GroupT b = GroupT::random_element();

    std::cout << "GroupT::a  =  " << a << std::endl;
    std::cout << "GroupT::2a =  " << a.dbl() << std::endl;

    assert(one != zero);
    assert(a != zero);
    assert(a != one);

    assert(b != zero);
    assert(b != one);

    assert(a.dbl() == a + a);
    assert(b.dbl() == b + b);
    assert(one.add(two) == three);
    assert(two.add(one) == three);
    assert(a + b == b + a);
    assert(a - a == zero);
    assert(a - b == a + (-b));
    assert(a - b == (-b) + a);

    // handle special cases
    assert(zero + (-a) == -a);
    assert(zero - a == -a);
    assert(a - zero == a);
    assert(a + zero == a);
    assert(zero + a == a);

    assert((a + b).dbl() == (a + b) + (b + a));
    assert(bigint<1>("2") * (a + b) == (a + b) + (b + a));

    assert((rand1 * a) + (rand2 * a) == (randsum * a));

    assert(GroupT::order() * a == zero);
    assert(GroupT::order() * one == zero);
    assert((GroupT::order() * a) - a != zero);
    assert((GroupT::order() * one) - one != zero);

    test_mixed_add<GroupT>();
}

/*
template<typename GroupT>
void test_mul_by_q()
{
    GroupT a = GroupT::random_element();
    assert((GroupT::base_field_char()*a) == a.mul_by_q());
}
*/

template<typename GroupT>
void test_output()
{
    GroupT g = GroupT::zero();

    for (size_t i = 0; i < 1000; ++i)
    {
        std::stringstream ss;
        ss << g;
        GroupT gg;
        ss >> gg;
        assert(g == gg);
        /* use a random point in next iteration */
        g = GroupT::random_element();
    }
}

int main(void)
{


/*
std::cout << "alt_bn128" << std::endl;
    alt_bn128_pp::init_public_params();
    test_group<G1<alt_bn128_pp> >();
    test_output<G1<alt_bn128_pp> >();
    test_group<G2<alt_bn128_pp> >();
    test_output<G2<alt_bn128_pp> >();
    //test_mul_by_q<G2<alt_bn128_pp> >();
*/

std::cout << "bls48" << std::endl;

    bls48_pp::init_public_params();
    std::cout << "check_test_group_G1" << std::endl;
    test_group<G1<bls48_pp> >();
    std::cout << "check_test_output_G1" << std::endl;
    test_output<G1<bls48_pp> >();
    std::cout << "check_test_group_G2" << std::endl;
    test_group<G2<bls48_pp> >();
    //std::cout << "check_test_output_G2" << std::endl;
    //test_output<G2<bls48_pp> >();
    std::cout << "check_end" << std::endl;
    //test_mul_by_q<G2<bls48_pp> >();

}
