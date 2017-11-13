/**
 *****************************************************************************
 * @author     This file is part of libsnark, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/
#include "common/profiling.hpp"

#include "algebra/fields/fp.hpp"
#include "algebra/fields/fp8.hpp"
#include "algebra/fields/fp24_2over2over2over3.hpp"
#include "algebra/fields/fp48_2over2over2over3over2.hpp"
#include "algebra/curves/bls48/bls48_pp.hpp"
#include "iostream"
///////
//#include "algebra/fields/fp24_2over2over2over3.hpp"
//#include "algebra/fields/fp48_2over2over2over3over2.hpp"
/////////

using namespace libsnark;

template<typename FieldT>
void test_field()
{
    bigint<1> rand1 = bigint<1>("76749407");
    bigint<1> rand2 = bigint<1>("44410867");
    bigint<1> randsum = bigint<1>("121160274");



std::cout<< "Fp_"<< FieldT::extension_degree() << std::endl;


    FieldT zero = FieldT::zero();
    std::cout << "zero = " << FieldT::zero() << std::endl;
    FieldT one = FieldT::one();
    std::cout << "one = " << FieldT::one() << std::endl;
    FieldT a = FieldT::random_element();
    std::cout << "a = " << a << std::endl;
    FieldT a_ser;
    std::cout << "a_ser = " << a_ser << std::endl;
    a_ser = reserialize<FieldT>(a);
    assert(a_ser == a);



    FieldT b = FieldT::random_element();
    FieldT c = FieldT::random_element();
    FieldT d = FieldT::random_element();

    assert(a != zero);
    assert(a != one);


    assert(a * a == a.squared());
std::cout << "a*a = " << a*a << std::endl;
std::cout << "a^2 = " << a.squared() << std::endl;

    assert((a + b).squared() == a.squared() + a*b + b*a + b.squared());
    assert((a + b)*(c + d) == a*c + a*d + b*c + b*d);
    assert(a - b == a + (-b));
    assert(a - b == (-b) + a);


    assert((a ^ rand1) * (a ^ rand2) == (a^randsum));
    assert(a * a.inverse() == one);
    assert((a + b) * c.inverse() == a * c.inverse() + (b.inverse() * c).inverse());

}

template<typename FieldT>
void test_sqrt()
{
    for (size_t i = 0; i < 100; ++i)
    {
        FieldT a = FieldT::random_element();
        FieldT asq = a.squared();
        assert(asq.sqrt() == a || asq.sqrt() == -a);
    }
}

template<typename FieldT>
void test_two_squarings()
{
    FieldT a = FieldT::random_element();
    assert(a.squared() == a * a);
    assert(a.squared() == a.squared_complex());
    assert(a.squared() == a.squared_karatsuba());
}

template<typename FieldT>
void test_Frobenius()
{

    //std::cout << "Frobenius_map " << std::endl;

    //std::cout << "base_field_char = " << FieldT::base_field_char() << std::endl;

    FieldT a = FieldT::random_element();
    assert(a.Frobenius_map(0) == a);
    //std::cout << "a = " << a << std::endl;
    FieldT a_q = a ^ FieldT::base_field_char();
    for (size_t power = 1; power < 10; ++power)
    {

        const FieldT a_qi = a.Frobenius_map(power);
        /*
        if(power == 1 || power == 8){
          std::cout << "power = " << power << std::endl;
          std::cout << "a_qi = " << a_qi << std::endl;
          std::cout << "a_q = " << a_q << std::endl;
        }
        */
        assert(a_qi == a_q);

        a_q = a_q ^ FieldT::base_field_char();
    }
}





template<typename FieldT>
void test_unitary_inverse()
{
    assert(FieldT::extension_degree() % 2 == 0);
    FieldT a = FieldT::random_element();
    FieldT aqcubed_minus1 = a.Frobenius_map(FieldT::extension_degree()/2) * a.inverse();
    assert(aqcubed_minus1.inverse() == aqcubed_minus1.unitary_inverse());
}


/*
template<typename FieldT>
void test_cyclotomic_squaring();

template<>
void test_cyclotomic_squaring<Fqk<edwards_pp> >()
{
    typedef Fqk<edwards_pp> FieldT;
    assert(FieldT::extension_degree() % 2 == 0);
    FieldT a = FieldT::random_element();
    FieldT a_unitary = a.Frobenius_map(FieldT::extension_degree()/2) * a.inverse();
    // beta = a^((q^(k/2)-1)*(q+1))
    FieldT beta = a_unitary.Frobenius_map(1) * a_unitary;
    assert(beta.cyclotomic_squared() == beta.squared());
}

template<>
void test_cyclotomic_squaring<Fqk<mnt4_pp> >()
{
    typedef Fqk<mnt4_pp> FieldT;
    assert(FieldT::extension_degree() % 2 == 0);
    FieldT a = FieldT::random_element();
    FieldT a_unitary = a.Frobenius_map(FieldT::extension_degree()/2) * a.inverse();
    // beta = a^(q^(k/2)-1)
    FieldT beta = a_unitary;
    assert(beta.cyclotomic_squared() == beta.squared());
}

template<>
void test_cyclotomic_squaring<Fqk<mnt6_pp> >()
{
    typedef Fqk<mnt6_pp> FieldT;
    assert(FieldT::extension_degree() % 2 == 0);
    FieldT a = FieldT::random_element();
    FieldT a_unitary = a.Frobenius_map(FieldT::extension_degree()/2) * a.inverse();
    // beta = a^((q^(k/2)-1)*(q+1))
    FieldT beta = a_unitary.Frobenius_map(1) * a_unitary;
    assert(beta.cyclotomic_squared() == beta.squared());
}
*/

template<typename ppT>
void test_all_fields()
{


    test_field<Fr<ppT> >();

    test_field<Fq<ppT> >();

    test_field<Fqe<ppT> >();

    test_field<Fqk<ppT> >();



    test_Frobenius<Fqe<ppT> >();
    test_Frobenius<Fqk<ppT> >();
    test_Frobenius<Fqt<ppT> >();



    //test_sqrt<Fr<ppT> >();

    //test_sqrt<Fq<ppT> >();
    //test_sqrt<Fqe<ppT> >();






    //test_two_squarings<Fqe<ppT> >();


    //test_unitary_inverse<Fqk<ppT> >();

}



int main(void)
{


    bls48_pp::init_public_params();

    test_all_fields<bls48_pp>();





}
