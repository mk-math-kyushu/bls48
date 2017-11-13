/** @file
 *****************************************************************************
 * @author     This file is part of libsnark, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#include "algebra/curves/bls48/bls48_pairing.hpp"
#include "algebra/curves/bls48/bls48_init.hpp"
#include "algebra/curves/bls48/bls48_g1.hpp"
#include "algebra/curves/bls48/bls48_g2.hpp"
#include <cassert>
#include "common/profiling.hpp"
#include <iostream>

namespace libsnark {


bool bls48_ate_G1_precomp::operator==(const bls48_ate_G1_precomp &other) const
{
    return (this->PX == other.PX &&
            this->PY == other.PY);
}

std::ostream& operator<<(std::ostream &out, const bls48_ate_G1_precomp &prec_P)
{
    out << prec_P.PX << OUTPUT_SEPARATOR << prec_P.PY;

    return out;
}

std::istream& operator>>(std::istream &in, bls48_ate_G1_precomp &prec_P)
{
    in >> prec_P.PX;
    consume_OUTPUT_SEPARATOR(in);
    in >> prec_P.PY;

    return in;
}

/*
bool  bls48_ate_ell_coeffs::operator==(const bls48_ate_ell_coeffs &other) const
{
    return (this->ell_0 == other.ell_0 &&
            this->ell_VW == other.ell_VW &&
            this->ell_VV == other.ell_VV);
}

std::ostream& operator<<(std::ostream &out, const bls48_ate_ell_coeffs &c)
{
    out << c.ell_0 << OUTPUT_SEPARATOR << c.ell_VW << OUTPUT_SEPARATOR << c.ell_VV;
    return out;
}

std::istream& operator>>(std::istream &in, bls48_ate_ell_coeffs &c)
{
    in >> c.ell_0;
    consume_OUTPUT_SEPARATOR(in);
    in >> c.ell_VW;
    consume_OUTPUT_SEPARATOR(in);
    in >> c.ell_VV;

    return in;
}
*/

bool bls48_ate_G2_precomp::operator==(const bls48_ate_G2_precomp &other) const
{
    return (this->QX == other.QX &&
            this->QY == other.QY);
             //&&this->coeffs == other.coeffs);
}

/*
std::ostream& operator<<(std::ostream& out, const bls48_ate_G2_precomp &prec_Q)
{
    out << prec_Q.QX << OUTPUT_SEPARATOR << prec_Q.QY << "\n";
    out << prec_Q.coeffs.size() << "\n";
    for (const bls48_ate_ell_coeffs &c : prec_Q.coeffs)
    {
        out << c << OUTPUT_NEWLINE;
    }
    return out;
}

std::istream& operator>>(std::istream& in, bls48_ate_G2_precomp &prec_Q)
{
    in >> prec_Q.QX;
    consume_OUTPUT_SEPARATOR(in);
    in >> prec_Q.QY;
    consume_newline(in);

    prec_Q.coeffs.clear();
    size_t s;
    in >> s;

    consume_newline(in);

    prec_Q.coeffs.reserve(s);

    for (size_t i = 0; i < s; ++i)
    {
        bls48_ate_ell_coeffs c;
        in >> c;
        consume_OUTPUT_NEWLINE(in);
        prec_Q.coeffs.emplace_back(c);
    }

    return in;
}
*/

/* final exponentiations */

bls48_Fq48 bls48_final_exponentiation_first_chunk(const bls48_Fq48 &elt)//easy part
{
    enter_block("Call to bls48_final_exponentiation_first_chunk");

    /*
      Computes result = elt^((q^24-1)*(q^8+1)).
      Follows, e.g., Beuchat et al page 9, by computing result as follows:
         elt^((q^24-1)*(q^8+1)) = (conj(elt) * elt^(-1))^(q^8+1)
      More precisely:
      A = conj(elt)
      B = elt.inverse()
      C = A * B
      D = C.Frobenius_map(8)
      result = D * C
    */

    const bls48_Fq48 A = bls48_Fq48(elt.c0,-elt.c1);
    const bls48_Fq48 B = elt.inverse();
    const bls48_Fq48 C = A * B;
    const bls48_Fq48 D = C.Frobenius_map(8);
    const bls48_Fq48 result = D * C;

    leave_block("Call to bls48_final_exponentiation_first_chunk");

    return result;
}

bls48_Fq48 bls48_exp_by_neg_z(const bls48_Fq48 &elt)
{
    enter_block("Call to bls48_exp_by_neg_z");

    bls48_Fq48 result = elt^bls48_final_exponent_z;
    if (!bls48_final_exponent_is_z_neg)
    {
        result = result.unitary_inverse();
    }

    leave_block("Call to bls48_exp_by_neg_z");

    return result;
}

bls48_Fq48 bls48_final_exponentiation_last_chunk(const bls48_Fq48 &elt)//hard part
{
    enter_block("Call to bls48_final_exponentiation_last_chunk");

    /*

      result =
      which equals

      result = elt^((q^16 - q^8 + 1)/r).

      Using the following addition chain:

      m_1 =  exp_by_neg_z(elt) // elt^(-z)

      result = V

    */


    const bls48_Fq48 m_1_temp = bls48_exp_by_neg_z(elt);
    const bls48_Fq48 m_2 = bls48_exp_by_neg_z(m_1_temp);
    const bls48_Fq48 m_1_temp_squared = m_1_temp.squared();
    const bls48_Fq48 m_1 = m_1_temp_squared.unitary_inverse();
    const bls48_Fq48 mu_15 = m_2*m_1*elt;
    const bls48_Fq48 mu_14 = bls48_exp_by_neg_z(mu_15);
    const bls48_Fq48 mu_13 = bls48_exp_by_neg_z(mu_14);
    const bls48_Fq48 mu_12 = bls48_exp_by_neg_z(mu_13);
    const bls48_Fq48 mu_11 = bls48_exp_by_neg_z(mu_12);
    const bls48_Fq48 mu_10 = bls48_exp_by_neg_z(mu_11);
    const bls48_Fq48 mu_9  = bls48_exp_by_neg_z(mu_10);
    const bls48_Fq48 mu_8  = bls48_exp_by_neg_z(mu_9);
    const bls48_Fq48 mu_15_inv = mu_15.unitary_inverse();
    const bls48_Fq48 mu_7  = bls48_exp_by_neg_z(mu_8)*mu_15_inv;
    const bls48_Fq48 mu_6 = bls48_exp_by_neg_z(mu_7);
    const bls48_Fq48 mu_5 = bls48_exp_by_neg_z(mu_6);
    const bls48_Fq48 mu_4 = bls48_exp_by_neg_z(mu_5);
    const bls48_Fq48 mu_3 = bls48_exp_by_neg_z(mu_4);
    const bls48_Fq48 mu_2 = bls48_exp_by_neg_z(mu_3);
    const bls48_Fq48 mu_1 = bls48_exp_by_neg_z(mu_2);
    const bls48_Fq48 mu_0_temp = bls48_exp_by_neg_z(mu_1);//
    const bls48_Fq48 m_cubed  = elt*elt*elt;//m_1_temp_squared*m_1_temp;//elt^3
    const bls48_Fq48 mu_0 = mu_0_temp*m_cubed;
    const bls48_Fq48 f_15_frob = mu_15.Frobenius_map(1);
    const bls48_Fq48 f_14 = f_15_frob*mu_14;
    const bls48_Fq48 f_14_frob = f_14.Frobenius_map(1);
    const bls48_Fq48 f_13 = f_14_frob*mu_13;
    const bls48_Fq48 f_13_frob = f_13.Frobenius_map(1);
    const bls48_Fq48 f_12 = f_13_frob*mu_12;
    const bls48_Fq48 f_12_frob = f_12.Frobenius_map(1);
    const bls48_Fq48 f_11 = f_12_frob*mu_11;
    const bls48_Fq48 f_11_frob = f_11.Frobenius_map(1);
    const bls48_Fq48 f_10 = f_11_frob*mu_10;
    const bls48_Fq48 f_10_frob = f_10.Frobenius_map(1);
    const bls48_Fq48 f_9 = f_10_frob*mu_9;
    const bls48_Fq48 f_9_frob = f_9.Frobenius_map(1);
    const bls48_Fq48 f_8 = f_9_frob*mu_8;
    const bls48_Fq48 f_8_frob = f_8.Frobenius_map(1);
    const bls48_Fq48 f_7 = f_8_frob*mu_7;
    const bls48_Fq48 f_7_frob = f_7.Frobenius_map(1);
    const bls48_Fq48 f_6 = f_7_frob*mu_6;
    const bls48_Fq48 f_6_frob = f_6.Frobenius_map(1);
    const bls48_Fq48 f_5 = f_6_frob*mu_5;
    const bls48_Fq48 f_5_frob = f_5.Frobenius_map(1);
    const bls48_Fq48 f_4 = f_5_frob*mu_4;
    const bls48_Fq48 f_4_frob = f_4.Frobenius_map(1);
    const bls48_Fq48 f_3 = f_4_frob*mu_3;
    const bls48_Fq48 f_3_frob = f_3.Frobenius_map(1);
    const bls48_Fq48 f_2 = f_3_frob*mu_2;
    const bls48_Fq48 f_2_frob = f_2.Frobenius_map(1);
    const bls48_Fq48 f_1 = f_2_frob*mu_1;
    const bls48_Fq48 f_1_frob = f_1.Frobenius_map(1);
    const bls48_Fq48 f = f_1_frob*mu_0;

    const bls48_Fq48 result = f;

    leave_block("Call to bls48_final_exponentiation_last_chunk");

    return result;
}

bls48_GT bls48_final_exponentiation(const bls48_Fq48 &elt)
{
    enter_block("Call to bls48_final_exponentiation");

    // OLD naive version:
        //bls48_GT result = elt^bls48_final_exponent;

    bls48_Fq48 A = bls48_final_exponentiation_first_chunk(elt);
    bls48_GT result = bls48_final_exponentiation_last_chunk(A);

    leave_block("Call to bls48_final_exponentiation");
    return result;
}

/* ate pairing */

bls48_ate_G1_precomp bls48_ate_precompute_G1(const bls48_G1& P)
{
    enter_block("Call to bls48_ate_precompute_G1");

    bls48_G1 Pcopy = P;
    Pcopy.to_affine_coordinates();
    //assert(Pcopy.Y.squared() == Pcopy.X.squared()*Pcopy.X + bls48_Fq::one());

    bls48_ate_G1_precomp result;
    result.PX = Pcopy.X;
    result.PY = Pcopy.Y;

    leave_block("Call to bls48_ate_precompute_G1");
    return result;
}


bls48_ate_G2_precomp bls48_ate_precompute_G2(const bls48_G2& Q)
{
    enter_block("Call to bls48_ate_precompute_G2");

    bls48_G2 Qcopy(Q);//bls48_G2 Qcopy = Q;
    Qcopy.to_affine_coordinates();
    //assert(Qcopy.Y.squared() == Qcopy.X.squared()*Qcopy.X + bls48_twist_coeff_b);

    //bls48_Fq two_inv = (bls48_Fq("2").inverse()); // could add to global params if needed

    bls48_ate_G2_precomp result;
    result.QX = Qcopy.X;
    result.QY = Qcopy.Y;



    leave_block("Call to bls48_ate_precompute_G2");
    return result;
}


bls48_Fq48 bls48_ate_miller_loop(const bls48_G1_precomp &prec_P,
                                     const bls48_G2_precomp &prec_Q)
{
    enter_block("Call to bls48_ate_miller_loop");

    bls48_Fq48 f = bls48_Fq48::one();//f <- 1
    bls48_G2 Q;//Q <- Q
    Q.X = prec_Q.QX;
    Q.Y = prec_Q.QY;
    Q.Z = bls48_Fq8::one();

    bls48_G2 T;//T <- Q
    T.X = prec_Q.QX;
    T.Y = prec_Q.QY;
    T.Z = bls48_Fq8::one();
    //std::cout << "T.Y = " << T.Y << std::endl;
    //std::cout << "T.Y + Q.Y = " << T.Y + Q.Y << std::endl;


    //assert(T.Y.squared() - T.X.squared()*T.X == bls48_twist_coeff_b);
    //long long j = -1;

    bls48_Fq48 l;
    l.c0.c0 = prec_P.PY*bls48_Fq8::one();
    bls48_Fq8 ramda;
    bls48_Fq8 c;
    bls48_Fq two = bls48_Fq::one() + bls48_Fq::one();
    bls48_Fq three = two + bls48_Fq::one();
    bls48_Fq8 u = bls48_Fq8(bls48_Fq4(bls48_Fq2(bls48_Fq::zero(),bls48_Fq::one()),bls48_Fq2::zero()),bls48_Fq4::zero());//bls48_Fq8 :: ( 1, u, v, uv, w, u*w, v*w, u*v*w )
    bls48_Fq8 xPu = prec_P.PX * u;


    //bool found_nonzero = false;

    //const bigint<bls48_Fr::num_limbs> &loop_count = bls48_ate_loop_count;//error
    //const bigint<bls48_Fq::num_limbs> &loop_count = bls48_ate_loop_count;

    //for (long i = loop_count.max_bits(); i >= 0; --i)
    for (int i = 31; i >= 0; --i)
    {
/*
        const bool bit = loop_count.test_bit(i);
        if (!found_nonzero)
        {
            // this skips the MSB itself
            found_nonzero |= bit;
            continue;
        }
*/

        /* code below gets executed for all bits (EXCEPT the MSB itself) of
           bls48_param_p (skipping leading zeros) in MSB to LSB
           order */


        //j = 2*j;

//////////
// l function
//////////
        ramda = two*T.Y;//2y
        ramda = ramda.inverse();//1/(2y)
        //assert(ramda*two*T.Y == bls48_Fq8::one());
        ramda = three*(T.X.squared())*ramda;//3*x^2/(2*y)
        c = ramda*T.X - T.Y;//y = ramda*x + -c
        l.c1.c0 = (-ramda)*xPu;//(-ramda*xP)*u where u in Fq8 and ramda is slope
        l.c1.c1 = c*u;//c*u where u is y-segment
        //std::cout << "l_TT(P)" << std::endl;
        //l.print();
        //std::cout << "end of print l_TT(P)" << std::endl;
//////////

        f = f.squared();//f <- f^2
        f = f*l;//f^2*l_{T,T}(P)
        T = T.dbl();
        T.to_affine_coordinates();//T <- 2T


        if (i==30 || i==10 )
        {
            //std::cout << "i = " << i << std::endl;
            ramda = Q.X - T.X;
            ramda = ramda.inverse();
            ramda = (Q.Y - T.Y)*ramda;//(yQ - yT)/(xQ - XT)
            c = ramda*Q.X - Q.Y;
            l.c1.c0 = (-ramda)*xPu;//(-ramda*xP)*u where u in Fq8 and ramda is slope
            l.c1.c1 = c*u;//c*u where c is - y-segment
            f = f*l;//f^2*l_{T,-Q}(P)
          T = T + Q;
          T.to_affine_coordinates();


          //j = j - 1;


        }



        if(i==0){
          //std::cout << "i = " << i << std::endl;
          ramda = Q.X - T.X;
          ramda = ramda.inverse();
          ramda = (Q.Y - T.Y)*ramda;//(yQ - yT)/(xQ - XT)
          c = ramda*Q.X - Q.Y;
          l.c1.c0 = (-ramda)*xPu;//(-ramda*xP)*u where u in Fq8 and ramda is slope
          l.c1.c1 = c*u;
          f = f*l;//f*l_{T,-Q}(P)


          //T = T + Q;
          //T.to_affine_coordinates();
          //assert(T == bls48_ate_loop_count*Q);

        }

        if(i==7){

          //std::cout << "i = " << i << std::endl;
          ramda = Q.X - T.X;
          ramda = ramda.inverse();
          ramda = (- Q.Y - T.Y)*ramda;//(y-Q - yT)/(x-Q - XT)
          c = ramda*Q.X + Q.Y;
          l.c1.c0 = (-ramda)*xPu;
          l.c1.c1 = c*u;
          f = f*l;//f*l_{T,Q}(P)
        T = T - Q;
        T.to_affine_coordinates();




        //j = j + 1;

        }

        //assert(T.Y.squared() - T.X.squared()*T.X == bls48_twist_coeff_b);

    }

    //std::cout << "j = " << j << std::endl;


    if (!bls48_ate_is_loop_count_neg)
    {
      //std::cout << "inverse" << std::endl;
    	f = f.unitary_inverse();
    }

    //std::cout << "f = " << f << std::endl;



    leave_block("Call to bls48_ate_miller_loop");
    return f;
}


bls48_Fq48 bls48_ate_pairing(const bls48_G1& P, const bls48_G2 &Q)
{
    enter_block("Call to bls48_ate_pairing");
    bls48_ate_G1_precomp prec_P = bls48_ate_precompute_G1(P);
    //assert(prec_P.PY.squared() - prec_P.PX.squared()*prec_P.PX == bls48_Fq::one());
    bls48_ate_G2_precomp prec_Q = bls48_ate_precompute_G2(Q);
    //assert(prec_Q.QY.squared() - prec_Q.QX.squared()*prec_Q.QX == bls48_twist_coeff_b);
    bls48_Fq48 result = bls48_ate_miller_loop(prec_P, prec_Q);
    leave_block("Call to bls48_ate_pairing");
    //std::cout << "after miller_loop result = " << result << std::endl;
    return result;
}

bls48_GT bls48_ate_reduced_pairing(const bls48_G1 &P, const bls48_G2 &Q)
{
    enter_block("Call to bls48_ate_reduced_pairing");
    const bls48_Fq48 f = bls48_ate_pairing(P, Q);
    const bls48_GT result = bls48_final_exponentiation(f);
    leave_block("Call to bls48_ate_reduced_pairing");
    //std::cout << "reduced result = " << result << std::endl;
    return result;
}

/* choice of pairing */

bls48_G1_precomp bls48_precompute_G1(const bls48_G1& P)
{
    return bls48_ate_precompute_G1(P);
}

bls48_G2_precomp bls48_precompute_G2(const bls48_G2& Q)
{
    return bls48_ate_precompute_G2(Q);
}

bls48_Fq48 bls48_miller_loop(const bls48_G1_precomp &prec_P,
                          const bls48_G2_precomp &prec_Q)
{
    return bls48_ate_miller_loop(prec_P, prec_Q);
}


bls48_Fq48 bls48_pairing(const bls48_G1& P,
                      const bls48_G2 &Q)
{
    return bls48_ate_pairing(P, Q);
}

bls48_GT bls48_reduced_pairing(const bls48_G1 &P,
                             const bls48_G2 &Q)
{
    return bls48_ate_reduced_pairing(P, Q);
}
} // libsnark
