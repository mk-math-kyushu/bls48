/** @file
 *****************************************************************************

 Implementation of interfaces for the (extension) field Fp8.

 See fp8.hpp .

 *****************************************************************************
 * @author     This file is part of libsnark, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef FP8_TCC_
#define FP8_TCC_

#include "algebra/fields/field_utils.hpp"
#include "algebra/scalar_multiplication/wnaf.hpp"

namespace libsnark {


//mul_non_residue is used for comp multiplication by Karatsuba
template<mp_size_t n, const bigint<n>& modulus>
Fp4_model<n, modulus> Fp8_model<n, modulus>::mul_by_non_residue(const Fp4_model<n, modulus> &elt)//2->4, 4->8, 2->4
{
    return Fp4_model<n, modulus>(non_residue * elt);//
}
//////////////////////////////////////


template<mp_size_t n, const bigint<n>& modulus>
Fp8_model<n, modulus> Fp8_model<n, modulus>::zero()//define zero in Fp8
{
    return Fp8_model<n,modulus>(my_Fp4::zero(),
                                my_Fp4::zero());
}

template<mp_size_t n, const bigint<n>& modulus>
Fp8_model<n, modulus> Fp8_model<n, modulus>::one()//define one in Fp8
{
    return Fp8_model<n,modulus>(my_Fp4::one(),
                                my_Fp4::zero());
}

template<mp_size_t n, const bigint<n>& modulus>
Fp8_model<n,modulus> Fp8_model<n,modulus>::random_element()
{
    Fp8_model<n, modulus> r;
    r.c0 = my_Fp4::random_element();
    r.c1 = my_Fp4::random_element();

    return r;
}

template<mp_size_t n, const bigint<n>& modulus>
bool Fp8_model<n,modulus>::operator==(const Fp8_model<n,modulus> &other) const
{
    return (this->c0 == other.c0 && this->c1 == other.c1);
}

template<mp_size_t n, const bigint<n>& modulus>
bool Fp8_model<n,modulus>::operator!=(const Fp8_model<n,modulus> &other) const
{
    return !(operator==(other));
}

template<mp_size_t n, const bigint<n>& modulus>
Fp8_model<n,modulus> Fp8_model<n,modulus>::operator+(const Fp8_model<n,modulus> &other) const//define + in Fp8
{
    return Fp8_model<n,modulus>(this->c0 + other.c0,
                                this->c1 + other.c1);
}

template<mp_size_t n, const bigint<n>& modulus>
Fp8_model<n,modulus> Fp8_model<n,modulus>::operator-(const Fp8_model<n,modulus> &other) const//define - in Fp8
{
    return Fp8_model<n,modulus>(this->c0 - other.c0,
                                this->c1 - other.c1);
}

template<mp_size_t n, const bigint<n>& modulus>
Fp8_model<n, modulus> operator*(const Fp_model<n, modulus> &lhs, const Fp8_model<n, modulus> &rhs)//define Fp*Fp8 in Fp8
{
    return Fp8_model<n,modulus>(lhs*rhs.c0,
                                lhs*rhs.c1);
}

template<mp_size_t n, const bigint<n>& modulus>
Fp8_model<n, modulus> operator*(const Fp2_model<n, modulus> &lhs, const Fp8_model<n, modulus> &rhs)//define Fp2*Fp8 in Fp8
{
    return Fp8_model<n,modulus>(lhs*rhs.c0,
                                lhs*rhs.c1);
}

//add
template<mp_size_t n, const bigint<n>& modulus>
Fp8_model<n, modulus> operator*(const Fp4_model<n, modulus> &lhs, const Fp8_model<n, modulus> &rhs)//define Fp4*Fp8 in Fp8
{
    return Fp8_model<n,modulus>(lhs*rhs.c0,
                                lhs*rhs.c1);
}
/////



//
template<mp_size_t n, const bigint<n>& modulus>
Fp8_model<n,modulus> Fp8_model<n,modulus>::operator*(const Fp8_model<n,modulus> &other) const//define * in Fp8 (Karatsuba)
{
    const my_Fp4 &B = other.c1, &A = other.c0,
        &b = this->c1, &a = this->c0;
    const my_Fp4 aA = a*A;
    const my_Fp4 bB = b*B;

    const my_Fp4 beta_bB = Fp8_model<n,modulus>::mul_by_non_residue(bB);
    return Fp8_model<n,modulus>(aA + beta_bB,
                                (a+b)*(A+B) - aA  - bB);


}
/////



//?
template<mp_size_t n, const bigint<n>& modulus>
Fp8_model<n,modulus> Fp8_model<n,modulus>::mul_by_023(const Fp8_model<n,modulus> &other) const
{
    //Devegili OhEig Scott Dahab --- Multiplication and Squaring on Pairing-Friendly Fields.pdf; Section 3 (Karatsuba)
    assert(other.c0.c1.is_zero());

    const my_Fp4 &B = other.c1, &A = other.c0,
        &b = this->c1, &a = this->c0;
    const my_Fp4 aA = my_Fp2(a.c0 * A.c0, a.c1 * A.c0);
    const my_Fp4 bB = b*B;

    const my_Fp4 beta_bB = Fp4_model<n,modulus>::mul_by_non_residue(bB);
    return Fp8_model<n,modulus>(aA + beta_bB,
                                (a+b)*(A+B) - aA  - bB);
}


template<mp_size_t n, const bigint<n>& modulus>
Fp8_model<n,modulus> Fp8_model<n,modulus>::operator-() const//change the sign of el in Fp8
{
    return Fp8_model<n,modulus>(-this->c0,
                                -this->c1);
}


///
template<mp_size_t n, const bigint<n>& modulus>
Fp8_model<n,modulus> Fp8_model<n,modulus>::squared() const
{
    /* Devegili OhEig Scott Dahab --- Multiplication and Squaring on Pairing-Friendly Fields.pdf; Section 3 (Complex) */

/*
    const my_Fp4 &a = this->c0, &b = this->c1;
    const my_Fp4 asq = a.squared();
    const my_Fp4 bsq = b.squared();

    return Fp8_model<n,modulus>(asq + Fp8_model<n, modulus>::mul_by_non_residue(bsq),
                                             (a + b).squared() - asq - bsq);
*/

    const my_Fp4 &b = this->c1, &a = this->c0;
    const my_Fp4 ab = a * b;

    return Fp8_model<n,modulus>((a+b)*(a+Fp8_model<n,modulus>::mul_by_non_residue(b))-ab-Fp8_model<n,modulus>::mul_by_non_residue(ab),
                                ab + ab);

}


template<mp_size_t n, const bigint<n>& modulus>
Fp8_model<n,modulus> Fp8_model<n,modulus>::inverse() const
{
    /* From "High-Speed Software Implementation of the Optimal Ate Pairing over Barreto-Naehrig Curves"; Algorithm 8 */
    const my_Fp4 &b = this->c1, &a = this->c0;
    const my_Fp4 t1 = b.squared();
    const my_Fp4 t0 = a.squared() - Fp8_model<n,modulus>::mul_by_non_residue(t1);
    const my_Fp4 new_t1 = t0.inverse();

    return Fp8_model<n,modulus>(a * new_t1, - (b * new_t1));
}



template<mp_size_t n, const bigint<n>& modulus>
Fp8_model<n,modulus> Fp8_model<n,modulus>::Frobenius_map(unsigned long power) const
{
    return Fp8_model<n,modulus>(c0.Frobenius_map(power),
                                Frobenius_coeffs_c1[power % 8] * c1.Frobenius_map(power));
}



template<mp_size_t n, const bigint<n>& modulus>
Fp8_model<n,modulus> Fp8_model<n,modulus>::unitary_inverse() const
{
    return Fp8_model<n,modulus>(this->c0,
                                -this->c1);
}


template<mp_size_t n, const bigint<n>& modulus>
Fp8_model<n,modulus> Fp8_model<n,modulus>::sqrt() const
{
    Fp8_model<n,modulus> one = Fp8_model<n,modulus>::one();

    size_t v = Fp8_model<n,modulus>::s;
    Fp8_model<n,modulus> z = Fp8_model<n,modulus>::nqr_to_t;
    Fp8_model<n,modulus> w = (*this)^Fp8_model<n,modulus>::t_minus_1_over_2;
    Fp8_model<n,modulus> x = (*this) * w;
    Fp8_model<n,modulus> b = x * w; // b = (*this)^t

#if DEBUG
    // check if square with euler's criterion
    Fp8_model<n,modulus> check = b;
    for (size_t i = 0; i < v-1; ++i)
    {
        check = check.squared();
    }
    if (check != one)
    {
        assert(0);
    }
#endif

    // compute square root with Tonelli--Shanks
    // (does not terminate if not a square!)

    while (b != one)
    {
        size_t m = 0;
        Fp8_model<n,modulus> b2m = b;
        while (b2m != one)
        {
            /* invariant: b2m = b^(2^m) after entering this loop */
            b2m = b2m.squared();
            m += 1;
        }

        int j = v-m-1;
        w = z;
        while (j > 0)
        {
            w = w.squared();
            --j;
        } // w = z^2^(v-m-1)

        z = w.squared();
        b = b * z;
        x = x * w;
        v = m;
    }

    return x;
}

/*
template<mp_size_t n, const bigint<n>& modulus>
Fp4_model<n,modulus> Fp4_model<n,modulus>::cyclotomic_squared() const
{
    const my_Fp2 A = this->c1.squared();
    const my_Fp2 B = this->c1 + this->c0;
    const my_Fp2 C = B.squared() - A;
    const my_Fp2 D = Fp4_model<n,modulus>::mul_by_non_residue(A); // Fp2(A.c1 * non_residue, A.c0)
    const my_Fp2 E = C - D;
    const my_Fp2 F = D + D + my_Fp2::one();
    const my_Fp2 G = E - my_Fp2::one();

    return Fp4_model<n,modulus>(F, G);
}
*/

/*
template<mp_size_t n, const bigint<n>& modulus>
template<mp_size_t m>
Fp4_model<n, modulus> Fp4_model<n,modulus>::cyclotomic_exp(const bigint<m> &exponent) const
{
    Fp4_model<n,modulus> res = Fp4_model<n,modulus>::one();
    Fp4_model<n,modulus> this_inverse = this->unitary_inverse();

    bool found_nonzero = false;
    std::vector<long> NAF = find_wnaf(1, exponent);

    for (long i = NAF.size() - 1; i >= 0; --i)
    {
        if (found_nonzero)
        {
            res = res.cyclotomic_squared();
        }

        if (NAF[i] != 0)
        {
            found_nonzero = true;

            if (NAF[i] > 0)
            {
                res = res * (*this);
            }
            else
            {
                res = res * this_inverse;
            }
        }
    }

    return res;
}
*/

template<mp_size_t n, const bigint<n>& modulus, mp_size_t m>
Fp8_model<n, modulus> operator^(const Fp8_model<n, modulus> &self, const bigint<m> &exponent)
{
    return power<Fp8_model<n, modulus> >(self, exponent);
}

template<mp_size_t n, const bigint<n>& modulus, mp_size_t m, const bigint<m>& modulus_p>
Fp8_model<n, modulus> operator^(const Fp8_model<n, modulus> &self, const Fp_model<m, modulus_p> &exponent)
{
    return self^(exponent.as_bigint());
}

template<mp_size_t n, const bigint<n>& modulus>
std::ostream& operator<<(std::ostream &out, const Fp8_model<n, modulus> &el)
{
    out << el.c0 << OUTPUT_SEPARATOR << el.c1;
    return out;
}

template<mp_size_t n, const bigint<n>& modulus>
std::istream& operator>>(std::istream &in, Fp8_model<n, modulus> &el)
{
    in >> el.c0 >> el.c1;
    return in;
}

} // libsnark

#endif // FP8_TCC_
