/** @file
 *****************************************************************************
 Implementation of arithmetic in the finite field F[(p^2)^3].
 *****************************************************************************
 * @author     This file is part of libsnark, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef FP24_2OVER2OVER2OVER3_TCC_
#define FP24_2OVER2OVER2OVER3_TCC_
#include "algebra/fields/field_utils.hpp"

namespace libsnark {

template<mp_size_t n, const bigint<n>& modulus>
Fp8_model<n, modulus> Fp24_2over2over2over3_model<n,modulus>::mul_by_non_residue(const Fp8_model<n, modulus> &elt)
{
    return Fp8_model<n, modulus>(non_residue * elt);
}

template<mp_size_t n, const bigint<n>& modulus>
Fp24_2over2over2over3_model<n,modulus> Fp24_2over2over2over3_model<n,modulus>::zero()
{
    return Fp24_2over2over2over3_model<n, modulus>(my_Fp8::zero(), my_Fp8::zero(), my_Fp8::zero());
}

template<mp_size_t n, const bigint<n>& modulus>
Fp24_2over2over2over3_model<n,modulus> Fp24_2over2over2over3_model<n,modulus>::one()
{
    return Fp24_2over2over2over3_model<n, modulus>(my_Fp8::one(), my_Fp8::zero(), my_Fp8::zero());
}

template<mp_size_t n, const bigint<n>& modulus>
Fp24_2over2over2over3_model<n,modulus> Fp24_2over2over2over3_model<n,modulus>::random_element()
{
    Fp24_2over2over2over3_model<n, modulus> r;
    r.c0 = my_Fp8::random_element();
    r.c1 = my_Fp8::random_element();
    r.c2 = my_Fp8::random_element();

    return r;
}

template<mp_size_t n, const bigint<n>& modulus>
bool Fp24_2over2over2over3_model<n,modulus>::operator==(const Fp24_2over2over2over3_model<n,modulus> &other) const
{
    return (this->c0 == other.c0 && this->c1 == other.c1 && this->c2 == other.c2);
}

template<mp_size_t n, const bigint<n>& modulus>
bool Fp24_2over2over2over3_model<n,modulus>::operator!=(const Fp24_2over2over2over3_model<n,modulus> &other) const
{
    return !(operator==(other));
}

template<mp_size_t n, const bigint<n>& modulus>
Fp24_2over2over2over3_model<n,modulus> Fp24_2over2over2over3_model<n,modulus>::operator+(const Fp24_2over2over2over3_model<n,modulus> &other) const
{
    return Fp24_2over2over2over3_model<n,modulus>(this->c0 + other.c0,
                                       this->c1 + other.c1,
                                       this->c2 + other.c2);
}

template<mp_size_t n, const bigint<n>& modulus>
Fp24_2over2over2over3_model<n,modulus> Fp24_2over2over2over3_model<n,modulus>::operator-(const Fp24_2over2over2over3_model<n,modulus> &other) const
{
    return Fp24_2over2over2over3_model<n,modulus>(this->c0 - other.c0,
                                       this->c1 - other.c1,
                                       this->c2 - other.c2);
}

template<mp_size_t n, const bigint<n>& modulus>
Fp24_2over2over2over3_model<n, modulus> operator*(const Fp_model<n, modulus> &lhs, const Fp24_2over2over2over3_model<n, modulus> &rhs)
{
    return Fp24_2over2over2over3_model<n,modulus>(lhs*rhs.c0,
                                       lhs*rhs.c1,
                                       lhs*rhs.c2);
}

template<mp_size_t n, const bigint<n>& modulus>
Fp24_2over2over2over3_model<n, modulus> operator*(const Fp2_model<n, modulus> &lhs, const Fp24_2over2over2over3_model<n, modulus> &rhs)
{
    return Fp24_2over2over2over3_model<n,modulus>(lhs*rhs.c0,
                                       lhs*rhs.c1,
                                       lhs*rhs.c2);
}

template<mp_size_t n, const bigint<n>& modulus>
Fp24_2over2over2over3_model<n, modulus> operator*(const Fp4_model<n, modulus> &lhs, const Fp24_2over2over2over3_model<n, modulus> &rhs)
{
    return Fp24_2over2over2over3_model<n,modulus>(lhs*rhs.c0,
                                       lhs*rhs.c1,
                                       lhs*rhs.c2);
}

template<mp_size_t n, const bigint<n>& modulus>
Fp24_2over2over2over3_model<n, modulus> operator*(const Fp8_model<n, modulus> &lhs, const Fp24_2over2over2over3_model<n, modulus> &rhs)
{
    return Fp24_2over2over2over3_model<n,modulus>(lhs*rhs.c0,
                                       lhs*rhs.c1,
                                       lhs*rhs.c2);
}

template<mp_size_t n, const bigint<n>& modulus>
Fp24_2over2over2over3_model<n,modulus> Fp24_2over2over2over3_model<n,modulus>::operator*(const Fp24_2over2over2over3_model<n,modulus> &other) const
{
    /* Devegili OhEig Scott Dahab --- Multiplication and Squaring on Pairing-Friendly Fields.pdf; Section 4 (Karatsuba) */

    const my_Fp8 &A = other.c0, &B = other.c1, &C = other.c2,
                 &a = this->c0, &b = this->c1, &c = this->c2;
    const my_Fp8 aA = a*A;
    const my_Fp8 bB = b*B;
    const my_Fp8 cC = c*C;

    return Fp24_2over2over2over3_model<n,modulus>(aA + Fp24_2over2over2over3_model<n,modulus>::mul_by_non_residue((b+c)*(B+C)-bB-cC),
                                       (a+b)*(A+B)-aA-bB+Fp24_2over2over2over3_model<n,modulus>::mul_by_non_residue(cC),
                                       (a+c)*(A+C)-aA+bB-cC);
}

template<mp_size_t n, const bigint<n>& modulus>
Fp24_2over2over2over3_model<n,modulus> Fp24_2over2over2over3_model<n,modulus>::operator-() const
{
    return Fp24_2over2over2over3_model<n,modulus>(-this->c0,
                                       -this->c1,
                                       -this->c2);
}

template<mp_size_t n, const bigint<n>& modulus>
Fp24_2over2over2over3_model<n,modulus> Fp24_2over2over2over3_model<n,modulus>::squared() const
{
    /* Devegili OhEig Scott Dahab --- Multiplication and Squaring on Pairing-Friendly Fields.pdf; Section 4 (CH-SQR2) */

    const my_Fp8 &a = this->c0, &b = this->c1, &c = this->c2;
    const my_Fp8 s0 = a.squared();
    const my_Fp8 ab = a*b;
    const my_Fp8 s1 = ab + ab;
    const my_Fp8 s2 = (a - b + c).squared();
    const my_Fp8 bc = b*c;
    const my_Fp8 s3 = bc + bc;
    const my_Fp8 s4 = c.squared();

    return Fp24_2over2over2over3_model<n,modulus>(s0 + Fp24_2over2over2over3_model<n,modulus>::mul_by_non_residue(s3),
                                       s1 + Fp24_2over2over2over3_model<n,modulus>::mul_by_non_residue(s4),
                                       s1 + s2 + s3 - s0 - s4);
}

template<mp_size_t n, const bigint<n>& modulus>
Fp24_2over2over2over3_model<n,modulus> Fp24_2over2over2over3_model<n,modulus>::inverse() const
{
    /* From "High-Speed Software Implementation of the Optimal Ate Pairing over Barreto-Naehrig Curves"; Algorithm 17 */

    const my_Fp8 &a = this->c0, &b = this->c1, &c = this->c2;
    const my_Fp8 t0 = a.squared();
    const my_Fp8 t1 = b.squared();
    const my_Fp8 t2 = c.squared();
    const my_Fp8 t3 = a*b;
    const my_Fp8 t4 = a*c;
    const my_Fp8 t5 = b*c;
    const my_Fp8 c0 = t0 - Fp24_2over2over2over3_model<n,modulus>::mul_by_non_residue(t5);
    const my_Fp8 c1 = Fp24_2over2over2over3_model<n,modulus>::mul_by_non_residue(t2) - t3;
    const my_Fp8 c2 = t1 - t4; // typo in paper referenced above. should be "-" as per Scott, but is "*"
    const my_Fp8 t6 = (a * c0 + Fp24_2over2over2over3_model<n,modulus>::mul_by_non_residue((c * c1 + b * c2))).inverse();
    return Fp24_2over2over2over3_model<n,modulus>(t6 * c0, t6 * c1, t6 * c2);
}

template<mp_size_t n, const bigint<n>& modulus>
Fp24_2over2over2over3_model<n,modulus> Fp24_2over2over2over3_model<n,modulus>::Frobenius_map(unsigned long power) const
{
    return Fp24_2over2over2over3_model<n,modulus>(c0.Frobenius_map(power),
                                       Frobenius_coeffs_c1[power % 24] * c1.Frobenius_map(power),
                                       Frobenius_coeffs_c2[power % 24] * c2.Frobenius_map(power));
}

template<mp_size_t n, const bigint<n>& modulus>
template<mp_size_t m>
Fp24_2over2over2over3_model<n,modulus> Fp24_2over2over2over3_model<n,modulus>::operator^(const bigint<m> &pow) const
{
    return power<Fp24_2over2over2over3_model<n, modulus>, m>(*this, pow);
}

template<mp_size_t n, const bigint<n>& modulus>
std::ostream& operator<<(std::ostream &out, const Fp24_2over2over2over3_model<n, modulus> &el)
{
    out << el.c0 << OUTPUT_SEPARATOR << el.c1 << OUTPUT_SEPARATOR << el.c2;
    return out;
}

template<mp_size_t n, const bigint<n>& modulus>
std::istream& operator>>(std::istream &in, Fp24_2over2over2over3_model<n, modulus> &el)
{
    in >> el.c0 >> el.c1 >> el.c2;
    return in;
}

template<mp_size_t n, const bigint<n>& modulus>
std::ostream& operator<<(std::ostream& out, const std::vector<Fp24_2over2over2over3_model<n, modulus> > &v)
{
    out << v.size() << "\n";
    for (const Fp24_2over2over2over3_model<n, modulus>& t : v)
    {
        out << t << OUTPUT_NEWLINE;
    }

    return out;
}

template<mp_size_t n, const bigint<n>& modulus>
std::istream& operator>>(std::istream& in, std::vector<Fp24_2over2over2over3_model<n, modulus> > &v)
{
    v.clear();

    size_t s;
    in >> s;

    char b;
    in.read(&b, 1);

    v.reserve(s);

    for (size_t i = 0; i < s; ++i)
    {
        Fp24_2over2over2over3_model<n, modulus> el;
        in >> el;
        v.emplace_back(el);
    }

    return in;
}

} // libsnark
#endif // FP24_2OVER2OVER2OVER3_TCC_
