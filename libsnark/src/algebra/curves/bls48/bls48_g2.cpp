/** @file
 *****************************************************************************
 * @author     This file is part of libsnark, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#include "algebra/curves/bls48/bls48_g2.hpp"

namespace libsnark {

#ifdef PROFILE_OP_COUNTS
long long bls48_G2::add_cnt = 0;
long long bls48_G2::dbl_cnt = 0;
#endif

std::vector<size_t> bls48_G2::wnaf_window_table;
std::vector<size_t> bls48_G2::fixed_base_exp_window_table;
bls48_G2 bls48_G2::G2_zero;
bls48_G2 bls48_G2::G2_one;


bls48_G2::bls48_G2()
{
    this->X = G2_zero.X;
    this->Y = G2_zero.Y;
    this->Z = G2_zero.Z;
}

/*
bls48_Fq8 bls48_G2::mul_by_b(const bls48_Fq8 &elt)
{
    return bls48_Fq8(bls48_twist_mul_by_b_c0 * elt.c0, bls48_twist_mul_by_b_c1 * elt.c1);
}
*/

void bls48_G2::print() const 
{
    if (this->is_zero())
    {
        printf("O\n");
    }
    else
    {
        bls48_G2 copy(*this);
        copy.to_affine_coordinates();
        gmp_printf("(((%Nd*u + %Nd)*v + %Nd*u + %Nd)*w + ((%Nd*u + %Nd)v + %Nd*u + %Nd) , ((%Nd*u + %Nd)*v + %Nd*u + %Nd)*w + (%Nd*u + %Nd)v + %Nd*u + %Nd))\n",
                   copy.X.c1.c1.c1.as_bigint().data, bls48_Fq::num_limbs,
                   copy.X.c1.c1.c0.as_bigint().data, bls48_Fq::num_limbs,
                   copy.X.c1.c0.c1.as_bigint().data, bls48_Fq::num_limbs,
                   copy.X.c1.c0.c0.as_bigint().data, bls48_Fq::num_limbs,
                   copy.X.c0.c1.c1.as_bigint().data, bls48_Fq::num_limbs,
                   copy.X.c0.c1.c0.as_bigint().data, bls48_Fq::num_limbs,
                   copy.X.c0.c0.c1.as_bigint().data, bls48_Fq::num_limbs,
                   copy.X.c0.c0.c0.as_bigint().data, bls48_Fq::num_limbs,
                   copy.Y.c1.c1.c1.as_bigint().data, bls48_Fq::num_limbs,
                   copy.Y.c1.c1.c0.as_bigint().data, bls48_Fq::num_limbs,
                   copy.Y.c1.c0.c1.as_bigint().data, bls48_Fq::num_limbs,
                   copy.Y.c1.c0.c0.as_bigint().data, bls48_Fq::num_limbs,
                   copy.Y.c0.c1.c1.as_bigint().data, bls48_Fq::num_limbs,
                   copy.Y.c0.c1.c0.as_bigint().data, bls48_Fq::num_limbs,
                   copy.Y.c0.c0.c1.as_bigint().data, bls48_Fq::num_limbs,
                   copy.Y.c0.c0.c0.as_bigint().data, bls48_Fq::num_limbs);
    }
}

void bls48_G2::print_coordinates() const
{
    if (this->is_zero())
    {
        printf("O\n");
    }
    else
    {
        gmp_printf("(((%Nd*u + %Nd)*v + %Nd*u + %Nd)*w + (%Nd*u + %Nd)v + %Nd*u + %Nd) : ((%Nd*u + %Nd)*v + %Nd*u + %Nd)*w + (%Nd*u + %Nd)v + %Nd*u + %Nd) : ((%Nd*u + %Nd)*v + %Nd*u + %Nd)*w + (%Nd*u + %Nd)v + %Nd*u + %Nd))\n",
                   this->X.c1.c1.c1.as_bigint().data, bls48_Fq::num_limbs,
                   this->X.c1.c1.c0.as_bigint().data, bls48_Fq::num_limbs,
                   this->X.c1.c0.c1.as_bigint().data, bls48_Fq::num_limbs,
                   this->X.c1.c0.c0.as_bigint().data, bls48_Fq::num_limbs,
                   this->X.c0.c1.c1.as_bigint().data, bls48_Fq::num_limbs,
                   this->X.c0.c1.c0.as_bigint().data, bls48_Fq::num_limbs,
                   this->X.c0.c0.c1.as_bigint().data, bls48_Fq::num_limbs,
                   this->X.c0.c0.c0.as_bigint().data, bls48_Fq::num_limbs,

                   this->Y.c1.c1.c1.as_bigint().data, bls48_Fq::num_limbs,
                   this->Y.c1.c1.c0.as_bigint().data, bls48_Fq::num_limbs,
                   this->Y.c1.c0.c1.as_bigint().data, bls48_Fq::num_limbs,
                   this->Y.c1.c0.c0.as_bigint().data, bls48_Fq::num_limbs,
                   this->Y.c0.c1.c1.as_bigint().data, bls48_Fq::num_limbs,
                   this->Y.c0.c1.c0.as_bigint().data, bls48_Fq::num_limbs,
                   this->Y.c0.c0.c1.as_bigint().data, bls48_Fq::num_limbs,
                   this->Y.c0.c0.c0.as_bigint().data, bls48_Fq::num_limbs,

                   this->Z.c1.c1.c1.as_bigint().data, bls48_Fq::num_limbs,
                   this->Z.c1.c1.c0.as_bigint().data, bls48_Fq::num_limbs,
                   this->Z.c1.c0.c1.as_bigint().data, bls48_Fq::num_limbs,
                   this->Z.c1.c0.c0.as_bigint().data, bls48_Fq::num_limbs,
                   this->Z.c0.c1.c1.as_bigint().data, bls48_Fq::num_limbs,
                   this->Z.c0.c1.c0.as_bigint().data, bls48_Fq::num_limbs,
                   this->Z.c0.c0.c1.as_bigint().data, bls48_Fq::num_limbs,
                   this->Z.c0.c0.c0.as_bigint().data, bls48_Fq::num_limbs);
    }
}

void bls48_G2::to_affine_coordinates()
{
    if (this->is_zero())
    {
        this->X = bls48_Fq8::zero();
        this->Y = bls48_Fq8::one();
        this->Z = bls48_Fq8::zero();
    }
    else
    {
        bls48_Fq8 Z_inv = Z.inverse();
        bls48_Fq8 Z2_inv = Z_inv.squared();
        bls48_Fq8 Z3_inv = Z2_inv * Z_inv;
        this->X = this->X * Z2_inv;
        this->Y = this->Y * Z3_inv;
        this->Z = bls48_Fq8::one();
    }
}

void bls48_G2::to_special()
{
    this->to_affine_coordinates();
}

bool bls48_G2::is_special() const
{
    return (this->is_zero() || this->Z == bls48_Fq8::one());
}

bool bls48_G2::is_zero() const
{
    return (this->Z.is_zero());
}

bool bls48_G2::operator==(const bls48_G2 &other) const
{
    if (this->is_zero())
    {
        return other.is_zero();
    }

    if (other.is_zero())
    {
        return false;
    }

    /* now neither is O */

    // using Jacobian coordinates so:
    // (X1:Y1:Z1) = (X2:Y2:Z2)
    // iff
    // X1/Z1^2 == X2/Z2^2 and Y1/Z1^3 == Y2/Z2^3
    // iff
    // X1 * Z2^2 == X2 * Z1^2 and Y1 * Z2^3 == Y2 * Z1^3

    bls48_Fq8 Z1_squared = (this->Z).squared();
    bls48_Fq8 Z2_squared = (other.Z).squared();

    if ((this->X * Z2_squared) != (other.X * Z1_squared))
    {
        return false;
    }

    bls48_Fq8 Z1_cubed = (this->Z) * Z1_squared;
    bls48_Fq8 Z2_cubed = (other.Z) * Z2_squared;

    if ((this->Y * Z2_cubed) != (other.Y * Z1_cubed))
    {
        return false;
    }

    return true;
}

bool bls48_G2::operator!=(const bls48_G2& other) const
{
    return !(operator==(other));
}

bls48_G2 bls48_G2::operator+(const bls48_G2 &other) const
{
    // handle special cases having to do with O
    if (this->is_zero())
    {
        return other;
    }

    if (other.is_zero())
    {
        return *this;
    }

    // no need to handle points of order 2,4
    // (they cannot exist in a prime-order subgroup)

    // check for doubling case

    // using Jacobian coordinates so:
    // (X1:Y1:Z1) = (X2:Y2:Z2)
    // iff
    // X1/Z1^2 == X2/Z2^2 and Y1/Z1^3 == Y2/Z2^3
    // iff
    // X1 * Z2^2 == X2 * Z1^2 and Y1 * Z2^3 == Y2 * Z1^3

    bls48_Fq8 Z1Z1 = (this->Z).squared();
    bls48_Fq8 Z2Z2 = (other.Z).squared();

    bls48_Fq8 U1 = this->X * Z2Z2;
    bls48_Fq8 U2 = other.X * Z1Z1;

    bls48_Fq8 Z1_cubed = (this->Z) * Z1Z1;
    bls48_Fq8 Z2_cubed = (other.Z) * Z2Z2;

    bls48_Fq8 S1 = (this->Y) * Z2_cubed;      // S1 = Y1 * Z2 * Z2Z2
    bls48_Fq8 S2 = (other.Y) * Z1_cubed;      // S2 = Y2 * Z1 * Z1Z1

    if (U1 == U2 && S1 == S2)
    {
        // dbl case; nothing of above can be reused
        return this->dbl();
    }

    // rest of add case
    bls48_Fq8 H = U2 - U1;                            // H = U2-U1
    bls48_Fq8 S2_minus_S1 = S2-S1;
    bls48_Fq8 I = (H+H).squared();                    // I = (2 * H)^2
    bls48_Fq8 J = H * I;                              // J = H * I
    bls48_Fq8 r = S2_minus_S1 + S2_minus_S1;          // r = 2 * (S2-S1)
    bls48_Fq8 V = U1 * I;                             // V = U1 * I
    bls48_Fq8 X3 = r.squared() - J - (V+V);           // X3 = r^2 - J - 2 * V
    bls48_Fq8 S1_J = S1 * J;
    bls48_Fq8 Y3 = r * (V-X3) - (S1_J+S1_J);          // Y3 = r * (V-X3)-2 S1 J
    bls48_Fq8 Z3 = ((this->Z+other.Z).squared()-Z1Z1-Z2Z2) * H; // Z3 = ((Z1+Z2)^2-Z1Z1-Z2Z2) * H

    return bls48_G2(X3, Y3, Z3);
}

bls48_G2 bls48_G2::operator-() const
{
    return bls48_G2(this->X, -(this->Y), this->Z);
}


bls48_G2 bls48_G2::operator-(const bls48_G2 &other) const
{
    return (*this) + (-other);
}

bls48_G2 bls48_G2::add(const bls48_G2 &other) const
{
    // handle special cases having to do with O
    if (this->is_zero())
    {
        return other;
    }

    if (other.is_zero())
    {
        return *this;
    }

    // no need to handle points of order 2,4
    // (they cannot exist in a prime-order subgroup)

    // handle double case
    if (this->operator==(other))
    {
        return this->dbl();
    }

#ifdef PROFILE_OP_COUNTS
    this->add_cnt++;
#endif
    // NOTE: does not handle O and pts of order 2,4
    // http://www.hyperelliptic.org/EFD/g1p/auto-shortw-projective.html#addition-add-1998-cmo-2

    bls48_Fq8 Z1Z1 = (this->Z).squared();             // Z1Z1 = Z1^2
    bls48_Fq8 Z2Z2 = (other.Z).squared();             // Z2Z2 = Z2^2
    bls48_Fq8 U1 = (this->X) * Z2Z2;                  // U1 = X1 * Z2Z2
    bls48_Fq8 U2 = (other.X) * Z1Z1;                  // U2 = X2 * Z1Z1
    bls48_Fq8 S1 = (this->Y) * (other.Z) * Z2Z2;      // S1 = Y1 * Z2 * Z2Z2
    bls48_Fq8 S2 = (other.Y) * (this->Z) * Z1Z1;      // S2 = Y2 * Z1 * Z1Z1
    bls48_Fq8 H = U2 - U1;                            // H = U2-U1
    bls48_Fq8 S2_minus_S1 = S2-S1;
    bls48_Fq8 I = (H+H).squared();                    // I = (2 * H)^2
    bls48_Fq8 J = H * I;                              // J = H * I
    bls48_Fq8 r = S2_minus_S1 + S2_minus_S1;          // r = 2 * (S2-S1)
    bls48_Fq8 V = U1 * I;                             // V = U1 * I
    bls48_Fq8 X3 = r.squared() - J - (V+V);           // X3 = r^2 - J - 2 * V
    bls48_Fq8 S1_J = S1 * J;
    bls48_Fq8 Y3 = r * (V-X3) - (S1_J+S1_J);          // Y3 = r * (V-X3)-2 S1 J
    bls48_Fq8 Z3 = ((this->Z+other.Z).squared()-Z1Z1-Z2Z2) * H; // Z3 = ((Z1+Z2)^2-Z1Z1-Z2Z2) * H

    return bls48_G2(X3, Y3, Z3);
}

bls48_G2 bls48_G2::mixed_add(const bls48_G2 &other) const
{
#ifdef DEBUG
    assert(other.is_special());
#endif

    // handle special cases having to do with O
    if (this->is_zero())
    {
        return other;
    }

    if (other.is_zero())
    {
        return *this;
    }

    // no need to handle points of order 2,4
    // (they cannot exist in a prime-order subgroup)

    // check for doubling case

    // using Jacobian coordinates so:
    // (X1:Y1:Z1) = (X2:Y2:Z2)
    // iff
    // X1/Z1^2 == X2/Z2^2 and Y1/Z1^3 == Y2/Z2^3
    // iff
    // X1 * Z2^2 == X2 * Z1^2 and Y1 * Z2^3 == Y2 * Z1^3

    // we know that Z2 = 1

    const bls48_Fq8 Z1Z1 = (this->Z).squared();

    const bls48_Fq8 &U1 = this->X;
    const bls48_Fq8 U2 = other.X * Z1Z1;

    const bls48_Fq8 Z1_cubed = (this->Z) * Z1Z1;

    const bls48_Fq8 &S1 = (this->Y);                // S1 = Y1 * Z2 * Z2Z2
    const bls48_Fq8 S2 = (other.Y) * Z1_cubed;      // S2 = Y2 * Z1 * Z1Z1

    if (U1 == U2 && S1 == S2)
    {
        // dbl case; nothing of above can be reused
        return this->dbl();
    }

#ifdef PROFILE_OP_COUNTS
    this->add_cnt++;
#endif

    // NOTE: does not handle O and pts of order 2,4
    // http://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#addition-madd-2007-bl
    bls48_Fq8 H = U2-(this->X);                         // H = U2-X1
    bls48_Fq8 HH = H.squared() ;                        // HH = H&2
    bls48_Fq8 I = HH+HH;                                // I = 4*HH
    I = I + I;
    bls48_Fq8 J = H*I;                                  // J = H*I
    bls48_Fq8 r = S2-(this->Y);                         // r = 2*(S2-Y1)
    r = r + r;
    bls48_Fq8 V = (this->X) * I ;                       // V = X1*I
    bls48_Fq8 X3 = r.squared()-J-V-V;                   // X3 = r^2-J-2*V
    bls48_Fq8 Y3 = (this->Y)*J;                         // Y3 = r*(V-X3)-2*Y1*J
    Y3 = r*(V-X3) - Y3 - Y3;
    bls48_Fq8 Z3 = ((this->Z)+H).squared() - Z1Z1 - HH; // Z3 = (Z1+H)^2-Z1Z1-HH

    return bls48_G2(X3, Y3, Z3);
}

bls48_G2 bls48_G2::dbl() const
{
#ifdef PROFILE_OP_COUNTS
    this->dbl_cnt++;
#endif
    // handle point at infinity
    if (this->is_zero())
    {
        return (*this);
    }

    // NOTE: does not handle O and pts of order 2,4
    // http://www.hyperelliptic.org/EFD/g1p/auto-shortw-projective.html#doubling-dbl-2007-bl

    bls48_Fq8 A = (this->X).squared();         // A = X1^2
    bls48_Fq8 B = (this->Y).squared();        // B = Y1^2
    bls48_Fq8 C = B.squared();                // C = B^2
    bls48_Fq8 D = (this->X + B).squared() - A - C;
    D = D+D;                        // D = 2 * ((X1 + B)^2 - A - C)
    bls48_Fq8 E = A + A + A;                  // E = 3 * A
    bls48_Fq8 F = E.squared();                // F = E^2
    bls48_Fq8 X3 = F - (D+D);                 // X3 = F - 2 D
    bls48_Fq8 eightC = C+C;
    eightC = eightC + eightC;
    eightC = eightC + eightC;
    bls48_Fq8 Y3 = E * (D - X3) - eightC;     // Y3 = E * (D - X3) - 8 * C
    bls48_Fq8 Y1Z1 = (this->Y)*(this->Z);
    bls48_Fq8 Z3 = Y1Z1 + Y1Z1;               // Z3 = 2 * Y1 * Z1

    return bls48_G2(X3, Y3, Z3);
}
/*
bls48_G2 bls48_G2::mul_by_q() const
{
    return bls48_G2(bls48_twist_mul_by_q_X * (this->X).Frobenius_map(1),
                      bls48_twist_mul_by_q_Y * (this->Y).Frobenius_map(1),
                      (this->Z).Frobenius_map(1));
}
*/

bool bls48_G2::is_well_formed() const
{
    if (this->is_zero())
    {
        return true;
    }
    else
    {
        /*
          y^2 = x^3 + b

          We are using Jacobian coordinates, so equation we need to check is actually

          (y/z^3)^2 = (x/z^2)^3 + b
          y^2 / z^6 = x^3 / z^6 + b
          y^2 = x^3 + b z^6
        */
        bls48_Fq8 X2 = this->X.squared();
        bls48_Fq8 Y2 = this->Y.squared();
        bls48_Fq8 Z2 = this->Z.squared();

        bls48_Fq8 X3 = this->X * X2;
        bls48_Fq8 Z3 = this->Z * Z2;
        bls48_Fq8 Z6 = Z3.squared();

        return (Y2 == X3 + bls48_twist_coeff_b * Z6);
    }
}

bls48_G2 bls48_G2::zero()
{
    return G2_zero;
}

bls48_G2 bls48_G2::one()
{
    return G2_one;
}

bls48_G2 bls48_G2::random_element()
{
    return (bls48_Fr::random_element().as_bigint()) * G2_one;
}

std::ostream& operator<<(std::ostream &out, const bls48_G2 &g)
{
    bls48_G2 copy(g);
    copy.to_affine_coordinates();
    out << (copy.is_zero() ? 1 : 0) << OUTPUT_SEPARATOR;
#ifdef NO_PT_COMPRESSION
    out << copy.X << OUTPUT_SEPARATOR << copy.Y;
#else
    /* storing LSB of Y */
    out << copy.X << OUTPUT_SEPARATOR << (copy.Y.c0.c0.c0.as_bigint().data[0] & 1);
#endif

    return out;
}

std::istream& operator>>(std::istream &in, bls48_G2 &g)
{
    char is_zero;
    bls48_Fq8 tX, tY;

#ifdef NO_PT_COMPRESSION
    in >> is_zero >> tX >> tY;
    is_zero -= '0';
#else
    in.read((char*)&is_zero, 1); // this reads is_zero;
    is_zero -= '0';
    consume_OUTPUT_SEPARATOR(in);

    unsigned char Y_lsb;
    in >> tX;
    consume_OUTPUT_SEPARATOR(in);
    in.read((char*)&Y_lsb, 1);
    Y_lsb -= '0';

    // y = +/- sqrt(x^3 + b)
    if (!is_zero)
    {
        bls48_Fq8 tX2 = tX.squared();
        bls48_Fq8 tY2 = tX2 * tX + bls48_twist_coeff_b;
        tY = tY2.sqrt();

        if ((tY.c0.c0.c0.as_bigint().data[0] & 1) != Y_lsb)
        {
            tY = -tY;
        }
    }
#endif
    // using projective coordinates
    if (!is_zero)
    {
        g.X = tX;
        g.Y = tY;
        g.Z = bls48_Fq8::one();
    }
    else
    {
        g = bls48_G2::zero();
    }

    return in;
}

template<>
void batch_to_special_all_non_zeros<bls48_G2>(std::vector<bls48_G2> &vec)
{
    std::vector<bls48_Fq8> Z_vec;
    Z_vec.reserve(vec.size());

    for (auto &el: vec)
    {
        Z_vec.emplace_back(el.Z);
    }
    batch_invert<bls48_Fq8>(Z_vec);

    const bls48_Fq8 one = bls48_Fq8::one();

    for (size_t i = 0; i < vec.size(); ++i)
    {
        bls48_Fq8 Z2 = Z_vec[i].squared();
        bls48_Fq8 Z3 = Z_vec[i] * Z2;

        vec[i].X = vec[i].X * Z2;
        vec[i].Y = vec[i].Y * Z3;
        vec[i].Z = one;
    }
}

} // libsnark
