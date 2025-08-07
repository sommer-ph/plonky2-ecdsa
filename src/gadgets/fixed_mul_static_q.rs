#![cfg(feature = "static-issuer")]

use alloc::vec::Vec;
use num::BigUint;
use plonky2::field::extension::Extendable;
use plonky2::hash::hash_types::RichField;
use plonky2::hash::keccak::KeccakHash;
use plonky2::plonk::circuit_builder::CircuitBuilder;
use plonky2::plonk::config::{GenericHashOut, Hasher};
use plonky2::field::types::Field;

use crate::curve::curve_types::{Curve, CurveScalar};
use crate::curve::p256::P256;
use crate::field::p256_scalar::P256Scalar;
use crate::gadgets::curve::{AffinePointTarget, CircuitBuilderCurve};
use crate::gadgets::lookup_tables::Q_TABLE;
use crate::gadgets::nonnative::NonNativeTarget;
use crate::gadgets::split_nonnative::CircuitBuilderSplit;

/// Number of bits per window.
const WINDOW_BITS: usize = 4;
/// Number of windows used for the fixed-base multiplication.
const NUM_WINDOWS: usize = 16;

/// Multiply the fixed issuer public key `Q` by a scalar using lookup tables.
pub(crate) fn fixed_mul_static_q<F: RichField + Extendable<D>, const D: usize>(
    builder: &mut CircuitBuilder<F, D>,
    scalar: &NonNativeTarget<P256Scalar>,
) -> AffinePointTarget<P256> {
    // Decompose the scalar into 4-bit windows.
    let limbs = builder.split_nonnative_to_4_bit_limbs(scalar);
    let zero = builder.zero();

    // Random starting point to avoid exceptional cases in incomplete addition.
    let hash_0 = KeccakHash::<32>::hash_no_pad(&[F::ZERO]);
    let hash_0_scalar = P256Scalar::from_noncanonical_biguint(BigUint::from_bytes_le(
        &GenericHashOut::<F>::to_bytes(&hash_0),
    ));
    let mut cur = CurveScalar(hash_0_scalar) * P256::GENERATOR_PROJECTIVE;
    let starting_point = cur;
    for _ in 0..NUM_WINDOWS * WINDOW_BITS {
        cur = cur.double();
    }
    let starting_point_multiplied = cur;

    // Initialize result with the random starting point.
    let mut result = builder.constant_affine_point(starting_point.to_affine());

    for i in (0..NUM_WINDOWS).rev() {
        result = builder.curve_repeated_double(&result, WINDOW_BITS);
        let limb = limbs[i];
        // Prepare lookup table for this window. Replace the identity with a dummy point.
        let mut points: Vec<_> = Q_TABLE[i]
            .iter()
            .map(|p| if p.zero { Q_TABLE[i][1] } else { *p })
            .collect();
        points[0] = points[1];
        let to_add = builder.random_access_affine_point(limb, &points);
        let is_zero = builder.is_equal(limb, zero);
        let should_add = builder.not(is_zero);
        result = builder.curve_conditional_add(&result, &to_add, should_add);
    }

    // Subtract the starting point multiplied by 2^{WINDOW_BITS * NUM_WINDOWS}.
    let to_subtract = builder.constant_affine_point(starting_point_multiplied.to_affine());
    let to_add = builder.curve_neg(&to_subtract);
    builder.curve_add(&result, &to_add)
}
