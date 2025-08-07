#![cfg(feature = "static-issuer")]

//! ECDSA verification gadget using a fixed issuer public key `Q`.

use core::marker::PhantomData;

use plonky2::field::extension::Extendable;
use plonky2::hash::hash_types::RichField;
use plonky2::plonk::circuit_builder::CircuitBuilder;

use crate::curve::curve_types::Curve;
use crate::curve::p256::P256;
use crate::field::p256_scalar::P256Scalar;
use crate::gadgets::curve::CircuitBuilderCurve;
use crate::gadgets::curve_fixed_base::fixed_base_curve_mul_circuit;
use crate::gadgets::fixed_mul_static_q::fixed_mul_static_q;
use crate::gadgets::nonnative::{CircuitBuilderNonNative, NonNativeTarget};
use crate::gadgets::ecdsa::ECDSASignatureTarget;

/// Verify an ECDSA signature with the fixed issuer public key `Q`.
///
/// * `msg` - The hashed message as a non-native field element.
/// * `sig` - The signature `(r, s)`.
pub fn add_static_pk_ecdsa_verify_constraints<F: RichField + Extendable<D>, const D: usize>(
    builder: &mut CircuitBuilder<F, D>,
    msg: NonNativeTarget<P256Scalar>,
    sig: ECDSASignatureTarget<P256>,
) {
    let ECDSASignatureTarget { r, s } = sig;

    let c = builder.inv_nonnative(&s);
    let u1 = builder.mul_nonnative(&msg, &c);
    let u2 = builder.mul_nonnative(&r, &c);

    let g_mul = fixed_base_curve_mul_circuit(builder, P256::GENERATOR_AFFINE, &u1);
    let q_mul = fixed_mul_static_q(builder, &u2);
    let sum = builder.curve_add(&g_mul, &q_mul);

    let x = NonNativeTarget::<P256Scalar> {
        value: sum.x.value,
        _phantom: PhantomData,
    };
    builder.connect_nonnative(&r, &x);
}
