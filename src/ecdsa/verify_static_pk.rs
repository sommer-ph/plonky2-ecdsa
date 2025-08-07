#![cfg(feature = "static-issuer")]

//! ECDSA verification gadget using a fixed issuer public key `Q`.

use plonky2::field::extension::Extendable;
use plonky2::field::types::Field;
use plonky2::hash::hash_types::RichField;
use plonky2::plonk::circuit_builder::CircuitBuilder;

use crate::curve::curve_types::Curve;
use crate::curve::p256::P256;
use crate::field::p256_scalar::P256Scalar;
use crate::gadgets::curve::CircuitBuilderCurve;
use crate::gadgets::curve_fixed_base::fixed_base_curve_mul_circuit;
use crate::gadgets::fixed_mul_static_q::fixed_mul_static_q;
use crate::gadgets::nonnative::{CircuitBuilderNonNative, NonNativeTarget};
use crate::gadgets::biguint::CircuitBuilderBiguint;
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

    // Range check: ensure 1 <= r < n and 1 <= s < n
    let n_biguint = builder.constant_biguint(&P256Scalar::order());
    
    let r_biguint = builder.nonnative_to_canonical_biguint(&r);
    let s_biguint = builder.nonnative_to_canonical_biguint(&s);
    
    // Check r < n
    let r_lt_n = builder.cmp_biguint(&r_biguint, &n_biguint);
    builder.assert_one(r_lt_n.target);
    
    // Check s < n  
    let s_lt_n = builder.cmp_biguint(&s_biguint, &n_biguint);
    builder.assert_one(s_lt_n.target);
    
    // Check r != 0 and s != 0 by using BigUint comparison
    let zero_biguint = builder.zero_biguint();
    
    // Check r != 0: cmp_biguint(a, 0) returns true iff a == 0 (for BigUint, a <= 0 means a == 0)
    let r_eq_zero = builder.cmp_biguint(&r_biguint, &zero_biguint);  // true iff r == 0
    let r_neq_zero = builder.not(r_eq_zero);
    builder.assert_one(r_neq_zero.target);
    
    // Check s != 0: same logic
    let s_eq_zero = builder.cmp_biguint(&s_biguint, &zero_biguint);  // true iff s == 0
    let s_neq_zero = builder.not(s_eq_zero);
    builder.assert_one(s_neq_zero.target);

    let c = builder.inv_nonnative(&s);
    let u1 = builder.mul_nonnative(&msg, &c);
    let u2 = builder.mul_nonnative(&r, &c);

    let g_mul = fixed_base_curve_mul_circuit(builder, P256::GENERATOR_AFFINE, &u1);
    let q_mul = fixed_mul_static_q(builder, &u2);
    let sum = builder.curve_add(&g_mul, &q_mul);

    // Convert sum.x (P256Base field element) to BigUint
    let x_biguint = builder.nonnative_to_canonical_biguint(&sum.x);
    
    // Get the scalar field order n as a constant
    let n_biguint = builder.constant_biguint(&P256Scalar::order());
    
    // Compute x mod n
    let x_mod_n_biguint = builder.rem_biguint(&x_biguint, &n_biguint);
    
    // Convert r to BigUint for comparison
    let r_biguint = builder.nonnative_to_canonical_biguint(&r);

    // Now we can compare r with x(sum) mod n
    builder.connect_biguint(&r_biguint, &x_mod_n_biguint);
}

#[cfg(test)]
mod tests {
    use super::*;
    use num::BigUint;
    use plonky2::iop::witness::PartialWitness;
    use plonky2::plonk::circuit_builder::CircuitBuilder;
    use plonky2::plonk::circuit_data::CircuitConfig;
    use plonky2::plonk::config::{GenericConfig, PoseidonGoldilocksConfig};
    use crate::gadgets::nonnative::CircuitBuilderNonNative;
    use crate::gadgets::lookup_tables::Q_TABLE;

    #[test]
    fn test_precompute_endianness() {
        // Test that our precomputed Q matches the expected value
        const D: usize = 2;
        type C = PoseidonGoldilocksConfig;
        type F = <C as GenericConfig<D>>::F;
        
        let config = CircuitConfig::standard_recursion_config();
        let pw = PartialWitness::new();
        let mut builder = CircuitBuilder::<F, D>::new(config);

        // Load Q.x from the first precompute entry (should be Q itself: 1 * Q)
        let q_entry = Q_TABLE[0][1]; // table[0][1] = 1 * (2^0) * Q = Q
        let q_x_target = builder.constant_nonnative(q_entry.x);
        
        // Expected x coordinate as decimal string
        let expected_x_dec = "66432692286261411630769223098970693805397596870633670159153355502222145619968";
        let expected_x = BigUint::parse_bytes(expected_x_dec.as_bytes(), 10).unwrap();
        
        // Convert the circuit target to BigUint and compare
        let q_x_biguint = builder.nonnative_to_canonical_biguint(&q_x_target);
        let expected_x_target = builder.constant_biguint(&expected_x);
        builder.connect_biguint(&q_x_biguint, &expected_x_target);

        let data = builder.build::<C>();
        let proof = data.prove(pw).unwrap();
        data.verify(proof).expect("Precompute endianness test failed - Q.x does not match expected value");
    }

    #[test] 
    fn test_precompute_window_1() {
        // Test that table[1][1] = 16*Q (verify windowing base multiplication)
        const D: usize = 2;
        type C = PoseidonGoldilocksConfig;
        type F = <C as GenericConfig<D>>::F;
        
        let mut config = CircuitConfig::standard_recursion_config();
        config.num_wires = 136;
        let pw = PartialWitness::new();
        let mut builder = CircuitBuilder::<F, D>::new(config);

        // Get both entries  
        let q_entry = Q_TABLE[0][1]; // Q
        let q16_entry = Q_TABLE[1][1]; // Should be 16*Q
        
        // Build targets
        let q_target = builder.constant_affine_point(q_entry);
        let q16_from_table = builder.constant_affine_point(q16_entry);
        
        // Compute 16*Q using repeated doubling
        let q16_computed = builder.curve_repeated_double(&q_target, 4); // 2^4 = 16
        
        // They should be equal
        builder.connect_affine_point(&q16_from_table, &q16_computed);

        let data = builder.build::<C>();
        let proof = data.prove(pw).unwrap();
        data.verify(proof).expect("Window 1 test failed - table[1][1] != 16*Q");
    }

    #[test]
    fn test_fixed_mul_static_q() {
        // Test fixed_mul_static_q against host p256 reference
        use p256::{AffinePoint as P256AffinePoint, ProjectivePoint, Scalar};
        use p256::elliptic_curve::sec1::{EncodedPoint, FromEncodedPoint, ToEncodedPoint};
        use p256::elliptic_curve::PrimeField;
        
        const D: usize = 2;
        type C = PoseidonGoldilocksConfig;
        type F = <C as GenericConfig<D>>::F;
        
        let config = CircuitConfig::standard_ecc_config();
        let pw = PartialWitness::new();
        let mut builder = CircuitBuilder::<F, D>::new(config);

        // Test scalars: 1, 2, and n-1 (where n is the curve order)
        // P-256 curve order n = 0xFFFFFFFF00000000FFFFFFFFFFFFFFFFBCE6FAADA7179E84F3B9CAC2FC632551
        let n_minus_1_bytes = [
            0xFF, 0xFF, 0xFF, 0xFF, 0x00, 0x00, 0x00, 0x00,
            0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
            0xBC, 0xE6, 0xFA, 0xAD, 0xA7, 0x17, 0x9E, 0x84,
            0xF3, 0xB9, 0xCA, 0xC2, 0xFC, 0x63, 0x25, 0x50,  // n-1
        ];
        let test_scalars = [
            BigUint::from(1u32),
            BigUint::from(2u32),
            BigUint::from_bytes_be(&n_minus_1_bytes),  // n-1 (largest valid scalar)
        ];

        // Q coordinates from build.rs constants
        let q_x_dec = "66432692286261411630769223098970693805397596870633670159153355502222145619968";
        let q_y_dec = "63182586149833488067701290985084360701345487374231728189741684364091950142361";
        let q_x_big = BigUint::parse_bytes(q_x_dec.as_bytes(), 10).unwrap();
        let q_y_big = BigUint::parse_bytes(q_y_dec.as_bytes(), 10).unwrap();
        let q_x_bytes = q_x_big.to_bytes_be();
        let q_y_bytes = q_y_big.to_bytes_be();
        let mut q_x_arr = [0u8; 32];
        let mut q_y_arr = [0u8; 32];
        q_x_arr[32 - q_x_bytes.len()..].copy_from_slice(&q_x_bytes);
        q_y_arr[32 - q_y_bytes.len()..].copy_from_slice(&q_y_bytes);
        let q_encoded = EncodedPoint::<p256::NistP256>::from_affine_coordinates(&q_x_arr.into(), &q_y_arr.into(), false);
        let q_host = P256AffinePoint::from_encoded_point(&q_encoded).unwrap();

        for scalar_big in test_scalars.iter() {
            // Host computation: scalar * Q
            let scalar_bytes = scalar_big.to_bytes_be();
            let mut scalar_arr = [0u8; 32];
            let offset = 32 - scalar_bytes.len();
            scalar_arr[offset..].copy_from_slice(&scalar_bytes);
            let scalar_host = Scalar::from_repr(scalar_arr.into()).unwrap();
            let expected_point = (ProjectivePoint::from(q_host) * scalar_host).to_affine();

            // Circuit computation using fixed_mul_static_q
            let scalar_target = builder.constant_biguint(scalar_big);
            let scalar_nonnative = builder.biguint_to_nonnative(&scalar_target);
            let result_point = fixed_mul_static_q(&mut builder, &scalar_nonnative);

            // Expected result as circuit constants
            let expected_encoded = expected_point.to_encoded_point(false);
            let expected_x_bytes = expected_encoded.x().unwrap();
            let expected_y_bytes = expected_encoded.y().unwrap();
            let expected_x_big = BigUint::from_bytes_be(expected_x_bytes);
            let expected_y_big = BigUint::from_bytes_be(expected_y_bytes);
            let expected_x_target = builder.constant_biguint(&expected_x_big);
            let expected_y_target = builder.constant_biguint(&expected_y_big);
            let expected_x_nonnative = builder.biguint_to_nonnative(&expected_x_target);
            let expected_y_nonnative = builder.biguint_to_nonnative(&expected_y_target);

            // Connect circuit result to expected
            builder.connect_nonnative(&result_point.x, &expected_x_nonnative);
            builder.connect_nonnative(&result_point.y, &expected_y_nonnative);
        }

        let data = builder.build::<C>();
        let proof = data.prove(pw).unwrap();
        data.verify(proof).expect("fixed_mul_static_q test failed - result != host reference");
    }
}
