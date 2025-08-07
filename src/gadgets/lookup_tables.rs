#![cfg(feature = "static-issuer")]

//! Precomputed lookup tables for the static issuer public key `Q`.

use crate::curve::curve_types::AffinePoint;
use crate::curve::p256::P256;

/// Lookup table for `Q`, organized as 16 windows of 16 points each.
/// Each window contains the multiples `k * 16^i * Q` for `k = 0..15`.
#[allow(dead_code)]
pub(crate) const Q_TABLE: [[AffinePoint<P256>; 16]; 16] =
    include!(concat!(env!("OUT_DIR"), "/static_q_table.in"));
