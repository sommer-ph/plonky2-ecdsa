#![allow(clippy::needless_range_loop)]
#![cfg_attr(not(test), no_std)]

extern crate alloc;

pub mod curve;
pub mod field;
pub mod gadgets;
#[cfg(feature = "static-issuer")]
pub mod ecdsa;
#[cfg(feature = "static-issuer")]
pub use ecdsa::verify_static_pk::add_static_pk_ecdsa_verify_constraints;
