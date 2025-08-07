use std::env;
use std::fs::File;
use std::io::Write;
use std::path::Path;

use num_bigint::BigUint;
use p256::elliptic_curve::group::Group;
use p256::elliptic_curve::sec1::{EncodedPoint, FromEncodedPoint, ToEncodedPoint};
use p256::{AffinePoint, ProjectivePoint};

const X_DEC: &str = "66432692286261411630769223098970693805397596870633670159153355502222145619968";
const Y_DEC: &str = "63182586149833488067701290985084360701345487374231728189741684364091950142361";

fn main() {
    // Parse decimal coordinates of Q.
    let x = BigUint::parse_bytes(X_DEC.as_bytes(), 10).expect("invalid x");
    let y = BigUint::parse_bytes(Y_DEC.as_bytes(), 10).expect("invalid y");

    // Convert to 32-byte big-endian arrays.
    let mut x_bytes = x.to_bytes_be();
    let mut y_bytes = y.to_bytes_be();
    assert!(x_bytes.len() <= 32 && y_bytes.len() <= 32);
    if x_bytes.len() < 32 {
        let mut padded = vec![0u8; 32 - x_bytes.len()];
        padded.extend_from_slice(&x_bytes);
        x_bytes = padded;
    }
    if y_bytes.len() < 32 {
        let mut padded = vec![0u8; 32 - y_bytes.len()];
        padded.extend_from_slice(&y_bytes);
        y_bytes = padded;
    }
    let x_arr: [u8; 32] = x_bytes.try_into().unwrap();
    let y_arr: [u8; 32] = y_bytes.try_into().unwrap();

    let encoded = EncodedPoint::<p256::NistP256>::from_affine_coordinates(&x_arr.into(), &y_arr.into(), false);
    let q_aff = AffinePoint::from_encoded_point(&encoded).expect("invalid point");
    let mut base = ProjectivePoint::from(q_aff);

    // Precompute the lookup table: 64 windows x 16 points each (covers full 256 bits).
    let mut table: Vec<Vec<AffinePoint>> = Vec::new();
    for _ in 0..64 {
        let mut multiples: Vec<AffinePoint> = Vec::new();
        let mut acc = ProjectivePoint::IDENTITY;
        for _ in 0..16 {
            multiples.push(acc.to_affine());
            acc += base;
        }
        table.push(multiples);
        // Move to the next window: multiply base by 16 (i.e., double four times).
        for _ in 0..4 {
            base = base.double();
        }
    }

    let out_dir = env::var("OUT_DIR").unwrap();
    let dest_path = Path::new(&out_dir).join("static_q_table.in");
    let mut f = File::create(dest_path).unwrap();

    writeln!(
        f,
        "{}",
        "["
    ).unwrap();

    for window in table {
        writeln!(f, "    [").unwrap();
        for point in window {
            if bool::from(point.is_identity()) {
                writeln!(f, "        crate::curve::curve_types::AffinePoint::ZERO,").unwrap();
                continue;
            }
            let encoded = point.to_encoded_point(false);
            let x_vals = bytes_to_u64s(encoded.x().unwrap());
            let y_vals = bytes_to_u64s(encoded.y().unwrap());
            writeln!(
                f,
                "        crate::curve::curve_types::AffinePoint {{ x: crate::field::p256_base::P256Base([{:#x}, {:#x}, {:#x}, {:#x}]), y: crate::field::p256_base::P256Base([{:#x}, {:#x}, {:#x}, {:#x}]), zero: false }},",
                x_vals[0], x_vals[1], x_vals[2], x_vals[3],
                y_vals[0], y_vals[1], y_vals[2], y_vals[3]
            ).unwrap();
        }
        writeln!(f, "    ],").unwrap();
    }
    writeln!(f, "]").unwrap();
}

fn bytes_to_u64s(bytes: &[u8]) -> [u64; 4] {
    let mut res = [0u64; 4];
    // bytes are big-endian. Convert to LSB-first limb order for plonky2.
    // res[0] should contain the least significant 64 bits, res[3] the most significant.
    for i in 0..4 {
        let start = 8 * i;
        let mut chunk = [0u8; 8];
        chunk.copy_from_slice(&bytes[start..start + 8]);
        // Since input is big-endian, we need to reverse the order
        res[3 - i] = u64::from_be_bytes(chunk);
    }
    res
}
