use std::ops::Mul;

use ark_ec::AffineCurve;
use ark_ff::{BigInteger, Field, PrimeField, ToBytes};
use rayon::{prelude::IndexedParallelIterator, slice::ParallelSliceMut};

#[cfg(feature = "ark-bls12-377")]
pub use ark_bls12_377 as curve;

#[cfg(feature = "ark-bn254")]
pub use ark_bn254 as curve;

use self::curve::{Fq, Fq2, G1Affine, G2Affine};

use rust_rw_device::rw_msm_to_dram as device_g1; //TODO: unify to one crate
use rust_rw_device_g2::rw_msm_to_dram as device_g2; //TODO: unify to one crate

use std::{any::TypeId, io::Cursor};

use std::mem::size_of;

use rayon::iter::ParallelIterator;

fn points_to_bytes<G: AffineCurve>(points: &[G]) -> Vec<u8> {
    let mut buff = Vec::<u8>::new();

    let field_size = size_of::<G::BaseField>();

    if TypeId::of::<G>() == TypeId::of::<G1Affine>() {
        let coord_size = field_size;

        buff.resize(coord_size * 2 * points.len(), 0u8);
        buff.par_chunks_mut(coord_size * 2)
            .zip(points)
            .for_each(|(f, g)| {
                let mut pbuff = Cursor::new(Vec::<u8>::new());
                g.write(&mut pbuff).unwrap();
                let pnt = pbuff.get_ref().to_vec();
                (*f)[..coord_size].copy_from_slice(&pnt[coord_size..coord_size * 2]);
                (*f)[coord_size..coord_size * 2].copy_from_slice(&pnt[..coord_size]);
            });
    } else if TypeId::of::<G>() == TypeId::of::<G2Affine>() {
        let coord_size = field_size / 2;
        buff.resize(coord_size * 4 * points.len(), 0u8);
        buff.par_chunks_mut(coord_size * 4)
            .zip(points)
            .for_each(|(f, g)| {
                let mut pbuff = Cursor::new(Vec::<u8>::new());
                g.write(&mut pbuff).unwrap();
                let pnt = pbuff.get_ref().to_vec();
                (*f)[..coord_size].copy_from_slice(&pnt[coord_size * 3..coord_size * 4]); //point.y.c1
                (*f)[coord_size..coord_size * 2].copy_from_slice(&pnt[coord_size..coord_size * 2]); //point.x.c1
                (*f)[coord_size * 2..coord_size * 3]
                    .copy_from_slice(&pnt[coord_size * 2..coord_size * 3]); //point.y.c0
                (*f)[coord_size * 3..].copy_from_slice(&pnt[..coord_size]); //point.x.c0
            });
    } else {
        panic!("unsupported curve")
    }
    buff
}

fn scalars_to_bytes<G: AffineCurve>(scalars: &[<G::ScalarField as PrimeField>::BigInt]) -> Vec<u8> {
    let scalar_size = size_of::<<G::ScalarField as PrimeField>::BigInt>();
    let mut buff = vec![0u8; scalar_size * scalars.len()];
    buff.par_chunks_mut(scalar_size)
        .zip(scalars)
        .for_each(|(f, g)| {
            g.write(&mut Cursor::new(f)).unwrap();
        });
    buff
}

pub struct HardwareVariableBaseMSM;

impl HardwareVariableBaseMSM {
    pub fn multi_scalar_mul<G: AffineCurve>(
        bases: &[G],
        scalars: &[<G::ScalarField as PrimeField>::BigInt],
    ) -> G::Projective {
        let size = std::cmp::min(bases.len(), scalars.len());
        let scalars = &scalars[..size];
        let bases = &bases[..size];

        let mut bases_filtered = Vec::<G>::new();
        let mut scalars_filtered = Vec::<<G::ScalarField as PrimeField>::BigInt>::new();

        for i in 0..size {
            //filtering zero scalars *AND* zero points is key for Ingo MSM compatibility for current version
            if !scalars[i].is_zero() && !bases[i].is_zero() {
                bases_filtered.push(bases[i]);
                scalars_filtered.push(scalars[i]);
            }
        }
        let points_bytes = points_to_bytes(&bases_filtered);
        let scalars_bytes = scalars_to_bytes::<G>(&scalars_filtered);
        let mut buff: Vec<u8>;
        let scalar_size = size_of::<G::ScalarField>();
        let size = scalars_bytes.len() / scalar_size;

        if TypeId::of::<G>() == TypeId::of::<G1Affine>() {
            let (result, _, _) = device_g1::msm_calc(&points_bytes, &scalars_bytes, size);
            let proj_x_field = Fq::from_random_bytes(&result[0]).unwrap();
            let proj_y_field = Fq::from_random_bytes(&result[1]).unwrap();
            let proj_z_field = Fq::from_random_bytes(&result[2]).unwrap();
            let aff_x = proj_x_field.mul(proj_z_field.inverse().unwrap());
            let aff_y = proj_y_field.mul(proj_z_field.inverse().unwrap());
            buff = Vec::<u8>::with_capacity(size_of::<G1Affine>());
            aff_x.write(&mut buff).unwrap();
            aff_y.write(&mut buff).unwrap();
        } else if TypeId::of::<G>() == TypeId::of::<G2Affine>() {
            let (result, _, _) = device_g2::msm_calc(&points_bytes, &scalars_bytes, size);
            let proj_x_field =
                Fq2::from_random_bytes(&[result[5].to_vec(), result[2].to_vec()].concat()).unwrap();
            let proj_y_field =
                Fq2::from_random_bytes(&[result[4].to_vec(), result[1].to_vec()].concat()).unwrap();
            let proj_z_field =
                Fq2::from_random_bytes(&[result[3].to_vec(), result[0].to_vec()].concat()).unwrap();

            let aff_x = proj_x_field.mul(proj_z_field.inverse().unwrap());
            let aff_y = proj_y_field.mul(proj_z_field.inverse().unwrap());
            buff = Vec::<u8>::with_capacity(size_of::<G2Affine>());
            aff_x.write(&mut buff).unwrap();
            aff_y.write(&mut buff).unwrap();
        } else {
            todo!("unsupported curve type")
        }

        buff.push(0);
        G::read(buff.as_slice()).unwrap().into_projective()
    }
}

#[cfg(test)]
mod test {
    use ark_ec::{msm::VariableBaseMSM, AffineCurve, ProjectiveCurve};
    use ark_ff::{PrimeField, UniformRand};
    use ark_std::test_rng;

    use std::str::FromStr;

    use super::curve::*;

    use super::HardwareVariableBaseMSM;

    fn generate_points_scalars<G: AffineCurve>(
        len: usize,
    ) -> (Vec<G>, Vec<<G::ScalarField as PrimeField>::BigInt>) {
        let rng = &mut test_rng();
        (
            <G::Projective as ProjectiveCurve>::batch_normalization_into_affine(
                &(0..len)
                    .map(|_| G::Projective::rand(rng))
                    .collect::<Vec<_>>(),
            ),
            (0..len)
                .map(|_| G::ScalarField::rand(rng).into_repr())
                .collect::<Vec<_>>(),
        )
    }

    fn msm_correctness<G: AffineCurve>() {
        let test_npow = std::env::var("TEST_NPOW").unwrap_or("11".to_string());
        let n_points = i32::from_str(&test_npow).unwrap();

        let len = 1 << n_points;
        let (points, scalars) = generate_points_scalars::<G>(len);

        let msm_ark_projective = VariableBaseMSM::multi_scalar_mul(&points.to_vec(), &scalars);

        let ret = HardwareVariableBaseMSM::multi_scalar_mul(&points, &scalars);
        assert_eq!(ret, msm_ark_projective);
    }

    #[test]
    fn msm_correctness_g1() {
        msm_correctness::<G1Affine>();
    }

    #[test]
    fn msm_correctness_g2() {
        msm_correctness::<G2Affine>();
    }
}
