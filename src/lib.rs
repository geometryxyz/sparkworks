use ark_ec::{msm::VariableBaseMSM, AffineCurve};
use ark_ff::{BigInteger, PrimeField, ToBytes};

use std::{any::TypeId, io::Cursor};

use ingo_x::curve::*;

use std::mem::size_of;

use ark_ff::Zero;

const BYTE_SIZE_SCALAR: usize = 32;

fn get_formatted_unified_points_from_g1_affine<G: AffineCurve>(points: &[G]) -> (Vec<u8>, usize) {
    /*
       In order to determine point size we can write zero to buffer.
       It writes: (x, y, is_zero_byte) so point size is (buff_size - 1) / 2, or just integer / 2 since 1 will be a reminder
    */
    let mut buff = Cursor::new(Vec::<u8>::new());
    G::zero().write(&mut buff).unwrap();
    let point_size = buff.get_ref().len() / 2;

    let mut points_buffer: Vec<u8> = vec![0; points.len() * 2 * point_size];

    for (i, base) in points.iter().enumerate() {
        // reset buffer in each iteration
        buff.set_position(0);
        base.write(&mut buff).unwrap();

        // write y
        points_buffer[2 * i * point_size..(2 * i + 1) * point_size]
            .copy_from_slice(&buff.get_ref()[point_size..2 * point_size]);
        // write x
        points_buffer[(2 * i + 1) * point_size..(2 * i + 2) * point_size]
            .copy_from_slice(&buff.get_ref()[0..point_size]);
    }

    (points_buffer, buff.get_ref().len())
}

fn get_formatted_unified_scalars_from_bigint<G: AffineCurve>(
    scalars: &[<G::ScalarField as PrimeField>::BigInt],
) -> Vec<u8> {
    let mut scalars_bytes: Vec<u8> = Vec::new();
    for i in 0..scalars.len() {
        let mut bytes_array = scalars[i].to_bytes_le();
        bytes_array.extend(std::iter::repeat(0).take(BYTE_SIZE_SCALAR - bytes_array.len()));
        scalars_bytes.extend(bytes_array);
    }
    scalars_bytes
}

fn points_to_bytes<G: AffineCurve>(points: &[G]) -> Vec<u8> { //TODO: this semi-generic method needs/can to be fully generic?
    let mut buff = Vec::<u8>::new();

    let field_size = size_of::<G::BaseField>();

    if TypeId::of::<G>() == TypeId::of::<G1Affine>() {
        let coord_size = field_size;
        points.into_iter().for_each(|f| {
            let mut pbuff = Cursor::new(Vec::<u8>::new());
            f.write(&mut pbuff).unwrap();
            let pnt = pbuff.get_ref().to_vec();
            buff.extend_from_slice(&pnt[coord_size..coord_size * 2]);
            buff.extend_from_slice(&pnt[..coord_size]);
        });
    } else if TypeId::of::<G>() == TypeId::of::<G2Affine>() {
        let coord_size = field_size / 2;
        points.into_iter().for_each(|f| {
            //point.y.c1.into_repr().into(),
            //point.x.c1.into_repr().into(),
            //point.y.c0.into_repr().into(),
            //point.x.c0.into_repr().into(),
            let mut pbuff = Cursor::new(Vec::<u8>::new());
            f.write(&mut pbuff).unwrap();
            let pnt = pbuff.get_ref().to_vec();
            buff.extend_from_slice(&pnt[coord_size * 3..coord_size * 4]);
            buff.extend_from_slice(&pnt[coord_size..coord_size * 2]);
            buff.extend_from_slice(&pnt[coord_size * 2..coord_size * 3]);
            buff.extend_from_slice(&pnt[..coord_size]);
        });
    } else {
        panic!("unsupported curve")
    }
    buff
}

fn scalars_to_bytes<G: AffineCurve>(scalars: &[<G::ScalarField as PrimeField>::BigInt]) -> Vec<u8> { //this is generic
    let mut buff = Cursor::new(Vec::<u8>::new());
    scalars.into_iter().for_each(|f| {
        f.write(&mut buff).unwrap();
    });
    buff.get_ref().to_vec()
}

pub struct FpgaVariableBaseMSM;

impl FpgaVariableBaseMSM {
    pub fn multi_scalar_mul<G: AffineCurve>(
        bases: &[G],
        scalars: &[<G::ScalarField as PrimeField>::BigInt],
    ) -> G::Projective {

        let size = std::cmp::min(bases.len(), scalars.len());
        let scalars = &scalars[..size];
        let bases = &bases[..size];

        let mut bases_filtered = Vec::<G>::new();
        let mut scalars_filtered = Vec::<<G::ScalarField as PrimeField>::BigInt>::new();

        for i in 0..size { //filtering zero scalars *AND* zero points is key for Ingo MSM compatibility for current version
            if !scalars[i].is_zero() && !bases[i].is_zero() {
                bases_filtered.push(bases[i]);
                scalars_filtered.push(scalars[i]);
            }
        }

        let points_bytes = points_to_bytes(&bases_filtered);
        let scalars_bytes = scalars_to_bytes::<G>(&scalars_filtered);
        let ret = ingo_x::msm_cloud_generic::<G>(&points_bytes, &scalars_bytes);

        ret.0
    }
}

#[cfg(test)]
mod test {
    use crate::{get_formatted_unified_scalars_from_bigint, points_to_bytes, scalars_to_bytes};

    use super::get_formatted_unified_points_from_g1_affine;
    use ark_bn254::{G1Affine, G2Affine};
    use ark_ec::{AffineCurve, ProjectiveCurve};
    use ark_ff::{BigInteger, BigInteger256, PrimeField, UniformRand, Zero};
    use ark_std::{rand::Rng, test_rng};
    use num_bigint::BigUint;

    use std::{ops::Add, str::FromStr};

    use ingo_x::curve::*;

    const BYTE_SIZE_POINT_COORD: usize = 32; // for BLS

    // ingonyama's implementation for asserting equality
    fn get_formatted_unified_points_from_biguint(points: &Vec<BigUint>) -> Vec<u8> {
        let mut points_bytes: Vec<u8> = Vec::new();
        for i in 0..points.len() {
            let mut bytes_array = points[i].to_bytes_le();
            bytes_array
                .extend(std::iter::repeat(0).take(BYTE_SIZE_POINT_COORD - bytes_array.len()));
            points_bytes.extend(bytes_array);
        }
        points_bytes
    }

    fn generate_points_scalars<G: AffineCurve, R: Rng>(len: usize, rng: &mut R) -> Vec<G> {
        <G::Projective as ProjectiveCurve>::batch_normalization_into_affine(
            &(0..len)
                .map(|_| G::Projective::rand(rng))
                .collect::<Vec<_>>(),
        )
    }

    #[test]
    fn test_affine_to_bytes() {
        let mut rng = test_rng();
        let len = 100;

        let points: Vec<G1Affine> = generate_points_scalars(len, &mut rng);

        let points_as_big_int = points
            .iter()
            .map(|point| [point.y.into_repr().into(), point.x.into_repr().into()])
            .flatten()
            .collect::<Vec<BigUint>>();

        let point_bytes_biguint = get_formatted_unified_points_from_biguint(&points_as_big_int);
        let (point_bytes_affine, _) = get_formatted_unified_points_from_g1_affine(&points);

        assert_eq!(point_bytes_biguint, point_bytes_affine);
    }

    #[test]
    fn test_ark_to_bytes() {
        let len = 1;

        let (points_g1, scalars) = ingo_x::util::generate_points_scalars::<G1Affine>(len);
        let (points_g2, _) = ingo_x::util::generate_points_scalars::<G2Affine>(len);

        let scalars = &scalars
            .to_vec()
            .into_iter()
            .map(|s| BigInteger256::try_from(s).unwrap())
            .collect::<Vec<BigInteger256>>(); //this is safe but slow conversion

        let point_g1_bytes = points_to_bytes(&points_g1);
        //let point_g2_bytes = points_to_bytes(&points_g2);
        let scalar_bytes = scalars_to_bytes::<G1Affine>(&scalars);

        let b = get_formatted_unified_points_from_g1_affine(&points_g1);
        let scalars = get_formatted_unified_scalars_from_bigint::<G1Affine>(&scalars);

        assert_eq!(scalar_bytes, scalars);
        assert!(
            point_g1_bytes == b.0,
            "1:::{:02X?}\n2:::{:02X?}",
            point_g1_bytes,
            b.0
        );
        // assert!(
        //     point_g2_bytes == b.0,
        //     "1:::{:02X?}\n2:::{:02X?}",
        //     point_g2_bytes,
        //     b.0
        // );
    }

    #[test]
    pub fn msm_correctness_g1() {
        let test_npow = std::env::var("TEST_NPOW").unwrap_or("11".to_string());
        let n_points = i32::from_str(&test_npow).unwrap();

        //TODO: conversion of inputs/outputs can be much much simplified as done for Sppark GPU and Ingo FPGA MSM

        let len = 1 << n_points;
        let (points, scalars) = ingo_x::util::generate_points_scalars::<G1Affine>(len);

        let msm_ark_projective = ingo_x::msm_ark(
            &points.to_vec(),
            &scalars
                .to_vec()
                .into_iter()
                .map(|s| BigInteger256::try_from(s).unwrap())
                .collect::<Vec<BigInteger256>>(), //this is safe but slow conversion
        );

        let mut msm_result_cpu_ingo_ref = G1Projective::zero(); //TODO: same as G1Affine::prime_subgroup_generator().mul(0);
        let mut msm_result_cpu_ref1 = G1Projective::zero();
        for i in 0..len {
            msm_result_cpu_ingo_ref = msm_result_cpu_ingo_ref.add(points[i].mul(scalars[i]));
            msm_result_cpu_ref1 =
                msm_result_cpu_ref1.add_mixed(&points[i].mul(scalars[i]).into_affine());
        }

        assert_eq!(msm_result_cpu_ingo_ref, msm_result_cpu_ref1);
        assert_eq!(msm_result_cpu_ingo_ref, msm_ark_projective);

        let points_bytes = points
            .to_vec()
            .into_iter()
            .map(|point| {
                [
                    point.y.into_repr().to_bytes_le(),
                    point.x.into_repr().to_bytes_le(),
                ]
            })
            .flatten()
            .flatten()
            .collect::<Vec<_>>();

        let scalar_bytes = scalars
            .to_vec()
            .into_iter()
            .map(|scalar| scalar.into_repr().to_bytes_le())
            .flatten()
            .collect::<Vec<u8>>();

        let msm_cloud_res = ingo_x::msm_cloud_generic::<G1Affine>(&points_bytes, &scalar_bytes);

        assert_eq!(msm_cloud_res.0, msm_ark_projective); //raw vec comparison isn't always meaningful

        let points_bytes_gen = points_to_bytes::<G1Affine>(&points.to_vec());
        let scalar_bytes_gen = scalars_to_bytes::<G1Affine>(
            &scalars
                .to_vec()
                .into_iter()
                .map(|s| BigInteger256::try_from(s).unwrap())
                .collect::<Vec<BigInteger256>>(),
        );

        assert!(
            points_bytes == points_bytes_gen,
            "\n1:::{:02X?}\n2:::{:02X?}",
            points_bytes,
            points_bytes_gen
        );

        assert!(
            scalar_bytes == scalar_bytes_gen,
            "{:02X?} {:02X?}",
            scalar_bytes,
            scalar_bytes_gen
        );

        let ret = ingo_x::msm_cloud_generic::<G1Affine>(&points_bytes_gen, &scalar_bytes_gen);
        assert_eq!(ret.0, msm_ark_projective);
    }

    #[test]
    pub fn msm_correctness_g2() {
        let test_npow = std::env::var("TEST_NPOW").unwrap_or("11".to_string());
        let n_points = i32::from_str(&test_npow).unwrap();

        //TODO: conversion of inputs/outputs can be much much simplified as done for Sppark GPU and Ingo FPGA MSM

        let len = 1 << n_points;
        let (points, scalars) = ingo_x::util::generate_points_scalars::<G2Affine>(len);

        let msm_ark_projective = ingo_x::msm_ark(
            &points,
            &scalars
                .to_vec()
                .into_iter()
                .map(|s| BigInteger256::try_from(s).unwrap())
                .collect::<Vec<BigInteger256>>(), //this is safe but slow conversion
        );

        let mut msm_result_cpu_ingo_ref = G2Projective::zero(); //TODO: same as G1Affine::prime_subgroup_generator().mul(0);
        let mut msm_result_cpu_ref1 = G2Projective::zero();
        for i in 0..len {
            msm_result_cpu_ingo_ref = msm_result_cpu_ingo_ref.add(points[i].mul(scalars[i]));
            msm_result_cpu_ref1 =
                msm_result_cpu_ref1.add_mixed(&points[i].mul(scalars[i]).into_affine());
        }

        assert_eq!(msm_result_cpu_ingo_ref, msm_result_cpu_ref1);
        assert_eq!(msm_result_cpu_ingo_ref, msm_ark_projective);

        let points_bytes = points
            .to_vec()
            .into_iter()
            .map(|point| {
                [
                    point.y.c1.into_repr().to_bytes_le(),
                    point.x.c1.into_repr().to_bytes_le(),
                    point.y.c0.into_repr().to_bytes_le(),
                    point.x.c0.into_repr().to_bytes_le(),
                ]
            })
            .flatten()
            .flatten()
            .collect::<Vec<u8>>();

        let scalar_bytes = scalars
            .to_vec()
            .into_iter()
            .map(|scalar| scalar.into_repr().to_bytes_le())
            .flatten()
            .collect::<Vec<u8>>();

        let msm_cloud_res = ingo_x::msm_cloud_generic::<G2Affine>(&points_bytes, &scalar_bytes);

        assert_eq!(msm_cloud_res.0, msm_ark_projective); //raw vec comparison isn't always meaningful

        let points_bytes_gen = points_to_bytes::<G2Affine>(&points.to_vec());
        let scalar_bytes_gen = scalars_to_bytes::<G2Affine>(
            &scalars
                .to_vec()
                .into_iter()
                .map(|s| BigInteger256::try_from(s).unwrap())
                .collect::<Vec<BigInteger256>>(),
        );

        assert!(
            points_bytes == points_bytes_gen,
            "\n1:::{:02X?}\n2:::{:02X?}",
            points_bytes,
            points_bytes_gen
        );

        assert!(
            scalar_bytes == scalar_bytes_gen,
            "{:02X?} {:02X?}",
            scalar_bytes,
            scalar_bytes_gen
        );

        let ret = ingo_x::msm_cloud_generic::<G2Affine>(&points_bytes_gen, &scalar_bytes_gen);
        assert_eq!(ret.0, msm_ark_projective);
    }
}
