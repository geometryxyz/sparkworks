use ark_ec::AffineCurve;
use ark_ff::{PrimeField, BigInteger, FromBytes, ToBytes, Field, BigInteger384, BigInteger256};
use rust_rw_device::rw_msm_to_dram::msm_core;
use core::num;
use std::io::Cursor;
use std::ops::Mul;

const BYTE_SIZE_SCALAR: usize = 32;

fn get_formatted_unified_points_from_affine<G: AffineCurve>(points: &[G]) -> Vec<u8> {
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
        points_buffer[2*i*point_size..(2*i+1)*point_size].copy_from_slice(&buff.get_ref()[point_size..2*point_size]);
        // write x
        points_buffer[(2*i+1)*point_size..(2*i+2)*point_size].copy_from_slice(&buff.get_ref()[0..point_size]);
    }

    points_buffer
}

fn get_formatted_unified_scalars_from_bigint<G: AffineCurve>(scalars: &[<G::ScalarField as PrimeField>::BigInt]) -> Vec<u8> {
    let mut scalars_bytes: Vec<u8> = Vec::new();
    for i in 0..scalars.len(){
        let mut bytes_array = scalars[i].to_bytes_le();
        bytes_array.extend(std::iter::repeat(0).take(BYTE_SIZE_SCALAR - bytes_array.len()));
        scalars_bytes.extend(bytes_array);
    }
    scalars_bytes
}

pub struct FpgaVariableBaseMSM;

impl FpgaVariableBaseMSM {
    pub fn multi_scalar_mul<G: AffineCurve>(
        bases: &[G],
        scalars: &[<G::ScalarField as PrimeField>::BigInt],
    ) -> G::Projective {
        let size = std::cmp::min(bases.len(), scalars.len());
        let scalars: &[<G::ScalarField as PrimeField>::BigInt] = &scalars[..size];
        let bases = &bases[..size];

        let points_bytes = get_formatted_unified_points_from_affine(bases);
        let scalars_bytes = get_formatted_unified_scalars_from_bigint::<G>(scalars);

        let (z_chunk, y_chunk, x_chunk, _, _) = msm_core(points_bytes.clone(), scalars_bytes.clone(), scalars.len());

        let proj_x_field = <<<G as AffineCurve>::BaseField as Field>::BasePrimeField as PrimeField>::from_le_bytes_mod_order(&x_chunk);
        let proj_y_field = <<<G as AffineCurve>::BaseField as Field>::BasePrimeField as PrimeField>::from_le_bytes_mod_order(&y_chunk);
        let proj_z_field = <<<G as AffineCurve>::BaseField as Field>::BasePrimeField as PrimeField>::from_le_bytes_mod_order(&z_chunk);

        let aff_x = proj_x_field.mul(proj_z_field.inverse().unwrap());
        let aff_y = proj_y_field.mul(proj_z_field.inverse().unwrap());

        let mut buff = Vec::<u8>::with_capacity(97);
        aff_x.write(&mut buff).unwrap();
        aff_y.write(&mut buff).unwrap();
        buff.push(0);

        G::read(buff.as_slice()).unwrap().into_projective()
        // (points_bytes, scalars_bytes, z_chunk, y_chunk, x_chunk)
    }
}

#[cfg(test)]
mod test {
    use ark_bls12_377::{G1Affine, G1Projective, Fq};
    use ark_ec::{AffineCurve, ProjectiveCurve};
    use ark_ff::{UniformRand, PrimeField, ToBytes, Zero};
    use ark_std::{test_rng, rand::Rng};
    use num_bigint::BigUint;
    use super::get_formatted_unified_points_from_affine;

    const BYTE_SIZE_POINT_COORD: usize = 48; // for BLS

    // ingonyama's implementation for asserting equality 
    fn get_formatted_unified_points_from_biguint(points: &Vec<BigUint>) -> Vec<u8> {
        let mut points_bytes: Vec<u8> = Vec::new();
        for i in 0..points.len(){
            let mut bytes_array = points[i].to_bytes_le();
            bytes_array.extend(std::iter::repeat(0).take(BYTE_SIZE_POINT_COORD - bytes_array.len()));
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

        let points_as_big_int = points.iter()
        .map(|point| [point.y.into_repr().into(), point.x.into_repr().into()])
        .flatten()
        .collect::<Vec<BigUint>>();

        let point_bytes_biguint = get_formatted_unified_points_from_biguint(&points_as_big_int);
        let point_bytes_affine = get_formatted_unified_points_from_affine(&points);

        assert_eq!(point_bytes_biguint, point_bytes_affine);
    }

    #[test]
    fn test_projective_read() {
        let mut rng = test_rng();

        let point = G1Projective::rand(&mut rng); 
        let point = G1Projective::zero(); 

        let mut buff = vec![];

        point.write(&mut buff).unwrap();

        println!("x: {}", point.x);
        println!("y: {}", point.y);
        println!("z: {}", point.z);

        // let buff = vec![Fq::zero(); 144];

        println!("{:?}, {}", buff.iter().take(48), buff.len());
    }

    #[test]
    fn write_affine() {
        let mut rng = test_rng();

        let point = G1Projective::rand(&mut rng).into_affine(); 

        let mut buff = vec![];

        // write by coord
        point.x.write(&mut buff).unwrap();
        point.y.write(&mut buff).unwrap();
        buff.push(0);

        // write full 
        let mut buff_full = vec![];
        point.write(&mut buff_full).unwrap();
        
        assert_eq!(buff, buff_full);
        // println!("{:?}, {}", buff.iter(), buff.len());

    }
}