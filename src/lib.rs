use ark_ec::AffineCurve;
use ark_ff::{PrimeField, BigInteger, FromBytes};
use rust_rw_device::rw_msm_to_dram::msm_core;

const BYTE_SIZE_POINT_COORD: usize = 48;
const BYTE_SIZE_SCALAR: usize = 32;

fn get_formatted_unified_points_from_affine<G: AffineCurve>(points: &[G]) -> Vec<u8> {
    let mut buff = vec![]; 
    let point_size = buff.len() / 2;
    G::zero().write(&mut buff).unwrap();

    let mut points_buffer = Vec::<u8>::with_capacity(points.len() * BYTE_SIZE_POINT_COORD);

    for (i, base) in points.iter().enumerate() {
        base.write(&mut buff).unwrap();
        points_buffer[2*i*point_size..(2*i+1)*point_size].copy_from_slice(&buff[point_size..2*point_size]);
        points_buffer[(2*i+1)*point_size..(2*i+2)*point_size].copy_from_slice(&buff[0..point_size]);
        points_buffer.extend(std::iter::repeat(0).take(BYTE_SIZE_POINT_COORD - points_buffer.len()));
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
        let points_bytes = get_formatted_unified_points_from_affine(bases);
        let scalars_bytes = get_formatted_unified_scalars_from_bigint::<G>(scalars);

        let (z_chunk, y_chunk, x_chunk, _, _) = msm_core(points_bytes, scalars_bytes, scalars.len());
        let mut result_buffer = Vec::new(); 
        result_buffer.extend_from_slice(&x_chunk);
        result_buffer.extend_from_slice(&y_chunk);
        result_buffer.extend_from_slice(&z_chunk);

        G::Projective::read(result_buffer.as_slice()).unwrap()
    }
}