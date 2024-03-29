use feos::ideal_gas::{DipprRecord, JobackRecord};
use serde::{Deserialize, Serialize};

pub const JOBACK_PARAMETER_OFFSETS: [f64; 5] = [-37.93, 0.21, -3.91e-4, 2.06e-7, 0.0];

#[derive(Clone, Serialize, Deserialize)]
#[serde(untagged)]
pub enum IdealGasRecord {
    Joback(JobackRecord),
    Dippr(DipprRecord),
}

impl IdealGasRecord {
    pub fn joback_params(&self) -> [f64; 5] {
        match &self {
            IdealGasRecord::Joback(JobackRecord { a, b, c, d, e }) => [*a, *b, *c, *d, *e],
            IdealGasRecord::Dippr(dippr) => {
                let DipprRecord::DIPPR100(v) = &dippr else {
                    panic!("All components have to use Eq. 100!")
                };
                let mut v = v.to_vec();
                if v.len() < 5 {
                    v.extend(vec![0.0; 5 - v.len()]);
                }
                let &[a, b, c, d, e, ..] = &v[..] else {
                    panic!("Eq. 100 needs to be provided with up to 5 coefficients!")
                };
                let [a0, b0, c0, d0, e0] = JOBACK_PARAMETER_OFFSETS;
                [
                    a * 0.001 - a0,
                    b * 0.001 - b0,
                    c * 0.001 - c0,
                    d * 0.001 - d0,
                    e * 0.001 - e0,
                ]
            }
        }
    }
}
