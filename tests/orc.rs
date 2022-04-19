use feos_campd::OrganicRankineCycleSuperStructure;
use feos_core::parameter::{IdentifierOption, Parameter};
use feos_core::EosResult;
use feos_pcsaft::{PcSaft, PcSaftParameters};
use std::fs::File;
use std::io::BufReader;
use std::rc::Rc;

#[test]
pub fn test_orc() -> EosResult<()> {
    let orc: OrganicRankineCycleSuperStructure =
        serde_json::from_reader(BufReader::new(File::open("tests/orc.json").unwrap())).unwrap();
    let params = PcSaftParameters::from_json(
        vec!["propane"],
        "tests/gross2001.json",
        None,
        IdentifierOption::Name,
    )?;
    let eos = Rc::new(PcSaft::new(Rc::new(params)));
    orc.solve(&eos, &[1.15, -4.5, -1.2, 0.1, 0.1, 0.5])?;
    Ok(())
}
