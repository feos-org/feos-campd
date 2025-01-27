use crate::{
    process::{ContinuousVariable, GeneralConstraint, ProcessModel},
    ChemicalRecord,
};
use feos_ad::{HelmholtzEnergyWrapper, PhaseEquilibriumAD, StateAD, TotalHelmholtzEnergy};
use feos_core::{DensityInitialization, EosResult, ReferenceSystem};
use nalgebra::SVector;
use num_dual::DualNum;
use quantity::{
    HeatCapacityRate, MoleFlowRate, Power, Pressure, Temperature, BAR, CELSIUS, KELVIN, KILO, MEGA,
    WATT,
};

struct PressureChanger<'a, E: TotalHelmholtzEnergy<1>, D: DualNum<f64> + Copy> {
    inlet: &'a StateAD<'a, E, D, 1>,
    outlet: StateAD<'a, E, D, 1>,
}

impl<'a, E: TotalHelmholtzEnergy<1>, D: DualNum<f64> + Copy> PressureChanger<'a, E, D> {
    fn pump(
        inlet: &'a StateAD<'a, E, D, 1>,
        pressure: Pressure<D>,
        efficiency: f64,
    ) -> EosResult<Self> {
        // calculate isentropic state
        let state2s = StateAD::new_ps(
            inlet.eos,
            pressure,
            inlet.molar_entropy(),
            inlet.molefracs,
            DensityInitialization::Liquid,
            Some(inlet.temperature.re()),
        )?;

        // calculate real state
        let h1 = inlet.molar_enthalpy();
        let h2s = state2s.molar_enthalpy();
        let h2 = h1 + (h2s - h1) / efficiency;
        let outlet = StateAD::new_ph(
            inlet.eos,
            pressure,
            h2,
            inlet.molefracs,
            DensityInitialization::Liquid,
            Some(state2s.temperature.re()),
        )?;

        Ok(Self { inlet, outlet })
    }

    fn turbine(
        inlet: &'a StateAD<'a, E, D, 1>,
        pressure: Pressure<D>,
        vle_condenser: &PhaseEquilibriumAD<'a, E, D, 1>,
        efficiency: f64,
    ) -> EosResult<Self> {
        // calculate isentropic state
        let s1 = inlet.molar_entropy();
        let s_v = vle_condenser.vapor.molar_entropy();
        let h2s = if s1.re() < s_v.re() {
            let s_l = vle_condenser.liquid.molar_entropy();
            let h_l = vle_condenser.liquid.molar_enthalpy();
            let h_v = vle_condenser.vapor.molar_enthalpy();
            let x = (s1 - s_l) / (s_v - s_l);
            h_l + x * (h_v - h_l)
        } else {
            StateAD::new_ps(
                inlet.eos,
                pressure,
                s1,
                inlet.molefracs,
                DensityInitialization::Vapor,
                Some(inlet.temperature.re()),
            )?
            .molar_enthalpy()
        };

        // calculate real state
        let h1 = inlet.molar_enthalpy();
        let h2 = h1 + (h2s - h1) * efficiency;
        let outlet = StateAD::new_ph(
            inlet.eos,
            pressure,
            h2,
            inlet.molefracs,
            DensityInitialization::Vapor,
            Some(inlet.temperature.re()),
        )?;

        Ok(Self { inlet, outlet })
    }

    fn power(&self, flow_rate: MoleFlowRate<D>) -> Power<D> {
        flow_rate * (self.outlet.molar_enthalpy() - self.inlet.molar_enthalpy())
    }
}

pub struct OrganicRankineCycle {
    c_p_hs: HeatCapacityRate,
    t_hs: Temperature,
    dt_hs: Temperature,
    eta_st: f64,
    eta_sp: f64,
    p_min: Pressure,
    p_max: Pressure,
    p_min_r: f64,
    p_max_r: f64,
    t_cool_in: Temperature,
    t_cool_out: Temperature,
    dt_cool: Temperature,
}

impl<E: TotalHelmholtzEnergy<1>> ProcessModel<E, 3, 1> for OrganicRankineCycle {
    fn variables(&self) -> [ContinuousVariable; 3] {
        [
            ContinuousVariable::new(298.15 / 300.0, 400.15 / 300.0, 1.0),
            ContinuousVariable::new(298.15 / 300.0, 400.15 / 300.0, 1.2),
            ContinuousVariable::new(0.0, 1.0, 0.1),
        ]
    }

    fn constraints(&self) -> Vec<GeneralConstraint> {
        vec![
            GeneralConstraint::Inequality(Some(1.0), None),
            GeneralConstraint::Inequality(Some(1.0), None),
            GeneralConstraint::Inequality(Some(self.p_min_r), Some(self.p_max_r)),
            GeneralConstraint::Inequality(Some(self.p_min_r), Some(self.p_max_r)),
            GeneralConstraint::Inequality(
                Some(self.p_min.convert_into(BAR)),
                Some(self.p_max.convert_into(BAR)),
            ),
            GeneralConstraint::Inequality(
                Some(self.p_min.convert_into(BAR)),
                Some(self.p_max.convert_into(BAR)),
            ),
            GeneralConstraint::Inequality(Some(1.0), None),
            GeneralConstraint::Inequality(Some(1.0), None),
        ]
    }

    fn evaluate<D: DualNum<f64> + Copy>(
        &self,
        eos: &HelmholtzEnergyWrapper<E, D, 1>,
        _: [&ChemicalRecord<D>; 1],
        x: [D; 3],
    ) -> EosResult<(D, Vec<D>)> {
        // unpack process variables
        let [t_cond, t_evap, dt_sh] = x;
        let t_cond = Temperature::from_reduced(t_cond * 300.0);
        let t_evap = Temperature::from_reduced(t_evap * 300.0);
        let dt_sh = Temperature::from_reduced(dt_sh * 50.0);
        let molefracs = SVector::from([D::one()]);

        // calculate isobars
        let (vle_cond, p_cond) = PhaseEquilibriumAD::new_t(eos, t_cond)?;
        let (vle_evap, p_evap) = PhaseEquilibriumAD::new_t(eos, t_evap)?;

        // calculate pump
        let pump = PressureChanger::pump(&vle_cond.liquid, p_evap, self.eta_sp)?;

        // calculate superheating
        let turbine_in = StateAD::new_tp(
            eos,
            t_evap + dt_sh,
            p_evap,
            molefracs,
            DensityInitialization::Vapor,
        )?;

        // calculate turbine
        let turbine = PressureChanger::turbine(&turbine_in, p_cond, &vle_cond, self.eta_st)?;

        // calculate mass flow rate
        let t_hs_pinch = vle_evap.liquid.temperature + self.dt_hs;
        let m_wf = (t_hs_pinch - self.t_hs) * self.c_p_hs
            / (vle_evap.liquid.molar_enthalpy() - turbine_in.molar_enthalpy());

        // target
        let target = (pump.power(m_wf) + turbine.power(m_wf)).convert_into(MEGA * WATT);

        // pinch constraint heat source
        let pinch_hs = ((-turbine_in.temperature + self.t_hs) / self.dt_hs).into_value();

        // pinch constraint condenser
        let h1 = vle_cond.liquid.molar_enthalpy();
        let h2 = vle_cond.vapor.molar_enthalpy();
        let h3 = turbine.outlet.molar_enthalpy();
        let t_cool_pinch =
            (h2 - h1) / (h3 - h1) * (self.t_cool_out - self.t_cool_in) + self.t_cool_in;
        let pinch_cond = ((t_cool_pinch - t_cond) / self.dt_cool).into_value();

        // critical pressure
        let p_crit = StateAD::critical_point_pure(eos)?.pressure();

        // reduced pressure constraints
        let pr_cond = p_cond.convert_into(p_crit);
        let pr_evap = p_evap.convert_into(p_crit);

        // absolute pressure constraints
        let pa_cond = p_cond.convert_into(BAR);
        let pa_evap = p_evap.convert_into(BAR);

        // pressure constraint
        let p_const = p_evap.convert_into(p_cond);

        // turbine outlet constraint
        let turb_out = turbine.outlet.temperature.convert_into(t_cond);

        let constraints = vec![
            pinch_hs, pinch_cond, pr_cond, pr_evap, pa_cond, pa_evap, p_const, turb_out,
        ];

        Ok((target, constraints))
    }
}

impl Default for OrganicRankineCycle {
    fn default() -> Self {
        Self {
            c_p_hs: 65.0 * KILO * WATT / KELVIN,
            t_hs: 175.0 * CELSIUS,
            dt_hs: 20.0 * KELVIN,
            eta_st: 0.65,
            eta_sp: 0.8,
            p_min: BAR,
            p_max: 50.0 * BAR,
            p_min_r: 1e-3,
            p_max_r: 0.8,
            t_cool_in: 25. * CELSIUS,
            t_cool_out: 40.0 * CELSIUS,
            dt_cool: -10.0 * KELVIN,
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use approx::assert_relative_eq;
    use feos_ad::eos::{Joback, PcSaftPure};
    use feos_ad::{EquationOfStateAD, ParametersAD};
    use feos_core::EosResult;
    use std::collections::HashMap;

    #[test]
    fn test_orc() -> EosResult<()> {
        let pcsaft = PcSaftPure([
            1.93988,
            3.229270435095473,
            225.21099088603418,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ]);
        let joback = Joback([
            13.939999999999998,
            0.19491999999999998,
            -0.00013533000000000002,
            0.000000042889999999999996,
            0.0,
        ]);
        let eos = EquationOfStateAD::new([joback], pcsaft).wrap();
        let x0 = [1.0, 1.2, 0.1];

        let orc = OrganicRankineCycle::default();
        let cr = ChemicalRecord::new(HashMap::new(), HashMap::new());
        let (obj, con) = orc.evaluate(&eos, [&cr], x0)?;
        println!("{obj} {con:?}");

        let ref_obj = -0.601275463470434;
        let ref_con = [
            4.157499999999999,
            -1.275381108337291,
            0.13413664979435497,
            0.5086603251144651,
            9.56756117623695,
            36.281201192353244,
            3.7921054826871905,
            1.0230776382730467,
        ];

        assert_relative_eq!(obj, ref_obj, max_relative = 1e-9);
        for (c, r) in con.into_iter().zip(ref_con) {
            assert_relative_eq!(c, r, max_relative = 1e-9);
        }

        Ok(())
    }
}
