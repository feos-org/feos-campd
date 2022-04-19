from .feos_campd import *
from .feos_campd import (
    SuperMolecule,
    OptimizationResult,
    OptimizationProblem,
    MetaOptimizationProblem,
    OrganicRankineCycle,
    OrganicRankineCycleSuperStructure,
    FixedMolecule,
)
from knitro import *


def setup_knitro_supermolecule(self, kc, solutions):
    n_y = self.variables

    indexVars = KN_add_vars(kc, n_y)
    KN_set_var_types(kc, indexVars, [KN_VARTYPE_BINARY] * n_y)

    # maximum size
    indexCon = KN_add_con(kc)
    coeffs, size = self.size_constraint
    KN_add_con_linear_struct(kc, [indexCon] * n_y, indexVars, coeffs)
    KN_set_con_upbnds(kc, indexCon, size)

    # integer cuts
    for solution in solutions:
        indexCon = KN_add_con(kc)
        coefs = [2 * yi - 1 for yi in solution.y]
        KN_add_con_linear_struct(kc, [indexCon] * n_y, indexVars, coefs)
        KN_set_con_upbnds(kc, indexCon, sum(solution.y) - 1)

    # functional group
    indexVars = self.functional_group_constraint
    if len(indexVars) > 0:
        indexCon = KN_add_con(kc)
        KN_add_con_linear_struct(
            kc, [indexCon] * len(indexVars), indexVars, [1] * len(indexVars)
        )
        KN_set_con_eqbnds(kc, indexCon, 1)

    # minimum size
    indexCon = KN_add_con(kc)
    KN_add_con_linear_struct(kc, indexCon, len(indexVars), 1)
    KN_set_con_eqbnds(kc, indexCon, 1)

    # bond constraints
    for bond in self.bond_constraints:
        indexCon = KN_add_con(kc)
        KN_add_con_linear_struct(kc, [indexCon] * 2, bond, [1, -1])
        KN_set_con_lobnds(kc, indexCon, 0)

    # symmetry constraints
    for indexVars, coefs in self.symmetry_constraints:
        indexCon = KN_add_con(kc)
        KN_add_con_linear_struct(kc, [indexCon] * len(coefs), indexVars, coefs)
        KN_set_con_lobnds(kc, indexCon, 0)


def setup_knitro_fixed_molecule(self, kc, solutions):
    n_y = self.variables

    indexVars = KN_add_vars(kc, n_y)
    KN_set_var_types(kc, indexVars, [KN_VARTYPE_BINARY] * n_y)

    # single component constraint
    indexCon = KN_add_con(kc)
    KN_add_con_linear_struct(kc, [indexCon] * n_y, indexVars, [1] * n_y)
    KN_set_con_eqbnds(kc, indexCon, 1)

    # integer cuts
    for solution in solutions:
        indexCon = KN_add_con(kc)
        coefs = [2 * yi - 1 for yi in solution.y]
        KN_add_con_linear_struct(kc, [indexCon] * n_y, indexVars, coefs)
        KN_set_con_upbnds(kc, indexCon, sum(solution.y) - 1)


def setup_knitro_process(self, kc, x0):
    # declare continuous variables
    n_x = len(self.variables)
    indexVars = KN_add_vars(kc, n_x)
    for i, (l, u) in zip(indexVars, self.variables):
        if l is not None:
            KN_set_var_lobnds(kc, i, l)
        if u is not None:
            KN_set_var_upbnds(kc, i, u)

    # set initial values
    KN_set_var_primal_init_values(kc, indexVars, x0[:n_x])

    # declare binary variables
    if self.binary_variables > 0:
        indexVars = KN_add_vars(kc, self.binary_variables)
        KN_set_var_types(kc, indexVars, [KN_VARTYPE_BINARY] * len(indexVars))

    # add constraints
    n_cons = len(self.constraints)
    indexCons = KN_add_cons(kc, n_cons)
    for i, (l, u, e) in zip(indexCons, self.constraints):
        if l is not None:
            KN_set_con_lobnds(kc, i, l)
        if u is not None:
            KN_set_con_upbnds(kc, i, u)
        if e is not None:
            KN_set_con_eqbnds(kc, i, e)

    return indexCons


SuperMolecule.setup_knitro = setup_knitro_supermolecule
FixedMolecule.setup_knitro = setup_knitro_fixed_molecule
OrganicRankineCycle.setup_knitro = setup_knitro_process
OrganicRankineCycleSuperStructure.setup_knitro = setup_knitro_process


def solve_knitro_once(
    molecule,
    property_model=None,
    process=None,
    solutions=[],
    x0=None,
    options=None,
    name=None,
):
    n_y = molecule.variables

    def callback(kc, cb, evalRequest, evalResult, userParams):
        if process is None:
            evalResult.obj = 0
            return 0
        y = evalRequest.x[:n_y]
        x = evalRequest.x[n_y:]
        eos = property_model.build_eos(molecule.build(y))
        try:
            _, target, constraints = process.solve(eos, x)
            evalResult.obj = target
            evalResult.c = constraints
        except:
            evalResult.obj = 0
            evalResult.c = [-1] * len(process.constraints)
        return 0

    kc = KN_new()

    molecule.setup_knitro(kc, solutions)
    if process is None:
        indexCons = []
    else:
        indexCons = process.setup_knitro(kc, x0)
    KN_add_eval_callback(kc, True, indexCons, callback)
    if options is not None:
        KN_load_param_file(kc, options)
    if KN_solve(kc) != 0:
        return
    _, t, x, _ = KN_get_solution(kc)
    y = [int(i) for i in x[:n_y]]
    x = x[n_y:]
    if name is None:
        print(t, y, x)
    else:
        print(f"{name:12} {t:.5f}", y, x)
    return OptimizationResult(t, y, x)


def solve_knitro(self, x0, n_solutions=1, options=None):
    for _ in range(n_solutions):
        result = solve_knitro_once(
            self.molecule,
            self.property_model,
            self.process,
            self.solutions,
            x0,
            options,
        )
        self.add_solution(result)


OptimizationProblem.solve_knitro = solve_knitro


def solve_knitro_meta(self, x0, n_solutions=1, options=None):
    candidates = {}
    for name, molecule in self.molecules.items():
        candidates[name] = solve_knitro_once(
            molecule, self.property_model, self.process, [], x0, options, name
        )
    self.initialize_candidates(candidates)

    for _ in range(n_solutions):
        chemical = self.best_candidate()
        solutions = self.candidates[chemical]
        result = solve_knitro_once(
            self.molecules[chemical],
            self.property_model,
            self.process,
            solutions,
            x0,
            options,
            chemical,
        )
        self.update_candidates(chemical, result)


MetaOptimizationProblem.solve_knitro = solve_knitro_meta


def _repr_png_(self):
    from rdkit import Chem

    return self.molecule([1] * self.variables)._repr_png_()


def canonical_smiles(self, y):
    from rdkit import Chem

    return Chem.MolToSmiles(Chem.MolFromSmiles(self.smiles(y)))


def molecule(self, y):
    from rdkit import Chem

    return Chem.MolFromSmiles(self.canonical_smiles(y))


SuperMolecule._repr_png_ = _repr_png_
SuperMolecule.canonical_smiles = canonical_smiles
SuperMolecule.molecule = molecule
