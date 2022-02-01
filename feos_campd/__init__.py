from operator import index
from .feos_campd import *
from .feos_campd import (
    SuperMolecule,
    OptimizationResult,
    OptimizationProblem,
    OrganicRankineCycle,
    HomoGc,
    FixedMolecule,
)
from knitro import *


def setup_knitro_supermolecule(self, kc, solutions):
    n_y = self.variables

    indexVars = KN_add_vars(kc, n_y)
    KN_set_var_types(kc, indexVars, [KN_VARTYPE_BINARY] * n_y)

    # maximum size
    indexCon = KN_add_con(kc)
    KN_add_con_linear_struct(kc, [indexCon] * n_y, indexVars, [1] * n_y)
    KN_set_con_upbnds(kc, indexCon, self.size + 1)

    # integer cuts
    for solution in solutions:
        indexCon = KN_add_con(kc)
        coefs = [2 * yi - 1 for yi in solution.y]
        KN_add_con_linear_struct(kc, [indexCon] * n_y, indexVars, coefs)
        KN_set_con_upbnds(kc, indexCon, sum(solution.y) - 1)

    # functional group
    indexVars = self.functional_group_constraint
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


def setup_knitro_homo_gc(self, kc, solutions):
    # alkyl groups
    indexVars_c = KN_add_vars(kc, 4)
    KN_set_var_types(kc, indexVars_c, [KN_VARTYPE_INTEGER] * 4)
    KN_set_var_lobnds(kc, indexVars_c, [0, 0, 0, 0])
    KN_set_var_upbnds(kc, indexVars_c, [self.size] * 4)

    # maximum size
    indexCon = KN_add_con(kc)
    KN_add_con_linear_struct(kc, [indexCon] * 4, indexVars_c, [1, 1, 1, 1])
    KN_set_con_upbnds(kc, indexCon, self.size)

    # octett rule
    indexCon = KN_add_con(kc)
    KN_add_con_linear_struct(kc, [indexCon] * 4, indexVars_c, [-2, -1, 0, 1])
    KN_set_con_eqbnds(kc, indexCon, 1)

    # functional groups
    indexVars_f = KN_add_vars(kc, self.variables - 4)
    KN_set_var_types(kc, indexVars_f, [KN_VARTYPE_BINARY] * len(indexVars_f))
    indexCon = KN_add_con(kc)
    KN_add_con_linear_struct(
        kc, [indexCon] * len(indexVars_f), indexVars_f, [1] * len(indexVars_f)
    )
    KN_set_con_eqbnds(kc, indexCon, 1)

    # integer cuts
    indexVars = indexVars_c + indexVars_f
    for solution in solutions:
        indexCon = KN_add_con(kc)
        qcoefs = [1] * len(solution.y)
        lcoefs = [-2 * yi for yi in solution.y]
        KN_add_con_quadratic_struct(
            kc, [indexCon] * len(indexVars), indexVars, indexVars, qcoefs
        )
        KN_add_con_linear_struct(kc, [indexCon] * len(indexVars), indexVars, lcoefs)
        KN_set_con_lobnds(kc, indexCon, 0.5 - sum(yi ** 2 for yi in solution.y))


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
    n_x = len(self.variables)

    # set bounds
    indexVars = KN_add_vars(kc, n_x)
    for i, (l, u) in zip(indexVars, self.variables):
        if l is not None:
            KN_set_var_lobnds(kc, i, l)
        if u is not None:
            KN_set_var_upbnds(kc, i, u)

    # set initial values
    KN_set_var_primal_init_values(kc, indexVars, x0)

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
HomoGc.setup_knitro = setup_knitro_homo_gc
FixedMolecule.setup_knitro = setup_knitro_fixed_molecule
OrganicRankineCycle.setup_knitro = setup_knitro_process


def solve_knitro(self, x0, n_solutions=1, options=None):
    n_y = self.molecule.variables

    def callback(kc, cb, evalRequest, evalResult, userParams):
        y = evalRequest.x[:n_y]
        x = evalRequest.x[n_y:]
        eos = self.property_model.build_eos(self.molecule.build(y))
        try:
            _, target, constraints = self.process.solve(eos, x)
            evalResult.obj = target
            evalResult.c = constraints
        except:
            evalResult.obj = 0
            evalResult.c = [-1] * len(self.process.constraints)
        return 0

    for _ in range(n_solutions):
        kc = KN_new()

        self.molecule.setup_knitro(kc, self.solutions)
        indexCons = self.process.setup_knitro(kc, x0)
        KN_add_eval_callback(kc, True, indexCons, callback)
        if options is not None:
            KN_load_param_file(kc, options)
        if KN_solve(kc) != 0:
            return
        _, t, x, _ = KN_get_solution(kc)
        y = [int(i) for i in x[:n_y]]
        x = x[n_y:]
        print(t, y, x)
        self.add_solution(OptimizationResult(t, y, x))


OptimizationProblem.solve_knitro = solve_knitro
