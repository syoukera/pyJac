""" Cache Optimizer

    Reorders loads of rate and species subs to optimize cache hits, etc.
"""

# Python 2 compatibility
from __future__ import division
from __future__ import print_function
from CUDAParams import Jacob_Unroll, ResetOnJacUnroll
from CParams import C_Jacob_Unroll
import utils
import multiprocessing
import pickle
import timeit

try:
    from Numberjack import *
    Model().load('SCIP')
    HAVE_NJ = True
except Exception, e:
    HAVE_NJ = False
    print(e)
    print('Cache-optimization support disabled...')

def get_mappings(specs, reacs, load_non_participating=False, consider_thd=False):
    r_to_s = [set() for i in range(len(reacs))]
    s_to_r = [set() for i in range(len(specs))]
    for sind, sp in enumerate(specs):
        for rind, rxn in enumerate(reacs):
            the_set = set(rxn.reac + rxn.prod)
            if consider_thd:
                thd_sp = [thd[0] for thd in rxn.thd_body_eff if thd[1] != 0]
                if thd_sp:
                    the_set = the_set.union(set(thd_sp))
            if any(sp.name == x for x in the_set):
                nu = utils.get_nu(sp, rxn)
                if (nu is not None and nu != 0) or load_non_participating:
                    r_to_s[rind].add(sind)
                    s_to_r[sind].add(rind)

    return r_to_s, s_to_r


def greedy_optimizer(lang, specs, reacs, multi_thread, force_optimize, build_path, time_lim=60):
    """
    An optimization method that reorders the species and reactions in a method to attempt to keep data in cache as
    long as possible

    Notes
    -----
    This method optimizes based on Jacobian matrix generation, as this is the most important and time consuming step
    of the reaction rate subroutines
    Species and reactions are reordered to optimize the cache rates for this.

    Next orderings for the evaluation of reactions, pressure dependent reactions and species rates are determined in
    order to optimize the various rate routines

    Parameters
    ----------
    lang : str
        The language
    specs : list of SpecInfo
        List of species in the mechanism.
    reacs : list of ReacInfo
        List of reactions in the mechanism.
    multi_thread : int
        The number of threads to use during optimization
    force_optimize : bool
        If true, reoptimize even if past data is available
    build_path : str
        The path to the build directory
    time_lim : int
        The time limit for optimization operations in minutes

    Returns
    _______

    specs : list of SpecInfo
        The reordered list of species in the mechanism

    reacs : list of ReacInfo
        the reordered list of reacs in the mechanism

    rxn_rate_order : list of int
        A list indicating the order to evaluate reaction rates in rxn_rates subroutine

    pdep_rate_order : list of int
        A list indicating the order to evaluate pressure dependent rates in rxn_rates_pres_mod subroutine

    spec_rate_order : list of int
        A list indicated the order to evaluate species rates in spec_rates subroutine

    old_species_order : list of int
        A list indicating the positioning of the species in the original mechanism, used in rate testing

     old_rxn_order : list of int
        A list indicating the positioning of the reactions in the original mechanism, used in rate testing
    """

    if not HAVE_NJ:
        print("Cache-optimization disabled, returning original mechanism")
        return None

    unroll_len = C_Jacob_Unroll if lang == 'c' else Jacob_Unroll
    # first try to load past data
    if not force_optimize:
        try:
            same_mech = False
            with open(build_path + 'optimized.pickle', 'rb') as file:
                old_specs = pickle.load(file)
                old_reacs = pickle.load(file)
                rxn_rate_order = pickle.load(file)
                spec_rate_order = pickle.load(file)
                pdep_rate_order = pickle.load(file)
                spec_ordering = pickle.load(file)
                rxn_ordering = pickle.load(file)
            same_mech = all(any(s == sp for sp in specs) for s in old_specs) and \
                        len(specs) == len(old_specs) and \
                        all(any(r == rxn for rxn in reacs) for r in old_reacs) and \
                        len(reacs) == len(old_reacs)
        except Exception, e:
            same_mech = False
        if same_mech:
            # we have to do the spec_rate_order each time
            return old_specs, old_reacs, rxn_rate_order, pdep_rate_order, spec_rate_order, spec_ordering,\
                   rxn_ordering

    #get the mappings
    r_to_s, s_to_r = get_mappings(specs, reacs, consider_thd=False, load_non_participating=True)

    #set up the Numberjack constraint problem to optimize reaction order
    model = Model()

    rxn_order = [Variable(0, len(reacs), "rxn_{}".format(i)) for i in range(len(reacs))]

    #set up the constraints to ensure unique rxn placement
    model.add(AllDiff(rxn_order))

    #now set up score function
    score_matrix = Matrix(len(specs), len(reacs), 2)

    class rxn_comp(Expression):
        def __init__ (self, reac1, reac2):
            Expression.__init__(self, "rxn_comp")
            self.set_children(reac1, reac2)
            self.cols = [reac1, reac2]

        def decompose(self):
            return [ Sum ([v1 < v2 for v1, v2 in zip(self.cols[0], self.cols[1])]) ]

    #set up the r_to_s constraints
    for i, rxn in enumerate(rxn_order):
        for sp in range(len(specs)):
            val = 1 if sp in r_to_s[i] else 0
            model.add(score_matrix[sp, rxn] == val)

    score_list = []
    for i in range(len(reacs)):
        if i + 1 == len(reacs):
            continue
        score_list.append(rxn_comp(score_matrix.col[i], score_matrix.col[i + 1]))

    score = Sum(score_list)
    model.add(Minimize(score))
    solver = model.load('SCIP')
    solver.setThreadCount(multi_thread)
    solver.setTimeLimit(time_lim * 60)
    solver.setVerbosity(2)

    def __rxn_score_check(r_to_s, ordering):
        the_sum = 0
        for i in range(len(r_to_s) - 1):
            the_sum += len(r_to_s[ordering[i + 1]].difference(r_to_s[ordering[i]]))

        return the_sum

    solved = solver.solve()
    solns = [x for x in solver.solutions()]
    print(solved, solver.is_opt(), len(solns))
    if solved and solver.is_opt():
        ordering = [x.get_value() for x in rxn_order]
    elif solved:
        bestVal = None
        for solution in solns:
            val = __rxn_score_check(r_to_s, solution)
            if bestVal is None or val < bestVal:
                ordering[:] = solution
                bestVal = val
    else:
        raise Exception('No solution found, try a longer timelimit')

    print([x for x in enumerate(ordering)])
    pre = __rxn_score_check(r_to_s, range(len(reacs)))
    post = __rxn_score_check(r_to_s, ordering)
    print('Reaction Cache Locality Heuristic changed from {} to {}'.format(pre, post))

    if post >= pre:
        print('Using newly optimized reaction order')
        rxn_ordering = ordering[:]
    else:
        print('Using original reaction order')
        rxn_ordering = range(len(reacs))

    #now set up the species optimization problem
    #set up the Numberjack constraint problem to optimize reaction order
    model = Model()

    sp_order = [Variable(0, len(specs), "sp_{}".format(i)) for i in range(len(specs))]

    #set up the constraints to ensure unique rxn placement
    model.add(AllDiff(sp_order))

    #now set up score function
    score_matrix = Matrix(len(specs), len(reacs), 2)

    #set up the r_to_s constraints
    for i, sp in enumerate(sp_order):
        for rxn in range(len(reacs)):
            val = 1 if rxn in s_to_r[i] else 0
            model.add(score_matrix[sp, rxn] == val)

    score_list = []
    for i in range(len(specs) - 1):
        sp = sp_order[i]
        for j in range(len(reacs)):
            score_list.append(score_matrix[sp, j] < score_matrix[sp + 1, j])

    score = Sum(score_list)
    model.add(Minimize(score))
    solver = model.load('SCIP')
    solver.setThreadCount(multi_thread)
    solver.setVerbosity(2)
    solver.setTimeLimit(time_lim * 60)
    print(solver.solve(), solver.is_opt())
    solver.printStatistics()

    def __sp_score_check(s_to_r, ordering):
        the_sum = 0
        for i in range(len(s_to_r) - 1):
            the_sum += len(s_to_r[ordering[i + 1]].difference(s_to_r[ordering[i]]))

        return the_sum

    ordering = [x.get_value() for x in sp_order]
    assert len(set(ordering)) == len(ordering)

    pre = __sp_score_check(s_to_r, range(len(specs)))
    post = __sp_score_check(s_to_r, ordering)
    print('Species Cache Locality Heuristic changed from {} to {}'.format(pre, post))

    if post >= pre:
        print('Using newly optimized species order')
        spec_ordering[:] = ordering
    else:
        print('Using original species order')
        spec_ordering = range(len(specs))

    print_spec_order = []
    # finally reorder the spec and rxn orderings to fix for printing
    for spec_ind in range(len(spec_ordering)):
        print_spec_order.append(
            spec_ordering.index(spec_ind)
        )

    print_rxn_order = []
    for rxn_ind in range(len(rxn_ordering)):
        print_rxn_order.append(
            rxn_ordering.index(rxn_ind)
        )

    # save to avoid reoptimization if possible
    with open(build_path + 'optimized.pickle', 'wb') as file:
        pickle.dump(specs, file)
        pickle.dump(reacs, file)
        pickle.dump(rxn_rate_order, file)
        pickle.dump(spec_rate_order, file)
        pickle.dump(pdep_rate_order, file)
        pickle.dump(print_spec_order, file)
        pickle.dump(print_rxn_order, file)

    # complete, so now return
    return specs, reacs, rxn_rate_order, pdep_rate_order, spec_rate_order, print_spec_order, print_rxn_order
