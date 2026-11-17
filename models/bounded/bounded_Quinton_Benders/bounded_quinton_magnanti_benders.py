from combinatorial_optimization_tools import *
from .bounded_quinton_master import *
from .bounded_quinton_dual_sub import *
from .bounded_quinton_primal_sub import *
from .bounded_quinton_auxilary_sub import *
import copy
from dataclasses import replace

class BoundedQuintonMagnantiWongBenders(MagnantiWongPapadakosBendersDecomposition):
    def   __init__(self, name, model_data, paper_faithful=True):
        # Store the model data
        self.data = model_data

        # paper_faithful=True (default): correction 1 (Papadakos core point,
        # Sec. VII-D) is active. paper_faithful=False restores the authors'
        # release-1.0.2 core-point logic. epsilon stays 1e-9 in BOTH modes by
        # the recorded project decision (tau = 1/mu ~ 1e-4; an absolute gap of
        # 1e-6 would allow period errors of tens of units).
        self.paper_faithful = paper_faithful

        # Initialize master problem
        master_problem = BoundedQuintonMaster(name + "_master", model_data)
        primal_sub_problem = BoundedQuintonPrimalSub(name + "_primal_sub", model_data)
        dual_sub_problem = BoundedQuintonDualSub(name + "_dual_sub", model_data)
        auxiliary_sub_problem = BoundedQuintonMagnantiWongSub(name + "_auxiliary_sub", model_data)

        super().__init__(name, master_problem, primal_sub_problem, dual_sub_problem, auxiliary_sub_problem)

    def _MagnantiWongPapadakosBendersDecomposition__extract_master_problem_results(self):
        """
        Extracts the relevant trail values from the master problem solution
        """

        self.master_problem_extract = copy.copy(self.master_problem.decision_vars)

        self.master_problem_extract.K = {x: round(y.solution_value)  for x, y in self.master_problem_extract.K.items()}
        self.master_problem_extract.m = {x: round(y.solution_value)  for x, y in self.master_problem_extract.m.items()}

        self.master_problem_extract.z = self.master_problem_extract.z.solution_value

    def _MagnantiWongPapadakosBendersDecomposition__calculate_core_point(self, iteration):
        """
        Obtains the core point for auxiliary sub problems
        """
        if self.paper_faithful:
            self._calculate_core_point_paper(iteration)
        else:
            self._calculate_core_point_legacy(iteration)

    def _calculate_core_point_paper(self, iteration):
        """
        Papadakos core point, paper Sec. VII-D (correction 1):
        v'_0 = 1/2 v_0 for k = 0; v'_k = 1/2 v'_{k-1} + 1/2 v_k for k > 0.
        The authors' version (kept in _calculate_core_point_legacy) collapses
        to the trial point from iteration 1 on and has a `0.5 +` typo in z.
        """
        trial = self.master_problem_extract
        self.core_point = copy.copy(trial)

        if iteration == 0:
            # v'_0 = 1/2 v_0  (paper Sec. VII-D, k = 0)
            self.core_point.K = {i: 0.5 * trial.K[i] for i in trial.K}
            self.core_point.m = {i: 0.5 * trial.m[i] for i in trial.m}
            self.core_point.z = 0.5 * trial.z
        else:
            # v'_k = 1/2 v'_{k-1} + 1/2 v_k  (paper Sec. VII-D, k > 0)
            prev = self.previous_core_point
            self.core_point.K = {i: 0.5 * prev.K[i] + 0.5 * trial.K[i] for i in trial.K}
            self.core_point.m = {i: 0.5 * prev.m[i] + 0.5 * trial.m[i] for i in trial.m}
            self.core_point.z = 0.5 * prev.z + 0.5 * trial.z

        # Keep v'_k as the running average for the next iteration.
        self.previous_core_point = copy.copy(self.core_point)

    def _calculate_core_point_legacy(self, iteration):
        """
        The authors' original core-point update (release 1.0.2), kept verbatim
        for the paper_faithful=False mode added in Stage 3.
        """
        self.core_point = copy.copy(self.master_problem_extract)

        if iteration == 0:
            self.previous_core_point = copy.copy(self.core_point)

            self.core_point.K = {x: 0.5 * y for x, y in self.core_point.K.items()}
            self.core_point.m = {x: 0.5 * y for x, y in self.core_point.m.items()}
            self.core_point.z = 0.5 * self.core_point.z
        else:
            self.previous_core_point.K = copy.deepcopy(self.core_point.K)
            self.previous_core_point.m = copy.deepcopy(self.core_point.m)
            self.previous_core_point.z = copy.deepcopy(self.core_point.z)

            self.core_point.K = {x: 0.5 * self.previous_core_point.K[x] + 0.5 * self.core_point.K[x] for x in self.core_point.K.keys()}
            self.core_point.m = {x: 0.5 * self.previous_core_point.m[x] + 0.5 * self.core_point.m[x]  for x in self.core_point.m.keys()}
            self.core_point.z = 0.5 + self.previous_core_point.z + 0.5 * self.core_point.z


    def _MagnantiWongPapadakosBendersDecomposition__apply_cuts(self):
        """
        Construct and apply Benders Cuts
        """
        mdl = self.master_problem.model

        # Obtain variable solutions
        sub_problem_extract = copy.copy(self.auxiliary_sub_problem.decision_vars)

        pi_d = {x: y.solution_value for x, y in sub_problem_extract.pi_d.items()}
        pi_nr = {x: y.solution_value for x, y in sub_problem_extract.pi_nr.items()}
        pi_p = {x: y.solution_value for x, y in sub_problem_extract.pi_p.items()}
        pi_y1 = {x: y.solution_value for x, y in sub_problem_extract.pi_y1.items()}
        pi_y2 = {x: y.solution_value for x, y in sub_problem_extract.pi_y2.items()}


        K = self.master_problem.decision_vars.K
        m = self.master_problem.decision_vars.m
        z = self.master_problem.decision_vars.z

        # P1 as in Quinton et al. 2020
        P1 = max([self.master_problem.data.TASKS[a][r] for a in self.data.TASKS for r in self.data.TASKS[a]])

        # Traditional Benders optimality cut
        cut = (
                z <= (mdl.sum(pi_nr[i] * mdl.sum(m[(i, r)] / self.data.TASKS[i][r] for r in self.data.TASKS[i]) for i in
                              self.data.TASKS) +
                      mdl.sum(
                    pi_d[(ij[0], ij[1], r)] * (P1 * (2 - m[(ij[0], r)] - m[(ij[1], r)]) + K[(ij[0], ij[1])])
                    for ij in self.data.ALLOCATION_OVERLAP.keys() for r in self.data.ALLOCATION_OVERLAP[ij]) +
                      mdl.sum(pi_p[c] * self.data.DEPENDENCIES[c] for c in self.data.DEPENDENCIES) +
                      mdl.sum(pi_y2[(i, r)] * m[(i, r)] for i in self.data.TASKS for r in self.data.TASKS[i]))
        )

        # Add cut to master problem
        self.master_problem.model.add(cut)
