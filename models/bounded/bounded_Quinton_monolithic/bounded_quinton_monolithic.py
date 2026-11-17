from combinatorial_optimization_tools import *
from models.configuration import *

class BoundedQuintonMonolithic(CplexModel):

    def setup_variables(self, *args):
        """
        Setup decision variables and solver parameters
        """
        mdl = self.model
        # Set solver parameters
        APSPParameters(mdl)
        # Set solver precision
        mdl.parameters.simplex.tolerances.optimality = 1e-6

        self.data.BOUND = self.data.NR_RESOURCES + self.data.NR_TASKS

        # Create decision vars
        self.decision_vars.u = mdl.continuous_var_dict([(a) for a in self.data.TASKS], name='u', lb=0)  # Linearized start time of each task
        self.decision_vars.tau = mdl.continuous_var(name='tau', lb=0) # Linearized period to minimize
        self.decision_vars.K = mdl.integer_var_dict([(a1, a2) for a1 in self.data.TASKS for a2 in self.data.TASKS if a1 != a2], name='K',
                                                    lb=-self.data.BOUND, ub=self.data.BOUND)   # Repetition distance between two tasks
        self.decision_vars.m = mdl.binary_var_dict([(a, r) for a in self.data.TASKS for r in self.data.TASKS[a]], name='m')  # Task allocated on a resource
        self.decision_vars.y = mdl.continuous_var_dict([(a, r) for a in self.data.TASKS for r in self.data.TASKS[a]], name='y', lb=0)  # Task allocated on a resource equals tau

        pass

    def setup_constraints(self, *args):
        """
        Setup constraints for solving
        """
        mdl = self.model

        u = self.decision_vars.u
        tau = self.decision_vars.tau
        K = self.decision_vars.K
        m = self.decision_vars.m
        y =self.decision_vars.y

        for a in self.data.TASKS:
            if 0 not in self.data.TASKS[a].values():
                # Constraint (6a) in Quinton et al. (2020)
                mdl.add(tau <= mdl.sum((m[(a, r)] / self.data.TASKS[a][r]) for r in self.data.TASKS[a]))

            # Constraint (6d) in Quinton et al. (2020)
            mdl.add(mdl.sum(y[(a, r)] for r in self.data.TASKS[a]) == tau)

            # Constraint (6g) in Quinton et al. (2020)
            mdl.add(mdl.sum(m[(a, r)] for r in self.data.TASKS[a]) == 1)

            for r in self.data.TASKS[a]:
                # Constraint (6e) in Quinton et al. (2020)
                mdl.add(y[(a, r)] <= m[(a, r)])
                # Constraint (6f) in Quinton et al. (2020)
                mdl.add(y[(a, r)] >= 0)

        for c in self.data.DEPENDENCIES:
            # Constraint (6b) in Quinton et al. (2020)
            mdl.add(u[c[1]] + self.data.DEPENDENCIES[c] >= u[c[0]] + mdl.sum(y[(c[0], r)] * self.data.TASKS[c[0]][r] for r in self.data.TASKS[c[0]]))


        P1 = max([self.data.TASKS[a][r] for a in self.data.TASKS for r in self.data.TASKS[a]])

        for ij in self.data.ALLOCATION_OVERLAP.keys():
            for r in self.data.ALLOCATION_OVERLAP[ij]:
                # Constraint (6c) in Quinton et al. (2020)
                mdl.add(P1 * (2 - m[(ij[0], r)] - m[(ij[1],r)]) + u[ij[1]] + K[ij[0], ij[1]] >= u[ij[0]] + y[(ij[0], r)] * self.data.TASKS[ij[0]][r])

            # Constraint (6h) in Quinton et al. (2020)
            mdl.add(K[ij[0], ij[1]] + K[ij[1], ij[0]] == 1)


        mdl.maximize(tau)

        pass

    def _OptimizationModel__additional_output_information(self):
        """
        Prints additional solution information
        """

        # Convert linearized start times back to regular start times
        for task, start_time in self.decision_vars.u.items():
            print(f"Task: {task}, Start Time: {start_time.solution_value * (1 / self.decision_vars.tau.solution_value)}")
        print(f"Period: {(1 / self.decision_vars.tau.solution_value)}")
