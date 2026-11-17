from combinatorial_optimization_tools import *
from models.configuration import *

class BoundedQuintonPrimalSub(CplexModel):

    def setup_variables(self):
        """
        Setup decision variables and solver parameters
        """
        mdl = self.model
        # Set solver parameters
        APSPParameters(mdl)
        # Set solver precision
        mdl.parameters.simplex.tolerances.optimality = 1e-9

        # Create decision vars
        self.decision_vars.u = mdl.continuous_var_dict([(a) for a in self.data.TASKS], name='u', lb=0)  # Linearized start time of each task
        self.decision_vars.tau = mdl.continuous_var(name='tau')  # Linearized period to minimize
        self.decision_vars.y = mdl.continuous_var_dict([(a, r) for a in self.data.TASKS for r in self.data.TASKS[a]], name='y', lb=0)  # Task allocated on a resource equals tau

        pass

    def setup_constraints(self, data):
        mdl = self.model

        u = self.decision_vars.u
        tau = self.decision_vars.tau
        y = self.decision_vars.y

        # Obtain master problem solution
        m_bar = data.m
        K_bar = data.K

        for a in self.data.TASKS:
            if 0 not in self.data.TASKS[a].values():
                # Constraint (9a) in Quinton et al. 2020
                mdl.add(self.decision_vars.tau <= mdl.sum((m_bar[(a, r)] / self.data.TASKS[a][r]) for r in self.data.TASKS[a]))

            # Constraint (9d) in Quinton et al. 2020
            mdl.add(mdl.sum(y[(a, r)] for r in self.data.TASKS[a]) == tau)

            for r in self.data.TASKS[a]:
                # Constraint (9e) in Quinton et al. 2020
                mdl.add(y[(a, r)] <= m_bar[(a, r)])

        for c in self.data.DEPENDENCIES:
            # Constraint (9b) in Quinton et al. 2020
            mdl.add(u[c[1]] + self.data.DEPENDENCIES[c] >= u[c[0]] + mdl.sum(y[(c[0], r)] * self.data.TASKS[c[0]][r] for r in self.data.TASKS[c[0]]))

        # P1 as in Quinton et al. 2020
        P1 = max([self.data.TASKS[a][r] for a in self.data.TASKS for r in self.data.TASKS[a]])

        for ij in self.data.ALLOCATION_OVERLAP.keys():
            for r in self.data.ALLOCATION_OVERLAP[ij]:
                # Constraint (9c) in Quinton et al. 2020
                mdl.add(P1 * (2 - m_bar[(ij[0], r)] - m_bar[(ij[1],r)]) + u[ij[1]] + K_bar[ij[0], ij[1]] >= u[ij[0]] + y[(ij[0], r)] * self.data.TASKS[ij[0]][r])

        # Objective as in Quinton et al. 2020
        mdl.maximize(tau)

        self.m_bar = m_bar
        self.K_bar = K_bar
        pass

    def additional_output_information(self):
        if self.decision_vars.tau.solution_value != 0:
            for actor, start_time in self.decision_vars.u.items():
                print(f"Actor: {actor}, Start Time: {start_time.solution_value * (1 / self.decision_vars.tau.solution_value)}")
            print(f"Period: {(1 / self.decision_vars.tau.solution_value)}")
