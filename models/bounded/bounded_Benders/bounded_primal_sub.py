from models.configuration import *

class BoundedPrimalSub(CplexModel):

    def setup_variables(self):
        """
        Setup decision variables and solver parameters
        """
        mdl = self.model
        APSPParameters(mdl)
        # Set solver precision
        mdl.parameters.simplex.tolerances.optimality = 1e-9

        self.data.BOUND = self.data.NR_RESOURCES + self.data.NR_TASKS

        # Create decision vars
        self.decision_vars.s = mdl.continuous_var_dict([(a) for a in self.data.TASKS.keys()], name='s', lb=0) # Start time of each actor
        self.decision_vars.mu = mdl.continuous_var(name='mu', lb=0)
        pass

    def setup_constraints(self, data):
        """
        Setup constraints for solving
        """
        mdl = self.model
        NR_PERIODS = self.data.NR_RESOURCES + self.data.NR_TASKS

        s = self.decision_vars.s
        mu = self.decision_vars.mu

        x_bar = data.x
        y_bar = data.y

        M = max([self.data.TASKS[a][r] for a in self.data.TASKS for r in self.data.TASKS[a]]) * self.data.NR_TASKS * self.data.BOUND

        # Begin model
        for c in self.data.DEPENDENCIES:
            mdl.add_constraint(
                s[c[1]] >= s[c[0]] - self.data.DEPENDENCIES[c] * mu + mdl.sum(x_bar[c[0], p] * self.data.TASKS[c[0]][p] for p in self.data.TASKS[c[0]]))

        for a1 in self.data.TASKS.keys():
            # Eq. (30): holds for ALL t in T (correction 3; the authors skipped
            # tasks with a zero execution time on any resource).
            mdl.add(s[a1] + mdl.sum(x_bar[a1, p] * self.data.TASKS[a1][p] for p in self.data.TASKS[a1]) <= s[a1] + mu)
            mdl.add(s[a1] <= self.data.BOUND * mu)

        for ij in self.data.ALLOCATION_OVERLAP.keys():
            for p in self.data.ALLOCATION_OVERLAP[ij]:
                for m1 in range(self.data.BOUND + 1):
                        mdl.add_constraint(
                            s[ij[0]] + self.data.TASKS[ij[0]][p] + m1 * mu <= s[ij[1]] + M * (
                                    3 - x_bar[(ij[0], p)] - x_bar[(ij[1], p)] - y_bar[(ij[0], ij[1], m1)]),
                        )
                        mdl.add_constraint(
                            s[ij[1]] + self.data.TASKS[ij[1]][p] <= s[ij[0]] + m1 * mu + M * (
                                    2 - x_bar[(ij[0], p)] - x_bar[(ij[1], p)] + y_bar[(ij[0], ij[1], m1)]),
                        )

        mdl.minimize(mu)

        self.x_bar = x_bar
        self.y_bar = y_bar
        pass
