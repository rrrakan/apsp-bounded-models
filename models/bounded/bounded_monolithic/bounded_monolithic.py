from models.configuration import *

class BoundedMonolithic(CplexModel):

    def setup_variables(self, *args):
        """
        Setup decision variables and solver parameters
        """
        mdl = self.model
        # Set solver parameters
        APSPParameters(mdl)
        # Set solver precision
        mdl.parameters.simplex.tolerances.optimality = 1e-6

        # Define bound from Theorem 1
        self.data.BOUND = self.data.NR_RESOURCES + self.data.NR_TASKS

        # Create decision vars
        self.decision_vars.s = mdl.continuous_var_dict([(t) for t in self.data.TASKS.keys()], name='s', lb=0) # Start time of each task
        self.decision_vars.mu = mdl.continuous_var(name='mu') # Period to minimize
        self.decision_vars.x = mdl.binary_var_dict([(t, r) for t in self.data.TASKS.keys() for r in self.data.TASKS[t]], name='x') # Task t is on resource r
        self.decision_vars.y = mdl.binary_var_dict([(t1, t2, k) for t1 in self.data.TASKS.keys() for t2 in self.data.TASKS.keys() if t1 != t2 for k in range(self.data.BOUND + 1)], name='y') # Repeition 1 of task t1 is before repetition m of t2


        pass

    def setup_constraints(self, *args):
        """
        Setup constraints for solving
        """
        mdl = self.model

        s = self.decision_vars.s
        mu = self.decision_vars.mu
        x = self.decision_vars.x
        y = self.decision_vars.y

        M = max([self.data.TASKS[t][r] for t in self.data.TASKS for r in self.data.TASKS[t]]) * self.data.NR_TASKS * self.data.BOUND

        # Begin constraints
        for d in self.data.DEPENDENCIES:
            # Add task dependencies constraints
            mdl.add_constraint(
                s[d[1]] >= s[d[0]] - self.data.DEPENDENCIES[d] * mu + mdl.sum(x[d[0], r] * self.data.TASKS[d[0]][r] for r in self.data.TASKS[d[0]]))

        # Add resource allocation constraint
        mdl.add_constraints(mdl.sum(x[t, p] for p in self.data.TASKS[t]) == 1 for t in self.data.TASKS)

        for t1 in self.data.TASKS.keys():
            # Add no-overlap constraint of t single tasks with another repetition.
            # Eq. (17): holds for ALL t in T (correction 3; the authors skipped
            # tasks with a zero execution time on any resource).
            mdl.add(s[t1] + mdl.sum(x[t1, r] * self.data.TASKS[t1][r] for r in self.data.TASKS[t1]) <= s[t1] + mu)
            # Add start time bound
            mdl.add(s[t1] <= self.data.BOUND * mu)

            # Add no-overlap constraints between various tasks and various repetitions
            for t2 in self.data.TASKS:
                if t1 != t2:
                    if (t1, t2) in self.data.ALLOCATION_OVERLAP:
                        for r in self.data.ALLOCATION_OVERLAP[(t1, t2)]:
                            for m1 in range(self.data.BOUND + 1):
                                mdl.add_constraint(
                                    s[t1] + self.data.TASKS[t1][r] + m1 * mu <= s[t2] + M * (
                                            3 - x[(t1, r)] - x[(t2, r)] - y[(t1, t2, m1)]),
                                )
                                mdl.add_constraint(
                                    s[t2] + self.data.TASKS[t2][r] <= s[t1] + m1 * mu + M * (
                                            2 - x[(t1, r)] - x[(t2, r)] + y[(t1, t2, m1)]),
                                )

        # Minimize mu
        mdl.minimize(mu)
