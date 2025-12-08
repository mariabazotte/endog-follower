import numpy as np
import hopsy

print("hopsy version:", hopsy.__version__)

# Define polytope Ax <= b.
A = np.array([
    [-1.0, 0.0, 0.0],   # x1 >= 0
    [ 0.0,-1.0, 0.0],   # x2 >= 0
    [ 0.0, 0.0,-1.0],   # x3 >= 0
    [ 1.0, 1.0, 1.0]    # x1+x2+x3 <= 1
], dtype=float)
b = np.array([0.0, 0.0, 0.0, 1.0], dtype=float)

# ---- Models (subclass hopsy.Model) ----
class LinearModel():
    def __init__(self, c, d):
        # store params as numpy arrays
        self.c = float(c)
        self.d = np.asarray(d, dtype=float)

    def log_density(self, x):
        # x will be numpy array in sampling space
        val = self.c + float(np.dot(self.d, x))
        if val <= 0.0:
            return -np.inf
        return np.log(val)

class ExponentialModel():
    def __init__(self, a, b):
        self.a = float(a)
        self.b = np.asarray(b, dtype=float)

    def log_density(self, x):
        # log π(x) = a + b·x
        return self.a + float(np.dot(self.b, x))

# ---- choose model and create Problem(A,b,model) ----
# Example: linear model (ensure c + d·x > 0 on domain or handle -inf)
lin_model = LinearModel(c=1.0, d=[1.0, 2.0, 3.0])
# Or exponential:
exp_model = ExponentialModel(a=0.0, b=[1.0, -1.0, 0.5])

# Attach the model to the Problem (this is the key step)
# You can pass the model as third argument or with keyword, depending on your hopsy version:
problem = hopsy.Problem(A, b, model=lin_model)

# Round (optional / recommended)
problem = hopsy.round(problem)

# Proposal: you can explicitly pick one or let MarkovChain choose a good default.
proposal = hopsy.UniformHitAndRunProposal  # class or instance - either is OK in constructor below

# Starting point:
starting_point = hopsy.compute_chebyshev_center(problem)

# Markov chain definitions (note we pass the proposal to MarkovChain)
markov_chains = [
    hopsy.MarkovChain(problem, proposal, starting_point=starting_point)
    for _ in range(4)
]

# RNGs
rng = [hopsy.RandomNumberGenerator(seed=i*42) for i in range(4)]

# Sampling (same call as before)
acceptance_rates, samples = hopsy.sample(markov_chains, rng, n_samples=10000, thinning=50, n_procs=1)

# back-transform to original space if desired
for i in range(samples.shape[0]):
    samples[i] = hopsy.back_transform(problem, points=samples[i])

print("acceptance_rates:", acceptance_rates)
print("rhat:", hopsy.rhat(samples))
print("ess:", hopsy.ess(samples))
