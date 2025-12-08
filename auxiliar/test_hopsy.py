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

# Problem definition and rounding.
problem = hopsy.Problem(A, b)
problem = hopsy.round(problem)

# Proposal definition.
proposal = hopsy.UniformHitAndRunProposal

# Starting point definition.
starting_point = hopsy.compute_chebyshev_center(problem)

# Markov chain definitions.
markov_chains = [hopsy.MarkovChain(problem, proposal, starting_point=starting_point) for i in range(4)]

# Random number generation.
rng = [hopsy.RandomNumberGenerator(seed=i*42) for i in range(4)]

# Sampling.
acceptance_rates, samples = hopsy.sample(markov_chains, rng, n_samples=10000, thinning=50, n_procs=1)

print(samples.shape)

for i in range(samples.shape[0]):  # loop over chains
    samples[i] = hopsy.back_transform(problem, points=samples[i])

# Print.
print(acceptance_rates)

# checks convergence and ESS
rhat = hopsy.rhat(samples)
print(rhat)
assert((rhat<=1.01).all()) # asserts that convergence has been reached. Here we use a strict limit of 1.01

ess = hopsy.ess(samples)
print(ess)
assert((ess >= 400).all()) # asserts that we have reached the minimum number of effective samples we want, in this case 400.