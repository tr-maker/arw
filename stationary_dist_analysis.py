"""
Read in a pickle containing the stationary distribution for the ARW on a connected simple graph with one sink vertex.
Print and analyze the stationary distribution.

External dependencies: a pickle, typically from stationary_dist.py
The pickle containing a list of stable states and a list of probabilities for each state.
"""
import sympy
import pickle
from itertools import combinations


def joint_intensities(k, sd):
    """
    Compute the k-point joint intensities (marginals) of the stationary distribution sd. Typically, sd is the output of
    stationary_dist() in the program stationary_dist.py.
    :param k: The joint intensity parameter.
    :type k: A nonnegative integer, typically between 0 and the number of non-sink vertices (otherwise an empty list is
    returned).
    :param sd: The stationary distribution.
    :type sd: A tuple with two elements.
    The first element is a list of stable states. Each element of this list is a list whose elements are either 0 or
    's', indicating no particle or a sleeping particle, respectively, at each vertex.
    (It is not required that 's' indicate sleeping particles, as long as it is different from 0.)
    The second element is a row vector of symbolic expressions representing probabilities that correspond to the stable
    states.
    :return: The k-point joint intensities.
    :rtype: A list with n choose k elements, where n is the number of non-sink vertices. Each element of the list
    is a symbolic expression representing the probability that a specific k-vertex subset of the non-sink vertices has
    all sleeping particles. The subsets are ordered in alphabetical order according to the ordering of the non-sink
    vertices in the list of stable states in sd.
    """
    (states, dist) = sd
    if len(states) != len(dist):
        raise ValueError('The number of states does not match the number of entries in the distribution vector.')
    t = len(states)  # number of states
    n = len(states[0])  # number of non-sink vertices
    joint_int = []
    for c in combinations(range(n), k):
        intensity = 0
        for s in range(t):
            c_in_s = True
            for idx in c:
                if states[s][idx] == 0:
                    c_in_s = False
                    break
            if c_in_s:
                intensity += dist[s]
        # sympy.simplify(intensity)
        joint_int.append(intensity)
    return joint_int


# Read in the computation of the stationary distribution from a pickle.
# Change the name of the pickle as necessary!
with open('4-clique.pickle', 'rb') as inf:  # Change the name of the pickle as necessary
    sd = pickle.loads(inf.read())
(states, dist) = sd

# Print the stationary distribution.
sympy.pprint(states)
for s in dist:
    # num, denom = sympy.fraction(s)
    # sympy.pprint(sympy.factor(num))
    # sympy.pprint(sympy.factor(denom))
    # print(sympy.latex(sympy.factor(num)))
    # print(sympy.latex(sympy.factor(denom)))
    print(sympy.latex(s))
    print()
print()

# The following example shows that the stationary distribution is not a determinantal process in the states.
# It requires a graph with 4 vertices.
dist = dist.subs([(sympy.symbols('q_0'), 0.5), (sympy.symbols('q_1'), 0.5), (sympy.symbols('q_2'), 0.5)])
sd = (states, dist)
dpp_mat = sympy.zeros(3)
# The diagonal entries of the kernel are determined by the 1-point joint intensities.
[dpp_mat[0, 0], dpp_mat[1, 1], dpp_mat[2, 2]] = joint_intensities(1, sd)
# The off-diagonal entries of the kernel are determined up to sign by the 2-point joint intensities.
# We pick all positive signs first.
joint_int_2 = joint_intensities(2, sd)
[dpp_mat[0, 1], dpp_mat[0, 2], dpp_mat[1, 2]] = [sympy.sqrt(dpp_mat[0, 0] * dpp_mat[1, 1] - joint_int_2[0]),
                                                 sympy.sqrt(dpp_mat[0, 0] * dpp_mat[2, 2] - joint_int_2[1]),
                                                 sympy.sqrt(dpp_mat[1, 1] * dpp_mat[2, 2] - joint_int_2[2])]
[dpp_mat[1, 0], dpp_mat[2, 0], dpp_mat[2, 1]] = [dpp_mat[0, 1], dpp_mat[0, 2], dpp_mat[1, 2]]
# sympy.pprint(dpp_mat[0, 0] * dpp_mat[1, 1] * dpp_mat[2, 2] - dpp_mat[0, 0] * dpp_mat[1, 2] * dpp_mat[2, 1]
#             - dpp_mat[0, 1] * dpp_mat[1, 0] * dpp_mat[2, 2] + dpp_mat[0, 1] * dpp_mat[1, 2] * dpp_mat[2, 0]
#             + dpp_mat[0, 2] * dpp_mat[1, 0] * dpp_mat[2, 1] - dpp_mat[0, 2] * dpp_mat[1, 1] * dpp_mat[2, 0])
# For the stationary distribution to be determinantal, the determinant of the 3x3 kernel must equal the 3-point joint
# intensity.
sympy.pprint(sympy.N(dpp_mat.det()))
sympy.pprint(sympy.N(joint_intensities(3, sd)[0]))
# We pick a different choice of signs for the off-diagonal entries.
# For determinant of a symmetric 3x3 matrix, it does not matter which off-diagonal entry (and its symmetric
# counterpart) we make negative. Furthermore, making any two off-diagonal entries (and their symmetric counterparts)
# negative is equivalent to the original matrix, and making all off-diagonal entries negative is equivalent to making
# just one off-diagonal entry (and its symmetric counterpart) negative.
[dpp_mat[0, 1], dpp_mat[0, 2], dpp_mat[1, 2]] = [-sympy.sqrt(dpp_mat[0, 0] * dpp_mat[1, 1] - joint_int_2[0]),
                                                 sympy.sqrt(dpp_mat[0, 0] * dpp_mat[2, 2] - joint_int_2[1]),
                                                 sympy.sqrt(dpp_mat[1, 1] * dpp_mat[2, 2] - joint_int_2[2])]
# For the stationary distribution to be determinantal, the determinant of the 3x3 kernel must equal the 3-point joint
# intensity.
sympy.pprint(sympy.N(dpp_mat.det()))
sympy.pprint(sympy.N(joint_intensities(3, sd)[0]))
