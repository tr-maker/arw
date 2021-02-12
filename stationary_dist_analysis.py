"""
Read in a pickle containing the stationary distribution for the ARW on a connected simple graph with one sink vertex.
Analyze the stationary distribution.

External dependencies: a pickle, typically from stationary_dist.py
The pickle contains a list of stable states and a list of probabilities for each state.
"""
import sympy
import pickle
from itertools import combinations
import os


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


def my_pretty(frac):
    """
    Write a symbolic rational function as a pretty string that can be printed.
    :param frac: A rational function.
    :type frac: A symbolic rational function.
    :return: A pretty version of the input.
    :rtype: A string.
    """
    num, denom = sympy.fraction(sympy.factor(frac))
    return sympy.pretty(num) + "\n" + "/\n" + sympy.pretty(denom) + "\n\n\n\n\n"


if __name__ == "__main__":
    graph_name = "4-clique"  # Change the name as necessary
    in_path = os.path.join(os.path.dirname(__file__), 'data/')
    with open(in_path + graph_name + '.pickle', 'rb') as in_file:
        sd = pickle.loads(inf.read())
    (states, dist) = sd

    out_path = os.path.join(os.path.dirname(__file__), 'data/')

    # Output the stationary distribution when all sleep rates are the same.
    univar = {k: sympy.symbols('q') for k in sympy.symbols(['q_{}'.format(x) for x in range(len(states[0]))])}
    dist_univar = dist.subs(univar)
    with open(out_path + graph_name + '-distribution-univar.txt', 'w') as out_file:
        for prob in dist_univar:
            out_file.write(my_pretty(prob))

    # Output the one-point joint intensities (marginals) and pair correlations.
    marginals = joint_intensities(1, sd)
    with open(out_path + graph_name + '-marginals.txt', 'w') as out_file:
        for marginal in marginals:
            marginal = sympy.factor(marginal)
            out_file.write(my_pretty(marginal))
    joints = joint_intensities(2, sd)
    correlations = [None] * len(joints)
    with open(out_path + graph_name + '-correlations.txt', 'w') as out_file:
        k = 0
        for i in range(len(states[0])):
            for j in range(i + 1, len(states[0])):
                correlations[k] = sympy.factor(joints[k] - marginals[i] * marginals[j])
                out_file.write(my_pretty(correlations[k]))
                print("Correlations: finished entry " + str(k))
                k += 1
    # Output the one-point joint intensities (marginals) and pair correlations when all sleep rates are the same.
    marginals_univar = [None] * len(marginals)
    for i in range(len(marginals)):
        marginals_univar[i] = marginals[i].subs(univar)
    with open(out_path + graph_name + '-marginals-univar.txt', 'w') as out_file:
        for marginal in marginals_univar:
            out_file.write(my_pretty(marginal))
    correlations_univar = [None] * len(correlations)
    for i in range(len(correlations)):
        correlations_univar[i] = correlations[i].subs(univar)
    with open(out_path + graph_name + '-correlations-univar.txt', 'w') as out_file:
        for correlation in correlations_univar:
            out_file.write(my_pretty(correlation))
