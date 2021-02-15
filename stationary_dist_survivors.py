"""
Read in a pickle from the 'data' folder containing the stationary distribution for the ARW on a connected simple graph
with one sink vertex.
Analyze the stationary distribution in terms of the probability of least k particles surviving and exactly k particles
surviving.
Output the results to the 'data' folder in pretty plaintext files.

External dependencies: a pickle, typically from stationary_dist.py
The pickle contains a list of stable states and a list of probabilities for each state.
"""
import sympy
import pickle
import os


def survivors(k, sd):
    """
    Compute the probability that at least k particles survive from the stationary distribution sd. Typically, sd is the
    output of stationary_dist() in the program stationary_dist.py.
    :param k: A number of particles.
    :type k: A nonnegative integer, typically between 0 and the number of non-sink vertices (otherwise an empty list is
    returned).
    :param sd: The stationary distribution.
    :type sd: A tuple with two elements.
    The first element is a list of stable states. Each element of this list is a list whose elements are either 0 or
    's', indicating no particle or a sleeping particle, respectively, at each vertex.
    (It is not required that 's' indicate sleeping particles, as long as it is different from 0.)
    The second element is a row vector of symbolic expressions representing probabilities that correspond to the stable
    states.
    :return: The probability that at least k particles survive.
    :rtype: A probability, determined from the stationary distribution sd.
    """
    (states, dist) = sd
    if len(states) != len(dist):
        raise ValueError('The number of states does not match the number of entries in the distribution vector.')
    t = len(states)  # number of states
    n = len(states[0])  # number of non-sink vertices

    prob = 0

    for s in range(t):
        if states[s].count(0) <= n - k:
            prob += dist[s]
    return prob


def exact_survivors(k, sd):
    """
    Compute the probability that exactly k particles survive from the stationary distribution sd. Typically, sd is the
    output of stationary_dist() in the program stationary_dist.py.
    :param k: A number of particles.
    :type k: A nonnegative integer, typically between 0 and the number of non-sink vertices (otherwise an empty list is
    returned).
    :param sd: The stationary distribution.
    :type sd: A tuple with two elements.
    The first element is a list of stable states. Each element of this list is a list whose elements are either 0 or
    's', indicating no particle or a sleeping particle, respectively, at each vertex.
    (It is not required that 's' indicate sleeping particles, as long as it is different from 0.)
    The second element is a row vector of symbolic expressions representing probabilities that correspond to the stable
    states.
    :return: The probability that exactly k particles survive.
    :rtype: A probability, determined from the stationary distribution sd.
    """
    (states, dist) = sd
    if len(states) != len(dist):
        raise ValueError('The number of states does not match the number of entries in the distribution vector.')
    t = len(states)  # number of states
    n = len(states[0])  # number of non-sink vertices

    prob = 0

    for s in range(t):
        if states[s].count(0) == n - k:
            prob += dist[s]
    return prob


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
        sd = pickle.loads(in_file.read())
    (states, dist) = sd

    out_path = os.path.join(os.path.dirname(__file__), 'data/')

    # Output the one-point joint intensities (marginals) and pair correlations.
    survivor_probs = [None] * (len(states[0]) + 1)
    with open(out_path + graph_name + '-survivors.txt', 'w') as out_file:
        for k in range(len(states[0]), -1, -1):
            survivor_probs[k] = sympy.factor(survivors(k, sd))
            out_file.write(my_pretty(survivor_probs[k]))
            print("Survivors: " + str(k))
    exact_survivor_probs = [None] * (len(states[0]) + 1)
    with open(out_path + graph_name + '-exact-survivors.txt', 'w') as out_file:
        for k in range(len(states[0]), -1, -1):
            exact_survivor_probs[k] = sympy.factor(exact_survivors(k, sd))
            out_file.write(my_pretty(exact_survivor_probs[k]))
            print("Exact survivors: " + str(k))
    univar = {k: sympy.symbols('q') for k in sympy.symbols(['q_{}'.format(x) for x in range(len(states[0]))])}
    with open(out_path + graph_name + '-survivors-univar.txt', 'w') as out_file:
        for k in range(len(states[0]) + 1):
            out_file.write(my_pretty(survivor_probs[k].subs(univar)))
            print("Survivors (univariate): " + str(k))
    with open(out_path + graph_name + '-exact-survivors-univar.txt', 'w') as out_file:
        for k in range(len(states[0]) + 1):
            out_file.write(my_pretty(exact_survivor_probs[k].subs(univar)))
            print("Exact survivors (univariate): " + str(k))
