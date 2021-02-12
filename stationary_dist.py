"""
Compute the stationary distribution for the ARW on a connected simple graph with one sink vertex.
Save the result in a pickle file and in 3 plaintext files in the 'data' folder.
One plaintext file ends in '-states.txt' and contains the stable states.
Another plaintext file ends in '-distribution-latex.txt' and contains the LaTeX expressions for the probabilities that
correspond to the stable states.
The last plaintext file ends in '-distribution.txt' and contains a pretty version of the probabilities that
correspond to the stable states.

External dependencies: solver.py.
"""
import sympy
from collections import deque
import copy
import solver
import time
import pickle
import os


def stationary_dist(a, sleep_probs):
    """
    Compute the stationary distribution of the ARW on a connected simple graph with one sink vertex.
    The initial state of the ARW consists of one active particle at each non-sink vertex.
    :param a: Adjacency list of the graph.
    :type a: List of lists of integers. The vth list (counting from 0) contains the neighbors of vertex v. The last
    vertex is assumed to be the sink vertex.
    :param sleep_probs: The sleep probabilities at each non-sink vertex.
    :type sleep_probs: List. The elements of the list are symbolic variables representing the sleep probabilities. The
    size of the list should be exactly one fewer than the number of vertices of the graph.
    :return: A subset of the absorbing states and the probability of ending up at each one.
    :rtype: A tuple with two elements.
    The first element is a list of stable states. Each element of this list is a list whose elements are either 0 or
    's', indicating no particle or a sleeping particle, respectively, at each vertex.
    The second element is a row vector of symbolic expressions representing probabilities that correspond to the stable
    states.
    """
    if len(a) - 1 != len(sleep_probs):
        raise ValueError('There should be exactly one sink vertex.')

    # n is the number of non-sink vertices
    n = len(a) - 1
    # deg is a list of degrees of the non-sink vertices
    deg = [len(a[v]) for v in range(n)]

    # t is a list of states in the transition matrix.
    # Initialize t with the state where all non-sink vertices have 1 active particle.
    t = [[1] * n]
    # print(t)
    t_absorb_idx = []  # list of indices of t that correspond to absorbing states
    # m is a Markov transition matrix between states, implemented as a list of lists in row-major order.
    m = [[0]]
    # q is a queue of current states to check.
    # Initialize q with the state where non-sink vertices have 1 active particle.
    # Each state in q includes at the end its index in t.
    init_state_q = [1] * n
    init_state_q.append(0)
    q = deque([init_state_q])
    while q:
        state = q.popleft()
        # Pick the vertex v to fire.
        # Fire the (only) vertex with 2 active particles; otherwise, fire the first vertex with 1 active particle.
        v = 0
        while v < n:
            if state[v] == 2:
                break
            v += 1
        if v == n:
            v = 0
            while v < n:
                if state[v] == 1:
                    break
                v += 1
        if v == n:
            # The state is absorbing.
            t_absorb_idx.append(state[n])
            m[state[n]][state[n]] = 1
        else:
            if state[v] == 2:
                # The vertex v has 2 active particles.
                temp_state = copy.copy(state)
                t_idx = temp_state.pop()

                # Consider the active particle trying to fall asleep (which leaves the state unchanged).
                m[t_idx][t_idx] = sleep_probs[v]

                temp_state[v] = 1
                # Consider the active particle at v jumping to a neighbor.
                for i in a[v]:
                    if i != n:
                        # The active particle at v jumps to a non-sink vertex.
                        temp_state_i_old = temp_state[i]
                        if temp_state[i] == 's':
                            temp_state[i] = 2
                        else:
                            temp_state[i] += 1
                        try:
                            # Check if we transition to a state that we have seen before (in t).
                            new_idx = t.index(temp_state)
                        except ValueError:
                            # We transition to a new state.
                            new_state = copy.copy(temp_state)
                            t.append(new_state)
                            new_idx = len(t) - 1
                            new_state_q = copy.copy(new_state)
                            new_state_q.append(new_idx)
                            q.append(new_state_q)
                            for row in m:
                                row.append(0)
                            m.append([0] * len(t))
                        m[t_idx][new_idx] = (1 - sleep_probs[v]) / deg[v]
                        temp_state[i] = temp_state_i_old
                    else:
                        # The active particle at v jumps to the sink.
                        try:
                            new_idx = t.index(temp_state)
                        except ValueError:
                            new_state = copy.copy(temp_state)
                            t.append(new_state)
                            new_idx = len(t) - 1
                            new_state_q = copy.copy(new_state)
                            new_state_q.append(new_idx)
                            q.append(new_state_q)
                            for row in m:
                                row.append(0)
                            m.append([0] * len(t))
                        m[t_idx][new_idx] = (1 - sleep_probs[v]) / deg[v]

            if state[v] == 1:
                # The vertex v has 1 active particle.
                temp_state = copy.copy(state)
                t_idx = temp_state.pop()

                # Consider the active particle at v falling asleep.
                temp_state[v] = 's'
                try:
                    new_idx = t.index(temp_state)
                except ValueError:
                    new_state = copy.copy(temp_state)
                    t.append(new_state)
                    new_idx = len(t) - 1
                    new_state_q = copy.copy(new_state)
                    new_state_q.append(new_idx)
                    q.append(new_state_q)
                    for row in m:
                        row.append(0)
                    m.append([0] * len(t))
                m[t_idx][new_idx] = sleep_probs[v]

                temp_state[v] = 0
                # Consider the active particle at v jumping to a neighbor.
                for i in a[v]:
                    if i != n:
                        temp_state_i_old = temp_state[i]
                        if temp_state[i] == 's':
                            temp_state[i] = 2
                        else:
                            temp_state[i] += 1
                        try:
                            new_idx = t.index(temp_state)
                        except ValueError:
                            new_state = copy.copy(temp_state)
                            t.append(new_state)
                            new_idx = len(t) - 1
                            new_state_q = copy.copy(new_state)
                            new_state_q.append(new_idx)
                            q.append(new_state_q)
                            for row in m:
                                row.append(0)
                            m.append([0] * len(t))
                        m[t_idx][new_idx] = (1 - sleep_probs[v]) / deg[v]
                        temp_state[i] = temp_state_i_old
                    else:
                        try:
                            new_idx = t.index(temp_state)
                        except ValueError:
                            new_state = copy.copy(temp_state)
                            t.append(new_state)
                            new_idx = len(t) - 1
                            new_state_q = copy.copy(new_state)
                            new_state_q.append(new_idx)
                            q.append(new_state_q)
                            for row in m:
                                row.append(0)
                            m.append([0] * len(t))
                        m[t_idx][new_idx] = (1 - sleep_probs[v]) / deg[v]

    # From m, t, and t_absorb_idx, we can calculate the probabilities of ending up at each absorbing state.
    # the transition matrix between transient states
    m_trans_trans = \
        sympy.Matrix([[row[j] for j in range(len(m)) if j not in t_absorb_idx]
                      for row in [m[i] for i in range(len(m)) if i not in t_absorb_idx]])
    # the transition matrix from transient states to absorbing states
    m_trans_absorb = \
        sympy.Matrix([[row[j] for j in t_absorb_idx]
                      for row in [m[i] for i in range(len(m)) if i not in t_absorb_idx]])

    t0 = time.process_time()
    print("Time to compute transition matrix: " + str(t0))

    # The probabilities of ending up at each absorbing state are given by the row vector
    # ((I - m_trans_trans)^T \ e1)^T * m_trans_absorb
    # The linear solve is the most time-consuming part of the program.
    # We use solver() (which seems to be faster than sympy.linsolve() and sympy.solve()):
    ell = len(t) - len(t_absorb_idx)
    mat = sympy.SparseMatrix(ell, ell, {})
    for i in range(ell):
        mat[i, i] = 1
    for i in range(ell):
        for j in range(ell):
            if m_trans_trans[i, j] != 0:
                mat[i, j] -= m_trans_trans[i, j]
    dist = solver.inverse(mat, [0], list(range(ell))) * m_trans_absorb

    t1 = time.process_time() - t0
    print("Time to compute final answer: " + str(t1))

    return [t[i] for i in t_absorb_idx], dist


if __name__ == "__main__":
    graph_name = "4-clique"  # Change the name as necessary!
    a = [[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]]
    # Possible graphs include:
    # a = [[1], [0]]  # 2-clique, a.k.a. 2-path, a.k.a. 2-cycle
    # a = [[1], [0, 2], [1]]  # 3-path
    # a = [[1, 2], [0, 2], [0, 1]]  # 3-clique, a.k.a. 3-cycle
    # a = [[1], [0, 2], [1, 3], [2]]  # 4-path
    # a = [[1, 3], [0, 2], [1, 3], [0, 2]]  # 4-cycle
    # a = [[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]]  # 4-clique
    sleep_probs = sympy.symbols(['q_{}'.format(x) for x in range(len(a)-1)])

    sd = stationary_dist(a, sleep_probs)
    out_path = os.path.join(os.path.dirname(__file__), 'data/')
    with open(out_path + graph_name + '.pickle', 'wb') as out_file:
        out_file.write(pickle.dumps(sd))
    with open(out_path + graph_name + '-states.txt', 'w') as out_file:
        for state in sd[0]:
            out_file.write(str(state) + "\n")
    with open(out_path + graph_name + '-distribution-latex.txt', 'w') as out_file:
        for prob in sd[1]:
            out_file.write(sympy.latex(prob) + "\n")
    with open(out_path + graph_name + '-distribution.txt', 'w') as out_file:
        for prob in sd[1]:
            num, denom = sympy.fraction(sympy.factor(prob))
            out_file.write(sympy.pretty(num) + "\n")
            out_file.write("/\n")
            out_file.write(sympy.pretty(denom) + "\n\n\n\n\n")
