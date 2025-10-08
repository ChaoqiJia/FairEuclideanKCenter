import csv
import math
import random
import networkx as nx
import numpy as np
import copy
from scipy.spatial import distance

global SQRT3, ONE_SQRT_3
SQRT3 = math.sqrt(3)
ONE_SQRT_3 = 1 + math.sqrt(3)


########################################################################################################################
# independent center set
def ICS(Gamma_init_X, X, classTable, constraints, m, epsilon, R_dict, left, right, j_bound, metric='euclidean'):
    Gamma_x = tuple(X)  # currently point
    para = 1 + epsilon

    j = j_bound + 1
    right = max(right, np.linalg.norm(np.array(Gamma_x) - np.array(Gamma_init_X[classTable])))

    while left < (para - 1) * (para ** j) <= right:
        R = round((para - 1) * (para ** j), 3)

        if R not in R_dict:
            R_dict[R] = [[] for _ in range(m)]
        Gamma_para = R_dict[R]
        if len(Gamma_para[classTable]) == 0:
            Gamma_para[classTable].append(tuple(Gamma_init_X[classTable]))

        # independent_center_set
        if min_metric(Gamma_x, Gamma_para[classTable], metric) > R * SQRT3:
            Gamma_para[classTable].append(Gamma_x)

        if any(len(g) > 2 * sum(constraints) for g in Gamma_para):
            left = R
            j_bound = j
        j += 1

    return left, right, j_bound, R_dict


########################################################################################################################
def streamingAlg(constraints, epsilon, PATH, PATH_G, metric='euclidean'):
    # r-cover
    left = 0
    right = 0
    j = -1
    m = len(constraints)
    Gamma_init_X = [-1] * m
    R_dict = dict()
    save_center = [[] for _ in range(m)]
    csv_filename = PATH
    # group_file = open(PATH_G, 'r')
    # next(group_file)
    with open(csv_filename, 'r') as csvfile:
        csv_reader = csv.reader(csvfile)
        next(csv_reader)
        for row in csv_reader:
            Y = int(row[-1])
            # Y = int(group_file.readline().strip())
            X = np.array(row[:-1], dtype=float)
            # X = np.array(row, dtype=float)
            if isinstance(Gamma_init_X[Y], int):
                Gamma_init_X[Y] = X
                print()
                # j = math.log(left / epsilon, (1 + epsilon))

            if len(save_center[Y]) < 2 * constraints[Y]:
                save_center[Y].append(tuple(X))
            left, right, j, R_dict = ICS(Gamma_init_X, X, Y, constraints, m, epsilon, R_dict, left, right, j,
                                         metric='euclidean')

    return left, right, {k: v for k, v in R_dict.items() if left < k <= right}, save_center


########################################################################################################################

def min_metric(x, X, metric='euclidean'):
    if x and X:
        x = np.array(x)
        X = np.array(X)

        distances = np.linalg.norm(X - x, axis=1)
        return np.min(distances)
    else:
        return math.inf
    # return np.min(distance_matrix)


def max_degree_on_side(G, V):
    """
    side: 'left'/'right' or 0/1
    """

    node, deg = max(G.degree(V), key=lambda x: x[1])
    return node, deg


def cover_point(C, fd, fd_value, constraints, ll, metric = "Euclidean"):
    m = len(constraints)
    U = []
    # V = [x for x in fd_value[ll] if x not in set(C[ll])]

    V = [
        x
        for m_r in range(m) if m_r != ll
        for x in fd_value[m_r] if x not in set(C[m_r])
    ]

    # V = fd_value[ll] - {C[0] + C[1]}
    # phase II-1
    for f in fd_value[ll]:
        if min_metric(f, sum(C,[]), metric) > fd * (1 + math.sqrt(3)):
            U.append(f)

    G = nx.Graph()
    G.add_nodes_from(V, bipartite=0) 
    G.add_nodes_from(U, bipartite=1) 

    U_prime = copy.copy(U)
    for u in U:
        for v in U + V:
            u_arr = np.array(u)
            v_arr = np.array(v)

            if np.linalg.norm(v_arr - u_arr) <= ((1 + math.sqrt(3)) * fd):
                G.add_edge(v, u)
                if v in U_prime:
                    U_prime.remove(u)
    U_prime_p = [x for x in U if x not in U_prime]
    C_u = []
    for c_u in U_prime_p:
        if min_metric(c_u, C_u, metric) > fd * 2:
            C_u.append(c_u)
    C[ll] += C_u
    U = [x for x in U if all(x not in G.neighbors(node) for node in C_u)]

    if len(U) == 0 and all([len(C[i]) <= constraints[i] for i in range(len(constraints))]):
        return C, True
    print("wrong")
    return C, False


########################################################################################################################
def one_two_sqrt_three_m_Approx(constraints, epsilon, PATH, PATH_G, metric='euclidean'):
    random.seed(42)
    m = len(constraints)
    left, right, stream_rst, save_center = streamingAlg(constraints, epsilon, PATH, PATH_G)

    # post-processing
    for fd, fd_value in stream_rst.items():
        if left < fd < right:
            # ll = -1
            C = [[] for _ in range(m)]
            # phase 1
            if np.all([len(fd_value[i]) <= constraints[i] for i in range(len(constraints))]):
                C = fd_value.copy()
                break
            else:
                fd_value_copy = copy.deepcopy(fd_value)
                while sum(bool(g) for g in fd_value_copy) != 1:
                    Psi = [[] for _ in range(m)]
                    Flag_last = True
                    sorted_indices = sorted(range(m), key=lambda i: constraints[i] - len(C[i]), reverse=True)
                    for i in sorted_indices:
                        if Flag_last and fd_value_copy[i] and len(C[i]) < constraints[i]:
                            Psi[i].append(random.choice(fd_value_copy[i]))
                            Flag_last = False
                        for item in fd_value_copy[i]:
                            if len(C[i]) < constraints[i] and min_metric(item, sum(Psi, []), metric) > fd * 2:
                                Psi[i].append(item)
                                break

                    for i in range(m):
                        print(f"!: {i} ")
                        print(len(fd_value_copy[i]))
                        for item in fd_value_copy[i]:
                            if min_metric(item, sum(Psi, []), metric) > fd * 2:
                                Psi[i].append(item)

                    if np.all([len(Psi[i]) + len(C[i]) <= constraints[i] for i in range(len(constraints))]):
                        for i in range(m):
                            C[i] += Psi[i]
                        right = fd
                        print(fd)
                        print("a")
                        break
                    elif sum(len(x) for x in Psi) + sum(len(x) for x in C)> sum(constraints):
                        left = fd
                        print("b")
                        break
                    else:
                        violations = [i for i in range(len(constraints)) if len(Psi[i]) + len(C[i]) < constraints[i]]
                        for i in violations:
                            C[i] += Psi[i]
                            for ii in range(m):
                                for phi in fd_value_copy[ii]:
                                    if min_metric(phi, sum(C, []), metric) <= fd * ONE_SQRT_3:
                                        fd_value_copy[ii].remove(phi)
                if sum(bool(g) for g in fd_value_copy) == 1:
                    print([len(g) for g in fd_value_copy])
                    non_empty_indices = [i for i, g in enumerate(fd_value_copy) if bool(g)][0]
                    # phase 2
                    if non_empty_indices:
                        C,flag = cover_point(C, fd, fd_value, constraints, non_empty_indices)
                        print("Aaa")
                        if flag:
                            right = fd
             
    for num in range(m):
        C_list = [x for x in save_center[num] if x not in C[num]]
        print(constraints[num] - len(C[num]))
        C[num] += random.sample(C_list, constraints[num] - len(C[num]))
    C = sum(C, [])
    print(len(C))
    return C
