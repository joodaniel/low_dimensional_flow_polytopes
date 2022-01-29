︠0cb5cd15-2d6f-4f3c-9f38-eb948bcdf409︠
'''
This file contains scripts and an example calculation of Mátyás Domokos and Dániel Joó created for the paper "LOW DIMENSIONAL FLOW POLYTOPES AND THEIR TORIC IDEALS" (https://arxiv.org/abs/2105.04004), which contains explanations for the terminology used in the docstrings here.
For questions contact the authors at: jood.mollia@gmail.com
'''

def make_eqs(A, theta):
    '''
    Calculates the equations and inequalities that define a quiver polyehdron, and formats them such that they are appropriate for the constructor of the Polyhedron() class.

    Parameters:
    A : the signed adjacency matrix of the quiver
    theta : integer weight

    Returns:

    list of tuples reperesenting equations as tuple[0] = tuple[1:] * x
    list of tuples reperesenting inequalities as tuple[0] <= tuple[1:] * x
    '''
    eqns = [tuple([-1 * theta[idx]] + A[idx]) for idx in range(len(A))]
    def unitv(idx, n):
        return [1 if i == idx else 0 for i in range(n)]
    ieqs = [tuple([0] + unitv(idx, len(A[0]))) for idx in range(len(A[0]))]
    return ieqs, eqns
def edges_to_quiver(edges):
    '''
    From the list of edges of an unsigned multi-graph computes the signed adjecency matrix of the quiver obtained by placing a valency 2 sink on each of the edges. Assumes that each vertex of the graph is incident to at least one edge and that the vertices are numbered 0,...,n.

    Parameters:
    edges (list): list of 2-tuples consisting of the endpoints of each edge
    Returns:
    integer matrix : signed adjecency matrix of the quiver, with rows corresponding to vertices, columns to arrows, values are -1 on the source, +1 on the target and 0 elsewhere
    '''
    num_vertices = max(sum(edges, [])) + 1
    A = [[0 for _ in range(2 * len(edges))] for _ in range(num_vertices + len(edges))]
    for idx in range(len(edges)):
        A[edges[idx][0]][idx] = -1
        A[edges[idx][1]][idx + len(edges)] = -1
        A[num_vertices + idx][idx] = 1
        A[num_vertices + idx][idx + len(edges)] = 1
    return A
def pol(A, w):
    '''Calculates the quiver polyhedron of a quiver and a weight
    Parameters:
    A (integer matrix) : the signed adjecency matrix of the quiver
    w (list of ints): integer weight

    Returns:
    Polyhedron: the quiver polyhedron
    '''
    ie, eq = make_eqs(A, w)
    return Polyhedron(ieqs=ie, eqns=eq)

def tight_list(A, w):
    '''
    Given a quiver and a weight returns the list of arrows that are not contractible nor removable. The idea is that an arrow is contractible if and only if the dimension of the quiver polytope drops after reversing its direction.
    Parameters:
    A (integer matrix) : the signed adjecency matrix of the quiver
    w (list of ints): integer weight

    Returns:
    list of int: indices of non contractible arrows
    '''
    dim = pol(A, w).dim()
    tight = []
    eq = [tuple([-1 * w[idx]] + A[idx]) for idx in range(len(A))]
    def unitv(idx, n, val=1):
        return [val if i == idx else 0 for i in range(n)]
    ieqs = [tuple([0] + unitv(idx, len(A[0]))) for idx in range(len(A[0]))]
    for idx in range(len(A[0])):
        ie = [tuple([0] + unitv(idx2, len(A[0]), val = -1 if idx==idx2 else 1)) for idx2 in range(len(A[0]))]

        P = Polyhedron(ieqs=ie, eqns=eq)
        if P.dim() == dim:
            tight.append(idx)
    return tight

def contract(A, w, idx):
    ''' Contracts an arrow in a quiver and calculates the new weight of the contracted quiver
    Parameters:
    A (integer matrix) : the signed adjecency matrix of the quiver
    w (list of ints): integer weight
    idx (int): index of arrow to be contracted

    Returns:
    integer matrix:  the signed adjecency matrix of the new quiver
    list of ints: the new weight
    '''
    A = deepcopy(A)
    w = list(deepcopy(w))
    source = next((j for j, x in enumerate([A[i][idx] for i in range(len(A))]) if x == -1), None)
    target = next((j for j, x in enumerate([A[i][idx] for i in range(len(A))]) if x == 1), None)
    if source == None:
        print(A)
        print(w)
        print(idx)
    for i in range(len(A)):
        A[i].pop(idx)
    for i in range(len(A[0])):
        A[target][i]+= A[source][i]

    w[target] += w[source]
    A.pop(source)
    w.pop(source)
    return A, w

def list_arrows(A):
    return [[next((j for j, x in enumerate([A[i][idx] for i in range(len(A))]) if x == -1), None), next((j for j, x in enumerate([A[i][idx] for i in range(len(A))]) if x == 1), None)] for idx in range(len(A[0]))]

def tighten(A, w, dim=None):
    '''
    Given a quiver and a weight returns a tight quiver with the same quiver polytope. The idea is that an arrow is contractible if and only if the dimension of the quiver polytope drops after reversing its direction. It is assumed that the quiver contains no removable arrows.
    Parameters:
    A (integer matrix) : the signed adjecency matrix of the quiver
    w (list of ints): integer weight

    Returns:
    integer matrix:  the signed adjecency matrix of the new quiver
    list of ints: the new weight
    '''
    if dim == None:
        dim = pol(A, w).dim()
    eq = [tuple([-1 * w[idx]] + A[idx]) for idx in range(len(A))]
    def unitv(idx, n, val=1):
        return [val if i == idx else 0 for i in range(n)]
    ieqs = [tuple([0] + unitv(idx, len(A[0]))) for idx in range(len(A[0]))]
    for idx in range(len(A[0])):
        ie = [tuple([0] + unitv(idx2, len(A[0]), val = -1 if idx==idx2 else 1)) for idx2 in range(len(A[0]))]

        P = Polyhedron(ieqs=ie, eqns=eq)
        if P.dim() != dim:
            new_A, new_w = contract(A, w, idx)
            return tighten(new_A, new_w, dim=dim)
    return A, w


def calculate_cells2(A):
    '''
    Calculates maximal dimensional cells that can be obtained from a quiver. Calculates the number of lattice points and the number of faces in each dimension for every cell.
    Parameters:
    A (integer matrix) : the signed adjecency matrix of the quiver

    Returns:
    list: list of tuples corresponding to each cell: (weight, number of lattice points, (number of faces in each dimension), Polytope)
    '''
    dim = len(A[0]) - len(A) + 1
    def val(idx):
        return sum(A[idx])
    from itertools import product
    def weight_list(idx):
        v = val(idx)
        if v == 2:
            return [1]
        elif v < 0:
            return [-w for w in list(range(1, -v))]
        else:
            raise TypeError()
    possible_weights = list(filter(lambda x: sum(x) == 0, product(*[weight_list(idx) for idx in range(len(A))])))
    result = [] #tuple (weight, num_lattice_points, num_faces, polynom)
    for w in possible_weights:
        P = pol(A, w)
        if P.dim() == dim:
            num_faces = tuple(len(P.faces(i)) for i in range(0, dim))
            result.append([w, len(P.integral_points()), num_faces, P])
    return result

def toric_ideal(P):
    '''Calculates the toric ideal of a polytope'''
    pts = P.integral_points()
    M = [[1] * len(pts)]
    for i in range(len(pts[0])):
        M.append([x[i] for x in pts])
    return ToricIdeal(M)

def ideal_equal(I1, I2):
    return all(g in I2 for g in I1.gens()) and all(g in I1 for g in I2.gens())

def similar_cell(c1, c2):
    return c1[1] == c2[1] and c1[2] == c2[2]

def group_cells(cells):
    result = []
    for c in cells:
        found = False
        for r in result:
            if similar_cell(c, r[0]):
                r.append(c)
                found = True
                break
        if not found:
            result.append([c])
    return result

def ideal_perm_equal(I1, I2):
    from itertools import permutations
    if I1.ring().ngens() != I2.ring().ngens():
        return False
    maps = (I1.ring().hom(p, I2.ring()) for p in permutations(I2.ring().gens()))
    return any(ideal_equal(m(I1), I2) for m in maps)
︡33d87ffc-9bbb-400c-b08f-7bf05b63a329︡{"done":true}
︠54e02ed8-2cd5-4572-a794-e5a77cd4b788s︠
# First we list the edges of the complete bipartite graph K3,3
edges = [[0, 3], [0, 4], [0, 5], [1, 3], [1, 4], [1, 5], [2, 3], [2, 4], [2, 5]]
# We calculate the signed adjecency matrix of the quiver obtained by putting a valency 2 sink on each edge
Q = edges_to_quiver(edges)
# First we calculate a complete but redundant list of maximal dimensional quiver cells that can be obtained from Q, the list below shows the corresponding weights, the number of lattice points and the number of 0..3 dimensional faces of eachy polytope
cells = calculate_cells2(Q)
for c in cells:
    print(c[:-1])
︡10293f09-2395-4bf6-a2eb-96017a9159f0︡{"stdout":"[(-1, -1, -1, -2, -2, -2, 1, 1, 1, 1, 1, 1, 1, 1, 1), 6, (6, 15, 18, 9)]\n[(-1, -1, -2, -1, -2, -2, 1, 1, 1, 1, 1, 1, 1, 1, 1), 5, (5, 10, 10, 5)]\n[(-1, -1, -2, -2, -1, -2, 1, 1, 1, 1, 1, 1, 1, 1, 1), 5, (5, 10, 10, 5)]\n[(-1, -1, -2, -2, -2, -1, 1, 1, 1, 1, 1, 1, 1, 1, 1), 5, (5, 10, 10, 5)]\n[(-1, -2, -1, -1, -2, -2, 1, 1, 1, 1, 1, 1, 1, 1, 1), 5, (5, 10, 10, 5)]\n[(-1, -2, -1, -2, -1, -2, 1, 1, 1, 1, 1, 1, 1, 1, 1), 5, (5, 10, 10, 5)]\n[(-1, -2, -1, -2, -2, -1, 1, 1, 1, 1, 1, 1, 1, 1, 1), 5, (5, 10, 10, 5)]\n[(-1, -2, -2, -1, -1, -2, 1, 1, 1, 1, 1, 1, 1, 1, 1), 5, (5, 10, 10, 5)]\n[(-1, -2, -2, -1, -2, -1, 1, 1, 1, 1, 1, 1, 1, 1, 1), 5, (5, 10, 10, 5)]\n[(-1, -2, -2, -2, -1, -1, 1, 1, 1, 1, 1, 1, 1, 1, 1), 5, (5, 10, 10, 5)]\n[(-2, -1, -1, -1, -2, -2, 1, 1, 1, 1, 1, 1, 1, 1, 1), 5, (5, 10, 10, 5)]\n[(-2, -1, -1, -2, -1, -2, 1, 1, 1, 1, 1, 1, 1, 1, 1), 5, (5, 10, 10, 5)]\n[(-2, -1, -1, -2, -2, -1, 1, 1, 1, 1, 1, 1, 1, 1, 1), 5, (5, 10, 10, 5)]\n[(-2, -1, -2, -1, -1, -2, 1, 1, 1, 1, 1, 1, 1, 1, 1), 5, (5, 10, 10, 5)]\n[(-2, -1, -2, -1, -2, -1, 1, 1, 1, 1, 1, 1, 1, 1, 1), 5, (5, 10, 10, 5)]\n[(-2, -1, -2, -2, -1, -1, 1, 1, 1, 1, 1, 1, 1, 1, 1), 5, (5, 10, 10, 5)]\n[(-2, -2, -1, -1, -1, -2, 1, 1, 1, 1, 1, 1, 1, 1, 1), 5, (5, 10, 10, 5)]\n[(-2, -2, -1, -1, -2, -1, 1, 1, 1, 1, 1, 1, 1, 1, 1), 5, (5, 10, 10, 5)]\n[(-2, -2, -1, -2, -1, -1, 1, 1, 1, 1, 1, 1, 1, 1, 1), 5, (5, 10, 10, 5)]\n[(-2, -2, -2, -1, -1, -1, 1, 1, 1, 1, 1, 1, 1, 1, 1), 6, (6, 15, 18, 9)]\n"}︡{"done":true}
︠8ebd113e-39e2-4747-ba69-248684eec42cs︠
# While there is no theoretical reason for this to be case, it turns out that the combinatorial data in the above table uniquely identifies each quiver cell (up to isomoprhism) in dimension 4 or less. This follows from the classification results in Propositions 5.2, 8.3. The short script below groups the cells
# by this data and verifies that in each group the toric ideals are isomorphic via a map that permutes the variables (this also needs not to work in a general case, but is sufficient for every quiver up to dim 4).
grouped_cells = group_cells(cells)
for group in grouped_cells:
    for g in group:
        g.append(toric_ideal(g[-1]))
for group in grouped_cells:
    print(all(ideal_perm_equal(group[0][4], g[4]) for g in group))
︡32d99e1a-31b1-4add-a428-c71efaa1b170︡{"stdout":"True\nTrue\n"}︡{"done":true}
︠d8f7a6c8-5b6f-4eb2-acc3-b80f0b913ba8s︠
for group in grouped_cells:
    print(group[0])
︡cc96a4c6-ce19-43f0-8c2e-ea35d8ff33ab︡{"stdout":"[(-1, -1, -1, -2, -2, -2, 1, 1, 1, 1, 1, 1, 1, 1, 1), 6, (6, 15, 18, 9), A 4-dimensional polyhedron in QQ^18 defined as the convex hull of 6 vertices, Ideal (z0*z2*z4 - z1*z3*z5) of Multivariate Polynomial Ring in z0, z1, z2, z3, z4, z5 over Rational Field]\n[(-1, -1, -2, -1, -2, -2, 1, 1, 1, 1, 1, 1, 1, 1, 1), 5, (5, 10, 10, 5), A 4-dimensional polyhedron in QQ^18 defined as the convex hull of 5 vertices, Ideal (0) of Multivariate Polynomial Ring in z0, z1, z2, z3, z4 over Rational Field]\n"}︡{"done":true}
︠ff979453-3975-4595-b802-308c16094910s︠
# In accordance with Proposition 5.3 we obtain that there are 2 quiver cells with this chassis, the Birkhoff polytope and the simplex. We take a weight representing each and calculate the tightened quiver, weight pairs
birkhoff_weight = grouped_cells[0][0][0]
simplex_weight = grouped_cells[1][0][0]
birkhoff_weight
simplex_weight
︡509bc479-ccc0-4394-afeb-efbd3b42e392︡{"stdout":"(-1, -1, -1, -2, -2, -2, 1, 1, 1, 1, 1, 1, 1, 1, 1)\n"}︡{"stdout":"(-1, -1, -2, -1, -2, -2, 1, 1, 1, 1, 1, 1, 1, 1, 1)\n"}︡{"done":true}
︠3dba9f03-b39b-4095-8eff-708043340da0︠
#the list of non contractible arrows
tight_list(Q, birkhoff_weight)
︡32410632-1816-47a6-bf9b-67d7b05cee2e︡{"stdout":"[0, 1, 2, 3, 4, 5, 6, 7, 8]"}︡{"stdout":"\n"}︡{"done":true}
︠29a431f5-acb2-44e8-8d30-5f3aa9c08061s︠
#We recover the "standard" way to obtain the Birkhoff polytope as a quiver cell, as a complete bipartite graph K3,3 with arrows pointing from one bipartition from the other and the weight is -1 on each source and 1 on each sink
tighten(Q, birkhoff_weight)
︡fe0838bc-6a89-4a01-9671-793c3aa64ae1︡{"stdout":"([[-1, -1, -1, 0, 0, 0, 0, 0, 0], [0, 0, 0, -1, -1, -1, 0, 0, 0], [0, 0, 0, 0, 0, 0, -1, -1, -1], [1, 0, 0, 1, 0, 0, 1, 0, 0], [0, 1, 0, 0, 1, 0, 0, 1, 0], [0, 0, 1, 0, 0, 1, 0, 0, 1]], [-1, -1, -1, 1, 1, 1])"}︡{"stdout":"\n"}︡{"done":true}
︠5ca88457-6160-4188-b844-ed69f87e3a7as︠
tight_list(Q, simplex_weight)
︡bd05cee0-7dd9-4681-a403-63a8432db931︡{"stdout":"[1, 2, 4, 5, 15]"}︡{"stdout":"\n"}︡{"done":true}
︠d0ee309c-a60f-40f1-ac4f-17db01b61930︠
#For the simplex we obtain the kronecker quiver after tightening
tighten(Q, simplex_weight)
︡7980d17d-6b74-4b13-9ead-b27909ad0dbb︡{"stdout":"([[-1, -1, -1, -1, -1], [1, 1, 1, 1, 1]], [-1, 1])"}︡{"stdout":"\n"}︡{"done":true}









