from sage.all import *
import random
import itertools
import json
import re
'''
This file contains scripts and an example calculation of Mátyás Domokos and Dániel Joó created for the paper "LOW DIMENSIONAL FLOW POLYTOPES AND THEIR TORIC IDEALS" (https://arxiv.org/abs/2105.04004), which contains explanations for the terminology used in the docstrings here.  
The file is intended to be used with the offline version of Sage and Python 7.2 or newer. As the purpose of making this file was to study some low dimensional examples, the scipts here are likely badly optimized. 
Obtaining full list of pulling triangulations of larger (vertices > 10) polytopes is lengthy and the resulting datafiles will be very large.
For questions contact the authors at: jood.mollia@gmail.com
'''

def toric_ideal(P):
    pts = P.integral_points()
    M = [[1] * len(pts)]
    for i in range(len(pts[0])):
        M.append([x[i] for x in pts])
    return ToricIdeal(M)

def cleanup_gens(ideal):
    ring = ideal.ring()
    gens = ideal.gens()
    degrees = [p.degree() for p in gens]
    max_deg = max(degrees)
    sorted = [[] for i in range(max_deg - 1)]
    for i in range(len(sorted)):
        if i == 0:
            sorted[0] = list(filter(lambda p: p.degree() == 2, gens))
        else:
            current_degree = list(filter(lambda p: p.degree() == i + 2, gens))
            lower_gens = sum(sorted[:i], [])
            lower_ideal = ring.ideal(lower_gens)
            for p in current_degree:
                if p not in lower_ideal:
                    sorted[i].append(p)
    return sorted

def max_non_empty(lists_of_things):
    l = []
    for idx in range(len(lists_of_things)):
        if len(lists_of_things[idx]) > 0:
            l.append(idx)
    return max(l)   

def is_compressed(P):
    return all(len(set(ineq.eval(v) for v in P.vertices())) <= 2 for ineq in P.Hrepresentation())

class PullingTriangulation:

    def __init__(self, triangulation, pull_order=[]):
        self.pull_order = pull_order
        self.triangulation = triangulation
    
    def extend(self, v):
        self.pull_order.insert(0, v)
        self.triangulation = [tri + [v] for tri in self.triangulation]

class PulledFacePoset:

    def __init__(self, face_poset, pulled_vertices=[]):
        self.face_poset = face_poset
        self.pulled_vertices = pulled_vertices
        self.top = face_poset.top()
        self.vertices = self.top.vertices()




def index_triangulation(triangulation, vertex_list):

    def replace_with_index(t):
        return [vertex_list.index(v) for v in t]
    return [replace_with_index(t) for t in triangulation]

def index_list(subset, fullset):
    return [fullset.index(x) for x in subset]

def pull_polyhedron(vertex, polyhedron):
    result = []
    for facet in polyhedron.facets():
        if vertex not in facet.vertices():
            result.append(list(facet.vertices()) + [vertex])
    return result

def pulling_refinement(pulling_subdivision, vertex):
    # pulling_subdivision[0]: pull order, [1]: [PulledFaceLattice, done]
    result = [[], []]
    result[0] = pulling_subdivision[0] + [vertex]
    
    for s in pulling_subdivision[1]:

        done = s[1]
        if done:
            result[1].append(s)
            continue
        if vertex not in s[0].vertices:
            result[1].append(s)
            continue

        for f in s[0].face_poset.lower_covers(s[0].top):
            if vertex in f.vertices():
                continue
            if f.dim() == len(f.vertices()) - 1:
                result[1].append([list(f.vertices()) + s[0].pulled_vertices + [vertex], True])
            else:
                new_pulled_poset = PulledFacePoset(s[0].face_poset.subposet(s[0].face_poset.order_ideal([f])), s[0].pulled_vertices + [vertex])
                result[1].append([new_pulled_poset, False])

    return result


def all_pulling_triangulations(P):

    if P.is_simplex():
        return [PullingTriangulation([list(P.vertices())])]    
    all_refinements = all_pulling_refinements([[], [[PulledFacePoset(P.face_lattice(), []), False]]])
    all_triang = [PullingTriangulation([x[0] for x in r[1]], r[0]) for r in all_refinements]
    return all_triang
def all_pulling_refinements(pulling_subdivision):

    if all(x[1] for x in pulling_subdivision[1]):
        return [pulling_subdivision]
    pullable_vertices = []
    for s in pulling_subdivision[1]:
        if not s[1]:
            for v in s[0].vertices:
                if v not in pulling_subdivision[0] and v not in pullable_vertices:
                    pullable_vertices.append(v)
    return sum((all_pulling_refinements(pulling_refinement(pulling_subdivision, v)) for v in pullable_vertices), [])

def get_polytope_data(P):
    '''Given a polytope P calculates every pulling triangulation of P. The dictionary returned contains a number of additional data about the polytope and the triangulations. 
    Parameters:
    P(Polytope): a lattice polytope

    Returns:
    dictionary: 'vertices': the vertices of P
                'compressed': is P a compressed polytope
                'generated_in_degree': the lowest degree that generates the toric ideal of P
                'triangulations': a full list of pulling triangulations. Simplices in the triangulation are given by the indices of the vertices
                'best_triangulation': the triangulation with the smallest maximal non-face
                'worst_triangulation': the triangulation with the largest maximal non-face
    '''
    data = {}
    vertices = P.vertices()
    data['vertices'] = [[int(x) for x in list(v)] for v in vertices]
    data['compressed'] = is_compressed(P)
    data['generated_in_degree'] = max_non_empty(cleanup_gens(toric_ideal(P))) + 2
    data['triangulations'] = []
    all_pulling = all_pulling_triangulations(P)
    for pull in all_pulling:
        S = SimplicialComplex(index_triangulation(pull.triangulation, vertices))
        D = S.alexander_dual()
        def invert(x):
            return [a for a in range(len(vertices)) if a not in x]
        non_faces = [invert(x) for x in D.facets()]  
        data['triangulations'].append([index_list(pull.pull_order, vertices), index_triangulation(pull.triangulation, vertices), non_faces, max(len(f) for f in non_faces)])
    data['best_triangulation'] = min(x[3] for x in data['triangulations'])
    data['worst_triangulation'] = max(x[3] for x in data['triangulations'])
    return data

#For an example calculation we use the quiver polytope VI_C as defined in Section 8 of the paper
VI_C = [[0, 0, 0, 0], [0, 0, 1, 0], [1, 0, 1, 0], [1, 0, 0, 1], [0, 1, 0, 0], [0, 1, 1, 0], [0, 1, 0, 1]]
c = [1, 2, 3, 4]

P = Polyhedron(VI_C)
data = get_polytope_data(P)

with open(".\\polytope_data_filtered\\VI_C_data.json", "w") as f:
    json.dump(data, f)

