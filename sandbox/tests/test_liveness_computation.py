from __future__ import absolute_import, division, print_function

import sys
sys.path.insert(0, '../')

from compiler import *
from constructs import *

from schedule import compute_liveness

class Node:
    def __init__(self, _name, _parents=None):
        self._name = _name
        self._parents = _parents

    @property
    def name(self):
        return self._name
    @property
    def parents(self):
        return self._parents

def get_children_map(nodes):
    children_map = {}
    for node in nodes:
        parents = node.parents
        if parents == None:
            continue
        for parent in parents:
            assert isinstance(parent, Node)
            if parent not in children_map:
                children_map[parent] = [node]
            else:
                children_map[parent] += [node]
    return children_map

def compare(map1, map2):
    is_equal = True

    keys1 = set([key for key in map1])
    keys2 = set([key for key in map2])

    is_equal = not keys1.difference(keys2)

    for key in map1:
        val1 = set(map1[key])
        val2 = set(map2[key])
        is_equal = is_equal and not val1.difference(val2)

    return is_equal

def test_graph1():
    # Test Graph 1 : A chain
    # a -> b -> c -> d -> e -> f
    a = Node('a')
    b = Node('b', [a])
    c = Node('c', [b])
    d = Node('d', [c])
    e = Node('e', [d])
    f = Node('f', [e])

    nodes = [a, b, c, d, e, f]

    children_map = get_children_map(nodes)

    schedule = { \
    a:0, b:1, c:2, d:3, e:4, f:5 \
    }

    liveness_map = compute_liveness(children_map, schedule)

    verify_map = { \
    1:[a], 2:[b], 3:[c], 4:[d], 5:[e] \
    }

    assert compare(liveness_map, verify_map)

def test_graph2():
    # Test Graph 2 : Chain-like with multiple children
    # a -> b -> c -> d -> e -> f
    # a -> c -> e
    # b -> d -> f
    a = Node('a')
    b = Node('b', [a])
    c = Node('c', [a, b])
    d = Node('d', [b, c])
    e = Node('e', [c, d])
    f = Node('f', [d, e])

    nodes = [a, b, c, d, e, f]

    children_map = get_children_map(nodes)

    schedule = { \
    a:0, b:1, c:2, d:3, e:4, f:5 \
    }

    liveness_map = compute_liveness(children_map, schedule)

    verify_map = { \
    2:[a], 3:[b], 4:[c], 5:[d, e] \
    }

    assert compare(liveness_map, verify_map)

def test_graph3():
    # Test Graph 3 : Tree
    # a -> b, c
    # b -> d, e
    # c -> f, g
    a = Node('a')
    b = Node('b', [a])
    c = Node('c', [a])
    d = Node('d', [b])
    e = Node('e', [b])
    f = Node('f', [c])
    g = Node('g', [c])

    nodes = [a, b, c, d, e, f, g]

    children_map = get_children_map(nodes)

    schedule = { \
    a:0, b:1, c:2, d:3, e:4, f:5, g:6 \
    }

    liveness_map = compute_liveness(children_map, schedule)

    verify_map = { \
    2:[a], 4:[b], 6:[c] \
    }

    assert compare(liveness_map, verify_map)

def test_graph4():
    # Test Graph 4 : Harris
    # img -> Ix, Iy
    # Ix -> Ixx, Ixy
    # Iy -> Ixy, Iyy
    # Ixx -> Sxx
    # Ixy -> Sxy
    # Iyy -> Syy
    # Sxx -> det, trace
    # Syy -> det, trace
    # Sxy -> det
    # det -> harris
    # trace -> harris
    img = Node('img')
    Ix = Node('Ix', [img])
    Iy = Node('Iy', [img])
    Ixx = Node('Ixx', [Ix])
    Ixy = Node('Ixy', [Ix, Iy])
    Iyy = Node('Iyy', [Iy])
    Sxx = Node('Sxx', [Ixx])
    Sxy = Node('Sxy', [Ixy])
    Syy = Node('Syy', [Iyy])
    det = Node('det', [Sxx, Sxy, Syy])
    trace = Node('trace', [Sxx, Syy])
    harris = Node('harris', [det, trace])

    nodes = [img, Ix, Iy, Ixx, Ixy, Iyy, Sxx, Sxy, Syy, det, trace, harris]

    children_map = get_children_map(nodes)

    def test_sched1():
        # schedule - 1
        schedule1 = { \
        img:0, Ix:1, Iy:2, Ixx:3, Ixy:4, Iyy:5, \
        Sxx:6, Sxy:7, Syy:8, det:9, trace:10, harris:11 \
        }

        liveness_map = compute_liveness(children_map, schedule1)

        verify_map = { \
        2:[img], 4:[Ix], 5:[Iy], 6:[Ixx], 7:[Ixy], 8:[Iyy], \
        9:[Sxy], 10:[Sxx, Syy], 11:[det, trace] \
        }

        assert compare(liveness_map, verify_map)

    def test_sched2():
        # schedule - 2
        schedule2 = { \
        img:0, Iy:1, Iyy:2, Syy:3, Ix:4, Ixy:5, \
        Sxy:6, Ixx:7, Sxx:8, trace:9, det:10, harris:11 \
        }

        liveness_map = compute_liveness(children_map, schedule2)

        verify_map = { \
        3:[Iyy], 4:[img], 5:[Iy], 6:[Ixy], 7:[Ix], 8:[Ixx], \
        10:[Sxx, Sxy, Syy], 11:[det, trace] \
        }

        assert compare(liveness_map, verify_map)

    test_sched1()
    test_sched2()

    return
