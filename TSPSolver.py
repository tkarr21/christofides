#!/usr/bin/python3

from which_pyqt import PYQT_VER

if PYQT_VER == 'PYQT5':
    from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
    from PyQt4.QtCore import QLineF, QPointF
else:
    raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

import time
import numpy as np
from TSPClasses import *
from scipy.sparse.csgraph import minimum_spanning_tree as min_tree
import heapq
import itertools


class TSPSolver:
    def __init__(self, gui_view):
        self._scenario = None

    def setupWithScenario(self, scenario):
        self._scenario = scenario

    def defaultRandomTour(self, time_allowance=60.0):
        results = {}
        cities = self._scenario.getCities()
        ncities = len(cities)
        foundTour = False
        count = 0
        bssf = None
        start_time = time.time()
        while not foundTour and time.time() - start_time < time_allowance:
            # create a random permutation
            perm = np.random.permutation(ncities)
            route = []
            # Now build the route using the random permutation
            for i in range(ncities):
                route.append(cities[perm[i]])
            bssf = TSPSolution(route)
            count += 1
            if bssf.cost < np.inf:
                # Found a valid route
                foundTour = True
        end_time = time.time()
        results['cost'] = bssf.cost if foundTour else math.inf
        results['time'] = end_time - start_time
        results['count'] = count
        results['soln'] = bssf
        results['max'] = None
        results['total'] = None
        results['pruned'] = None
        return results

    def greedy(self, time_allowance=60.0):
        # Greedy algorithm to find a first tour fast that is near the optimal
        results = {}
        cities = self._scenario.getCities()
        ncities = len(cities)
        foundTour = False
        bssf = None
        start_time = time.time()
        start_node = 0
        # Loop until a tour is found or we time out
        while not foundTour and time.time() - start_time < time_allowance and start_node < ncities:
            badTour = False
            current_matrix = self.generateInitialMatrix()
            route = []
            # Start at the current start node, this will loop through every possible node as a start until a valid tour
            # is found or we time out
            current_node = start_node
            route.append(cities[current_node])
            current_matrix = blockCol(current_matrix, current_node)
            # We loop through enough times to create the length of a tour
            for i in range(ncities - 1):
                # From our current node grab the index for the smallest cost
                current_node = self.findMinIndex(current_matrix, current_node)
                # if our current_node is infinite then that means the lowest cost was infinite so this won't be a valid tour
                if current_node == np.inf:
                    badTour = True
                    break
                # append the node to the route and update the matrix so we don't revisit
                route.append(cities[current_node])
                current_matrix = blockCol(current_matrix, current_node)
            # create a TSPSolution based on our tour, if we had a bad tour or if the cost is infinite then it is not valid
            # so throw it out and try the next node as our starting position
            bssf = TSPSolution(route)
            if badTour:
                bssf.cost = np.inf
            if bssf.cost < np.inf:
                foundTour = True
                self._global_bssf = bssf
            start_node += 1
        end_time = time.time()
        results['cost'] = bssf.cost if foundTour else np.inf
        results['time'] = end_time - start_time
        results['count'] = None
        results['soln'] = bssf
        results['max'] = None
        results['total'] = None
        results['pruned'] = None
        return results

    def branchAndBound(self, time_allowance=60.0):
        pass

    def fancy(self, time_allowance=60.0):
        results = {}
        initial_matrix = self.generateInitialMatrix()
        min_tree = self.minTree(initial_matrix)
        odd_verts = self.getOddVerts(min_tree)
        print (min_tree)
        print (odd_verts)

    def generateInitialMatrix(self):
        i = 0
        j = 0
        cities = self._scenario.getCities()
        ncities = len(cities)
        matrix = np.empty([ncities, ncities])
        for i in range(ncities):
            for j in range(ncities):
                matrix[i, j] = cities[i].costTo(cities[j])
        return matrix

    def getOddVerts(self, matrix):
        odds = []
        for i in range(matrix.shape[0]):
            size = 0
            for j in range(matrix.shape[0]):
                if matrix[i, j] > 0:
                    size += 1
            for k in range(matrix.shape[0]):
                if matrix[k, i] > 0:
                    size += 1
            if (size % 2) != 0:
                print(size)
                odds.append(i)
        return odds



    def minTree(self, matrix):
        min_matrix = min_tree(matrix)
        min_matrix = min_matrix.toarray().astype(float)
        return min_matrix
