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
        perfect = self.perfectMatch(odd_verts,initial_matrix.copy(), min_tree)
        multigraph = self.multigraph(min_tree,perfect)
        print(multigraph)

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
                odds.append(i)
        return odds

    def minTree(self, matrix):
        min_matrix = min_tree(matrix)
        min_matrix = min_matrix.toarray().astype(float)
        return min_matrix

    def hierholzer(self, graph):
        # Initialize variables
        start_vertex = 0
        circuit = [start_vertex]
        edges_visited = {}
        # Loop through all edges that connect to the starting vertex
        for v in graph[start_vertex]:
            # If an edge exists and it hasn't been visited
            if graph[start_vertex][v] != np.inf and (start_vertex, v) not in edges_visited:
                # Mark as visited
                edges_visited.append((start_vertex, v))
                edges_visited.append((v, start_vertex))
                # Add it to the circuit
                circuit.append(v)
                # Initialize current path to be updated from following the edge to the next vertex
                curr_path = []
                self.search_new_vertex(graph, v, curr_path, edges_visited, start_vertex)
                # add the new path to the current circuit
                for i in range(len(curr_path)):
                    circuit.append(curr_path[i])

        return circuit

    def search_new_vertex(self, graph, u, curr_path, edges_visited, starting_vertex):
        # Loop through all edges that connect to the current vertex (u)
        for v in graph[u]:
            # If an edge exists and it hasn't been visited
            if graph[u][v] != np.inf and (u, v) not in edges_visited:
                # Mark as visited
                edges_visited.append((u, v))
                edges_visited.append((v, u))
                # Add it to the current path
                curr_path.append(v)
                # If we have completed the circuit, return; else, keep searching until the circuit is completed
                if v == starting_vertex:
                    return
                else:
                    self.search_new_vertex(graph, v, curr_path, edges_visited, starting_vertex)

    def perfectMatch(self, vertices, matrix, minMatrix):
        newmatrix = np.zeros(matrix.shape)
        numvertices = len(vertices)
        for i in range(matrix.shape[0]):
            if i not in vertices:
                matrix[i] = math.inf
                for j in range(matrix.shape[1]):
                    matrix[j][i] = math.inf
        while len(vertices) != 0:
            if len(vertices) == 1:
                print("this should never happen")
            else:
                pos = np.argmin(matrix)
                cols = matrix.shape[0]
                y = np.mod(pos, cols)
                x = pos // matrix.shape[0]
                if x in vertices and y in vertices:
                    if minMatrix[x][y] != matrix[x][y]:
                        vertices.remove(x)
                        vertices.remove(y)
                        newmatrix[x][y] = matrix[x][y]
                    elif minMatrix[y][x] != matrix[y][x]:
                        vertices.remove(x)
                        vertices.remove(y)
                        newmatrix[y][x] = matrix[y][x]
                    matrix[x][y] = math.inf
                    matrix[y][x] = math.inf
                    if self.checkPerfect(newmatrix, numvertices):
                        return newmatrix
                else:
                    matrix[x][y] = math.inf
                    matrix[y][x] = math.inf
                    continue
        return newmatrix

    def checkPerfect(self, matrix, numvertices):
        min = matrix.max(1)
        if np.count_nonzero(min) == numvertices // 2:
            return True
        else:
            return False

    def multigraph(self, matrix, perfectMatrix):
        return matrix + perfectMatrix
