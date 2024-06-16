import heapq
import dnaParser
import offset
import dnaParser
from dnaParser import Cell
import time



def sort_by_pos(data: list[Cell]) -> list[Cell]:
    # sort by getPosH() and getPosL()
    return sorted(data, key=lambda x: (x.getPosH(), x.getPosL()))


# map of oligonucleotides and connections between nodes
class Graph:
    def __init__(self, data: dnaParser.Probe, start: str, length: int) -> None:
        self.__data = list(sort_by_pos(data.getCells()))
        self.__start = start
        self.__graph = []
        self.__length = int(length)
        self.__start_id = -1
        self.__diff_matrix = self.calculateDiffMatrix(self.__data)
        self.createGraph()

    def getGraph(self) -> dict:
        return self.__graph

    def printGraph(self) -> None:
        # print all adjacencies
        for i in range(len(self.__graph)):
            sequence = self.__data[i].getSequence()
            posl = self.__data[i].getPosL()
            posh = self.__data[i].getPosH()
            #list of adjacencies and their diff matrix values
            adjacencies = self.__graph[i]
            diffs = [self.__diff_matrix[j][i] for j in adjacencies]

            xd = [f"{adjacencies[x]}({diffs[x]})" for x in range(len(adjacencies))]

            print(f"{i} ({sequence}) [ {posl},{posh} ]: [ {', '.join(xd)} ]")
            #print (f"{i} ({self.__data[i].getSequence()}) [{self.__data[i].getPosL()},{self.__data[i].getPosH()}]: {self.__graph[i]}")
            #print(f"{i} ({self.__data[i].getSequence()}) [{self.__data[i].getPosL()},{self.__data[i].getPosH()}]: {self.__graph[i][j] self.__diff_matrix[i][j]  for j in self.__graph[i]}")

    def createGraph(self) -> None:
        #add edge 0 as starting point
        self.__graph = [[] for _ in range(len(self.__data))]

        for i in range(len(self.__data)):
            if 1 >= self.__data[i].getPosL() and self.__start == self.__data[i].getSequence():
                self.__start_id = i
                self.__data[i].setPosH(0)
                self.__data[i].setPosL(0)
                break

        for i in range(len(self.__data)):
            for j in range(len(self.__data)):
                if i == j: continue

                # PosL[i] = 1, PosH[i] = 3
                # PosL[j] = 6, PosH[j] = 8
                # 3+1 >= 6, 1 <= 8
                # false   , true

                # posL[i] = 11, posH[i] = 13
                # posL[j] = 6, posH[j] = 8
                # 13+1 >= 6, 11 <= 8
                # true     , false

                # posL[i] = 11, posH[i] = 13
                # posL[j] = 14, posH[j] = 16
                # 13+1 >= 14, 11 <= 16
                # true     , true

                if (self.__data[i].getPosH() + 1 >= self.__data[j].getPosL()) and (self.__data[i].getPosL() <= self.__data[j].getPosH()):
                    self.__graph[i].append(j)

        return None


    def diff(self, cell1: Cell, cell2: Cell) -> int:
        if len(cell1.getSequence()) != len(cell2.getSequence()):
            raise ValueError("Sequences have different lengths")

        diff = offset.find_offset_naive(cell1.getSequence(), cell2.getSequence())

        if diff < 0:
            raise ValueError("Negative diff")

        return diff if diff != 0 else 2**30

    def costAt(self, i: int, j: int) -> int:
        return self.__diff_matrix[i][j]

    def calculateDiffMatrix(self, data: list[Cell]):
        n = len(data)

        matrix = [[self.diff(data[i], data[j]) for i in range(n)] for j in range(n)]

        return matrix

    def sortGraphAdjacencies(self):
        for i in range(len(self.__graph)):
            self.__graph[i] = sorted(self.__graph[i], key=lambda x: self.__diff_matrix[x][i])

    def reccursive(self, current: int, visited: list[bool], path: list[int], cost: int, best_path) -> None:
        if len(path) == self.__length:
            return path

        for i in self.__graph[current]:
            if visited[i] == False:
                visited[i] = True
                path.append(i)
                cost += self.costAt(current, i)
                self.reccursive(i, visited, path, cost)
                if len(path) == self.__length:
                    best_path = path
                    return

                visited[i] = False
                path.pop()



        return path

    def total_cost(self, path: list[int]) -> int:
        cost = 0
        for i in range(len(path) - 1):
            cost += self.costAt(path[i], path[i+1])

        return cost

    def getSequence(self, path: list[int]) -> None:
            #based on cell Sequence and diff matrix
        print()
        cost = 0
        for i in range(len(path) - 1):
            print(" " * cost, self.__data[path[i]].getSequence())
            count = self.costAt(path[i+1], path[i])
            cost += count
            #string += self.__data[path[i]].getSequence()[:count]


        #string += self.__data[path[-1]].getSequence()

        return None



    def modified_dijkstra(self):
        # Initialize distances to infinity and the start node distance to 0

        distances = {node: float('inf') for node in range(len(self.__data))}
        distances[self.__start_id] = len(self.__data[self.__start_id].getSequence())

        max_nodes = len(self.__graph)

        # Priority queue to store nodes with their tentative distances
        priority_queue = [(0, self.__start_id)]

        # Visited set to keep track of visited nodes
        visited = set()

        # Keep track of total cost and visited nodes
        total_cost = 0
        visited_nodes = 0

        while priority_queue:
            # Get the node with the smallest tentative distance
            current_cost, current_node = heapq.heappop(priority_queue)

            # Check if the node has already been visited
            if current_node in visited:
                continue

            # Update total cost and visited nodes
            total_cost += current_cost
            visited_nodes += 1

            # Mark the current node as visited
            visited.add(current_node)

            # Check if the maximum traversal cost has been reached
            if total_cost > self.__length:
                break

            # Check if the maximum number of visited nodes has been reached
            if visited_nodes > max_nodes:
                break

            # Explore neighbors of the current node
            for neighbor in self.__graph[current_node]:
                if neighbor not in visited:
                    new_distance = distances[current_node] + self.costAt(current_node, neighbor)
                    if new_distance < distances[neighbor]:
                        distances[neighbor] = new_distance
                        heapq.heappush(priority_queue, (new_distance, neighbor))

        # Return the visited nodes and total cost
        return visited, total_cost



def main() -> None:

    dna = dnaParser.DNA().loadFile("input_easy.xml")
    #print(dna.getStart())
    #solve(dna.getProbes()[0].getCells())

    graph = Graph(dna.getProbes()[0], dna.getStart(), dna.getLength())
    graph.sortGraphAdjacencies()
    graph.printGraph()
    path, cost = graph.modified_dijkstra()
    print("Path: ", path, ", Cost: ", cost)

    graph.getSequence(path)


    end = time.time()
    print("End")
    print("Time:", end-start)   


if __name__ == '__main__':

    main()