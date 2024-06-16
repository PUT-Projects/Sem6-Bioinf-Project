from dnaParser import Cell
import dnaParser
import offset
import time

def sort_by_pos(data: list[Cell]) -> list[Cell]:
    # sort by getPosH() and getPosL()
    return sorted(data, key=lambda x: (x.getPosH(), x.getPosL()))

class PathFinder:
    def __init__(self, graph, diff_matrix, start_id, length):
        self.__graph = graph
        self.__diff_matrix = diff_matrix
        self.__length = length
        self.__start_id = start_id
        self.best_path = []
        self.max_vertices = 0

    def run(self):
        # Start the search from the specified start_id
        visited = [False] * len(self.__graph)
        self.__backtrack(self.__start_id, [], visited, 0)
        return self.best_path

    def __backtrack(self, current_vertex, current_path, visited, current_cost):
        # Include current vertex in path
        visited[current_vertex] = True
        current_path.append(current_vertex)

        # Check if the current cost equals the desired length
        if current_cost == self.__length:
            # Check if the current path is longer (in terms of vertices) than the previously found paths
            if len(current_path) > self.max_vertices:
                self.max_vertices = len(current_path)
                self.best_path = current_path.copy()
            # Since we found a perfect cost match, we can return here for optimization
            visited[current_vertex] = False
            current_path.pop()
            return

        # Explore adjacent vertices
        for neighbor in self.__graph[current_vertex]:
            if not visited[neighbor]:
                edge_cost = self.__diff_matrix[neighbor][current_vertex]
                if current_cost + edge_cost <= self.__length:  # Only consider if within total cost
                    self.__backtrack(neighbor, current_path, visited, current_cost + edge_cost)

        # Backtrack
        visited[current_vertex] = False
        current_path.pop()

# map of oligonucleotides and connections between nodes
class Graph:
    def __init__(self, data: dnaParser.Probe, start: str, length: int) -> None:
        data.getCells().append(Cell(0, 0, start))

        self.__data = sort_by_pos(data.getCells())
        self.__start = start
        self.__graph = []
        self.__length = int(length)
        self.__start_id = 0
        self.__diff_matrix = self.calculateDiffMatrix(self.__data)
        self.__pattern_length = len(self.__data[0].getSequence())
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
            for j in range(len(self.__data)):
                if i == j: continue

                if (self.__data[i].getPosH() + 1 >= self.__data[j].getPosL()) and (self.__data[i].getPosL() + self.__diff_matrix[j][i] <= self.__data[j].getPosH() + 1):
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
            cost += self.costAt(path[i+1], path[i])

        return cost

    def sts(self):
        #Selective traveling salesman
        #Maximize the number of nodes visited
        #cost of the path cannot exceed the self.__length
        # if len(path) == self.__length, return path
        #Start at the start_id node
        # use Bellman-Held-Karp algorithm

        finder = PathFinder(self.__graph, self.__diff_matrix, self.__start_id, self.__length - self.__pattern_length)
        return finder.run()


        # for minLen in range(len(self.__graph), 0, -1):
        #   while len(path) < minLen:
        #     # TOD
        #   if cost <= self.__length:
        #     return path


    def getSequence(self, path: list[int]) -> str:
            #based on cell Sequence and diff matrix
        string = ""
        cost = 0
        for i in range(len(path) - 1):
            #print(" " * cost, self.__data[path[i]].getSequence())
            count = self.costAt(path[i+1], path[i])
            cost += count
            string += self.__data[path[i]].getSequence()[:count]

        #print(" " * cost, self.__data[path[-1]].getSequence())

        string += self.__data[path[-1]].getSequence()

        return string


def solve(data: list[Cell]):
    sorted_data = sort_by_pos(data)
    result = []
    #print(sorted_data)
    print('\n'.join(str(s) for s in sorted_data))


def main() -> None:
    #dna = dnaParser.DNA().loadFile("input_easy2.xml")
    #dna = dnaParser.DNA().loadXML(dnaParser.getInputFromWeb(18, 5, sqpe=4, pose=4))
    dna = dnaParser.DNA().loadXML(dnaParser.getInputFromFile('input.xml'))
    #print(dna.getStart())
    #solve(dna.getProbes()[0].getCells())

    print("Start")
    start = time.time()

    graph = Graph(dna.getProbes()[0], dna.getStart(), dna.getLength())
    graph.sortGraphAdjacencies()
    graph.printGraph()
    sts = graph.sts()
    print(sts)

    seq = graph.getSequence(sts)
    print(seq, len(seq))

    end = time.time()
    print("End")
    print("Time:", end-start)




if __name__ == "__main__":
    main()