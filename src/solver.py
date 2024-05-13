from dnaParser import Cell
import dnaParser

def sort_by_pos(data: list[Cell]) -> list[Cell]:
    # sort by getPosH() and getPosL()
    return sorted(data, key=lambda x: (x.getPosH(), x.getPosL()))

# map of oligonucleotides and connections between nodes
class Graph:
    def __init__(self, data: dnaParser.Probe):
        self.__probe = data
        self.graph = {}
        self.createGraph()

    def createGraph(self):
        sort_by_pos(self.__probe.getCells())

        for cell in self.__probe.getCells():
            #PosH_i + 1 >= PosL_{j}
            for cell2 in self.__probe.getCells():
                if cell.getPosH() + 1 >= cell2.getPosL():
                    self.addEdge(cell, cell2)
                else:
                    break

    def addEdge(self, cell1: Cell, cell2: Cell):
        if cell1 not in self.graph:
            self.graph[cell1] = []
        self.graph[cell1].append((cell2, self.getDiff(cell1, cell2)))

    def find_offset_optimized(self, str1, str2):
        # Assign integer values to characters
        char_to_int = {'A': 1, 'C': 2, 'T': 3, 'G': 4}
        base = 5
        max_len = len(str2)

        # Compute hash for str2
        target_hash = 0
        for char in str2:
            target_hash = target_hash * base + char_to_int[char]

        current_hash = 0
        highest_base_pow = base ** (max_len - 1)

        # Compute the hash of substrings in str1 and compare
        for i in range(len(str1)):
            current_hash = current_hash * base + char_to_int[str1[i]]

            if i >= max_len:
                current_hash -= char_to_int[str1[i - max_len]] * highest_base_pow

            if i >= max_len - 1:
                if current_hash == target_hash:
                    if str1[i - max_len + 1:i + 1] == str2:
                        return i - max_len + 1

        return len(str1)  # If no match found, return length of str1 as specified

    def getDiff(self, cell1: Cell, cell2: Cell):
        #get diff in letters between two cells
        #find substrings in cell1 and cell2
        #like: AACTGC and ACTCCG -> 1
        #like: AACTGC and TGCCTA -> 3
        #like: AACTGC and GAATTC -> 6
        #like: AACTGC and AACTGC -> 0

        if len(cell1.getSequence()) != len(cell2.getSequence()):
            raise ValueError("Sequences have different lengths")

        diff = 0

        #AACTGC and ACTCCG
        #A|ACTGC ACTCC|G -> 1

        #like: AACTGC and TGCCTA
        # AAC|TGC TGC|CTA -> 3

        #like: AACTGC and GAATTC -> 6
        #like: AACTGC and AACTGC -> 0

        #worst case (no substrings)
        return len(cell1.getSequence())

def solve(data: list[Cell]):
    sorted_data = sort_by_pos(data)
    result = []
    #print(sorted_data)
    print('\n'.join(str(s) for s in sorted_data))

def main() -> None:


    dna = dnaParser.DNA().loadFile("input.xml")
    print(dna.getStart())
    solve(dna.getProbes()[0].getCells())


if __name__ == "__main__":
    main()