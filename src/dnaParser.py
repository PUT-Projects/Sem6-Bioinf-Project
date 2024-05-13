## File for parsing input xml file

## <dna key="22408638" length="500" start="GAGTATTTCT">
# <probe pattern="NNNNNNNNNN">
# <cell posL="250" posH="350">AAAAAGTGTG</cell>
# <cell posL="391" posH="491">AAAAATTATC</cell>
# <cell posL="370" posH="470">AAAAATTGTT</cell>

# store it in a list of objects
import xml.etree.ElementTree as ET

import requests


class Cell:
    def __init__(self, posL, posH, sequence):
        self.__posL = posL
        self.__posH = posH
        self.__sequence = sequence

    def __str__(self) -> str:
        return f"({self.__posL}, {self.__posH}, {self.__sequence})"

    def __repr__(self) -> str:
        return f"Cell({self.__posL}, {self.__posH}, {self.__sequence})"

# Getters
    def getPosL(self) -> int:
        return self.__posL

    def getPosH(self) -> int:
        return self.__posH

    def getSequence(self) -> str:
        return self.__sequence


# The `Probe` class has attributes for pattern and cells, with methods to retrieve the probe, pattern,
# and cells.
class Probe:
    def __init__(self, pattern, cells):
        self.__pattern = pattern
        self.__cells = cells

# Getters
    def getPattern(self) -> str:
        return self.__pattern

    def getCells(self) -> list[Cell]:
        return self.__cells

class DNA:
    def __init__(self, key=None, length=None, start=None, probes=None):
        self.__key = key
        self.__length = length
        self.__start = start
        self.__probes = probes


    def loadXML(self, xmlString):
        """
        The `loadXML` function parses an XML string, extracts specific attributes and elements, and creates
        Probe objects with Cell objects based on the XML structure.

        :param xmlString: The `loadXML` method you provided is used to parse an XML string and extract
        information to initialize the attributes of an object. The XML string should contain data in a
        specific format that the method expects
        :return: The `loadXML` method is returning the instance of the class itself after parsing the XML
        string and populating its attributes with the data extracted from the XML.
        """
        #parse xml string
        root = ET.fromstring(xmlString)
        self.__key = root.attrib['key']
        self.__length = root.attrib['length']
        self.__start = root.attrib['start']
        self.__probes = []

        for probe in root.findall('probe'):
            pattern = probe.attrib['pattern']
            cells = []
            for cell in probe.findall('cell'):
                posL = cell.attrib['posL']
                posH = cell.attrib['posH']
                sequence = cell.text
                cells.append(Cell(posL, posH, sequence))
            self.__probes.append(Probe(pattern, cells))
        return self


# Getters
    def getKey(self) -> str:
        return self.__key

    def getLength(self) -> int:
        return self.__length

    def getStart(self) -> str:
        return self.__start

    def getProbes(self) -> list[Probe]:
        return self.__probes

    def getProbe(self, index) -> Probe:
        return self.__probes[index]

    def loadFile(self, inputFile) -> None:
    #load xml file as string
        with open(inputFile, 'r') as file:
            data = file.read()
        self.loadXML(data)


    def printDNA(self) -> None:
        print("Key: ", self.getKey())
        print("Length: ", self.getLength())
        print("Start: ", self.getStart())
        for probe in self.__probes:
            print("Pattern: ", probe.getPattern())
            for cell in probe.getCells():
                print("PosL: ", cell.getPosL(), " PosH: ", cell.getPosH(), " Sequence: ", cell.getSequence())



def getInputFromWeb(n=500, k=10, mode="basic", intensity=0, position=1, sqpe=100, sqne=0, pose=100):
    url = "https://www.cs.put.poznan.pl/pwawrzyniak/bio/bio.php?n=" + str(n) + "&k=" + str(k) + "&mode=" + mode + "&intensity=" + str(intensity) + "&position=" + str(position) + "&sqpe=" + str(sqpe) + "&sqne=" + str(sqne) + "&pose=" + str(pose)
    ## get the data from the website and save to string
    xml = requests.get(url).text
    return xml


if __name__ == "__main__":
    inputFile = "input.xml"
    dna = DNA().loadFile(inputFile)
    dna.printDNA()
