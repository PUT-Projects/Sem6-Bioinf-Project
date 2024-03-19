## File for parsing input xml file

## <dna key="22408638" length="500" start="GAGTATTTCT">
# <probe pattern="NNNNNNNNNN">
# <cell posL="250" posH="350">AAAAAGTGTG</cell>
# <cell posL="391" posH="491">AAAAATTATC</cell>
# <cell posL="370" posH="470">AAAAATTGTT</cell>

# store it in a list of objects
import xml.etree.ElementTree as ET

import requests

# The `Probe` class has attributes for pattern and cells, with methods to retrieve the probe, pattern,
# and cells.
class Probe:
    def __init__(self, pattern, cells):
        self.__pattern = pattern
        self.__cells = cells

# Getters
    def getProbe(self):
        # The `return self` statement in the methods of the classes (`Probe`, `Cell`, and `DNA`) is
        # used to return the instance of the class itself.
        return self

    def getPattern(self):
        return self.__pattern

    def getCells(self):
        return self.__cells


class Cell:
    def __init__(self, posL, posH, sequence):
        self.__posL = posL
        self.__posH = posH
        self.__sequence = sequence

# Getters
    def getCell(self):
        return self
    def getPosL(self):
        return self.__posL
    def getPosH(self):
        return self.__posH
    def getSequence(self):
        return self.__sequence


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

    def loadFile(self, inputFile):
    #load xml file as string
        with open(inputFile, 'r') as file:
            data = file.read()
        self.loadXML(data)


    def printDNA(self):
        print("Key: ", self.key)
        print("Length: ", self.length)
        print("Start: ", self.start)
        for probe in self.probes:
            print("Pattern: ", probe.pattern)
            for cell in probe.cells:
                print("PosL: ", cell.posL, " PosH: ", cell.posH, " Sequence: ", cell.sequence)


# Getters
    def getDNA(self):
        return self

    def getKey(self):
        return self.__key

    def getLength(self):
        return self.__length

    def getStart(self):
        return self.__start


def getInputFromWeb(n=500, k=10, mode="basic", intensity=0, position=1, sqpe=100, sqne=0, pose=100):
    url = "https://www.cs.put.poznan.pl/pwawrzyniak/bio/bio.php?n=" + str(n) + "&k=" + str(k) + "&mode=" + mode + "&intensity=" + str(intensity) + "&position=" + str(position) + "&sqpe=" + str(sqpe) + "&sqne=" + str(sqne) + "&pose=" + str(pose)
    ## get the data from the website and save to string
    xml = requests.get(url).text
    return xml


if __name__ == "__main__":
    inputFile = "input.xml"
    dna = DNA().loadFile(inputFile)
    dna.printDNA()
