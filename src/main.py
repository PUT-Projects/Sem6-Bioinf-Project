import requests
import dnaParser

if __name__ == "__main__":
    #dna = dnaParser.DNA().loadFile("input.xml")
    dna = dnaParser.DNA().loadXML(dnaParser.getInputFromWeb())
    dna.printDNA()