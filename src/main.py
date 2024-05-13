import requests
import dnaParser

def main():
        #dna = dnaParser.DNA().loadFile("input.xml")
    dna = dnaParser.DNA().loadXML(dnaParser.getInputFromWeb())
    dna.printDNA()

if __name__ == "__main__":
    main()