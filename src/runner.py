#run analysis

import dnaParser
import solver
import solverHeuristic
import time



def runHeuristic(n,k,sqpep,posep):
    try:
        webxml = dnaParser.getInputFromWeb(n, k, sqpep=sqpep, posep=posep)
        dna = dnaParser.DNA().loadXML(webxml)
    except Exception as e:
        print(e)
        print(n,k,sqpep,posep)
        print(webxml)
        return None, None

    start = time.time()
    graph = solverHeuristic.Graph(dna.getProbes()[0], dna.getStart(), dna.getLength())
    graph.sortGraphAdjacencies()
    #graph.printGraph()
    sts = graph.sts()

    seq = graph.getSequence(sts)
    end = time.time()

    time_taken = end-start

    return seq, time_taken


def runExact(n,k,sqpep,posep):

    try:
        webxml = dnaParser.getInputFromWeb(n, k, sqpep=sqpep, posep=posep)
        dna = dnaParser.DNA().loadXML(webxml)
    except Exception as e:
        print(e)
        print(n,k,sqpep,posep)
        print(webxml)
        return None, None


    start = time.time()

    graph = solver.Graph(dna.getProbes()[0], dna.getStart(), dna.getLength())
    graph.sortGraphAdjacencies()
    #graph.printGraph()
    sts = graph.sts()

    seq = graph.getSequence(sts)

    end = time.time()
    time_taken = end-start

    return seq, time_taken

def writeToFile(filename, param, time, correct):
    mstime = round(time*1000, 2)
    with open('output/' + filename, 'a') as file:
        file.write(str(param) + " " + str(mstime) + " " + str(correct) + "\n")



def heuristicLoop():
    #Default values
    n = 100
    k = 8
    posep = 10
    sqpep = 10

    #For variable n
    n_table = [20,25,30,40, 50 ,100, 500]
    for i in n_table:
        s,t = runHeuristic(i, k, sqpep, posep)
        #save to file n,t and if len(s) = n
        writeToFile("heuristic_results_n.txt", i, t, len(s) == i)

    #For variable k
    k_table = [6,7,8,9,10]
    for i in k_table:
        s,t = runHeuristic(n, i, sqpep, posep)
        writeToFile("heuristic_results_k.txt", i, t, len(s) == n)

    #For variable sqpep
    sqpep_table = [5,10,15,20,25]
    for i in sqpep_table:
        s,t = runHeuristic(n, k, i, posep)
        writeToFile("heuristic_results_sqpep.txt", i, t, len(s) == n)

    #For variable posep
    posep_table = [5,10,15,20,25, 30,40,50]
    for i in posep_table:
        s,t= runHeuristic(n, k, sqpep, i)
        writeToFile("heuristic_results_posep.txt", i, t, len(s) == n)


def exactLoop():
    #Default values
    n = 25
    k = 8
    posep = 10
    sqpep = 10

    #For variable n
    n_table = [ 20,25,30, 35]
    for i in n_table:
        s,t = runExact(i, k, sqpep, posep)
        #save to file n,t and if len(s) = n
        writeToFile("exact_results_n.txt", i, t, len(s) == i)

    #For variable k
    k_table = [4,5,6,7,8,9,10]
    for i in k_table:
        s,t = runExact(n, i, sqpep, posep)
        writeToFile("exact_results_k.txt", i, t, len(s) == n)

    #For variable sqpep
    sqpep_table = [5,10,15,20,25]
    for i in sqpep_table:
        s,t = runExact(n, k, i, posep)
        writeToFile("exact_results_sqpep.txt", i, t, len(s) == n)

    #For variable posep
    posep_table = [5,10,15,20,25,30,40,50]
    for i in posep_table:
        s,t= runExact(n, k, sqpep, i)
        writeToFile("exact_results_posep.txt", i, t, len(s) == n)


def main():
    heuristicLoop()
    exactLoop()

if __name__ == "__main__":
    main()