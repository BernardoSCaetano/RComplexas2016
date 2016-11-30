from networkx import *
from random import randint

import random 
import numpy
import matplotlib.pyplot as plt

globalBAedges = 3
globalProb = 0.1
EQUILIBRIUM = 0.05
EXPERIENCESNR = 10
GLOBALRATES = [0.05,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.6,0.7]


def our_barabasi_albert_graph(n, m):
    G = barabasi_albert_graph(n,m)
    addInfectionAttributeToGraph(G)
    return G        

def our_erdos_renyi_graph(n,p):
    G = erdos_renyi_graph(n,p)
    addInfectionAttributeToGraph(G)
    return G

def setInfection(G,nodeID,value):
    G.node[nodeID]['infected'] = value
    
def isInfected(G,nodeID):
    return G.node[nodeID]['infected'] == True
    
def infectNode(G,nodeID):
    setInfection(G,nodeID,True)
    
def healNode(G,nodeID):
    setInfection(G,nodeID,False)


def addInfectionAttributeToGraph(G):
    for x in range(len(G.nodes())):
        healNode(G,x)
        
def minimal_graph(n):
    if(n<=2):
        raise nx.NetworkXError("minimal network must have n>2")
    else:
        random.seed()
        G=complete_graph(2)
    
        nodeID=2
        while nodeID<n:
            G.add_node(nodeID)
            edges=list(G.edges())
            chosenOnes=random.choice(edges)
            G.add_edge(nodeID,chosenOnes[0])
            G.add_edge(nodeID,chosenOnes[1])
            nodeID+=1        
    
    addInfectionAttributeToGraph(G)    
    return G
    
        
def evalGraph(graph):
    
    #iterates over 1000 trials, choose a node at random, choose two of its neighbors at random, and check if they are connected
    clustCoef = average_clustering(graph) #TODO use average_clustering or clustering function?
    
    #degree centrality for a node v is the fraction of nodes it is connected to
    #degreeCentrality=degree_centrality(graph)
    
    
    #diametro = diameter(graph)                
    #avgspl = average_shortest_path_length(graph) #TODO same question avg or function
   
    #Careful with disconnected graphs
    
    dictionary = dict([('cc', clustCoef), ('diameter', 0), ('avgSPL', 0)])
    
    return dictionary



        
def clusterCoefByDegreeLogLog(graph):
    graus=degree(graph)
    coefs=clustering(graph)
    
    graphic=plt.plot(graus.values(),coefs.values(), 'ro')
    plt.ylabel('log C(k)')
    plt.xlabel('log k')
    plt.yscale('log')
    plt.xscale('log')
    plt.show()
    
    
def clusterCoefByDegree(graph):
    graus=degree(graph)
    coefs=clustering(graph)
    
    graphic=plt.plot(graus.values(),coefs.values(), 'ro')
    plt.ylabel('C(k)')
    plt.xlabel('k')
    
    plt.show()    
    
    
    
def experimentation(graphType,numberOfGraphs,numberOfNodes):
    i=0
    finalDict = dict([('cc',0), ('diameter', 0), ('avgSPL', 0)])
    
    while i < numberOfGraphs:
        graph = createGraph(graphType,numberOfNodes)
        newDict = evalGraph(graph)
        
        newClustCoef = newDict['cc']
        newDiametro = newDict['diameter']
        newAvgspl = newDict['avgSPL']
        
        oldClustCoef = finalDict['cc']
        oldDiametro = finalDict['diameter']
        oldAvgspl = finalDict['avgSPL']
        
        sumClustCoef = newClustCoef + oldClustCoef
        sumDiametro = newDiametro + oldDiametro
        sumAvgspl =  newAvgspl +  oldAvgspl 
        
        finalDict['cc'] = sumClustCoef
        finalDict['diameter'] = sumDiametro
        finalDict['avgSPL'] = sumAvgspl
        
        i=i+1
        
    if finalDict['cc'] != 0:
        finalDict['cc'] = finalDict['cc']/numberOfGraphs
        
    if finalDict['diameter'] != 0:
        finalDict['diameter'] = finalDict['diameter']/numberOfGraphs
        
    if finalDict['avgSPL'] != 0:
        finalDict['avgSPL'] = finalDict['avgSPL']/numberOfGraphs

    
    print "AVERAGES::\nClustering coefficient: %f\nDiameter: %f\nAvgShortestPathLength: %f"%(finalDict['cc'],finalDict['diameter'],finalDict['avgSPL'])
    
    return [finalDict['cc'],finalDict['diameter'],finalDict['avgSPL']]
 
        
def createGraph(graphType,nodesNr):
    if graphType == 'minimal':
        return minimal_graph(nodesNr)
    if graphType == 'barabasi-albert':
        return our_barabasi_albert_graph(nodesNr,globalBAedges)
    if graphType == 'erdos-renyi':
        return our_erdos_renyi_graph(nodesNr,globalProb)

#given a graph, returns a graph with 1 infected node
def startRandomInfection(graph):
    nodes = graph.nodes()
    patientZero = random.choice(nodes)
    infectNode(graph,patientZero)
    return graph


#given a graph, returns a list with the ids of infected nodes
def getInfectedNodes(graph):
    dictionary = get_node_attributes(graph,'infected')
    infectedNodes = list()
    for k,v in dictionary.iteritems():
        if v == True:
            infectedNodes.append(k)
    return infectedNodes
        
#calculates graph of disease for instant t+1 using SI model
def spreadInfectionSI(graph,TransmissionRate):
    infectedNodes = getInfectedNodes(graph)
    for nodeNr in infectedNodes:
        for neighbour in all_neighbors(graph,nodeNr):
            randomNr = random.uniform(0,1)
            if randomNr < TransmissionRate:
                infectNode(graph,neighbour)
    return graph



#calculates graph of disease for instant t+1 using SIS model
#we can optimize: by adding attribute "hasChanged" to each node
#have to set it to 1, for each node who changes from S->I (or vice versa), on  each run of spreadInfectionSIS 
#the above removes the 2 "for" loops required to infect/heal the nodes selected
def spreadInfectionSIS(graph,TransmissionRate,RecoveryRate):
    
    RecoveryRate=1
    
    nodesToInfect = list()
    nodesToHeal = list()
    
    for nodeNr in graph.nodes():
        
        if not(isInfected(graph,nodeNr)):
            for neighbour in all_neighbors(graph,nodeNr):
                if isInfected(graph,neighbour):
                    randomNr1 = random.uniform(0,1)
                    if randomNr1 < TransmissionRate:
                        nodesToInfect.append(nodeNr)
                    break
        if isInfected(graph,nodeNr):
            randomNr2 = random.uniform(0,1)
            if randomNr2 < RecoveryRate:
                nodesToHeal.append(nodeNr)
    
    for nodeSusceptible in nodesToInfect:
        infectNode(graph,nodeSusceptible)
    
    for nodeInfected in nodesToHeal:
        healNode(graph,nodeInfected)
        
    return graph
                
                
#returns the fraction of infected nodes, given a graph
def getFractionOfInfected(graph):
    numberOfNodes=len(graph.nodes())
    numberOfInfected=len(getInfectedNodes(graph))
    
    return float(numberOfInfected)/numberOfNodes
    
    
#True if 2 numbers are close    
def areNumbersClose(n,m):
    if abs(n-m) < EQUILIBRIUM:
        return True
    return False


#True if the n-last elements of lista are close to the first one in the list
def hasStabilized(lista,n):
    if len(lista) <= n:
        return False
    
    lista=lista[-n:]
    for elem in lista[1:]:
        if not(areNumbersClose(lista[0],elem)):
            return False     
    return True
   
   
#returns the average fraction of infected after the spreading stabilizes    
def fractionInfectedAfterStabilizing(graph_model,TransmissionRate):
    G = createGraph(graph_model,1000)
    
    infectedFractions = list()
    infectedFractions.append(0) 
    G = startRandomInfection(G)
    infectedFractions.append(getFractionOfInfected(G))
    
    timeline=list()
    timeline.append(0)
    timeline.append(1)
    t=2
    while not(hasStabilized(infectedFractions,100)):
        G = spreadInfectionSIS(G,TransmissionRate,1)
        infectedFractions.append(getFractionOfInfected(G))
        timeline.append(t)
        t=t+1
        if t > 800:
            break
    
    last100 = infectedFractions[-100:]
    averageInfected = numpy.mean(last100)
    
    
    infectedFractions = infectedFractions[:-100]
    infectedFractions.append(averageInfected)
    
    #graphic=plt.plot(timeline[:-99],infectedFractions, 'ro')
    #plt.ylabel('Infected Fraction')
    #plt.xlabel('t')    
    #plt.show()
    return averageInfected


def guessBetween(n,m):
    return round(random.uniform(n,m),4)


#True if the number is really close to 0
def collapses(n):
    if round(n,3)<=0.001:
        return True
    return False
    
    
#deprecated
def calculateTimeRequired():
    G = createGraph('minimal',100)
    G = startRandomInfection(G)
    i=0
    sizes=list()
    while(i<100):
        i=i+1
        size = len(getInfectedNodes(G))
        sizes.append(size)
        G = spreadInfectionSIS(G,0.5,0.5)    
        
    average = numpy.mean(sizes)
    
    for index,elem in enumerate(sizes):#sorted?
        if abs(elem-average)<3:
            firstIndex = index
    
    print average
    print index
    print sizes
    print sizes[99]
   
   
#deprecated    
def main1():
    ccs=[]
    diameters=[]
    avgSPLs=[]
    
    #clusterCoefByDegreeLogLog(createGraph('minimal',100))
    #clusterCoefByDegreeLogLog(createGraph('minimal',1000))
    #clusterCoefByDegreeLogLog(createGraph('minimal',10000))

        
        #diameters.extend(experimentation('minimal',10,100)[1])
        #avgSPLs.extend(experimentation('minimal',10,100)[1])


def testOneStepSIS():
    G = createGraph('erdos-renyi',10)
    startRandomInfection(G)
    print "patient zero: " + str(getInfectedNodes(G))
    print "edges: " + str(G.edges())
    spreadInfectionSIS(G,1,1) #testing with rate = 1
    print "after one step infection: "+str(getInfectedNodes(G))          
    

def testOneStepSI():
    G = createGraph('erdos-renyi',10)
    startRandomInfection(G)
    print "patient zero: " + str(getInfectedNodes(G))
    print "edges: " + str(G.edges())
    spreadInfectionSI(G,1) #testing with rate = 1
    print "after one step infection: "+str(getInfectedNodes(G))    


#returns a pair of lists, cointaining the rates and the experimental(experienceNR times) fraction of infected for each rate
def infectedFractionByTransmissionRate(graph_model):
    rates=GLOBALRATES
    fractionByRate=list()
    averages=list()
    experiments=0
    
    for rate in rates:       
        while experiments < EXPERIENCESNR:
            experiments=experiments+1
            averages.append(fractionInfectedAfterStabilizing(graph_model,rate))
            pointInPlot=numpy.mean(averages)
        
        fractionByRate.append(pointInPlot)
        experiments=0
        averages=[]
        print "experienced for "+str(EXPERIENCESNR)+ " instances with transmission rate "+ str(rate)
        
    graphic=plt.plot(rates,fractionByRate, 'ro')
    plt.ylabel('Infected Fraction')
    plt.xlabel('TransmissionRate')    
    plt.show()
    return rates, fractionByRate


#returns the experimental approximation of the threshold for a given model
def calcThreshold(graph_model):
    rates, fractions = infectedFractionByTransmissionRate(graph_model)
    limsup=0
    liminf=0
    
    for i,fraction in enumerate(fractions):
        if fraction < 0.1 and fraction > 0.008:
            limsup=rates[i]
            break
    
    while abs(liminf-limsup) > 0.0005:
        rateHypothesis = guessBetween(liminf,limsup)
        fraction = fractionInfectedAfterStabilizing(graph_model,rateHypothesis)
        if collapses(fraction):
            liminf=rateHypothesis
        if not(collapses(fraction)):
            limsup=rateHypothesis
        print "[ "+str(liminf)+" , "+str(limsup)+" ]"
    
    threshold=float(limsup+liminf)/2
    print "threshold: "+ str(threshold)
    return threshold
        
    
calcThreshold('barabasi-albert')