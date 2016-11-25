from networkx import *
from random import randint

import random 
import numpy
import matplotlib.pyplot as plt

globalBAedges=3
globalProb=0.5

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
def spreadInfectionSIS(graph,TransmissionRate,RecoveryRate):
    infectedNodes = getInfectedNodes(graph)
    for nodeNr in infectedNodes:
        for neighbour in all_neighbors(graph,nodeNr):
            randomNr1 = random.uniform(0,1)
            if randomNr1 < TransmissionRate:
                infectNode(graph,neighbour)
            
        randomNr2 = random.uniform(0,1)
        if randomNr2 < RecoveryRate:
            healNode(graph,nodeNr)
            
    return graph
                
                
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
    G = createGraph('erdos-renyi',7)
    startRandomInfection(G)
    print "patient zero: " + str(getInfectedNodes(G))
    print "edges: " + str(G.edges())
    spreadInfectionSIS(G,1,1) #testing with rate = 1
    print "after infection: "+str(getInfectedNodes(G))          
    

def testOneStepSI():
    G = createGraph('erdos-renyi',7)
    startRandomInfection(G)
    print "patient zero: " + str(getInfectedNodes(G))
    print "edges: " + str(G.edges())
    spreadInfectionSI(G,1) #testing with rate = 1
    print "after infection: "+str(getInfectedNodes(G))    

testOneStepSIS()

#cc_by_node[0.7265637348008186, 0.7365163316339398, 0.7359093143372706, 0.7358060601107528, 0.7402482954234243, 0.7354239628654052, 0.7350685263852069, 0.7378760587088745, 0.7346786156776204, 0.7390508706152371, 0.7367485248779614, 0.7388354161471171, 0.7376105016939264, 0.7396181868942342, 0.738036211678682, 0.7380555784813824, 0.7393140437301213, 0.7383356255055944, 0.7370154626660501, 0.7384302696892411, 0.7387004232057347, 0.7384176587694157, 0.7380755343022549, 0.7383797304451774, 0.7377326878099811, 0.7390769433793558, 0.7389729130056194, 0.7395671400080576, 0.7385502870518751, 0.7401749713501835]

#nodes=(100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,2100,2200,2300,2400,2500,2600,2700,2800,2900,3000)