from networkx import *
from random import randint

import random 
import numpy
import matplotlib.pyplot as plt

globalBAedges=3

def _random_subset(seq,m):
    """ Return m unique elements from seq.

    This differs from random.sample which can return repeated
    elements if seq holds repeated elements.
    """
    targets=set()
    while len(targets)<m:
        x=random.choice(seq)
        targets.add(x)
    return targets


def barabasi_albert_graph(n, m, seed=None):
    """Returns a random graph according to the Barabási–Albert preferential
    attachment model.

    A graph of ``n`` nodes is grown by attaching new nodes each with ``m``
    edges that are preferentially attached to existing nodes with high degree.

    Parameters
    ----------
    n : int
        Number of nodes
    m : int
        Number of edges to attach from a new node to existing nodes
    seed : int, optional
        Seed for random number generator (default=None).

    Returns
    -------
    G : Graph

    Raises
    ------
    NetworkXError
        If ``m`` does not satisfy ``1 <= m < n``.

    References
    ----------
    .. [1] A. L. Barabási and R. Albert "Emergence of scaling in
       random networks", Science 286, pp 509-512, 1999.
    """
    if m < 1 or  m >=n:
        raise nx.NetworkXError("Barabási–Albert network must have m >= 1"
                               " and m < n, m = %d, n = %d" % (m, n))
    if seed is not None:
        random.seed(seed)
    
    # Add m initial nodes (m0 in barabasi-speak)
    G=empty_graph(m)
    G.name="barabasi_albert_graph(%s,%s)"%(n,m)
    # Target nodes for new edges
    targets=list(range(m))
    # List of existing nodes, with nodes repeated once for each adjacent edge
    repeated_nodes=[]
    # Start adding the other n-m nodes. The first node is m.
    source=m
    while source<n:
        # Add edges to m nodes from the source.
        G.add_edges_from(zip([source]*m,targets))
        # Add one node to the list for each new edge just created.
        repeated_nodes.extend(targets)
        # And the new node "source" has m edges to add to the list.
        repeated_nodes.extend([source]*m)
        # Now choose m unique nodes from the existing nodes
        # Pick uniformly from repeated_nodes (preferential attachement)
        targets = _random_subset(repeated_nodes,m)
        source += 1
    addInfectionToGraph(G)
    return G

        

def setInfection(G,nodeID,value):
    G.node[nodeID]['infected'] = value
    
def infectNode(G,nodeID):
    setInfection(G,nodeID,True)
    
def healNode(G,nodeID):
    setInfection(G,nodeID,False)


def addInfectionToGraph(G):
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
            edges=G.edges()
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

def createGraph(graphType,nodesNr):
    if graphType == 'minimal':
        return minimal_graph(nodesNr)
    if graphType == 'barabasi-albert':
        return barabasi_albert_graph(nodesNr,globalBAedges)

        
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
        
def main():
    ccs=[]
    diameters=[]
    avgSPLs=[]
    
    #clusterCoefByDegreeLogLog(createGraph('minimal',100))
    #clusterCoefByDegreeLogLog(createGraph('minimal',1000))
    #clusterCoefByDegreeLogLog(createGraph('minimal',10000))

        
        #diameters.extend(experimentation('minimal',10,100)[1])
        #avgSPLs.extend(experimentation('minimal',10,100)[1])
        
    

main()

#cc_by_node[0.7265637348008186, 0.7365163316339398, 0.7359093143372706, 0.7358060601107528, 0.7402482954234243, 0.7354239628654052, 0.7350685263852069, 0.7378760587088745, 0.7346786156776204, 0.7390508706152371, 0.7367485248779614, 0.7388354161471171, 0.7376105016939264, 0.7396181868942342, 0.738036211678682, 0.7380555784813824, 0.7393140437301213, 0.7383356255055944, 0.7370154626660501, 0.7384302696892411, 0.7387004232057347, 0.7384176587694157, 0.7380755343022549, 0.7383797304451774, 0.7377326878099811, 0.7390769433793558, 0.7389729130056194, 0.7395671400080576, 0.7385502870518751, 0.7401749713501835]

#nodes=(100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,2100,2200,2300,2400,2500,2600,2700,2800,2900,3000)