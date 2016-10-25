from networkx import *
from random import randint

import random 
import numpy

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
    
    
    diametro = diameter(graph)
                     
    avgspl = average_shortest_path_length(graph) #TODO same question avg or function
    #Careful with disconnected graphs
    
    dictionary = dict([('cc', clustCoef), ('diameter', diametro), ('avgSPL', avgspl)])
    
    return dictionary

def createGraph(graphType,nodesNr):
    if graphType == 'minimal':
        return minimal_graph(nodesNr)
    if graphType == 'barabasi-albert':
        return barabasi_albert_graph(nodesNr,globalBAedges)

        
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
        