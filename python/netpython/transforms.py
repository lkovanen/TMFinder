import pynet,netext,percolator
import random
import numpy as np

def mst(net,maximum=False):
    """Find a minimum/maximum spanning tree

    """
    return mst_kruskal(net,True,maximum)

def mst_kruskal(net,randomize=True,maximum=False):
    """Find a minimum/maximum spanning tree using Kruskal's algorithm

    If random is set to true and the mst is not unique, a random
    mst is chosen.

    >>> t=pynet.SymmNet()
    >>> t[1,2]=1
    >>> t[2,3]=2
    >>> t[3,1]=3
    >>> m=mst_kruskal(t)
    >>> print m.edges
    [[1, 2, 1], [2, 3, 2]]
    """
    edges=list(net.edges)
    if randomize:
        random.shuffle(edges) #the sort has been stable since python version 2.3
    edges.sort(lambda x,y:cmp(x[2],y[2]),reverse=maximum)
    mst=pynet.SymmNet()
    numberOfNodes=len(net)
    ktree=percolator.Ktree() #just use dict
    addedEdges=0

    #First add the nodes
    for node in net:
        mst.addNode(node)
    
    #Add the edges
    for edge in edges:
        if ktree.getParent(edge[0])!=ktree.getParent(edge[1]):
            mst[edge[0],edge[1]]=edge[2]
            ktree.setParent(edge[0],edge[1])
            addedEdges+=1

        if addedEdges==numberOfNodes-1:
            #the mst is a tree
            netext.copyNodeProperties(net,mst)
            return mst

    # else it is a forest
    netext.copyNodeProperties(net,mst)
        
    return mst


def snowball(net, seed, depth, includeLeafEdges=False):
    """Snowball sampling

    Works for both directed and undirected networks. For directed
    networks all edges all followed during the sampling (as opposed to
    following only outbound edges).

    Parameters
    ----------
    net : pynet.SymmNet or pynet.Net object
        The network to be sampled.
    seed : int or a sequence of ints
        The seed of the snowball, either a single node index or
        several indices.
    depth : int
        The depth of the snowball. Depth 1 corresponds to first
        neighbors of the seed only.
    includeLeafEdges : bool (default: False)
        If True, then the edges between the leaves (i.e. the nodes at
        final depth) will also be included in the snowball network. By
        default these edges are not included.

    Return
    ------
    snowball : pynet.SymmNet or pynet.Net object
        The snowball sample, will be of the same type as `net`.
    """
    if isinstance(seed, int):
        seed = [seed]
    toVisit=set(seed)

    # Create a network for the sample with the same type as `net`.
    newNet=type(net)()
    visited=set()

    for d in range(1,depth+1):
        #print "Depth: ",d," visited ", len(visited)," to visit ", len(toVisit)
        visited=visited|toVisit
        newToVisit=set()

        if len(toVisit) == 0:
            break
        
        for nodeIndex in toVisit:
            node = net[nodeIndex]

            # Go through outbound edges (this equals all neighbors in
            # an undirected network.
            for outIndex in node.iterOut():
                newNet[nodeIndex][outIndex] = net[nodeIndex][outIndex]
                if outIndex not in visited:
                    newToVisit.add(outIndex)

            # If we are dealing with a directed network, then we must
            # also go through the inbound edges.
            if isinstance(net, pynet.Net):
                for inIndex in node.iterIn():
                    newNet[inIndex][nodeIndex] = net[inIndex][nodeIndex]
                    if inIndex not in visited:
                        newToVisit.add(inIndex)

        # If this is the last depth and `includeLeafEdges` is
        # True, we add the edges between the most recently added
        # nodes, that is, those currently in the set `newToVisit`.
        if d == depth and includeLeafEdges:
            for nodeIndex in newToVisit:
                node = net[nodeIndex]

                for outIndex in node.iterOut():
                    if outIndex in newToVisit:
                        newNet[nodeIndex][outIndex] = net[nodeIndex][outIndex]

                if isinstance(net, pynet.Net):
                    for inIndex in node.iterIn():
                        if inIndex in newToVisit:
                            newNet[inIndex][nodeIndex] = net[inIndex][nodeIndex]

        # The nodes to be visited on the next round are the leaves
        # found in the current round.
        toVisit=newToVisit

    netext.copyNodeProperties(net,newNet)

    return newNet


def collapseIndices(net, returnIndexMap=False):
    """Changes the indices of net to run from 0 to len(net)-1.
    """

    newNet = type(net)()
    indexmap = {}
    index = 0

    for i in net:
        newNet.addNode(index)
        indexmap[i] = index;
        index += 1

    for edge in net.edges:
        i,j,w=edge
        newNet[indexmap[i]][indexmap[j]] = w

    netext.copyNodeProperties(net,newNet)

    if returnIndexMap:
        return newNet, indexmap
    else:
        return newNet


def threshold_by_value(net,threshold,accept="<",keepIsolatedNodes=False):
    '''Generates a new network by thresholding the input network. 
       If using option keepIsolatedNodes=True, all nodes in the
       original network will be included in the thresholded network;
       otherwise only those nodes which have links will remain (this
       is the default). 
    
       Inputs: net = network, threshold = threshold value,
       accept = "foobar": accept weights foobar threshold (e.g accept = "<": accept weights < threshold)
       Returns a network.'''

    newnet=pynet.SymmNet()
    edges=list(net.edges)
    if accept == "<":
        for edge in edges:
            if (edge[2] < threshold):
                newnet[edge[0],edge[1]]=edge[2]
    elif accept == ">":
        for edge in edges:
            if (edge[2] > threshold):
                newnet[edge[0],edge[1]]=edge[2] 
    elif accept == ">=":
        for edge in edges:
            if (edge[2] >= threshold):
                newnet[edge[0],edge[1]]=edge[2] 
    elif accept == "<=":
        for edge in edges:
            if (edge[2] <= threshold):
                newnet[edge[0],edge[1]]=edge[2] 
    else:
        raise Exception("Parameter 'accept' must be either '<', '>', '<=' or '>='.")

    # Add isolated nodes to the network.
    if keepIsolatedNodes==True:
        for node in net:
            if not newnet.__contains__(node):
                newnet.addNode(node)
            
    netext.copyNodeProperties(net,newnet)
                
    return newnet


def dist_to_weights(net,epsilon=0.001):
    '''Transforms a distance matrix / network to a weight
    matrix / network using the formula W = 1 - D / max(D)+epsilon.
    Returns a matrix/network'''

    N=len(net._nodes)

    if (isinstance(net,pynet.SymmFullNet)):
        newmat=pynet.SymmFullNet(N)
    else:
        newmat=pynet.SymmNet()    

    edges=list(net.edges)

    maxd=0.0
    for edge in edges:
        if edge[2]>maxd:
            maxd=edge[2]

    # epsilon trick; lowest weight will be almost but
    # not entirely zero    
    for edge in edges:
        newmat[edge[0]][edge[1]]=1-edge[2]/maxd+epsilon

    netext.copyNodeProperties(net,newmat)

        
    return newmat

def filterNet(net,keep_these_nodes):
    return getSubnet(net,keep_these_nodes)

def getSubnet(net,nodes):
    """Get induced subgraph.

    Parameters
    ----------
    net: pynet.Net, pynet.SymmNet or pynet.SymmFullNet
        The original network.
    nodes : sequence
        The nodes that span the induces subgraph.

    Return
    ------
    subnet : type(net)
        The induced subgraph that contains only nodes given in
        `nodes` and the edges between those nodes that are
        present in `net`. Node properties etc are left untouched.
    """
    # Handle both directed and undirected networks.
    newnet = type(net)() # Initialize to same type as `net`.
    degsum=0
    for node in nodes:        
        degsum += net[node].deg()
        newnet.addNode(node)
    if degsum >= len(nodes)*(len(nodes)-1)/2:
        othernodes=set(nodes)
        for node in nodes:
            if net.isSymmetric():
                othernodes.remove(node)
            for othernode in othernodes:
                if net[node,othernode]!=0:
                    newnet[node,othernode]=net[node,othernode]
    else:
        for node in nodes:
            for neigh in net[node]:
                if neigh in nodes:
                    newnet[node,neigh]=net[node,neigh]

    netext.copyNodeProperties(net, newnet)

    return newnet


    


def collapseBipartiteNet(net,nodesToRemove):
    """
    Returns an unipartite projection of a bipartite network.
    """
    newNet=pynet.SymmNet()
    for node in nodesToRemove:
        degree=float(net[node].deg())
        for node1 in net[node]:
            for node2 in net[node]:
                if node1.__hash__()>node2.__hash__():
                    newNet[node1,node2]=newNet[node1,node2]+1.0/degree

    netext.copyNodeProperties(net,newNet)
    return newNet

def local_threshold_by_value(net,threshold):
    '''Generates a new network by thresholding the input network.
       Inputs: net = network, threshold = threshold value,
       mode = 0 (accept weights < threshold), 1 (accept weights > threshold)
       Returns a network. Note! threshold is really alpha which is defined in
       "Extracting the multiscale backbone of complex weighted networks"
       http://www.pnas.org/content/106/16/6483.full.pdf'''

    newnet=pynet.SymmNet()
    for node in net:
        s=net[node].strength()
        k=net[node].deg()
        for neigh in net[node]:
            w=net[node,neigh]
            if (1-w/s)**(k-1)<threshold:
                newnet[node,neigh]=w
    netext.copyNodeProperties(net,newnet)

    return newnet

def getLineGraph(net, useWeights=False, output=None, format='edg'):
    """Return a line graph constructed from `net`.

    The nodes in the line graph correspond to edges in the original
    graph, and there is an edge between two nodes if they have a
    common incident node in the original graph. 

    If weights are not used (`useWeights = False`), the resulting
    network will be undirected and the weight of each new edge will be
    1/(k_i-1), where k_i is the degree of the common node in `net`.

    If weights are used (`useWeights = True`), the resulting network
    will be directed and the weight of edge (e_ij, e_jk) will be
    w_jk/sum_{x != i} w_jx, where the indices i, j and k refer to
    nodes in `net`.

    Parameters
    ----------
    net : pynet.SymmNet object
        The original graph that is used for constructing the line
        graph.
    useWeights : boolean
        If True, the edge weights will be used when constructing the
        line graph.
    output : file object
        If given, the edges will be written to output in edg-format
        instead of returning a pynet.Net() or pynet.SymmNet() object.
    format : str, 'edg' or 'net'
        If `output` is specified, `format` specifies how the output is
        written. 'edg' is the standard edge format (FROM TO WEIGHT)
        and 'net' gives the Pajek format.

    Return
    ------
    IF `output` is None:
        linegraph : pynet.SymmNet or pynet.Net object
            The weighted line graph.
    id_array : numpy.array with shape (len(net.edges), 2)
        Array for converting the nodes in the line graph back into the
        edges of the original graph. id_array[EDGE_ID] contains the
        two end nodes of given edge, where EDGE_ID is the same as used
        in `linegraph`.
    """

    if output is None:
        if useWeights:
            linegraph = pynet.Net()
        else:
            linegraph = pynet.SymmNet()
    edge_map = dict() # edge_map[sorted([n_i, n_j])] = new_node_ID

    if output is not None and format == 'net':
        # Print Pajek file header.
        N_edges = len(list(net.edges))
        output.write("*Vertices %d\n" % N_edges)
        for i in range(N_edges):
            output.write('%d "%d"\n' % (i, i))
        N_edge_links = 0
        for n in net:
            degree = len(list(net[n]))
            N_edge_links += (degree*(degree-1))/2
        if useWeights:
            output.write("*Arcs %d\n" % (2*N_edge_links,))
        else:
            output.write("*Edges %d\n" % N_edge_links)

    # Go through all nodes (n_c = center node), and for each node, go
    # through all pairs of neighbours (n_i and n_j). The edges
    # e_i=(n_c,n_i) and e_j=(n_c,n_j) are nodes in the line graph, so
    # we add a link between them.
    for n_c in net:
        strength = net[n_c].strength()
        nb = list(net[n_c]) # List of neighbours
        for i, n_i in enumerate(nb):
            e_i = edge_map.setdefault(tuple(sorted([n_c,n_i])), len(edge_map))
            other_nb = (nb[:i]+nb[i+1:] if useWeights else nb[i+1:])
            for n_j in other_nb:
                e_j = edge_map.setdefault(tuple(sorted([n_c,n_j])), len(edge_map))
                if useWeights:
                    w = net[n_c][n_j]/(strength - net[n_c][n_i])
                else:
                    w = 1.0/(len(nb)-1)
                if output is None:
                    linegraph[e_i][e_j] = w
                else:
                    output.write(" ".join(map(str, [e_i, e_j, w])) + "\n")

    # Construct id_array from edge_map
    id_array = np.zeros((len(edge_map), 2), int)
    for node_pair, edgeID in edge_map.iteritems():
        id_array[edgeID] = list(node_pair)

    if output is None:
        return linegraph, id_array
    else:
        return id_array


def netConfiguration(net, keepsOrigNet=False, seed=None):
    """Generate configuration network
    
    This function generates a configuration network from any arbitrary
    net. It retains the degree of each node but randomize the edges
    between them.
		
    Parameters
    ----------
    net : pynet.SymmNet object
        The network to be used as the basis for the configuration
        model.
    keepsOrigNet : bool (default: False)
        If False, the input network, `net`, will be overwritten by the
        configuration network.
    seed : int (default: None)
        A seed for the random number generator. If None, the RNG is
        not be re-initialized but the current state is used.

    Return
    ------
    configuration_net : pynet.SymmNet object
        The shuffled network. Note that if `keepsOrigNet` is False,
        the returned value will be identical to `net`.
    """
    if seed is not None:
        random.seed(int(seed))

    newNet = pynet.SymmNet()
    if keepsOrigNet:
        testNet = pynet.SymmNet()
        for edge in net.edges:
            testNet[edge[0],edge[1]] = edge[2]	
    else:
        testNet=net

    edgeList = list(net.edges)

    for i in range(len(edgeList)):

        j=i
        while j==i:
            j=random.randint(0,len(edgeList)-1)

        if ((edgeList[i][1]==edgeList[j][0]) 
            or (edgeList[j][1]==edgeList[i][0])):
            continue
        if ((edgeList[i][1]==edgeList[j][0]) 
            and (edgeList[j][1]==edgeList[i][0])):
            continue
        if ((edgeList[i][0]==edgeList[j][0]) 
            or (edgeList[i][1]==edgeList[j][1])):
            continue
        if ((newNet[edgeList[i][0],edgeList[j][1]]>0.0) 
            or (newNet[edgeList[j][0],edgeList[i][1]]>0.0)):
            continue
        if ((testNet[edgeList[i][0],edgeList[j][1]]>0.0) 
            or (testNet[edgeList[j][0],edgeList[i][1]]>0.0)):
            continue

        edgeList[i][1]+=edgeList[j][1]
        edgeList[j][1]=edgeList[i][1]-edgeList[j][1]
        edgeList[i][1]=edgeList[i][1]-edgeList[j][1]

        newNet[edgeList[i][0],edgeList[j][1]]=0.0
        newNet[edgeList[j][0],edgeList[i][1]]=0.0

        testNet[edgeList[i][0],edgeList[j][1]]=0.0
        testNet[edgeList[j][0],edgeList[i][1]]=0.0

        newNet[edgeList[i][0],edgeList[i][1]]=1.0
        newNet[edgeList[j][0],edgeList[j][1]]=1.0
        testNet[edgeList[i][0],edgeList[i][1]]=1.0
        testNet[edgeList[j][0],edgeList[j][1]]=1.0


    return newNet


def copyNet(net):
    """ Copy a network. Also the node properties are copied.

    Parameters
    ----------
    net : any pynet object
        The network to be copied.

    Return
    ------
    newNet : net.__class__ object
        A copy of the network.

    """
    newNet=net.__copy__()
    netext.copyNodeProperties(net,newNet)
    return newNet

if __name__ == '__main__':
    """Run unit tests if called."""
    from tests.test_transforms import *
    unittest.main()
