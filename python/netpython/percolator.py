"""
This module contains percolator functions and some related data structures.
This contains general framework for edge and node percolation studies, but
also functions for k-clique percolation.
"""

import pynet,netext,array,math,netio,communities,numpy,transforms
from operator import mul



class KtreeInteger:
    def __init__(self,size=0):
        """
        Parameters
        ----------
        size : int
            Size of the Kruskal tree. If size is set to 0, the Kruskal tree is
            considered of having a dynamic size, and it can be made larger
            by using the setSize methods. For positive values of size, the
            Kruskal tree size is static. 
        """
        if size>0:
            self.ktree=numpy.ones(size,dtype="uint")
            self.subTreeWeight=numpy.ones(size,dtype="uint")
        else:
            self.ktree=[]
            self.subTreeWeight=[]

        self.mappingOn=False
        self.sizeDistribution={}
        self.size=size
        if size!=0:
            self.sizeDistribution[1]=size        
            for index in xrange(0,size):
                self.ktree[index]=index;
            self.giantSize=1
        else:
            self.giantSize=0
        
        self.suscSum=size #The sum of the squared sizes in for the susceptibility

    def __getRealParent(self,node):
        """
        Private method. Reads elements directly from the tree.        
        """
        try:
            return self.ktree[node]
        except IndexError:
            self.setSize(node)
            return node

    def __setRealParent(self,node,newParent):
        """
        Private.
        """
        try:
            self.ktree[node]=newParent
        except IndexError:
            self.setSize(node)
            self.ktree[node]=newParent        

    def getSetIndex(self,node):
        parent=self.__getRealParent(node)
        if node!=parent:
            self.__setRealParent(node,self.getSetIndex(parent))
        return self.__getRealParent(node)

            
    def mergeSets(self,node1,node2):
        set_of_node1=self.getSetIndex(node1)
        set_of_node2=self.getSetIndex(node2)
        if set_of_node1!=set_of_node2:
            if self.subTreeWeight[set_of_node1]>self.subTreeWeight[set_of_node2]:
                large_set,small_set=set_of_node1,set_of_node2
            else:
                large_set,small_set=set_of_node2,set_of_node1 

            small_set_size=self.subTreeWeight[small_set]
            large_set_size=self.subTreeWeight[large_set]
            self.sizeDistribution[small_set_size]-=1
            self.sizeDistribution[large_set_size]-=1
            #we remove empty elements
            if self.sizeDistribution[small_set_size]==0:
                self.sizeDistribution.__delitem__(small_set_size)
            if small_set_size!=large_set_size: 
                # If sets are equal size, "larger" set size is already
                # removed.
                if self.sizeDistribution[large_set_size]==0:
                    self.sizeDistribution.__delitem__(large_set_size)
            
            self.sizeDistribution[small_set_size+large_set_size] = self.sizeDistribution.get(small_set_size+large_set_size,0)+1
            self.subTreeWeight[large_set]+=self.subTreeWeight[small_set]

            new_set_size=small_set_size+large_set_size
            if new_set_size>self.giantSize:
                self.giantSize=new_set_size
            self.suscSum-=small_set_size*small_set_size
            self.suscSum-=large_set_size*large_set_size
            self.suscSum+=new_set_size*new_set_size
            assert self.suscSum>1

            self.__setRealParent(small_set,large_set)

    def getSusceptibility(self):
        if self.size!=self.giantSize:
            return (self.suscSum-self.giantSize*self.giantSize)/float(self.size-self.giantSize)
        else:
            return 0.0
        
    def mergeSetsWithElements(self,elements):
        first=elements[0]
        for i in range(1,len(elements)):
            self.setParent(elements[i],first)

    def __iter__(self):
        for i in self.ktree:
            yield i

    def getCommStruct(self,separateElements=True):
        communityMap={}
        if self.mappingOn:
            nodes=self.ktree
        else:
            nodes=range(0,len(self.ktree))
        
        for node in nodes:
            communityKey=self.getSetIndex(node)
            if separateElements or communityKey!=node:
                if communityKey not in communityMap:
                    communityMap[communityKey]=[node]
                else:
                    communityMap[communityKey].append(node)

        return communities.NodePartition(communityMap)

    def __len__(self):
        return len(self.ktree)


    def addEdge(self,edge):
        self.setParent(edge[0],edge[1])
       
    def setSize(self,newSize):
        if newSize<self.size:
            raise Exception("The size cannot be decreased.")
        if self.size==0 and newSize>0:
            self.giantSize=1
        self.size=newSize
        for index in range( len(self.ktree),newSize):
           self.ktree.append(index)
           self.subTreeWeight.append(1)
           self.suscSum+=1
           self.sizeDistribution[1]=self.sizeDistribution.get(1,0)+1

    #this is for legacy support
    def setParent(self,node,newParent):
        self.mergeSets(node,newParent)
    def getParent(self,node):
        return self.getSetIndex(node)

class KtreeInteger_old:
    def __init__(self,size=0):
        self.ktree=[]
        self.mappingOn=False
        if size!=0:
            for index in range(0,size+1):
                self.ktree.append(index);

        
    def __getRealParent(self,node):
        """
        Private method. Reads elements directly from the tree.        
        """
        try:
            return self.ktree[node]
        except IndexError:
            self.setSize(node)
            return node

    def __setRealParent(self,node,newParent):
        """
        Private.
        """
        try:
            self.ktree[node]=newParent
        except IndexError:
            self.setSize(node)
            self.ktree[node]=newParent

    def getParent(self,node):
        parent=self.__getRealParent(node)
        if node!=parent:
            self.__setRealParent(node,self.getParent(parent))
        return self.__getRealParent(node)
            
    def setParent(self,node,newParent):
        self.__setRealParent(self.getParent(node),self.getParent(newParent))

    def __iter__(self):
        for i in self.ktree:
            yield i

    def getCommStruct(self,separateElements=True):
        communityMap={}
        if self.mappingOn:
            nodes=self.ktree
        else:
            nodes=range(0,len(self.ktree))
        
        for node in nodes:
            communityKey=self.getParent(node)
            if separateElements or communityKey!=node:
                if communityKey not in communityMap:
                    communityMap[communityKey]=[node]
                else:
                    communityMap[communityKey].append(node)

        return communities.NodeCover(communityMap)

    def __len__(self):
        return len(self.ktree)

    def mergeSetsWithElements(self,elements):
        first=elements[0]
        for i in range(1,len(elements)):
            self.setParent(elements[i],first)

    def addEdge(self,edge):
        self.setParent(edge[0],edge[1])
       
    def setSize(self,newSize):
        for index in range( len(self.ktree),newSize+1):
           self.ktree.append(index);



#class KtreeMapping(KtreeInteger):
class Ktree: 
    """
    A Kruskal tree with mapping frontend. This means that node names can be any
    hashable objects.
    """
    def __init__(self,size=0,nodeNames=None):
        self.nodeIndex=netext.Enumerator()
        self.mappingOn=True

        if nodeNames!=None:
            if size<len(nodeNames):
                size=len(nodeNames)
            for nodeName in nodeNames:
                self.nodeIndex[nodeName]

        self.ktree=KtreeInteger(size)


    def getParent(self,node):
        if self.ktree.size<(self.nodeIndex[node]+1):
            self.ktree.setSize(self.nodeIndex[node]+1)
        return self.nodeIndex.getReverse(self.ktree.getParent(self.nodeIndex[node]))
    def setParent(self,node,newParent):
        requiredSize=max(self.nodeIndex[node],self.nodeIndex[newParent])+1        
        if self.ktree.size<requiredSize:
            self.ktree.setSize(requiredSize)
        self.ktree.setParent(self.nodeIndex[node],self.nodeIndex[newParent])

    def addEdge(self,edge):
        self.setParent(edge[0],edge[1])

    def __iter__(self):
        return self.nodeIndex.__iter__()

    def getCommStruct(self):
        cs=self.ktree.getCommStruct()
        newcs=communities.NodeCover()
        for c in cs:
            newc=[]
            for node in c:
                newc.append(self.nodeIndex.getReverse(node))
            newcs._addCommunity(newc)
            #newcs.comm.append(newc)
        return newcs

class Percolator:
    def __init__(self,edgesAndEvaluations,buildNet=True,symmetricNet=True,nodes=None,returnKtree=False):
        self.edges=edgesAndEvaluations
        self.buildNet=buildNet
        self.symmetricNet=symmetricNet
        self.nodes=nodes
        self.returnKtree=returnKtree

    def __iter__(self):
        if self.buildNet:
            if self.symmetricNet:
                net=pynet.SymmNet()
            else:
                net=pynet.Net()

        if self.nodes==None:
            ktree=Ktree()
        else:
            ktree=Ktree(size=len(self.nodes),nodeNames=self.nodes)
            
        for edge in self.edges:
            if isinstance(edge,EvaluationEvent):
                if self.returnKtree:
                    ktree.threshold=edge.threshold
                    ktree.addedEdges=edge.addedElements
                    yield ktree
                else:
                    cs=ktree.getCommStruct()
                    cs.threshold=edge.threshold
                    cs.addedEdges=edge.addedElements
                    if self.buildNet:
                        cs.net=net
                    yield cs
            else:
                ktree.addEdge(edge)
                if self.buildNet:
                    net[edge[0],edge[1]]=edge[2]


class EvaluationList:
    """
    EvaluationList object is an iterable object that iterates through a list returning
    EvaluationEvent objects according to given rules.
    Todo: better behavior when stacking EvaluationList objects.

    """
    def __init__(self,thelist,evaluationType="none",evaluations=None,weightFunction=lambda x:x[2],listlen=None):
        """
        Parameters
        ----------
        thelist : iterable object
            The list which is iterated when this EvaluationList is iterated. Possible EvaluationEvents
            are added between elements in thelist.
        evaluationType : string
            The method for adding evaluation events. See section 'Evaluation Types' for more information.
        evaluation : None or list of numbers
            The evalation point data. Meaning of this argument depends on the evaluationType parameter.
        weightFunction : function
            If evaluationType uses weights, this parameter defines how the weight is calculated from the 
            elements in the list.
        listlen : integer
            Lenght of the list given as thelist argument. This is useful if thelist does not have __len__
            function but the length is still known.


        Evaluation Types
        ----------------
        The evaluationType parameter defines how the data given as parameter 'evaluations'
        should be used.

        'none' : No evaluations.

        'fraclist' : Evaluations is a list of fractional number of elements to be added before 
                     EvaluationEvent. The the total number of elements must be known either
                     given in listlen parameter or by len function in thelist parameter.

        'abslist' : Evaluations is a list of absolute numbers of elements.

        'fraclinear' : Evaluation events are generated to be linearly. Evaluations is a triplet 
                       (tuple or a list), where the first element is the first evaluation, the
                       second element is the last evaluation and the third element is the total
                       number of evaluations. The first and last element is given as fraction of
                       elements in comparison to the total number of elements in the thelist parameter.
                       The the total number of elements must be known either
                       given in listlen parameter or by len function in thelist parameter.                       

        'abslinear' : Evaluation events are generated to be linearly. Evaluations is a triplet 
                      (tuple or a list), where the first element is the first evaluation, the
                      second element is the last evaluation and the third element is the total
                      number of evaluations. The first and last element is given as absolute number
                      of elements.

        'weights' : Evaluations parameter is not used. Instead EvaluationEvents are produced
                    every time when weight of an next element is different than weight of the
                    previous element.
        """
        
        self.thelist=thelist
        self.weightFunction=weightFunction
        self.strengthEvaluations=False
        self.evaluationPoints=[]
        self.lastEvaluation=False
        self.listlen=listlen

        if evaluationType=="fraclist":
            nElements=len(self)
            flist=map(lambda x:int(self.listlen*x),evaluations)
            self.setEvaluations(flist)
        elif evaluationType=="abslist":
            self.setEvaluations(evaluations)
        elif evaluationType=="weights":
            self.setStrengthEvaluations()
        elif evaluationType=="abslinear":
            self.setLinearEvaluations(evaluations[0],evaluations[1],evaluations[2])
        elif evaluationType=="fraclinear":
            nElements=len(self)
            self.setLinearEvaluations(int(evaluations[0]*nElements),int(evaluations[1]*nElements),evaluations[2])
        elif evaluationType=="none":
            pass
        else:
            raise Exception("Invalid evaluationType: "+str(evaluationType))



    def __len__(self):
        if self.listlen==None:
            self.listlen=len(self.thelist)
        return self.listlen
        
    def setEvaluations(self,evaluationPoints):
        self.evaluationPoints=sorted(evaluationPoints)
    def setLinearEvaluations(self,first,last,numberOfEvaluationPoints):
        self.strenghtEvaluations=False
        if last<=first:
            raise Exception("last<=first")
        if numberOfEvaluationPoints<2:
            raise Exception("Need 2 or more evaluation points")
        last=last-1
        self.setEvaluations(map(lambda x:int(first+(last-first)*x/float(numberOfEvaluationPoints-1)),range(0,numberOfEvaluationPoints)))
        self.lastEvaluation=False
    def setStrengthEvaluations(self):
        self.strengthEvaluations=True
        self.lastEvaluation=False
    def setLastEvaluation(self):
        self.lastEvaluation=True
    def __iter__(self):
        if not self.strengthEvaluations and not self.lastEvaluation:
            index=0
            evalIter=self.evaluationPoints.__iter__()
            try:
                nextEvaluationPoint=evalIter.next()
            except StopIteration: #thelist is empty
                nextEvaluationPoint=None
            for element in self.thelist:
                yield element
                if index==nextEvaluationPoint:
                    yield EvaluationEvent(self.weightFunction(element),index+1)
                    nextEvaluationPoint=evalIter.next()
                index+=1
        elif not self.lastEvaluation:
            last=None
            numberOfElements=0
            for element in self.thelist:
                numberOfElements+=1
                #print str(element)+"          \r", #debug purposes only
                if last!=self.weightFunction(element) and last!=None:
                    yield EvaluationEvent(last,numberOfElements-1)
                last=self.weightFunction(element)
                yield element
            yield EvaluationEvent(last,numberOfElements)

        else:
            for element in self.thelist:
                yield element
            yield EvaluationEvent()
        
def getComponents(net):
    """Get connected components of a network
    
    Parameters
    ----------
    net : A network object
    
    Returns
    -------
    Connected components as a NodePartition object.

    """
    edges=net.edges
    ee=EvaluationList(edges)
    ee.setLastEvaluation()
    p=Percolator(ee,buildNet=False,nodes=net)
    for cs in p:
        return cs

def getKCliqueComponents(net,k):
    """
    Returns community structure calculated with unweighted k-clique percolation.

    Parameters
    ----------
    net : A network object
    k : integer (larger than 2)
        The clique size.

    Returns
    -------
    The community structure as a NodeCover object.

    Examples
    --------
    >>> n=netio.loadNet('nets/co-authorship_graph_cond-mat_small.edg')
    >>> getKCliqueComponents(n,5).getSizeDist()=={5: 36, 6: 9, 7: 4, 8: 2}
    True
    """
    def evaluateAtEnd(edges):
        for edge in edges:
            yield edge
        yield EvaluationEvent()
    edgesAndEvaluations=evaluateAtEnd(net.edges)

    kcliques=kcliquesByEdges(edgesAndEvaluations,k) #unweighted clique percolation        
    for community in communitiesByKCliques(kcliques):
        return community


class KClique(object):
    """
    A class for presenting cliques of size k. Realizations
    of this class just hold a sorted list of nodes in the clique.
    """
    def __init__(self,nodelist,notSorted=True):
        self.nodes=nodelist
        if notSorted:
            self.nodes.sort()
        self.hash=None
    def __hash__(self):
        #this function is very important for speed
        if self.hash==None:
            #self.hash=hash(reduce(mul,self.nodes))
            self.hash=hash(reduce(mul,map(self.nodes[0].__class__.__hash__,self.nodes)))
            #self.hash=str.__hash__(reduce(lambda x,y:x+y,map(str,self.nodes)))
            #self.hash=int.__hash__(sum(map(self.nodes[0].__class__.__hash__,self.nodes)))
        return self.hash
    def __iter__(self):
        for node in self.nodes:
            yield node
    def __cmp__(self,kclique):
        if self.nodes==kclique.nodes:
            return 0
        else:
            return 1
    def __add__(self,kclique):
        return KClique(self.nodes+kclique.nodes)
    def getSubcliques(self):
        for i in range(0,len(self.nodes)):
            yield KClique(self.nodes[:i]+self.nodes[(i+1):],notSorted=False)
    def __str__(self):
        return str(self.nodes)
    def iterEdges(self):
	for node in self.nodes:
	    for othernode in self.nodes:
		if node!= othernode:
		   yield (node,othernode)
    getEdges = iterEdges
    def getK(self):
	return len(self.nodes)

def getIntensity(kclique,net):
    intensity=1
    for edge in kclique.getEdges():
	intensity*=net[edge[0],edge[1]]
    return pow(intensity,1.0/float(kclique.getK()*((kclique.getK()-1)/2)))

class EvaluationEvent:
    def __init__(self,threshold=None,addedElements=None):
        self.threshold=threshold
        self.addedElements=addedElements

def kcliquesAtSubnet(nodes,net,k):
    """List all k-cliques in a subnet of `net` induced by `nodes`.

    Any implementation is fine, but as this routine is a part of a
    clique percolator anyway we will use itself to find cliques larger
    than 2. Cliques of size 1 and 2 are trivial.

    Parameters
    ----------
    nodes : iterable
        Nodes that induces the subnet.
    net : pynet.SymmNet object
        The full networks where we look for the subnet induced by
        `nodes`.
    k : int, >= 1
        The number of nodes in the cliques.

    Yield
    -----
    kclique : KClique object
        The k-cliques found in the subnet of `net` induced by `nodes`.
    """
    if len(nodes)>=k:
        if k==1:
            for node in nodes:
                yield KClique([node])
        elif k==2:
            subnet=netext.getSubnet(net,nodes)
            for edge in subnet.edges:
                yield KClique([edge[0],edge[1]])
        else:
	    subnet=netext.getSubnet(net,nodes)
	    for kclique in kcliquesByEdges(subnet.edges,k):
		yield kclique

def kcliquesByEdges(edges, k):
    """Generate k-cliques from edges.

    Generator function that generates a list of cliques of size k in
    the order they are formed when edges are added in the order
    defined by the `edges` argument.  If many cliques are formed by
    adding one edge, the order of the cliques is arbitrary.
    
    This generator will pass through any EvaluationEvent objects that
    are passed to it in the `edges`.

    Parameters
    ----------
    edges : iterable with elements (node_1, node_2, weight)
        The edges that form the network. `edges` may also contain
        EvaluationEvent objects, which are simply passed through.
    k : int, >= 3
        The function returns `k`-cliques, that is, induced full subnets
        of `k` nodes.

    Yield
    -----
    kclique : KClique object
        When a new k-clique is formed, it is returned as a KClique
        object.

    Notes
    -----
    If an edge is included in `edges` multiple times, all k-cliques in
    the network constructed so far will be returned every time. Most
    of the time this is not what is wanted, so take care not the
    supply multiple edges. (LK 31.7.2009)
    """
    newNet=pynet.SymmNet() # Edges are added to an empty network one by one
    for edge in edges:
        if isinstance(edge, EvaluationEvent):
            yield edge
        else:
            # First we find all new triangles that are born when the new edge is added
            triangleEnds=set() # We keep track of the tip nodes of the new triangles
            for adjacentNode in newNet[edge[0]]: # Neighbor of one node of the edge ...
                if newNet[adjacentNode,edge[1]] != 0: #...is a neighbor of the other node of the edge...
                    triangleEnds.add(adjacentNode) #...then the neighbor is a tip of a new triangle

            # New k-cliques are now (k-2)-cliques at the triangle end points plus
            # the two nodes at the tips of the edge we are adding to the network
            for kclique in kcliquesAtSubnet(triangleEnds,newNet,k-2):
                yield kclique+KClique([edge[0],edge[1]])

            # Finally we add the new edge to the network.
            newNet[edge[0],edge[1]] = edge[2] 

def kcliquesWeight(net,k,weightFunction):
    kcliques=list(kcliquesByEdges(net.edges,k))
    kcliques.sort(lambda x,y: cmp(weightFunction(x,net),weightFunction(y,net)))
    return kcliques
    #for kclique in kcliques:
    #    yield kclique

def communitiesByKCliques(kcliques):
    # Calculate the neighboring relations
    krTree=Ktree()
    for kcliqueNumber,kclique in enumerate(kcliques):
        if isinstance(kclique,EvaluationEvent):
            communityStructure=krTree.getCommStruct().getCollapsed()
            communityStructure.threshold=kclique.threshold
            communityStructure.numberOfEdges=kclique.addedElements
            communityStructure.numberOfKCliques=kcliqueNumber+1
            yield communityStructure
        else:
            #for fewer operations at ktree, names of new cliques should be saved
            #and given to ktree when merging the sets
            krcliques=list(kclique.getSubcliques()) #list all k-1 cliques that are subcliques
            krTree.mergeSetsWithElements(krcliques) #merge the sets of k-1 cliques at the list 


def cliquePercolator(net,k,evaluations,weightFunction="minweight",evaluationType="fraclist",reverse=False):    
    """
    Weighted and unweighted clique percolation with sequential clique percolation algorithm.    

    Parameters
    ----------
    net : Network object
    k : Size of the clique to be used.
    evaluations : The evaluation data.
    weightFunction : The weight function for the cliques.
    evaluationType : How to use the evaluation data.
    reverse : True if the order of the percolation process should be reversed. That is, if small weights
              are added before large weights.

    Weight Functions
    ----------------
    Weight function which should be used to sort the cliques.

    'minweight' : Do sequential clique percolation by sorting the edges according to their weigths.

    'intensity' : Sort the cliques according to their intensities.

    Evaluation Types
    ----------------
    The evaluationType parameter defines how the data given as parameter 'evaluations'
    should be used.

    'fraclist' : Evaluations is a list of fractions of cliques (or edges for unweighted CP) 
                 from the total number of cliques (edges).

    'abslist' : Evaluations is a list of absolute numbers of cliques (or edges for unweighted CP).

    'fraclinear' : Evaluation events are generated to be linearly. Evaluations is a triplet 
                   (tuple or a list), where the first element is the first evaluation, the
                   second element is the last evaluation and the third element is the total
                   number of evaluations. The first and last element is given as fraction of
                   cliques (edges) in comparison to the total number of cliques (edges) in the thelist parameter.

    'abslinear' : Evaluation events are generated to be linearly. Evaluations is a triplet 
                  (tuple or a list), where the first element is the first evaluation, the
                  second element is the last evaluation and the third element is the total
                  number of evaluations. The first and last element is given as absolute number
                  of cliques (edges).

    'weights' : Evaluations parameter is not used. Instead EvaluationEvents are produced
                every time when weight of an next clique (edge) is different than weight of the
                previous clique (edge).

    Returns
    -------
    An iterable object generating the sequence of clique community structures.

    Scaling
    -------
    Time : Number of k-cliques * k.
    Memory : For unweighted case the memory consumption is dominated by the number of (k-1)-cliques
             in the network. For weighted case also all the k-cliques of the network are needed to be 
             kept in the memory.
    """

    #TODO: add sanity checks for the parameters.

    #Phase I: find k-cliques
    if weightFunction=="minweight":
        edges=list(net.edges)
        edges.sort(lambda x, y: cmp(x[2],y[2]),reverse=reverse)
        edgesAndEvaluations=EvaluationList(edges,evaluations=evaluations,evaluationType=evaluationType)
        kcliques=kcliquesByEdges(edgesAndEvaluations,k) #unweighted clique percolation        
    elif weightFunction=="intensity":
        kcliqueList=kcliquesWeight(net,k,getIntensity)
        kcliques=EvaluationList(kcliqueList,evaluations=evaluations,evaluationType=evaluationType,weightFunction=lambda x:getIntensity(x,net))
    else:
        raise Exception("No such weight function: "+str(weightFunction))

    #Phase II: using k-cliques, find clique communities
    for community in communitiesByKCliques(kcliques):
        yield community
        

def getKCliqueBipartiteNet(net,k):
    """
    Returns a bipartite network where to partitions are k-cliques and 
    (k-1)-cliques in the net that is given as a parameter. There is a link
    between a k-clique and (k-1)-clique if the (k-1)-clique is a subclique of 
    the k-clique.
    """

    kcliques=set()
    krcliques=set()
    kbinet=pynet.SymmNet()
    for kclique in kcliquesByEdges(net.edges,k):
        kcliques.add(kclique)
        for krclique in kclique.getSubcliques():
            krcliques.add(krclique)
            kbinet[kclique,krclique]=1
    return kbinet,kcliques,krcliques

def getKCliqueNet(net,k):
    """
    Returns a network of k-cliques in the network given as a parameter.
    Two k-cliques are adjacent if they share a (k-1)-clique.
    """
    kbinet,kcliques,krcliques=getKCliqueBipartiteNet(net,k)
    return transforms.collapseBipartiteNet(kbinet,krcliques)

def getKRCliqueNet(net,k):
    """
    Returns a network of (k-1)-cliques, which are subcliques of some k-clique in the
    network given as a parameter.
    Two (k-1)-cliques are adjacent if they are subcliques of the same k-clique.
    """
    krnet=pynet.SymmNet()
    for kclique in kcliquesByEdges(net.edges,k):
        krcliques=list(kclique.getSubcliques())
        for krclique in krcliques:
            for krclique2 in krcliques:
                krnet[krclique][krclique2]=1
    return krnet
                
