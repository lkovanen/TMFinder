import numpy
import scipy.sparse
import itertools


class VirtualNet(object):
	"""
	Virtual network class. Node names are sparse in the
	interface, but dense in the backend.
	Backend must implement following functions:
	...
	"""
	def __init__(self,sizeLimit=0):
		self._removedNodes=[]
		self._nodes={}
		self._indexToName=[]
		self.sizeLimit=sizeLimit
	def __getitem__(self,args):
		if isinstance(args, tuple):
			if len(args) != 2:
				raise KeyError, "Don't be silly. One or two indices"
			try:
				n1=self._nodes[args[0]]
			except KeyError:
				return 0.0
				#raise KeyError, "No such node: "+str(args[0])
			try:
				n2=self._nodes[args[1]]
			except KeyError:
				return 0.0
				#raise KeyError, "No such node: "+str(args[1])
			if n1==n2:
				return 0
			#assert n1!=n2, "No self-edges."
			return self._getEdge(n1,n2)
		else:
			return Node(self,args)
	def __setitem__(self,key,val):
		if isinstance(key, tuple):
			if len(key) != 2:
				raise KeyError, "Don't be silly. Give two indices"
		if key[0] not in self:
			self.addNode(key[0])
		if key[1] not in self:
			self.addNode(key[1])
		assert key[0]!=key[1], "No self-edges."		
		self._setEdge(self._nodes[key[0]],self._nodes[key[1]],val)
		return val

	def __delitem__(self,args):
		if isinstance(args, tuple):
			if len(args) != 2:
				raise KeyError, "Don't be silly. One or two indices"
			self[args[0], args[1]]=0
		else:
			self.delNode(args)

	def __contains__(self,nodeName):
		return nodeName in self._nodes
	def __iter__(self):
		return self._nodes.iterkeys()
	def __len__(self):
		return len(self._nodes)

	def __copy__(self):
		copyNet=self.__class__(self.sizeLimit)

		#add nodes to the new net
		for node in self:
			copyNet.addNode(node)
			
		#add the edges to the backend
		if self.symmetric:
			for node1 in self:			
				for node2 in self[node1]:
					if node1.__hash__()<node2.__hash__():
						copyNet[node1,node2]=self[node1,node2]
		else:
			for node1 in self:			
				for node2 in self[node1]:
					copyNet[node1,node2]=self[node1,node2]			

		return copyNet

	def addNode(self,nodeName):
		"""
		Adds an empty node to the network.
		"""
		if nodeName not in self:
			if len(self._removedNodes)>0:
				newIndex=self._removedNodes.pop()
				self._indexToName[newIndex]=nodeName
			else:
				#assert len(self._nodes)==len(self._indexToName)
				newIndex=len(self._nodes)
				#check for size limit
				if self.sizeLimit!=0 and newIndex>self.sizeLimit:
					raise Exception,"Size limit for the network reached." 
				self._indexToName.append(nodeName)
				self._addNode()
			self._nodes[nodeName]=newIndex

	def delNode(self,nodeName): #override for directed
		"""
		Removes a node from the network.	     
		"""
		#Check that node is in the net
		if nodeName in self:
		        #First, remove all the edges
			self._removeNodeEdges(nodeName)
			
			#then remove the empty node
			removedIndex=self._nodes[nodeName]
			self._removedNodes.append(removedIndex)
			del self._nodes[nodeName]

		
	       
	def isSymmetric(self):
		return self.symmetric

	def deg(self,nodeName):
		return self._degIndex(self._nodes[nodeName])

	def inDeg(self,nodeName):
		return self._inDegIndex(self._nodes[nodeName])

	def outDeg(self,nodeName):
		return self._outDegIndex(self._nodes[nodeName])

	#--- Symmetric net:
	def _inDegIndex(self,nodeIndex):
		return self._degIndex(nodeIndex)
	def _outDegIndex(self,nodeIndex):
		return self._degIndex(nodeIndex)
	def _iterNodeIn(self,nodeIndex):
		return self._iterNode(nodeIndex)
	def _iterNodeOut(self,nodeIndex):
		return self._iterNode(nodeIndex)
	def _removeNodeEdges(self,nodeName):
		for neigh in list(self[nodeName]):
			self[nodeName,neigh]=0
	symmetric=True

	#--- Virtual functions to override
	def _addNode(self):
		pass
	def _degIndex(self,nodeIndex):
		raise NotImplemented
	def _getEdge(self,src,dest):
		raise NotImplemented
	def _setEdge(self,src,dest,val):
		raise NotImplemented
	def _iterNode(self,nodeIndex):
		raise NotImplemented
	

class Node(object):
	def __init__(self,net,name):
		self.net=net
		self.name=name
	def _indexToNameIter(self,iterator):
		for index in iterator:
			yield self.net._indexToName[index]
	def __iter__(self):
		if not self.name in self.net:
			return [].__iter__()
			#raise KeyError, "Node not in the network."
		#name -> index -> backend iterator -> index -> name
		return self._indexToNameIter(self.net._iterNode(self.net._nodes[self.name]))

	def __getitem__(self, name):
		return self.net[self.name, name]

	def __setitem__(self, name, val):
		 self.net[self.name, name]=val
		 return val

	def __contains__(self,nodeName):
		return self[nodeName]!=0

	def deg(self):
		return self.net.deg(self.name)

	def inDeg(self):
		return self.net.inDeg(self.name)

	def outDeg(self):
		return self.net.outDeg(self.name)


	def iterIn(self):
		if not self.name in self.net:
			raise KeyError, "Node not in the network."
		#name -> index -> backend iterator -> index -> name
		return self._indexToNameIter(self.net._iterNodeIn(self.net._nodes[self.name]))
	def iterOut(self):
		if not self.name in self.net:
			raise KeyError, "Node not in the network."
		#name -> index -> backend iterator -> index -> name
		return self._indexToNameIter(self.net._iterNodeOut(self.net._nodes[self.name]))
	

class VirtualDirNet(VirtualNet):
	def _removeNodeEdges(self,nodeName):
		for neigh in list(self[nodeName]):
			self[nodeName,neigh]=0
			self[neigh,nodeName]=0

	symmetric=False

	#--- Virtual functions to override
	def _iterNodeIn(self,nodeIndex):
		raise NotImplemented
	def _iterNodeOut(self,nodeIndex):
		raise NotImplemented
	def _inDegIndex(self,nodeIndex):
		raise NotImplemented
	def _outDegIndex(self,nodeIndex):
		raise NotImplemented

class ScipySparseSymmNet(VirtualNet):
	def __init__(self,sizeLimit=0):
		VirtualNet.__init__(self,sizeLimit=sizeLimit)
		self._nodeList=[]
	def _addNode(self):
		newNode=scipy.sparse.dok_matrix((10**10,1),dtype='float')
		self._nodeList.append(newNode)

	#--- Methods used by VirtualNet:
	def _degIndex(self,nodeIndex):
		return len(self._nodeList[nodeIndex])
	def _getEdge(self,src,dest):
		return self._nodeList[src][dest,0]
	def _setEdge(self,src,dest,val):
		nNodes=len(self._nodeList)
		while nNodes<=src or nNodes<=dest:
			self._addNode()
			nNodes+=1
		if val==0:
			try:
				del self._nodeList[src][dest,0]
			except KeyError:
				pass
			try:
				del self._nodeList[dest][src,0]
			except KeyError:
				pass
		else:
			self._nodeList[src][dest,0]=val
			self._nodeList[dest][src,0]=val
	def _iterNode(self,nodeIndex):
		return itertools.imap(lambda x:x[0],self._nodeList[nodeIndex].iterkeys())


class ScipySparseDirNet(VirtualDirNet):
	def __init__(self,sizeLimit=0):
		VirtualNet.__init__(self,sizeLimit=sizeLimit)
		self._nodeList=[] #directed edges and weights
		self._backNodeList=[] #reversed directed edges 
		self._totalDeg=[] #total degree of nodes
	def _addNode(self):
		newNode=scipy.sparse.dok_matrix((10**10,1),dtype='float')
		self._nodeList.append(newNode)
		newBackNode=scipy.sparse.dok_matrix((10**10,1),dtype='bool')
		self._backNodeList.append(newBackNode)
		self._totalDeg.append(0)

	#--- Methods used by VirtualDirNet:
	def _degIndex(self,nodeIndex):
		return self._totalDeg[nodeIndex]
	def _getEdge(self,src,dest):
		return self._nodeList[src][dest,0]
	def _setEdge(self,src,dest,val):
		nNodes=len(self._nodeList)
		if src>=nNodes: 
			self._addNode()
		if dest>=nNodes: 
			self._addNode()
		if val==0: #removing an edge
			#nothing happens if the link does not exist:
			if (dest,0) in self._nodeList[src]:
				if not (src,0) in self._nodeList[dest]:
					self._totalDeg[src]+=-1
					self._totalDeg[dest]+=-1
				self._nodeList[src][dest,0]=0
				self._backNodeList[dest][src,0]=0
		else:
			if not (dest,0) in self._nodeList[src] and not (src,0) in self._nodeList[dest]:
				self._totalDeg[src]+=1
				self._totalDeg[dest]+=1

			self._nodeList[src][dest,0]=val
			self._backNodeList[dest][src,0]=1


	def _iterNode(self,nodeIndex):
		#First iter through all outgoing neighbors:
		for neigh in self._nodeList[nodeIndex].iterkeys():
			yield neigh[0]
		#Then iter through only incoming neighbors
		for neigh in self._backNodeList[nodeIndex].iterkeys():
			if not neigh in self._nodeList[nodeIndex]:
				yield neigh[0]
	def _iterNodeIn(self,nodeIndex):
		return itertools.imap(lambda x:x[0],self._backNodeList[nodeIndex].iterkeys())
	def _iterNodeOut(self,nodeIndex):
		return itertools.imap(lambda x:x[0],self._nodeList[nodeIndex].iterkeys())

	def _inDegIndex(self,nodeIndex):
		return len(self._backNodeList[nodeIndex])
	def _outDegIndex(self,nodeIndex):
		return len(self._nodeList[nodeIndex])

class NumpyFullSymmNet(VirtualNet):
	def __init__(self,sizeLimit):
		VirtualNet.__init__(self,sizeLimit=sizeLimit)
		self._adjMatrix=numpy.zeros([sizeLimit,sizeLimit])
		self._degree=numpy.zeros(self.sizeLimit,dtype='uint')

	#--- Methods used by VirtualNet:
	def _degIndex(self,nodeIndex):
		return int(self._degree[nodeIndex])
	def _getEdge(self,src,dest):
		return self._adjMatrix[src,dest]
	def _setEdge(self,src,dest,val):
		if val==0 and self._adjMatrix[src,dest]!=0:
			self._degree[src]+=-1
			self._degree[dest]+=-1
		elif val!=0 and self._adjMatrix[src,dest]==0:
			self._degree[src]+=1
			self._degree[dest]+=1
		self._adjMatrix[src][dest]=val
		self._adjMatrix[dest][src]=val
	def _iterNode(self,nodeIndex):
		for i in range(0,self.sizeLimit):
                        if self._adjMatrix[nodeIndex,i]!=0:
                                yield i

class NumpyFullDirNet(VirtualDirNet):
	def __init__(self,sizeLimit):
		VirtualNet.__init__(self,sizeLimit=sizeLimit)
		self._adjMatrix=numpy.zeros([sizeLimit,sizeLimit])
		self._degree=numpy.zeros(self.sizeLimit)
		self._inDegree=numpy.zeros(self.sizeLimit)
		self._outDegree=numpy.zeros(self.sizeLimit)

	#--- Methods used by VirtualDirNet:
	def _degIndex(self,nodeIndex):
		return int(self._degree[nodeIndex])
	def _getEdge(self,src,dest):
		return self._adjMatrix[src,dest]
	def _setEdge(self,src,dest,val):
		if val==0:
			if self._adjMatrix[src,dest]!=0:
				self._outDegree[src]+=-1
				self._inDegree[dest]+=-1
				if self._adjMatrix[dest,src]==0:
					self._degree[src]+=-1
					self._degree[dest]+=-1
		elif val!=0:
			if self._adjMatrix[src,dest]==0:
				self._outDegree[src]+=1
				self._inDegree[dest]+=1
				if self._adjMatrix[dest,src]==0:
					self._degree[src]+=1
					self._degree[dest]+=1
		self._adjMatrix[src][dest]=val
	def _iterNode(self,nodeIndex):
		for i in range(0,self.sizeLimit):
                        if self._adjMatrix[nodeIndex,i]!=0 or self._adjMatrix[i,nodeIndex]!=0:
                                yield i
	def _iterNodeIn(self,nodeIndex):
		for i in range(0,self.sizeLimit):
                        if self._adjMatrix[i,nodeIndex]!=0:
                                yield i
	def _iterNodeOut(self,nodeIndex):
		for i in range(0,self.sizeLimit):
                        if self._adjMatrix[nodeIndex,i]!=0:
                                yield i
	def _inDegIndex(self,nodeIndex):
		return int(self._inDegree[nodeIndex])
	def _outDegIndex(self,nodeIndex):
		return int(self._outDegree[nodeIndex])



class LCELibSparseSymmNet(VirtualNet):
	def __init__(self,sizeLimit=0):
		VirtualNet.__init__(self,sizeLimit=sizeLimit)
		self._net=_cnet.new_Sn(0)
	def __del__(self):
		_cnet.delete_Sn(self._net)
		

	#--- Methods used by VirtualNet:
	def _degIndex(self,nodeIndex):
		return _cnet.Sn_getDegree(self._net,nodeIndex)
	def _getEdge(self,src,dest):
		return _cnet.Sn_getEdge(self._net,src,dest)
	def _setEdge(self,src,dest,val):
		_cnet.Sn_setEdge(self._net,src,dest, float(val))
	def _iterNode(self,nodeIndex):
		citerator=_cnet.Sn_getNeighborIterator(self._net,nodeIndex)
		next=_cnet.NeighborIterator_getNext(citerator)
		while next!=-1:
			yield next
			next=_cnet.NeighborIterator_getNext(citerator)
		_cnet.delete_NeighborIterator(citerator)

class LCELibSparseDirNet(VirtualDirNet):
	def __init__(self,sizeLimit=0):
		VirtualNet.__init__(self,sizeLimit=sizeLimit)
		self._net=_cnet.new_Dn(0)
	def __del__(self):
		_cnet.delete_Dn(self._net)
		

	#--- Methods used by VirtualNet:
	def _degIndex(self,nodeIndex):
		return _cnet.Dn_getDegree(self._net,nodeIndex)
	def _getEdge(self,src,dest):
		return _cnet.Dn_getEdge(self._net,src,dest)
	def _setEdge(self,src,dest,val):
		_cnet.Dn_setEdge(self._net,src,dest, float(val))
	def _iterNode(self,nodeIndex):
		citerator=_cnet.Dn_getNeighborIteratorAll(self._net,nodeIndex)
		next=_cnet.NeighborIteratorAll_getNext(citerator)
		while next!=-1:
			yield next
			next=_cnet.NeighborIteratorAll_getNext(citerator)
		_cnet.delete_NeighborIteratorAll(citerator)

	#--- Methods used by VirtualDirNet:
	def _iterNodeIn(self,nodeIndex):
		citerator=_cnet.Dn_getNeighborIteratorIn(self._net,nodeIndex)
		next=_cnet.NeighborIteratorIn_getNext(citerator)
		while next!=-1:
			yield next
			next=_cnet.NeighborIteratorIn_getNext(citerator)
		_cnet.delete_NeighborIteratorIn(citerator)
	def _iterNodeOut(self,nodeIndex):
		citerator=_cnet.Dn_getNeighborIteratorOut(self._net,nodeIndex)
		next=_cnet.NeighborIteratorOut_getNext(citerator)
		while next!=-1:
			yield next
			next=_cnet.NeighborIteratorOut_getNext(citerator)
		_cnet.delete_NeighborIteratorOut(citerator)
	def _inDegIndex(self,nodeIndex):
		return _cnet.Dn_getInDegree(self._net,nodeIndex)
	def _outDegIndex(self,nodeIndex):
		return _cnet.Dn_getOutDegree(self._net,nodeIndex)




#--- Implementation lists
SymmBackends=[LCELibSparseSymmNet,ScipySparseSymmNet,NumpyFullSymmNet]
DirBackends=[LCELibSparseDirNet,ScipySparseDirNet,NumpyFullDirNet]

#--- Default implementations
DirNet=ScipySparseDirNet #LCELibSparseDirNet overrides this
Net=DirNet
SymmNet=ScipySparseSymmNet #LCELibSparseSymmNet overrides this
SymmFullNet=NumpyFullSymmNet
FullNet=NumpyFullDirNet


#--- Try to import the C++-implementation.
#try:
#        from cnet import _cnet
#	SymmNet=LCELibSparseSymmNet
#	DirNet=LCELibSparseDirNet
#	Net=DirNet
#	SymmBackends.append(LCELibSparseSymmNet)
#except ImportError:
#        print "Importing the LCELib network data structure failed."





if __name__ == '__main__':
    """Run unit tests if called."""
    from tests.test_pynet import *
    test_pynet()
