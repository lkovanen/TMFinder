"""Input and output functions for pynet
  
  Current status
  --------------
  This module contains functions for reading and writing files
  containing networks.  loadNet and writeNet are general functions for
  loading and writing networks in different file types. This means
  that they try to recognize the filetype of the network file from the
  filename and then use appropriate function to read/write it.

  The user can also use loadNet_[format] and writeNet_[format]
  directly to force loading or writing to be done in given format.

  Currently only the edg, gml and matrix format has been implemented 

  Future additions
  ----------------
  - Important: make loadNet_? work as loadNet, so that no that
    filename is input
  - Support for metadata in network files
  - graphXML-format and others
  - make gml-io stable
"""

import pynet,netext,warnings
import sys
knownFiletypes=["edg","gml","mat","net","adj"]

def getFiletype(fileName):
    """Infer the type of a file.

    (Current behaviour is to just return the file name suffix after
    the last dot in fileName.)

    Parameters
    ----------
    filename : string
        The filename whose type you want to know.

    Return
    ------
    filetype : string
        A string literal depicting the file type.
    """

    # Return the file identifier after the last dot.
    # Examples: mynet.edg     ==>   edg
    #           mynet.old.mat ==>   mat
    return fileName.split('.')[-1]


def loadNet_gml(input):
    """
    Reads a networks data from input in gml format.

    Note: This is not a complete gml-parser, because gml format can
    hold almost anykind of data. Instead this parser tries to find
    edges of one graph from the given input. Use at you own risk with
    complicated gml-files.
    """

    source=None
    target=None
    value=None
    for line in input:
        line=line.strip()
        if line.startswith("directed"):
            if line[9:10]=="0":
                net=pynet.SymmNet()
            elif line[9:10]=="1":
                net=pynet.Net()
        elif line.startswith("source"):
            source=line.split()[1]
        elif line.startswith("target"):
            target=line.split()[1]
        elif line.startswith("value"):
            value=line.split()[1]
        elif line.startswith("edge"):
            if source!=None and target!=None:
                if value!=None:
                    net[source][target]=float(value)
                else:
                    net[source][target]=1
                source=None
                target=None
                value=None
    if source!=None and target!=None:
        if value!=None:
            net[source][target]=float(value)
        else:
            net[source][target]=1


    return net

def loadNet_edg(input, mutualEdges=False, splitterChar=None, symmetricNet=True,
                numerical=None, allowSelfEdges=True, hasHeaderLine=False):
    """Read network data from input in edg format.

    If `mutualEdges` is set to True, an edge is added between nodes i
    and j only if both edges (i,j) and (j,i) are listed. The weight of
    the edge is the sum of the weights of the original edges.

    If the same edge is encountered multiple times, the edge weight
    will be the sum of all weights.

    If 'allowSelfEdges', the self edges are translated as nodes with
    no edges. Otherwise those edges are just thrown out.

    If hasHeaderLine is True the first line is skipped, if it is False
    the first line is read normally, and if it is None the first line
    is skanned for "from to weight" string to decide if it is a header
    line.
    """
    def isNumerical(input):
	try:
	   for line in input:
		int(line.split(splitterChar)[0])
		int(line.split(splitterChar)[1])
	except ValueError:
	   input.seek(0)
	   return False
	input.seek(0)
	return True

    if numerical is None:
        numerical = isNumerical(input)
    
    if symmetricNet:
        newNet=pynet.SymmNet()
    else:
        newNet=pynet.Net()

    #Check for headers
    possibleHeaders=[["from","to","weight"],["head","tail","weight"]]
    if hasHeaderLine != False:
        firstLine=input.readline().strip().lower()
        fields=firstLine.split(splitterChar)
        if hasHeaderLine==None and fields not in possibleHeaders:
            input.seek(0)
        if hasHeaderLine==True:
            input.seek(0)

    nodeMap = {} # Used only if mutualEdges = True.

    for line in input:
        fields=line.split(splitterChar)
        if len(fields)>1:
            if len(fields)==2: #if weight is missing:
                fields[1]=fields[1].strip('\n') #strip the endline from the node name
                fields.append(1) #add one as the weight
            if fields[0] != fields[1] or allowSelfEdges:
                if numerical:
                    fields[0]=int(fields[0])
                    fields[1]=int(fields[1])
                if fields[0]!=fields[1]:
                    if mutualEdges:
                        if nodeMap.has_key( (fields[1], fields[0]) ):
                            newNet[fields[0]][fields[1]] += nodeMap[(fields[1], fields[0])]
                            newNet[fields[0]][fields[1]] += float(fields[2])
                            nodeMap[(fields[1], fields[0])] = 0
                            nodeMap[(fields[0], fields[1])] = 0
                        else:
                            nodeMap[(fields[0], fields[1])] = nodeMap.get((fields[0], fields[1]), 0) + float(fields[2])
                    else:
                        newNet[fields[0]][fields[1]] += float(fields[2])
                else:
                    newNet.addNode(fields[0])

    return newNet

def loadNet_adj(input, mutualEdges=False, splitterChar=None, symmetricNet=True,
                numerical=None):
    raise Exception("Reading adj file format is not implemented.")

def loadNet_mat(input, mutualEdges=False, splitterChar=None,symmetricNet=True,nodeNames=[],type="square"):
    """
    Loads a network from a file which is in a weight matrix format. Nodes are ordered in ascending order by their names.
    
    Parameters
    ----------
    type : string
        "square": Input is in a square matrix format. Weight of an edge between node i an j 
        is in row i and column j of that row. 
        "upperdiag": Same as the square matrix, but now only the elements at the diagonal and above
        are present. 
        "lowerdiag": See upperdiag
        "supperdiag": Strictly upper diagonal. The diagonal elements are not present.
        "slowerdiag": see supperdiag.

    Returns
    -------
    The network that is loaded in a FullNet or in SymmFullNet format.
    """


    usenodelist=False	
    if len(nodeNames)>0:
	usenodelist=True

    #Get the network size
    if usenodelist:
        netSize=len(nodeNames)
    else: #we have to infer the size from the file
        if type=="square" or type=="upperdiag" or type=="supperdiag":
            netSize=len(input.readline().strip().split(splitterChar))
            input.seek(0)
            if type=="supperdiag": netSize+=1
        elif type=="lowerdiag" or type=="slowerdiag":
            for line in input: pass #it would be faster to read backwards
            input.seek(0)
            netSize=len(line.strip().split(splitterChar))
            if type=="slowerdiag": netSize+=1
        else:
            raise Exception("Invalid type for the matrix: "+str(type))

    if symmetricNet:
        newNet=pynet.SymmFullNet(netSize)
    else:
        newNet=pynet.FullNet(netSize)

    if type=="square":
        for rowIndex,line in enumerate(input):
            fields=line.strip().split(splitterChar)
            
            #Check that the matrix is of right form
            if netSize!=(len(fields)):
                if usenodelabels and rowIndex==0:
                    raise Exception("The length of the node label list does not macth the number of columns.")
                else:
                    raise Exception("The length of row "+str(rowIndex)+" does not match the previous rows.")

            #Now fill the row of the matrix
            for columnIndex,element in enumerate(fields):
                if columnIndex!=rowIndex:
                    if usenodelist:
                        newNet[nodeNames[rowIndex],nodeNames[columnIndex]]=float(element)
                    else:
                        newNet[rowIndex,columnIndex]=float(element)
        if (rowIndex+1)!=netSize:
            raise Exception("Invalid number of rows: There are "+str(rowIndex+1)+" rows and "+str(netSize)+" columns.")            

    elif type=="upperdiag":
        for rowIndex,line in enumerate(input):
            fields=line.strip().split(splitterChar)

            #Check that the matrix is of right form
            if netSize!=(len(fields)+rowIndex):
                if usenodelabels and rowIndex==0:
                    raise Exception("The length of the node label list does not macth the number of columns.")
                else:
                    raise Exception("The length of row "+str(rowIndex)+" does not match the previous rows.")

            #Now fill the row of the matrix
            for columnIndex,element in enumerate(fields[1:]):
                columnIndex+=1
                if usenodelist:
                    newNet[nodeNames[rowIndex],nodeNames[columnIndex+rowIndex]]=float(element)
                else:
                    newNet[rowIndex,columnIndex+rowIndex]=float(element)
        if (rowIndex+1)!=netSize:
            raise Exception("Invalid number of rows: There are "+str(rowIndex+1)+" rows and "+str(netSize)+" columns.")

    elif type=="supperdiag":
        for rowIndex,line in enumerate(input):
            fields=line.strip().split(splitterChar)

            #Check that the matrix is of right form
            if netSize!=(len(fields)+rowIndex+1):
                if usenodelabels and rowIndex==0:
                    raise Exception("The length of the node label list does not macth the number of columns.")
                else:
                    raise Exception("The length of row "+str(rowIndex)+" does not match the previous rows.")

            #Now fill the row of the matrix
            for columnIndex,element in enumerate(fields):
                columnIndex+=1
                if usenodelist:
                    newNet[nodeNames[rowIndex],nodeNames[columnIndex+rowIndex]]=float(element)
                else:
                    newNet[rowIndex,columnIndex+rowIndex]=float(element)
        if (rowIndex+2)!=netSize:
            raise Exception("Invalid number of rows: There are "+str(rowIndex+1)+" rows and "+str(netSize)+" columns.")

    elif type=="lowerdiag":
        for rowIndex,line in enumerate(input):
            fields=line.strip().split(splitterChar)

            #Check that the matrix is of right form
            if len(fields)!=(rowIndex+1):
                raise Exception("The length of row "+str(rowIndex)+" does not match the previous rows.")

            #Now fill the row of the matrix
            for columnIndex,element in enumerate(fields[:-1]):
                if usenodelist:
                    newNet[nodeNames[rowIndex],nodeNames[columnIndex]]=float(element)
                else:
                    newNet[rowIndex,columnIndex]=float(element)
        if (rowIndex+1)!=netSize:
            raise Exception("Invalid number of rows: There are "+str(rowIndex+1)+" rows and "+str(netSize)+" columns.")

    elif type=="slowerdiag":
        for rowIndex,line in enumerate(input):
            fields=line.strip().split(splitterChar)
            rowIndex+=1

            #Check that the matrix is of right form
            if len(fields)!=rowIndex:
                raise Exception("The length of row "+str(rowIndex)+" does not match the previous rows.")

            #Now fill the row of the matrix
            for columnIndex,element in enumerate(fields):
                if usenodelist:
                    newNet[nodeNames[rowIndex],nodeNames[columnIndex]]=float(element)
                else:
                    newNet[rowIndex,columnIndex]=float(element)
        if (rowIndex+1)!=netSize:
            raise Exception("Invalid number of rows: There are "+str(rowIndex+1)+" rows and "+str(netSize)+" columns.")


    return newNet

def loadNet_net(input):
    raise Exception("Reading Pajek file format is not implemented.")

def writeNet_gml(net,outputFile):
    if not hasattr(outputFile, 'write'):
        raise ValueError("Parameter 'outputFile' must be a file object.")
    outputFile.write("graph [\n")
    if net.__class__==pynet.SymmNet:
        outputFile.write("directed 0\n")
    else:
        outputFile.write("directed 1\n")

    nodeIndex=netext.Enumerator()
    for node in net:
        outputFile.write("node [\n")
        outputFile.write("id "+str(nodeIndex[node])+"\n")
        outputFile.write("label "+str(node))
        outputFile.write("]\n")

    for edge in net.edges:
        outputFile.write("edge [\n")
        outputFile.write("source " + str(nodeIndex[edge[0]])+"\n")
        outputFile.write("target " + str(nodeIndex[edge[1]])+"\n")
        outputFile.write("value " + str(edge[2])+"\n")
        outputFile.write("]\n")

    outputFile.write("]\n")
    outputFile.close()

def writeNet_edg(net, outputFile, headers=False):
    if not hasattr(outputFile, 'write'):
        raise ValueError("Parameter 'outputFile' must be a file object.")
    #edges=netext.getEdges(net)
    edges=net.edges
    if headers==True:
        outputFile.write("HEAD\tTAIL\tWEIGHT\n")
    for edge in edges:
        outputFile.write("\t".join(map(str, edge)) + "\n")

def writeNet_net(net, outputFile):
    """
    Write network files in Pajek format.

    Todo: add writing metadata to the vertices rows
    """
    if not hasattr(outputFile, 'write'):
        raise ValueError("Parameter 'outputFile' must be a file object.")
        
    #Writing vertices to the disk.
    numberOfNodes = len(net)
    nodeNameToIndex = {}
    outputFile.write("*Vertices "+str(numberOfNodes)+"\n")
    for index,node in enumerate(net):
        outputFile.write(str(index+1)+' "'+str(node)+'"\n')
        nodeNameToIndex[node]=index+1

    #Writing edges to the disk
    outputFile.write("*Edges\n")
    for edge in net.edges:
        outputFile.write(str(nodeNameToIndex[edge[0]]) + "\t" 
                         + str(nodeNameToIndex[edge[1]]) + "\t"
                         + str(edge[2]) + "\n")

    del nodeNameToIndex

def writeNet_mat(net, outputFile,type="square"):
    """
    Writes network into file in weight matrix format. Nodes are ordered in ascending order by their names.
    
    Parameters
    ----------
    type : string
        "square": Outputs the matrix in a square format. Weight of an edge between node i an j 
        is in row i and column j of that row. 
        "upperdiag": Same as the square matrix, but now only the elements at the diagonal and above
        are writen. 
        "lowerdiag": See upperdiag
        "supperdiag": Strictly upper diagonal. The diagonal elements are not written.
        "slowerdiag": see supperdiag.

    The list of the node names in the order they are used in the weight matrix is 
    returned.
    """
    if not hasattr(outputFile, 'write'):
        raise ValueError("Parameter 'outputFile' must be a file object.")

    nodes=list(net)
    nodes.sort()
    if type=="square":
        for i in nodes:
            first=True
            for j in nodes:
                if first:first=False
                else:outputFile.write("\t")
                outputFile.write(str(net[i,j]))
            outputFile.write("\n")
    elif type=="upperdiag":
        for iindex,i in enumerate(nodes):
            first=True
            for jindex,j in enumerate(nodes):
                if first:first=False
                else:outputFile.write("\t")
                if iindex<=jindex:
                    outputFile.write(str(net[i,j]))
            outputFile.write("\n")
    elif type=="supperdiag":
        for iindex,i in enumerate(nodes):
            if iindex!=len(nodes)-1:
                for jindex,j in enumerate(nodes):
                    if jindex>1:outputFile.write("\t")
                    if iindex<jindex:
                        outputFile.write(str(net[i,j]))
                outputFile.write("\n")
    elif type=="lowerdiag":
        for iindex,i in enumerate(nodes):
            for jindex,j in enumerate(nodes):
                if jindex!=0 and jindex<=iindex: outputFile.write("\t")
                if iindex>=jindex:
                    outputFile.write(str(net[i,j]))
            outputFile.write("\n")    
    elif type=="slowerdiag":
        for iindex,i in enumerate(nodes):
            if iindex!=0:
                for jindex,j in enumerate(nodes):
                    if jindex!=0 and jindex<iindex:outputFile.write("\t")
                    if iindex>jindex:
                        outputFile.write(str(net[i,j]))
                outputFile.write("\n")
    else:
        raise ValueError("Invalid value for parameter 'type'.")
    return nodes

def writeNet_adj(net,outputFile,splitterChar=" "):
    """
    Writes the network into disk as a directed adjacency list.
    Each line in the output file corresponds to node and its neighbors. That is,
    first node in each line is connected to the consequent nodes in the same line.
    """

    if not hasattr(outputFile, 'write'):
        raise ValueError("Parameter 'outputFile' must be a file object.")
    for node in net:
        nodes=list(net[node])
        nodes.append(node)
        nodes.reverse()
        outputFile.write(splitterChar.join(map(str,nodes)))
        outputFile.write("\n")        

def writeNet(net, output, headers=False, fileType=None):
    """Write network to disk.

    Parameters
    ----------
    net : pynet network object
        The network to write.
    output : str or file
        Name of the file to be opened.
    headers : bool
        If true, print headers before the actual network data (affects
        only edg format).
    fileType : str
        Type of the output file. In None, the suffix of fileName will
        be used to guess the file type.

    Exceptions
    ----------
    ValueError : If file type is unknown or unable to write to
                 `output`.
    """
    # If `output` is a string, we assume it is a file name and open
    # it. Otherwise if it implements 'write'-method we assume it is a
    # file object.
    fileOpened = False
    if isinstance(output, str):
        outputFile = open(output, 'w')
        fileOpened = True
    elif not hasattr(output, 'write'):
        raise ValueError("'output' must be a string or an object "
                         "with a 'write'-method.")
    else:
        outputFile = output

    try:
        # Infer file type if not explicitely given.
        if fileType is None and hasattr(outputFile, 'name'):
            fileType = getFiletype(outputFile.name)

        # Write out the network.
        if fileType == 'edg':
            writeNet_edg(net, outputFile, headers)
        #elif fileType in ('gml', 'mat', 'net'):
        elif fileType in knownFiletypes:
            eval("writeNet_%s(net,outputFile)" % fileType)
        else:
            raise ValueError("Unknown file type, try writeNet_[filetype].")
    finally:
        if fileOpened:
            outputFile.close()

def loadNet(input, fileType=None, **keywords):
    """Read network from disk.

    Parameters
    ----------
    input : str or file
        Name of the file to be opened or a file object.
    fileType : str
        Type of the output file. In None, the suffix of fileName will
        be used to guess the file type.

    Exceptions
    ----------
    ValueError : If file type is unknown or unable to read from
                 `input`.
    """
    #mutualEdges=False, splitterChar=None, symmetricNet=True,
    #numerical=None, fileType=None, hasHeaderLine=False, allowSelfEdges=True):


    # If `input` is a string, we assume it is a file name and open
    # it. Otherwise if it implements 'write'-method we assume it is a
    # file object.
    fileOpened = False
    if isinstance(input, str):
        inputFile = open(input, 'r')
        fileOpened = True
    elif not isinstance(input, file):
        raise ValueError("'input' must be a string or a file object.")
    else:
        inputFile = input
    
    # Infer file type if not explicitely given.
    if fileType is None and hasattr(inputFile, 'name'):
        fileType = getFiletype(inputFile.name)

    # Read in the network.
    try:
        # edg-files need different behaviour.
        if fileType == 'edg':
            newNet = loadNet_edg(inputFile, **keywords)
        #elif fileType in ('gml', 'mat', 'net'):
        elif fileType in knownFiletypes:
            newNet = eval("loadNet_%s(inputFile)" % fileType)
        else:
            raise ValueError("Unknown file type '%s', try loadNet_[filetype]."
                             % fileType)
    finally:
        if fileOpened:
            inputFile.close()

    return newNet
    
def loadNodeProperties(net,input,splitterChar=None,propertyNames=None,
                       allowMissingData=False,allowExtraData=False):
    """Read metadata (properties for nodes) from a file.

    The metadata file can contain any number of columns. The first
    column should contain names of nodes contained in 'net', and the
    other columns contain user-defined properties.
    
    If a list 'propertyNames' is not given, the first row must contain
    headers. The first column header should be node_label, and the
    other column headers are names of the user-defined properties.
    They are automatically appended to the property list in 'net'.
    Alternatively, you can provide a list 'propertyNames' containing a
    label for each column. In this case, your file should not contain
    a header. The function 'loadNodeProperties' checks whether
    'propertyNames' contains 'node_label' as the first element, and
    adds it if it doesn't, so you do not need to give it explicitly.

    Example input file format:
    node_label node_color node_class
    node1      blue       class1

    """

    #todo: default properties as argument for missing property lines

    #tested for i) SymmNet with string node labels, ii) -"- with integer labels,
    # iii) FullNet, iv) the above cases with too many lines / non-existing nodes
    # in input metadata file
    # <--- you should have written a unit test for that, please do that next time

    def isanum(str):
        """
        Checks if a string contains only digits, decimal points "." or minus signs "-"        
        """
        from string import digits
        for c in str:
            if not c in digits and c!="." and c!="-": return False
        return True


    def isint(str):
        """
        Checks if a string contains only digits or minus signs "-"        
        """
        from string import digits
        for c in str:
            if not c in digits and c!="-": return False
        return True

    def getNumberOfLines(input,nfields):
        """
        Returns length of the given file in lines. Throws IOError if the file does not 
        exist.
        """
        if isinstance(input, str):
            theFile = open(input, 'rU')
        else:
            theFile = input
        
        i=0
        for line in theFile:
            fields=line.split(splitterChar)
            if len(fields)!=nfields:
                f.close()
                raise Exception("Invalid number of fields on row: "+str(i+1))
            i+=1
        #theFile.close()
        return i

    def checkNodes(input,net,fieldNames,hasHeader,splitterChar):
        """
        Returns
        -------
        A tuple where: 
        The first element is True if each node in a network is in the property file, and otherwise False
        The second element is True if each node in the property file is in the network, and otherwise False
        """

        if isinstance(input, str):
            f = open(input, 'rU')
        else:
            f = input


        nfields=len(fieldNames)
        nNodesFound=0
        netHasAllNodes=True
        nodeLabelField=fieldNames.index('node_label')
        if hasHeader:
            f.readline()
        for i,line in enumerate(f):
            fields=line.split(splitterChar)
            if len(fields)!=nfields:
                f.close()
                raise Exception("Invalid number of fields on row: "+str(i+1))

            if isint(fields[nodeLabelField]):
                nodeLabel=int(fields[nodeLabelField]) # if node name/label is an integer, convert to int
            else:
                nodeLabel=fields[nodeLabelField]

            if nodeLabel in net:
                nNodesFound+=1
            else:
                netHasAllNodes=False
        f.close()
        
        fileHasAllNodes=(nNodesFound==len(net))

        return fileHasAllNodes,netHasAllNodes

    def addProperty(net,node,propertyName,theString):
        if isanum(theString):  # if it is a number 
            net.nodeProperty[propertyName][node]=float(theString)
            if isint(theString):  # if it is integer 
                net.nodeProperty[propertyName][node]=int(theString)
        else: # if it is a string
            net.nodeProperty[propertyName][node]=theString 

        
    if isinstance(input, str):
        f = open(input, 'rU')
    else:
        f = input

    #Read in the field names
    if propertyNames==None:
        line=f.readline() 
        fieldNames=line.strip().split(splitterChar)  # field names are taken from first line (header)
    else:
        fieldNames=list(propertyNames)    # fieldNames are given, copy field names to a new list
        if type(fieldNames)==str:
            fieldNames=[fieldNames] # if only a single field name string was given, 
                                    # convert it into a list (this makes it easier for the user 
                                    # if only one property is added, as he only needs to type 'property' 
                                    # and not ['property'] )

    if "node_label" not in fieldNames: # if node_label is not a field 
            fieldNames=["node_label"]+ fieldNames[:]    # add node_label as the first element  

    nfields=len(fieldNames)
    nodeLabelField=fieldNames.index("node_label")

    #enforce the rules of having no missing or extra nodes:
    fileHasAllNodes,netHasAllNodes=checkNodes(input,net,fieldNames,propertyNames==None,splitterChar)
    if fileHasAllNodes==False and not allowMissingData:
        f.close()
        raise Exception("The property file is missing some nodes that are in the network.")
    if netHasAllNodes==False and not allowExtraData:
        f.close()
        raise Exception("The property file has some extra nodes that are not in the network.")

    #Add the property names to the net
    for field in range(1,nfields):
        netext.addNodeProperty(net,fieldNames[field])

    #Add the properties for each node
    someNodesNotInNetwork=False
    for i,line in enumerate(f):        
        fields=line.strip().split(splitterChar)        

        assert len(fields)==nfields, "The number of fields on a row does not match"

        if isint(fields[nodeLabelField]):
            nodeName=int(fields[nodeLabelField]) # if node name/label is an integer, convert to int
        else:
            nodeName=fields[nodeLabelField]

        if nodeName in net:
            for field in range(1,nfields): #assumes that node label on the first field

                addProperty(net,nodeName,fieldNames[field],fields[field])

    f.close()

        
def saveNodeProperties(net,filename):
    plist=list(net.nodeProperty)
    f=open(filename,'w')

    #Write headers
    f.write("name")
    for p in plist:
        f.write(" "+p)
    f.write("\n")

    #Write values
    for node in net:
        f.write(str(node))
        for p in plist:
            f.write(" "+str(net.nodeProperty[p][node]))
        f.write("\n")


if __name__ == '__main__':
    """Run unit tests if called."""
    from tests.test_netio import *
    unittest.main()
