from __future__ import with_statement

# -*- coding: latin-1 -*-
import pynet
import netext,percolator,netio,transforms
import os
from pylab import *
import numpy as np
import copy
import random
import shutil
import operator
import matplotlib.colors
import pygraphviz

def normalizeValue(value,valueLimits):
    """Transforms a numerical value to the range (0,1).

    It is intended that the user should set valueLimits such that the
    true values fall between the limits. If this is not the case,
    values above given maxval or below given minval are truncated. The
    rest of the values are transformed linearly, such that the range
    (given minval, given maxval) becomes (0,1).
     
    normalizedValue= (true_val-given minval)/(given_maxval-given_minval)
    """ 
    if (valueLimits[0]-valueLimits[1]) == 0: 
        # If given minval and maxval are the same, all values will be
        # equal.
        normalizedValue=1
    elif value < valueLimits[0]:
        # If value is smaller than given minval
        normalizedValue=0
    elif value>valueLimits[1]:
        # If value is larger than given maxval
        normalizedValue=1
    else:
        normalizedValue=(value-valueLimits[0])/float(valueLimits[1] -
                                                     valueLimits[0])
    return normalizedValue 


def setColorMap(colorMap):
    """Set a colormap for edges.

    Two options of our own ('orange' and 'primary') are available in
    addition to the 150 pylab readymade colormaps (which can be listed
    with help matplotlib.cm ).

    Usage:
        myMap = setColorMap('bone')
    """
    if hasattr(colorMap, '_segmentdata'):
        return colorMap

    known_colormaps = ('primary', 'orange', 'bluered')
    if colorMap in known_colormaps:
        if colorMap == 'primary':
            # Jari's map: yellow->blue->red 
            segmentdata={'red': ( (0,1,1),(0.5,0,0), (1,1,1)  ),
                         'green': ( (0,1,1), (0.5,0.5,0.5), (1,0,0) ),
                         'blue': ( (0,0,0), (0.5,1,1), (1,0,0) )}
        elif colorMap=='orange':
            # Riitta's color map from white through yellow and orange to red 
            segmentdata = { 'red'  : ( (0.,.99,.99), 
                                       (0.2,.98,.98), 
                                       (0.4,.99,.99), 
                                       (0.6,.99,.99), 
                                       (0.8,.99,.99), 
                                       (1.0,.92,.92) ),
                            'green': ( (0,0.99,0.99), 
                                       (0.2,.89,.89),  
                                       (0.4,.80,.80), 
                                       (0.6,.50,.50), 
                                       (0.8,.33,.33), 
                                       (1.0,.10,.10) ),
                            'blue' : ( (0,.99,.99), 
                                       (0.2,.59,.59), 
                                       (0.4,.20,.20), 
                                       (0.6,0.0,0.0), 
                                       (0.8,0.0,0.0), 
                                       (1.0,.03,.03) )  }
        elif colorMap=='bluered':
            segmentdata={'red':  ( (0,0,0), 
                                   (0.17,0.25,0.25), 
                                   (0.33,0.7,0.7), 
                                   (0.5,.87,.87), 
                                   (0.67,.97,.97),  
                                   (0.83,.93,.93), 
                                   (1,.85,.85) ),
                         'green': ( (0,0,0), 
                                    (0.1667,0.53,0.53), 
                                    (0.3333,.8,.8), 
                                    (0.5,.9,.9), 
                                    (0.6667,.7,.7),
                                    (0.8333,.32,.32), 
                                    (1,.07,.07) ),
                         'blue': ( (0,.6,.6),  
                                   (0.1667,.8,.8),    
                                   (0.3333,1,1),    
                                   (0.5,.8,.8),    
                                   (0.6667,.33,.33),    
                                   (0.8333,.12,.12),
                                   (1,.05,.05) ) }
        myMap = matplotlib.colors.LinearSegmentedColormap(colorMap, segmentdata)
    else:
        try:
            myMap=get_cmap(colorMap)
        except AssertionError:
            comment = "Could not recognize given colorMap name '%s'" % colorMap
            raise AssertionError(comment)
    return myMap


def getConstantColorMap(rgb=(0,0,0)):
    """Return a colormap with constant color.

    Parameters
    ----------
    rgb : tuple (r, g, b)
        The color as RGB tuple. Each value must be between 0 and 1.

    Return
    ------
    cm : colorMap
        The colormap that has just one constant color.
    """
    cm={
        'red':  ( (0,rgb[0],rgb[0]), 
                  (1,rgb[0],rgb[0]) ),
        'green': ( (0,rgb[1],rgb[1]), 
                   (1,rgb[1],rgb[1]) ),
        'blue': ( (0,rgb[2],rgb[2]), 
                  (1,rgb[2],rgb[2]) ) }

    return matplotlib.colors.LinearSegmentedColormap("constant colormap",cm)  


# ---------------------------------------

def isListOfColors(theList):
    """
    Returns True if each element in the given list can be converted to a color
    by Matplotlib. Otherwise returns False.
    """
    cc=matplotlib.colors.ColorConverter()
    for element in theList:
        try:
            cc.to_rgb(element)
        except ValueError:
            return False
    return True
        

def getNodeColors(net,colorwith="strength",useColorMap="orange",parentnet=[]):
    """Returns a dictionary {node:color}. The colors are set based
    on either node strengh (colorwith="strength", default) 
    or any nodeProperty. For cases where e.g. nodes which have been thresholded
    out (k=0), the input parameter parentnet can be used - parentnet should contain the original
    network *before* thresholding, i.e. containing all original nodes and
    their attributes. IF parentnet is given, i) if strength is used, its nodes
    which are NOT in net colored gray, ii) if properties
    are used, its nodes are colored similarly to those nodes in net. Also the
    dictionary which is returned contains then all nodes in parentnet"""

    myNodeColors=setColorMap(useColorMap)

    nodeColors={}

    if colorwith=="strength":

        if hasattr(net,'matrixtype'):
            if net.matrixtype==0:        
                net=transforms.dist_to_weights(net)

        strengths = netext.strengths(net)
        max_value = max(strengths.values())
        min_value = min(strengths.values())

        if len(parentnet)>0:        # if we want the dict to include nodes not in net
            for node in parentnet:
                if node in net:     # if the node is in net, use its strength for color
                    nodeColors[node]=setColor(strengths[node],(min_value,max_value),myNodeColors)
                else:               # otherwise color it gray
                    nodeColors[node]=(0.5,0.5,0.5)
        else:
            for node in net:        # if parentnet not given, just color nodes by strength
                nodeColors[node]=setColor(strengths[node],(min_value,max_value),myNodeColors)
    else:

        numeric_props=netext.getNumericProperties(net)
        # first check if colorwith is a numeric property
        if colorwith in numeric_props:
            values=[]
            if len(parentnet)>0:    # again if we want to include nodes not in net
                for node in parentnet:  # first get min and max value of property
                    values.append(parentnet.nodeProperty[colorwith][node])

                min_value=min(values)
                max_value=max(values)
                for node in parentnet: # then set colors according to property
                    nodeColors[node]=setColor(parentnet.nodeProperty[colorwith][node],(min_value,max_value),myNodeColors)
            else:                   # otherwise do the above for nodes in net
                for node in net:
                    values.append(net.nodeProperty[colorwith][node])
                
                min_value=min(values)
                max_value=max(values)

                for node in net:
                    nodeColors[node]=setColor(net.nodeProperty[colorwith][node],(min_value,max_value),myNodeColors)
        
        else:
            # colorwith is not a numeric property, so look up unique values
            # and give them integer numbers

            values={} # empty dict for values
           
            if len(parentnet)>0:# if there are nodes not in net
                props=list(set(parentnet.nodeProperty[colorwith].values()))
            else:
                props=list(set(net.nodeProperty[colorwith].values()))

            #Check if properties can be converted to colors:
            if isListOfColors(props):
                propToColor={}
                cc=matplotlib.colors.ColorConverter()
                for p in props:
                    propToColor[p]=cc.to_rgb(p)
                if len(parentnet)>0:
                    for node in parentnet:
                        nodeColors[node]=propToColor[parentnet.nodeProperty[colorwith][node]]
                else:
                    for node in net:
                        nodeColors[node]=propToColor[parentnet.nodeProperty[colorwith][node]]                    
            else:
                for i,prop in enumerate(props):
                    values[prop]=i+1
                # now all property strings have a numerical value
                min_value=1
                max_value=max(values.values())
                if len(parentnet)>0:

                    for node in parentnet:
                        nodeColors[node]=setColor(values[parentnet.nodeProperty[colorwith][node]],(min_value,max_value),myNodeColors)
                else:
                    for node in net:
                        nodeColors[node]=setColor(values[net.nodeProperty[colorwith][node]],(min_value,max_value),myNodeColors)



    if len(nodeColors)==0:  # finally if for whatever reason no nodes were colored, just set them gray
        if len(parentnet)>0:
            for node in parentnet:
                nodeColors[node]=(0.5,0.5,0.5)
        else:
            for node in net:
                nodeColors[node]=(0.5, 0.5, 0.5)

    return nodeColors

# ------------------------------------------

def getNodeSizes(net,size_by="strength",minsize=2.0,maxsize=6.0):
    """Returns a dictionary {node:size} for visualizations. The sizes
    are set using either node strength"""

    nodeSizes={}

    if size_by=="strength":

        if hasattr(net,'matrixtype'):
            if net.matrixtype==0:        
                net=transforms.dist_to_weights(net)

        strengths = netext.strengths(net)
        maxs = max(strengths.values())
        mins = min(strengths.values())           

        if maxs==mins:
            A=0
        else:
            A=(maxsize-minsize)/(maxs-mins)
        B=maxsize-A*maxs

        for node in net:
            nodeSizes[node]=A*strengths[node]+B

    elif size_by=="fixed":
        for node in net:
            nodeSizes[node]=maxsize
    else:
        numeric_props=netext.getNumericProperties(net)
        if size_by in numeric_props:
            values=[]
            for node in net:
                values.append(net.nodeProperty[size_by][node])

            minval=min(values)
            maxval=max(values)

            if maxval==minval:
                A=0
            else:
                A=(maxsize-minsize)/(maxval-minval)

            B=maxsize-A*maxval
            for node in net:
                nodeSizes[node]=A*net.nodeProperty[size_by][node]+B

    return nodeSizes
          

def setColor(value,valueLimits,colorMap):
    """Converts a numerical value to a color.

    The value is scaled linearly to the range (0...1) using the
    function normalizeValue and the limits valueLimits. This scaled
    value is used to pick a color from the given colormap. The
    colormap should take in values in the range (0...1) and produce a
    three-tuple containing an RGB color, as in (r,g,b).
    """
    if valueLimits[0] < valueLimits[1]: 
        normalizedValue = normalizeValue(value,valueLimits)
        color = colorMap(normalizedValue)
    else:
        color=(0.5,0.5,0.5)  # gray if all values are equal
    return color


def setEdgeWidth(value,weightLimits,minwidth,maxwidth):
    """Transforms edge weights to widths in the range (minwidth,
    maxwidth). If given minwidth and maxwidth are the same, simply use
    that given width.
    """
    if not(weightLimits[0]-weightLimits[1])==0:
        # Normalizes the weight linearly to the range (0,1)
        normalizedWeight=normalizeValue(value,weightLimits)  
        # Transforms the normalized weight linearly to the range
        # (minwidth,maxwidth)
        width=minwidth+normalizedWeight*(maxwidth-minwidth)   
    else:
        # If given minwidth and maxwidth are the same, simply use that width.
        width=minwidth 
    return width

def get_edge_angle(xcoords, ycoords):
    """Return the edge angle in [-pi/2, 3*pi/2]."""
    dx, dy = xcoords[1]-xcoords[0], ycoords[1]-ycoords[0]
    if dx == 0:
        if dy > 0: theta = np.pi/2
        else: theta = -np.pi/2
    elif dx > 0:
        theta = np.arctan(dy/dx)
    else:
        theta = np.arctan(dy/dx) + np.pi
    return theta

def points_to_data(ax, x_points):
    """Converts point length `x_points` (1/72th of inch) to length
    in data coordinates. Note that this only works (in general) if
    the axis are equal; otherwise converting x and y-coordinates
    gives different results."""
    ax_pos = ax.get_position()
    V = ax.axis()
    width_inches = ax_pos.width*ax.get_figure().get_figwidth()
    return (V[1]-V[0])*x_points/(72*width_inches)

def draw_edge(axes, xcoords, ycoords, width, color, 
              symmetric, zorder, nodesize):
    """Used by visualizeNet to draw edges."""
    if symmetric:
        line_obj, = axes.plot(xcoords, ycoords, '-', lw=width,
                              color=color, zorder=zorder)
    else:
        dx, dy = xcoords[1]-xcoords[0], ycoords[1]-ycoords[0]
        theta = get_edge_angle(xcoords, ycoords)

        x_diff = points_to_data(axes, 0.5*nodesize*np.cos(theta))
        y_diff = points_to_data(axes, 0.5*nodesize*np.sin(theta))
        x, y = xcoords[0]+x_diff, ycoords[0]+y_diff
        dx, dy = dx-2*x_diff, dy-2*y_diff

        #print ("Start (%.4f, %.4f), End (%.4f, %.4f)" 
        #       % (x, y, xcoords[0]+dx, ycoords[0]+dy))

        arrow_width = points_to_data(axes, width)

        # See pylab.matplotlib.patches.FancyArrow for
        # documentation of the arrow command.
        line_obj = axes.arrow(x, y, dx, dy,
                              color=color, width=arrow_width,
                              head_width=10*arrow_width,
                              head_length=14*arrow_width,
                              shape='left', overhang=0.05,
                              length_includes_head=True, 
                              head_starts_at_zero=False,
                              zorder=zorder)
    return line_obj

def draw_node(axes, x, y, shape, color, size, edgecolor, edgewidth, zorder):
    """Used by visualizeNet to draw nodes."""
    #if size == 0: return
    node_obj, = axes.plot([x], [y], shape,
                          markerfacecolor=color,
                          markeredgecolor=edgecolor,
                          markeredgewidth=edgewidth,
                          markersize=size,
                          zorder=zorder)
    return node_obj

from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection

class NetNodes(object):
    """Helper for drawing nodes in visualizeNet."""
    def __init__(self, ax):
        self.nodes = []
        self.ax = ax

    def add(self, x, y, color, size, edgecolor, edgewidth, zorder):
        self.nodes.append(Circle((x,y), 0.5*points_to_data(self.ax, size),
                                 facecolor=color, edgecolor=edgecolor,
                                 linewidth=edgewidth, zorder=zorder))

    def __del__(self):
        self.ax.add_collection(PatchCollection(self.nodes, match_original=True))

from matplotlib.path import Path
from matplotlib.patches import PathPatch, FancyArrow

class NetEdges(object):
    """Helper for drawing edges in visualizeNet."""
    def __init__(self, ax, symmetric):
        self.edges = []
        self.ax = ax
        self.symmetric = symmetric

    def add(self, xcoords, ycoords, width, color, zorder, nodesize):
        if self.symmetric:
            self.edges.append(PathPatch(Path(zip(xcoords,ycoords)),
                                        linewidth=width, color=color, zorder=zorder))
        else:
            dx, dy = xcoords[1]-xcoords[0], ycoords[1]-ycoords[0]
            theta = get_edge_angle(xcoords, ycoords)

            x_diff = points_to_data(self.ax, 0.5*nodesize*np.cos(theta))
            y_diff = points_to_data(self.ax, 0.5*nodesize*np.sin(theta))
            x, y = xcoords[0]+x_diff, ycoords[0]+y_diff
            dx, dy = dx-2*x_diff, dy-2*y_diff

            #print ("Start (%.4f, %.4f), End (%.4f, %.4f)" 
            #       % (x, y, xcoords[0]+dx, ycoords[0]+dy))

            arrow_width = points_to_data(self.ax, width)

            # See pylab.matplotlib.patches.FancyArrow for
            # documentation of the arrow command.
            self.edges.append(FancyArrow(x, y, dx, dy,
                                         color=color, width=arrow_width,
                                         head_width=10*arrow_width,
                                         head_length=14*arrow_width,
                                         shape='left', overhang=0.05,
                                         length_includes_head=True, 
                                         head_starts_at_zero=False,
                                         zorder=zorder))
    def __del__(self):
        self.ax.add_collection(PatchCollection(self.edges, match_original=True))

def visualizeNet(net, coords=None, axes=None, frame=False, animated=False,
                 scaling=True, margin=0.025,
                 nodeShapes=None, defaultNodeShape='o',
                 nodeSizes=None, defaultNodeSize=None,
                 nodeColors=None, defaultNodeColor=None,
                 nodeEdgeColors=None, defaultNodeEdgeColor='black',
                 nodeLabels=None, nodeLabelSize=None, labelAllNodes=False,
                 labelPositions=None, defaultLabelPosition='out',
                 edgeColors=None, defaultEdgeColor=None,
                 edgeWidths=None, defaultEdgeWidth=None,
                 nodeEdgeWidths=None, defaultNodeEdgeWidth=0.2,
                 edgeLabels=None, edgeLabelSize=None, labelAllEdges=False,
                 nodePlotOrders=None, defaultNodePlotOrder=1,
                 edgePlotOrders=None, defaultEdgePlotOrder=0):
    """Visualize a network.
    
    Note that all sizes (node size, link width, etc.) are given in
    points: one point equals 1/72th of an inch.

    Basic parameters
    ----------------
    net : pynet.SymmNet
        The network to visualize
    coords : dictionary of tuples {node_ID: (x,y)}
        Coordinates of all nodes. If None, the coordinates will be
        calculated. The x and y coordinates are assumed to
        have the same scale, so for example an edge between nodes 0
        and 1 with `coords[0]=(0,0)` and `coords[1]=(2,2)` will be at
        an angle of 45 degrees.
    axes : pylab.axes object
        If given, the network will be drawn in this axis. Otherwise a
        new figure is created for the plot and the figure handle is
        then returned.
    frame : bool
        If False, the frame will be not be shown in the plot. You can
        still use axis labels.
    scaling : bool
        If True, the coordinate axes will be scaled for best fit. If
        false, the coordinate axes will not be altered.
    margin : float (>= 0)
        The relative size of empty margin around the network. Margin
        of 0.0 means that some nodes touch the edge of the plot,
        margin of 0.2 adds 20 % on all sides etc. This parameter has
        an effect only if scaling is True.

    Defining node and edge colors
    -----------------------------

    Colors for nodes and edges are defined similarly. The following
    explains the procedure for nodes; to control edge coloring simply
    replace the word 'node' (or 'Node') with the word 'edge' (or
    'Edge') in the parameter names.

    The color of a node is defined with the dictionary `nodeColors`:
    key is the node index and the value is any valid coloring scheme
    (see below under 'Coloring schemes'). If a node index is not in
    `nodeColors`, it is colored according to `defaultNodeColor`. This
    variable can have the same values as the values in `nodeColors`.

    Coloring schemes
    ----------------

    A constant color can be defined in any way allowed by pylab. For
    example 'k', 'black' and (0,0,0) all give black color.

    Alternatively the color can be based on the node strength, degree
    or any node property. In this case the coloring definition is a
    dictionary. The following examples illustrate the idea:
    
    color_scheme = {'by':'weight', 'scale':'log', 'cmap':'winter'}
    color_scheme = {'by':'degree', 'scale':'lin', 'min':1, 'max':10}
    color_scheme = {'by':'property:myProperty', 'scale':'log'}

    The possible keys and their default values are
    
        KEY      DEFAULT VALUE       OTHER POSSIBLE VALUES
        'by'     'strength'/'weight' 'degree', 'property:<property_name>'
        'scale'  'log'               'lin'
        'cmap'   'jet'               Any colormap
        'min'    (Min value in data) Any integer x,     1 <= x <= 'max'
        'max'    (Max value in data) Any integer x, 'min' <= x         

    Any keys that are omitted are filled in with the default
    value. Note the syntax for using node properties, where the word
    'property' is followed by a semicolon and the property name.

    Node size
    ---------

    The node size is controlled with a syntax similar to that used
    with colors. Node size is defined by dictionary `nodeSizes`, and
    if a node is not in it, the default value given by
    `defaultNodeSize` is used. 

    The value can be a single integer, which gives the node size in
    pixels. Alternatively the node size can be controlled by node
    strength, degree or any property:
    
    node_size_scheme = {'by':'strength', 'scale':'log', 'min':2, 'max':10}
    node_size_scheme = {'by':'degree', 'scale':'lin'}
    node_size_scheme = {'by':'property:myProperty', 'scale':'log'}

    The possible keys and their default values are
        KEY        DEFAULT VALUE       OTHER POSSIBLE VALUES
        'by'       'strength'          'degree', 'property:<property_name>'
        'scale'    'log'               'lin'
        'min'      (Min value in data) int;          1 <= x <= 'max'
        'max'      (Max value in data) int;      'min' <= x         
        'min_size'   1                 int;          1 <= x <= 'max_size'
        'max_size'   6                 int; 'min_size' <= x         
    Again, keys that are omitted are filled with default values.

    Edge width
    ----------

    Edge width is defined by dictionary `edgeWidths`, and if an edge
    is not in it, the default value given by `defaultEdgeWidth` is
    used.

    The value can be a single integer, which gives the edge width in
    pixels. Alternatively the edge width can be controlled by edge
    weight:
    
    edge_width_scheme = {'by':'weight', 'scale':'log', 'min':1, 'max':5}

    The possible keys and their default values are
    
        KEY         DEFAULT VALUE       OTHER POSSIBLE VALUES
        'by'        'weight'          
        'scale'     'log'               'lin'
        'min'       (Min value in data) int;   1 <= x <= 'max'
        'max'       (Max value in data) int;   'min' <= x         
        'min_size'    0.2               float; 1 <= x <= 'max_width'
        'max_size'    2.0               float; 'min_width' <= x         

    Note that the 'by'-key can always be omitted since it has only one
    possible value.

    Node labels
    -----------

    Node labels can be given in `nodeLabels` dictionary, where the key
    is node index and the value is the corresponding labels. If
    `labelAllNodes` is True, also nodes not in `nodeLabels` will
    receive a label, which is the node index.

    The values in `nodeLabels` are converted to string with str().

    There are two possibilities for the positioning of node labels:
    'in' and 'out' (default). 'in' means that the label is printed at
    the exact position of the node; if the node is hollow, this
    effectively prints the node labes inside the nodes (make sure the
    node size and font sizes are compatible). 'out' prints the label
    next to the node.

    Edge labels
    -----------

    Also edges can have labels, given in `edgeLabels` dictionary,
    where the key is a tuple (i,j) of end node indices. The edge
    labels are always printed on the right side of each edge
    (with direction defined from i to j).

    If `labelAllEdges` is True, also the edges not listed in
    `edgeLabels` will be given a label. In this case the label is
    '(i,j)', where i and j are the indices of the end nodes.

    Return
    ------
    fig : pylab.Figure (None if `axes` is given.)
        Figure object with one axes containing the plotted network
        figure.

    Examples
    --------
    >>> # Construct an example network.
    >>> from netpython import pynet, visuals
    >>> net = pynet.SymmNet()
    >>> net[0][1] = 1.0
    >>> net[1][2] = 3.5
    >>> net[0][2] = 5.0

    >>> # Simplest case: get coordinates, plot into a
    >>> # new figure and save it to disk
    >>> fig = visuals.visualizeNet(net)
    >>> fig.savefig('myNet.eps')

    >>> # Draw the figure in the upper left subfigure, with predefined
    >>> # coordinates. Note that drawNet does not return anything.
    >>> import pylab
    >>> coords = {0:(0,0), 1:(4,0), 2:(2,3)}
    >>> fig = pylab.figure()
    >>> ax = fig.add_subplot(2,2,1)
    >>> visuals.visualizeNet(net, coords=coords, axes=ax)
    """

    #
    # DEFAULT VALUES. These will be used whenever the user has not
    # defined a given value for defaultNodeColor etc.
    # 
    internal_defaultNodeColor = {'by':'strength', 'scale':'log', 'cmap':'jet'}
    internal_defaultEdgeColor = {'by':'weight', 'scale':'log', 'cmap':'jet'}

    internal_defaultNodeSize = {'by':'strength', 'scale':'log',
                                'min_size':2, 'max_size':6}
    internal_defaultEdgeWidth = {'by':'weight', 'scale':'log',
                                 'min_size':0.2, 'max_size':2.0}

    node_label_font_color = 'k'
    node_label_font_size = (8 if nodeLabelSize==None else nodeLabelSize)
    edge_label_font_color = None #'k'
    edge_label_font_size = (5 if edgeLabelSize==None else edgeLabelSize)

    #
    # PROCESS INPUT PARAMETERS
    #
    
    if coords is None:
        coords = calculateCoordinates(net)

    fig = None
    if axes is None:
        fig = figure()
        axes = fig.add_subplot(111)

    nodeShapes = (nodeShapes or {})

    nodeColors = (nodeColors or {})
    defaultNodeColor = (defaultNodeColor or {})
    if isinstance(defaultNodeColor, dict):
        for k,v in internal_defaultNodeColor.iteritems():
            if k not in defaultNodeColor:
                defaultNodeColor[k] = v

    nodeEdgeColors = (nodeEdgeColors or {})
                
    edgeColors = (edgeColors or {})
    defaultEdgeColor = (defaultEdgeColor or {})
    if isinstance(defaultEdgeColor, dict):
        for k,v in internal_defaultEdgeColor.iteritems():
            if k not in defaultEdgeColor:
                defaultEdgeColor[k] = v
    
    nodeSizes = (nodeSizes or {})
    if defaultNodeSize is None:
        defaultNodeSize = {}
    if isinstance(defaultNodeSize, dict):
        for k,v in internal_defaultNodeSize.iteritems():
            if k not in defaultNodeSize:
                defaultNodeSize[k] = v

    edgeWidths = (edgeWidths or {})
    if defaultEdgeWidth is None:
        defaultEdgeWidth = {}
    if isinstance(defaultEdgeWidth, dict):
        for k,v in internal_defaultEdgeWidth.iteritems():
            if k not in defaultEdgeWidth:
                defaultEdgeWidth[k] = v

    nodeEdgeWidths = (nodeEdgeWidths or {})

    nodeLabels = (nodeLabels or {})
    labelPositions = (labelPositions or {})

    edgeLabels = (edgeLabels or {})

    nodePlotOrders = (nodePlotOrders or {})
    edgePlotOrders = (edgePlotOrders or {})

    if margin < 0: margin = 0.0

    # Initialize dictionary where the returned plotted artist objects will be collected.
    obj = dict()

    #
    # AUXILIARY FUNCTIONS
    #

    def scaled(scaling_type, value, value_limits, final_limits):

        def lin_scaling(value, value_limits, final_limits):
            value_span = value_limits[1] - value_limits[0]
            final_span = final_limits[1] - final_limits[0]
            if final_span == 0:
                return final_limits[0]
            if value_span == 0:
                p = 0.5
            else:
                p = float(value - value_limits[0])/value_span
            return final_limits[0]+p*final_span

        if value <= value_limits[0]:
            return final_limits[0]
        if value >= value_limits[1]:
            return final_limits[1]

        if scaling_type == 'log' or scaling_type == 'logarithmic':
            return lin_scaling(np.log(value),
                               np.log(value_limits),
                               final_limits)
        else:
            return lin_scaling(value, value_limits, final_limits)

    def determine_size(scheme, i, net, values, limits, defaults):
        if not isinstance(scheme, dict):
            return scheme
        else:
            # Determine what defines the size. Calculate the limits
            # for this property if not yet done.
            size_by = scheme.get('by', defaults['by'])
            if size_by not in limits:
                property_name = "".join(size_by.split(':')[1:])
                np_ = sorted(net.nodeProperty[property_name].values())
                limits[size_by] = (np_[0], np_[-1])
            if size_by not in values:
                property_name = "".join(size_by.split(':')[1:])
                values[size_by] = net.nodeProperty[property_name][i]
                
            scale = scheme.get('scale', defaults['scale'])
            val_min = scheme.get('min', limits[size_by][0])
            val_max = scheme.get('max', limits[size_by][1])
            size_min = scheme.get('min_size', defaults['min_size'])
            size_max = scheme.get('max_size', defaults['max_size'])

            #print size_by, scale, val_min, val_max, size_min, size_max
            return scaled(scale, values[size_by], [val_min, val_max],
                          [size_min, size_max])
                    
    def determine_color(scheme, i, net, values, limits, defaults):
        if not isinstance(scheme, dict):
            return scheme
        else:
            color_by = scheme.get('by', defaults['by'])
            if color_by not in limits:
                property_name = "".join(color_by.split(':')[1:])
                np_ = sorted(net.nodeProperty[property_name].values())
                limits[color_by] = (np_[0], np_[-1])
            if color_by not in values:
                property_name = "".join(color_by.split(':')[1:])
                values[color_by] = net.nodeProperty[property_name][i]
                
            scale = scheme.get('scale', defaults['scale'])
            cmap = scheme.get('cmap', defaults['cmap'])
            val_min = scheme.get('min', limits[color_by][0])
            val_max = scheme.get('max', limits[color_by][1])

            cval = scaled(scale, values[color_by],
                          [val_min, val_max], [0.0,1.0])
            cm = setColorMap(cmap)
            return cm(float(cval))

    def luminance(c):
        """Return luminance of color `c`."""
        c_vec = matplotlib.colors.colorConverter.to_rgb(c)
        return 0.2126*c_vec[0]+0.7152*c_vec[1]+0.0722*c_vec[2]

    def edge_label_pos(axes, xcoords, ycoords, edge_width, label_size, offset=1.5):
        """Return the baseline position and label rotation (in angles)
        for an edge label. The label will be on the right side of the
        edge, in proper orientation for reading, and located `offset`
        points from the edge.
        """
        theta = get_edge_angle(xcoords, ycoords)
        if theta > -np.pi/2 and theta < np.pi/2:
            # Edge goes from left to right.
            label_rotation = theta
        else:
            # Edge goes from right to left.
            label_rotation = theta - np.pi

        offset_points = offset+0.5*edge_width+0.5*label_size
        label_position = (0.5*sum(xcoords), 0.5*sum(ycoords))
        offset_dir = theta - np.pi/2
        label_offset = (offset_points*np.cos(offset_dir),
                        offset_points*np.sin(offset_dir))
        return label_position, label_offset, 180*label_rotation/np.pi


    #
    # INITIALIZE SOME DATA STRUCTURES
    #
        
    # Find out the minimum and maximum value for strength and degree.
    strengths = netext.strengths(net)
    smin, smax = min(strengths.values()), max(strengths.values())
    degrees = netext.deg(net)
    dmin, dmax = min(degrees.values()), max(degrees.values())

    limits = {"strength":(smin, smax), "degree":(dmin,dmax)}

    #
    # INITIALIZE AXES SIZE
    #

    node_diameters = {};
    for nodeIndex in net:
        values = {"strength": strengths[nodeIndex],
                  "degree": degrees[nodeIndex]}
        node_diameters[nodeIndex] = (determine_size(nodeSizes.get(nodeIndex, defaultNodeSize),
                                                    nodeIndex, net, values, limits,
                                                    defaultNodeSize) +
                                     determine_size(nodeEdgeWidths.get(nodeIndex, 
                                                                       defaultNodeEdgeWidth),
                                                    nodeIndex, net, values, limits,
                                                    defaultNodeEdgeWidth))

    if scaling:
        # Make axis equal making sure the nodes on the edges are not
        # clipped. We cannot use `axis('equal')` because nothing has been
        # drawn yet, but we still need to do this so the arrows will be
        # draw properly in the plotting phase. This is a bit tricky
        # because the axes might not be square.

        max_node_diameter = max(node_diameters.values())
        y_coords = sorted(map(operator.itemgetter(1), coords.values()))
        x_coords = sorted(map(operator.itemgetter(0), coords.values()))
        ax_pos = axes.get_position()
        x_span = x_coords[-1] - x_coords[0]
        y_span = y_coords[-1] - y_coords[0]
        ax_width_inches = ax_pos.width*axes.get_figure().get_figwidth()
        ax_height_inches = ax_pos.height*axes.get_figure().get_figheight()

        if (x_span*ax_height_inches > y_span*ax_width_inches):
            # The x-span dictates the coordinates. Calculate the margin
            # necessary to fit in the nodes on the edges.
            rad_frac = 0.5*max_node_diameter/(72*ax_width_inches)
            x_margin = x_span*rad_frac/(1-2*rad_frac)
            x_min, x_max = x_coords[0]-x_margin, x_coords[-1]+x_margin
            y_mid = 0.5*(y_coords[-1] + y_coords[0])
            y_axis_span = (x_max-x_min)*ax_height_inches/ax_width_inches
            y_min, y_max = y_mid - 0.5*y_axis_span, y_mid + 0.5*y_axis_span,
        else:
            # The y-span dictates the coordinates. Calculate the margin
            # necessary to fit in the nodes on the edges.
            rad_frac = 0.5*max_node_diameter/(72*ax_height_inches)
            y_margin = y_span*rad_frac/(1-2*rad_frac)
            y_min, y_max = y_coords[0]-y_margin, y_coords[-1]+y_margin
            x_mid = 0.5*(x_coords[-1] + x_coords[0])
            x_axis_span = (y_max-y_min)*ax_width_inches/ax_height_inches
            x_min, x_max = x_mid - 0.5*x_axis_span, x_mid + 0.5*x_axis_span,

        y_span, x_span = y_max - y_min, x_max - x_min
        axes.set_ylim(ymin=y_min-margin*y_span, ymax=y_max+margin*y_span)
        axes.set_xlim(xmin=x_min-margin*x_span, xmax=x_max+margin*x_span)
    prev_autoscale = axes.get_autoscale_on()
    axes.set_autoscale_on(False)

    point_trans = lambda x_: points_to_data(ax,x_) # Transform points to data coordinates.


    #
    # DRAW EDGES
    #

    edges = list(net.edges)
    if edges:
        # Sort by edge weight.
        edges.sort(key=operator.itemgetter(2))
        limits['weight'] = (edges[0][2], edges[-1][2])

        # DEBUGGING: Print statistics of edge weights.
        #import data_utils
        #weights_ = map(operator.itemgetter(2), edges) 
        #print "Edges:", limits['weight'], " 10/25/50/75/90:",
        #mp_ = len(edges)/2
        #print "%d, %d, %d, %d, %d" % (int(data_utils.percentile(weights_, 0.1)),
        #                              int(data_utils.percentile(weights_, 0.25)),
        #                              int(data_utils.percentile(weights_, 0.5)),
        #                              int(data_utils.percentile(weights_, 0.75)),
        #                              int(data_utils.percentile(weights_, 0.9)))
        net_edges = NetEdges(axes, net.isSymmetric())
        for i,j,w in edges:
            values = {"weight": w,
                      "strength": strengths[j],
                      "degree": degrees[j]}
            
            # Determine edge width.
            if (i,j) in edgeWidths:
                width = determine_size(edgeWidths[(i,j)], (i,j), net,
                                       values, limits, defaultEdgeWidth)
            elif (j,i) in edgeWidths:
                width = determine_size(edgeWidths[(j,i)], (j,i), net,
                                       values, limits, defaultEdgeWidth)
            else:
                width = determine_size(defaultEdgeWidth, (i,j), net,
                                       values, limits, defaultEdgeWidth)

            # Determine edge color.
            if (i,j) in edgeColors:
                color = determine_color(edgeColors[(i,j)], (i,j), net,
                                       values, limits, defaultEdgeColor)
            elif (j,i) in edgeColors:
                color = determine_color(edgeColors[(j,i)], (j,i), net,
                                       values, limits, defaultEdgeColor)
            else:
                color = determine_color(defaultEdgeColor, (j,i), net,
                                       values, limits, defaultEdgeColor)

            if (i,j) in edgePlotOrders:
                zorder = edgePlotOrders[(i,j)]
            elif (j,i) in edgePlotOrders:
                zorder = edgePlotOrders[(j,i)]
            else:
                zorder = defaultEdgePlotOrder

            # FOR DEBUGGING:
            #print "Edge (%d,%d / %.2f) : %.1f %s %f" % (i,j,w,width,str(color),zorder)
            xcoords, ycoords = (coords[i][0], coords[j][0]), (coords[i][1], coords[j][1])
            #obj[(i,j)] = draw_edge(axes, xcoords, ycoords, width, color,
            #                       net.isSymmetric(), zorder, node_diameters[j])
            net_edges.add(xcoords, ycoords, width, color, zorder, node_diameters[j])

            # Add edge label.
            if (labelAllEdges or (i,j) in edgeLabels or 
                (net.isSymmetric() and (j,i) in edgeLabels)):
                if (i,j) in edgeLabels:
                    label = str(edgeLabels[(i,j)])
                elif (net.isSymmetric() and (j,i) in edgeLabels):
                    label = str(edgeLabels[(j,i)])
                else:
                    label = "(%d,%d)" % (i,j)

                lpos, loffset, lrot = edge_label_pos(axes, xcoords, ycoords,
                                                     width, edge_label_font_size)
                axes.annotate(label, lpos, xytext=loffset,
                              textcoords='offset points',
                              color=edge_label_font_color,
                              size=edge_label_font_size,
                              horizontalalignment='center',
                              verticalalignment='center',
                              rotation=lrot,
                              zorder=zorder+0.5)

    #
    # DRAW NODES
    #

    max_node_diameter = 0;
    net_nodes = NetNodes(axes)
    node_internal_color = {}
    for nodeIndex in net:
        values = {"strength": strengths[nodeIndex],
                  "degree": degrees[nodeIndex]}

        # Determine node shape.
        shape = nodeShapes.get(nodeIndex, defaultNodeShape)

        # Determine node size (and update max size).
        size = determine_size(nodeSizes.get(nodeIndex, defaultNodeSize),
                              nodeIndex, net, values, limits,
                              defaultNodeSize)

        # Determine node edge width.
        edgewidth = determine_size(nodeEdgeWidths.get(nodeIndex, defaultNodeEdgeWidth),
                                   nodeIndex, net, values, limits,
                                   defaultNodeEdgeWidth)
        max_node_diameter = max(max_node_diameter, size+edgewidth)

        # Determine node color
        color = determine_color(nodeColors.get(nodeIndex, defaultNodeColor),
                                nodeIndex, net, values, limits,
                                defaultNodeColor)
        node_internal_color[nodeIndex] = color

        # Determine node edge color
        edgecolor = determine_color(nodeEdgeColors.get(nodeIndex, defaultNodeEdgeColor),
                                    nodeIndex, net, values, limits,
                                    defaultNodeEdgeColor)

        # Determine z-order.
        zorder = nodePlotOrders.get(nodeIndex, defaultNodePlotOrder)
    
        # FOR DEBUGGING:
        #print "Node %d : %f %s %f %s" % (nodeIndex, size, str(color), edgewidth, str(edgecolor))

        #obj[nodeIndex] = draw_node(axes, coords[nodeIndex][0], coords[nodeIndex][1],
        #                           shape, color, size, edgecolor, edgewidth, zorder)
        net_nodes.add(coords[nodeIndex][0], coords[nodeIndex][1],
                      color, size, edgecolor, edgewidth, zorder)

        # Add node label.
        if nodeIndex in nodeLabels or labelAllNodes:
            if nodeIndex in nodeLabels:
                label = str(nodeLabels[nodeIndex])
            else:
                label = str(nodeIndex)

            label_pos = labelPositions.get(nodeIndex, defaultLabelPosition)
            if label_pos == 'out':
                nodeLabel_offset = int(np.ceil(float(size)/2))+1
                axes.annotate(label,
                              (coords[nodeIndex][0],coords[nodeIndex][1]),
                              textcoords='offset points',
                              xytext=(nodeLabel_offset, nodeLabel_offset),
                              color=node_label_font_color,
                              size=node_label_font_size,
                              zorder=zorder+0.5)
            elif label_pos == 'in':
                axes.annotate(label,
                              (coords[nodeIndex][0],coords[nodeIndex][1]),
                              color=("black" if luminance(node_internal_color[nodeIndex]) > 0.5 else "white"),
                              size=max(5, min(size-1, 0.7*size)),
                              horizontalalignment='center',
                              verticalalignment='center',
                              zorder=zorder+0.5)
 
    # Remove frame.
    if not frame:
        # Using 'axes.set_axis_off()' would also turn of the axis
        # labels, which is too much. The following lines are required
        # to turn off the frame and tick labels while keeping axis
        # labels.
        axes.set_frame_on(False)
        axes.set_xticklabels([])
        axes.xaxis.set_ticks_position('none')
        axes.set_yticklabels([])
        axes.yaxis.set_ticks_position('none')

    # Return autoscaling to the original value.
    axes.set_autoscale_on(prev_autoscale)

    # Return figure (or objects drawn).
    return (obj if animated else axes.get_figure())
    

def calculateCoordinates(net):
    """Calculate coordinates for a network."""
    # Construct a pygraphviz graph.
    G = pygraphviz.AGraph(strict=True, directed=False)
    for i,j,w in net.edges:
        G.add_edge(i,j)
    ## Calculate coordinates and extract them.
    G.layout()
    coords = dict((int(i), tuple(map(float, G.get_node(i).attr["pos"].split(","))))
                  for i in G.nodes_iter())
    return coords
