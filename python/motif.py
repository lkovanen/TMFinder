"""General functions for representing and plotting motifs."""

import sys
from netpython import pynet
from netpython import visuals
from netpython import transforms
import numpy as np
from PyBliss import Graph
from digraph import Digraph
import operator
import data_utils
import pickle
import itertools
import matplotlib
import pygraphviz as pgv

def nof_lines(fileName):
    c = 0
    with open(fileName,'r') as f:
        for line in f:
            c += 1
    return c

# Do not use the C++ interface so we can use pickle.
DirNet = pynet.ScipySparseDirNet

class Bunch(object):
    """Empty class for holding named values."""
    def __init__(self, **kwds):
        self.__dict__.update(kwds)


#def copynet(g):
#    g2 = DirNet()
#    for n1 in g:
#        g2.addNode(n1)
#        for n2 in g[n1].iterOut():
#            g2[n1][n2] = g[n1][n2]
#    return g2
copynet = transforms.copyNet
       
def get_components(g):
    """Find connected components of graph `g`."""
    unseen_nodes = set(g._nodes.keys())
    while unseen_nodes:
        to_process = set([unseen_nodes.pop()])
        component_nodes = []
        while to_process:
            component_nodes.append(to_process.pop())
            for neighbor in g[component_nodes[-1]]:
                if neighbor in unseen_nodes:
                    unseen_nodes.discard(neighbor)
                    to_process.add(neighbor)

        yield transforms.getSubnet(g, component_nodes)

def is_connected(g, first_node=None):
    """Check if graph `g` is connected."""
    unseen_nodes = set(g._nodes.keys())
    if first_node is not None:
        to_process = set([first_node])
    else:
        to_process = set([unseen_nodes.pop()])
    while to_process:
        node = to_process.pop()
        for neighbor in g[node]:
            if neighbor in unseen_nodes:
                unseen_nodes.discard(neighbor)
                to_process.add(neighbor)

    return not bool(unseen_nodes)

def construct_directed_hash(edges):
    """Shortcut for construct_hash(edges, 'digraph')."""
    return construct_hash(edges, 'digraph')

def construct_hash(events, graph_type):
    """Construct a graph and return canonical hash.

    A motif 'm' is isomorphic to a edge list 'edges' if
    'm.get_directed_hash(m.basic_type_map)' equals
    'construct_directed_hash(edges)'.

    Parameters
    ----------
    events : list of pairs
       List of events (node pairs) that define the topology of the
       motif. For example [(0,1),(1,2),(2,0)] would correspond to a
       3-cycle. The labels of the nodes (i.e. 0,1, and 2 in this
       example) have no relevance, they are only used to describe the
       motif.
    graph_type : str
       The graph type for which the hash is calculated, one of
       "motif", "multigraph", "digraph" or "graph".

    Return
    ------
    hash : int
       The hash of the canonical form of the graph defined by the edge
       list.
    """
    # To make sure everything works exactly as with actual motifs the
    # hash is calculated by constructing a motif object.
    events = map(tuple,events)
    all_nodes = set(reduce(operator.add, events))
    node_map = dict(zip(all_nodes,range(len(all_nodes))))
    event_ids = zip(events, range(len(node_map),len(node_map)+len(events)))
    
    node_id_str = ",".join(["%d:0" % i for i in node_map.itervalues()])
    event_id_str = ",".join(["%d:1" % i for (n1,n2),i in event_ids])
    edges_str=""
    i_ev_prev = None
    for (n1,n2),i_ev in event_ids:
        edges_str += " %d,%d %d,%d" % (node_map[n1], i_ev, i_ev, node_map[n2])
        if i_ev_prev is not None:
            edges_str += " %d,%d" % (i_ev_prev, i_ev)
        i_ev_prev = i_ev

    motif_line = "0 0 %d [%s,%s]%s" % (len(node_map)+len(events), node_id_str, event_id_str, edges_str)
    #print motif_line
    m = Motif(motif_line, node_types={0:Bunch(name='node')}, event_types={1:Bunch(name='event')})
    return m.get_hash(graph_type=graph_type, type_map=m.basic_type_map)

class Motif(object):
    """A temporal motif.

    After initialization the following variables can be accessed:
    self.count : int
        The total number of this motif found.
    self.ref_count : float
        The total count in the reference.
    self.n_lt_ref : int
        Number of times the count is smaller in the data than in the
        reference.
    self.nof_nodes : int
        Number of nodes in this motif.
    self.nof_events : int
        Number of events in this motif.
    self.nof_edges : int
        Number of edges in the undirected underlying graph.
    self.nof_dir_edges : int
        Number of edges in the directed underlying graph.
    self.nodes : {int: int}
        Dictionary of nodes. The keys are the node ids obtained for
        the canonical graph of the underlying undirected and untyped
        graph, and values are the node colors.
    self.events : {int: (int, int, int)}
        Dictionary of events. The keys are event ids in
        [0,N_events], and the value is a tuple (node_1, node_2,
        color), where the node ids correspond to those in the
        dictionary `nodes` and color is the event color.
    self.node_net : pynet.Net
        A weighted, directed graph where vertices correspond to nodes
        and edge weights tell the number of events between them.
    self.event_net : pynet.Net
        An unweighted, directed graph where vertices correspond to
        events and the edges denote their temporal order.
    self.topological_hash : int
        The canonical hash of the underlying undirected simple graph
        when node types are irrelevant. This can be used to plot
        motifs so that similar motifs look the same.

    The user node ids are identical in all motifs with the same
    topological_hash.
    """

    # Possible graph types when calculating isomorphism. These have
    # been sorted from the most to least specific.
    graph_types = ('motif','multigraph','digraph','graph')

    # Possible node roles in each graph type. These are updated in
    # __read_vertices() if it turns out that there is a motif with
    # more nodes than expected.
    N_nodes_max = 2
    node_roles = {'motif': ("C","R","T","N"),
                  'multigraph': ("C","R","T"),
                  'digraph': ("C","R","T"),
                  'graph': range(N_nodes_max)}

    def __init__(self, line, node_types=None, event_types=None, header_line=None):
        """Initialize motif.

        Parameters
        ----------
        line : str
          The line that describes the motif. The first field must be
          the count, and the line must end with format
          "[0:c_0,..,n:c_n,...] i,j ..." where n-1 is the total number
          of nodes and events, 'c_i' is the color of node/event i and
          the part after ']' describes all directed edges.
        node_types : dict {int: Bunch}
          A dictionary that tells which colors correspond to actual
          nodes. The key is the color and node_types[color].name is a
          description of the node type.
        event_types : dict {int: Bunch}
          A dictionary that tells which colors correspond to events
          nodes. The key is the color and event_types[color].name a
          description of the event type.
        """
        self.node_types = node_types
        self.event_types = event_types
        self.__line = line
        self.__canonical_labels = False
        self.__read_properties(line, header_line)
        # The full initialization will require running
        # self.__read_vertices(), self.__read_edges() and
        # self.__canonize(), and in this order. The class is laid out
        # so that these will be run only when needed. This allows
        # reading in the motifs quickly when some are discarded.

    def __read_properties(self, line, header_line):
        """Read additional properties from line."""
        self.__property_names = []
        if header_line:
            fields = line.split()
            headers = header_line.split()
            i = 2
            while i < len(headers) and headers[i] != "N":
                try:
                    val = int(fields[i])
                except ValueError:
                    val = float(fields[i])
                self.__dict__[headers[i].replace('-','_')] = val
                self.__property_names.append(headers[i].replace('-','_'))
                i += 1;

    def copy(self):
        """Create a copy of the motif."""

        # Note that the copy does not use the variable self.__line,
        # because the motif might have been altered after
        # initialization and these changes wouldn't show up in
        # self.__line. (And if current object is a copy, self.__line
        # won't exist anyway.)

        new_motif = Motif(None, self.node_types, self.event_types)
        new_motif.__count = self.count
        new_motif.__ref_count = self.ref_count

        for property_name in self.__property_names:
            new_motif.__dict__[property_name] = self.__dict__[property_name]
        new_motif.__property_names = self.__property_names

        # Copy event net. Make sure self.__read_edges() is called for
        # self before continuing with copying.
        new_motif.__event_net = DirNet(self.nof_events)
        try:
            edges = list(self.__event_net.edges)
        except AttributeError:
            self.__read_edges()
            edges = list(self.__event_net.edges)          
        for i,j,w in edges:
            new_motif.__event_net[i][j] = w

        # Because self.__read_edges() has been called, self.__events
        # has been fully initialized and we can go on with copying.
        new_motif.__events = self.__events.copy()
        new_motif.__nodes = self.__nodes.copy()
        
        # Other variables in new_motif will constructed when needed
        # (by an implicit call of __canonize()).
        return new_motif

    @property
    def nof_nodes(self):
        try:
            return len(self.__nodes)
        except AttributeError:
            self.__read_vertices()
            return len(self.__nodes)

    @property
    def nof_events(self):
        try:
            return len(self.__events)
        except AttributeError:
            self.__read_vertices()
            return len(self.__events)

    @property
    def nof_edges(self):
        return len(set( frozenset((x[0],x[1])) for x in self.events.itervalues() ))

    @property
    def nof_dir_edges(self):
        return len(set( (x[0],x[1]) for x in self.events.itervalues() ))

    def _get_count(self):
        try:
            return self.__count
        except AttributeError:
            self.__read_vertices()
            return self.__count
    def _set_count(self, x):
        try:
            self.__count = x
        except AttributeError:
            self.__read_vertices()
            self.__count = x
    count = property(_get_count, _set_count)

    def get_ref_count(self):
        try:
            return self.__ref_count
        except AttributeError:
            self.__read_vertices()
            return self.__ref_count
    def set_ref_count(self, x):
        try:
            self.__ref_count = x
        except AttributeError:
            self.__read_vertices()
            self.__ref_count = x
    ref_count = property(get_ref_count, set_ref_count)

    @property
    def nodes(self):
        if not self.__canonical_labels:
            self.__canonize()
        return self.__nodes

    @property
    def events(self):
        if not self.__canonical_labels:
            self.__canonize()
        return self.__events

    @property
    def event_net(self):
        if not self.__canonical_labels:
            self.__canonize()
        return self.__event_net

    @property
    def topological_hash(self):
        try:
            return self.__topological_hash
        except AttributeError:
            self.__canonize()
            return self.__topological_hash

    @property
    def node_net(self):
        try:
            return self.__node_net
        except AttributeError:
            self.__canonize()
            return self.__node_net

    @property
    def basic_type_map(self):
        """Map all node types to 0 and event types to 1."""
        return dict([(x,0) for x in self.node_types]+[(x,1) for x in self.event_types])        

    def full_init(self):
        """Ensure the motif is fully initialized."""
        if not hasattr(self, "__node_net"):
            self.__canonize()
        #self.__line = ""

    def __read_vertices(self):
        """Read nodes, events and their colors from string data.

        This method allows using the variables nof_nodes and
        nof_events.
        """
        fields = self.__line.split()
        self.__count = int(fields[0])
        self.__ref_count = float(fields[1])
        #self.__n_lt_ref = int(fields[5])

        i_field = 1
        while fields[i_field][0] != "[":
            i_field += 1

        # Read in the vertices (identities of nodes and events). The
        # node and event indices are not final. The final indices are
        # calculated only when needed.
        self.__nodes, self.__events = {}, {}
        for vertex_color_str in fields[i_field][1:-1].split(","):
            i, color = map(int, vertex_color_str.split(":"))
            if color in self.node_types:
                self.__nodes[i] = color
            elif color in self.event_types:
                self.__events[i] = [None, None, color]
            else:
                sys.stderr.write("Motif.__init__: Unidentified color %d\n"
                                 % (color,))
                exit(1)

        self.__edges_str = fields[i_field+1:]

        # If there are more nodes in this node that previously seen,
        # update the max node count and node roles for graph.
        if len(self.__nodes) > Motif.N_nodes_max:
            Motif.N_nodes_max = len(self.__nodes)
            Motif.node_roles['graph'] = range(Motif.N_nodes_max)

    def __read_edges(self):
        """Read edges from string data."""

        # Using self.nof_events makes sure self.__read_vertices() is
        # called before the rest of this method is run.
        self.__event_net = DirNet(self.nof_events)

        # Read in the edges. Also construct the directed graph induced
        # by event vertices (event_graph).
        for edge_pair_str in self.__edges_str:
            i,j = map(int, edge_pair_str.split(","))
            if i in self.__nodes:
                self.__events[j][0] = i
            elif j in self.__nodes:
                self.__events[i][1] = j
            else:
                # An edge between two event nodes.
                self.__event_net[i][j] = 1

    def __canonize(self):
        """Find canonical labels.

        Uses the current values of self.__nodes, self.__events and
        self.__event_net to calculate new indices for nodes and
        edges. The new node indices will be based on the underlying
        undirected canonical graph (of nodes), while the events will
        be relabeled from 0 to n_events-1.
        """

        # Make sure self.__read_edges() is called.
        try:
            prev_edges = list(self.__event_net.edges)
        except AttributeError:
            self.__read_edges()
            prev_edges = list(self.__event_net.edges)            

        # Map event indices to [0,N_events-1] and update event_net to
        # use the same indices.
        event_map = dict(zip(self.__events.keys(), range(self.nof_events)))
        self.__event_net = DirNet(self.nof_events)
        for i,j,w in prev_edges:
            self.__event_net[event_map[i]][event_map[j]] = w
        
        # Recalculate node indices so that similar nodes will be
        # plotted similarly. This is done by calculating the untyped
        # hashs, but needs to be done in three steps: 
        #
        #    1. Remap indices based on untyped, directed, and ordered
        #    graph. After this motifs that differ only in types will
        #    have identical indices.
        #
        #    2. Remap indices based on untyped and directed underlying
        #    graph. After this motifs that have identical underlying
        #    graph will have identical indices.
        #
        #    3. Remap indices based on untyped and undirected underlying
        #    graph. After this motifs that differ in the direction of
        #    events will have identical indices.

        # Remap based on untyped, directed and fully ordered graph.
        dg = Digraph()
        edge_set = set([(x[0],x[1]) for x in self.__events.itervalues()])
        for i_ev,(i,j,c) in enumerate(self.__events.itervalues()):
            dg.add_event(i_ev, i, j)
            if i_ev > 0:
                dg.add_event_order(i_ev-1, i_ev)
        node_map = dg.canonical_labeling()

        # Update self.__nodes and self.__events.
        self.__nodes = dict([(node_map[i],self.__nodes[i]) for i in self.__nodes])
        self.__events = dict([(event_map[i],(node_map[j],node_map[k],c)) 
                            for i,(j,k,c) in self.__events.iteritems()])

        # Recalculate node indices by using the directed, untyped
        # underlying graph. to calculate new node and event
        # indices. This is needed to make sure all motifs that differ
        # only in the order/number of events or node/event types will
        # be drawn using the same coordinates.
        dg = Digraph()
        edge_set = set([(x[0],x[1]) for x in self.__events.itervalues()])
        for i,j in edge_set:
            dg.add_edge(i,j)
        node_map = dg.canonical_labeling()

        # Update self.__nodes and self.__events.
        self.__nodes = dict([(node_map[i],self.__nodes[i]) for i in self.__nodes])
        self.__events = dict([(i,(node_map[j],node_map[k],c)) 
                              for i,(j,k,c) in self.__events.iteritems()])

        # Construct the undirected and untyped underlying graph to
        # obtain the undirected hash (this is used for example to draw
        # all triangles with the same coordinates)
        g = Graph()
        edge_set = set([(min(x[0],x[1]),max(x[0],x[1])) for x in self.__events.itervalues()])
        for i,j in edge_set:
            g.add_edge(i,j)
        node_map = g.canonical_labeling()
        self.__topological_hash = hash(g.relabel(node_map))

        # Update self.__nodes and self.__events.
        self.__nodes = dict([(node_map[i],self.__nodes[i]) for i in self.__nodes])
        self.__events = dict([(i,(node_map[j],node_map[k],c)) 
                              for i,(j,k,c) in self.__events.iteritems()])

        self.__canonical_labels = True

        # Create underlying directed graph.
        self.__node_net = DirNet()
        for i,j,c in self.__events.itervalues():
            self.__node_net[i][j] += 1

    def get_user_coords(self):
        """Calculate coordinates for nodes."""
        return visuals.calculateCoordinates(self.node_net)
        #return visuals.Himmeli(self.node_net).getCoordinates()

    def event_iter(self):
        """Iterate event ids in temporal order."""
        g = copynet(self.event_net)
        in_degs = sorted([(g[j].inDeg(), j) for j in g])
        while in_degs:
            yield in_degs[0][1]
            g.delNode(in_degs[0][1])
            in_degs = sorted([(g[j].inDeg(), j) for j in g])
    
    def get_edge_labels(self, use_colors=False):
        """Create edge labels.

        An edge label will consist of ranks for all events taking
        place on that edge, sorted in ascenting order. The elements of
        the sets defining components will be denoted by letters.

        Parameters
        ----------
        use_colors : bool
          If True, the labels of individual events will be tagged with
          a latex color. Note that this only works with the ps-backend
          in matplotlib. In addition '\usepackage{color}' must be
          defined in the LaTex preamble.
        """
        if (self.nof_events == 1):
            return {}

        # Find all directed paths between event nodes.
        g = copynet(self.event_net)
        paths_to = dict([node,set()] for node in g)
        in_degs = sorted([(g[j].inDeg(), j) for j in g])
        while in_degs and in_degs[0][0] == 0:
            node = in_degs[0][1]
            # Push this node and its set of paths to all children.
            for child in g[node].iterOut():
                paths_to[child].update(paths_to[node], set([node]))
            g.delNode(node)
            in_degs = sorted([(g[j].inDeg(), j) for j in g])
        # Now 'paths_to[i]' is a set that contains all starting nodes
        # of directed paths that lead to node 'i'.

        # Find nodes with zero in-degree, and initialize each one's
        # set with a unique element (unless there is only one root, in
        # which case an empty set is fine).
        g = copynet(self.event_net)
        event_ranks = dict([(node,[1,set()]) for node in g])
        s_max = 0
        in_degs = sorted([(g[j].inDeg(), j) for j in g])
        if len(in_degs) > 1 and in_degs[1][0] == 0:
            while in_degs and in_degs[0][0] == 0:
                event_ranks[in_degs[0][1]][1] = set([s_max])
                s_max += 1
                in_degs = in_degs[1:]
                
        for i_ in range(len(self.event_net)):
            # There must always be at least one node with zero
            # in-degree; pick any if there are several.
            in_degs = sorted([(g[j].inDeg(), j) for j in g])
            node = in_degs[0][1] # Current parent
            assert(in_degs[0][0] == 0)

            # Find those children who do not have incoming paths from
            # other children.
            children = set([child for child in g[node].iterOut()])
            good_children = filter(lambda x: set.isdisjoint(children, paths_to[x]), children)
            for child in good_children:
                event_ranks[child][0] = max(1+event_ranks[node][0],
                                            event_ranks[child][0])
                event_ranks[child][1].update(event_ranks[node][1])
                if len(good_children) > 1:
                    event_ranks[child][1].add(s_max)
                    s_max += 1
            g.delNode(node)
        # event_ranks[event_id] = [int, set()]
            
        # Event ranks done; combine the ranks of different events on
        # the same edges.
        edge_labels_vec = dict([((i,j),[]) for i,j,c in self.events.itervalues()])
        for e,(i,j,c) in self.events.iteritems():
            edge_labels_vec[(i,j)].append((event_ranks[e],c))
        # edge_labels_vec[(node_id,node_id)] = [ ([int, set()],event_color), ...]

        # Construct RGB color labels for different event types if colors are used.
        if use_colors:
            event_colors = dict()
            for c, et in self.event_types.iteritems():
                color_rgb = "%.4f,%.4f,%.4f" % matplotlib.colors.colorConverter.to_rgb(et.color)
                event_colors[c] = "\color[rgb]{%s}" % (color_rgb,)

        edge_labels_str = {}
        for k,v in edge_labels_vec.iteritems():
            edge_labels_str[k] = ""
            for ir, (r,c) in enumerate(sorted(v)):
                # r = [int, set()]
                # Get number
                event_label_str = str(r[0])
                # Get possible letter string to denote sets.
                if r[1]: event_label_str += "".join(map(lambda x:chr(97+x), r[1]))
                # Wrap in color if colored labels are used.
                if use_colors:
                    event_label_str = r"{%s%s}" % (event_colors[c],event_label_str)
                # Add comma if there are still more events on this edge.
                edge_labels_str[k] += event_label_str + ("," if len(v) > ir + 1 else "")

        return edge_labels_str

    def get_node_colors(self):
        return dict([(i,self.node_types[c].color) for i,c in self.nodes.iteritems()])

    def get_node_edge_colors(self):
        return dict([(i,self.node_types[c].edge_color) for i,c in self.nodes.iteritems()])
    
    def get_edge_colors(self):
        return dict([((j,k),self.event_types[c].color) for j,k,c in self.events.itervalues()])

    def get_label(self):
        """Return label describing the whole motif."""
        return r"$%d$" % (self.count,)

    def update(self, other_motif):
        """Update occurrence count from isomorphic motif."""
        try:
            self.__ref_count += other_motif.ref_count
        except AttributeError:
            self.__read_vertices()
            self.__ref_count += other_motif.ref_count

        self.count += other_motif.count
        

    def remap_types(self, type_map):
        """Remap the types of nodes and/or events.

        Parameters
        ----------
        type_map : {int: int}
          If you do not want to differentiate between some node/event
          types you can use this to map different types into the same
          type. If for example `node_types` has keys 0 and 1 but you
          want to treat them as equal, you could give `type_map =
          {1:0}` to map type 1 also to 0.
        """
        # Update nodes and events to use the new types.
        self.__nodes = dict([(i,type_map.get(c,c)) 
                           for i,c in self.nodes.iteritems()])
        self.__events = dict([(i,(j,k,type_map.get(c,c))) 
                            for i,(j,k,c) in self.events.iteritems()])


    def __get_directed_multigraph(self, type_map=None):
        """Return directed multigraph for calculating hashes.
    
        Arguments
        ---------
        type_map : {int: int}
          Temporarily change some node/event types when calculating
          the hash; see remap_types() for format. Cannot include other
          types than those originally used.
                  
        Returns
        -------
        dg : Digraph
          A directed graph.
        """
        type_map = (type_map or {})

        dg = Digraph()
        # Add vertices.
        for i,c in self.nodes.iteritems():
            dg.add_vertex(i, type_map.get(c,c))
        # Add events. 
        for e,(i,j,c) in self.events.iteritems():
            dg.add_event(e, i, j, type_map.get(c,c))
        return dg

    def __get_directed_graph(self, type_map=None):
        """Return directed graph for calculating hashes.
    
        Arguments
        ---------
        type_map : {int: int}
          Temporarily change some node/event types when calculating
          the hash; see remap_types() for format. Cannot include other
          types than those originally used. Note that because we use
          underlying (simple) edges instead of events, event types
          have no effect.
                  
        Returns
        -------
        dg : Digraph
          A directed graph.
        """
        type_map = (type_map or {})

        dg = Digraph()
        # Add vertices.
        for i,c in self.nodes.iteritems():
            dg.add_vertex(i, type_map.get(c,c))
        # Add edges (no types!)
        for e,(i,j,c) in self.events.iteritems():
            dg.add_edge(i, j)
        return dg

    def __get_undirected_graph(self, type_map=None):
        """Return undirected graph for calculating hashes.
    
        Arguments
        ---------
        type_map : {int: int}
          Temporarily change some node/event types when calculating
          the hash; see remap_types() for format. Cannot include other
          types than those originally used. Note that because we use
          underlying (simple) edges instead of events, event types
          have no effect.
                  
        Returns
        -------
        g : Graph
          An undirected graph.
        """
        type_map = (type_map or {})

        g = Graph()
        # Add vertices.
        for i,c in self.nodes.iteritems():
            g.add_vertex(i, type_map.get(c,c))
        # Add edges (no types).
        edge_set = set([(min(x[0],x[1]),max(x[0],x[1])) for x in self.__events.itervalues()])
        for i,j in edge_set:
            g.add_edge(i,j)
        return g

    def get_typed_hash(self, type_map=None):
        """Shortcut for get_hash('motif', type_map)."""
        return self.get_hash('motif', type_map)

    def get_directed_hash(self, type_map=None):
        """Shortcut for get_hash('multigraph', type_map)."""
        return self.get_hash('multigraph', type_map)

    def get_hash(self, graph_type, type_map=None, canonical_labels=False):
        """Return the canonical hash.
        
        Arguments
        ---------
        graph_type : str
          Defines what kind of graph is used to calculate the
          hash. The possible values are:
           'motif':      Directed multigraph with event order
           'multigraph': Directed multigraph
           'digraph':    Directed underlying graph
           'graph':      Undirected underlying graph
        type_map : {int: int}
          Temporarily change some node/event types when calculating
          the hash; see remap_types() for format. The default is to
          use the original types. For untyped hashed give
          `self.basic_type_map` as `type_map`.
        canonical_labels : bool
          If True, returns also a dictionary whose keys correspond to
          node indices in the motif and the values to node indices in
          the undirected untyped canonical graph.
          
        Returns
        -------
        hash : int
          The hash of the given graph.
        canonical_labels : {int: int}
          The canonical labels, a dictionary whose keys correspond to
          node indices in the motif and the values to node indices in
          the undirected untyped canonical graph. Returned only if
          canonical_labels is True.
        """
        if graph_type == 'motif' or graph_type == 'multigraph':
            g = self.__get_directed_multigraph(type_map)
            if graph_type == 'motif':
                # Add event order.
                for i,j,w in self.event_net.edges:
                    g.add_event_order(i,j)
        elif graph_type == 'digraph':
            g = self.__get_directed_graph(type_map)
        elif graph_type == 'graph':
            g = self.__get_undirected_graph(type_map)
        else:
            raise TypeError

        if canonical_labels:
            clabels = self.__get_undirected_graph(self.basic_type_map).canonical_labeling()
            return hash(g), clabels
        return hash(g)

    def get_node_roles(self, graph_type, latex=False):
        """Return the node roles.
        
        The values for different graph types are
           'motif':      {R,C,T,N}
           'multigraph': {R,C,T}
           'digraph':    {R,C,T}
           'graph':      {0,1,2,...}
        where the values stand for
           R:  Receiver, out-degree is zero
           C:  Caller, in-degree is zero
           T:  Transmitter, non-zero in- and out-degree
           T:  (motif) Causal transmitter, first receivers then sends
           N:  Non-causal transmitter, first sends then receives (at
               least once)

        With 'graph' the node role corresponds to the canonical labels
        of the undirected untyped graph.

        Arguments
        ---------
        graph_type : str
          Defines what kind of graph is used to calculate the
          hash. The possible values are:
           'motif':      Directed multigraph with event order
           'multigraph': Directed multigraph
           'digraph':    Directed underlying graph
           'graph':      Undirected underlying graph
        latex : bool
          If True, return the node labels in LaTeX format (but without
          the '$'-signs)'.

        Returns
        -------
        node_roles : {int: str}
          The node roles for given topology.
        """
        roles = {}
        if graph_type == 'graph':
            roles = self.__get_undirected_graph(self.basic_type_map).canonical_labeling()
        elif graph_type == 'multigraph' or graph_type == 'digraph':
            for i in self.nodes:
                if self.node_net.outDeg(i):
                    roles[i] = ("T" if self.node_net.inDeg(i) else "C")
                else:
                    roles[i] = "R"
        else: # motif
            for event_id in self.event_iter():
                i,j,c = self.events[event_id]
                # If the caller had no role, set it to "C".
                # If the previous type was "R", set it to "T".
                # Otherwise no change.
                if roles.setdefault(i,"C") == "R":
                    roles[i] = "T"
                # If the receiver had no role, set it to "R". 
                # If previous role was other than "R", set it to "N"
                if roles.setdefault(j,"R") != "R":
                    roles[j] = "N"
        #if latex:
        #    latex_map = {"Tc":"T_c", "Tn":"T_n"}
        #    return dict((i,latex_map.get(r,r)) for i,r in roles.iteritems())
        return roles

  
    def __str__(self):
        """Print motif as a string."""
        vec = (["%d:%d" % (k,v) for k,v in self.nodes.iteritems()]
               +["%d:%d" % (self.nof_nodes+k,v[2]) for k,v in self.events.iteritems()])
        s = str(self.nof_nodes+self.nof_events)+" ["+",".join(vec)+"] "
        vec = ["%d,%d %d,%d" % (v[0],self.nof_nodes+k,self.nof_nodes+k,v[1]) 
               for k,v in self.events.iteritems()]
        vec += ["%d,%d" % (self.nof_nodes+i,self.nof_nodes+j)
                for i,j,w in self.event_net.edges]
        s += " ".join(vec)
        return s
                        

def count_to_str(count):
    if isinstance(count, int):
        if count < 10000:
            return "%d" % count
    else:
        if count < 100:
            return "%.2f" % count
        elif count < 1000:
            return "%.1f" % count
        elif count < 10000:
            return "%d" % int(count)

    a = int(np.floor(np.log10(count)))
    b = np.power(10,a)
    return r"%.2f\textrm{e}%d" % (float(count)/(10**a),a)

def motif_reader(motifFile, node_types=None, event_types=None):
    if isinstance(motifFile, str):
        f = open(motifFile, 'r')
    else:
        f = motifFile

    # If the first line contains headers, use it to get names for the
    # variables.
    first_line = f.next()
    header_line = None
    try:
        tmp = float(first_line.split()[0])
        # The first line did not contain headers, so initialize the motif.
        m = Motif(first_line, node_types, event_types)
    except ValueError:
        header_line = first_line

    for line in f:
        m = Motif(line, node_types, event_types, header_line)
        yield m

    if isinstance(motifFile, str):
        f.close()


def read_motifs(fileName, node_types, event_types):
    """Read motifs from file.

    Parameters
    ----------
    fileName: string
       The file to read the motifs from. First line is skipped.
    node_types, event_types: dict
       The keys of each tell the vertex types that correspond to nodes
       and events.

    Returns
    -------
    motifs: {motif_hash: Motif object}
       The motifs addressed by their full typed hash.
    """
    # Get line count.
    lineCount = nof_lines(fileName)
    sys.stderr.write("Reading motifs (%d lines)\n" % (lineCount,))
    motifs = {}

    for line in data_utils.progress_reader(fileName):
        # Make sure we are not at header, then read motif.
        if line.split()[0].lower() == 'count':
            continue
        m = Motif(line, node_types, event_types)
        h = m.get_typed_hash()
        motifs[h] = m
    return motifs

def filter_motifs(motifs, N_events=None, N_nodes=None, type_map=None,
                  use_node_types=None, node_type_policy=None, incl_other_nodes=False,
                  use_event_types=None, event_type_policy=None, incl_other_events=False,
                  accepted_hashes=None, min_count=None, verbose=False):
    """Read all motifs in a smart way.

    Parameters
    ----------
    motifs: {motif_hash: Motif object}
       The motifs addressed by their full typed hash.
    N_events: int or set(int)
       If given, include only motifs with the given number of events.
    N_nodes: int or set(int)
       If given, include only motifs with the given number of nodes.
    type_map: {old_type: new_type}
       Map some node/event types into something else (use this to make
       some motifs appear identical, i.e. to disregard the types).
    use_node_types: set
       The node types to use. The meaning of this set is defined by the
       arguments `node_type_policy` and `other_nodes`. Note that the checking
       of types is done after mapping types with the `type_map`.
    node_type_policy : 'one', 'any', or 'all'
       'one': Include motif only if a single type from `use_node_types`
              is present.
       'any': Include motif if there is at least one node with 
              type included in `use_node_types`.
       'all': Include motif if a node of each type given in `use_node_types`
              is present.
    incl_other_nodes : bool
       If True, the motif can also include other types than those in 
       `use_node_types`. If `use_node_types` is not given or empty, this
       parameter has no effect.
    accepted_hashes : {str, set(int)}
       Sets of hashes of motifs that are accepted. The possible keys
       in the dict are "motif", "multigraph", "digraph" or "graph",
       and the corresponding value gives those hashes that are
       accepted. This can be used to filter motifs by the topology.
       Use 'construct_hash' to easily create a hash to compare with.
    min_count : int
       The minimum count for motif. Motifs with count strictly smaller
       than `min_count` will be excluded.
    verbose : bool
       If True, print info on the number of motifs processed.

    The selection of event types is identical to node types.
       
    Returns
    -------
    valid_motifs: {motif_hash: Motif object}
       The motifs that pass the filter.
    motifs_of_type: {type: set of motif_hash}
       The keys will be all types in 'use_node_types' and
       'use_event_types', and the value is a set of all motif
       hashes where all nodes/events are of given type.
    """
    type_map = (type_map or {})
    use_node_types = (use_node_types or set())
    use_event_types = (use_event_types or set())
    N_events = (set([N_events]) if type(N_events) == int else N_events)
    N_nodes = (set([N_nodes]) if type(N_nodes) == int else N_nodes)
    node_type_policy = (node_type_policy or 'any')
    event_type_policy = (event_type_policy or 'any')
    if node_type_policy not in ('one','any','all'):
        sys.stderr.write("Error: Invalid node type policy '%s'.\n" % node_type_policy)
        return False
    if event_type_policy not in ('one','any','all'):
        sys.stderr.write("Error: Invalid event type policy '%s'.\n" % event_type_policy)
        return False
        
    valid_motifs = {}
    motifs_of_type = dict([(x,set()) for x in set.union(use_node_types, 
                                                        use_event_types)])

    # Get line count.
    n_motifs = len(motifs)
    if verbose: sys.stderr.write("Filtering %d motifs ..." % (n_motifs,))
    
    for h,m in motifs.iteritems():
        # Create a motif object and check that it has the correct
        # number of events.
        if N_events and m.nof_events not in N_events:
            #sys.stderr.write("Failed at event number check! (%d)\n" % m.nof_events)
            continue
        
        if N_nodes and m.nof_nodes not in N_nodes:
            #sys.stderr.write("Failed at node number check! (%d)\n" % m.nof_nodes)
            continue
        
        # Check count.
        if min_count and m.count < min_count:
            continue
        
        # Check directed hash. Note that we need to go through the
        # hashes starting from the highest level ('graph'). A match at
        # the higher level is enough.
        if accepted_hashes:
            accepted = False
            for graph_type in ("graph", "digraph", "multigraph", "motif"):
                if accepted_hashes.get(graph_type):
                    if m.get_hash(graph_type, m.basic_type_map) in accepted_hashes[graph_type]:
                        #sys.stderr.write("Failed at topology check.\n")
                        accepted = True
                        break
            if not accepted:
                continue

        # Map types before checking them. Because this alters the
        # motif object, we take a copy of the motif first. 
        m_cp = m.copy()
        if type_map:
            m_cp.remap_types(type_map)
            h = m_cp.get_typed_hash()
                                    
        def check_types(curr_types, use_types, type_policy, incl_other):
            """Make sure the node/event types are ok."""
            if not use_types: return True
            
            # Either all types are ok OR those included are ok.
            if not (incl_other or curr_types.issubset(use_types)):
                return False
                
            if type_policy == 'one':
                # Only one of the used types is in the motif.
                return len(curr_types.intersection(use_types)) == 1
            elif type_policy == 'any':
                # At least one of the used types is present.
                return not curr_types.isdisjoint(use_types)
            else:
                # All of the types must be in the motif.
                return use_types.issubset(curr_types)
                
        # Check that node and event types fill the requirements.
        if use_node_types:
            curr_node_types =  set(m_cp.nodes.itervalues())
            if not check_types(curr_node_types, use_node_types, node_type_policy, incl_other_nodes):
                #sys.stderr.write("Failed at node type check.\n")
                continue
        if use_event_types:
            curr_event_types = set(map(operator.itemgetter(2), m_cp.events.itervalues()))
            if not check_types(curr_event_types, use_event_types, event_type_policy, incl_other_events):
                #sys.stderr.write("Failed at event type check.\n")
                continue

        # This motif is good. Take it.
        if h in valid_motifs:
            valid_motifs[h].update(m_cp)
        else:
            valid_motifs[h] = m_cp
            if use_node_types:
                for nt in curr_node_types:
                    if nt in motifs_of_type:
                        motifs_of_type[nt].add(h)
            if use_event_types:
                for et in curr_event_types:
                    if et in motifs_of_type:
                        motifs_of_type[et].add(h)

    if verbose: sys.stderr.write(" %d motifs found.\n" % (len(valid_motifs),))
    return valid_motifs, motifs_of_type

def get_motifs(fileName, node_types, event_types,
               N_events=None, N_nodes=None, type_map=None,
               use_node_types=None, node_type_policy=None, incl_other_nodes=False,
               use_event_types=None, event_type_policy=None, incl_other_events=False,
               accepted_hashes=None):
    """Read and filter motifs.
    
    The first three parameters are given to 'read_motifs', the rest to
    'filter_motifs'.
    """
    motifs = read_motifs(fileName, node_types, event_types)
    return filter_motifs(motifs, N_events, N_nodes, type_map,
                         use_node_types, node_type_policy, incl_other_nodes,
                         use_event_types, event_type_policy, incl_other_events,
                         accepted_hashes)

def init_prev_coords():
    """Coordinates for some simple directed motifs."""
    prev_coords = {}
    g = Graph()

    # Two nodes only, put them side by side.
    g.add_edge(0,1)
    h = hash(g.relabel(g.canonical_labeling()))
    prev_coords[h] = {0:(0.,0.),1:(1.,0.)}

    # Chain of two edges, put them into a triangular shape.
    g.add_edge(0,2)
    h = hash(g.relabel(g.canonical_labeling()))
    prev_coords[h] = {0:(0.,0.),1:(1.,0.), 2:(0.5,-0.6)}

    # A triangle, draw as an equilateral triangle.
    g.add_edge(1,2)
    h = hash(g.relabel(g.canonical_labeling()))
    prev_coords[h] = {0:(0.,0.),1:(1.,0.), 2:(0.5,-np.sqrt(0.75))}

    # A 3-star, draw as an equilateral triangle.
    g.del_edge(1,2)
    g.add_edge(0,3)
    h = hash(g.relabel(g.canonical_labeling()))
    prev_coords[h] = {0:(0.,0.),1:(1.,0.), 2:(0.5,np.sqrt(0.75)), 3:(0.5, np.sqrt(0.75)/3.0)}

    # A chain of 3 edges.
    g.del_edge(0,3)
    g.add_edge(2,3)
    h = hash(g.relabel(g.canonical_labeling()))
    prev_coords[h] = {0:(0.,0.), 1:(1.,0.), 2:(0.15,1.), 3:(0.85,1.)}

    # A 4-star, draw symmetric.
    g.del_edge(2,3)
    g.add_edge(0,3)
    g.add_edge(0,4)
    h = hash(g.relabel(g.canonical_labeling()))
    prev_coords[h] = {0:(0.,0.),1:(1.,0.), 2:(0.,1.), 3:(1.,1.), 4:(0.5,0.5)}

    # A 5-chain, draw as a saw-pattern.
    g.del_edge(0,4)
    g.add_edge(2,4)
    g.del_edge(0,1)
    g.add_edge(1,3)
    h = hash(g.relabel(g.canonical_labeling()))
    theta = 2.0*np.pi/5.0
    beta = np.pi/10.0
    labels = (4,2,3,0,1)
    prev_coords[h] = dict([(labels[i],(np.cos(i*theta+beta),np.sin(i*theta+beta))) for i in range(5)])
    #prev_coords[h] = {0:(0.,0.),4:(1.5,0.7),2:(1.,0.),1:(2.,0.),3:(0.5,0.7)}

    return prev_coords

# Initialize motif coordinates that can be used globally.
coords = init_prev_coords()

def plot_motif(ax, motif, prev_coords=None,
               node_color=None, node_edge_color=None, nodeShape='o',
               labelFontSize=4, orderFontSize=None,
               label_fun=None, label=None, 
               node_size=5, node_edge_width=1.0,
               edge_width=0.4, colored_edges=False, colored_edge_labels=False,
               node_labels = None, node_label_size=None,
               topology="motif"):
    """Draw a single temporal motif as a directed graph."""
    # Use global if nothing else given.
    prev_coords = (prev_coords or coords)

    # First find out the coordinates.
    user_coords = prev_coords.setdefault(motif.topological_hash, motif.get_user_coords())
    node_colors = (None if node_color else motif.get_node_colors())
    node_edge_colors = (None if node_edge_color else motif.get_node_edge_colors())
    edge_colors = (motif.get_edge_colors() if colored_edges else None)

    # Create node labels.
    node_labels = (node_labels or {})
    if node_label_size and not node_labels:
        all_labels = tuple('abcdefghijklmnopqrstuvwxyz'.upper())
        node_labels = dict(zip(motif.nodes,all_labels[:motif.nof_nodes]))

    edge_labels = motif.get_edge_labels(use_colors=colored_edge_labels)
    if topology=="multigraph":
        if colored_edge_labels:
            edge_labels = motif.get_edge_labels(use_colors=False)
        edge_labels = dict((e_id,"*"*len(s.split(","))) for e_id,s in edge_labels.iteritems())
    elif topology=='digraph' or topology=='graph':
        edge_labels={}

    net = motif.node_net
    if topology=='graph':
        edge_width=4*edge_width
        net = pynet.SymmNet()
        for i,j,w in motif.node_net.edges:
            net[i][j] = 1
        
    # Draw first the user nodes.
    margin=0.2
    visuals.visualizeNet(net, coords=user_coords,
                         scaling=True, margin=margin, axes=ax, 
                         frame=False,
                         defaultNodeShape=nodeShape,
                         defaultNodeColor=(node_color or 'w'),
                         nodeColors=node_colors,
                         defaultNodeEdgeColor=(node_edge_color or 'k'),
                         nodeEdgeColors=node_edge_colors,
                         defaultEdgeColor='k', edgeColors=edge_colors,
                         defaultNodeSize=node_size,
                         defaultNodeEdgeWidth=node_edge_width,
                         defaultEdgeWidth=edge_width,
                         edgeLabels=edge_labels,
                         edgeLabelSize=orderFontSize,
                         nodeLabels=node_labels,
                         nodeLabelSize=node_label_size,
                         defaultLabelPosition='in')

    if label is not None:
        label_str = label
    elif label_fun is not None:
        label_str = label_fun(motif)
    else:
        label_str = motif.get_label()
    
    if label_str:
        ax.text(0.5, 0.75*(margin/2.0), label_str,
                horizontalalignment='center',
                verticalalignment='top',
                transform = ax.transAxes,
                fontsize=labelFontSize)


def type_reader(f):
    """Read node or event types from file.
    
    Parameters
    ----------
    f : file    
       File that contains the types of nodes or events. The format is
           0    Description_0    ShortName_0    color    [edge_color]
           1	Description_1    ShortName_1    color    [edge_color]
           ...
       The integers that denote the types must all be distinct, i.e. the
       same integer can not be used for the type of both a node and an
       event. The descriptions cannot contain whitespace. The color given
       is used for plotting; if edge_color is not given, it
       will be the same as color.

       The reading ends when a line cannot be interpreted to represent
       a type.
       
    Returns
    -------
    (int, str, str, str)
       The values corresponds to type id, type description, color and
       edge color. The edge color corresponds to color if not given in
       the file
    """
    for line in f:
        try:
            fields = iter(line.split())
            type_id = int(fields.next())
            type_name = fields.next()
            first_char = type_name[0]
            if first_char in ('"', "'"):
                type_name = type_name[1:]
                if type_name[-1] == first_char:
                    type_name = type_name[:-1]
                else:
                    for s in fields:
                        if s[-1] == first_char:
                            type_name += " "+s[:-1]
                            break
                        type_name += " "+s
            try:
                short_name = fields.next()
                node_color = fields.next()
            except StopIteration:
                sys.stderr.write("Error: Invalid type description file format. "
                                 "At least 4 columns required.\n")
                exit(1)
            try:
                edge_color = fields.next()
            except StopIteration:
                edge_color = node_color
        except (IndexError, ValueError):
            break
        yield type_id, type_name, short_name, node_color, edge_color

def read_types(typeFile):
    """Read motif and event types.
    
    Parameters
    ----------
    typeFile: str OR file
       The file that contains the types of nodes and events. The format is
           #NODES
           0    Description_of_node_0    ShortName_0    color    [edge_color]
           1	Description_of_node_1    ShortName_0    color    [edge_color]
           #EVENTS
           2	Description_of_event_2   ShortName_0    color
           5	Description_of_event_5   ShortName_0    color
       The integers that denote the types must all be distinct, i.e. the
       same integer can not be used for the type of both a node and an
       event. The descriptions cannot contain whitespace. The color given
       is used for plotting; if edge_color is not given for nodes, it
       will be the same as color.
       
       The reading ends when a line cannot be interpreted to represent
       a type.

    Returns
    -------
    node_types: {int: Bunch}
       The key is the node type id, and the value has three
       properties: name, color and edge_color. Use as
       'node_types[ID].name', 'node_types[ID].color', etc.
    event_types: {int : Bunch}
       The key is the event type id, and the value has two
       properties: name and color.
    """
    if type(typeFile) == str:
        f = open(typeFile, 'r')
    elif hasattr(typeFile, '__iter__'):
        f = typeFile
    else:
        return None
    
    node_types, event_types = {}, {}
    first_line = f.next()
    if first_line.strip() != "#NODES":
        sys.stderr.write("Invalid type file.\n")
    else:
        for type_id, type_desc, short_name, color, edge_color in type_reader(f):
            node_types[type_id] = Bunch(name=type_desc, short_name=short_name, color=color, edge_color=edge_color)
        for type_id, type_desc, short_name, color, edge_color in type_reader(f):
            event_types[type_id] = Bunch(name=type_desc, short_name=short_name, color=color)

    if type(typeFile) == str:
        f.close()

    return node_types, event_types


def read_tmf(motifFileName, verbose=False, progress_callback=None):
    """Read a tmf file that has both types and motifs.

    Parameters
    ----------
    motifFileName: str    
       The type file that contains the types of nodes and events. The format is
           #NODES
           0    Description_of_node_0    color    [edge_color]
           1	Description_of_node_1    color    [edge_color]
           #EVENTS
           2	Description_of_event_2   color
           5	Description_of_event_5   color
           #MOTIFS
           ...
       The integers that denote the types must all be distinct, i.e. the
       same integer can not be used for the type of both a node and an
       event. The descriptions cannot contain whitespace. The color given
       is used for plotting; if edge_color is not given for nodes, it
       will be the same as color.

    Returns
    -------
    node_types: {int : (str, str, str) }
       A dictionary where the keys correspond to node types and the value is
       a tuple with node description, node color and node edge color.
    event_types: {int : (str, str) }
       A dictionary where the keys correspond to event types and the value is
       a tuple with event description and edge color.
    motifs: {motif_hash: Motif object}
       The motifs addressed by their full typed hash.
    """
    motifs = {}
    if verbose:
        f = data_utils.progress_reader(motifFileName, callback=progress_callback)
    else:
        f = open(motifFileName,'r')
    node_types, event_types = read_types(f)
    for m in motif_reader(f, node_types=node_types, event_types=event_types):
        h = m.get_typed_hash()
        motifs[h] = m
    return node_types, event_types, motifs

def save_npy(motifs, npyFileName):
    """Save motifs as npy file.

    Parameters
    ----------
    motifs : {motif_hash: Motif object}
       The motifs addressed by their full typed hash. All motifs
       must use the same node and event types!
    npyFileName : str
       The file name where the motifs will be saved.
    """
    # To save space, skim the motif objects a bit: make sure they are
    # fully initialized, then remove the input line. Also remove the
    # node and event type dictionaries from all others except the
    # first motifs.
    m = motifs.itervalues().next()
    node_types, event_types = m.node_types.copy(), m.event_types.copy()
    motifs_copy = {}
    for h,m in data_utils.progress(motifs.iteritems(), len(motifs)):
        m_cp = m.copy()
        m_cp.full_init()
        m_cp.node_types = {}
        m_cp.event_types = {}
        motifs_copy[h] = m_cp
    with open(npyFileName, 'wb') as f:
        pickle.dump((node_types, event_types, motifs_copy), f, protocol=1)

def load_npy(npyFileName):
    """Load motifs from npy file.

    Parameters
    ----------
    npyFileName : str
       The file name where the motifs will be loaded from.

    Returns
    -------
    node_types: {int : (str, str, str) }
       A dictionary where the keys correspond to node types and the value is
       a tuple with node description, node color and node edge color.
    event_types: {int : (str, str) }
       A dictionary where the keys correspond to event types and the value is
       a tuple with event description and edge color.
    motifs: {motif_hash: Motif object}
       The motifs addressed by their full typed hash.
    """
    with open(npyFileName, 'rb') as f:
        node_types, event_types, motifs = pickle.load(f)

    # Restore node and event types.
    for m in data_utils.progress(motifs.itervalues(), len(motifs)):
        m.node_types = node_types
        m.event_types = event_types

    return node_types, event_types, motifs


if __name__ == '__main__':

    # Test Motif class.
    node_types =  {0: Bunch(name='Postpaid, Female'),
                   1: Bunch(name='Prepaid, Female'),
                   3: Bunch(name='Postpaid, Male'),
                   4: Bunch(name='Prepaid, Male'),
                   6: Bunch(name='Postpaid, Unknown'),
                   7: Bunch(name='Prepaid, Unknown')}
    event_types = {2: Bunch(name='Call'),
                   5: Bunch(name='SMS')}
    type_map = {3:0,4:0}

    line_1 = "5638733 628976.50 8.9649 430.5000 2 11637.0650 4 [0:2,1:2,2:3,3:3] 0,2 1,0 1,2 3,0 3,1"
    motif_1 = Motif(line_1, node_types, event_types)
    print motif_1.nof_nodes, "== 2"
    print motif_1.nof_events, "== 2"
    print motif_1.get_edge_labels()
    print

    line_2 = "5638733 628976.50 8.9649 430.5000 2 11637.0650 4 [0:2,1:2,2:3,3:4] 0,2 0,1 1,2 3,0 3,1"
    motif_2 = Motif(line_2, node_types, event_types)
    print motif_1.topological_hash, "==", motif_2.topological_hash
    print motif_1.get_typed_hash(), "!=", motif_2.get_typed_hash()
    print motif_1.get_typed_hash(type_map), "==", motif_2.get_typed_hash(type_map)
    print motif_1.get_directed_hash(), "!=", motif_2.get_directed_hash()
    print motif_1.get_directed_hash(type_map), "==", motif_2.get_directed_hash(type_map)
    print

    motif_2 = Motif(line_2, node_types, event_types)    
    motif_2_copy = motif_2.copy()
    print motif_2.get_typed_hash(), "==", motif_2_copy.get_typed_hash()
    print

    line_3 = "5638733 628976.50 8.9649 430.5000 2 11637.0650 4 [0:2,1:2,2:2,3:3,4:3,5:4,6:4] 4,1 1,3 4,0 0,5 5,2 2,6 0,1 1,2"
    motif_3 = Motif(line_3, node_types, event_types)    
    print motif_3.nof_nodes, "== 4"
    print motif_3.nof_events, "== 3"
    print motif_3.events
    print motif_3.get_node_roles(graph_type="motif")
    print motif_3.get_node_roles(graph_type="motif")
    print motif_3.get_edge_labels()
    print motif_3.get_edge_labels()
    print 

    #for i, (mf_ev, mf_sub) in enumerate(motif_3.submotifs()):
    #    print "Remove event",i
    #    print "   |   ".join(map(str,[mf_ev,mf_sub]))
    #print

    line_4 = "5638733 628976.50 8.9649 430.5000 2 11637.0650 4 [0:2,1:2,2:3,3:4,5:2,6:4] 0,2 0,1 1,2 3,0 3,1 3,5 5,6 1,5"
    motif_4 = Motif(line_4, node_types, event_types)
    print motif_4.nof_nodes, "== 3"
    print motif_4.nof_events, "== 3"
    print motif_4.events
    
    #print "Count={count:d}".format(**(motif_1.__dict__))
    #for i, (mf_ev, mf_sub) in enumerate(motif_4.submotifs()):
    #    print "Remove event",i
    #    print "   |   ".join(map(str,[mf_ev,mf_sub]))
    #print
