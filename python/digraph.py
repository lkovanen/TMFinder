"""Directed graph for calculating canonical labels.

PyBliss only supports undirected graphs. To calculate canonical form
of a directed graphs we therefore map it to an undirected graph by
splitting each node into two and giving them different colors to
denote whether their edges correspond to incoming or outgoing edges.
"""

from PyBliss import Graph

class Enumerator(object):
    def __init__(self, init_val):
        self.val = init_val
    def __call__(self):
        self.val += 1
        return self.val-1

class Digraph(object):
    """Directed graph for calculating canonical forms.
    
    Digraph can have either edges or events; the difference is that
    only one edge can exist between two nodes, whereas multiple events
    are possible. The events are also ordered in time.

    For any Digraph use either `add_edge` or `add_event`, do not mix
    the two. If events are used, add information about their mutual
    order with `add_event_order`.
    """
    def __init__(self):
        """Initialize an empty digraph."""
        self.g = Graph()
        self.c_n = 0 # Color used internally for nodes.
        self.c_e = 1 # Color used internally for edges.
        self.node_map = {}
        self.edge_map = {}
        self.event_map = {}
        self.enum = Enumerator(0)

    def add_vertex(self, v, color=0):
        """Add a new vertex with given id and color.

        Does nothing if a vertex with the given id already
        exists. Note that missing vertices are added automatically
        when new edges or events are added.
        """
        if v in self.node_map:
            return False
        n_0, n_1 = self.enum(), self.enum()
        self.g.add_vertex(n_0, self.c_n)
        self.g.add_vertex(n_1, color+2)
        self.g.add_edge(n_0, n_1)
        self.node_map[v] = (n_0,n_1,color)
        return True

    def __add_link(self, v1, v2, color):
        self.add_vertex(v1)
        self.add_vertex(v2)
        e_0, e_1 = self.enum(), self.enum()
        self.g.add_vertex(e_0, self.c_e)
        self.g.add_vertex(e_1, color+2)
        self.g.add_edge(e_0, e_1)
        self.g.add_edge(self.node_map[v1][1], e_1)
        self.g.add_edge(e_1, self.node_map[v2][0])
        return e_0, e_1

    def add_edge(self, v1, v2, color=1):
        """Add an edge with a given color.

        There can be only one edge per ordered node pair. Trying to
        create an edge that already exists will do nothing.
        """
        if (v1,v2) in self.edge_map:
            return
        e_0, e_1 = self.__add_link(v1,v2,color)
        self.edge_map[(v1,v2)] = (e_0,e_1,color)

    def add_event(self, e_id, v1, v2, color=1):
        """Add an event with a given id and color.

        There can be multiple events between any pair of nodes.
        """
        e_0, e_1 = self.__add_link(v1,v2,color)
        self.event_map[e_id] = (e_0,e_1,color)

    def add_event_order(self, e_id_1, e_id_2):
        """State that event `e_id_1` takes place before `e_id_2`."""
        self.g.add_edge(self.event_map[e_id_1][1], self.event_map[e_id_2][0])

    def get_vertices(self):
        """Returns a list of node ids."""
        return self.node_map.keys()

    def get_edge(self):
        """Returns a list of edge ids."""
        return self.edge_map.keys()

    def get_events(self):
        """Returns a list of event ids."""
        return self.event_map.keys()

    def canonical_labeling(self):
        """Create canonical labeling of nodes.

        The returned dictionary has keys corresponding to the node
        labels given as input when constructing the Digraph and values
        corresponding to canonical labeling [0 .. n-1] that is
        guaranteed to be identical for isomorphic digraphs or
        multigraphs.

        NB! Do not use this labeling to calculate the hash, use
        hash(Digraph) instead.
        """
        # The implementation is slighty more complicated than just
        # using the canonical labeling of the Graph object because
        # there are many other vertices than just those corresponding
        # to nodes. However, by starting from the canonical labeling
        # of the Graph object we can be sure the returned labels are
        # the same for any isomorphic graphs.
        labels = self.g.canonical_labeling()
        real_nodes = dict((v[0],k) for k,v in self.node_map.iteritems())
        inv_labels = {}
        for n,v in labels.iteritems():
            if n in real_nodes:
                inv_labels[v] = real_nodes[n]
        new_labels = {}
        i = 0
        for v in sorted(inv_labels.keys()):
            new_labels[inv_labels[v]] = i
            i += 1
        return new_labels

    def __hash__(self):
        """Return canonical hash."""
        return hash(self.g.relabel(self.g.canonical_labeling()))

if __name__ == '__main__':

    g1 = Digraph()
    g1.add_edge(2,0)
    g1.add_vertex(1,100)
    g1.add_edge(0,1,77)
    g1.add_edge(2,1)
    print g1.canonical_labeling()
    print hash(g1)

    g2 = Digraph()
    g2.add_vertex(1,100)
    g2.add_edge(22,1)
    g2.add_edge(22,17)
    g2.add_edge(17,1,77)
    print g2.canonical_labeling()
    print hash(g2)

    assert(hash(g1) == hash(g2))
