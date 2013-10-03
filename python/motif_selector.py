"""Common command line interface for reading motifs.
"""

import os
import sys
import motif as mf
import operator
import argparse
import collections

class Error(Exception):
    """Base class for exceptions in this module."""
    pass

class FileError(Error):
    """Something wrong with reading input file.

    Attributes:
        expr -- input expression in which the error occurred
        msg  -- explanation of the error
    """
    def __init__(self, expr, msg):
        self.expr = expr
        self.msg = msg

def dist_median(dist):
    """Return median index of a distribution."""
    return len(dist) - (dist.cumsum() >= 0.5).sum()

class MotifSelector(object):
    """Read, filter, group and order motifs.

    Create your own argparse.ArgumentParser and supply
    MotifSelector.parser as a parent. Process arguments, then given
    the processed arguments as parameter to a new MotifSelector
    object:
        
      my_parser = argparse.ArgumentParser(parents=[MotifSelector.parser])
      arg = my_parser.parse_args()
      motif_selector = MotifSelector(arg)

    Initializing motif_selector reads in the motifs, then filters,
    groups and orders them according to command line parameters. The
    results can then be accessed using the following iterators:
    
       MotifSelector.__iter__
          Iterate through all motif objects. If the motifs have been
          grouped (using the argument '--group_by'), the argument
          specifies which group to iterate.

       MotifSelector.iter_groups
          Iterator through all groups. This simply
          MotifSelector.__iter__ at every iteration. Note that if the
          argument '--group_by' was not used, there is only one group.

      The following structures can also be accessed directly:

      motifs : {int: Motif}
        A map of motif objects. The key corresponds to full typed hash. 

      motifs_of_type : {int: [int]}
        Key is a motif type given in '--node_types' or
        '--event_types', and the value is a list of motif hashes that
        have at least one node with that type.
          
      motif_groups : [[int]]
        Each element is a single group, and each group consists of a
        list of motif hashes in that group.

      datafile : str
        The currently used motif data file.  

    MotifSelector supports multiple data files, as long as they all
    use the same parameters. Use next_data() to update motifs and
    motif_groups based on the next file; if this returns false, there
    are no more data files available.
    
    If you change the arguments in variable 'arg', call process_arg()
    to reprocess arguments or reset() to first read in the motifs from
    file and then reprocess.
    """

    typefile_desc="""
    'TYPEFILE' is a text file that contains the descriptions of the types of
    nodes and events. The format is
       #NODES
       0	User of type 0     Node_color_0     [Node_edge_color_0]
       1	User of type 1     Node_color_0     [Node_edge_color_1]
       #EVENTS
       2	Event of type 2    Edge_color_2
       5	Event of type 5    Edge_color_5
    The integers that denote the types must all be distinct, i.e. the same
    integer can not be used for the type of both a node and an
    event. Colors are any color strings understood by matplotlib.
    """
    parser = argparse.ArgumentParser(add_help=False, epilog=typefile_desc,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    datafile_help="""
        Motif distribution files. If this is a tmf or npy file (that
        contains motif and event types), the argument --typefile may
        be omitted. If however --typefile is also given, the types
        will be read from that file instead."""
    parser.add_argument('-df', '--datafile', nargs='*', help=datafile_help, required=True)
    parser.add_argument('-tf', '--typefile', help='Node and event type file.')
    parser.add_argument('-e', '--events', nargs='*', type=int, help="Number of events in the motifs.")
    parser.add_argument('-n', '--nodes', nargs='*', type=int, help="Number of nodes in the motifs.")
    node_types_help="""
        List of node types. The parameter '--type_policy' defines how
        these are treated.  Empty value corresponds to using no types
        (default)."""
    parser.add_argument('-nt', '--node_types', nargs='*', type=int, help=node_types_help)
    event_types_help="""
        List of event types. The parameter '--type_policy' defines how
        these are treated.  Empty value corresponds to using no types
        (default)."""
    parser.add_argument('-et', '--event_types', nargs='*', type=int, help=event_types_help)
    like_motif_help="""
        Output only motifs that are isomorphic to a given motif. The
        motif is given as a list of events, where each event is a
        comma separated pair of integers (nodes) and the order of
        these pairs denotes their temporal order. For example '0,1 1,2
        2,0' would match any 3-cycle. If multiple motifs are given,
        a motif must be isomorphic to one of them. If the other
        '--like' arguments are used, a motif must be isomorphic at any
        level to be accepted."""
    parser.add_argument('--like_motif', nargs='*', action='append', help=like_motif_help)
    like_digraph_help="""
        As '--like_motif', but matches the underlying directed graph
        of the motif; the order of edges in the input has no
        significance.
        """
    parser.add_argument('--like_digraph', nargs='*', action='append', help=like_digraph_help)
    like_graph_help="""
        As '--like_motif', but matches the underlying undirected graph
        of the motif; neither the order of edges in the input or the
        ordering of the pairs describing each edge has no
        significance.
        """
    parser.add_argument('--like_graph', nargs='*', action='append', help=like_graph_help)
    type_help="""
        Defines how the types should be treated. If 'one', the motif
        can have only single node type (of those given in '--types'),
        if 'any', at least one node of input types must be present,
        and if 'all' all input types must be present. Append '+' to
        allow other nodes that those in '--types' to be present."""
    parser.add_argument('-tp', '--type_policy', 
                        choices=('one','one+','any','any+','all','all+'),
                        default='any', help=type_help)
    order_help="""
        Defines the order in which the motifs are plotted. The value
        is any property available for the motif; 'rel_count' refers to
        count relative to motif's group (see parameter
        '--group_by'). Note that if '--group_by' is also used, the
        motifs are sorted inside each group, but the groups themselves
        are not ordered."""
    parser.add_argument('-o', '--order_by', type=str, default='count', help=order_help)
    parser.add_argument('-or', '--order_reversed', action='store_true', 
                        default=False, help="Order largest first.")
    group_help="""
        Defines how to group motifs by their topology. The value
        consists of two strings: [typed|untyped]
        [motif|multigraph|digraph|graph|nodes]. 'nodes' refers to the
        set of node types in the motif.
        """
    parser.add_argument('-g', '--group_by', type=str, nargs=2, default=None, help=group_help)
    mc_help="""
        The minimum count for motif. Motifs with smaller count are not read."""
    parser.add_argument('-mc', '--min_count', type=int, default=1, help=mc_help)

    npy_help="""
        If given, try to load the data from an npy file instead of the
        one given. If there is no npy file but the data was
        succesfully read from a tmf/dat file, saves the data as npy
        file."""
    parser.add_argument('--npy', action='store_true', default=False, help=npy_help)

    parser.add_argument('--verbose', action='store_true', default=False,
                        help="Print info on processing steps to stderr.")
    
    def __init__(self, arg, filter_motifs=True):
        """Initialize MotifSelector from arguments.

        Parameters
        ----------
        arg : argparse.ParseArguments
          Contains the input arguments. Note that this object must be
          one where MotifSelector.parser was used are parent.
        """
        self.arg = arg
        self.i_data = 0
        self.all_motifs = None
        self.check_data_files()
        self.__read_all_motifs()
        if filter_motifs:
            self.process_arg()

    def reset(self, filter_motifs=True):
        self.i_data = 0
        self.i_group = 0
        self.all_motifs = None
        self.__read_all_motifs()
        if filter_motifs:
            self.process_arg()

    @property
    def datafile(self):
        return self.arg.datafile[self.i_data]

    def check_data_files(self):
        """Make sure all data files exist.

        Raises a FileError exception if there is a problem with data
        files, so can also be used during initialization.
        """
        i_data_prev = self.i_data # Reset back to this in the end.
        for i in range(len(self.arg.datafile)):
            self.i_data = i
            if not self.is_valid():
                raise FileError("MotifSelector.check_data_files", "Error with input data files.")
        self.i_data = i_data_prev

    def reset(self):
        self.i_data = 0
        self.__read_all_motifs()
        self.process_arg()
    
    def has_data(self):
        return self.i_data < len(self.arg.datafile)

    def nof_files(self):
        return len(self.arg.datafile)

    def next_data(self):
        """Move to next good data file."""
        self.i_data += 1
        self.all_motifs = None
        while self.has_data():
            if self.__read_all_motifs() and self.process_arg():
                # Read succesfully, return data.
                return True
            # Couldn't read data, try next one.
            self.i_data += 1
        # Could not read any data.
        return False

    def __iter__(self, i_group=0):
        """Iterate through motifs.

        If there are multiple groups, the parameter can be used to
        select the group.
        """
        return (self.motifs[h] for h in self.motif_groups[i_group][1])

    def iter_groups(self):
        """Iterate through motif groups.

        Yields
        ------
        group_id, motif_iterator : int/str, iter(self)
           The group key is whatever value that defines this
           group. The iterator is obtained by calling iter(self).
        """
        for i_group in range(len(self.motif_groups)):
            yield self.motif_groups[i_group][0], self.__iter__(i_group)

    def iter_all(self):
        """Iterate through all motif by groups.

        Yields
        ------
        motif: Motif object
        """
        for i_group in range(len(self.motif_groups)):
            for m in self.__iter__(i_group):
                yield m

    def nof_groups(self):
        """Number of groups."""
        return len(self.motif_groups)

    def get_description(self, show_node_types=True, show_event_types=True, short=False):
        """String description of parameters."""
        policy_str = (self.arg.type_policy if self.arg.node_types else "")
        events_str = ("_E%s" % ("-".join(map(str,sorted(self.arg.events))),) if self.arg.events else "")
        nodes_str = ("_N%s" % ("-".join(map(str,sorted(self.arg.nodes))),) if self.arg.nodes else "")
        if short: return '%s%s%s' % (policy_str, events_str, nodes_str)

        node_types_str = ("_NT%s" % ("-".join(map(str,sorted(self.arg.node_types))),) 
                          if self.arg.node_types and show_node_types else "")
        event_types_str = ("_ET%s" % ("-".join(map(str,sorted(self.arg.event_types))),) 
                           if self.arg.event_types and show_event_types else "")
        return '%s%s%s%s%s' % (node_types_str, event_types_str, policy_str, events_str, nodes_str)

    def is_valid(self):
        """Make sure arguments are ok."""
        retval = True
        if not self.has_data():
            sys.stderr.write("No datafile given.\n")
            return False
        if not os.path.isfile(self.datafile):
            sys.stderr.write("File '%s' does not exist!\n" % self.datafile)
            retval = False
        if self.arg.typefile and not os.path.isfile(self.arg.typefile):
            sys.stderr.write("File '%s' does not exist!\n" % self.arg.typefile)
            retval = False
        if self.datafile[-4:] not in ('.npy','.tmf') and not self.arg.typefile:
            sys.stderr.write("Type file not given!\n")
            retval = False
        return retval

    def __read_all_motifs(self):
        """Read in all types and motifs."""

        if not self.is_valid():
            raise FileError("Invalid file.\n")

        npyFileName = self.datafile[:-4] + ".npy"
        if self.datafile[-4:] == '.npy' or (self.arg.npy and os.path.isfile(npyFileName)):
            if self.arg.verbose: sys.stderr.write("Reading motifs from file '%s'\n" % npyFileName)
            self.node_types, self.event_types, self.all_motifs = mf.load_npy(npyFileName)
        else:
            if self.arg.verbose: sys.stderr.write("Reading motifs from file '%s'\n" % self.datafile)
            if self.datafile[-4:] == '.tmf':
                self.node_types, self.event_types, self.all_motifs = mf.read_tmf(self.datafile)
            if self.arg.typefile:
                self.node_types, self.event_types = mf.read_types(self.arg.typefile)
            if self.datafile[-4:] == '.dat':
                self.all_motifs = mf.read_motifs(self.datafile, self.node_types, self.event_types)

            # Save the npy-file if the data was read from something else
            # than npy and the npy file does not yet exist.
            if self.arg.npy:
                mf.save_npy(self.all_motifs, npyFileName)
            
        return True

    def process_arg(self):
        """Filter, group and order motifs."""
        if self.arg.verbose: sys.stderr.write("Processing motifs ...\n")

        # Accepted numbers of nodes and events.
        N_events = set(self.arg.events or [])
        N_nodes = set(self.arg.nodes or [])

        # Create directed hashes if any graphs are given as input.
        accepted_hashes = {}
        if self.arg.like_motif:
            accepted_hashes['motif'] = set()
            for edges_str in self.arg.like_motif:
                edges = map(lambda s: map(int,s.split(',')), edges_str)
                accepted_hashes['motif'].add(mf.construct_hash(edges, 'motif'))
        if self.arg.like_digraph:
            accepted_hashes['digraph'] = set()
            for edges_str in self.arg.like_digraph:
                edges = map(lambda s: map(int,s.split(',')), edges_str)
                accepted_hashes['digraph'].add(mf.construct_hash(edges, 'digraph'))
        if self.arg.like_graph:
            accepted_hashes['graph'] = set()
            for edges_str in self.arg.like_graph:
                edges = map(lambda s: map(int,s.split(',')), edges_str)
                accepted_hashes['graph'].add(mf.construct_hash(edges, 'graph'))

        # Figure out what to do with node types.
        type_map = None
        if not self.arg.node_types:
            # Map everything to the smallest node type.
            common_node_type = min(self.node_types.keys())
            use_node_types=set([common_node_type])
            type_map = dict((k,common_node_type) for k in self.node_types)
        else:
            use_node_types = set(self.arg.node_types)
            if not use_node_types.issubset(set(self.node_types.keys())):
                sys.stderr.write('Error: Unidentified node types, possible types are %s\n' 
                                 % str(sorted(self.node_types.keys())))
                exit(1)

        # Figure out what to do with event types.
        if not self.arg.event_types:
            # Map everything to the smallest event type.
            common_event_type = min(self.event_types.keys())
            use_event_types=set([common_event_type])
            type_map = dict((k,common_event_type) for k in self.event_types)
        else:
            use_event_types = set(self.arg.event_types)
            if not use_event_types.issubset(set(self.event_types.keys())):
                sys.stderr.write('Error: Unidentified event types, possible types are %s\n' 
                                 % str(sorted(self.event_types.keys())))
                exit(1)

        # Check type policy.
        if self.arg.type_policy[-1] == '+':
            type_policy = self.arg.type_policy[:-1]
            incl_other_nodes = True
            incl_other_events = True
        else:
            type_policy = self.arg.type_policy
            incl_other_nodes = False
            incl_other_events = False

        # Filter all motifs with the given number of events.
        # motifs:         key = typed hash 
        #                 value = motif object
        # motifs_of_type: key = node type AND event type
        #                 value = list of typed motif hashes where the
        #                         motif has this node/event.
        self.motifs, self.motifs_of_type = mf.filter_motifs(self.all_motifs, N_events=N_events, N_nodes=N_nodes,
                                                            type_map=type_map,
                                                            use_node_types=use_node_types, 
                                                            node_type_policy=type_policy, incl_other_nodes=incl_other_nodes,
                                                            use_event_types=use_event_types,
                                                            event_type_policy=type_policy, incl_other_events=incl_other_events,
                                                            accepted_hashes=accepted_hashes,
                                                            min_count=self.arg.min_count,
                                                            verbose=self.arg.verbose)

        # Check the motif grouping argument.
        if self.arg.group_by is not None:
            # Make sure the format is correct.
            if not (len(self.arg.group_by) == 2 and
                    self.arg.group_by[0] in set(['typed','untyped']) and
                    self.arg.group_by[1] in set(['motif','multigraph','digraph','graph','nodes'])):
                sys.stderr.write("  Wrong format for '--group_by': %s\n" % self.arg.group_by)
        
        if self.arg.group_by == 'nodes':
            # Group by the set of node types in the motif.
            hash_groups = collections.defaultdict(list)
            for h,m in self.motifs.iteritems():
                hash_groups[set(m.nodes.values())].append(h)
        elif self.arg.group_by is not None:
            # Group by topology.
            hash_groups = collections.defaultdict(list)
            for h,m in self.motifs.iteritems():
                h_topo = m.get_hash(self.arg.group_by[1], (None if self.arg.group_by[0] == 'typed' else m.basic_type_map))
                hash_groups[h_topo].append(h)
        else:
            hash_groups = {"all": self.motifs.keys()}

        if len(hash_groups) > 1:
            sys.stderr.write("Motif group sizes (%s): %s\n" % (" ".join(self.arg.group_by), str(map(len, hash_groups.values()))),)

        # Create motif groups and order by the defined criterion. Ordering is done separately for each group.
        self.motif_groups = []
        for group_key, motif_hashes in hash_groups.iteritems():
            # If the motifs don't have the variable they should be ordered by, add that variable.
            if not hasattr(self.motifs.itervalues().next(), self.arg.order_by):
                if self.arg.order_by == 'median':
                    for h in motif_hashes:
                        self.motifs[h].median = dist_median(self.motifs[h].distribution)
                elif self.arg.order_by == 'rel_count':
                    total_count = float(sum(self.motifs[h].count for h in motif_hashes))
                    for h in motif_hashes:
                        self.motifs[h].rel_count = self.motifs[h].count/total_count
                elif self.arg.order_by == 'ratio':
                    for h in motif_hashes:
                        self.motifs[h].ratio = (float(self.motifs[h].count)/self.motifs[h].ref_count if self.motifs[h].ref_count > 0 else 0)
                        
            # Sort motifs and take only the motif hash.
            motifs_tmp = sorted([(getattr(self.motifs[h],self.arg.order_by),h) for h in motif_hashes],
                                key=operator.itemgetter(0), reverse=self.arg.order_reversed)
            self.motif_groups.append((group_key, map(operator.itemgetter(1), motifs_tmp)))

        # Motifs are now ready: filtered, grouped and ordered.
        return True

