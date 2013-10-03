"""Motif topology chart."""

import motif as mf
import fig_utils
import pylab
import operator

def plot_topology_chart(motifs, chartFileName, plot_single_motifs=False):
    """Create topology chart of motifs.

    Parameters
    ----------
    motifs : iterable
       The typed motifs for which the chart is produced. The order is
       not important (the names of the motifs are independent of the
       order).
    chartFileName : str
       File name where the chart is saved. The suffix should be some
       image file format. If empty, no plot is created.
    plot_single_motifs : bool
       If True, every untyped motif is plotted in a separate figure.

    Returns
    -------
    labels : {int, str}
       A dictionary where the key is the hash of the untyped motif and
       the value gives a textual name for the motif that is shown in
       the motif chart.
    
    """

    def get_label(h, v_label, labels, lvl):
        v_label[lvl] += 1 
        v_label[lvl+1:] = [0]*(3-lvl)
        label = "-".join(map(str, v_label[:lvl+1]))
        labels[h] = label
        return label

    graph_data = {}
    hash_types = ('graph','digraph','multigraph','motif')
    # First find out the topology tree. We go through hashes in
    # increasing detail: first undirected graph, then digraph, then
    # multigraph, and finally motif that includes also the timings.
    for m in motifs:
        data = graph_data
        for hash_type in hash_types:
            h = m.get_hash(hash_type, type_map=m.basic_type_map)#, canonical_labels=True)
            clabels = m.get_node_roles(hash_type, latex=True)
            if h not in data:
                m_cp = m.copy()
                m_cp.remap_types(m.basic_type_map)
                node_labels = dict((i,r"$%s$" % str(r)) for i,r in clabels.iteritems())
                data[h] = (m_cp,node_labels,{})
            data = data[h][-1]

    # Get max size of a row.
    N_v, N_h = 0, 0
    for m, tmp, digraph_data in graph_data.itervalues():
        #N_v += 1
        for m, tmp, multigraph_data in digraph_data.itervalues():
            #N_v += 1
            for m, tmp, motif_data in multigraph_data.itervalues():
                N_v += 1
                N_h = max(N_h, 4 + len(motif_data))

    # Physical size of the output.
    w_cm_motif, h_cm_motif = 2.05, 1.85 # width and height of motif in cm
    t_cm, b_cm = 1.0, 0.5 # top and bottom margin in cm
    l_cm, r_cm = 0.30, 0.30 # left and right margin in cm
    ws_cm, hs_cm = 0.13, 0.36 # wspace and hspace in cm

    # Derived physical sizes of the figure.
    w_cm = l_cm+r_cm+N_h*w_cm_motif+(N_h-1)*ws_cm
    h_cm = t_cm+b_cm+N_v*h_cm_motif+(N_v-1)*hs_cm 

    # Set plotting properties
    font_sizes = {'text': 14, 'title':13}
    params = fig_utils.get_rcParams(w_cm, h_cm/w_cm, font_sizes)
    pylab.rcParams.update(params)
    subplot_specs = {'left':l_cm/w_cm, 'right':1-r_cm/w_cm,
                     'top':1-t_cm/h_cm, 'bottom':b_cm/h_cm,
                     'wspace':ws_cm/w_cm_motif, 'hspace':hs_cm/h_cm_motif}
    labelFontSize=9
    orderFontSize=7
    node_size=10
    node_label_size=0.8*node_size

    # Create figure and subplot helper.
    fig = pylab.figure()
    fig.subplots_adjust(**subplot_specs)
    prev_coords = mf.init_prev_coords()
    i_ax = fig_utils.SubplotHelper(N_v,N_h)

    # Go through the data depth first, plotting the motifs in each branch.
    v_label = [0,0,0,0]
    labels = {} # To map untyped hashes to their labels.

    # Sort graphs by the number of nodes and edges.
    graph_hashes = [(m.nof_nodes, m.nof_edges, h) for h,(m,tmp1,tmp2) in graph_data.iteritems()]
    for h in map(operator.itemgetter(-1), sorted(graph_hashes)):
        m, node_labels, digraph_data = graph_data[h]
        # Update label
        label = get_label(h, v_label, labels, lvl=0)
        # Plot motif m as graph.
        ax = i_ax.add_subplot(fig)
        mf.plot_motif(ax, m, prev_coords, node_color="white", node_edge_color="black",
                      label=label, labelFontSize=labelFontSize, orderFontSize=orderFontSize,
                      topology="graph", node_size=node_size,
                      node_labels=node_labels, node_label_size=node_label_size)

        if ax.is_first_row():
            ax.set_title("Graph")

        # Sort digraphs by the number of directed edges (they all have
        # the same number of nodes because the undirected underlying
        # graph is the same).
        digraph_hashes = [(m.nof_dir_edges, h) for h,(m,tmp1,tmp2) in digraph_data.iteritems()]
        for h in map(operator.itemgetter(-1), sorted(digraph_hashes)):
            m, node_labels, multigraph_data = digraph_data[h]
            i_ax.set_index_at(col=2)
            label = get_label(h, v_label, labels, lvl=1)
            # Plot motif m as digraph.
            ax = i_ax.add_subplot(fig)
            mf.plot_motif(ax, m, prev_coords, node_color="white", node_edge_color="black",
                          label=label, labelFontSize=labelFontSize, orderFontSize=orderFontSize,
                          topology="digraph", node_size=node_size,
                          node_labels=node_labels, node_label_size=node_label_size)

            if ax.is_first_row():
                ax.set_title("Digraph")

            # Sort digraphs by the total number of events.
            multigraph_hashes = [(m.nof_events, h) for h,(m,tmp1,tmp2) in multigraph_data.iteritems()]
            for h in map(operator.itemgetter(-1), sorted(multigraph_hashes)):
                m, node_labels, motif_data = multigraph_data[h]
                i_ax.set_index_at(col=3)
                label = get_label(h, v_label, labels, lvl=2)
                # Plot motif m as multigraph.
                ax = i_ax.add_subplot(fig)
                mf.plot_motif(ax, m, prev_coords, node_color="white", node_edge_color="black",
                              label=label, labelFontSize=labelFontSize, orderFontSize=orderFontSize,
                              topology="multigraph", node_size=node_size,
                              node_labels=node_labels, node_label_size=node_label_size)

                if ax.is_first_row():
                    ax.set_title("Multigraph")

                i_ax.right()

                # Sort motifs by their hash.
                for h in sorted(motif_data):
                    m, node_labels, tmp = motif_data[h]
                    i_ax.right()
                    label = get_label(h, v_label, labels, lvl=3)
                    # Plot motif m as motif.
                    ax = i_ax.add_subplot(fig)
                    mf.plot_motif(ax, m, prev_coords, node_color="white", node_edge_color="black",
                                  label=label, labelFontSize=labelFontSize, orderFontSize=orderFontSize,
                                  topology="motif", node_size=node_size,
                                  node_labels=node_labels, node_label_size=node_label_size)
                    
                    if ax.is_first_row() and i_ax.col==5:
                        ax.set_title("Motifs ...")

                    if plot_single_motifs and chartFileName:
                        # Plot the motif alone in a separate figure.
                        params = fig_utils.get_rcParams(w_cm_motif, h_cm_motif/w_cm_motif, font_sizes)
                        pylab.rcParams.update(params)
                        fig_motif = pylab.figure()
                        mf.plot_motif(fig_motif.add_axes([0.,0.,1.,1.]), m, prev_coords, node_color="white", node_edge_color="black",
                                      label="",orderFontSize=orderFontSize,
                                      topology="motif", node_size=node_size,
                                      node_labels=node_labels, node_label_size=node_label_size)
                        motifFileName = "%s/motif_%s.pdf" % ("/".join(chartFileName.split("/")[:-1]), label)
                        fig_utils.savefig(fig_motif, motifFileName, verbose=True)

                i_ax.next_row()

    if chartFileName:
        fig_utils.savefig(fig, chartFileName, verbose=True)
    return labels



if __name__ == "__main__":
    pass
