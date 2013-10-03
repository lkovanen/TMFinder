"""Plot motifs grouped by untyped topology.
"""
import sys
import pylab
import fig_utils
import numpy as np
import motif as mf
from motif_chart import plot_topology_chart
import argparse
from motif_selector import MotifSelector
import collections
import operator
import scipy.misc
import re

class LabelFun(object):
    """Call to get motif label from motif."""

    def __init__(self, label_def, sample_motif=None):
        """Construct function for motif labels."""
        self.label_def = label_def

    def __call__(self, m):
        motif_attributes = dict((name, getattr(m, name)) for name in dir(m) if not name.startswith('_'))
        return self.label_def.format(**motif_attributes)


def plot_single_motif(ax, m, label_fun, type_label_dict=None):

    node_labels = None
    if type_label_dict:
        node_labels = {}
        for i,t in m.nodes.iteritems():
            node_labels[i] = type_label_dict[t]

    mf.plot_motif(ax, m, 
                  label_fun=label_fun,
                  labelFontSize=9, 
                  orderFontSize=7,
                  node_size=(10 if node_labels else 6),
                  node_labels=node_labels,
                  node_label_size=(7 if node_labels else None),
                  colored_edges=True)


def create_motif_figure(N_h, N_v):
    """Create a figure for plotting motifs.

    `event_types` is needed to construct colors used in the labels of
    events.
    """
    # Physical size of the output.
    #t_cm_scatter, l_cm_scatter = 1.0, 1.2 # scatter left margin in cm
    #h_cm_scatter, w_cm_scatter = 8.0, 8.0 # scatter plot width in cm
    w_cm, h_cm = 2.05, 1.85 # width and height of motif in cm
    t_cm, b_cm = 0.3, 0.6   # top and bottom margin in cm
    l_cm, r_cm = 0.30, 0.30 # left and right margin in cm
    ws_cm, hs_cm = 0.16, 0.36 # wspace and hspace in cm

    # Derived physical sizes of the figure.
    w_cm = l_cm+r_cm+N_h*w_cm+(N_h-1)*ws_cm
    h_cm = t_cm+b_cm+N_v*h_cm+(N_v-1)*hs_cm 

    # Set plotting properties
    font_sizes = {'text': 14, 'title':14}
    params = fig_utils.get_rcParams(w_cm, h_cm/w_cm, font_sizes)
    pylab.rcParams.update(params)
    subplot_specs = {'left':l_cm/w_cm, 'right':1-r_cm/w_cm,
                     'top':1-t_cm/h_cm, 'bottom':b_cm/h_cm,
                     'wspace':ws_cm/w_cm, 'hspace':hs_cm/h_cm}
    #motif_colors = {2: "#990033", 3: "#0000ff", 4:"#C850BE", 5:"#64C55F"}
    
    fig = pylab.figure()
    fig.subplots_adjust(**subplot_specs)
    i_ax = fig_utils.SubplotHelper(N_v,N_h)
    return fig, i_ax

def draw_motifs(fig, i_ax, motifs, label_fun, N_max=None):
    """Draw as many motifs as fits in the figure."""

    #label_fun = lambda x: r"%.3f" % x.ratio
    for i, m in enumerate(motifs):
        ax = i_ax.add_subplot(fig)
        plot_single_motif(ax, m, label_fun)
        if (N_max and i == N_max-1) or not i_ax.right(n=1,loop=True):
            break

def plot_most_common_motifs(motifs, motifsFileName, N_h, N_v, label_fun):
    sys.stderr.write("Plotting most common motifs.\n")
    
    # Order motifs.
    #motif_sort = sorted([(m.ratio, m) for m in motifs], reverse=True)
    #ratios_sorted, motifs_sorted = zip(*motif_sort)
    motifs = list(motifs)

    if N_v is None:
        N_v = int(np.ceil(float(len(motifs))/N_h))

    if len(motifs) < N_h*(N_v-1):
        # All rows are not needed.
        N_v = int(np.ceil(len(motifs)/float(N_h)))
    if N_v == 1:
        N_h = len(motifs)

    # Draw the most common motifs.
    fig, i_ax = create_motif_figure(N_h, N_v)
    draw_motifs(fig, i_ax, motifs, label_fun)
    fig_utils.savefig(fig, motifsFileName, verbose=True)


def plot_most_common_groups(motif_group_iter, N_groups, motifsFileName,
                            N_h, N_v, label_fun):
    sys.stderr.write("Plotting most common motifs by group.\n")

    motif_groups = {}
    for group_id, motif_iter in motif_group_iter:
        motif_groups[group_id] = list(motif_iter)

    # Number of motif rows and columns to plot in each group.
    N_h = min(N_h, max(map(len, motif_groups.itervalues())))
    N_v_all = [int(np.ceil(float(len(motifs))/N_h)) for motifs in motif_groups.values()]
    if N_v is not None:
        N_v_all = min(N_v_all, [N_v]*len(motif_groups))
    N_v_tot = sum(N_v_all)

    # Draw the most common motifs in each group.
    fig, i_ax = create_motif_figure(N_h, N_v_tot)
    for i, motifs in enumerate(motif_groups.itervalues()):
        i_ax.next_row() # New group, so start from the beginning of a row.
        draw_motifs(fig, i_ax, motifs, label_fun, N_max=N_h*N_v_all[i])
    fig_utils.savefig(fig, motifsFileName, verbose=True)


if __name__ == '__main__':

    # Read input arguments.
    parser = argparse.ArgumentParser(description=__doc__, parents=[MotifSelector.parser])
    parser.add_argument('--chart', default=None, help="Name for the motif chart.")
    parser.add_argument('-pre', '--prefix', default=None, help="Output file name prefix.", required=True)
    #parser.add_argument('-suf', '--output_suffix', default="", help="Output file name suffix.")
    parser.add_argument('-Nh', type=int, default=8, help="Maximum number of columns in the plot.")
    Nv_help="""Maximum number of rows in the motif plot. If omitted,
    all motifs are plotted."""
    parser.add_argument('-Nv', type=int, default=None, help=Nv_help)
    parser.add_argument('--most_common', action='store_true', default=False,
                        help="Plot most common motifs in a single figure.")
    parser.add_argument('--most_common_grouped', action='store_true', default=False,
                        help="Plot most common motifs in groups in a single figure.")
    parser.add_argument('--most_common_grouped_single', action='store_true', default=False,
                        help="Plot most common motifs in groups, each group in a single figure.")
    label_help="""Motif label format. For example
        '{count:d}/{ref_count:.2f}' would output motif count followed
        by '/'-character and the reference count with two
        decimals. The default value is '{ratio:.4f}' which plots the
        ratio with four decimals."""
    parser.add_argument('-lb', '--label', type=str, default='{ratio:.4f}', help=label_help)

    arg = parser.parse_args()
    if not (arg.most_common or arg.most_common_grouped or arg.most_common_grouped_single):
        sys.stderr.write("Nothing selected to be plotted.\n")
        sys.exit(0)

    #print arg
    motif_selector = MotifSelector(arg, False)

    # Get function for labeling motifs.
    label_fun = LabelFun(arg.label)

    # Remove types that are not used.
    node_types = motif_selector.node_types.copy()
    event_types = motif_selector.event_types.copy()
    used_node_types = set(arg.node_types or [])
    for t in list(node_types.keys()):
        if t not in used_node_types:
            node_types.pop(t)

    # Filter motifs.
    motif_selector.process_arg()

    # Create and plot untyped motif chart.
    #suffix = ("_"+arg.output_suffix if arg.output_suffix else "")
    #chartFileName = "%s_chart.pdf" % (arg.prefix, motif_selector.get_description(short=True), suffix)
    topo_labels = plot_topology_chart(motif_selector.motifs.itervalues(), arg.chart, plot_single_motifs=False)

    fileNameBase = arg.prefix

    # All motifs sorted by ratio.
    if arg.most_common:
        plot_most_common_motifs(motif_selector.iter_all(), 
                                fileNameBase+"_all.pdf",
                                arg.Nh, arg.Nv, label_fun)

    # Grouped by topology and sorted by ratio (single figure)
    if arg.most_common_grouped:
        if (arg.group_by is None):
            sys.stderr.write("Error: Argument '--group_by' must be given "
                             "when plotting motifs in groups.\n")
            exit(0)
        plot_most_common_groups(motif_selector.iter_groups(), motif_selector.nof_groups(),
                                fileNameBase+"_grouped.pdf", arg.Nh, arg.Nv, label_fun)

    # Grouped by topology and sorted by ratio (multiple figures)
    if arg.most_common_grouped_single:
        for h, motif_iter in motif_selector.iter_groups():
            fileName = '%s_M%s.pdf' % (fileNameBase, topo_labels[h])
            plot_most_common_motifs(motif_iter, fileName, arg.Nh, arg.Nv, label_fun)
