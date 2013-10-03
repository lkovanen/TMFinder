# Different command for making figures.

import sys
import os
import pylab
import matplotlib.ticker as ticker
import numpy as np

def savefig(fig, name, extensions=None, verbose=False, dpi=300):
    """Save figure.

    Save matplotlib.figure object `fig` as `name`.EXT, where EXT are
    given in `extensions`. If only one save type is used, the full
    name (including the extension) can also be given as `name`.

    Note! Saving as 'svg' requires that the program 'pdf2svg' is
    installed on your system.
    """

    if len(name) == 0:
        raise ValueError("File name can not be empty.")

    if extensions is None:
        fields = name.split(".")
        if len(fields) == 1:
            raise ValueError("File name must contain an extension if"
                             " extensions are not given explicitely.")
        extensions = fields[-1]
        name = ".".join(fields[:-1])

    if isinstance(extensions, str):
        extensions = (extensions,)

    if verbose:
        if len(extensions) == 1:
            sys.stderr.write("Saving figure %s.%s ...\n" % (name,extensions[0]))
        else:
            sys.stderr.write("Saving figures:\n")
            for extension in extensions:
                sys.stderr.write("              %s.%s ...\n" % (name,extension))
    
    # Check if pdf should be generated (both eps and svg will be
    # created from pdf) and generate if necessary.
    pdf_generated = False
    pdf_tmp = "%s_tmp_%d.pdf" % (name, np.random.randint(100000))
    if set(['pdf','svg','eps']).intersection(extensions):
        fig.savefig(pdf_tmp, dpi=dpi)
        pdf_generated = True

    for ext in extensions:
        if not isinstance(ext, str):
            raise ValueError("'extensions' must be a list of strings.")
        if ext[0] == '.':
            ext = ext[1:]
            
        if ext == 'eps':
            pipe = os.popen("pdftops -eps %s %s.eps" % (pdf_tmp, name))
            exit_status = pipe.close()
            if exit_status:
                if os.WEXITSTATUS(exit_status) == 127:
                    sys.stderr.write("%s could not be created because program "
                                     "'pdftoeps' could not be found.\n" % (ext,))
                else:
                    sys.stderr.write("Problem saving '%s'.\n" % (ext,))
        elif ext == 'svg':
            pipe = os.popen("pdf2svg %s %s.svg" % (pdf_tmp, name))
            exit_status = pipe.close()
            if exit_status:
                if os.WEXITSTATUS(exit_status) == 127:
                    sys.stderr.write("%s could not be created because program "
                                     "'pdf2svg' could not be found.\n" % (ext,))
                else:
                    sys.stderr.write("Problem saving '%s'.\n" % (ext,))
        elif ext != 'pdf':
            # fig.savefig raises a ValueError if the extension is not identified.
            #
            # According to "http://stackoverflow.com/questions/4581504/how-to-set-
            # opacity-of-background-colour-of-graph-wit-matplotlib" it is necessary
            # to set the face- and edgecolor again!.
            fig.savefig(name + "." + ext, dpi=dpi)

    if pdf_generated:
        if 'pdf' in extensions:
            os.popen("mv %s %s.pdf" % (pdf_tmp,name))
        else:
            os.popen("rm %s" % (pdf_tmp,))

def get_rcParams(fig_width_cm, fig_ratio=0.8, font_sizes=None):
    """Set good parameters for LaTeX-figures.

    The idea idea is to set the figure width in centimeters to be the
    same as the final size in your LaTeX document. This way the font
    sizes will be correct also.

    Parameters
    ----------
    fig_width_cm: int or float
        The width of the final figure in centimeters.
    fig_ratio: float (between 0 and 1)
        The ratio height/width. < 1 is landscape, 1.0 is square and
        > 1.0 is portrait.
    font_sizes: dictionary
        The font sizes used in the figure. Default is size 10 for the
        title and 8 for everything else. Possible keys are 'default',
        'label', 'title', 'text', 'legend' and 'tick'.  'default' is
        used when the specific value is not defined, other keys should
        be self explanatory.
    """
    default_font_sizes = {'label':8, 'title':10, 'text':8, 'legend':8, 'tick':8}
    font_sizes = (font_sizes or {})
    for k in default_font_sizes:
        if k not in font_sizes:
            font_sizes[k] = font_sizes.get('default', default_font_sizes[k])

    latex_preamble = [r"\usepackage{amsmath}",
                      r"\usepackage{color}"]

    inches_per_cm = 1/2.54
    fig_width = 1.0*fig_width_cm*inches_per_cm  # width in inches
    fig_height = 1.0*fig_width*fig_ratio        # height in inches
    fig_size =  [fig_width,fig_height]
    params = {'font.family':'serif',
              'font.serif':'Computer Modern Roman',
              'axes.labelsize': font_sizes['label'],
              'axes.titlesize': font_sizes['title'],
              'text.fontsize': font_sizes['text'],
              'font.size': font_sizes['text'],
              'legend.fontsize': font_sizes['legend'],
              'xtick.labelsize': font_sizes['tick'],
              'ytick.labelsize': font_sizes['tick'],
              'text.usetex': True,
              'text.latex.preamble': latex_preamble,
              'figure.figsize': fig_size,
              'legend.labelspacing': 0.0,
              'lines.markersize': 3,
              'lines.linewidth': 0.5}
    return params


def pretty_number(x):
    """Turn integer into a pretty string.
    
    Make a pretty number by adding extra spaces between every three
    digits. The LaTeX environment must be active.
    
    Parameters
    ----------
    x : int
       The integer to print.

    Returns
    -------
    s : str
       String representation of the integer.
    """

    if isinstance(x,float):
        return r"$%.2f$" % x

    s_tmp = str(abs(x))
    s = ""
    while s_tmp:
        if len(s_tmp) > 3:
            s = r"\,"+s_tmp[-3:]+s
            s_tmp = s_tmp[:-3]
        else:
            s = ("-" if x<0 else "")+s_tmp+s
            break
    return r"$%s$" % s


class SubplotHelper(object):
    """Helps in iterating through the indices of subplot."""

    def __init__(self, rows, columns):
        self.N_v = rows
        self.N_h = columns
        self.i_max = self.N_v * self.N_h
        self.__i = 1

    def copy(self):
        ax_copy = SubplotHelper(self.N_v, self.N_h)
        ax_copy.set_index_at(self.row, self.col)
        return ax_copy

    @property
    def i(self):
        """Current index."""
        return self.__i

    def add_subplot(self, fig):
        return fig.add_subplot(self.N_v, self.N_h, self.__i)

    @property
    def col(self):
        """"Current column (1 based)."""
        return 1 + ((self.__i-1)%self.N_h)

    @property
    def row(self):
        """"Current row (1 based)."""
        return 1 + ((self.__i-1)/self.N_h)

    def set_index_at(self, row=None, col=None):
        """Set index corresponding to given row and column."""
        row = (row or self.row)
        col = (col or self.col)
        self.__i = self.N_h*(row-1)+col

    def move(self, down, right):
        """Move the index right `right` steps and down `down` steps.
        Negative values move left and up. Returns False if the
        operation would move index out of bounds.""" 
        new_col = self.col + right
        new_row = self.row + down
        if new_col < 1 or new_col > self.N_h or new_row < 1 or new_row > self.N_v:
            return False
        self.set_index_at(new_row, new_col)
        return True

    def down(self, n=1):
        """Shortcut for `move(n,0)`."""
        return self.move(n,0)

    def up(self, n=1):
        """Shortcut for `move(-n,0)`."""
        return self.move(-n,0)

    def right(self, n=1, loop=False):
        """Shortcut for `move(0,n)`."""
        ret_val = self.move(0,n)
        if loop and not ret_val:
            i_new = self.__i + n
            if i_new <= self.i_max:
                self.__i = i_new
                ret_val = True
        return ret_val

    def left(self, n=1):
        """Shortcut for `move(0,-n)`."""
        return self.move(0,-n)

    def next_row(self, forced=False):
        """Move the index at the beginning of the next row. If
        `forced` is False, the index is not moved if it is already at
        the beginning of a row. Returns False if this would move index
        out of bounds (index is not updated)."""
        if not forced and self.col == 1:
            return True
        if self.row == self.N_v:
            return False
        self.set_index_at(self.row+1, 1)
        return True

    def next_column(self, forced=False):
        """Move the index at the beginning of the next column. If
        `forced` is False, the index is not moved if it is already at
        the beginning of a row. Returns False if this would move index
        out of bounds (index is not updated)."""
        if not forced and self.row == 1:
            return True
        if self.col == self.N_h:
            return False
        self.set_index_at(1, self.col+1)
        return True
