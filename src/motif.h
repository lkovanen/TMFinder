/* Simple extension to the bliss digraph class to allow pretty
 * printing.
 *
 * Lauri Kovanen, 2/2012
 */

#include <ostream>
#include "graph.hh"

#ifndef MOTIF_H
#define MOTIF_H

class Motif : public bliss::Digraph
{
  /**
   * Print the graph on one line. Format is
   *    N [0:c_0,1:c_1,...] i,j k,l ...
   * where N is the number of vertices, the next part lists the vertex
   * colors and the last part gives the edges (only out-edges are
   * listed). (Added by Lauri Kovanen, BECS 6/2010).
   */
  friend std::ostream& operator<<(std::ostream& output, const Motif& g);
};

#endif
