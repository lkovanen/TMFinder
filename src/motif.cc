/* Motif is just a bliss::digraph object. The only thing added here is
 *  a method for printing it nicely. 
 */

#include "motif.h"

std::ostream& operator<<(std::ostream& output, const Motif& m) 
{
  output << m.get_nof_vertices() << " [0:" << m.vertices[0].color;
  for (unsigned int i = 1; i < m.get_nof_vertices(); ++i)
    output << "," << i << ":" << m.vertices[i].color;
  output << "] ";
  for (unsigned int i = 0; i < m.get_nof_vertices(); ++i)
    {
      for (unsigned int j = 0; j < m.vertices[i].nof_edges_out(); ++j)
	{
	  output << i << "," << m.vertices[i].edges_out[j] << " ";
	}
    }
  return output;
}
