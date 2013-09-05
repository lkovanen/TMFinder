#include "bin_limits.h"

std::ostream& operator<<(std::ostream& stream, const Interval& ival)
{
  stream << "[" << ival.left << "," << ival.right << "]:" << ival.count;
  return stream;
}

/* Get the weight frequency for the directed edges of a network. The
 * frequencies are saved into `weights` so that `weights[w]` gives the
 * number of edges with weight `w`. The function also returns the
 * number of directed edges with non-zero weights.
 *
 */
unsigned int weight_dist(const NetType& net, std::map<unsigned int, unsigned int>& weights)
{
  weights.clear();
  unsigned int total_edges = 0;
  for (unsigned int i = 0; i < net.size(); ++i)
    {
      for (NetType::const_edge_iterator j = net(i).begin(); !j.finished(); ++j) 
        {
	  unsigned int w = j.value().out();
	  if (w) {
	    weights[w]++;
	    total_edges++;
	  }
        }
    }
  //std::cout << "Weight distribution: " << weights << std::endl;
  //std::cout << "   Total number of edges: " << total_edges << std::endl;
  return total_edges;
}

unsigned int initialize_limit_seq(LimitSeq& limit_seq, const NetType& net, double p_leave_out)
{
  // Get the total weight distribution.
  std::map<unsigned int, unsigned int> weights;
  unsigned int total_edges = weight_dist(net, weights);
  unsigned int n_leave_out = (unsigned int)(p_leave_out*total_edges);

  // We leave out at least n_leave_out edges (could be more).
  std::map<unsigned int, unsigned int>::const_reverse_iterator rit = weights.rbegin();
  unsigned int cum_edges = 0;
  for (; rit != weights.rend(); ++rit)
    {
      if (cum_edges >= n_leave_out) break;
      cum_edges += rit->second;
    }

  std::cout << "   Highest weight edges to omit: " << cum_edges 
	    << " (those with w > " << rit->first << ")" << std::endl;

  // Create limit sequence where each weight has its own bin.
  limit_seq.clear();
  for (; rit != weights.rend(); ++rit) 
    {
      Interval ival;
      ival.left = rit->first;
      ival.right = rit->first;
      ival.count = rit->second;
      limit_seq.push_front(ival);
    }
  return total_edges;
}

bool bin_limits_from_limit_seq(const LimitSeq& limit_seq, std::set<unsigned int>& limits)
{
  limits.clear();
  LimitSeq::const_iterator lcit = limit_seq.begin();
  for (; lcit != limit_seq.end(); ++lcit) limits.insert(lcit->left);
  limits.insert(limit_seq.rbegin()->right + 1);
  return true;
}

// Create unbinned bin limits (each weight has its own bin).
bool get_limits_unbinned(const NetType& net, std::set<unsigned int>& limits, double p_leave_out)
{
  // Initialize limit sequence, then construct the bin limits.
  LimitSeq limit_seq;
  initialize_limit_seq(limit_seq, net, p_leave_out);
  bin_limits_from_limit_seq(limit_seq, limits);

  //std::cerr << "   Weight bin limits: " << limits << std::endl;
  return true;
}
