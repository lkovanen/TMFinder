
/* 
 * A collection of codes for calculating bin limits based on edge weights.
 */

#include <map>
#include <set>
#include <list>
#include "Nets.H"
#include "std_printers.h"

typedef DirNet<unsigned int> NetType;

struct Interval
{ 
  unsigned int left;
  unsigned int right;
  unsigned int count;
};
typedef std::list<Interval> LimitSeq;
std::ostream& operator<<(std::ostream& stream, const Interval& ival);

/* Get the weight frequency for the directed edges of a network. The
 * frequencies are saved into `weights` so that `weights[w]` gives the
 * number of edges with weight `w`. The function also returns the
 * number of directed edges with non-zero weights.
 */
unsigned int weight_dist(const NetType& net, std::map<unsigned int, unsigned int>& weights);

/* Other auxiliary functions. */
unsigned int initialize_limit_seq(LimitSeq& limit_seq, const NetType& net, double p_leave_out);
bool bin_limits_from_limit_seq(const LimitSeq& limit_seq, std::set<unsigned int>& limits);

/* Create bin limits for edge weights so that each bin has at least
 * min_count data points.
 */
bool get_limits_by_min_count(const NetType& net,
			     std::set<unsigned int>& limits,
			     unsigned int bin_min_count,
			     double p_leave_out);

/* Create bin limits so that there are N_bins bins, each of with
 * approximately the same number of data points, and exclude the
 * largest p_leave_out fraction of weights.
 */
bool get_limits_by_N_bin(const NetType& net,
			 std::set<unsigned int>& limits,
			 unsigned int N_bins,
			 double p_leave_out);

/* Create unbinned bin limits (each weight has its own bin).
 */
bool get_limits_unbinned(const NetType& net,
			 std::set<unsigned int>& limits,
			 double p_leave_out);
