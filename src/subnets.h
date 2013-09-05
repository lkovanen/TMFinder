/*
Iterate through ordered, connected sets of edges.

Lauri Kovanen (2/2012)
*/

#ifndef SUBNETS_H
#define SUBNETS_H

#include <stack>
#include <list>
#include <set>
#include <map>
#include "Nets.H"
#include "std_printers.h"
#include "edges.h"

typedef unsigned int node_id;
typedef DirNet<unsigned int> NetType;

class SubnetIterator
{
private:
  typedef std::pair<node_id, node_id> NodePair;

  const NetType& net;
  unsigned int N_max;
  std::stack<std::set<NodePair> > neighbors;
  std::list<std::set<NodePair>::const_iterator> nit;
  unsigned int curr_node;
  NodepairVector edges;

  /* Initialize neighbors with all out-edges of current node. */
  void fill_first_level();

  /* Go deeper: add all neighboring edges of the current edge. */
  void fill_next_level();

  /* Fill until full. */
  void fill_up();

public:
  SubnetIterator(const NetType& net, unsigned int N_max);
  SubnetIterator(const SubnetIterator& sgIt);

  // Reset state.
  void reset();

  // Advance state.
  SubnetIterator& operator++();

  // True if all subnets have been processed.
  inline bool finished() {return (curr_node >= net.size() && neighbors.size() == 0);}

  // Return a vector of edges at current state.
  NodepairVector& operator*();

  // Return a vector of edges at current state.
  void print_state() { 
    std::cerr << "Stack: " << neighbors << std::endl;
    std::cerr << "Value: " << **this << std::endl;
  };

};

#endif
