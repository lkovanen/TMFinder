#include "subnets.h"

void SubnetIterator::fill_first_level()
{
  // Initialize neighbors with all out-edges of current node.
  std::set<NodePair> curr_neighbors;
  while (curr_neighbors.empty() && curr_node < net.size())
    {
      for (NetType::const_edge_iterator j = net(curr_node).begin(); !j.finished(); ++j) 
	{
	  if (j.value().out()) curr_neighbors.insert(std::make_pair(curr_node,*j));
	}
      curr_node++; // Go to next node.
    }
  if (!curr_neighbors.empty())
    {
      neighbors.push(curr_neighbors);
      nit.push_back(neighbors.top().begin());
    }
}

void SubnetIterator::fill_next_level()
{
  node_id node = nit.back()->first;
  std::set<NodePair> curr_neighbors(neighbors.top());
  for (int i_iter = 0; i_iter < 2; i_iter++)
    {
      for (NetType::const_edge_iterator j = net(node).begin(); !j.finished(); ++j) 
	{
	  if (j.value().out()) curr_neighbors.insert(std::make_pair(node,*j));
	  if (j.value().in()) curr_neighbors.insert(std::make_pair(*j,node));
	}
      node = nit.back()->second; // Change to other node and repeat.
    }
  neighbors.push(curr_neighbors);
  nit.push_back(neighbors.top().begin());
}

void SubnetIterator::fill_up()
{
  if (neighbors.empty()) fill_first_level(); 
  if (!neighbors.empty())
    {
      while (neighbors.size() < N_max) fill_next_level();
    }
}

SubnetIterator::SubnetIterator(const NetType& net, unsigned int N_max) : net(net),
									 N_max(N_max),
									 neighbors(),
									 nit(),
									 curr_node(0),
									 edges() 
{ fill_up(); }

SubnetIterator::SubnetIterator(const SubnetIterator& sgIt) : net(sgIt.net),
							     N_max(sgIt.N_max),
							     neighbors(sgIt.neighbors),
							     nit(sgIt.nit),
							     curr_node(sgIt.curr_node),
							     edges()
{ std::cerr << "Copy!\n" << std::endl; }

void SubnetIterator::reset()
{ 
  while (!neighbors.empty()) neighbors.pop();
  nit.clear();
  curr_node = 0;
  edges.clear();
  fill_up();
}

SubnetIterator& SubnetIterator::operator++()
{
  // Next state.
  if (!nit.empty()) 
    {
      nit.back()++;

      // If we are at the end of current set, back out as far as
      // needed and call it done. The current state will be some
      // shorter edge list (if it is empty, we need to fill it again).
      if (nit.back() == neighbors.top().end())
	{
	  neighbors.pop();
	  nit.pop_back();
	}
      else if (neighbors.size() < N_max) fill_next_level();
    }
  
  if (nit.empty() && curr_node < net.size()) fill_up();

  return *this;
}

NodepairVector& SubnetIterator::operator*()
{ 
  edges.clear();
  std::list<std::set<NodePair>::const_iterator>::const_iterator it;
  for (it = nit.begin(); it != nit.end(); ++it) edges.push_back(**it);
  return edges;
}
