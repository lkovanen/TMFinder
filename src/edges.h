/* Data structures for dealing with edges.
 */

#ifndef EDGES_H
#define EDGES_H

#include <vector>
#include <set>
#include <map>
#include "events.h"

struct Edge
{ 
  node_id from;
  node_id to;
  short int type;
  Edge();
  Edge(const Event& event);
  Edge(node_id from, node_id to, short int type);
  friend bool operator<(const Edge& e1, const Edge& e2);
};
typedef std::vector<Edge> EdgeVector;
typedef std::set<Edge> EdgeSet;
typedef std::map<EdgeVector, unsigned int> EdgeVectorMap;
//typedef std::map<unsigned int, EdgeVectorMap> LocationMap;

bool operator<(const Edge& e1, const Edge& e2);
std::ostream& operator<<(std::ostream& output, const Edge& e);

typedef std::vector<std::pair<node_id, node_id> > NodepairVector;

// Quick hack for creating edge vector.
bool create_edges(EdgeVector& edges,
		  const NodepairVector& nodePairs, 
		  const std::vector<short int>& eventTypes);

#endif
