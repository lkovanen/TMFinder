#include "edges.h"

Edge::Edge():
  from(0),
  to(0),
  type(-1)
{}

Edge::Edge(const Event& event):
  from(event.from()),
  to(event.to()),
  type(event.type())
{}

Edge::Edge(node_id from, node_id to, short int type):
  from(from),
  to(to),
  type(type)
{}

bool operator<(const Edge& e1, const Edge& e2)
{
  if (e1.from == e2.from)
    {
      if (e1.to == e2.to) return e1.type < e2.type;
      return e1.to < e2.to;
    }
  return e1.from < e2.from; 
}

std::ostream& operator<<(std::ostream& output, const Edge& e)
{
  output << e.from << "->" << e.to << ":" << e.type;
  return output;
}

bool create_edges(EdgeVector& edges, const NodepairVector& nodePairs, const std::vector<short int>& eventTypes)
{
  for (unsigned int i = 0; i < nodePairs.size(); ++i)
    {
      edges.push_back(Edge(nodePairs[i].first, nodePairs[i].second, eventTypes[i]));
    }
  return true;
}
