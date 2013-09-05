/* Representation for temporal subraph and temporal motif.
 *
 * Lauri Kovanen, BECS (August 2011)
 */
#include <queue>
#include <iostream>
#include "events.h"
#include "tsubgraph.h"

TSubgraph::TSubgraph(Events const& events,
		     const EventSet& eventSet,
		     std::vector<unsigned short int> const& node_types,
		     unsigned int dt_max)
  :  edgeVector(eventSet.size()),
     node_types(node_types),
     motif_typed(NULL),
     motif_untyped(NULL),
     motif_static(NULL),
     nodeSet(),
     edgeSet(),
     __dt_max(dt_max),
     __is_valid(false)
{
  // Construct edgeVector and eventTypes.
  __is_valid = check_validity(events, eventSet);
}

TSubgraph::TSubgraph(const EdgeVector& edgeVector,
		     const std::vector<unsigned short int>& node_types)
  :  edgeVector(edgeVector),
     node_types(node_types),
     motif_typed(NULL),
     motif_untyped(NULL),
     motif_static(NULL),
     nodeSet(),
     edgeSet(),
     __dt_max(0),
     __is_valid(true)
{ }


TSubgraph::TSubgraph(TSubgraph const& sg)
  : edgeVector(sg.edgeVector),
    node_types(sg.node_types),
    motif_typed(NULL),
    motif_untyped(NULL),
    motif_static(NULL),
    nodeSet(),
    edgeSet(),
    __dt_max(sg.__dt_max),
    __is_valid(sg.__is_valid)
{}


unsigned int TSubgraph::add_event_to_motif(Motif& g, node_id fr, node_id to, int prev_ge_id,
					   std::map<node_id, unsigned int>& nodeMap,
					   unsigned int e_type, unsigned int fr_type, unsigned int to_type,
					   bool is_static) const
{
  // Add new node vertices if necessary, collecting the node labels.
  if (nodeMap.find(fr) == nodeMap.end()) nodeMap[fr] = g.add_vertex(fr_type);
  if (nodeMap.find(to) == nodeMap.end()) nodeMap[to] = g.add_vertex(to_type);
  // Add new event vertex.
  unsigned int ge_id = g.add_vertex(e_type);
  // Add edges between the node vertices and event vertex.
  g.add_edge(nodeMap[fr], ge_id);
  g.add_edge(ge_id, nodeMap[to]);
  // Add edge from previous event vertex if there is one (and we are
  // not constructing a static motif).
  if (prev_ge_id && !is_static) g.add_edge(prev_ge_id, ge_id);
  return ge_id;
}

/* Construct motif and return it. */
Motif* TSubgraph::get_motif(bool use_node_types, bool use_event_types, bool is_static) const
{
  // Create a map of nodes and a temporary graph.
  Motif g;
  g.set_splitting_heuristic(Motif::shs_f);
  g.set_component_recursion(false);
  std::map<node_id, unsigned int> nodeMap;
  unsigned int ge_id = 0;
	
  /* Build a directed graph that includes information about the
   * temporal order of all events (irregardless of whether they are
   * adjacent).
   */
  //std::cerr << edgeVector << " " << use_node_types << " " << use_event_types << " " << is_static << std::endl;
  for (EdgeVector::const_iterator it = edgeVector.begin(); it != edgeVector.end(); ++it)
    {
      // Add edge.
      ge_id = add_event_to_motif(g, it->from, it->to, ge_id, nodeMap, 
				 (use_event_types ? it->type : 1),
				 (use_node_types ? node_types[it->from] : 0), 
				 (use_node_types ? node_types[it->to] : 0),
				 is_static);
    }
	
  /* Calculate the canonical form. Note that permute() reserves memory
     for a new Motif object. */
  bliss::Stats stats;
  return static_cast<Motif*>(g.permute(g.canonical_form(stats,NULL,NULL)));
}

Motif* TSubgraph::get_static_motif() const
{
  if (motif_static == NULL) motif_static = get_motif(false,false,true);
  return motif_static;
}

Motif* TSubgraph::get_untyped_motif() const
{
  if (motif_untyped == NULL) motif_untyped = get_motif(false,false,false);
  return motif_untyped;
}

Motif* TSubgraph::get_typed_motif() const
{
  if (motif_typed == NULL) motif_typed = get_motif(true,true,false);
  return motif_typed;
}

/* Make sure a tw-connected subgraph is valid by making sure that the
 * tw-connected events of each node must be consecutive. If yes, the
 * return value is the canonical motif and the value dt_max will be
 * the largest time difference between two consecutive events of a any
 * node in the motif (i.e. the smallest time window with which this
 * motif is valid).
 */
bool TSubgraph::check_validity(const Events& events, const EventSet& eventSet)
{
  // prev_events[v] = iterator to the previous event of node in 'event_neighbors'
  std::map<node_id, node_iterator> prev_events;
  std::map<node_id, node_iterator>::iterator it;

  /* We go through the events of each node in temporal order. If two
   * consecutive events in the subgraph are also consecutive in the
   * full data, then everything is ok: whether they are tw-adjacent or
   * not makes no difference for the validity of the
   * subgraph. However, if they are not consecutive in the full data,
   * then the subgraph is valid only if the events in between are not
   * in the same maximal subgraph.
   */
  unsigned int i_ev = 0;
  for (EventSet::const_iterator event_iter = eventSet.begin(); 
       event_iter != eventSet.end(); ++event_iter)
    {
      event_id curr = *event_iter;
      Event const& e = events[curr];
		
      node_id node = e.from();
      for (int i_node = 0; i_node < 2; i_node++)
	{
	  it = prev_events.find(node);
	  if (it == prev_events.end())
	    {
	      prev_events[node] = events.find_node_event(node, curr);
	    }   
	  else
	    {
	      node_iterator uit(it->second);
	      ++uit;
				
	      while (*uit != curr)
		{
		  // At least one event in between. Make sure they are
		  // not in the same maximal subgraph. Note that if
		  // components have not been calculated this check
		  // will _always_ fail, and hence there can be no
		  // events in between (even if they in reality were
		  // in another component).
		  if (events[*uit].component() == events[curr].component()) return false;
		  ++uit;
		}
                    
	      // Update previous event of this node.
	      prev_events[node] = uit;
	    }
	  node = e.to(); // Change to the other node and repeat.
	}

      // Check ok so far; update edgeVector and eventTypes.
      edgeVector[i_ev] = Edge(e);
      ++i_ev;
    }
  return true;
}

void TSubgraph::create_node_and_edge_sets() const
{
  nodeSet.clear();
  edgeSet.clear();
  for (EdgeVector::const_iterator it = edgeVector.begin(); it != edgeVector.end(); ++it)
    {
      unsigned int v1 = it->to;
      unsigned int v2 = it->from;
      if (it->from < it->to)
	{
	  v1 = it->from;
	  v2 = it->to;
	}
      nodeSet.insert(v1);
      nodeSet.insert(v2);
      edgeSet.insert(Edge(v1,v2,it->type));
    }
}

TSubgraph::~TSubgraph()
{
  //std::cerr << "Deleting motif at " << motif << std::endl;
  if (motif_typed) delete motif_typed;
  if (motif_untyped) delete motif_untyped;
  if (motif_static) delete motif_static;
}



TSubgraphFinder::TSubgraphFinder(event_id root_event_id,
				 unsigned int time_window,
				 unsigned int max_subgraph_size,
				 Events const& events,
				 std::vector<unsigned short int> const& node_types)
 :root_event_id(root_event_id),
  tw(time_window),
  max_subgraph_size(max_subgraph_size),
  events(events),
  node_types(node_types),
  subgraphs() {}


void TSubgraphFinder::add_subgraph(const EventSet& eventSet, unsigned int dt_max)
{
  //std::cerr << "      Adding subgraph " << subgraphs.size() << ": " << eventSet << " (dt_max = " << dt_max << ")\n";
  subgraphs.push_back(new TSubgraph(events, eventSet, node_types, dt_max));
}

void TSubgraphFinder::create_subgraphs(const EventSet& eventSet,
				      const EventSet& excludedEvents,
				      const EventMMap& validNeighbors,
				      unsigned int dt_max)
{
  /*
  std::cerr << "    eventSet       : " << eventSet << std::endl;
  std::cerr << "    excludedEvents : " << excludedEvents << std::endl;
  std::cerr << "    validNeighbors : " << validNeighbors << std::endl;
  */

  // Exclude all events that were already excluded (this includes the
  // current event). This set will be updated with new events as they
  // are added.
  EventSet excludedEvents_new(excludedEvents);

  // Go through the valid neighbors in the order of smallest time
  // difference. Note that it is possible that an event is in
  // validNeighbors twice if it can be reached via two different
  // routes. The faster route will be used, and the event is added to
  // excludedEvents, so this does not cause problems.
  EventMMap::const_iterator it;
  for (it = validNeighbors.begin(); it != validNeighbors.end(); ++it)
    {
      // Make sure the event hasn't been added yet, and if not, add
      // the current event to the list of excluded events.
      if (excludedEvents_new.find(it->second) != excludedEvents_new.end()) continue;
      excludedEvents_new.insert(it->second);
      
      // Copy the previous event set, add the new event and construct
      // the corresponding motif.
      EventSet eventSet_new(eventSet);
      eventSet_new.insert(it->second);
      dt_max = (it->first > dt_max ? it->first : dt_max);
      add_subgraph(eventSet_new, dt_max);
	      
      // Continue recursion if the maximum subgraph size has not been reached.
      if (max_subgraph_size == 0 || eventSet_new.size() < max_subgraph_size)
	{
	  // Valid neighbors are immediate neighbors of the current
	  // event that either were valid before and have a time
	  // difference larger (or possibly equal) than for the
	  // current event,
	  EventMMap::const_iterator it2(it);
	  EventMMap validNeighbors_new(++it2, validNeighbors.end());
	  // ... or are valid neighbors of the newly added event, take
	  // place after the root event, have not been excluded so far
	  // and have a time difference smaller than the time window:
	  EventMMap potentialNeighbors;
	  events.prev_immediate_events(it->second, potentialNeighbors);
	  events.next_immediate_events(it->second, potentialNeighbors);
	  for (EventMMap::const_iterator pnit = potentialNeighbors.begin();
	       pnit != potentialNeighbors.end(); ++pnit)
	    {
	      if ((pnit->second > root_event_id) && 
		  (excludedEvents_new.find(pnit->second) == excludedEvents_new.end()) &&
		  pnit->first <= tw)
		{
		  validNeighbors_new.insert(*pnit);
		}
	    }

	  // Continue recursion.
	  create_subgraphs(eventSet_new, excludedEvents_new, validNeighbors_new, dt_max);
	}
    }	
}

void TSubgraphFinder::find_subgraphs()
{
  // Events in the subgraph; initially only root event.
  EventSet eventSet;
  eventSet.insert(root_event_id);

  // Construct the motif that consists of only this one event.
  add_subgraph(eventSet,0);
  
  // Put the root event to excluded events so it is no longer added.
  EventSet excludedEvents;
  excludedEvents.insert(root_event_id);

  // Place the next events of both nodes that occur within tw into
  // validneighbors (not previous events, because we only return those
  // subgraphs where the root event is the first event.)
  
  // First get all valid neighbors ...
  EventMMap validNeighbors;
  events.next_immediate_events(root_event_id, validNeighbors);
  // ... and then cut out the part where the time differences are too big.
  EventMMap::iterator it = validNeighbors.begin();
  for (; it != validNeighbors.end(); ++it) if (it->first > tw) break;
  if (it != validNeighbors.end()) validNeighbors.erase(it, validNeighbors.end());

  // Create subgraphs recursively.
  create_subgraphs(eventSet, excludedEvents, validNeighbors, 0);

}

TSubgraphFinder::~TSubgraphFinder()
{
  for (TSubgraphList::iterator it = subgraphs.begin(); it != subgraphs.end(); ++it) delete *it;
}
