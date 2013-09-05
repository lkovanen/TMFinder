/* Representation for temporal subraph and temporal motif.
 *
 * Lauri Kovanen, BECS (August 2011)
 */

#ifndef TSUBGRAPH_H
#define TSUBGRAPH_H

#include <assert.h>
#include <string>
#include <iostream>
#include <set>
#include <map>
#include <queue>
#include "motif.h"
#include "edges.h"

class Events;

typedef std::set<event_id> EventSet;
typedef std::set<node_id> NodeSet;

/* Temporal subgraph consists of an ordered sequence of edges. The
   edges may repeat (multiple events on the same edge)

 */
class TSubgraph
{
 private:
  EdgeVector edgeVector;   // Edge sequence for this motif.
  const std::vector<unsigned short int>& node_types; // Reference to the vector of node types.

  /* Different motifs. These are saved because they are owned by this
     class and need to be deleted when this object is deleted. */
  mutable Motif *motif_typed;   // Fully typed motif (uses both node types and event types)
  mutable Motif *motif_untyped; // Event type is always 1 and node type is 0.
  mutable Motif *motif_static;  // Types are as in motif_untyped, disregard order of the events.

  mutable NodeSet nodeSet; // The set of distinct nodes in this subgraph.
  mutable EdgeSet edgeSet;  // The set of distinct undirected edges in this subgraph.

  unsigned int __dt_max;   // Max time gap in this subgraph.
  bool __is_valid; // True after passing validity check.

  /* Methods for constructing the set of nodes and edges. */
  void create_node_and_edge_sets() const;

  // Make sure the subgraph is valid. Also constructs edgeVector.
  bool check_validity(const Events& events, const EventSet& eventSet);

  // Auxiliary method for constructing motifs.
  unsigned int add_event_to_motif(Motif& g, node_id fr, node_id to, int prev_ge_id,
				  std::map<node_id, unsigned int>& nodeMap,
				  unsigned int e_type, unsigned int fr_type, unsigned int to_type,
				  bool is_static) const;
  
 public:

  /* Construct the subgraph from a set of events. If this constructor
     is used, the validity of the subgraph is also checked (i.e. that
     the nodes do not have events between those given in the event
     set).
  */
  TSubgraph(const Events& events,
	    const EventSet& eventSet,
	    const std::vector<unsigned short int>& node_types,
	    unsigned int dt_max);

  /* Construct the subgraph from a sequence of edges. The subgraph is
     always valid. */
  TSubgraph(const EdgeVector& edgeVector,
	    const std::vector<unsigned short int>& node_types);

  // Copy constructor.
  TSubgraph(TSubgraph const& sg);

  /* Delete pointers to motifs. */
  ~TSubgraph();

  // Return motif.
  Motif* get_motif(bool use_node_types, bool use_event_types, bool is_static) const;
  // Shortcuts:
  Motif* get_typed_motif() const;
  Motif* get_untyped_motif() const;
  Motif* get_static_motif() const;

  /* Get the number of events in the most recently returned
     motif. */
  inline unsigned int nof_events() const { return edgeVector.size(); };
  inline unsigned int nof_edges() const { if (edgeSet.empty()) create_node_and_edge_sets(); return edgeSet.size(); };
  inline unsigned int nof_nodes() const { if (nodeSet.empty()) create_node_and_edge_sets(); return nodeSet.size(); };

  /* Methods for iterating through edge vector. */
  typedef EdgeVector::const_iterator events_iterator;
  inline events_iterator begin() const { return edgeVector.begin(); };
  inline events_iterator end() const { return edgeVector.end(); };

  /* Methods for iterating through nodes. */
  typedef NodeSet::iterator nodes_iterator;
  inline nodes_iterator nbegin() { if (nodeSet.empty()) create_node_and_edge_sets(); return nodeSet.begin(); };
  inline nodes_iterator nend() { return nodeSet.end(); };
	
  /* Methods for iterating through edge set. */
  typedef EdgeSet::iterator edges_iterator;
  inline edges_iterator ebegin() { if (edgeSet.empty()) create_node_and_edge_sets(); return edgeSet.begin(); };
  inline edges_iterator eend() { return edgeSet.end(); };

  /* Check if this is a valid subgraph. */
  inline bool is_valid() const { return __is_valid; };

  inline unsigned int dt_max() const { return __dt_max; };

  /* Build a set of event types. */
  //inline void build_event_type_set(std::set<short int>& evt) const { evt.insert(eventTypes.begin(), eventTypes.end());};
    
};

typedef std::list<TSubgraph*> TSubgraphList;

class TSubgraphFinder
{
 private:
  /* The event where the search is started. */
  const event_id root_event_id;
	
  /* Time window. */
  const unsigned int tw;
	
  /* The maximum number of events in submotifs. If this is set
     to 0, all submotifs will be searched (which might take a
     while).*/
  unsigned int max_subgraph_size;
	
  /* Reference to the events object where this motif is to be
     found. */
  const Events& events;
	
  /* Reference to the vector of node types. */
  const std::vector<unsigned short int>& node_types;
	
  /* Subgraphs (event sets) where root_event_id is the first
     one. */
  TSubgraphList subgraphs;
	
  /* Find all motifs corresponding to valid subgraphs up to size
   * `max_subgraph_size` where the root event is the first
   * event. The method find_submotifs() can only be accessed
   * through the public iterator interface; create_submotifs()
   * takes care of the recursive search.
   */
  void find_subgraphs();
  void create_subgraphs(const EventSet& eventMap,
			const EventSet& excludedEvents,
			const EventMMap& validNeighbors,
			unsigned int dt_max);
  void add_subgraph(const EventSet& eventMap, unsigned int dt_max);
	
 public:
  /* Simple constructor, only initializes parameters. */
  TSubgraphFinder(event_id root_event_id,
		 unsigned int time_window,
		 unsigned int max_submotif_size,
		 Events const& events,
		 std::vector<unsigned short int> const& node_types);
	
  /* Methods for iterating through subgraphs. This actually
     first finds all subgraphs and saves them into a list. The
     pointers to the subgraphs are owned by the subgraph finder,
     and they will be deleted when the finder object goes out of
     scope. */
  typedef TSubgraphList::iterator iterator;
  inline iterator begin() { find_subgraphs(); return subgraphs.begin(); };
  inline iterator end() { return subgraphs.end(); };

  /* Delete the pointers to subgraphs. */
  ~TSubgraphFinder();
	
};



#endif
