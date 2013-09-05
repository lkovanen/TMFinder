/* Representation for mobile phone event data.
 *
 * Lauri Kovanen, BECS (June 2010)
 */

#ifndef EVENTS_H
#define EVENTS_H

#include <assert.h>
#include <string>
#include <iostream>
#include <vector>
#include <list>
#include <set>
#include <sstream>
#include <algorithm>
#include <map>
#include <limits>
#include <cstdlib>
#include <stdint.h>
#include <math.h>
#include "fixed_tree.h"
#include "std_printers.h"

typedef uint32_t event_id;
typedef uint32_t node_id;
typedef FixedTree<event_id> event_tree;
typedef event_tree::iterator node_iterator;
typedef std::multimap<unsigned int, event_id> EventMMap;

class Events;

class Event
{
 public:
  /* Denote "no event" by the max unsigned int. If someone actually
     has that many events this will be a very nasty bug, but this is
     very unlikely and if this happens you'd be screwed anyway.*/
  static const event_id null_event;

 private:
  event_id _id;
  node_id _fr;
  node_id _to;
  unsigned int _start_time;
  unsigned int _end_time;
  short int _type;
  event_id _component_id;
  
 public:

  void Init(event_id id, node_id fr, node_id to,
	    unsigned int start_time, unsigned int duration,
	    short int type);

  void Reset(node_id fr, node_id to);

  inline event_id id() const {return _id;};
  inline node_id from() const {return _fr;};
  inline node_id to() const {return _to;};
  inline node_id other_node(node_id node) const {return (node==_fr?_to:_fr);};
  inline unsigned int start_time() const {return _start_time; };
  inline unsigned int duration() const {return _end_time-_start_time;};
  inline unsigned int end_time() const {return _end_time;};
  inline short int type() const {return _type;};
  inline void set_type(short int new_type) { _type = new_type;};
  inline event_id component() const {return _component_id;};
  inline bool has_component() const {return _component_id != null_event;};

  inline void set_component(event_id cid) {_component_id = cid;};

  void print() const
  {
    std::cerr << "Event " << _id 
	      << ": t = " << _start_time << "-" << _end_time << ", "
          << _fr << " -> " << _to << " [" << _type << "]" << std::endl;
  }
};

std::ostream& operator<<(std::ostream& output, const Event& e);

class Events
{
 private:
  /* Contains all events.
   */
  std::vector<Event> events;

  /* A set of events where a node is involved. This allows iterating
     over events of a single node.
   */
  std::vector<event_tree> node_events;

  /* The first and last time in data. */
  unsigned int t_first, t_last, t_last_start; 

  /* Switch the time of events i and j. This method does not really
     change the time, but instead all other data except time. This way
     the events will remain sorted.
   */
  void switch_times(event_id i, event_id j);

  /* Return true if the times of the input events overlap. Assumes
     that i_first < i_second and that neither is equal to
     Event::null_event.
   */
  bool check_overlap(event_id i_first, event_id i_second);

 public:

  inline unsigned int size() const {return events.size();};
  inline unsigned int get_nof_events() const {return size();};
  inline unsigned int get_nof_nodes() const {return node_events.size();};
  inline unsigned int first_time() const {return t_first;};
  inline unsigned int last_time() const {return t_last;};
  inline unsigned int last_start_time() const {return t_last_start;};

  /* Time difference between two events. */
  inline unsigned int dt(event_id i_1, event_id i_2) const 
  { 
    return events[i_2].start_time()-events[i_1].end_time();
  }

  /* Randomly shuffle event times. This method will also reset the
     component id of all events.
  */
  void shuffle();

  /* Randomly shuffle event types.
   */
  void shuffle_event_types();

  /* Randomly shuffle event types assuming all events on a given edge
     have the same type. Returns false if this assumption fails. Note
     that this method does not retain the number of events of each
     type, but the number of edges of each type.
   */
  bool shuffle_edge_types();

  /* At each time step two random events are selected for
     shuffling. The total number of (valid) selections is
     N_events*N_shuffle. The value of N_shuffle affects the goodness
     of the shuffling. The following list tells what fraction of data
     points will (on average) remain in the original position:

        N_shuffle    f(not shuffled) = e^(-2*N_shuffle)
	    2               1.83 %
	    3               0.248 %
	    4               0.0336 %
	    5               0.00454 %
	    6               0.000614 %
	    7               0.0000832 %
	    8               0.0000113 %
  */
  void shuffle_constrained(unsigned int N_shuffle);

  /* Shuffling with artificial correlation. At each of the
     N_events*N_shuffle iterations selects first one event i for
     switching, and then selects N_corr other events and picks the one
     that gives closest distance to other event of the nodes in
     events i. "N_corr = 1" corresponds to unbiased shuffling.
  */
  void shuffle_constrained_corr(unsigned int N_shuffle, unsigned int N_corr);

  /* Check that the events are properly constructed.
   */
  void check_events() const;

  // The constructor reads in the events from a file.
  Events(std::istream& event_file);
  ~Events() {};

  void print() const;

  inline Event& operator[](event_id event_id) { return events[event_id]; };
  inline const Event& operator[](event_id event_id) const { return events[event_id]; };

  /* Interface for iterating over the events of a single node.
   */
  inline node_iterator find_node_event(node_id node, event_id i) const { return node_events[node].find(i); };
  inline node_iterator begin(node_id node) const {return node_events[node].begin(); };
  inline node_iterator end(node_id node) const {return node_events[node].end(); };
  inline node_iterator rbegin(node_id node) const {return node_events[node].rbegin(); };
  inline node_iterator rend(node_id node) const {return node_events[node].rend(); };

  /* Interface for iterating through all events. (This is just a
     shortcut to iterating through the events list.) */
  typedef std::vector<Event>::iterator iterator;
  typedef std::vector<Event>::const_iterator const_iterator;
  inline iterator begin() {return events.begin();};
  inline iterator end() {return events.end();};
  inline const_iterator begin() const {return events.begin();};
  inline const_iterator end() const {return events.end();};

  /* Get the immediate next and previous events.
  */
  void next_immediate_events(event_id e_id, EventMMap& next_events) const;
  void prev_immediate_events(event_id e_id, EventMMap& prev_events) const;

  /* Identify maximal subgraphs with given time window. The ID of
     maximal subgraphs is set as the component id of each event.
  */
  void find_maximal_subgraphs(unsigned int tw);
};


#endif
