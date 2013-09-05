/* Representation for mobile phone event data.
 *
 * Lauri Kovanen, BECS (June 2010)
 */
#include <iostream>
#include <math.h>
#include "events.h"

const event_id Event::null_event = std::numeric_limits<event_id>::max();

void Event::Init(event_id id, node_id fr, node_id to,
		 unsigned int start_time, unsigned int duration,
		 short int type)
{
  _id = id;
  _start_time = start_time;
  _end_time = start_time+duration;
  _type = type;
  Reset(fr, to);
}

void Event::Reset(node_id fr, node_id to)
{
  _fr = fr;
  _to = to;
  _component_id = Event::null_event;
}

std::ostream& operator<<(std::ostream& output, const Event& e) 
{
  output << "id " << e.id() << " (t " << e.start_time() << "): "
	 << e.from() << "->" << e.to() << " [" << e.type() << "]"; 
  return output;
}

Events::Events(std::istream& event_file):events(1000),
					 node_events(1000),
					 t_first(0),
					 t_last(0),
					 t_last_start(0)
{
  // Set initial values for the number of nodes and the number and
  // events. These will be scaled up or down when necessary.
  unsigned int N_events = events.size();
  unsigned int N_nodes = node_events.size();

  //std::vector<std::list<event_id> > tmp_node_events(N_nodes_init);

  assert(sizeof(event_id) >= 4);
  event_id max_node_id = 0;
  event_id id = 0;

  // Read in the events.
  if (event_file.good())
    {
      std::string line;
      while (!event_file.eof())
	{
	  getline(event_file, line);
	  if (!line.empty())
	    {
	      // Read in the basic information about this event.
	      std::istringstream is(line);
	      //std::cerr << line << std::endl;
	      if (id >= N_events)
		{
		  //N_events = (node_id)ceil(1.2*(1+id));
		  N_events = 1+id;
		  events.resize(N_events);
		}
	      Event & e = events[id];
	      unsigned int start_time, duration, fr, to;
	      short int event_type = 1;
	      is >> start_time;
	      is >> duration;
	      is >> fr;
	      is >> to;
	      assert(fr != to);
	      if (!is.eof()) is >> event_type;
	      e.Init(id, fr, to, start_time, duration, event_type);

	      max_node_id = std::max(max_node_id, std::max(fr, to));
	      if ((fr >= N_nodes) || (to >= N_nodes))
		{
		  //std::cerr << "Increasing node_events size from " 
		  //	    << N_nodes << " to " << max_node_id+1 << "\n";
		  N_nodes = max_node_id+1;
		  //tmp_node_events.resize(N_nodes);
		  node_events.resize(N_nodes);
		}

	      if (e.end_time() > t_last) t_last = e.end_time();
	      //tmp_node_events[fr].push_back(id);
	      //tmp_node_events[to].push_back(id);
	      //std::cerr << id << ": for node " << fr << " (" << node_events[fr].size() << ")\n";
	      node_events[fr].add(id);

	      //std::cerr << id << ": for node " << to << " (" << node_events[to].size() << ")\n";
	      node_events[to].add(id);

	      //e.print();
	      ++id;
	    }
	}
    }

  /* Check if the number of events and nodes matches and resize the
   * event vector if necessary. Note that at this point the value of
   * 'id' is one larger than the largest value added.
   */
  //std::cerr << id << " " << N_events << std::endl;
  //std::cerr << max_node_id << " " << N_nodes << std::endl;
  if (id < N_events)
    {
      N_events = id;
      //std::cerr << "Resizing ...\n";
      events.resize(N_events);
      //std::cerr << "done.\n";
    }
  if (max_node_id + 1 < N_nodes)
    {
      N_nodes = max_node_id + 1;
      //tmp_node_events.resize(N_nodes);
      node_events.resize(N_nodes);
    }

  //node_events.resize(N_nodes);
  //std::cerr << "  sizeof(Event) = " << sizeof(Event) << "\n";
  //std::cerr << "  sizeof(event_tree) = " << sizeof(event_tree) << "\n";
  //std::cerr << "  sizeof(FixedNode<event_id>) = " << sizeof(FixedNode<event_id>) << "\n";

  // Get the starting times of the first and last events.
  t_first = events[0].start_time();
  t_last_start = events[N_events-1].start_time();
  
  // Initialize node_events.
  //std::cerr << "   Initialize node_events.\n";
  for (unsigned int i = 0; i < node_events.size(); ++i)
    {
        if (!node_events[i].empty()) node_events[i].Init();
    }

  std::cout << "   Events read, found "
	    << get_nof_nodes() << " nodes and " 
	    << get_nof_events() << " events.\n";

}

void Events::switch_times(event_id i, event_id j)
{
  node_id i_from = events[i].from();
  node_id i_to = events[i].to();
  node_id j_from = events[j].from();
  node_id j_to = events[j].to();
  
  // Switch all other data of events i and j except those
  // related to time and id.
  events[i].Reset(j_from, j_to);
  events[j].Reset(i_from, i_to);
  
  //std::cerr << "Remove events ...\n";
  node_events[i_from].replace(i,j);
  node_events[i_to].replace(i,j);
  node_events[j_from].replace(j,i);
  node_events[j_to].replace(j,i);
}

void Events::shuffle()
{
  unsigned int N_events = get_nof_events();

  //srand(time(NULL));

  for (event_id i = 0; i < N_events; ++i)
    {
      // Get random number from U(i,N_events-1).
      double diff = ((N_events-i)*((double)rand()))/(RAND_MAX+1.0);
      event_id j = i + (event_id)diff;
      std::cerr << "Shuffling events " << i << " and " << j << std::endl;
      /*
      std::cerr << "  Events of node fr("<<i<<")=" << events[i].from() << ": ";
      node_events[events[i].from()].print();
      std::cerr << "  Events of node fr("<<i<<")=" << events[i].to() << ": ";
      node_events[events[i].from()].print();
      std::cerr << "  Events of node fr("<<j<<")=" << events[i].from() << ": ";
      node_events[events[j].from()].print();
      std::cerr << "  Events of node fr("<<j<<")=" << events[i].to() << ": ";
      node_events[events[j].from()].print();
      */
      if (i != j) switch_times(i, j);
    }
};

void Events::shuffle_event_types()
{
  unsigned int N_events = get_nof_events();
  for (event_id i = 0; i < N_events; ++i)
    {
      // Get random number from U(i,N_events-1).
      double diff = ((N_events-i)*((double)rand()))/(RAND_MAX+1.0);
      event_id j = i + (event_id)diff;
      //std::cerr << "Shuffling types of events " << i << " and " << j << std::endl;
      if (i != j)
	{
	  short int i_type = events[i].type();
	  events[i].set_type(events[j].type());
	  events[j].set_type(i_type);
	}
    }
};

bool Events::shuffle_edge_types()
{
  // Get the type of each edge. Returns false if for some edge there
  // are two types of events. We collect the types in a separate
  // vector so they can be easily shuffled.
  unsigned int i_edge = 0;
  std::map<std::pair<node_id, node_id>, unsigned int> edges;
  std::vector<short int> edge_types;
  for (event_id i = 0; i < get_nof_events(); ++i)
    {
      std::pair<node_id, node_id> curr_edge(events[i].from(), events[i].to());
      std::map<std::pair<node_id, node_id>, unsigned int>::const_iterator ed_it = edges.find(curr_edge);
      if (ed_it == edges.end()) 
	{
	  edges[curr_edge] = i_edge;
	  edge_types.push_back(events[i].type());
	  i_edge++;
	}
      else if (edge_types[ed_it->second] != events[i].type()) return false;
    }

  // Shuffle event types.
  std::random_shuffle(edge_types.begin(), edge_types.end());

  // Re-assign randomized event types.
  for (event_id i = 0; i < get_nof_events(); ++i)
    {
      std::pair<node_id, node_id> curr_edge(events[i].from(), events[i].to());
      events[i].set_type(edge_types[edges[curr_edge]]);
    }
  return true;
};

bool Events::check_overlap(event_id i_first, event_id i_second)
{
  if (events[i_first].end_time() >= events[i_second].start_time())
    return true;
  return false;
}

void Events::shuffle_constrained(unsigned int N_shuffle)
{
  unsigned int N_nodes = get_nof_nodes();
  unsigned int N_events = get_nof_events();

  std::cerr << N_nodes << " nodes, " << N_events << " events." << std::endl << std::flush;

  // Randomize first the time stamps, then build up again the pointers
  // to next and previous events.
  unsigned int n_shuffles = 0;
  unsigned int shuffle_tries = 0;
  while (n_shuffles < N_events*N_shuffle)
    {
      ++shuffle_tries;

      if (shuffle_tries % 10000000 == 0) 
	{
	  float p_done = ((float)n_shuffles)/(N_events*N_shuffle);
	  std::cerr << "    Shuffled " << n_shuffles << " out of " 
		    << shuffle_tries << " tries (" 
		    << (int)(100*p_done) << "% done)"
		    << std::endl << std::flush;
	}

      event_id i,j;
      i = (event_id)(N_events*(rand()/(RAND_MAX+1.0)));
      __builtin_prefetch(node_events[events[i].from()].nodes, 0, 3);
      __builtin_prefetch(node_events[events[i].to()].nodes, 0, 3);

      do {
	j = (event_id)(N_events*(rand()/(RAND_MAX+1.0)));
      } while (events[i].type() != events[j].type() || i == j);
      __builtin_prefetch(node_events[events[j].from()].nodes, 0, 3);
      __builtin_prefetch(node_events[events[j].to()].nodes, 0, 3);

      //std::cerr << "Trying to shuffle " << i << " and " << j << std::endl;

      // The only two events that might overlap after the switch are
      // those with id exactly smaller and larger than j: the first
      // one might have started before j but end during j; the second
      // one might start during j. There can be other events that
      // overlap, but if these two do not overlap we know that no
      // other event will either.

      bool overlap_found = false;
      event_id i0 = i;
      event_id i1 = j;
      for (int e_ = 0; e_ < 2; ++e_)
	{
	  Event const& e = events[i0];
	  node_id tmp_node = e.from();
	  for (int u_ = 0; u_ < 2; ++u_)
	    {
	      event_id i2;
	      
	      i2 = node_events[tmp_node].find_prev(i1, Event::null_event);
	      if (i2 == i0) node_events[tmp_node].find_prev(i2, Event::null_event);
	      if (i2 != Event::null_event && check_overlap(i2, i1))
		{
		  overlap_found = true; break;
		}

	      i2 = node_events[tmp_node].find_next(i1, Event::null_event);
	      if (i2 == i0) node_events[tmp_node].find_next(i2, Event::null_event);
	      if (i2 != Event::null_event && check_overlap(i1, i2))
		{
		  overlap_found = true; break;
		}

	      tmp_node = e.to(); // Repeat for the other node.
	    }
	  if (overlap_found) break;
	  
	  i0 = j; i1 = i; // Repeat for the other event.
	}
      if (overlap_found) continue;
	  
      // If we got this far we know the switch is valid, i.e. after
      // making the switch there will be no overlapping events.
      //std::cerr << "Shuffling events " << i << " and " << j << std::endl;
      switch_times(i,j);

      // Increase the successful shuffles count.
      ++n_shuffles;

      //std::cerr << "Check events:" << std::endl;
      //check_events(); // FOR DEBUGGING ONLY!
    }

  // Restore order of the underlying data structures of node events after shuffling.
  std::cerr << "Restore order ...\n";
  std::vector<event_tree>::iterator uit;
  for (uit = node_events.begin(); uit != node_events.end(); ++uit) uit->restore_order();

  std::cerr << "Accepted " << n_shuffles << "/" << shuffle_tries 
	    << " switches during shuffling.\n";
};


void Events::shuffle_constrained_corr(unsigned int N_shuffle, unsigned int N_corr)
{
  unsigned int N_nodes = get_nof_nodes();
  unsigned int N_events = get_nof_events();

  std::cerr << N_nodes << " nodes, " << N_events << " events." << std::endl << std::flush;
  std::cerr << "Shuffling a total of " << N_events*N_shuffle << " times." << std::endl << std::flush;

  unsigned int n_shuffles = 0;
  unsigned int shuffle_tries = 0;
  while (n_shuffles < N_events*N_shuffle)
    {
      ++shuffle_tries;

      if (shuffle_tries % 10000000 == 0) 
	{
	  float p_done = ((float)n_shuffles)/(N_events*N_shuffle);
	  std::cerr << "    Shuffled " << n_shuffles << " out of " 
		    << shuffle_tries << " tries (" 
		    << (int)(100*p_done) << "% done)"
		    << std::endl << std::flush;
	}

      bool overlap_found = false;

      // Get the first event.
      event_id i;
      i = (event_id)(N_events*(rand()/(RAND_MAX+1.0)));
      __builtin_prefetch(node_events[events[i].from()].nodes, 0, 3);
      __builtin_prefetch(node_events[events[i].to()].nodes, 0, 3);
      Event const& e_i = events[i];

      //std::cerr << "Trying to shuffle " << i << " with ..." << std::endl;

      // Get N_corr other events and then select the one closest to i.
      event_id j = Event::null_event;
      unsigned int diff_best = t_last;
      for (event_id i_rnd = 0; i_rnd < N_corr; ++i_rnd)
	{
	  // Get another random event ...
	  event_id j_try;
	  do {
	    j_try = (event_id)(N_events*(rand()/(RAND_MAX+1.0)));
	  } while (events[i].type() != events[j_try].type() || i == j_try);

	  // ... and calculate how close it would be to other events
	  // of nodes in e_i after switching the times.
	  //unsigned int diff_min = t_last;
	  double diff_gmean = 0;
	  unsigned short N_adjacent_events = 0;

	  node_id tmp_node = e_i.from();
	  for (int u_ = 0; u_ < 2; ++u_)
	  {
	    event_id i_prev = node_events[tmp_node].find_prev(j_try, Event::null_event);
	    if (i_prev == i) i_prev = node_events[tmp_node].find_prev(i_prev, Event::null_event);
	    if (i_prev != Event::null_event)
	      {
		int diff = events[j_try].start_time() - events[i_prev].end_time();
		if (diff <= 0) {
		  overlap_found = true; break;
		}
		diff_gmean += log(diff); N_adjacent_events++;
	      }

	    event_id i_next = node_events[tmp_node].find_next(j_try, Event::null_event);
	    if (i_next == i) i_next = node_events[tmp_node].find_next(i_next, Event::null_event);
	    if (i_next != Event::null_event)
	      {
		int diff = events[i_next].start_time() - events[j_try].end_time();
		if (diff <= 0) {
		  overlap_found = true; break;
		}
		diff_gmean += log(diff); N_adjacent_events++;
	      }
	    
	    tmp_node = e_i.to(); // Repeat for the other node.
	  }
	  if (overlap_found) continue;
	  
	  //std::cerr << "   " << j_try << " (" << diff_min << ")" << std::endl;	  

	  // To get the geometric mean we should take the exponent
	  // here, but because we are only interested in the order of
	  // the values we can skip that. If 'i' had no neighboring
	  // events we simply pick the first random event that does
	  // not overlap.
	  if (N_adjacent_events == 0)
	    {
	      j = j_try;
	      break;
	    }
	  diff_gmean = diff_gmean/(double)N_adjacent_events;
	  if (diff_gmean < diff_best)
	    {
	      j = j_try;
	      diff_best = diff_gmean;
	    }
	}

      // Make sure at least one non-overlapping other event was found.
      if (j == Event::null_event) continue;
      __builtin_prefetch(node_events[events[j].from()].nodes, 0, 3);
      __builtin_prefetch(node_events[events[j].to()].nodes, 0, 3);
      Event const& e_j = events[j];

      // Check the overlapping in the other direction (at the nodes of
      // j).
      node_id tmp_node = e_j.from();
      for (int u_ = 0; u_ < 2; ++u_)
	{
	  event_id j_prev = node_events[tmp_node].find_prev(i, Event::null_event);
	  if (j_prev == j) j_prev = node_events[tmp_node].find_prev(j_prev, Event::null_event);	  
	  if (j_prev != Event::null_event && (e_i.start_time() < events[j_prev].end_time()))
	    {
	      overlap_found = true; break;
	    }

	  event_id j_next = node_events[tmp_node].find_next(i, Event::null_event);
	  if (j_next == j) j_next = node_events[tmp_node].find_next(j_next, Event::null_event);	  
	  if (j_next != Event::null_event && (events[j_next].start_time() < e_i.end_time()))
	    {
	      overlap_found = true; break;
	    }
	  tmp_node = e_j.to(); // Repeat for the other node.
	}
      if (overlap_found) continue;

      // If we got this far we know the switch is valid, i.e. after
      // making the switch there will be no overlapping events.
      //std::cerr << "Shuffling events " << i << " and " << j << std::endl;
      switch_times(i,j);

      // Increase the successful shuffles count.
      ++n_shuffles;
    }

  // Restore order of the underlying data structures of node events after shuffling.
  std::cerr << "Restore order ...\n";
  std::vector<event_tree>::iterator uit;
  for (uit = node_events.begin(); uit != node_events.end(); ++uit) uit->restore_order();

  std::cerr << "Accepted " << n_shuffles << "/" << shuffle_tries 
	    << " switches during shuffling.\n";
};

void Events::print() const
{
  //std::cerr << "First events: " << first_events << std::endl;
  //std::cerr << "Last events: " << last_events << std::endl;
  Events::const_iterator it;
  for (it = events.begin(); it != events.end(); ++it)
      it->print();
  node_id node;
  for (node = 0; node < get_nof_nodes(); ++node)
    {
      //node_events[node].debug_print();
      std::cerr << "Node " << node << ": ";
      node_iterator uit = begin(node);
      for (; uit != end(node); ++uit) std::cerr << *uit << " ";
      std::cerr << std::endl;
    }
};



void Events::next_immediate_events(event_id e_id, EventMMap& next_events) const
{
  const Event& e = events[e_id];

  node_id node = e.from();
  node_iterator it = find_node_event(node, e_id); it++;
  if (it == end(node))
    {
      // Simple case: This node does not have a next event, so the
      // only possible next event is that of the other node, if it has
      // one.
      node = e.to();
      it = find_node_event(node, e_id); it++;
      if (it != end(node)) next_events.insert(std::make_pair(dt(e_id,*it),*it));
      return;
    }
  event_id e_fr = *it;

  // The first node has a next event. Now we need to find the next
  // event of the other node to figure out what to do next.
  node = e.to();
  it = find_node_event(node, e_id); it++;
  if (it == end(node))
    {
      // Simple case: The second node doesn't have a next event within
      // tw, so just return the next event of the first node.
      next_events.insert(std::make_pair(dt(e_id,e_fr),e_fr));
      return;
    }
  event_id e_to = *it;

  // Both nodes have a next event within tw. If its the same one, then
  // it takes place between the same two nodes and we can safely
  // return it.
  if (e_fr == e_to)
    {
      next_events.insert(std::make_pair(dt(e_id,e_fr),e_fr));
      return;
    }

  // The two nodes have different next events. Both cannot be on this
  // same edge (because then they would be the same event), but it is
  // possible that one of them is on the same edge (or neither is). If
  // one of them is on the same edge, it takes place after the other
  // one (if it were the first one, it would have been returned for
  // both nodes) so it is only the other one we return. Only when
  // neither is on the same edge we return both events. In summary, at
  // this point we can safely return events that are not on the same
  // edge.
  node_id third_node;
  third_node = events[e_fr].other_node(e.from());
  if (third_node != e.to()) next_events.insert(std::make_pair(dt(e_id,e_fr),e_fr));
  third_node = events[e_to].other_node(e.to());
  if (third_node != e.from()) next_events.insert(std::make_pair(dt(e_id,e_to),e_to));

  return;
}

void Events::prev_immediate_events(event_id e_id, EventMMap& prev_events) const
{
  const Event& e = events[e_id];
  node_id node = e.from();

  node_iterator it = find_node_event(node, e_id); it--;
  if (it == rend(node))
    {
      // Simple case: This node does not have a previous event, so the
      // only previous event is that of the other node, if it has one.
      node = e.to();
      it = find_node_event(node, e_id); it--;
      if (it != rend(node)) prev_events.insert(std::make_pair(dt(*it,e_id),*it));
      return;
    }
  event_id e_fr = *it;

  // The first node has a previous event. Now we need to find the
  // previous event of the other node to figure out what to do next.
  node = e.to();
  it = find_node_event(node, e_id); it--;
  if (it == rend(node))
    {
      // Simple case: The second node doesn't have a previous event
      // within tw, so just return the previous event of the first
      // node.
      prev_events.insert(std::make_pair(dt(e_fr,e_id),e_fr));
      return;
    }
  event_id e_to = *it;

  // Both nodes have a previous event within tw. If its the same one,
  // then it takes place between the same two nodes and we can safely
  // return it.
  if (e_fr == e_to)
    {
      prev_events.insert(std::make_pair(dt(e_fr,e_id),e_fr));
      return;
    }

  // The two nodes have different previous events. Both cannot be on
  // this same edge (because then they would be the same event), but
  // it is possible that one of them is on the same edge (or neither
  // is). If one of them is on the same edge, it takes place before
  // the other one (if it were the later one, it would have been the
  // previous event of both nodes) so it is only the other one we
  // return. Only when neither is on the same edge do we return both
  // events. In summary, at this point we can safely return events
  // that are not on the same edge.
  node_id third_node;
  third_node = events[e_fr].other_node(e.from());
  if (third_node != e.to()) prev_events.insert(std::make_pair(dt(e_fr,e_id),e_fr));
  third_node = events[e_to].other_node(e.to());
  if (third_node != e.from()) prev_events.insert(std::make_pair(dt(e_to,e_id),e_to));

  return;
}

void Events::check_events() const
{
  // Go through all events and make sure that the event is listed for
  // both nodes in node_events.

  std::vector<Event>::const_iterator it;
  for (it = events.begin(); it != events.end(); ++it)
    {
      node_id node = it->from();
      for (int i = 0; i < 2; ++i)
	{
	  if (node_events[node].empty())
	    {
	      std::cerr << "Error: Node " << node << " has no event listed but " 
			<< " involved in event " << it->id() << ".\n";
	    }
	  else if (find_node_event(node, it->id()) == end(node))
	    {
	      std::cerr << "Error: Event " << it->id() << " not list in events " 
			<< " of node " << node << ".\n";
	    }
	  node = it->to();
	}
    }

  // Go through all node events and make sure that the node is involved in the event.
  node_id node;
  for (node = 0; node < get_nof_nodes(); ++node)
    {
      node_iterator uit = begin(node);
      for (; uit != end(node); ++uit)
	{
	  Event const& e = events[*uit];
	  if (node != e.from() && node != e.to())
	    {
	      std::cerr << "Error: Node " << node << " not involved in " 
			<< " event " << e.id() << ".\n";
	    }
	}
    }
}

void Events::find_maximal_subgraphs(unsigned int tw)
{
  // Go through all events and recursively find the maximal temporal
  // subgraph. The maximal subgraph id will be the event id of the
  // earliest id in the subgraph (this happens implicitely because we
  // go through the events in temporal order).
  iterator e_it;
  for (e_it = begin(); e_it != end(); ++e_it)
    {
      // See if this event has already been assigned to a temporal
      // component.
      if (e_it->has_component()) continue;

      unsigned int component_id = e_it->id();
      std::set<event_id> to_process;
      to_process.insert(component_id);

      while (!to_process.empty())
	{
	  // Take an event from processing queue and set its component
	  // id.
	  std::set<event_id>::iterator it = to_process.begin();
	  Event& e = events[*it];
	  e.set_component(component_id);
	  to_process.erase(it);

	  // Add immediate neighbors of this event to queue, unless
	  // they already have a community id (which shows that they
	  // have already been processed).
	  EventMMap neighbors;
	  prev_immediate_events(e.id(), neighbors);
	  next_immediate_events(e.id(), neighbors);
	  for (EventMMap::const_iterator it = neighbors.begin(); it != neighbors.end(); it++)
	    {
	      if (it->first > tw) break;
	      if (!events[it->second].has_component()) to_process.insert(it->second);
	    }
	}
    }
}
