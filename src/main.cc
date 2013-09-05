/* TMFinder - Identify motifs in temporal networks

Copyright (C) 2013  Lauri Kovanen

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see {http://www.gnu.org/licenses/}.
*/

#include <utility>
#include <math.h>
#include <time.h>
#include <iterator>
#include "events.h"
#include "tsubgraph.h"
#include "subnets.h"
#include "binner.h"
#include "motif_counter.h"
#include "progress_counter.h"
#include "lcelib/Nets.H"
#include "edges.h"
#include "bin_limits.h"

// LocationMap[motif_hash][edge_id_list] = count

// weightsMap[motif_hash].add(weightVector, value)
// weightsMap[motif_hash].get_mean(weightVector, result)
typedef Binner<unsigned int> wBinner;
typedef std::vector<unsigned int> WeightVector;
typedef std::map<unsigned int, wBinner> WeightsMap;

typedef DirNet<unsigned int> NetType;

typedef std::vector<short int> TypeSeq;
typedef std::vector<TypeSeq> TypeSeqs;
typedef std::map<unsigned int, TypeSeqs> TypeSeqsMap;

// Get current date/time, format is YYYY-MM-DD.HH:mm:ss
const std::string currentDateTime() {
  time_t now = time(0);
  struct tm tstruct;
  char buf[80];
  tstruct = *localtime(&now);
  // Visit http://www.cplusplus.com/reference/clibrary/ctime/strftime/
  // for more information about date/time format
  strftime(buf, sizeof(buf), "%Y-%m-%d %X", &tstruct);
  return buf;
}

bool update_location_count(const TSubgraph& sg, 
			   EdgeVectorMap& locationMap)
{
  if (!sg.is_valid()) return false;

  // Note that if the motif hash or edge vector are not found in the
  // corresponding maps, they are automatically created by
  // map.operator[].
  EdgeVector edges(sg.begin(), sg.end());
  locationMap[edges]++;
  return true;
}

unsigned int get_location_count(const EdgeVectorMap& locationMap,
				const EdgeVector& edges)
{
  EdgeVectorMap::const_iterator ev_it = locationMap.find(edges);
  if (ev_it == locationMap.end()) return 0;
  return ev_it->second;
}

bool get_edge_weights(const EdgeVector& edges,
		      const std::map<short int, NetType*>& nets,
		      std::vector<unsigned int>& weights)
{
  weights.clear();
  for (EdgeVector::const_iterator ev_it = edges.begin();
       ev_it != edges.end(); ++ev_it) 
    {
      const NetType& net = *(nets.find(ev_it->type)->second);
      unsigned int w = net[ev_it->from][ev_it->to];
      if (w == 0) return false;
      weights.push_back(w);
    }
  return true;
}

unsigned int read_node_types(std::vector<unsigned short int>& node_types, std::string node_file_name)
{
  std::ifstream node_file(node_file_name.c_str(), std::ifstream::in);
  unsigned int node_count = 0;
  if (node_file.is_open())
    {
      std::cerr << "Reading node types ...\n";
      std::string line;
      while (node_file.good())
        {
	  getline(node_file, line);
	  if (!line.empty())
            {
	      std::istringstream is(line);
	      unsigned int node_id;
	      unsigned short int node_type;
	      is >> node_id;
	      is >> node_type;
	      if (node_id >= node_types.size()) node_types.resize(node_id+1);
	      node_types[node_id] = node_type;
	      node_count++;
            }
        }
    }
  node_file.close();
  return node_count;
}

bool fill_event_type_seq(TypeSeqsMap& ets, unsigned int n_max, const std::set<short int>& eventTypes, bool multiple_event_types)
{
  if (multiple_event_types)
    {
      // Initialize with an empty vector.
      TypeSeq c0;
      ets[0].push_back(c0);

      for (unsigned int n = 1; n <= n_max; ++n) {
	for (std::set<short int>::const_iterator et_it = eventTypes.begin();
	     et_it != eventTypes.end(); ++et_it)
	  {
	    for (TypeSeqs::const_iterator cp_it = ets[n-1].begin();
		 cp_it != ets[n-1].end(); ++cp_it)
	      {
		// Combine the previous combination and the event type.
		TypeSeq cn(cp_it->begin(), cp_it->end());
		cn.push_back(*et_it);
		// Add to map.
		ets[n].push_back(cn);
	      }
	  }
      }
    }
  else //Only single event type in a motif.
    {
      for (std::set<short int>::const_iterator et_it = eventTypes.begin();
	   et_it != eventTypes.end(); ++et_it)
        {
	  for (unsigned int n = 1; n <= n_max; ++n)
            {
	      TypeSeq c(n,*et_it);
	      ets[n].push_back(c);
            }
        }
    }
  return true;
}


// Class for reading the input parameters.
class Parameters
{
private:
  const static int N_required_param = 2;
  bool verbose;

  void print_licence()
  {
      std::cout << "TMFinder - Identify motifs in temporal networks\n\n"
          << "Copyright (C) 2013  Lauri Kovanen\n"
          << "\n"
          << "This program is free software: you can redistribute it and/or modify\n"
          << "it under the terms of the GNU General Public License as published by\n"
          << "the Free Software Foundation, either version 3 of the License, or\n"
          << "(at your option) any later version.\n"
          << "\n"
          << "This program is distributed in the hope that it will be useful,\n"
          << "but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
          << "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n"
          << "GNU General Public License for more details.\n"
          << "\n"
          << "You should have received a copy of the GNU General Public License\n"
          << "along with this program.  If not, see http://www.gnu.org/licenses/.\n"
          << "\n";
  };

  void print_help()
  {
    std::cout << "TMFinder  Copyright (C) 2013  Lauri Kovanen\n\n"
          << "This program comes with ABSOLUTELY NO WARRANTY. This is free software,\n"
          << "and you are welcome to redistribute it under certain conditions.\n"
          << "Type `./tmf --licence' for details.\n"
          << "\n"
          << "Usage:\n"
	      << "   ./tmf TW OUTPUTNAME < EVENTFILE\n\n"
	      << "The input event data has one event per line, with the first four columns corresponding to\n"
	      << "starting time, duration, and the id's of the two nodes involved. A fifth column may also be\n"
	      << "used to denote event type by an integer. When omitted, the type is assumed to be 1.\n\n"
	      << "There are two required parameters:\n\n"
	      << "  TW is the time window.\n\n"
	      << "  OUTPUTNAME is the beginning of the output file name.\n\n"
	      << "The optional parameters are:\n\n"  
	      << "-m INT | --max_size INT\n"
	      << "  The maximum number of events in valid subgraphs that are used to create motifs. If 0,\n"
	      << "  detect all subgraphs. This can take a very long time if the time window is large.\n\n"
	      << "--maximal\n"
	      << "  If given, detect only maximal subgraphs with at most '--max_size' events. If '--max_size'\n"
	      << "  is 0, detects all maximal subgraphs.\n\n"
	      << "-r INT | --references INT\n"
	      << "  The number of independent references to create. The references are created by generating\n"
	      << "  random motif counts at each location.\n\n"
	      << "-nf STR | --node_file STR\n"
	      << "  The file that contains the node types. It has one line per node, giving the node id\n"
	      << "  and type (both integers) separated by whitespace. The nodes that are not listed are\n"
	      << "  assumed to have type 0, which is the type of every node in case the file does not\n"
	      << "  exist. Note that node types must use different integers than the event types given\n"
	      << "  in the input data. Also, try to avoid using value 0 in the file as these nodes could\n"
	      << "  get mixed up with missing values.\n\n"
	      << "-wo FLOAT | --weight_omit FLOAT\n"
	      << "  The fraction of largest weights to exclude from the analysis. The purpose of excluding\n"
	      << "  largest weights is to reduce the possible bias caused by the very largest weight.\n"
	      << "  A value of 0.001 or below is typically good.\n\n"
	      << "--allow_multiple_event_types\n"
	      << "  By default only those motifs are detected where all events have the same type. If this\n"
	      << "  parameter is given, also includes motifs that have multiple event types. If this argument\n"
	      << "  is not used the run time increases by factor `N` compared to having only single event type\n"
	      << "  (here `N` is the number of event types). When this argument is given, the increase is\n"
	      << "  approximately `N^m` where `m` is the maximum number of events in a motif.\n\n"
	      << "-h INT | --hypothesis INT\n"
	      << "  The null hypothesis to test. There are three possible values:\n"
	      << "     0 : Motif count does not depend on node types. (default)\n"
	      << "     1 : Motif count does not depend on event types.\n"
	      << "     2 : Motif count does not depend on node or event types.\n"
	      << "  In the implementation these choices correspond to different motif hash that defines the\n"
	      << "  reference distribution. For '0' the hash is calculated after excluding node types,\n"
	      << "  for '1' after excluding event types, and for '2' after excluding both. Note that if the\n"
	      << "  data does not contain multiple event types, '0' and '2' give identical result and\n"
	      << "  '1' is non-sensical; if there is only one node type, '1' and '2' are identical and using \n"
	      << "  '0' makes no sense.\n\n"
	      << "-t INT | --time_gap INT\n"
	      << "  The time gap to leave in the beginning and end of the data: only those subgraphs\n"
	      << "  will be used that contain at least one event in the valid region. The purpose is to\n"
	      << "  avoid using data in the ends where maximal motifs are necessarily incomplete. The\n"
	      << "  default value is no gap.\n\n"
	      << "-st INT | --shuffle_type INT\n"
	      << "  If given, shuffled data is used instead of empirical one. The value is either\n"
	      << "      -1 : shuffle edge types\n"
	      << "       0 : shuffle node types\n"
	      << "       1 : shuffle event times (uniform)\n"
	      << "      >1 : shuffle event times (with bias corresponding to value)\n\n"
	      << "-s INT | --seed INT\n"
	      << "  The seed for the random number generator. If omitted the system time is used.\n"
	      << "\n"
	      << "The columns in the output files are\n"
	      << "  0 count    : count of temporal motif in input data\n"
	      << "  1 ref_avg  : mean count in the reference\n"
	      << "  2 ratio    : count/ref_avg (-1 if undefined)\n"
	      << "  3 ref_std  : standard deviation in the reference\n"
	      << "  4 N_ref    : number of references where this motif occurs at least once\n"
	      << "  5 N_lt_ref : number of references where count is smaller\n"
	      << "               than the ref count\n"      
	      << "  6 z-score  : (count-ref_avg)/ref_std (0 also if undefined)\n"
	      << "  7 N_loc_tot: number of locations in the aggregated network\n"
	      << "               where the motif could occur\n"
	      << "  8 N_loc    : number of locations in the aggregated network\n"
	      << "               where the motif occurs in the input data\n"
	      << "  9 N_loc_ref: mean number of locations in the aggregated network\n"
	      << "               where the motif occur in the reference\n"
	      << " 10 motif    : the motif itself\n"
	      << "\n";
  };

  bool get_parameter(int& i, int argc, char *argv[])
  {
    std::string name(argv[i]);
    if ((name.compare("-m") == 0) || (name.compare("--max_size") == 0))
      {
	i++; if (i > argc) return false;
	max_size = atoi(argv[i]);
      }
    else if (name.compare("--maximal") == 0)
      {
        maximal = true;
      }
    else if ((name.compare("-r") == 0) || (name.compare("--references") == 0))
      {
	i++; if (i > argc) return false;
	references = atoi(argv[i]);
      }
    else if ((name.compare("-nf") == 0) || (name.compare("--node_file") == 0))
      {
	i++; if (i > argc) return false;
	node_file_name = argv[i];
      }
    else if ((name.compare("-t") == 0) || (name.compare("--time_gap") == 0))
      {
	i++; if (i > argc) return false;
	if (atoi(argv[i]) < 0) return false;
	time_gap = atoi(argv[i]);
      } 
    else if ((name.compare("-wo") == 0) || (name.compare("--weight_omit") == 0))
      {
	i++; if (i > argc) return false;
	weight_omit = atof(argv[i]);
      } 
    else if (name.compare("--allow_multiple_event_types") == 0)
      {
	allow_multiple_event_types = true;
      }
    else if ((name.compare("-h") == 0) || (name.compare("--hypothesis") == 0))
      {
	i++; if (i > argc) return false;
	hypothesis = atoi(argv[i]);
	if ((hypothesis < 0) || (hypothesis > 2)) return false;
      }
    else if ((name.compare("-st") == 0) || (name.compare("--shuffle_type") == 0))
      {
	i++; if (i > argc) return false;
	int val = atoi(argv[i]);
	if (val == -1) edge_type_shuffling = true;
	else if (val == 0) node_type_shuffling = true;
	else if (val >= 1) {
	  time_shuffling = true;
	  bias_strength = val;
	}
	else return false;
      }
    else if ((name.compare("-s") == 0) || (name.compare("--seed") == 0))
      {
	i++; if (i > argc) return false;
	rng_seed = atoi(argv[i]);
      } 
    else
      {
	if (verbose) std::cout << "   Unidentified parameter '" << name << "'.\n";
	return false;
      }
    return true;
  };

  bool post_process()
  {
    if (verbose) 
      {
	std::cout << "Input parameters read (" << currentDateTime() << "):\n";
	std::cout << "   Time window : " << tw << std::endl;
	std::cout << "   Skipping " << time_gap << " time units.\n";
	if (weight_omit > 0.0) std::cout << "   Omitting highest " << weight_omit << " of edge weights." << std::endl;
      }

    // Construct file names. The value of max_size determines
    // whether only maximal motifs are detected or all motifs up to a
    // given size.
    output_file_name = output_file_trunk + ".dat";

    if (verbose) 
      {
	std::cout << "   Output file: " << output_file_name << std::endl;
	if (maximal)
	  {
	    if (max_size) std::cout << "   Finding maximal motifs with up to " << max_size << " events.\n";
	    else std::cout << "   Finding all maximal motifs.\n";
	  }
	else
	  {
	    if (max_size) std::cout << "   Finding motifs with up to " << max_size << " events.\n";
	    else std::cout << "   Finding all motifs with any number of events.\n";
	  }

	if (!node_file_name.empty())
	  {
	    std::cout << "   Reading node types from '" << node_file_name << "' (if the file exists)." << std::endl;
	  }
	else
	  {
	    std::cout << "   No node types!" << std::endl;
	    return false;
	  }

	if (references) std::cout << "   Creating "<<references<<" references (seed " << rng_seed << ")\n";
	else std::cout << "   No references.\n";

	if (allow_multiple_event_types) std::cout << "   Multiple event type motifs included (assuming there are multiple event types).\n";
	else std::cout << "   Including only motifs with single event type.\n";

	if (hypothesis == 0) std::cout << "   H_0 : Motif count does not depend on node types.\n";
	else if (hypothesis == 1) std::cout << "   H_0 : Motif count does not depend on event types.\n";
	else std::cout << "   H_0 : Motif count does not depend on node or event types.\n";
      }

    // Print the value of the rng seed.
    if (verbose && time_shuffling) 
      {
	std::cout << "   Shuffling event times (seed " << rng_seed << ")\n";
	if (bias_strength > 1) std::cout << "      Shuffling with bias strength " << bias_strength << ".\n";
      }
    if (verbose && edge_type_shuffling) std::cout << "   Shuffling edge types (seed " << rng_seed << ")\n";
    if (verbose && node_type_shuffling) std::cout << "   Shuffling node types (seed " << rng_seed << ")\n";

    return true;
  };

public:
  // Required parameters.
  unsigned int tw;
  std::string output_file_trunk;

  // Parameters constructed from the required parameters.
  std::string output_file_name;

  // Optional parameters.
  unsigned int max_size;
  bool maximal;
  unsigned int references;
  std::string node_file_name;
  unsigned int time_gap;
  double weight_omit;
  bool allow_multiple_event_types;
  unsigned int hypothesis;
  bool time_shuffling;
  unsigned int bias_strength;
  bool edge_type_shuffling;
  bool node_type_shuffling;
  unsigned int rng_seed;

  // Constructor sets default values for optional parameters.
  Parameters(bool verbose):
    verbose(verbose),
    max_size(0),
    maximal(false),
    references(0),
    node_file_name(),
    time_gap(0),
    weight_omit(0.0),
    allow_multiple_event_types(false),
    hypothesis(0),
    time_shuffling(false),
    bias_strength(1),
    edge_type_shuffling(false),
    node_type_shuffling(false),
    rng_seed(time(NULL)) 
  {};

  bool Init(int argc, char *argv[])
  {
        
    // Print help text or licence if there are not enough parameters.
    if (argc < 1 + N_required_param)
      {
        if (argc == 2) {
            std::string name(argv[1]);
            if (name.compare("--licence") == 0) {
                print_licence();
                return false;
            }
        }
	    print_help();
	    return false;
      }

    // Read the required parameters (time window and max motif file name).
    tw = atoi(argv[1]);
    output_file_trunk = argv[2];

    // Read optional parameters.
    int i = 3;
    while(i < argc)
      {
	if (!get_parameter(i, argc, argv)) return false;
	i++;
      }
    return post_process();
  };

};


/* Get all motifs and use them to fill locationMap.
 */
bool get_motifs(EdgeVectorMap& locationMap, 
		const Events& events,
		const Parameters& param,
		std::vector<unsigned short int> const& node_types)
{
  unsigned int gap_0 = events.first_time() + param.time_gap;
  unsigned int gap_1 = events.last_start_time() - param.time_gap;

  ProgressCounter evCounter(std::cerr, events.size(), 10);
  for (Events::const_iterator e_it = events.begin(); e_it != events.end(); ++e_it)
    {
      // Run through the first time_gap events and break if
      // the last time gap has been reached.
      if (e_it->start_time() < gap_0) continue;
      if (e_it->start_time() > gap_1) break;

      // Print progress.
      evCounter.next(*e_it);

      // Iterate through all subgraphs where the current event is the
      // first one, and update the count of the corresponding motif at
      // that location.
      TSubgraphFinder sgf(e_it->id(), param.tw, param.max_size, events, node_types);
      for (TSubgraphFinder::iterator sit = sgf.begin(); sit != sgf.end(); ++sit)
	{
	  const TSubgraph& sg = **sit;
	  if (sg.is_valid()) update_location_count(sg, locationMap);
	}
    }
  return true;
}

/* Get maximal motifs and use them to fill locationMap.
 */
bool get_maximal_motifs(EdgeVectorMap& locationMap, 
			const Events& events,
			const Parameters& param,
			std::vector<unsigned short int> const& node_types)
{
  unsigned int gap_0 = events.first_time() + param.time_gap;
  unsigned int gap_1 = events.last_start_time() - param.time_gap;

  // Maximal subgraphs are detected by inverting the maximal motif ids
  // so that we get the set of events that correspond to each maximal
  // subgraph id.
  std::map<event_id, EventSet> maximal_subgraphs;
  ProgressCounter evCounter(std::cerr, events.size(), 10);
  for (Events::const_iterator e_it = events.begin(); e_it != events.end(); ++e_it)
    {
      // Run through the first time_gap events and break if
      // the last time gap has been reached.
      if (e_it->start_time() < gap_0) continue;
      if (e_it->start_time() > gap_1) break;

      // Print progress.
      evCounter.next(*e_it);

      // Add event to the maximal subgraph.
      maximal_subgraphs[e_it->component()].insert(e_it->id());
    }

  // Now we know the exact set of events in each maximal
  // subgraph, so we just need to construct the corresponding
  // motifs.
  std::map<event_id, EventSet>::const_iterator ms_it;
  for (ms_it = maximal_subgraphs.begin(); ms_it != maximal_subgraphs.end(); ms_it++)
    {
      // Skip maximal subgraphs that are too large.
      if (param.max_size && ms_it->second.size() > param.max_size) continue;
      // Create temporal sugraph and update the location count.
      TSubgraph sg(events, ms_it->second, node_types, param.tw);
      if (sg.is_valid()) update_location_count(sg, locationMap);
    }
  return true;
}

int main(int argc, char *argv[])
{
  // Read command line parameters.
  Parameters param(true);
  if (!param.Init(argc, argv)) exit(1);
  std::cout << std::endl;

  // Initialize RNG.
  srand(param.rng_seed);

  // Read in the events.
  std::cerr << "Reading events from stdin ...\n";
  Events events(std::cin);

  // Try to read in the node types.
  std::vector<unsigned short int> node_types(events.get_nof_nodes());
  if (!param.node_file_name.empty())
    {
      unsigned int types_read = read_node_types(node_types, param.node_file_name);
      if (types_read)
        {
	  unsigned int max_node_index = node_types.size()-1;
	  std::vector<unsigned short int>::const_reverse_iterator rit = node_types.rbegin();
	  while (rit != node_types.rend() && *rit == 0) 
            {
	      ++rit; max_node_index--;
            }
	  std::cout << "   Read the type of " << types_read << " nodes (max index " << max_node_index << ").\n";
        }
      else 
        {
	  std::cout << "   Failed to read node types from file '" << param.node_file_name << "'.\n";
	  exit(1);
        }
    }
  else std::cout << "Only one type (0) of nodes used.\n";

  // Shuffle event times and/or node types.
  if (param.time_shuffling)
    {
      unsigned int shuffle_multiplier = 10;
      std::cout << "Shuffling (" << shuffle_multiplier << " x N_events";
      if (param.bias_strength > 1) std::cout << " with bias " << param.bias_strength;
      std::cout << ")...\n" << std::flush;

      if (param.bias_strength > 1) events.shuffle_constrained_corr(shuffle_multiplier, param.bias_strength);
      else events.shuffle_constrained(shuffle_multiplier);
    }
  if (param.node_type_shuffling)
    {
      std::random_shuffle(node_types.begin(), node_types.end());
    }
  if (param.edge_type_shuffling) 
    {
      if (!events.shuffle_edge_types())
        {
	  std::cerr << "Error: Unable to shuffle edge types because there were multiple event types on some edge.\n";
	  exit(1);
        }
    }

  unsigned int gap_0 = events.first_time() + param.time_gap;
  unsigned int gap_1 = events.last_start_time() - param.time_gap;

  // EVENT TYPES: Create a set for saving event types in the data.
  std::set<short int> eventTypes;

  // Construct the weighted, directed aggregate network.
  std::cerr << "Constructing aggregate network.\n";
  std::cout << "Constructing aggregate network ("<< currentDateTime() <<").\n";
  NetType net;
  std::map<short int, NetType*> nets;
  for (Events::iterator e_it = events.begin(); e_it != events.end(); ++e_it)
    {
      // Run through the first time_gap events and break if
      // the last time gap has been reached.
      if (e_it->start_time() < gap_0) continue;
      if (e_it->start_time() > gap_1) break;

      // Print progress.
      //evCounter.next(*e_it);
      net[e_it->from()][e_it->to()] += 1;

      // Update the set of event types.
      eventTypes.insert(e_it->type());

      // Increase count at the typed network.
      std::map<short int, NetType*>::iterator m_it = nets.find(e_it->type());
      if (m_it == nets.end()) 
        {
	  nets[e_it->type()] = new NetType;
	  m_it = nets.find(e_it->type());
        }
      NetType& typed_net = *(m_it->second);
      typed_net[e_it->from()][e_it->to()] += 1;
    }

  // Make all typed nets large enough to contain all nodes.
  for (std::map<short int, NetType*>::iterator m_it = nets.begin();
       m_it != nets.end(); ++m_it)
    {
      NetType& typed_net = *(m_it->second);
      typed_net.resize(net.size());
    }

  // Find the maximal subgraph ids of each event.
  // This is needed to properly detect motifs.
  std::cerr << "Finding maximal subgraphs.\n";
  std::cout << "Finding maximal subgraphs ("<< currentDateTime() <<").\n"; 
  events.find_maximal_subgraphs(param.tw);

  // Create maps for counting motifs by location.
  EdgeVectorMap locationMap;

  // ***************************
  // *** FILL IN locationMap ***
  // ***************************
  // Find the number of motifs at each location where there is a motif.
  if (param.maximal)
    {
      std::cerr << "Finding maximal typed motifs in data.\n";
      get_maximal_motifs(locationMap, events, param, node_types);      
    }
  else
    {
      std::cerr << "Finding typed motifs in data.\n";
      get_motifs(locationMap, events, param, node_types);
    }
  // Now 'get_location_count(locationMap, edges)' gives the number
  // of motifs at location specified by 'edges'. Note that 'edges'
  // includes information about the event types. Note that the edge
  // sequence uniquely gives the node types at this location, so
  // there is a unique typed hash to which this count corresponds
  // to.

  //std::cout << locationMap << std::endl;

  // Construct binner limits.
  std::set<unsigned int> bin_limits;
  get_limits_unbinned(net, bin_limits, param.weight_omit);

  // Object for counting motifs, both empirical and expected. The empirical counts 
  // should be placed at position 0, the references at position from 1 to param.references.
  ReferenceMotifCounter<double> motif_counts(param.references);

  // Create maps for counting the number of motifs by edge weights.
  // weightsMap[untyped_hash] is a binner instance.
  WeightsMap weightsMap;

  // Get the event type sequences that we go through next.
  TypeSeqsMap event_type_seqs;
  fill_event_type_seq(event_type_seqs, param.max_size, eventTypes, param.allow_multiple_event_types);
  std::cout << "Combinations of event types to go through per number of events:" << std::endl;
  std::cout << "   " << event_type_seqs << std::endl;

  // Whether to use node or event types when calculating the hash.
  bool use_node_types = (param.hypothesis == 1);
  bool use_event_types = (param.hypothesis == 0);


  // ***************************
  // *** FILL IN weightsMap ***
  // ***************************

  std::cerr << "Finding all locations with at most "
	    << param.max_size <<" edges in aggregate network.\n";
  std::cout << "Finding all locations with at most "
	    << param.max_size <<" edges in aggregate network ("<< currentDateTime() <<").\n"; 
  ProgressCounter pcounter(std::cerr, 1000000);
  for (SubnetIterator sn_it(net, param.max_size); !sn_it.finished(); ++sn_it)
    {
      // Print progress.
      pcounter.next(*sn_it);

      const NodepairVector& nodePairs = *sn_it;

      // Iterate through all assignments of event types on these edges.
      for (TypeSeqs::const_iterator ets_it = event_type_seqs[nodePairs.size()].begin();
	   ets_it != event_type_seqs[nodePairs.size()].end(); ++ets_it)
        {
	  // Get the untyped motif of the temporal subgraph with
	  // events on the given edges.
	  EdgeVector edges;
	  create_edges(edges, nodePairs, *ets_it);

	  TSubgraph sg(edges, node_types);
	  Motif* motif_untyped = sg.get_motif(use_node_types, use_event_types, false);
	  unsigned int untyped_hash = motif_untyped->get_hash();
	  delete motif_untyped;

	  // Get the weight sequence of edges. Continue if some edge
	  // has zero weight (this is possible because we are
	  // detecting the subgraphs in a network where all event
	  // types have been aggregated.
	  std::vector<unsigned int> curr_weights;
	  if (!get_edge_weights(edges, nets, curr_weights)) continue;

	  // Get the binner for this motif, and initialize it if one didn't exist.
	  wBinner& curr_binner = weightsMap[untyped_hash];
	  if (!curr_binner.is_initialized()) curr_binner.Init(bin_limits, edges.size());

	  // Increase the binner at index given by weights by a value given
	  // by the number of this motif at this exact location.
	  unsigned int location_count = get_location_count(locationMap, edges);
	  if (curr_binner.add(curr_weights, location_count))
            {
	      // The count at this location was successfully added, which means that 
	      // the weights at this location are included in statistics. Increase the
	      // count of the typed motif also.
	      Motif* motif_typed = sg.get_typed_motif();
	      motif_counts.add_at(motif_typed, 0, location_count);
            }

	  /* // DEBUG
	     unsigned int sum, count;
	     weightsMap[untyped_hash].get_sum(curr_weights, sum);
	     weightsMap[untyped_hash].get_count(curr_weights, count);
	     std::cerr << *sn_it << " " << curr_weights << " (" << sum << ", " << count << ")" << std::endl;
	  */
	}
    }
  // weightsMap[motif_hash].get_random(edge_weights) now gives a
  // random sample from the distribution of motif counts at
  // locations with given weights sequence.

  // Note that the 'motif_hash' used as key in weightsMap defines
  // the reference system. If the hash corresponds to the untyped
  // motif (no event or node types), the null hypothesis is "Node
  // and event types do not affect motif counts." This is the only
  // option if the data has only node types or only event types.

  // If however the node has both node and event types, there are
  // more alternatives. We can then test against the null hypothesis
  // "Event types do not affect motif counts" (when motif hash is
  // obtained by omitting event types) or "Node types do not affect
  // motif counts" (when motif hash is obtained by omitting node
  // types).

  // ***************************
  // *** SAMPLE MOTIF COUNTS ***
  // ***************************
  // Iterate again over all connected sequences of at most
  // param.max_size edges. At each location get a sample from the
  // distribution corresponding to that untyped motif and weight
  // sequence.
  std::cerr << "Calculating expected number of each motif.\n";
  std::cout << "Calculating expected number of each motif ("<< currentDateTime() <<").\n"; 
  pcounter.reset();
  for (SubnetIterator sn_it(net, param.max_size); !sn_it.finished(); ++sn_it)
    {
      // Print progress.
      pcounter.next(*sn_it);

      const NodepairVector& nodePairs = *sn_it;

      // Iterate through all assignments of event types on these edges.
      for (TypeSeqs::const_iterator ets_it = event_type_seqs[nodePairs.size()].begin();
	   ets_it != event_type_seqs[nodePairs.size()].end(); ++ets_it)
	{
	  // Get the untyped motif of the temporal subgraph with
	  // events on the given edges.
	  EdgeVector edges;
	  create_edges(edges, nodePairs, *ets_it);

	  TSubgraph sg(edges, node_types);
	  Motif* motif_untyped = sg.get_motif(use_node_types, use_event_types, false);
	  unsigned int untyped_hash = motif_untyped->get_hash();
	  delete motif_untyped;

	  // Get the weight sequence of edges. Continue if some edge
	  // has zero weight (this is possible because we are
	  // detecting the subgraphs in a network where all event
	  // types have been aggregated.
	  std::vector<unsigned int> curr_weights;
	  if (!get_edge_weights(edges, nets, curr_weights)) continue;

	  // Get a random number of this motif given the edge weights at
	  // this location for each reference.
	  std::vector<unsigned int> ref_counts(param.references);
	  if (weightsMap[untyped_hash].get_random(curr_weights, ref_counts))
	    { 
	      // The weight sequence is included in the statistics.

	      /* // DEBUG
		 unsigned int sum, count;
		 weightsMap[motif_ut_hash].get_sum(curr_weights, sum);
		 weightsMap[motif_ut_hash].get_count(curr_weights, count);
		 std::cout << *sn_it << " " << curr_weights << " (" << sum << ", " << count << ") = " << mean_count << std::endl;
	      */

	      // Get the typed motif spanned by these events.
	      Motif* motif_typed = sg.get_typed_motif();

	      // Add counts to the reference value of the typed motif (if
	      // non-zero).
	      unsigned int i_ref = 1;
	      for (std::vector<unsigned int>::const_iterator ref_it = ref_counts.begin();
		   ref_it != ref_counts.end(); ++ref_it)
		{
		  if (*ref_it) motif_counts.add_at(motif_typed, i_ref, *ref_it);
		  ++i_ref;
		}
	    }
	}
    }

  // Print out the results.
  std::cout << "Calculations finished ("<< currentDateTime() <<")." << std::endl;
 if (motif_counts.print(param.output_file_name)) std::cout << "Results written ("<< currentDateTime() <<")." << std::endl;
  else std::cout << "Error writing results to file! ("<< currentDateTime() <<")." << std::endl;

  // Free nets.
  for (std::map<short int, NetType*>::iterator m_it = nets.begin();
       m_it != nets.end(); ++m_it) 
    {
      //std::cerr << "Delete net with event type " << m_it->first << std::endl;
      delete m_it->second;
    }

}
