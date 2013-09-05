/*
Binner for binning data.

2/2012: The bin limits must be integers; only some methods have been
implemented, the code is far from complete.

Lauri Kovanen (2/2012)
*/

#ifndef BINNER_H
#define BINNER_H

#include <map>
#include <set>
#include <vector>
#include <algorithm>
#include <assert.h>
#include "std_printers.h"

/* Class: Binner

   A multidimensional binner object. The bin limits are identical in
   each dimension. The template argument gives the data type.

   The bin limits must be unsigned int (for now).
*/
typedef unsigned int limit_type;
typedef unsigned int index_type;
template<typename value_type> class Binner
{
  typedef std::map<limit_type, index_type> BinLimits;
  typedef std::map<value_type, unsigned int> ValueDist; // value, count
  typedef std::map<std::vector<index_type>, std::pair<ValueDist, unsigned int> > DataMap;

 private:
  BinLimits bin_limits;
  unsigned int dim;
  limit_type min_val, max_val;
  DataMap data;

  /* Auxiliary functions used by find() to find indices for the
     position. Both return false if position is outside bin limits. */
  bool get_index(limit_type pos, index_type& index);
  bool get_indices(const std::vector<limit_type>& pos, std::vector<index_type>& indices);

  /* Find the distribution at given position. Returns false if the
     position is invalid (outside bin limits). If the value is not
     found, it = data.end(). */
  bool find(const std::vector<limit_type>& pos, typename DataMap::iterator& it);

  /* Extract sum from value distribution. */
  value_type get_sum(const ValueDist& valueDist);

  /* Convert distribution to cumulative distribution. */
  void create_cumulative();

  /* Return a single random value from distribution. */
  value_type get_random(const ValueDist& valueDist, unsigned int count);

 public:
  Binner();
  bool Init(const std::set<limit_type>& limits, unsigned int dim);
  bool is_initialized() { return dim > 0; };

  /* Add a value to distribution. */
  bool add(const std::vector<limit_type>& pos, value_type val);

  /* Both get_mean() and get_sum() return false if pos is not within
     bin limits or if no values corresponding to pos have been
     added. */
  bool get_mean(const std::vector<limit_type>& pos, double& res);
  bool get_sum(const std::vector<limit_type>& pos, value_type& res);

  /* get_count() returns false only if pos is outside bin limits. */
  bool get_count(const std::vector<limit_type>& pos, unsigned int& res);

  /* Get a random value(s) from distribution corresponding to
     pos. Returns false if pos is outside bin limits or the
     distribution is empty. */
  bool get_random(const std::vector<limit_type>& pos, value_type& res);
  bool get_random(const std::vector<limit_type>& pos, std::vector<value_type>& res);

  void print_data();
};

template<typename value_type> Binner<value_type>::Binner() 
: bin_limits(), dim(0), data() {}

template<typename value_type> bool Binner<value_type>::Init(const std::set<limit_type>& limits, unsigned int dim)
{
  this->dim = dim;
  index_type i = 0;
  for (std::set<limit_type>::const_iterator it = limits.begin(); it != limits.end(); ++it)
    {
      bin_limits[*it] = i;
      i++;
    }
  min_val = *limits.begin();
  max_val = *(--limits.end());
  //std::cerr << "bin_limits: " << bin_limits << " (" << min_val << ", " << max_val << ")\n";
  return true;
}

template<typename value_type> 
bool Binner<value_type>::add(const std::vector<limit_type>& pos, value_type val)
{
  std::vector<index_type> indices(dim);
  if (!get_indices(pos, indices)) return false;
  std::pair<ValueDist, unsigned int>& data_pair = data[indices];
  data_pair.first[val]++;
  data_pair.second++;
  return true;
}

template<typename value_type> 
bool Binner<value_type>::get_indices(const std::vector<limit_type>& pos, std::vector<index_type>& indices)
{
  //std::cerr << pos << std::endl;
  //std::cerr << indices << std::endl;
  assert(pos.size() == indices.size());
  std::vector<index_type>::iterator index = indices.begin();
  for (std::vector<limit_type>::const_iterator it = pos.begin(); it != pos.end(); ++it)
    {
      if (!get_index(*it, *index)) return false;
      ++index;
    }
  return true;
}

template<typename value_type> 
bool Binner<value_type>::get_index(limit_type pos, index_type& index)
{
  //std::cerr << "Looking for " << pos;
  if (pos < min_val || pos >= max_val) {
    //std::cerr << " (not found, "<< min_val << " " << max_val <<")\n";
    return false;
  }
  index = (bin_limits.upper_bound(pos))->second - 1;
  //std::cerr << " (found at "<<index<<")\n";
  return true;
}

template<typename value_type> 
bool Binner<value_type>::find(const std::vector<limit_type>& pos, typename DataMap::iterator& it)
{
  std::vector<index_type> indices(dim);
  if (!get_indices(pos, indices)) return false;
  it = data.find(indices);
  return true;
}

template<typename value_type> 
value_type Binner<value_type>::get_sum(const ValueDist& valueDist)
{
  value_type sum = 0;
  for (typename ValueDist::const_iterator vit = valueDist.begin(); vit != valueDist.end(); ++vit)
    {
      sum += vit->first * vit->second;
    }
  return sum;
}

template<typename value_type> 
bool Binner<value_type>::get_mean(const std::vector<limit_type>& pos, double& res)
{
  typename DataMap::iterator it;
  if (!find(pos, it) || it == data.end()) return false;
  res = ((double)get_sum(it->second.first))/((double)it->second.second);
  return true;
}

template<typename value_type> 
bool Binner<value_type>::get_count(const std::vector<limit_type>& pos, unsigned int& res)
{
  typename DataMap::iterator it;
  if (!find(pos, it)) return false;
  if (it == data.end()) res = 0.0;
  else res = it->second.second;
  return true;
}

template<typename value_type> 
bool Binner<value_type>::get_sum(const std::vector<limit_type>& pos, value_type& res)
{
  typename DataMap::iterator it;
  if (!find(pos, it) || it == data.end()) return false;
  res = get_sum(it->second.first);
  return true;
}

template<typename value_type> 
value_type Binner<value_type>::get_random(const ValueDist& valueDist, unsigned int count)
{
  unsigned int i = 1 + (unsigned int)(count*(rand()/(RAND_MAX+1.0)));
  unsigned int cum_count = 0;
  for (typename ValueDist::const_iterator vit = valueDist.begin(); vit != valueDist.end(); ++vit)
    {
      cum_count += vit->second;
      //std::cerr << "i "<<i<<", cum_count "<<cum_count<<", value "<<res<<std::endl;	  
      if (cum_count >= i) return vit->first;
    }
  return 0; // This should never happen.
}

template<typename value_type> 
bool Binner<value_type>::get_random(const std::vector<limit_type>& pos, value_type& res)
{
  typename DataMap::iterator it;
  if (!find(pos, it) || it == data.end()) return false;
  res = get_random(it->second.first, it->second.second);
  return true;
}

template<typename value_type> 
bool Binner<value_type>::get_random(const std::vector<limit_type>& pos, std::vector<value_type>& res)
{
  typename DataMap::iterator it;
  if (!find(pos, it) || it == data.end()) return false;
  for (typename std::vector<value_type>::iterator res_it = res.begin(); res_it != res.end(); ++res_it)
    {
      *res_it = get_random(it->second.first, it->second.second);
    }
  return true;
}

template<typename value_type> 
void Binner<value_type>::print_data()
{
  std::cerr << data << std::endl;
}


#endif
