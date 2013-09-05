/*
Save and print the number of motifs in data.

Lauri Kovanen (2/2012)
*/

#ifndef MOTIFCOUNTER_H
#define MOTIFCOUNTER_H

#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <iomanip>
#include "assert.h"
#include "std_printers.h"
#include "motif.h"

template <typename T> struct MotifCount
{
  unsigned int count;
  std::string* desc;
  std::vector<T> ref_counts;
  MotifCount();
  ~MotifCount();
};

template <typename T> MotifCount<T>::MotifCount() 
: count(0),desc(NULL),ref_counts() {}

template <typename T> MotifCount<T>::~MotifCount() 
{
  delete desc;
}

template <class T> 
std::ostream& operator<<(std::ostream& output, const MotifCount<T>& mc) 
{
  output << mc.count << "/" << mc.ref_counts;
  return output;
}

// Abstract base class for counting motifs.
template <typename T> class MotifCounter
{
 private:
  // Map for saving motif counts. The key is the motif hash, the value
  // is a pair that consists of motif counts and a string description of
  // the motifs.
  typedef std::map<unsigned int, MotifCount<T> > MotifCountMap;

  MotifCountMap mcm;

 protected:
  unsigned int N; // Number of other values in addition to empirical one.

  // Sort motifs by count in data.
  void sort_by_count(std::multimap<unsigned int, unsigned int>& sorted_counts) const;
  void sort_by_ref_count(std::multimap<T, unsigned int>& sorted_counts) const;

  // Access results.
  inline unsigned int count(unsigned int hash) const { return mcm.find(hash)->second.count; };
  inline std::vector<T> const& counts(unsigned int hash) const { return mcm.find(hash)->second.ref_counts; };
  inline std::string const& desc(unsigned int hash) const { return *(mcm.find(hash)->second.desc); };

 public:
  // Really simple constructor.
  MotifCounter(unsigned int N);

  // Increase count of motif m at position i by 1.
  virtual void increment_at(Motif *m, unsigned int i) { add_at(m,i,1); };

  // Increase count of motif m at position i by val.
  virtual void add_at(Motif *m, unsigned int i, T val);

  // Print output.
  virtual bool print(const std::string& fileName) const =0;

  void debug_print() const { std::cerr << mcm << std::endl; };
};

template <typename T> class ReferenceMotifCounter : public MotifCounter<T>
{
 protected:
  // Save the number of locations where the motif occurs (index 1 for 
  // empirical data and 2-> for references) and the total number of
  // locations in the aggregate network (index 0). Key is the motif hash.
  std::map<unsigned int, std::vector<unsigned int> > locationCounts;

 public:
  ReferenceMotifCounter(unsigned int N_ref);
  void add_at(Motif *m, unsigned int i, T val);
  bool print(const std::string& fileName) const;
};

template <typename T> class SingleRefMotifCounter : public MotifCounter<T>
{
 public:
  SingleRefMotifCounter();
  bool print(const std::string& fileName) const;
};

template <typename T> class DistributionMotifCounter : public MotifCounter<T>
{
 public:
  DistributionMotifCounter(unsigned int tw);
  void increment_at(Motif *m, unsigned int i) { MotifCounter<T>::increment_at(m, i+1); };
  bool print(const std::string& fileName) const;
};


template<typename T>
MotifCounter<T>::MotifCounter(unsigned int N) : mcm(), N(N) {}

template<typename T>
void MotifCounter<T>::add_at(Motif *m, unsigned int i, T value)
{
  MotifCount<T>& mc = mcm[m->get_hash()];
  if (mc.desc == NULL) 
    {
      mc.desc = new std::string(to_string(*m));
      mc.ref_counts.resize(N);
    }
  if (i == 0) mc.count += value;
  else mc.ref_counts[i-1] += value;
}

template<typename T>
void MotifCounter<T>::sort_by_count(std::multimap<unsigned int, unsigned int>& sorted_counts) const
{
  typename MotifCountMap::const_iterator m_it;
  for (m_it = mcm.begin(); m_it != mcm.end(); ++m_it)
    {
      sorted_counts.insert(std::make_pair(m_it->second.count,m_it->first));
    }
}

template<typename T>
void MotifCounter<T>::sort_by_ref_count(std::multimap<T, unsigned int>& sorted_counts) const
{
  typename MotifCountMap::const_iterator m_it;
  for (m_it = mcm.begin(); m_it != mcm.end(); ++m_it)
    {
      T total_count = 0;
      for (typename std::vector<T>::const_iterator v_it = m_it->second.ref_counts.begin();
	   v_it != m_it->second.ref_counts.end(); v_it++) total_count += *v_it;
      sorted_counts.insert(std::make_pair(total_count, m_it->first));
    }
}

template<typename T>
ReferenceMotifCounter<T>::ReferenceMotifCounter(unsigned int N_ref) 
: MotifCounter<T>(N_ref),locationCounts()
{}

template<typename T>
void ReferenceMotifCounter<T>::add_at(Motif *m, unsigned int i, T value)
{
  MotifCounter<T>::add_at(m,i,value);

  // Increment the location count.
  std::vector<unsigned int>& v = locationCounts[m->get_hash()];
  if (v.empty()) v.resize((this->N)+2);
  if (i == 0) v[0]++;
  if (value > 0) v[i+1]++;
}

template <typename T>
bool ReferenceMotifCounter<T>::print(const std::string& fileName) const
{
  /* Sort the map keys by count in the actual data so they can be
     easily printed in sorted order. (The multimap is automatically
     sorted by its key.) */
  std::multimap<unsigned int, unsigned int> sorted_counts;
  this->sort_by_count(sorted_counts);
  
  std::cout << "Writing results to file " << std::endl;
  std::cout << "   " << fileName.c_str() << std::endl;

  /* Open output stream. */
  std::ofstream output;
  output.open(fileName.c_str());
  if (output.fail())
    {
      perror("Failed to open output file");
      return false;
    }
	
  /* Print output. */
  unsigned int fw = 14; // field width
  unsigned int fp = 4; // field precision
  output << "count     "
	 << "ref_avg       "
	 << "ratio         "
	 << "ref_std       "
	 << "N_ref "
	 << "N_lt_ref "
	 << "z-score       "
	 << "N_loc_tot   "
	 << "N_loc       "
	 << "N_loc_ref     "
	 << "N [node:color ...] edges ..." << std::endl;
	
  std::multimap<unsigned int, unsigned int>::reverse_iterator s_it;
  for (s_it = sorted_counts.rbegin(); s_it != sorted_counts.rend(); ++s_it)
    {
      unsigned int h = s_it->second; // Motif hash.
		
      unsigned int data_count = this->count(h);
      const std::vector<T>& ref_counts = this->counts(h);
      unsigned int rnd_found_count = 0;
      unsigned int N_lt_ref = 0;
      double avg_count = 0, std_count = 0;
      for (typename std::vector<T>::const_iterator v_it = ref_counts.begin(); v_it != ref_counts.end(); ++v_it)
	{
	  if (*v_it > 0) {
	    double c = (double)(*v_it);
	    avg_count += c;
	    std_count += c*c;
	    rnd_found_count += 1;
	  }
	  if (data_count < *v_it) N_lt_ref += 1;
	}
		
      avg_count /= this->N;
      std_count = std_count/this->N - avg_count*avg_count;
      std_count = sqrt(std_count);
      double z_score = 0.0;
      if (std_count > 0) z_score = (((double)data_count) - avg_count)/std_count;
      double ratio = -1.0;
      if (avg_count > 0) ratio = ((double)data_count)/avg_count;

      const std::vector<unsigned int>& v = locationCounts.find(h)->second;
      std::vector<unsigned int>::const_iterator lv_it = v.begin();
      unsigned int N_loc_tot = *lv_it; ++lv_it;
      unsigned int N_loc = *lv_it; ++lv_it;
      double N_loc_ref = 0;
      for (; lv_it != v.end(); ++lv_it) N_loc_ref += *lv_it;
      N_loc_ref /= this->N;
		
      output << std::setiosflags(std::ios::left)
	     << std::setw(10) << data_count
	     << std::setw(fw) << std::setiosflags(std::ios::fixed) << std::setprecision(2) << avg_count
	     << std::setw(fw) << std::setiosflags(std::ios::fixed) << std::setprecision(fp) << ratio
	     << std::setw(fw) << std::setiosflags(std::ios::fixed) << std::setprecision(fp) << std_count
	     << std::setw(6) << rnd_found_count
	     << std::setw(9) << N_lt_ref
	     << std::setw(fw) << std::setiosflags(std::ios::fixed) << std::setprecision(fp) << z_score
	     << std::setw(12) << N_loc_tot
	     << std::setw(12) << N_loc
	     << std::setw(fw) << std::setiosflags(std::ios::fixed) << std::setprecision(2) << N_loc_ref
	     << this->desc(h) << std::endl;
    }
	
  output.close();
  if (output.fail())
    {
      perror("Failed to close output file");
      return false;
    }
  return true;
}


template<typename T>
SingleRefMotifCounter<T>::SingleRefMotifCounter() 
  : MotifCounter<T>(1)
{}

template <typename T>
bool SingleRefMotifCounter<T>::print(const std::string& fileName) const
{
  /* Sort the map keys by count in the actual data so they can be
     easily printed in sorted order. (The multimap is automatically
     sorted by its key.) */
  std::multimap<unsigned int, unsigned int> sorted_counts;
  this->sort_by_count(sorted_counts);
  
  /* Open output stream. */
  std::ofstream output;
  output.open(fileName.c_str());
  if (output.fail())
    {
      perror("Failed to open output file");
      return false;
    }
	
  /* Print output. */
  unsigned int fw = 14; // field width
  unsigned int fp = 4; // field precision
  output << "count     "
	 << "ref           "
	 << "ratio         "
	 << "N [node:color ...] edges ...\n";
	
  std::multimap<unsigned int, unsigned int>::reverse_iterator s_it;
  for (s_it = sorted_counts.rbegin(); s_it != sorted_counts.rend(); ++s_it)
    {
      unsigned int h = s_it->second; // Motif hash.
		
      unsigned int data_count = this->count(h);
      T ref_count = this->counts(h)[0];
      double ratio = -1.0;
      if (ref_count > 0) ratio = ((double)data_count)/ref_count;
		
      output << std::setiosflags(std::ios::left)
	     << std::setw(10) << data_count
	     << std::setw(fw) << std::setiosflags(std::ios::fixed) << std::setprecision(4) << ref_count
	     << std::setw(fw) << std::setiosflags(std::ios::fixed) << std::setprecision(fp) << ratio
	     << this->desc(h) << std::endl;
    }
	
  output.close();
  if (output.fail())
    {
      perror("Failed to close output file");
      return false;
    }
  return true;
}


template<typename T>
DistributionMotifCounter<T>::DistributionMotifCounter(unsigned int tw) 
  : MotifCounter<T>(tw+1)
{}

template<typename T>
bool DistributionMotifCounter<T>::print(const std::string& fileName) const 
{
  /* Sort the map keys by total count so they can be easily printed in
     sorted order. (The multimap is automatically sorted by its
     key.) */
  std::multimap<T, unsigned int> sorted_counts;
  this->sort_by_ref_count(sorted_counts);

  /* Open output stream. */
  std::ofstream output;
  output.open(fileName.c_str());
  if (output.fail())
    {
      perror("Failed to open output file");
      return false;
    }

  /* Print output. */
  typename std::multimap<T, unsigned int>::reverse_iterator s_it;
  for (s_it = sorted_counts.rbegin(); s_it != sorted_counts.rend(); ++s_it)
    {
      unsigned int h = s_it->second; // Motif hash.

      // Print counts separated by commas.
      const std::vector<T>& ref_counts = this->counts(h);      
      std::vector<unsigned int>::const_iterator v_it = ref_counts.begin();
      output << *v_it; v_it++;
      for (; v_it != ref_counts.end(); v_it++) output << "," << *v_it;

      // Print the motif string.
      output << " " << this->desc(h) << std::endl;
    }
	
  output.close();
  if (output.fail())
    {
      perror("Failed to close output file");
      return false;
    }
  return true;
}

#endif
