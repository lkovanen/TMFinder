#ifndef STD_PRINTERS
#define STD_PRINTERS

#include <ostream>
#include <sstream>
#include <vector>
#include <list>
#include <stack>
#include <set>
#include <map>
#include <math.h>

/* 
   A collection of methods that allow printing C++ data structures.
 */

template <class T> 
void invert_vector(std::vector<T>& v) 
{
  size_t n = v.size();
  size_t i_mid = (unsigned int)floor(0.5*((float)n));
  for (size_t i = 0; i < i_mid; ++i)
    {
      size_t j = n-i-1;
      T i_copy = v[i];
      v[i] = v[j];
      v[j] = i_copy;
    }
}

/* Turn any printable object into a string. */
template <class T>
inline std::string to_string(const T& t)
{
  std::stringstream ss;
  ss << t;
  return ss.str();
}

template <class T> 
std::ostream& operator<<(std::ostream& output, const std::vector<T>& v) 
{
  std::string separator;
  typename std::vector<T>::const_iterator it;
  output << "[";
  for (it = v.begin(); it != v.end(); ++it)
    {
      output << separator << *it;
      separator = ", ";
    }
  output << "]";
  return output;

}

template <class T> 
std::ostream& operator<<(std::ostream& output, const std::set<T>& s)
{
  output << "{";
  std::string separator;
  typename std::set<T>::const_iterator it;
  for (it = s.begin(); it != s.end(); ++it)
    {
      output << separator << *it;
      separator = ", ";
    }
  output << "}";
  return output;
}

template <class K, class V> 
std::ostream& operator<<(std::ostream& output, const std::map<K,V>& m) 
{
  std::string separator;
  typename std::map<K,V>::const_iterator it;
  output << "{";
  for (it = m.begin(); it != m.end(); ++it)
    {
      output << separator << it->first << ":" << it->second;
      separator = ", ";
    }
  output << "}";
  return output;
}

template <class K, class V> 
std::ostream& operator<<(std::ostream& output, const std::multimap<K,V>& m) 
{
  std::string separator;
  typename std::map<K,V>::const_iterator it;
  output << "{";
  for (it = m.begin(); it != m.end(); ++it)
    {
      output << separator << it->first << ":" << it->second;
      separator = ", ";
    }
  output << "}";
  return output;
}

template <class V> 
std::ostream& operator<<(std::ostream& output, const std::list<V>& l) 
{
  std::string separator;
  typename std::list<V>::const_iterator it;
  for (it = l.begin(); it != l.end(); ++it)
    {
      output << separator << *it;
      separator = ", ";
    }
  return output;
}

template <class V> 
std::ostream& operator<<(std::ostream& output, const std::stack<V>& s) 
{
  // Stacks don't have iterators; make a copy of the stack first,
  // because we will need to destroy it to print contents.
  std::stack<V> s_cp(s);
  std::string separator;
  while (!s_cp.empty())
    {
      output << separator << s_cp.top();
      separator = " | ";
      s_cp.pop();
    }
  return output;
}

template <class V1, class V2> 
  std::ostream& operator<<(std::ostream& output, const std::pair<V1,V2>& p) 
{
  output << "(" << p.first << ", " << p.second << ")";
  return output;
}

#endif
