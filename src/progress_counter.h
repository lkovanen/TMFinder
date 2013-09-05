/*
Print output during long iterations to follow progress.

Lauri Kovanen (3/2012)
*/

#ifndef PROGRESSCOUNTER_H
#define PROGRESSCOUNTER_H

#include <ostream>
#include <iomanip>

/* Progress counter */
class ProgressCounter
{
 private:
  const unsigned int N;
  unsigned int n;
  const double increment;
  double counter;
  std::ostream& output;

  bool increment_counter();
  void output_position();

 public:
  ProgressCounter(std::ostream& output, unsigned int iter_length, unsigned int N_output);
  ProgressCounter(std::ostream& output, unsigned int output_every);
  template <class T> void next(const T& t);
  template <class T1, class T2> void next(const T1& t1, const T2& t2);
  void next();
  void reset() { n=0; counter=increment; };
};

template <class T>
void ProgressCounter::next(const T& t)
{
  if (increment_counter())
    {
      output_position();
      output << " : " << t << std::endl;
    }
}

template <class T1, class T2>
void ProgressCounter::next(const T1& t1, const T2& t2)
{
  if (increment_counter())
    {
      output_position();
      output << " : " << t1 << ", " << t2 << std::endl;
    }
}

#endif
