/*
Print output during long iterations to follow progress.

Lauri Kovanen (3/2012)
*/

#include "progress_counter.h"

ProgressCounter::ProgressCounter(std::ostream& output, unsigned int iter_length, unsigned int N_output)
: N(iter_length),
  increment(((double)N_output)/((double)iter_length)),
  output(output)
{ reset(); }

ProgressCounter::ProgressCounter(std::ostream& output, unsigned int output_every)
: N(0),
  increment(1.0/((double)output_every)),
  output(output)
{ reset(); }

bool ProgressCounter::increment_counter()
{
  n++;
  counter += increment;
  if (counter > 1.0)
    {
      while (counter > 1.0) counter -= 1.0;
      return true;
    }
  return false;
}

void ProgressCounter::output_position()
{
  if (N)
    {
      output << std::setw(5) << std::setiosflags(std::ios::fixed) 
	     << std::setprecision(2) << 100*((double)(n))/N << "%";
    }
  else output << n;
}

void ProgressCounter::next()
{
  if (increment_counter()) output_position();
}


