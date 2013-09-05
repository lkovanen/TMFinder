#ifndef LCE_WEIGHT_POLICY
#define LCE_WEIGHT_POLICY

/** 
 * Policies for element weights in containers. Getting, and in the case
 * of implicit trees only storing sums, setting elements. 
 *
 * Don't worry about the setting part: it is only used with specific 
 * data structures. If you still do, just leave the correponding 
 * functions undefined. 
 *
 * You probably want to consider specializing the specialization
 * for container data types instead of this. 
 */

template<typename DataType> 
struct WeightPolicy {
  typedef DataType WeightType;
  /* The type of a reference to a weight */
  typedef DataType & weight_reference;
  
  static WeightType getWeight(const DataType & subj) {return subj;}
  
  static weight_reference refToWeight(DataType & subj) {
    return subj;
  }
};

#endif

