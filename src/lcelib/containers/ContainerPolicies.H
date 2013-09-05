#ifndef LCE_CONTAINER_POLICIES
#define LCE_CONTAINER_POLICIES

/**
 * The standard specializations for the policies for containers 
 * is defined here. This refers to the containers as whole, not to 
 * anytyhing contained. That is, if containers are contained, this is 
 * what is used. Uses CRTP.
 */

/**
 * This is the base for weight policies of weighed tables, that is, ones
 * for which you can get the weight out without going through all the 
 * elements.
 */

template<typename TableType>
struct WeightWeightPolicy {
  typedef typename TableType::WeightType WeightType;
  static WeightType getWeight(const TableType & subj) {
    return subj.weight();
  }
};

/**
 * For unweighed tables, this is what we probably mean...
 */

template<typename TableType>
struct CountWeightPolicy {
  /** 
   * If we do not store weights _explicitly_ we probably want to use
   * the number of elements as the weight. We do not give the direct
   * reference out, though. */
  typedef size_t WeightType;

  static WeightType getWeight(const TableType & subj) {
    return subj.elemsByCRTP();
  }
};





#endif
