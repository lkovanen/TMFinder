#ifndef LCE_KEY_POLICY
#define LCE_KEY_POLICY

/** 
 * Default hash calculator merely xors the whole object together.
 * A nice candidate for loop unrolling: this should in fact be good
 * enough for anything, although endianness could matter. 
 */
template<typename KeyType>
struct KeyPolicy {
  static const size_t numWords=sizeof(KeyType)/sizeof(size_t);
  static const size_t numBytes=sizeof(KeyType)-numWords*sizeof(size_t);
  static size_t getHashValue(const KeyType & subj) {
    size_t val=0;
    /* These should be calculated, resolved & unwrapped at 
     * compile-time 
     */
    /* First, the part that is a multiple of the machine word size */
    for (size_t i=0; i<numWords; ++i) {
      val^=((size_t *)(&subj))[i];
    }     
    /* Then, rest. */
    for (size_t j=0; j<numBytes; ++j) {
      ((char *)(&val))[j] ^=
	((char *)(&subj))[j+numWords*sizeof(size_t)];
    }
    return val;
  }
};

/** Direct cast for pointers */

template<typename KeyType> 
struct KeyPolicy<KeyType *> {
  static size_t getHashValue(const KeyType * & key) {
    return static_cast<size_t>(key);
  }
};

template<typename IntegralType>
struct IntKeyPolicy {
  static size_t getHashValue(const IntegralType & key) {return key;}
  static bool firstSmaller(const IntegralType & first,
			   const IntegralType & second) {
    return first<second;
  }
};

/** Specializations for small integral types: */

template<> 
struct KeyPolicy<char>: public IntKeyPolicy<char> {};

template<> 
struct KeyPolicy<unsigned char>: public IntKeyPolicy<unsigned char> {};

template<> 
struct KeyPolicy<short>: public IntKeyPolicy<short> {};

template<> 
struct KeyPolicy<unsigned short>: public IntKeyPolicy<unsigned short> {};

template<> 
struct KeyPolicy<int>: public IntKeyPolicy<int> {};

template<> 
struct KeyPolicy<unsigned int>: public IntKeyPolicy<unsigned int> {};

/* size_t is unsigned int. Bugger. */

//template<> 
//struct KeyPolicy<size_t>: public IntKeyPolicy<size_t> {};


#endif
