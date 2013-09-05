#ifndef LCE_PAIR
#define LCE_PAIR
#include<cstring>
#include<memory>

/**
 * An utility pair class for containers. Partial specializations for void
 * members exist: in order to simplify the design, they return the non-void
 * member if a reference for the void one is called. For safety, only 
 * const references are given out this way, however.
 * 
 * Includes typedefs in order to avoid refs to void. 
 */

template<typename FirstType, typename SecondType>
class Pair {
protected:
  FirstType _first;
  SecondType _second;
  
public:
  typedef FirstType & first_reference;
  typedef SecondType & second_reference;
  typedef const FirstType & const_first_reference;
  typedef const SecondType & const_second_reference;
  FirstType & first() {return _first;}
  SecondType & second() {return _second;}
  const FirstType & first() const {return _first;}
  const SecondType & second() const {return _second;}

  void removeSecond() {_second.~SecondType();}
  void clearSecond() {
    SecondType temp=SecondType();
    memcpy(&(_second),
	   &temp,
	   sizeof(SecondType));
  }
};

template<typename FirstType>
class Pair<FirstType, void> {
  FirstType _first;

public:
  typedef FirstType & first_reference;
  typedef const FirstType & second_reference;
  typedef const FirstType & const_first_reference;
  typedef const FirstType & const_second_reference;
  FirstType & first() {return _first;}
  const FirstType & first() const {return _first;}
  const FirstType & second() const {return _first;}

  void removeSecond() {}
};

template<typename SecondType>
class Pair<void, SecondType> {
  SecondType _second;

public:
  typedef const SecondType & first_reference;
  typedef SecondType & second_reference;
  typedef const SecondType & const_first_reference;
  typedef const SecondType & const_second_reference;
  SecondType & second() {return _second;}
  const SecondType & second() const {return _second;}
  const SecondType & first() const {return _second;}

  void removeSecond() {_second.~SecondType();}
  
  void clearSecond() {
    SecondType temp=SecondType();
    memcpy(&(_second),
	   &temp,
	   sizeof(SecondType));
  }
};
  

#endif
