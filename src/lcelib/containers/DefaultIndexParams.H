#ifndef LCE_DEFAULT_INDEX_PARAMS
#define LCE_DEFAULT_INDEX_PARAMS




template<template<typename, 
		  typename, 
		  typename, 
		  typename, 
		  template<typename, 
			   typename, 
			   typename,
			   typename,
			   typename> class Foo> class Index>
struct DefaultIndexParams {
			     /** This is not cared about if no status is used. Void means auto.*/		    
			     typedef typename EmbStatusPolicy UsagePolicy;
			   };

#endif
