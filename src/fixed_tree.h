/*
Multiset implementation with a binary tree that uses an array to store
the data. The implementation is based on
http://eternallyconfuzzled.com/tuts/datastructures/jsw_tut_avl.aspx

Lauri Kovanen (4/2011)
*/

#ifndef FIXEDTREE_H
#define FIXEDTREE_H

#include <stdlib.h>
#include <stdint.h>
#include <iterator>
#include <list>
#include <iostream>
#include <limits>

// Pre-declare Tree so it can be used by tree_iterator.
template<typename T> class FixedTree;

typedef uint32_t node_id;

template<typename T> class FixedNode
{
public:
  T value;
  node_id child[2];
  FixedNode();
};

template<typename T> FixedNode<T>::FixedNode():
value(0)
{
  child[0] = FixedTree<T>::null_node;
  child[1] = FixedTree<T>::null_node;
}

template <typename T> 
std::ostream& operator<<(std::ostream& output, const FixedNode<T>& node) 
{
  output << node.value << " (";
  if (node.child[0] == FixedTree<T>::null_node) output << "-,";
  else output << node.child[0] << ",";
  if (node.child[1] == FixedTree<T>::null_node) output << "-)";
  else output << node.child[1] << ")";
  return output;
}

template<typename T> inline T maximum(T a, T b)
{
  return (a > b ? a : b);
};

/* Class tree_iterator

   This is a very simple iteration interface to the array of the tree
   itself: contains just a pointer to the tree itself ('fixed_tree')
   and the current position ('pos').

   Note that while the iterator can be used at any time, the returned
   elements are ordered only when the underlying array structure in
   the tree is sorted.
 */
template<typename T>
class tree_iterator:public std::iterator<std::bidirectional_iterator_tag,T>
{
  friend class FixedTree<T>;
 private:
  const FixedTree<T> * fixed_tree;
  node_id pos;
  tree_iterator(const FixedTree<T> *at, node_id pos);

 public:
  tree_iterator();
  tree_iterator(const tree_iterator & ait);
  inline tree_iterator& operator++() { pos++; return *this;};
  inline void operator++(int) { pos++; };
  inline tree_iterator& operator--() { pos--; return *this;};
  inline void operator--(int) { pos--; };
  bool operator==(const tree_iterator & ait) const;
  bool operator!=(const tree_iterator & ait) const;
  inline const T& operator*() { return fixed_tree->nodes[pos].value; };
};

template <typename T>
tree_iterator<T>::tree_iterator():
fixed_tree(NULL),pos(0)
{}

template <typename T>
tree_iterator<T>::tree_iterator(const FixedTree<T> *at, node_id pos):
fixed_tree(at),pos(pos)
{}

template <typename T>
tree_iterator<T>::tree_iterator(const tree_iterator<T> & ait):
fixed_tree(ait.fixed_tree),pos(ait.pos)
{ }

template <typename T>
bool tree_iterator<T>::operator==(const tree_iterator<T> & ait) const
{
  return (fixed_tree == ait.fixed_tree && pos == ait.pos);
}

template <typename T>
bool tree_iterator<T>::operator!=(const tree_iterator<T> & ait) const
{
  
  return (fixed_tree != ait.fixed_tree || pos != ait.pos);
}

/* Class: FixedTree

   A binary tree with fixed size. FixedTree has been designed with the
   following usage scenario in mind:

     1) Create the tree with methods add() and Init().

     2) Use the tree as if it was a sorted array: log-time find() and
        constant time iteration of previous and next elements throught
	the iterator interface.

     3) Alter the tree: multiple calls of replace(), during which the
        iterator cannot be used but all tree-based methods run in
        log-time. Finally call restore_order() and return to 2.

   Because the size of the tree is fixed, we can use an array to store
   the elements which makes memory handling more efficient than using
   linked lists. At stage 2 the underlying array is ordered, which
   makes it possible to use the iterator to get the previous and next
   elements in constant time (the methods find_prev() and find_next()
   use the tree structure and run in log-time, but are usable also
   when the array is not sorted).
 */
template<typename T> class FixedTree
{
  friend class tree_iterator<T>;
private:
  unsigned int _size;
  node_id root;
 public:
  static const node_id null_node;
  typedef tree_iterator<T> iterator;

  FixedTree();
  FixedTree(const FixedTree& other);
  ~FixedTree();
  FixedTree& operator=(const FixedTree& other);

  /* These methods can be called at any time, even if the internal
     array is not sorted. All find methods are based on the tree
     structure, so they have log(n) running time; this is the expected
     running time if fast_replace has been called (without a call to
     restore_order), otherwise it is the worst case.
  */
  void Init(std::list<T> & values);
  void clear();
  bool empty() const { return (_size == 0); };
  int size() const { return _size; };
  void print() const;
  iterator find(T value) const;
  T find_min() const;
  T find_max() const;
  T find_prev(T value, T null_value) const;
  T find_next(T value, T null_value) const;

  /* An alternative to using Init(std::list<T>&) is to add values one
     by one with add(), and then call Init() after having added all
     elements. Even when using add() the values must still be added in
     correct order. The benefit is that no extra memory is needed to
     save the values twice; the downside is that this way of adding
     values is somewhat slower (though not that much).
   */
  void add(T value);
  void Init();

  /* Calling replace() will destroy the sorting of the internal
     array. Make sure to call restore_order() after (one or several)
     calls to replace(). The bidirectional iterator _cannot_ be used
     between these calls (or technically you can use it, but it
     returns absolute non-sense).
  */
  void replace(T old_value, T new_value);
  void restore_order();

  /* Bidirectional iterator interface. The iterator assumes the
     internal array is sorted. The only method that destroys sorting
     is replace(), so make sure to call restore_order() after calling
     replace(). The iterator will _not_ check if the internal array is
     sorted!  If it is not, the returned value can be anything.
  */
  inline iterator begin() const { return iterator(this,0); };
  inline iterator end() const { return iterator(this,_size); };
  inline iterator rbegin() const { return iterator(this,_size-1); };
  inline iterator rend() const { return iterator(this,null_node); };

  FixedNode<T> *nodes;

 private:

  node_id build_children(node_id first, node_id last);
  void sorted_array_copy(FixedNode<T> **new_nodes, node_id i);

  node_id _erase(T value);
  void _insert(T value, node_id pos);

  void debug_print() const { debug_print(nodes); };
  void debug_print(FixedNode<T> *_nodes) const;
  void _print(node_id v) const;

};

template<typename T> const node_id FixedTree<T>::null_node = std::numeric_limits<node_id>::max();

template<typename T>
  FixedTree<T>::FixedTree():_size(0),root(null_node),nodes(NULL)
{}

template<typename T> FixedTree<T>::FixedTree(const FixedTree<T>& other)
:_size(other._size),root(other.root),nodes(NULL)
{
  if (other.nodes != NULL)
    {
      nodes = new FixedNode<T>[_size];
      for (unsigned int i = 0; i < _size; ++i) nodes[i] = other.nodes[i];
    }
}

template<typename T>
FixedTree<T>& FixedTree<T>::operator=(const FixedTree<T>& other)
{
  if (this != &other)
    {
      root = other.root;
      _size = other._size;
      if (other.nodes != NULL)
	{
	  FixedNode<T> *new_nodes(NULL);
	  new_nodes = new FixedNode<T>[_size];
	  for (unsigned int i = 0; i < _size; ++i) new_nodes[i] = other.nodes[i];
	  delete[] nodes;
	  nodes = new_nodes;
	}
      else
	{
	  delete[] nodes;
	  nodes = NULL;
	}
    }
  return *this;
}

template<typename T>
void FixedTree<T>::Init()
{
  // Build the pointers to children so that the tree is balanced.
  root = build_children(0,_size-1);
}

template<typename T>
void FixedTree<T>::add(T value)
{
  // Create new array for nodes, copy contents to it, destroy the old
  // one and increase _size by one.
  FixedNode<T> *new_nodes = new FixedNode<T>[_size+1];
  for (unsigned int i = 0; i < _size; ++i) new_nodes[i] = nodes[i];
  new_nodes[_size].value = value;
  delete[] nodes;
  nodes = new_nodes;
  _size++;
}

template<typename T>
void FixedTree<T>::Init(std::list<T> & sorted_values)
{
  delete[] nodes;

  // Reserve space for the nodes.
  _size = sorted_values.size();
  if (_size > (unsigned int)null_node)
    {
      std::cerr << "Error in FixedTree::Init : Unable to create a tree with "<< _size <<" nodes.\n";
      exit(1);
    }
  nodes = new FixedNode<T>[_size];

  // Copy values to array.
  FixedNode<T> *curr_node = nodes;
  typename std::list<T>::const_iterator it;
  for (it = sorted_values.begin(); it != sorted_values.end(); ++it)
    {
      curr_node->value = *it;
      ++curr_node;
    }
  
  Init();
}

template<typename T>
void FixedTree<T>::restore_order()
{
  if (nodes == NULL) return;
  FixedNode<T> *new_nodes = new FixedNode<T>[_size];
  FixedNode<T> * const new_nodes_copy = new_nodes;
  sorted_array_copy(&new_nodes, root);
  delete[] nodes;
  nodes = new_nodes_copy;

  // Build the pointers to children so that the tree is balanced.
  root = build_children(0,_size-1);
  //std::cerr << "nodes " << nodes << std::endl;
  //debug_print(nodes);
}

template<typename T>
void FixedTree<T>::sorted_array_copy(FixedNode<T> **new_nodes, node_id i)
{
  if (nodes[i].child[0] != null_node) sorted_array_copy(new_nodes, nodes[i].child[0]);
  (**new_nodes).value = nodes[i].value;
  (*new_nodes)++;
  if (nodes[i].child[1] != null_node) sorted_array_copy(new_nodes, nodes[i].child[1]);
}

template<typename T>
node_id FixedTree<T>::build_children(node_id first, node_id last)
{
  if (first == last)
    {
      return first; // This node has no children, nothing needs to be done.
    }
  else if (last - first == 1)
    {
      // Only two nodes left in this branch, make the bigger one the
      // right child of the smaller one.
      nodes[first].child[1] = last;
      return first;
    }
  else
    {
      // At least three nodes left, continue recursively.
      node_id mid = first+(last-first)/2; // rounded down
      nodes[mid].child[0] = build_children(first, mid-1);
      nodes[mid].child[1] = build_children(mid+1,last);
      return mid;
    }
}

template<typename T>
FixedTree<T>::~FixedTree() 
{
  if (nodes != NULL) delete[] nodes;
}

template<typename T>
T FixedTree<T>::find_min() const
{
  node_id i = root;
  while (nodes[i].child[0] != null_node) i = nodes[i].child[0];
  return nodes[i].value;
}

template<typename T>
T FixedTree<T>::find_max() const
{
  node_id i = root;
  while (nodes[i].child[1] != null_node) i = nodes[i].child[1];
  return nodes[i].value;
}

template<typename T>
tree_iterator<T> FixedTree<T>::find(T value) const
{
  node_id i = root;
  while (i != null_node)
    {
      if (value == nodes[i].value) 
	{
	  return tree_iterator<T>(this, i);
	}
      i = nodes[i].child[nodes[i].value < value];
    }
  return tree_iterator<T>(this, _size);
}

template<typename T>
T FixedTree<T>::find_prev(T value, T null_value) const
{
  if(_size == 0) return null_value;
  T best = null_value;
  node_id i = root;
  while (i != null_node)
    {
      if(nodes[i].value < value)
        {
	  best = nodes[i].value;
	  i = nodes[i].child[1];
        }
      else i = nodes[i].child[0];
    }
  return best;
}

template<typename T>
T FixedTree<T>::find_next(T value, T null_value) const
{
    if(_size==0) return null_value;
    T best = null_value;
    node_id i = root;
    while (i != null_node)
    {
        if(nodes[i].value > value)
        {
            best = nodes[i].value;
            i = nodes[i].child[0];
        }
        else i = nodes[i].child[1];
    }
    return best;
}

template<typename T>
void FixedTree<T>::clear()
{
  if (nodes != NULL) delete[] nodes;
  nodes = NULL;
  _size = 0;
  root = null_node;
}

template<typename T>
void FixedTree<T>::replace(T old_value, T new_value)
{
  if (_size == 1) nodes[0].value = new_value;
  else
    {
      node_id pos = _erase(old_value);
      _insert(new_value, pos);
    }
}

template<typename T>
void FixedTree<T>::_insert(T value, node_id pos)
{
  node_id i = root;
  while(true)
    {
      int dir = (nodes[i].value < value);
      if (nodes[i].child[dir] == null_node)
	{
	  nodes[i].child[dir] = pos;
	  nodes[pos].value = value;
	  break;
	}
      i = nodes[i].child[dir];
    }
}

template<typename T>
node_id FixedTree<T>::_erase(T value)
{
  node_id p = null_node, i = root;
  while(true)
    {
      if (nodes[i].value == value) break;
      int dir = (nodes[i].value < value);
      p = i;
      i = nodes[i].child[dir];
    }

  if (nodes[i].child[0] != null_node && nodes[i].child[1] != null_node)
    {
      /* Find inorder successor of i. */
      p = i;
      node_id j = nodes[i].child[1];
      while (nodes[j].child[0] != null_node)
	{
	  p = j;
	  j = nodes[j].child[0];
	}
      
      /* Now i is the node to delete, j is its inorder successor
	 and p is the parent of j. We remove i by replacing the
	 contents of i by j. Note that j has no left child.
       */
      nodes[i].value = nodes[j].value;
      nodes[p].child[nodes[p].child[1] == j] = nodes[j].child[1];
      nodes[j].child[0] = null_node;
      nodes[j].child[1] = null_node;
      return j;
    }
  else
    {
      /* Node i has at most one child.
       */
      int dir = (nodes[i].child[0] == null_node);
      if (p == null_node) root = nodes[i].child[dir];
      else nodes[p].child[nodes[p].child[1] == i] = nodes[i].child[dir];
      nodes[i].child[0] = null_node;
      nodes[i].child[1] = null_node;
      return i;
    }
}


template<typename T>
void FixedTree<T>::print() const
{
  if (root != null_node) _print(root);
  std::cerr << std::endl;
}

template<typename T>
void FixedTree<T>::_print(node_id i) const
{
  if(nodes[i].child[0] != null_node) _print(nodes[i].child[0]);
  std::cerr << nodes[i].value << " ";
  if(nodes[i].child[1] != null_node) _print(nodes[i].child[1]);	
}

template<typename T>
void FixedTree<T>::debug_print(FixedNode<T>* _nodes) const
{
    for (node_id i = 0; i < _size; ++i)
        {
	  std::cerr << "nodes[" << i << "] = " << _nodes[i];
	  if (i == root) std::cerr << " <root>";
	  std::cerr << std::endl;
        }
}

#endif
