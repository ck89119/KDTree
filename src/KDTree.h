/**
 * File: KDTree.h
 * Author: Cao Kai
 * ------------------------
 * An interface representing a kd-tree in some number of dimensions. The tree
 * can be constructed from a set of data and then queried for membership and
 * nearest neighbors.
 */

#ifndef KDTREE_INCLUDED
#define KDTREE_INCLUDED

#include "Point.h"
#include "BoundedPQueue.h"
#include <stdexcept>
#include <cmath>
#include <map>
#include <iostream>

// "using namespace" in a header file is conventionally frowned upon, but I'm
// including it here so that you may use things like size_t without having to
// type std::size_t every time.
using namespace std;

template <size_t N, typename ElemType>
class KDTree {
private:
    struct Node {
      Node(const Point<N>& p): point(p) { left = right = NULL; }
      Point<N> point;
      Node *left, *right;
    };

public:
    // Constructor: KDTree();
    // Usage: KDTree<3, int> myTree;
    // ----------------------------------------------------
    // Constructs an empty KDTree.
    KDTree();
    
    // Destructor: ~KDTree()
    // Usage: (implicit)
    // ----------------------------------------------------
    // Cleans up all resources used by the KDTree.
    ~KDTree();
    
    // KDTree(const KDTree& rhs);
    // KDTree& operator=(const KDTree& rhs);
    // Usage: KDTree<3, int> one = two;
    // Usage: one = two;
    // -----------------------------------------------------
    // Deep-copies the contents of another KDTree into this one.
    KDTree(const KDTree& rhs);
    KDTree& operator=(const KDTree& rhs);
    
    // size_t dimension() const;
    // Usage: size_t dim = kd.dimension();
    // ----------------------------------------------------
    // Returns the dimension of the points stored in this KDTree.
    size_t dimension() const;
    
    // size_t size() const;
    // bool empty() const;
    // Usage: if (kd.empty())
    // ----------------------------------------------------
    // Returns the number of elements in the kd-tree and whether the tree is
    // empty.
    size_t size() const;
    bool empty() const;
    
    // bool contains(const Point<N>& pt) const;
    // Usage: if (kd.contains(pt))
    // ----------------------------------------------------
    // Returns whether the specified point is contained in the KDTree.
    bool contains(const Point<N>& pt) const;
    
    // void insert(const Point<N>& pt, const ElemType& value);
    // Usage: kd.insert(v, "This value is associated with v.");
    // ----------------------------------------------------
    // Inserts the point pt into the KDTree, associating it with the specified
    // value. If the element already existed in the tree, the new value will
    // overwrite the existing one.
    void insert(const Point<N>& pt, const ElemType& value);
    
    // ElemType& operator[](const Point<N>& pt);
    // Usage: kd[v] = "Some Value";
    // ----------------------------------------------------
    // Returns a reference to the value associated with point pt in the KDTree.
    // If the point does not exist, then it is added to the KDTree using the
    // default value of ElemType as its key.
    ElemType& operator[](const Point<N>& pt);
    
    // ElemType& at(const Point<N>& pt);
    // const ElemType& at(const Point<N>& pt) const;
    // Usage: cout << kd.at(v) << endl;
    // ----------------------------------------------------
    // Returns a reference to the key associated with the point pt. If the point
    // is not in the tree, this function throws an out_of_range exception.
    ElemType& at(const Point<N>& pt);
    const ElemType& at(const Point<N>& pt) const;
    
    // ElemType kNNValue(const Point<N>& key, size_t k) const
    // Usage: cout << kd.kNNValue(v, 3) << endl;
    // ----------------------------------------------------
    // Given a point v and an integer k, finds the k points in the KDTree
    // nearest to v and returns the most common value associated with those
    // points. In the event of a tie, one of the most frequent value will be
    // chosen.
    ElemType kNNValue(const Point<N>& key, size_t k) const;
    
    Node* CopyTree(const Node* node);

    void DestructTree(Node* node);

    void rec(const Point<N>& key, const Node* cur, BoundedPQueue<ElemType>& bpq, int d) const;

    void TestInsert(Node* cur);

    Node* GetRoot() const { return root; }

private:
    Node *root;
    map<Point<N>, ElemType> data;
};

/** KDTree class implementation details */

template <size_t N, typename ElemType>
KDTree<N, ElemType>::KDTree() {
  root = NULL;
  data.clear();
}

template <size_t N, typename ElemType>
KDTree<N, ElemType>::KDTree(const KDTree& rhs) {
  data = rhs.data;
  root = CopyTree(rhs.root);
}

template <size_t N, typename ElemType>
typename KDTree<N, ElemType>::Node* 
KDTree<N, ElemType>::CopyTree(const Node* node) {
  if (node == NULL) return NULL;

  Node *cur = new Node(node->point);
  cur->left = CopyTree(node->left);
  cur->right = CopyTree(node->right);
  return cur;
}

template <size_t N, typename ElemType>
KDTree<N, ElemType>::~KDTree() {
  data.clear();
  DestructTree(root);
}

template <size_t N, typename ElemType>
void KDTree<N, ElemType>::DestructTree(Node* node) {
  if (node == NULL) return ;
  DestructTree(node->left);
  DestructTree(node->right);
  delete node;
}

template <size_t N, typename ElemType>
ElemType& KDTree<N, ElemType>::at(const Point<N>& pt) {
  if (contains(pt)) return data[pt];
  else throw out_of_range("no such point");
}

template <size_t N, typename ElemType>
const ElemType& KDTree<N, ElemType>::at(const Point<N>& pt) const {
  auto it = data.find(pt);
  if (contains(pt)) return const_cast<ElemType&>(it->second);
  else throw out_of_range("no such point");
}

template <size_t N, typename ElemType>
bool KDTree<N, ElemType>::contains(const Point<N>& pt) const {
  return data.find(pt) != data.end();
}

template <size_t N, typename ElemType>
size_t KDTree<N, ElemType>::dimension() const {
  return N;
}

template <size_t N, typename ElemType>
bool KDTree<N, ElemType>::empty() const {
  return data.empty();
}

template <size_t N, typename ElemType>
void KDTree<N, ElemType>::insert(const Point<N>& pt, const ElemType& value) {
  if (contains(pt)) {
    data[pt] = value;
    return;
  }

  data[pt] = value;
  if (root == NULL) {
    root = new Node(pt);
    return;
  }

  Node* cur = root;
  int d = 0;
  while (true) {
    if (pt[d] < cur->point[d]) {
      if (cur->left == NULL) {
        cur->left = new Node(pt);
        return;
      }
      cur = cur->left;
    }
    else {
      if (cur->right == NULL) {
        cur->right = new Node(pt);
        return;
      }
      cur = cur->right;
    }
    d = (d + 1) % N;
  }
}

template <size_t N, typename ElemType>
void KDTree<N, ElemType>::TestInsert(Node* cur) {
  cout << "###test_insert###" << endl;
  if (cur == NULL) return;
  for (int i = 0; i < N; i++)
    cout << cur->point[i];
  cout << ": " << endl;

  if (cur->left != NULL) {
    cout << "left: ";
    for (int i = 0; i < N; i++)
      cout << cur->left->point[i];
    cout << endl;
    test_insert(cur->left);
  }

  if (cur->right != NULL) {
    cout << "right: ";
    for (int i = 0; i < N; i++)
      cout << cur->right->point[i];
    cout << endl;
    test_insert(cur->right);
  }
}

template <size_t N, typename ElemType>
ElemType KDTree<N, ElemType>::kNNValue(const Point<N>& key, size_t k) const {
  BoundedPQueue<ElemType> bpq(k);
  rec(key, root, bpq, 0);
  map<ElemType, int> elem_count;
  while (!bpq.empty()) {
    ElemType tmp = bpq.dequeueMin();
    elem_count[tmp] += 1;
  }

  ElemType ret;
  int max_count = -1;
  for (auto it: elem_count) {
    if (it.second > max_count) {
      max_count = it.second;
      ret = it.first;
    }
  }
  return ret;
}

template <size_t N, typename ElemType>
void KDTree<N, ElemType>::rec(const Point<N>& key, const Node* cur, BoundedPQueue<ElemType>& bpq, int d) const {
  if (cur == NULL) return;

  auto it = data.find(cur->point);
  bpq.enqueue(it->second, Distance(key, cur->point));

  Node* other_sub = NULL;
  if (key[d] < cur->point[d]) {
    rec(key, cur->left, bpq, (d + 1) % N);
    other_sub = cur->right; 
  } else {
    rec(key, cur->right, bpq, (d + 1) % N);
    other_sub = cur->left; 
  }
  if (bpq.size() < bpq.maxSize() || fabs(key[d] - cur->point[d]) + 1e-9 < bpq.worst()) {
    rec(key, other_sub, bpq, (d + 1) % N);
  }
}

template <size_t N, typename ElemType>
KDTree<N, ElemType>& KDTree<N, ElemType>::operator=(const KDTree& rhs) {
  if (this == &rhs) 
    return *this;

  data.clear();
  DestructTree(root);
  data = rhs.data;
  root = CopyTree(rhs.root);
  return *this;
}

template <size_t N, typename ElemType>
ElemType& KDTree<N, ElemType>::operator[](const Point<N>& pt) {
  if (!contains(pt)) {
    insert(pt, 0); 
  }
  return data[pt];          
}

template <size_t N, typename ElemType>
size_t KDTree<N, ElemType>::size() const {
  return data.size();
}


#endif // KDTREE_INCLUDED
