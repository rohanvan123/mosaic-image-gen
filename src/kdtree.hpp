/**
 * @file kdtree.cpp
 * Implementation of KDTree class.
 */

#include <utility>
#include <algorithm>

using namespace std;

template <int Dim>
bool KDTree<Dim>::smallerDimVal(const Point<Dim>& first,
                                const Point<Dim>& second, int curDim) const
{
    /**
     * @todo Implement this function!
     */

    // return true if the first point is smaller than the second given the dim
    if (first[curDim] < second[curDim]) {
      return true;
    } else if (first[curDim] > second[curDim]) {
      return false;
    } else {
      return first < second;
    }
}

template <int Dim>
bool KDTree<Dim>::shouldReplace(const Point<Dim>& target,
                                const Point<Dim>& currentBest,
                                const Point<Dim>& potential) const
{
    /**
     * @todo Implement this function!
     */

     // need to calculate the distance between currentBest and target and potential and Target, and compare
    int pot_dist = 0;
    int curr_dist = 0;
    for (int d = 0; d < Dim; d++) {
      pot_dist += (potential[d] - target[d]) * (potential[d] - target[d]);
      curr_dist += (currentBest[d] - target[d]) * (currentBest[d] - target[d]);
    }

    if (pot_dist < curr_dist) {
      return true;
    } else if (pot_dist > curr_dist) {
      return false;
    }

    return potential < currentBest;
}

template <int Dim>
KDTree<Dim>::KDTree(const vector<Point<Dim>>& newPoints)
{
    /**
     * @todo Implement this function!
     */
    if (newPoints.empty()) {
      root = NULL;
      size = 0;
    } else {
      vector<Point<Dim>> copy(newPoints.begin(), newPoints.end());
      root = buildNode(copy, 0, copy.size() - 1, 0);
      size = newPoints.size();
    }
}

template <int Dim>
KDTree<Dim>::KDTree(const KDTree<Dim>& other) {
  /**
   * @todo Implement this function!
   */
}

template <int Dim>
const KDTree<Dim>& KDTree<Dim>::operator=(const KDTree<Dim>& rhs) {
  /**
   * @todo Implement this function!
   */

  return *this;
}

template <int Dim>
KDTree<Dim>::~KDTree() {
  /**
   * @todo Implement this function!
   */
}

template <int Dim>
Point<Dim> KDTree<Dim>::findNearestNeighbor(const Point<Dim>& query) const
{
    /**
     * @todo Implement this function!
     */
    if (root == NULL) {
      return Point<Dim>();
    }
    return nearestNeighbor(query, root, 0);
}

// helper functions
template <int Dim>
typename KDTree<Dim>::KDTreeNode* KDTree<Dim>::buildNode(std::vector<Point<Dim>>& points, int l, int r, int d) {
  if (points.empty() || l > r) {
    return NULL;
  }
  int med = (l + r) / 2;
  Point<Dim> median = quickselect(points, l, r, med, d);
  KDTreeNode* node = new KDTreeNode(median);
  node->left = buildNode(points, l, med - 1, (d + 1) % Dim);
  node->right = buildNode(points, med + 1, r, (d + 1) % Dim);
  return node;
}

template <int Dim>
Point<Dim> KDTree<Dim>::quickselect(std::vector<Point<Dim>>& points, int l, int r, int k, int d) {
  if (l == r) {
      return points[l];
  }
  
  int pivotIndex = (l + r) / 2;
  pivotIndex = partition(points, l, r, pivotIndex, d);
  
  if (k == pivotIndex) {
      return points[k];
  } else if (k < pivotIndex) {
      return quickselect(points, l, pivotIndex - 1, k, d);
  } else {
      return quickselect(points, pivotIndex + 1, r, k, d);
  }
}

template <int Dim>
int KDTree<Dim>::partition(std::vector<Point<Dim>>& points, int l, int r, int pivotIndex, int d) {
  Point<Dim> pivotPoint = points[pivotIndex];
  swap(points, pivotIndex, r);
  int storeIndex = l;
  for (int i = l; i < r; i++) {
      if (smallerDimVal(points[i], pivotPoint, d)) {
          swap(points, storeIndex, i);
          storeIndex++;
      }
  }

  swap(points, storeIndex, r);
  return storeIndex;
}

template <int Dim> 
void KDTree<Dim>::swap(std::vector<Point<Dim>>& points, int b, int s) {
  Point<Dim> tmp = points[b];
  points[b] = points[s];
  points[s] = tmp;
}

template <int Dim>
Point<Dim> KDTree<Dim>::nearestNeighbor(const Point<Dim>& query, KDTreeNode* curr, int d) const {
  // Look in the left or right subtree depending on whether target is smaller or larger than our current root
  bool left = false;
  bool right = false;
  Point<Dim> childResult;

  if (smallerDimVal(query, curr->point, d)) {
      if (curr->left == NULL)
        return curr->point;
      childResult = nearestNeighbor(query, curr->left, (d + 1) % Dim);
      left = true;
  } else {
      if (curr->right == NULL)
        return curr->point;
      childResult = nearestNeighbor(query, curr->right, (d + 1) % Dim);
      right = true;
  }

  if (shouldReplace(query, childResult, curr->point)) {
    childResult = curr->point;
  } 
  // after replacing, then we just check if the split distance (at that dimension) is smaller then the current radius bc it being in the same subspace is dependent on split dim 
  int radius = 0;
  for (int i = 0; i < Dim; i++) {
    radius += (query[i] - childResult[i]) * (query[i] - childResult[i]);
  }
  int split = (query[d] - curr->point[d]) * (query[d] - curr->point[d]);

  if (split <= radius) { 
    if (left && curr->right) {
      Point<Dim> potential = nearestNeighbor(query, curr->right, (d + 1) % Dim);
      if (shouldReplace(query, childResult, potential)) {
        childResult = potential;
      }
    } else if (right && curr->left) {
      Point<Dim> potential = nearestNeighbor(query, curr->left, (d + 1) % Dim);
      if (shouldReplace(query, childResult, potential)) {
        childResult = potential;
      }
    }
  }

  return childResult;
}

template <int Dim>
bool KDTree<Dim>::sameDistance(const Point<Dim>& target,
                                const Point<Dim>& currentBest,
                                const Point<Dim>& potential) const {
                                  int pot_dist = 0;
    int curr_dist = 0;
    for (int d = 0; d < Dim; d++) {
      pot_dist += (potential[d] - target[d]) * (potential[d] - target[d]);
      curr_dist += (currentBest[d] - target[d]) * (currentBest[d] - target[d]);
    }
    return pot_dist == curr_dist;
}