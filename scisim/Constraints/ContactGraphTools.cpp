// ContactGraphTools.cpp
//
// Breannan Smith
// Last updated: 10/27/2015

#include "ContactGraphTools.h"

#include "Constraint.h"
#include "scisim/Math/MathDefines.h"

// TODO: Faster to reserve space in each entry c[i] of a[i].size + b[i].size ?
template <typename T>
static void mergeWithoutDuplicates(const std::vector<T> &a,
                                   const std::vector<T> &b, std::vector<T> &c) {
  c.clear();
  typename std::vector<T>::const_iterator aitr{a.begin()};
  typename std::vector<T>::const_iterator bitr{b.begin()};
  while (true) {
    // If we have exhausted all elements of a
    if (aitr == a.end()) {
      // Insert the remainder of b
      c.insert(c.end(), bitr, b.end());
      break;
    }
    // If we have exhausted all elements of b
    if (bitr == b.end()) {
      // Insert the remainder of a
      c.insert(c.end(), aitr, a.end());
      break;
    }
    // If the next element of a is less than the next element of b
    if (*aitr < *bitr) {
      // Insert the next element of a
      c.emplace_back(*aitr++);
    }
    // If the next element of b is greater than the next element of b
    else if (*aitr > *bitr) {
      // Insert the next element of b
      c.emplace_back(*bitr++);
    }
    // Otherwise, the next element of a and the next element of b must be equal
    else {
      // Only insert the equal element once, but advance the pointers to the
      // next entry for each array
      assert(*aitr == *bitr);
      c.emplace_back(*aitr++);
      ++bitr;
    }
  }
}

static void buildContactGraph(
    const unsigned nbodies,
    const std::vector<std::unique_ptr<Constraint>> &constraints,
    std::vector<
        std::vector<std::vector<std::unique_ptr<Constraint>>::size_type>> &c) {
  typedef std::vector<std::unique_ptr<Constraint>>::size_type stype;

  const stype ncons{constraints.size()};

  // Build a map that, for each body, contains the indices of constraints acting
  // on that body
  // TODO: Use a hash map here
  std::vector<std::vector<stype>> b(nbodies);
  for (stype cidx = 0; cidx < ncons; ++cidx) {
    std::pair<int, int> bodies;
    constraints[cidx]->getSimulatedBodyIndices(bodies);
    assert(bodies.first < int(nbodies));
    assert(bodies.first >= 0);
    assert(bodies.second < int(nbodies));
    assert(bodies.second >= -1);
    b[bodies.first].emplace_back(cidx);
    if (bodies.second >= 0) {
      b[bodies.second].emplace_back(cidx);
    }
  }

  // For each constraint, build a map of neighboring constraints
  c.resize(ncons);
  for (stype cidx = 0; cidx < ncons; ++cidx) {
    std::pair<int, int> bodies;
    constraints[cidx]->getSimulatedBodyIndices(bodies);
    assert(bodies.first < int(nbodies));
    assert(bodies.first >= 0);
    assert(bodies.second < int(nbodies));
    assert(bodies.second >= -1);
    if (bodies.second != -1) {
      assert(!b[bodies.first].empty());
      assert(!b[bodies.second].empty());
      mergeWithoutDuplicates(b[bodies.first], b[bodies.second], c[cidx]);
      assert(c[cidx].size() <=
             b[bodies.first].size() + b[bodies.second].size());
    } else {
      c[cidx] = b[bodies.first];
    }
  }

#ifndef NDEBUG
  for (stype cidx = 0; cidx < ncons; ++cidx) {
    // The neighbor list must contain at least the constraint itself
    assert(!c[cidx].empty());
    for (std::vector<stype>::size_type idx = 0; idx < c[cidx].size() - 1;
         ++idx) {
      // And the neighbors should be in sorted order
      assert(c[cidx][idx] < c[cidx][idx + 1]);
    }
  }
#endif
}

static VectorXu
degrees(const std::vector<
        std::vector<std::vector<std::unique_ptr<Constraint>>::size_type>> &c) {
  const std::vector<std::unique_ptr<Constraint>>::size_type nc{c.size()};
  VectorXu node_degrees{static_cast<VectorXu::Index>(nc)};
  for (std::vector<std::vector<std::vector<
           std::unique_ptr<Constraint>>::size_type>>::size_type idx = 0;
       idx < nc; ++idx) {
    // The degree of each node is simply the count of neighbors
    node_degrees(idx) = unsigned(c[idx].size());
  }
  return node_degrees;
}

// Reorders the given vector according to the permutation, but sorts the
// permutation in the process
static void
destructiveReorder(VectorXi &permutation,
                   std::vector<std::unique_ptr<Constraint>> &constraints) {
  assert(permutation.size() ==
         static_cast<VectorXi::Index>(constraints.size()));
  for (int idx = 0; idx < permutation.size(); ++idx) {
    while (idx != permutation[idx]) {
      const int fnl_idx{permutation[idx]};
      assert(fnl_idx < permutation.size());
      using std::swap;
      swap(permutation[idx], permutation[fnl_idx]);
      swap(constraints[idx], constraints[fnl_idx]);
    }
  }
#ifndef NDEBUG
  for (int idx = 0; idx < permutation.size(); ++idx) {
    assert(permutation(idx) == idx);
  }
#endif
}

// Argsorts a collection of indices based on the the values in 'vals'
static VectorXi generateSortedIndices(const VectorXu &vals) {
  VectorXi indices{vals.size()};
  for (int idx = 0; idx < indices.size(); ++idx) {
    indices(idx) = idx;
  }
  std::sort(indices.data(), indices.data() + indices.size(),
            [&vals](const size_t i0, const size_t i1) {
              return vals(i0) < vals(i1);
            });
  return indices;
}

static VectorXi generateInverseIndices(const VectorXi &indices) {
  VectorXi inverse{indices.size()};
  for (int idx = 0; idx < indices.size(); ++idx) {
    inverse[indices[idx]] = idx;
  }
  return inverse;
}

static void
cuthillMckee(const unsigned long num_bodies,
             const std::vector<std::unique_ptr<Constraint>> &constraints,
             VectorXi &permutation) {
  // Build the contact graph from the constraints
  std::vector<std::vector<std::vector<std::unique_ptr<Constraint>>::size_type>>
      contact_graph;
  buildContactGraph(unsigned(num_bodies), constraints, contact_graph);

  // Compute the degree of each node in the contact graph
  const VectorXu node_degrees{degrees(contact_graph)};

  assert(node_degrees.size() == permutation.size());

  // Indices corresponding to the DoFs in node_degrees, sorted in ascending
  // order. Also used to label visited nodes.
  VectorXi indices{generateSortedIndices(node_degrees)};
#ifndef NDEBUG
  for (int idx = 0; idx < indices.size() - 1; ++idx) {
    assert(node_degrees(indices(idx)) <= node_degrees(indices(idx + 1)));
  }
#endif

  // Maps an original index to its position in indices
  const VectorXi rev_indices{generateInverseIndices(indices)};
#ifndef NDEBUG
  for (int index = 0; index < indices.size(); ++index) {
    assert(rev_indices[indices[index]] == index);
    assert(indices[rev_indices[index]] == index);
  }
#endif

  // Temporary storage for the insertion sort
  VectorXu tmp_degrees{node_degrees.maxCoeff()};

  unsigned num_inserted_in_permute{0};

  // Step through the dofs one by one; this loop is required for non-connected
  // graphs
  for (unsigned outer_permute_idx = 0; outer_permute_idx < permutation.size();
       ++outer_permute_idx) {
    // If we haven't visited the current dof, start a breadth first search from
    // the current dof
    if (indices(outer_permute_idx) != -1) {
      // Insert the current node into its final position in the permutation
      assert(outer_permute_idx < indices.size());
      permutation(num_inserted_in_permute++) = indices(outer_permute_idx);
      // Mark the node as visitied
      indices(outer_permute_idx) = -1;

      // For the first iteration of this island, start the breadth first search
      // from a single vertex
      assert(num_inserted_in_permute >= 1);
      unsigned last_insert_p_start{num_inserted_in_permute - 1};
      unsigned last_insert_p_end{num_inserted_in_permute};

      // If we inserted a node into the permutation during the last iteration
      while (last_insert_p_start < last_insert_p_end) {
        // For each node that we inserted during the last iteration of breadth
        // first search
        for (unsigned cur_node_p_idx = last_insert_p_start;
             cur_node_p_idx < last_insert_p_end; ++cur_node_p_idx) {
          const unsigned num_inserted_prev{num_inserted_in_permute};

          // For each child node of the current node
          for (const size_t child_node_orig_idx :
               contact_graph[permutation(cur_node_p_idx)]) {
            // If the child node has not been visited
            if (indices(rev_indices(child_node_orig_idx)) != -1) {
              // Mark the child as visited
              assert(indices(rev_indices(child_node_orig_idx)) ==
                     int(child_node_orig_idx));
              indices(rev_indices(child_node_orig_idx)) = -1;
              // Assign the neighbor a position in the permuted state
              permutation(num_inserted_in_permute++) = int(child_node_orig_idx);
            }
          }

          // Copy the degrees of the nodes we inserted into contiguous storage
          for (unsigned inserted_idx = num_inserted_prev;
               inserted_idx < num_inserted_in_permute; ++inserted_idx) {
            assert(inserted_idx < permutation.size());
            assert(permutation(inserted_idx) < node_degrees.size());
            tmp_degrees(inserted_idx - num_inserted_prev) =
                node_degrees(permutation(inserted_idx));
          }

          // Sort the inserted nodes in the permutation by degree (by insertion
          // sort)
          for (unsigned node_deg_idx = 1;
               node_deg_idx < num_inserted_in_permute - num_inserted_prev;
               ++node_deg_idx) {
            const unsigned mvd_degree{tmp_degrees(node_deg_idx)};
            const int mvd_p_val{permutation(num_inserted_prev + node_deg_idx)};
            unsigned cur_deg_idx{node_deg_idx};
            while ((cur_deg_idx > 0) &&
                   (mvd_degree < tmp_degrees[cur_deg_idx - 1])) {
              tmp_degrees(cur_deg_idx) = tmp_degrees(cur_deg_idx - 1);
              permutation(num_inserted_prev + cur_deg_idx) =
                  permutation(num_inserted_prev + cur_deg_idx - 1);
              cur_deg_idx -= 1;
            }
            tmp_degrees(cur_deg_idx) = mvd_degree;
            permutation(num_inserted_prev + cur_deg_idx) = mvd_p_val;
          }

// Verify that the nodes inserted at this level are sorted by degree
#ifndef NDEBUG
          for (unsigned inserted_idx = num_inserted_prev;
               inserted_idx < num_inserted_in_permute - 1; ++inserted_idx) {
            assert(inserted_idx < permutation.size());
            assert(permutation(inserted_idx) < node_degrees.size());
            assert(node_degrees(permutation(inserted_idx)) <=
                   node_degrees(permutation(inserted_idx + 1)));
          }
#endif
        }

        // Update the range of nodes that we just inserted into the permutation
        last_insert_p_start = last_insert_p_end;
        last_insert_p_end = num_inserted_in_permute;
      }
    }
    // Possible early out, if this island exhausted all nodes
    if (num_inserted_in_permute == permutation.size()) {
      break;
    }
  }
// Verify that we have a valid permutation
#ifndef NDEBUG
  VectorXi sorted_permutation{permutation};
  std::sort(sorted_permutation.data(),
            sorted_permutation.data() + sorted_permutation.size());
  for (int i = 0; i < sorted_permutation.size(); ++i) {
    assert(sorted_permutation(i) == i);
  }
#endif
}

void ContactGraphTools::reduceBandwidth(
    const unsigned long num_bodies,
    std::vector<std::unique_ptr<Constraint>> &constraints) {
  VectorXi permutation{static_cast<VectorXi::Index>(constraints.size())};
  cuthillMckee(num_bodies, constraints, permutation);
  // TODO: Generate the inverted version directly in the above code
  VectorXi inverted_permutation{generateInverseIndices(permutation)};
  destructiveReorder(inverted_permutation, constraints);
}
