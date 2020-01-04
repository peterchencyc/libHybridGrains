// Utilities.h
//
// Breannan Smith
// Last updated: 09/03/2015

#ifndef UTILITIES_H
#define UTILITIES_H

#include <bitset>
#include <cassert>
#include <istream>
#include <memory>
#include <vector>

namespace Utilities {
template <typename Test, template <typename...> class Ref>
struct is_specialization : public std::false_type {};

template <template <typename...> class Ref, typename... Args>
struct is_specialization<Ref<Args...>, Ref> : public std::true_type {};

// Supresses unused variable warnings
template <class T> void ignoreUnusedVariable(const T &) {}

// Performs a deep copy of a vector of pointers
// Precondition: copy.empty() is true
template <class T>
void cloneVector(const std::vector<T *> &original, std::vector<T *> &copy) {
  assert(copy.empty());
  copy.resize(original.size());
  for (typename std::vector<T *>::size_type i = 0; i < original.size(); ++i) {
    copy[i] = original[i]->clone();
  }
}

// Performs a deep copy of a vector of unique_ptrs
template <class T>
void cloneVector(const std::vector<std::unique_ptr<T>> &original,
                 std::vector<std::unique_ptr<T>> &copy) {
  copy.clear();
  copy.resize(original.size());
  for (typename std::vector<std::unique_ptr<T>>::size_type i = 0;
       i < original.size(); ++i) {
    copy[i] = std::unique_ptr<T>(original[i]->clone());
  }
}

// Returns a 'deep copy' of a vector of unique_ptrs
template <class T>
std::vector<std::unique_ptr<T>>
cloneVector(const std::vector<std::unique_ptr<T>> &original) {
  std::vector<std::unique_ptr<T>> copy;
  cloneVector(original, copy);
  return copy;
}

// Computes the index of an element of a vector
// Precondition: element is contained in vec
template <class T> unsigned index(const std::vector<T> &vec, const T &element) {
  return unsigned(&element - vec.data());
}

// TODO: Don't use a string to serialize
template <std::size_t N>
void serializeBitset(const std::bitset<N> &bs, std::ostream &output_stream) {
  const std::string bitset_string{bs.to_string()};
  assert(bitset_string.size() == N);
  output_stream.write(const_cast<char *>(bitset_string.c_str()),
                      N * sizeof(char));
}

// TODO: Don't use a string to deserialize
template <std::size_t N>
std::bitset<N> deserializeBitset(std::istream &input_stream) {
  std::vector<char> bitset_string(N);
  input_stream.read(bitset_string.data(), N * sizeof(char));
  return std::bitset<N>{
      std::string{bitset_string.cbegin(), bitset_string.cend()}};
}

template <typename T>
void serializeBuiltInType(const T &var, std::ostream &output_stream) {
  static_assert(std::is_integral<T>::value ||
                    std::is_floating_point<T>::value || std::is_enum<T>::value,
                "Error, serializeBuiltInType only works with integral and "
                "floating point types.");
  assert(output_stream.good());
  T local_var = var;
  output_stream.write(reinterpret_cast<char *>(&local_var), sizeof(T));
}

template <typename T>
void deserializeBuiltInType(T &var, std::istream &input_stream) {
  static_assert(std::is_integral<T>::value ||
                    std::is_floating_point<T>::value || std::is_enum<T>::value,
                "Error, deserializeBuiltInType only works with integral and "
                "floating point types.");
  assert(input_stream.good());
  input_stream.read(reinterpret_cast<char *>(&var), sizeof(T));
}

template <typename T> T deserialize(std::istream &input_stream) {
  T type_to_read;
  deserializeBuiltInType(type_to_read, input_stream);
  return type_to_read;
}

template <class T>
void serializeVectorBuiltInType(const std::vector<T> &vector,
                                std::ostream &output_stream) {
  assert(output_stream.good());
  // Write out the length of the vector
  Utilities::serializeBuiltInType(vector.size(), output_stream);
  assert(output_stream.good());
  // Output each element of the vector
  for (typename std::vector<T>::size_type idx = 0; idx < vector.size(); ++idx) {
    Utilities::serializeBuiltInType(vector[idx], output_stream);
    assert(output_stream.good());
  }
}

template <>
void serializeVectorBuiltInType(const std::vector<bool> &vector,
                                std::ostream &output_stream);

template <class T>
void deserializeVectorBuiltInType(std::vector<T> &vector,
                                  std::istream &input_stream) {
  assert(input_stream.good());
  // Read the length of the vector
  const typename std::vector<T>::size_type length{
      deserialize<typename std::vector<T>::size_type>(input_stream)};
  assert(input_stream.good());
  // Allocate space for the new data
  vector.resize(length);
  // Read each element of the vector
  for (typename std::vector<T>::size_type idx = 0; idx < vector.size(); ++idx) {
    T loaded_val;
    deserializeBuiltInType(loaded_val, input_stream);
    vector[idx] = loaded_val;
    assert(input_stream.good());
  }
}

template <class T>
void serializeVectorCustomType(const std::vector<T> &vector,
                               std::ostream &output_stream) {
  assert(output_stream.good());
  // Write out the length of the vector
  Utilities::serializeBuiltInType(vector.size(), output_stream);
  assert(output_stream.good());
  // Output each element of the vector
  for (typename std::vector<T>::size_type idx = 0; idx < vector.size(); ++idx) {
    vector[idx].serialize(output_stream);
    assert(output_stream.good());
  }
}

template <class T, class U>
void serializeVectorCustomType(const std::vector<T, U> &vector,
                               std::ostream &output_stream) {
  assert(output_stream.good());
  // Write out the length of the vector
  Utilities::serializeBuiltInType(vector.size(), output_stream);
  assert(output_stream.good());
  // Output each element of the vector
  for (typename std::vector<T, U>::size_type idx = 0; idx < vector.size();
       ++idx) {
    vector[idx].serialize(output_stream);
    assert(output_stream.good());
  }
}

template <class T>
void serializeVectorCustomType(const std::vector<T> &vector,
                               void (*serializationRoutine)(const T &,
                                                            std::ostream &),
                               std::ostream &output_stream) {
  assert(output_stream.good());
  // Write out the length of the vector
  Utilities::serializeBuiltInType(vector.size(), output_stream);
  assert(output_stream.good());
  // Output each element of the vector
  for (typename std::vector<T>::size_type idx = 0; idx < vector.size(); ++idx) {
    serializationRoutine(vector[idx], output_stream);
    assert(output_stream.good());
  }
}

template <class T>
void deserializeVectorCustomType(std::vector<T> &vector,
                                 std::istream &input_stream) {
  assert(input_stream.good());
  // Read the length of the vector
  const typename std::vector<T>::size_type length{
      deserialize<typename std::vector<T>::size_type>(input_stream)};
  assert(input_stream.good());
  // Allocate space for the new data
  vector.reserve(length);
  // Read each element of the vector
  for (typename std::vector<T>::size_type idx = 0; idx < length; ++idx) {
    vector.emplace_back(input_stream);
    assert(input_stream.good());
  }
  assert(vector.size() == length);
}

template <class T, class U>
void deserializeVectorCustomType(std::vector<T, U> &vector,
                                 std::istream &input_stream) {
  assert(input_stream.good());
  // Read the length of the vector
  const typename std::vector<T, U>::size_type length{
      deserialize<typename std::vector<T, U>::size_type>(input_stream)};
  assert(input_stream.good());
  // Allocate space for the new data
  vector.reserve(length);
  // Read each element of the vector
  for (typename std::vector<T, U>::size_type idx = 0; idx < length; ++idx) {
    vector.emplace_back(input_stream);
    assert(input_stream.good());
  }
  assert(vector.size() == length);
}

template <class T>
void deserializeVectorCustomType(std::vector<T> &vector,
                                 T (*deserializationRoutine)(std::istream &),
                                 std::istream &input_stream) {
  assert(input_stream.good());
  // Read the length of the vector
  const typename std::vector<T>::size_type length{
      deserialize<typename std::vector<T>::size_type>(input_stream)};
  assert(input_stream.good());
  // Allocate space for the new data
  vector.resize(length);
  // Read each element of the vector
  for (typename std::vector<T>::size_type idx = 0; idx < vector.size(); ++idx) {
    vector[idx] = deserializationRoutine(input_stream);
    assert(input_stream.good());
  }
}

template <class T>
std::vector<T>
deserializeVectorCustomType(T (*deserializationRoutine)(std::istream &),
                            std::istream &input_stream) {
  assert(input_stream.good());
  // Read the length of the vector
  const typename std::vector<T>::size_type length{
      deserialize<typename std::vector<T>::size_type>(input_stream)};
  assert(input_stream.good());
  // Allocate space for the new data
  std::vector<T> vector(length);
  // Read each element of the vector
  for (typename std::vector<T>::size_type idx = 0; idx < vector.size(); ++idx) {
    vector[idx] = deserializationRoutine(input_stream);
    assert(input_stream.good());
  }
  return vector;
}

template <class T>
void serializeVectorCustomTypePointers(const std::vector<T> &vector,
                                       std::ostream &output_stream) {
  assert(output_stream.good());
  // Write out the length of the vector
  Utilities::serializeBuiltInType(vector.size(), output_stream);
  assert(output_stream.good());
  // Output each element of the vector
  for (typename std::vector<T>::size_type idx = 0; idx < vector.size(); ++idx) {
    assert(vector[idx] != nullptr);
    vector[idx]->serialize(output_stream);
    assert(output_stream.good());
  }
}

} // namespace Utilities

#endif
