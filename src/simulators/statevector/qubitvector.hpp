/**
 * This code is part of Qiskit.
 *
 * (C) Copyright IBM 2018, 2019.
 *
 * This code is licensed under the Apache License, Version 2.0. You may
 * obtain a copy of this license in the LICENSE.txt file in the root directory
 * of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
 *
 * Any modifications or derivative works of this code must retain this
 * copyright notice, and modified files need to carry a notice indicating
 * that they have been altered from the originals.
 */



#ifndef _qv_qubit_vector_hpp_
#define _qv_qubit_vector_hpp_

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "framework/json.hpp"

#ifdef SIMD_PPC64LE
#include<altivec.h>
#endif


namespace QV {

// Type aliases
using uint_t = uint64_t;
using int_t = int64_t;
using reg_t = std::vector<uint_t>;
using indexes_t = std::unique_ptr<uint_t[]>;
template <size_t N> using areg_t = std::array<uint_t, N>;
template <typename T> using cvector_t = std::vector<std::complex<T>>;

//============================================================================
// BIT MASKS and indexing
//============================================================================

/*
# Auto generate these values with following python snippet
import json

def cpp_init_list(int_lst):
    ret = json.dumps([str(i) + 'ULL' for i in int_lst])
    return ret.replace('"', '').replace('[','{{').replace(']','}};')

print('const std::array<uint_t, 64> BITS = ' + cpp_init_list([(1 << i) for i in range(64)]) + '\n')
print('const std::array<uint_t, 64> MASKS = ' + cpp_init_list([(1 << i) - 1 for i in range(64)]))
*/

const std::array<uint_t, 64> BITS {{
  1ULL, 2ULL, 4ULL, 8ULL,
  16ULL, 32ULL, 64ULL, 128ULL,
  256ULL, 512ULL, 1024ULL, 2048ULL,
  4096ULL, 8192ULL, 16384ULL, 32768ULL,
  65536ULL, 131072ULL, 262144ULL, 524288ULL,
  1048576ULL, 2097152ULL, 4194304ULL, 8388608ULL,
  16777216ULL, 33554432ULL, 67108864ULL, 134217728ULL,
  268435456ULL, 536870912ULL, 1073741824ULL, 2147483648ULL,
  4294967296ULL, 8589934592ULL, 17179869184ULL, 34359738368ULL,
  68719476736ULL, 137438953472ULL, 274877906944ULL, 549755813888ULL,
  1099511627776ULL, 2199023255552ULL, 4398046511104ULL, 8796093022208ULL,
  17592186044416ULL, 35184372088832ULL, 70368744177664ULL, 140737488355328ULL, 
  281474976710656ULL, 562949953421312ULL, 1125899906842624ULL, 2251799813685248ULL,
  4503599627370496ULL, 9007199254740992ULL, 18014398509481984ULL, 36028797018963968ULL,
  72057594037927936ULL, 144115188075855872ULL, 288230376151711744ULL, 576460752303423488ULL,
  1152921504606846976ULL, 2305843009213693952ULL, 4611686018427387904ULL, 9223372036854775808ULL
}};


const std::array<uint_t, 64> MASKS {{
  0ULL, 1ULL, 3ULL, 7ULL,
  15ULL, 31ULL, 63ULL, 127ULL,
  255ULL, 511ULL, 1023ULL, 2047ULL,
  4095ULL, 8191ULL, 16383ULL, 32767ULL,
  65535ULL, 131071ULL, 262143ULL, 524287ULL,
  1048575ULL, 2097151ULL, 4194303ULL, 8388607ULL,
  16777215ULL, 33554431ULL, 67108863ULL, 134217727ULL,
  268435455ULL, 536870911ULL, 1073741823ULL, 2147483647ULL,
  4294967295ULL, 8589934591ULL, 17179869183ULL, 34359738367ULL,
  68719476735ULL, 137438953471ULL, 274877906943ULL, 549755813887ULL,
  1099511627775ULL, 2199023255551ULL, 4398046511103ULL, 8796093022207ULL,
  17592186044415ULL, 35184372088831ULL, 70368744177663ULL, 140737488355327ULL,
  281474976710655ULL, 562949953421311ULL, 1125899906842623ULL, 2251799813685247ULL,
  4503599627370495ULL, 9007199254740991ULL, 18014398509481983ULL, 36028797018963967ULL,
  72057594037927935ULL, 144115188075855871ULL, 288230376151711743ULL, 576460752303423487ULL,
  1152921504606846975ULL, 2305843009213693951ULL, 4611686018427387903ULL, 9223372036854775807ULL
}};


//============================================================================
// QubitVector class
//============================================================================

// Template class for qubit vector.
// The arguement of the template must have an operator[] access method.
// The following methods may also need to be template specialized:
//   * set_num_qubits(size_t)
//   * initialize()
//   * initialize_from_vector(cvector_t<data_t>)
// If the template argument does not have these methods then template
// specialization must be used to override the default implementations.

template <typename data_t = double>
class QubitVector {

public:

  //-----------------------------------------------------------------------
  // Constructors and Destructor
  //-----------------------------------------------------------------------

  QubitVector();
  explicit QubitVector(size_t num_qubits);
  virtual ~QubitVector();
  QubitVector(const QubitVector& obj) = delete;
  QubitVector &operator=(const QubitVector& obj) = delete;

  //-----------------------------------------------------------------------
  // Data access
  //-----------------------------------------------------------------------

  // Element access
  std::complex<data_t> &operator[](uint_t element);
  std::complex<data_t> operator[](uint_t element) const;

  // Returns a reference to the underlying data_t data class
  std::complex<data_t>* &data() {return data_;}

  // Returns a copy of the underlying data_t data class
  std::complex<data_t>* data() const {return data_;}

  //-----------------------------------------------------------------------
  // Utility functions
  //-----------------------------------------------------------------------

  // Set the size of the vector in terms of qubit number
  void set_num_qubits(size_t num_qubits);

  // Returns the number of qubits for the current vector
  virtual uint_t num_qubits() const {return num_qubits_;}

  // Returns the size of the underlying n-qubit vector
  uint_t size() const {return data_size_;}

  // Returns required memory
  size_t required_memory_mb(uint_t num_qubits) const;

  // Returns a copy of the underlying data_t data as a complex vector
  cvector_t<data_t> vector() const;

  // Return JSON serialization of QubitVector;
  json_t json() const;

  // Set all entries in the vector to 0.
  void zero();

  // convert vector type to data type of this qubit vector
  cvector_t<data_t> convert(const cvector_t<double>& v) const;

  // index0 returns the integer representation of a number of bits set
  // to zero inserted into an arbitrary bit string.
  // Eg: for qubits 0,2 in a state k = ba ( ba = 00 => k=0, etc).
  // indexes0([1], k) -> int(b0a)
  // indexes0([1,3], k) -> int(0b0a)
  // Example: k = 77  = 1001101 , qubits_sorted = [1,4]
  // ==> output = 297 = 100101001 (with 0's put into places 1 and 4).
  template<typename list_t>
  uint_t index0(const list_t &qubits_sorted, const uint_t k) const;

  // Return a std::unique_ptr to an array of of 2^N in ints
  // each int corresponds to an N qubit bitstring for M-N qubit bits in state k,
  // and the specified N qubits in states [0, ..., 2^N - 1]
  // qubits_sorted must be sorted lowest to highest. Eg. {0, 1}.
  // qubits specifies the location of the qubits in the returned strings.
  // NOTE: since the return is a unique_ptr it cannot be copied.
  // indexes returns the array of all bit values for the specified qubits
  // (Eg: for qubits 0,2 in a state k = ba:
  // indexes([1], [1], k) = [int(b0a), int(b1a)],
  // if it were two qubits inserted say at 1,3 it would be:
  // indexes([1,3], [1,3], k) -> [int(0b0a), int(0b1a), int(1b0a), (1b1a)]
  // If the qubits were passed in reverse order it would swap qubit position in the list:
  // indexes([3,1], [1,3], k) -> [int(0b0a), int(1b0a), int(0b1a), (1b1a)]
  // Example: k=77, qubits=qubits_sorted=[1,4] ==> output=[297,299,313,315]
  // input: k = 77  = 1001101
  // output[0]: 297 = 100101001 (with 0's put into places 1 and 4).
  // output[1]: 299 = 100101011 (with 0 put into place 1, and 1 put into place 4).
  // output[2]: 313 = 100111001 (with 1 put into place 1, and 0 put into place 4).
  // output[3]: 313 = 100111011 (with 1's put into places 1 and 4).
  indexes_t indexes(const reg_t &qubits, const reg_t &qubits_sorted, const uint_t k) const;

  // As above but returns a fixed sized array of of 2^N in ints
  template<size_t N>
  areg_t<1ULL << N> indexes(const areg_t<N> &qs, const areg_t<N> &qubits_sorted, const uint_t k) const;

  // State initialization of a component
  // Initialize the specified qubits to a desired statevector
  // (leaving the other qubits in their current state)
  // assuming the qubits being initialized have already been reset to the zero state
  // (using apply_reset)
  void initialize_component(const reg_t &qubits, const cvector_t<double> &state);

  //-----------------------------------------------------------------------
  // Check point operations
  //-----------------------------------------------------------------------

  // Create a checkpoint of the current state
  void checkpoint();

  // Revert to the checkpoint
  void revert(bool keep);

  // Compute the inner product of current state with checkpoint state
  std::complex<double> inner_product() const;

  //-----------------------------------------------------------------------
  // Initialization
  //-----------------------------------------------------------------------

  // Initializes the current vector so that all qubits are in the |0> state.
  void initialize();

  // Initializes the vector to a custom initial state.
  // If the length of the data vector does not match the number of qubits
  // an exception is raised.
  void initialize_from_vector(const cvector_t<double> &data);

  // Initializes the vector to a custom initial state.
  // If num_states does not match the number of qubits an exception is raised.
  void initialize_from_data(const std::complex<data_t>* data, const size_t num_states);

  //-----------------------------------------------------------------------
  // Apply Matrices
  //-----------------------------------------------------------------------

  // Apply a 1-qubit matrix to the state vector.
  // The matrix is input as vector of the column-major vectorized 1-qubit matrix.
  void apply_matrix(const uint_t qubit, const cvector_t<double> &mat);

  // Apply a N-qubit matrix to the state vector.
  // The matrix is input as vector of the column-major vectorized N-qubit matrix.
  void apply_matrix(const reg_t &qubits, const cvector_t<double> &mat);

  // Apply a stacked set of 2^control_count target_count--qubit matrix to the state vector.
  // The matrix is input as vector of the column-major vectorized N-qubit matrix.
  void apply_multiplexer(const reg_t &control_qubits, const reg_t &target_qubits, const cvector_t<double> &mat);

  // Apply a 1-qubit diagonal matrix to the state vector.
  // The matrix is input as vector of the matrix diagonal.
  void apply_diagonal_matrix(const uint_t qubit, const cvector_t<double> &mat);

  // Apply a N-qubit diagonal matrix to the state vector.
  // The matrix is input as vector of the matrix diagonal.
  void apply_diagonal_matrix(const reg_t &qubits, const cvector_t<double> &mat);
  
  // Swap pairs of indicies in the underlying vector
  void apply_permutation_matrix(const reg_t &qubits,
                                const std::vector<std::pair<uint_t, uint_t>> &pairs);

#ifdef SIMD_PPC64LE
  // Apply a 1-qubit matrix to the state vector.
  // The matrix is input as vector of the column-major vectorized 1-qubit matrix.
  void apply_matrix_ppc64le(const uint_t qubit, const cvector_t<double> &mat);

  // Apply a N-qubit matrix to the state vector.
  // The matrix is input as vector of the column-major vectorized N-qubit matrix.
  void apply_matrix_ppc64le(const reg_t &qubits, const cvector_t<double> &mat);
  
  // Apply a N-qubit diagonal matrix to the state vector.
  // The matrix is input as vector of the matrix diagonal.
  void apply_diagonal_matrix_ppc64le(const reg_t &qubits, const cvector_t<double> &mat);
#endif //SIMD_PPC64LE
  //-----------------------------------------------------------------------
  // Apply Specialized Gates
  //-----------------------------------------------------------------------

  // Apply a general N-qubit multi-controlled X-gate
  // If N=1 this implements an optimized X gate
  // If N=2 this implements an optimized CX gate
  // If N=3 this implements an optimized Toffoli gate
  void apply_mcx(const reg_t &qubits);

  // Apply a general multi-controlled Y-gate
  // If N=1 this implements an optimized Y gate
  // If N=2 this implements an optimized CY gate
  // If N=3 this implements an optimized CCY gate
  void apply_mcy(const reg_t &qubits);
  
  // Apply a general multi-controlled single-qubit phase gate
  // with diagonal [1, ..., 1, phase]
  // If N=1 this implements an optimized single-qubit phase gate
  // If N=2 this implements an optimized CPhase gate
  // If N=3 this implements an optimized CCPhase gate
  // if phase = -1 this is a Z, CZ, CCZ gate
  void apply_mcphase(const reg_t &qubits, const std::complex<double> phase);

  // Apply a general multi-controlled single-qubit unitary gate
  // If N=1 this implements an optimized single-qubit U gate
  // If N=2 this implements an optimized CU gate
  // If N=3 this implements an optimized CCU gate
  void apply_mcu(const reg_t &qubits, const cvector_t<double> &mat);

  // Apply a general multi-controlled SWAP gate
  // If N=2 this implements an optimized SWAP  gate
  // If N=3 this implements an optimized Fredkin gate
  void apply_mcswap(const reg_t &qubits);

  //-----------------------------------------------------------------------
  // Z-measurement outcome probabilities
  //-----------------------------------------------------------------------

  // Return the Z-basis measurement outcome probability P(outcome) for
  // outcome in [0, 2^num_qubits - 1]
  virtual double probability(const uint_t outcome) const;

  // Return the probabilities for all measurement outcomes in the current vector
  // This is equivalent to returning a new vector with  new[i]=|orig[i]|^2.
  // Eg. For 2-qubits this is [P(00), P(01), P(010), P(11)]
  virtual std::vector<double> probabilities() const;

  // Return the Z-basis measurement outcome probabilities [P(0), ..., P(2^N-1)]
  // for measurement of N-qubits.
  virtual std::vector<double> probabilities(const reg_t &qubits) const;

  // Return M sampled outcomes for Z-basis measurement of all qubits
  // The input is a length M list of random reals between [0, 1) used for
  // generating samples.
  virtual reg_t sample_measure(const std::vector<double> &rnds) const;

  //-----------------------------------------------------------------------
  // Norms
  //-----------------------------------------------------------------------
  
  // Returns the norm of the current vector
  double norm() const;

  // These functions return the norm <psi|A^dagger.A|psi> obtained by
  // applying a matrix A to the vector. It is equivalent to returning the
  // expectation value of A^\dagger A, and could probably be removed because
  // of this.

  // Return the norm for of the vector obtained after apply the 1-qubit
  // matrix mat to the vector.
  // The matrix is input as vector of the column-major vectorized 1-qubit matrix.
  double norm(const uint_t qubit, const cvector_t<double> &mat) const;

  // Return the norm for of the vector obtained after apply the N-qubit
  // matrix mat to the vector.
  // The matrix is input as vector of the column-major vectorized N-qubit matrix.
  double norm(const reg_t &qubits, const cvector_t<double> &mat) const;

  // Return the norm for of the vector obtained after apply the 1-qubit
  // diagonal matrix mat to the vector.
  // The matrix is input as vector of the matrix diagonal.
  double norm_diagonal(const uint_t qubit, const cvector_t<double> &mat) const;

  // Return the norm for of the vector obtained after apply the N-qubit
  // diagonal matrix mat to the vector.
  // The matrix is input as vector of the matrix diagonal.
  double norm_diagonal(const reg_t &qubits, const cvector_t<double> &mat) const;

  //-----------------------------------------------------------------------
  // JSON configuration settings
  //-----------------------------------------------------------------------

  // Set the threshold for chopping values to 0 in JSON
  void set_json_chop_threshold(double threshold);

  // Set the threshold for chopping values to 0 in JSON
  double get_json_chop_threshold() {return json_chop_threshold_;}

  //-----------------------------------------------------------------------
  // OpenMP configuration settings
  //-----------------------------------------------------------------------

  // Set the maximum number of OpenMP thread for operations.
  void set_omp_threads(int n);

  // Get the maximum number of OpenMP thread for operations.
  uint_t get_omp_threads() {return omp_threads_;}

  // Set the qubit threshold for activating OpenMP.
  // If self.qubits() > threshold OpenMP will be activated.
  void set_omp_threshold(int n);

  // Get the qubit threshold for activating OpenMP.
  uint_t get_omp_threshold() {return omp_threshold_;}

  //-----------------------------------------------------------------------
  // Optimization configuration settings
  //-----------------------------------------------------------------------

  // Set the sample_measure index size
  void set_sample_measure_index_size(int n) {sample_measure_index_size_ = n;}

  // Get the sample_measure index size
  int get_sample_measure_index_size() {return sample_measure_index_size_;}

protected:

  //-----------------------------------------------------------------------
  // Protected data members
  //-----------------------------------------------------------------------
  size_t num_qubits_;
  size_t data_size_;
  std::complex<data_t>* data_;
  std::complex<data_t>* checkpoint_;

  // If data_t is double, SIMD function is enable
#ifdef SIMD_PPC64LE
  bool type_double = ((typeid(data_t) == typeid(double)) ? true : false);
#endif

  //-----------------------------------------------------------------------
  // Config settings
  //----------------------------------------------------------------------- 
  uint_t omp_threads_ = 1;     // Disable multithreading by default
  uint_t omp_threshold_ = 14;  // Qubit threshold for multithreading when enabled
  int sample_measure_index_size_ = 10; // Sample measure indexing qubit size
  double json_chop_threshold_ = 0;  // Threshold for choping small values
                                    // in JSON serialization

  //-----------------------------------------------------------------------
  // Error Messages
  //-----------------------------------------------------------------------

  void check_qubit(const uint_t qubit) const;
  void check_vector(const cvector_t<data_t> &diag, uint_t nqubits) const;
  void check_matrix(const cvector_t<data_t> &mat, uint_t nqubits) const;
  void check_dimension(const QubitVector &qv) const;
  void check_checkpoint() const;

  //-----------------------------------------------------------------------
  // Statevector update with Lambda function
  //-----------------------------------------------------------------------
  // Apply a lambda function to all entries of the statevector.
  // The function signature should be:
  //
  // [&](const int_t k)->void
  //
  // where k is the index of the vector
  template <typename Lambda>
  void apply_lambda(Lambda&& func);

  //-----------------------------------------------------------------------
  // Statevector block update with Lambda function
  //-----------------------------------------------------------------------
  // These functions loop through the indexes of the qubitvector data and
  // apply a lambda function to each block specified by the qubits argument.
  //
  // NOTE: The lambda functions can use the dynamic or static indexes
  // signature however if N is known at compile time the static case should
  // be preferred as it is significantly faster.

  // Apply a N-qubit lambda function to all blocks of the statevector
  // for the given qubits. The function signature should be either:
  //
  // (Static): [&](const areg_t<1ULL<<N> &inds)->void
  // (Dynamic): [&](const indexes_t &inds)->void
  //
  // where `inds` are the 2 ** N indexes for each N-qubit block returned by
  // the `indexes` function.
  template <typename Lambda, typename list_t>
  void apply_lambda(Lambda&& func, const list_t &qubits);

  // Apply an N-qubit parameterized lambda function to all blocks of the
  // statevector for the given qubits. The function signature should be:
  //
  // (Static): [&](const areg_t<1ULL<<N> &inds, const param_t &params)->void
  // (Dynamic): [&](const indexes_t &inds, const param_t &params)->void
  //
  // where `inds` are the 2 ** N indexes for each N-qubit block returned by
  // the `indexes` function and `param` is a templated parameter class.
  // (typically a complex vector).
  template <typename Lambda, typename list_t, typename param_t>
  void apply_lambda(Lambda&& func, const list_t &qubits, const param_t &par);

  //-----------------------------------------------------------------------
  // State reduction with Lambda functions
  //-----------------------------------------------------------------------
  // Apply a complex reduction lambda function to all entries of the
  // statevector and return the complex result.
  // The function signature should be:
  //
  // [&](const int_t k, double &val_re, double &val_im)->void
  //
  // where k is the index of the vector, val_re and val_im are the doubles
  // to store the reduction.
  // Returns std::complex<double>(val_re, val_im)
  template <typename Lambda>
  std::complex<double> apply_reduction_lambda(Lambda&& func) const;

  //-----------------------------------------------------------------------
  // Statevector block reduction with Lambda function
  //-----------------------------------------------------------------------
  // These functions loop through the indexes of the qubitvector data and
  // apply a reduction lambda function to each block specified by the qubits
  // argument. The reduction lambda stores the reduction in two doubles
  // (val_re, val_im) and returns the complex result std::complex<double>(val_re, val_im)
  //
  // NOTE: The lambda functions can use the dynamic or static indexes
  // signature however if N is known at compile time the static case should
  // be preferred as it is significantly faster.

  // Apply a N-qubit complex matrix reduction lambda function to all blocks
  // of the statevector for the given qubits.
  // The lambda function signature should be:
  //
  // (Static): [&](const areg_t<1ULL<<N> &inds, const param_t &mat,
  //               double &val_re, double &val_im)->void
  // (Dynamic): [&](const indexes_t &inds, const param_t &mat,
  //                double &val_re, double &val_im)->void
  //
  // where `inds` are the 2 ** N indexes for each N-qubit block returned by
  // the `indexes` function, `val_re` and `val_im` are the doubles to
  // store the reduction returned as std::complex<double>(val_re, val_im).
  template <typename Lambda, typename list_t>
  std::complex<double> apply_reduction_lambda(Lambda&& func,
                                              const list_t &qubits) const;

  // Apply a N-qubit complex matrix reduction lambda function to all blocks
  // of the statevector for the given qubits.
  // The lambda function signature should be:
  //
  // (Static): [&](const areg_t<1ULL<<N> &inds, const param_t &parms,
  //               double &val_re, double &val_im)->void
  // (Dynamic): [&](const indexes_t &inds, const param_t &params,
  //                double &val_re, double &val_im)->void
  //
  // where `inds` are the 2 ** N indexes for each N-qubit block returned by
  // the `indexe`s function, `params` is a templated parameter class
  // (typically a complex vector), `val_re` and `val_im` are the doubles to
  // store the reduction returned as std::complex<double>(val_re, val_im).
  template <typename Lambda, typename list_t, typename param_t>
  std::complex<double> apply_reduction_lambda(Lambda&& func,
                                              const list_t &qubits,
                                              const param_t &params) const;
};

/*******************************************************************************
 *
 * Implementations
 *
 ******************************************************************************/

//------------------------------------------------------------------------------
// JSON Serialization
//------------------------------------------------------------------------------

template <typename data_t>
inline void to_json(json_t &js, const QubitVector<data_t> &qv) {
  js = qv.json();
}

template <typename data_t>
json_t QubitVector<data_t>::json() const {
  const int_t END = data_size_;
  const json_t ZERO = std::complex<data_t>(0.0, 0.0);
  json_t js = json_t(data_size_, ZERO);

  if (json_chop_threshold_ > 0) {
    #pragma omp parallel for if (num_qubits_ > omp_threshold_ && omp_threads_ > 1) num_threads(omp_threads_)
    for (int_t j=0; j < END; j++) {
      if (std::abs(data_[j].real()) > json_chop_threshold_)
        js[j][0] = data_[j].real();
      if (std::abs(data_[j].imag()) > json_chop_threshold_)
        js[j][1] = data_[j].imag();
    }
  } else {
    #pragma omp parallel for if (num_qubits_ > omp_threshold_ && omp_threads_ > 1) num_threads(omp_threads_)
    for (int_t j=0; j < END; j++) {
      js[j][0] = data_[j].real();
      js[j][1] = data_[j].imag();
    }
  }
  return js;
}

//------------------------------------------------------------------------------
// Error Handling
//------------------------------------------------------------------------------

template <typename data_t>
void QubitVector<data_t>::check_qubit(const uint_t qubit) const {
  if (qubit + 1 > num_qubits_) {
    std::string error = "QubitVector: qubit index " + std::to_string(qubit) +
                        " > " + std::to_string(num_qubits_);
    throw std::runtime_error(error);
  }
}

template <typename data_t>
void QubitVector<data_t>::check_matrix(const cvector_t<data_t> &vec, uint_t nqubits) const {
  const size_t DIM = BITS[nqubits];
  const auto SIZE = vec.size();
  if (SIZE != DIM * DIM) {
    std::string error = "QubitVector: vector size is " + std::to_string(SIZE) +
                        " != " + std::to_string(DIM * DIM);
    throw std::runtime_error(error);
  }
}

template <typename data_t>
void QubitVector<data_t>::check_vector(const cvector_t<data_t> &vec, uint_t nqubits) const {
  const size_t DIM = BITS[nqubits];
  const auto SIZE = vec.size();
  if (SIZE != DIM) {
    std::string error = "QubitVector: vector size is " + std::to_string(SIZE) +
                        " != " + std::to_string(DIM);
    throw std::runtime_error(error);
  }
}

template <typename data_t>
void QubitVector<data_t>::check_dimension(const QubitVector &qv) const {
  if (data_size_ != qv.size_) {
    std::string error = "QubitVector: vectors are different shape " +
                         std::to_string(data_size_) + " != " +
                         std::to_string(qv.num_states_);
    throw std::runtime_error(error);
  }
}

template <typename data_t>
void QubitVector<data_t>::check_checkpoint() const {
  if (!checkpoint_) {
    throw std::runtime_error("QubitVector: checkpoint must exist for inner_product() or revert()");
  }
}

//------------------------------------------------------------------------------
// Constructors & Destructor
//------------------------------------------------------------------------------

template <typename data_t>
QubitVector<data_t>::QubitVector(size_t num_qubits) : num_qubits_(0), data_(nullptr), checkpoint_(0){
  set_num_qubits(num_qubits);
}

template <typename data_t>
QubitVector<data_t>::QubitVector() : QubitVector(0) {}

template <typename data_t>
QubitVector<data_t>::~QubitVector() {
  if (data_)
    free(data_);

  if (checkpoint_)
    free(checkpoint_);
}

//------------------------------------------------------------------------------
// Element access operators
//------------------------------------------------------------------------------

template <typename data_t>
std::complex<data_t> &QubitVector<data_t>::operator[](uint_t element) {
  // Error checking
  #ifdef DEBUG
  if (element > data_size_) {
    std::string error = "QubitVector: vector index " + std::to_string(element) +
                        " > " + std::to_string(data_size_);
    throw std::runtime_error(error);
  }
  #endif
  return data_[element];
}

template <typename data_t>
std::complex<data_t> QubitVector<data_t>::operator[](uint_t element) const {
  // Error checking
  #ifdef DEBUG
  if (element > data_size_) {
    std::string error = "QubitVector: vector index " + std::to_string(element) +
                        " > " + std::to_string(data_size_);
    throw std::runtime_error(error);
  }
  #endif
  return data_[element];
}

template <typename data_t>
cvector_t<data_t> QubitVector<data_t>::vector() const {
  cvector_t<data_t> ret(data_size_, 0.);
  const int_t END = data_size_;
  #pragma omp parallel for if (num_qubits_ > omp_threshold_ && omp_threads_ > 1) num_threads(omp_threads_)
  for (int_t j=0; j < END; j++) {
    ret[j] = data_[j];
  }
  return ret;
}

//------------------------------------------------------------------------------
// Indexing
//------------------------------------------------------------------------------

template <typename data_t>
template <typename list_t>
uint_t QubitVector<data_t>::index0(const list_t &qubits_sorted, const uint_t k) const {
  uint_t lowbits, retval = k;
  for (size_t j = 0; j < qubits_sorted.size(); j++) {
    lowbits = retval & MASKS[qubits_sorted[j]];
    retval >>= qubits_sorted[j];
    retval <<= qubits_sorted[j] + 1;
    retval |= lowbits;
  }
  return retval;
}

template <typename data_t>
template <size_t N>
areg_t<1ULL << N> QubitVector<data_t>::indexes(const areg_t<N> &qs,
                                               const areg_t<N> &qubits_sorted,
                                               const uint_t k) const {
  areg_t<1ULL << N> ret;
  ret[0] = index0(qubits_sorted, k);
  for (size_t i = 0; i < N; i++) {
    const auto n = BITS[i];
    const auto bit = BITS[qs[i]];
    for (size_t j = 0; j < n; j++)
      ret[n + j] = ret[j] | bit;
  }
  return ret;
}

template <typename data_t>
indexes_t QubitVector<data_t>::indexes(const reg_t& qubits,
                                       const reg_t& qubits_sorted,
                                       const uint_t k) const {
  const auto N = qubits_sorted.size();
  indexes_t ret(new uint_t[BITS[N]]);
  // Get index0
  ret[0] = index0(qubits_sorted, k);
  for (size_t i = 0; i < N; i++) {
    const auto n = BITS[i];
    const auto bit = BITS[qubits[i]];
    for (size_t j = 0; j < n; j++)
      ret[n + j] = ret[j] | bit;
  }
  return ret;
}

//------------------------------------------------------------------------------
// State initialize component
//------------------------------------------------------------------------------
template <typename data_t>
void QubitVector<data_t>::initialize_component(const reg_t &qubits, const cvector_t<double> &state0) {

  cvector_t<data_t> state = convert(state0);

  // Lambda function for initializing component
  const size_t N = qubits.size();
  auto lambda = [&](const indexes_t &inds, const cvector_t<data_t> &_state)->void {
    const uint_t DIM = 1ULL << N;
    std::complex<data_t> cache = data_[inds[0]];  // the k-th component of non-initialized vector
    for (size_t i = 0; i < DIM; i++) {
      data_[inds[i]] = cache * _state[i];  // set component to psi[k] * state[i]
    }    // (where psi is is the post-reset state of the non-initialized qubits)
   };
  // Use the lambda function
  apply_lambda(lambda, qubits, state);
}

//------------------------------------------------------------------------------
// Utility
//------------------------------------------------------------------------------

template <typename data_t>
void QubitVector<data_t>::zero() {
  const int_t END = data_size_;    // end for k loop

#pragma omp parallel for if (num_qubits_ > omp_threshold_ && omp_threads_ > 1) num_threads(omp_threads_)
  for (int_t k = 0; k < END; ++k) {
    data_[k] = 0.0;
  }
}

template <typename data_t>
cvector_t<data_t> QubitVector<data_t>::convert(const cvector_t<double>& v) const {
  cvector_t<data_t> ret(v.size());
  for (size_t i = 0; i < v.size(); ++i)
    ret[i] = v[i];
  return ret;
}


template <typename data_t>
void QubitVector<data_t>::set_num_qubits(size_t num_qubits) {

  size_t prev_num_qubits = num_qubits_;
  num_qubits_ = num_qubits;
  data_size_ = BITS[num_qubits];

  if (checkpoint_) {
    free(checkpoint_);
    checkpoint_ = nullptr;
  }

  // Free any currently assigned memory
  if (data_) {
    if (prev_num_qubits != num_qubits_) {
      free(data_);
      data_ = nullptr;
    }
  }

  // Allocate memory for new vector
  if (data_ == nullptr)
    data_ = reinterpret_cast<std::complex<data_t>*>(malloc(sizeof(std::complex<data_t>) * data_size_));
}

template <typename data_t>
size_t QubitVector<data_t>::required_memory_mb(uint_t num_qubits) const {

  size_t unit = std::log2(sizeof(std::complex<data_t>));
  size_t shift_mb = std::max<int_t>(0, num_qubits + unit - 20);
  size_t mem_mb = 1ULL << shift_mb;
  return mem_mb;
}


template <typename data_t>
void QubitVector<data_t>::checkpoint() {
  if (!checkpoint_)
    checkpoint_ = reinterpret_cast<std::complex<data_t>*>(malloc(sizeof(std::complex<data_t>) * data_size_));

  const int_t END = data_size_;    // end for k loop
  #pragma omp parallel for if (num_qubits_ > omp_threshold_ && omp_threads_ > 1) num_threads(omp_threads_)
  for (int_t k = 0; k < END; ++k)
    checkpoint_[k] = data_[k];
}


template <typename data_t>
void QubitVector<data_t>::revert(bool keep) {

  #ifdef DEBUG
  check_checkpoint();
  #endif

  // If we aren't keeping checkpoint we don't need to copy memory
  // we can simply swap the pointers and free discarded memory
  if (!keep) {
    free(data_);
    data_ = checkpoint_;
    checkpoint_ = nullptr;
  } else {
    // Otherwise we need to copy data
    const int_t END = data_size_;    // end for k loop
#pragma omp parallel for if (num_qubits_ > omp_threshold_ && omp_threads_ > 1) num_threads(omp_threads_)
    for (int_t k = 0; k < END; ++k)
      data_[k] = checkpoint_[k];
  }
}

template <typename data_t>
std::complex<double> QubitVector<data_t>::inner_product() const {

  #ifdef DEBUG
  check_checkpoint();
  #endif
  // Lambda function for inner product with checkpoint state
  auto lambda = [&](int_t k, double &val_re, double &val_im)->void {
    const std::complex<double> z = data_[k] * std::conj(checkpoint_[k]);
    val_re += std::real(z);
    val_im += std::imag(z);
  };
  return apply_reduction_lambda(lambda);
}

//------------------------------------------------------------------------------
// Initialization
//------------------------------------------------------------------------------

template <typename data_t>
void QubitVector<data_t>::initialize() {
  zero();
  data_[0] = 1.;
}

template <typename data_t>
void QubitVector<data_t>::initialize_from_vector(const cvector_t<double> &statevec) {
  if (data_size_ != statevec.size()) {
    std::string error = "QubitVector::initialize input vector is incorrect length (" + 
                        std::to_string(data_size_) + "!=" +
                        std::to_string(statevec.size()) + ")";
    throw std::runtime_error(error);
  }

  const int_t END = data_size_;    // end for k loop

#pragma omp parallel for if (num_qubits_ > omp_threshold_ && omp_threads_ > 1) num_threads(omp_threads_)
  for (int_t k = 0; k < END; ++k)
    data_[k] = statevec[k];
}

template <typename data_t>
void QubitVector<data_t>::initialize_from_data(const std::complex<data_t>* statevec, const size_t num_states) {
  if (data_size_ != num_states) {
    std::string error = "QubitVector::initialize input vector is incorrect length (" +
                        std::to_string(data_size_) + "!=" + std::to_string(num_states) + ")";
    throw std::runtime_error(error);
  }

  const int_t END = data_size_;    // end for k loop

#pragma omp parallel for if (num_qubits_ > omp_threshold_ && omp_threads_ > 1) num_threads(omp_threads_)
  for (int_t k = 0; k < END; ++k)
    data_[k] = statevec[k];
}


/*******************************************************************************
 *
 * CONFIG SETTINGS
 *
 ******************************************************************************/

template <typename data_t>
void QubitVector<data_t>::set_omp_threads(int n) {
  if (n > 0)
    omp_threads_ = n;
}

template <typename data_t>
void QubitVector<data_t>::set_omp_threshold(int n) {
  if (n > 0)
    omp_threshold_ = n;
}

template <typename data_t>
void QubitVector<data_t>::set_json_chop_threshold(double threshold) {
  json_chop_threshold_ = threshold;
}

/*******************************************************************************
 *
 * LAMBDA FUNCTION TEMPLATES
 *
 ******************************************************************************/


//------------------------------------------------------------------------------
// State update
//------------------------------------------------------------------------------

template <typename data_t>
template<typename Lambda>
void QubitVector<data_t>::apply_lambda(Lambda&& func) {
  const int_t END = data_size_;
#pragma omp parallel if (num_qubits_ > omp_threshold_ && omp_threads_ > 1) num_threads(omp_threads_)
  {
#pragma omp for
    for (int_t k = 0; k < END; k++) {
      std::forward<Lambda>(func)(k);
    }
  }
}

template <typename data_t>
template<typename Lambda, typename list_t>
void QubitVector<data_t>::apply_lambda(Lambda&& func, const list_t &qubits) {

  // Error checking
  #ifdef DEBUG
  for (const auto &qubit : qubits)
    check_qubit(qubit);
  #endif

  const auto NUM_QUBITS = qubits.size();
  const int_t END = data_size_ >> NUM_QUBITS;
  auto qubits_sorted = qubits;
  std::sort(qubits_sorted.begin(), qubits_sorted.end());
#pragma omp parallel if (num_qubits_ > omp_threshold_ && omp_threads_ > 1) num_threads(omp_threads_)
  {
#pragma omp for
    for (int_t k = 0; k < END; k++) {
      // store entries touched by U
      const auto inds = indexes(qubits, qubits_sorted, k);
      std::forward<Lambda>(func)(inds);
    }
  }
}

template <typename data_t>
template<typename Lambda, typename list_t, typename param_t>
void QubitVector<data_t>::apply_lambda(Lambda&& func,
                                       const list_t &qubits,
                                       const param_t &params) {

  // Error checking
  #ifdef DEBUG
  for (const auto &qubit : qubits)
    check_qubit(qubit);
  #endif

  const auto NUM_QUBITS = qubits.size();
  const int_t END = data_size_ >> NUM_QUBITS;
  auto qubits_sorted = qubits;
  std::sort(qubits_sorted.begin(), qubits_sorted.end());

#pragma omp parallel if (num_qubits_ > omp_threshold_ && omp_threads_ > 1) num_threads(omp_threads_)
  {
#pragma omp for
    for (int_t k = 0; k < END; k++) {
      const auto inds = indexes(qubits, qubits_sorted, k);
      std::forward<Lambda>(func)(inds, params);
    }
  }
}


//------------------------------------------------------------------------------
// Reduction Lambda
//------------------------------------------------------------------------------

template <typename data_t>
template<typename Lambda>
std::complex<double> QubitVector<data_t>::apply_reduction_lambda(Lambda &&func) const {
  // Reduction variables
  double val_re = 0.;
  double val_im = 0.;
  const int_t END = data_size_;
#pragma omp parallel reduction(+:val_re, val_im) if (num_qubits_ > omp_threshold_ && omp_threads_ > 1)         \
                                               num_threads(omp_threads_)
  {
#pragma omp for
    for (int_t k = 0; k < END; k++) {
        std::forward<Lambda>(func)(k, val_re, val_im);
      }
  } // end omp parallel
  return std::complex<double>(val_re, val_im);
}


template <typename data_t>
template<typename Lambda, typename list_t>
std::complex<double>
QubitVector<data_t>::apply_reduction_lambda(Lambda&& func,
                                            const list_t &qubits) const {

  // Error checking
  #ifdef DEBUG
  for (const auto &qubit : qubits)
    check_qubit(qubit);
  #endif

  const size_t NUM_QUBITS =  qubits.size();
  const int_t END = data_size_ >> NUM_QUBITS;
  auto qubits_sorted = qubits;
  std::sort(qubits_sorted.begin(), qubits_sorted.end());

  // Reduction variables
  double val_re = 0.;
  double val_im = 0.;
#pragma omp parallel reduction(+:val_re, val_im) if (num_qubits_ > omp_threshold_ && omp_threads_ > 1)         \
                                               num_threads(omp_threads_)
  {
#pragma omp for
    for (int_t k = 0; k < END; k++) {
      const auto inds = indexes(qubits, qubits_sorted, k);
      std::forward<Lambda>(func)(inds, val_re, val_im);
    }
  } // end omp parallel
  return std::complex<double>(val_re, val_im);
}


template <typename data_t>
template<typename Lambda, typename list_t, typename param_t>
std::complex<double>
QubitVector<data_t>::apply_reduction_lambda(Lambda&& func,
                                            const list_t &qubits,
                                            const param_t &params) const {

  const auto NUM_QUBITS = qubits.size();
  // Error checking
  #ifdef DEBUG
  for (const auto &qubit : qubits)
    check_qubit(qubit);
  #endif

  const int_t END = data_size_ >> NUM_QUBITS;
  auto qubits_sorted = qubits;
  std::sort(qubits_sorted.begin(), qubits_sorted.end());

  // Reduction variables
  double val_re = 0.;
  double val_im = 0.;
#pragma omp parallel reduction(+:val_re, val_im) if (num_qubits_ > omp_threshold_ && omp_threads_ > 1)         \
                                               num_threads(omp_threads_)
  {
#pragma omp for
    for (int_t k = 0; k < END; k++) {
      const auto inds = indexes(qubits, qubits_sorted, k);
      std::forward<Lambda>(func)(inds, params, val_re, val_im);
    }
  } // end omp parallel
  return std::complex<double>(val_re, val_im);
}


/*******************************************************************************
 *
 * MATRIX MULTIPLICATION
 *
 ******************************************************************************/
template <typename data_t>
void QubitVector<data_t>::apply_matrix(const reg_t &qubits,
                                       const cvector_t<double> &mat) {

#ifdef SIMD_PPC64LE
	apply_matrix_ppc64le(qubits, mat);
#else
  const size_t N = qubits.size();
  // Error checking
  #ifdef DEBUG
  check_vector(mat, 2 * N);
  #endif

  // Static array optimized lambda functions
  switch (N) {
    case 1:
      apply_matrix(qubits[0], mat);
      return;
    case 2: {
      // Lambda function for 2-qubit matrix multiplication
      auto lambda = [&](const areg_t<4> &inds, const cvector_t<data_t> &_mat)->void {
        std::array<std::complex<data_t>, 4> cache;
        for (size_t i = 0; i < 4; i++) {
          const auto ii = inds[i];
          cache[i] = data_[ii];
          data_[ii] = 0.;
        }
        // update state vector
        for (size_t i = 0; i < 4; i++)
          for (size_t j = 0; j < 4; j++)
            data_[inds[i]] += _mat[i + 4 * j] * cache[j];
      };
      apply_lambda(lambda, areg_t<2>({{qubits[0], qubits[1]}}), convert(mat));
      return;
    }
    case 3: {
      // Lambda function for 3-qubit matrix multiplication
      auto lambda = [&](const areg_t<8> &inds, const cvector_t<data_t> &_mat)->void {
        std::array<std::complex<data_t>, 8> cache;
        for (size_t i = 0; i < 8; i++) {
          const auto ii = inds[i];
          cache[i] = data_[ii];
          data_[ii] = 0.;
        }
        // update state vector
        for (size_t i = 0; i < 8; i++)
          for (size_t j = 0; j < 8; j++)
            data_[inds[i]] += _mat[i + 8 * j] * cache[j];
      };
      apply_lambda(lambda, areg_t<3>({{qubits[0], qubits[1], qubits[2]}}), convert(mat));
      return;
    }
    case 4: {
      // Lambda function for 4-qubit matrix multiplication
      auto lambda = [&](const areg_t<16> &inds, const cvector_t<data_t> &_mat)->void {
        std::array<std::complex<data_t>, 16> cache;
        for (size_t i = 0; i < 16; i++) {
          const auto ii = inds[i];
          cache[i] = data_[ii];
          data_[ii] = 0.;
        }
        // update state vector
        for (size_t i = 0; i < 16; i++)
          for (size_t j = 0; j < 16; j++)
            data_[inds[i]] += _mat[i + 16 * j] * cache[j];
      };
      apply_lambda(lambda, areg_t<4>({{qubits[0], qubits[1], qubits[2], qubits[3]}}), convert(mat));
      return;
    }
    default: {
      const uint_t DIM = BITS[N];
      // Lambda function for N-qubit matrix multiplication
      auto lambda = [&](const indexes_t &inds, const cvector_t<data_t> &_mat)->void {
        auto cache = std::make_unique<std::complex<data_t>[]>(DIM);
        for (size_t i = 0; i < DIM; i++) {
          const auto ii = inds[i];
          cache[i] = data_[ii];
          data_[ii] = 0.;
        }
        // update state vector
        for (size_t i = 0; i < DIM; i++)
          for (size_t j = 0; j < DIM; j++)
            data_[inds[i]] += _mat[i + DIM * j] * cache[j];
      };
      apply_lambda(lambda, qubits, convert(mat));
    }
  } // end switch
#endif //#ifdef SIMD_PPC64LE
}

template <typename data_t>
void QubitVector<data_t>::apply_multiplexer(const reg_t &control_qubits,
                                            const reg_t &target_qubits,
                                            const cvector_t<double>  &mat) {
  
  // General implementation
  const size_t control_count = control_qubits.size();
  const size_t target_count  = target_qubits.size();
  const uint_t DIM = BITS[(target_count+control_count)];
  const uint_t columns = BITS[target_count];
  const uint_t blocks = BITS[control_count];
  // Lambda function for stacked matrix multiplication
  auto lambda = [&](const indexes_t &inds, const cvector_t<data_t> &_mat)->void {
    auto cache = std::make_unique<std::complex<data_t>[]>(DIM);
    for (uint_t i = 0; i < DIM; i++) {
      const auto ii = inds[i];
      cache[i] = data_[ii];
      data_[ii] = 0.;
    }
    // update state vector
    for (uint_t b = 0; b < blocks; b++)
      for (uint_t i = 0; i < columns; i++)
        for (uint_t j = 0; j < columns; j++)
	{
	  data_[inds[i+b*columns]] += _mat[i+b*columns + DIM * j] * cache[b*columns+j];
	}
  };
  
  // Use the lambda function
  auto qubits = target_qubits;
  for (const auto &q : control_qubits) {qubits.push_back(q);}
  apply_lambda(lambda, qubits, convert(mat));
}

template <typename data_t>
void QubitVector<data_t>::apply_diagonal_matrix(const reg_t &qubits,
                                                const cvector_t<double> &diag) {
#ifdef SIMD_PPC64LE
  apply_diagonal_matrix_ppc64le(qubits, diag);
#else
  const int_t N = qubits.size();
  // Error checking
  #ifdef DEBUG
  check_vector(diag, N);
  #endif

  if (N == 1) {
    apply_diagonal_matrix(qubits[0], diag);
    return;
  }

  auto lambda = [&](const areg_t<2> &inds, const cvector_t<data_t> &_diag)->void {
    for (int_t i = 0; i < 2; ++i) {
      const int_t k = inds[i];
      int_t iv = 0;
      for (int_t j = 0; j < N; j++)
        if ((k & (1ULL << qubits[j])) != 0)
          iv += (1 << j);
      if (_diag[iv] != (data_t) 1.0)
        data_[k] *= _diag[iv];
    }
  };
  apply_lambda(lambda, areg_t<1>({{qubits[0]}}), convert(diag));
#endif
}

template <typename data_t>
void QubitVector<data_t>::apply_permutation_matrix(const reg_t& qubits,
                                                   const std::vector<std::pair<uint_t, uint_t>> &pairs) {
  const size_t N = qubits.size();

  // Error checking
  #ifdef DEBUG
  check_vector(diag, N);
  #endif

  switch (N) {
    case 1: {
      // Lambda function for permutation matrix
      auto lambda = [&](const areg_t<2> &inds)->void {
        for (const auto& p : pairs) {
          std::swap(data_[inds[p.first]], data_[inds[p.second]]);
        }
      };
      apply_lambda(lambda, areg_t<1>({{qubits[0]}}));
      return;
    }
    case 2: {
      // Lambda function for permutation matrix
      auto lambda = [&](const areg_t<4> &inds)->void {
        for (const auto& p : pairs) {
          std::swap(data_[inds[p.first]], data_[inds[p.second]]);
        }
      };
      apply_lambda(lambda, areg_t<2>({{qubits[0], qubits[1]}}));
      return;
    }
    case 3: {
      // Lambda function for permutation matrix
      auto lambda = [&](const areg_t<8> &inds)->void {
        for (const auto& p : pairs) {
          std::swap(data_[inds[p.first]], data_[inds[p.second]]);
        }
      };
      apply_lambda(lambda, areg_t<3>({{qubits[0], qubits[1], qubits[2]}}));
      return;
    }
    case 4: {
      // Lambda function for permutation matrix
      auto lambda = [&](const areg_t<16> &inds)->void {
        for (const auto& p : pairs) {
          std::swap(data_[inds[p.first]], data_[inds[p.second]]);
        }
      };
      apply_lambda(lambda, areg_t<4>({{qubits[0], qubits[1], qubits[2], qubits[3]}}));
      return;
    }
    case 5: {
      // Lambda function for permutation matrix
      auto lambda = [&](const areg_t<32> &inds)->void {
        for (const auto& p : pairs) {
          std::swap(data_[inds[p.first]], data_[inds[p.second]]);
        }
      };
      apply_lambda(lambda, areg_t<5>({{qubits[0], qubits[1], qubits[2],
                                       qubits[3], qubits[4]}}));
      return;
    }
    case 6: {
      // Lambda function for permutation matrix
      auto lambda = [&](const areg_t<64> &inds)->void {
        for (const auto& p : pairs) {
          std::swap(data_[inds[p.first]], data_[inds[p.second]]);
        }
      };
      apply_lambda(lambda, areg_t<6>({{qubits[0], qubits[1], qubits[2],
                                       qubits[3], qubits[4], qubits[5]}}));
      return;
    }
    default: {
      // Lambda function for permutation matrix
      auto lambda = [&](const indexes_t &inds)->void {
        for (const auto& p : pairs) {
          std::swap(data_[inds[p.first]], data_[inds[p.second]]);
        }
      };
      // Use the lambda function
      apply_lambda(lambda, qubits);
    }
  } // end switch
}


/*******************************************************************************
 *
 * APPLY OPTIMIZED GATES
 *
 ******************************************************************************/

//------------------------------------------------------------------------------
// Multi-controlled gates
//------------------------------------------------------------------------------

template <typename data_t>
void QubitVector<data_t>::apply_mcx(const reg_t &qubits) {
  // Calculate the permutation positions for the last qubit.
  const size_t N = qubits.size();
  const size_t pos0 = MASKS[N - 1];
  const size_t pos1 = MASKS[N];

  switch (N) {
    case 1: {
      // Lambda function for X gate
      auto lambda = [&](const areg_t<2> &inds)->void {
        std::swap(data_[inds[pos0]], data_[inds[pos1]]);
      };
      apply_lambda(lambda, areg_t<1>({{qubits[0]}}));
      return;
    }
    case 2: {
      // Lambda function for CX gate
      auto lambda = [&](const areg_t<4> &inds)->void {
        std::swap(data_[inds[pos0]], data_[inds[pos1]]);
      };
      apply_lambda(lambda, areg_t<2>({{qubits[0], qubits[1]}}));
      return;
    }
    case 3: {
      // Lambda function for Toffli gate
      auto lambda = [&](const areg_t<8> &inds)->void {
        std::swap(data_[inds[pos0]], data_[inds[pos1]]);
      };
      apply_lambda(lambda, areg_t<3>({{qubits[0], qubits[1], qubits[2]}}));
      return;
    }
    default: {
      // Lambda function for general multi-controlled X gate
      auto lambda = [&](const indexes_t &inds)->void {
        std::swap(data_[inds[pos0]], data_[inds[pos1]]);
      };
      apply_lambda(lambda, qubits);
    }
  } // end switch
}

template <typename data_t>
void QubitVector<data_t>::apply_mcy(const reg_t &qubits) {
  // Calculate the permutation positions for the last qubit.
  const size_t N = qubits.size();
  const size_t pos0 = MASKS[N - 1];
  const size_t pos1 = MASKS[N];
  const std::complex<data_t> I(0., 1.);

  switch (N) {
    case 1: {
      // Lambda function for Y gate
      auto lambda = [&](const areg_t<2> &inds)->void {
        const std::complex<data_t> cache = data_[inds[pos0]];
        data_[inds[pos0]] = -I * data_[inds[pos1]];
        data_[inds[pos1]] = I * cache;
      };
      apply_lambda(lambda, areg_t<1>({{qubits[0]}}));
      return;
    }
    case 2: {
      // Lambda function for CY gate
      auto lambda = [&](const areg_t<4> &inds)->void {
        const std::complex<data_t> cache = data_[inds[pos0]];
        data_[inds[pos0]] = -I * data_[inds[pos1]];
        data_[inds[pos1]] = I * cache;
      };
      apply_lambda(lambda, areg_t<2>({{qubits[0], qubits[1]}}));
      return;
    }
    case 3: {
      // Lambda function for CCY gate
      auto lambda = [&](const areg_t<8> &inds)->void {
        const std::complex<data_t> cache = data_[inds[pos0]];
        data_[inds[pos0]] = -I * data_[inds[pos1]];
        data_[inds[pos1]] = I * cache;
      };
      apply_lambda(lambda, areg_t<3>({{qubits[0], qubits[1], qubits[2]}}));
      return;
    }
    default: {
      // Lambda function for general multi-controlled Y gate
      auto lambda = [&](const indexes_t &inds)->void {
        const std::complex<data_t> cache = data_[inds[pos0]];
        data_[inds[pos0]] = -I * data_[inds[pos1]];
        data_[inds[pos1]] = I * cache;
      };
      apply_lambda(lambda, qubits);
    }
  } // end switch
}

template <typename data_t>
void QubitVector<data_t>::apply_mcswap(const reg_t &qubits) {
  // Calculate the swap positions for the last two qubits.
  // If N = 2 this is just a regular SWAP gate rather than a controlled-SWAP gate.
  const size_t N = qubits.size();
  const size_t pos0 = MASKS[N - 1];
  const size_t pos1 = pos0 + BITS[N - 2];

  switch (N) {
    case 2: {
      // Lambda function for SWAP gate
      auto lambda = [&](const areg_t<4> &inds)->void {
        std::swap(data_[inds[pos0]], data_[inds[pos1]]);
      };
      apply_lambda(lambda, areg_t<2>({{qubits[0], qubits[1]}}));
      return;
    }
    case 3: {
      // Lambda function for C-SWAP gate
      auto lambda = [&](const areg_t<8> &inds)->void {
        std::swap(data_[inds[pos0]], data_[inds[pos1]]);
      };
      apply_lambda(lambda, areg_t<3>({{qubits[0], qubits[1], qubits[2]}}));
      return;
    }
    default: {
      // Lambda function for general multi-controlled SWAP gate
      auto lambda = [&](const indexes_t &inds)->void {
        std::swap(data_[inds[pos0]], data_[inds[pos1]]);
      };
      apply_lambda(lambda, qubits);
    }
  } // end switch
}

template <typename data_t>
void QubitVector<data_t>::apply_mcphase(const reg_t &qubits, const std::complex<double> phase) {
  const size_t N = qubits.size();
  switch (N) {
    case 1: {
      // Lambda function for arbitrary Phase gate with diagonal [1, phase]
      auto lambda = [&](const areg_t<2> &inds)->void {
        data_[inds[1]] *= phase;
      };
      apply_lambda(lambda, areg_t<1>({{qubits[0]}}));
      return;
    }
    case 2: {
      // Lambda function for CPhase gate with diagonal [1, 1, 1, phase]
      auto lambda = [&](const areg_t<4> &inds)->void {
        data_[inds[3]] *= phase;
      };
      apply_lambda(lambda, areg_t<2>({{qubits[0], qubits[1]}}));
      return;
    }
    case 3: {
      auto lambda = [&](const areg_t<8> &inds)->void {
         data_[inds[7]] *= phase;
      };
      apply_lambda(lambda, areg_t<3>({{qubits[0], qubits[1], qubits[2]}}));
      return;
    }
    default: {
      // Lambda function for general multi-controlled Phase gate
      // with diagonal [1, ..., 1, phase]
      auto lambda = [&](const indexes_t &inds)->void {
         data_[inds[MASKS[N]]] *= phase;
      };
      apply_lambda(lambda, qubits);
    }
  } // end switch
}

template <typename data_t>
void QubitVector<data_t>::apply_mcu(const reg_t &qubits,
                                    const cvector_t<double> &mat) {

  // Calculate the permutation positions for the last qubit.
  const size_t N = qubits.size();
  const size_t pos0 = MASKS[N - 1];
  const size_t pos1 = MASKS[N];

  // Check if matrix is actually diagonal and if so use 
  // diagonal matrix lambda function
  // TODO: this should be changed to not check doubles with ==
  if (mat[1] == 0.0 && mat[2] == 0.0) {
    // Check if actually a phase gate
    if (mat[0] == 1.0) {
      apply_mcphase(qubits, mat[3]);
      return;
    }
    // Otherwise apply general diagonal gate
    const cvector_t<double> diag = {{mat[0], mat[3]}};
    // Diagonal version
    switch (N) {
      case 1: {
        // If N=1 this is just a single-qubit matrix
        apply_diagonal_matrix(qubits[0], diag);
        return;
      }
      case 2: {
        // Lambda function for CU gate
        auto lambda = [&](const areg_t<4> &inds,
                          const cvector_t<data_t> &_diag)->void {
          data_[inds[pos0]] = _diag[0] * data_[inds[pos0]];
          data_[inds[pos1]] = _diag[1] * data_[inds[pos1]];
        };
        apply_lambda(lambda, areg_t<2>({{qubits[0], qubits[1]}}), convert(diag));
        return;
      }
      case 3: {
        // Lambda function for CCU gate
        auto lambda = [&](const areg_t<8> &inds,
                          const cvector_t<data_t> &_diag)->void {
          data_[inds[pos0]] = _diag[0] * data_[inds[pos0]];
          data_[inds[pos1]] = _diag[1] * data_[inds[pos1]];
        };
        apply_lambda(lambda, areg_t<3>({{qubits[0], qubits[1], qubits[2]}}), convert(diag));
        return;
      }
      default: {
        // Lambda function for general multi-controlled U gate
        auto lambda = [&](const indexes_t &inds,
                          const cvector_t<data_t> &_diag)->void {
          data_[inds[pos0]] = _diag[0] * data_[inds[pos0]];
          data_[inds[pos1]] = _diag[1] * data_[inds[pos1]];
        };
        apply_lambda(lambda, qubits, convert(diag));
        return;
      }
    } // end switch
  }

  // Non-diagonal version
  switch (N) {
    case 1: {
      // If N=1 this is just a single-qubit matrix
      apply_matrix(qubits[0], mat);
      return;
    }
    case 2: {
      // Lambda function for CU gate
      auto lambda = [&](const areg_t<4> &inds,
                        const cvector_t<data_t> &_mat)->void {
      const auto cache = data_[inds[pos0]];
      data_[inds[pos0]] = _mat[0] * data_[inds[pos0]] + _mat[2] * data_[inds[pos1]];
      data_[inds[pos1]] = _mat[1] * cache + _mat[3] * data_[inds[pos1]];
      };
      apply_lambda(lambda, areg_t<2>({{qubits[0], qubits[1]}}), convert(mat));
      return;
    }
    case 3: {
      // Lambda function for CCU gate
      auto lambda = [&](const areg_t<8> &inds,
                        const cvector_t<data_t> &_mat)->void {
      const auto cache = data_[inds[pos0]];
      data_[inds[pos0]] = _mat[0] * data_[inds[pos0]] + _mat[2] * data_[inds[pos1]];
      data_[inds[pos1]] = _mat[1] * cache + _mat[3] * data_[inds[pos1]];
      };
      apply_lambda(lambda, areg_t<3>({{qubits[0], qubits[1], qubits[2]}}), convert(mat));
      return;
    }
    default: {
      // Lambda function for general multi-controlled U gate
      auto lambda = [&](const indexes_t &inds,
                        const cvector_t<data_t> &_mat)->void {
      const auto cache = data_[inds[pos0]];
      data_[inds[pos0]] = _mat[0] * data_[inds[pos0]] + _mat[2] * data_[inds[pos1]];
      data_[inds[pos1]] = _mat[1] * cache + _mat[3] * data_[inds[pos1]];
      };
      apply_lambda(lambda, qubits, convert(mat));
      return;
    }
  } // end switch
}

//------------------------------------------------------------------------------
// Single-qubit matrices
//------------------------------------------------------------------------------

template <typename data_t>
void QubitVector<data_t>::apply_matrix(const uint_t qubit,
                                       const cvector_t<double>& mat) {

#ifdef SIMD_PPC64LE
  apply_matrix_ppc64le(qubit, mat);
#else 
  // Check if matrix is diagonal and if so use optimized lambda
  if (mat[1] == 0.0 && mat[2] == 0.0) {
    const cvector_t<double> diag = {{mat[0], mat[3]}};
    apply_diagonal_matrix(qubit, diag);
    return;
  }
  
  // Convert qubit to array register for lambda functions
  areg_t<1> qubits = {{qubit}};

  // Check if anti-diagonal matrix and if so use optimized lambda
  if(mat[0] == 0.0 && mat[3] == 0.0) {
    if (mat[1] == 1.0 && mat[2] == 1.0) {
      // X-matrix
      auto lambda = [&](const areg_t<2> &inds)->void {
        std::swap(data_[inds[0]], data_[inds[1]]);
      };
      apply_lambda(lambda, qubits);
      return;
    }
    if (mat[2] == 0.0) {
      // Non-unitary projector
      // possibly used in measure/reset/kraus update
      auto lambda = [&](const areg_t<2> &inds,
                        const cvector_t<data_t> &_mat)->void {
        data_[inds[1]] = _mat[1] * data_[inds[0]];
        data_[inds[0]] = 0.0;
      };
      apply_lambda(lambda, qubits, convert(mat));
      return;
    }
    if (mat[1] == 0.0) {
      // Non-unitary projector
      // possibly used in measure/reset/kraus update
      auto lambda = [&](const areg_t<2> &inds,
                        const cvector_t<data_t> &_mat)->void {
        data_[inds[0]] = _mat[2] * data_[inds[1]];
        data_[inds[1]] = 0.0;
      };
      apply_lambda(lambda, qubits, convert(mat));
      return;
    }
    // else we have a general anti-diagonal matrix
    auto lambda = [&](const areg_t<2> &inds,
                      const cvector_t<data_t> &_mat)->void {
      const std::complex<data_t> cache = data_[inds[0]];
      data_[inds[0]] = _mat[2] * data_[inds[1]];
      data_[inds[1]] = _mat[1] * cache;
    };
    apply_lambda(lambda, qubits, convert(mat));
    return;
  }
  // Otherwise general single-qubit matrix multiplication
  auto lambda = [&](const areg_t<2> &inds, const cvector_t<data_t> &_mat)->void {
    const auto cache = data_[inds[0]];
    data_[inds[0]] = _mat[0] * cache + _mat[2] * data_[inds[1]];
    data_[inds[1]] = _mat[1] * cache + _mat[3] * data_[inds[1]];
  };
  apply_lambda(lambda, qubits, convert(mat));
#endif
}

template <typename data_t>
void QubitVector<data_t>::apply_diagonal_matrix(const uint_t qubit,
                                                const cvector_t<double>& diag) {

  // TODO: This should be changed so it isn't checking doubles with ==
  if (diag[0] == 1.0) {  // [[1, 0], [0, z]] matrix
    if (diag[1] == 1.0)
      return; // Identity

    if (diag[1] == std::complex<double>(0., -1.)) { // [[1, 0], [0, -i]]
      auto lambda = [&](const areg_t<2> &inds,
                        const cvector_t<data_t> &_mat)->void {
        const auto k = inds[1];
        double cache = data_[k].imag();
        data_[k].imag(data_[k].real() * -1.);
        data_[k].real(cache);
      };
      apply_lambda(lambda, areg_t<1>({{qubit}}), convert(diag));
      return;
    }
    if (diag[1] == std::complex<double>(0., 1.)) {
      // [[1, 0], [0, i]]
      auto lambda = [&](const areg_t<2> &inds,
                        const cvector_t<data_t> &_mat)->void {
        const auto k = inds[1];
        double cache = data_[k].imag();
        data_[k].imag(data_[k].real());
        data_[k].real(cache * -1.);
      };
      apply_lambda(lambda, areg_t<1>({{qubit}}), convert(diag));
      return;
    } 
    if (diag[0] == 0.0) {
      // [[1, 0], [0, 0]]
      auto lambda = [&](const areg_t<2> &inds,
                        const cvector_t<data_t> &_mat)->void {
        data_[inds[1]] = 0.0;
      };
      apply_lambda(lambda, areg_t<1>({{qubit}}), convert(diag));
      return;
    } 
    // general [[1, 0], [0, z]]
    auto lambda = [&](const areg_t<2> &inds,
                      const cvector_t<data_t> &_mat)->void {
      const auto k = inds[1];
      data_[k] *= _mat[1];
    };
    apply_lambda(lambda, areg_t<1>({{qubit}}), convert(diag));
    return;
  } else if (diag[1] == 1.0) {
    // [[z, 0], [0, 1]] matrix
    if (diag[0] == std::complex<double>(0., -1.)) {
      // [[-i, 0], [0, 1]]
      auto lambda = [&](const areg_t<2> &inds,
                        const cvector_t<data_t> &_mat)->void {
        const auto k = inds[1];
        double cache = data_[k].imag();
        data_[k].imag(data_[k].real() * -1.);
        data_[k].real(cache);
      };
      apply_lambda(lambda, areg_t<1>({{qubit}}), convert(diag));
      return;
    } 
    if (diag[0] == std::complex<double>(0., 1.)) {
      // [[i, 0], [0, 1]]
      auto lambda = [&](const areg_t<2> &inds,
                        const cvector_t<data_t> &_mat)->void {
        const auto k = inds[1];
        double cache = data_[k].imag();
        data_[k].imag(data_[k].real());
        data_[k].real(cache * -1.);
      };
      apply_lambda(lambda, areg_t<1>({{qubit}}), convert(diag));
      return;
    } 
    if (diag[0] == 0.0) {
      // [[0, 0], [0, 1]]
      auto lambda = [&](const areg_t<2> &inds,
                        const cvector_t<data_t> &_mat)->void {
        data_[inds[0]] = 0.0;
      };
      apply_lambda(lambda, areg_t<1>({{qubit}}), convert(diag));
      return;
    } 
    // general [[z, 0], [0, 1]]
    auto lambda = [&](const areg_t<2> &inds,
                      const cvector_t<data_t> &_mat)->void {
      const auto k = inds[0];
      data_[k] *= _mat[0];
    };
    apply_lambda(lambda, areg_t<1>({{qubit}}), convert(diag));
    return;
  } else {
    // Lambda function for diagonal matrix multiplication
    auto lambda = [&](const areg_t<2> &inds,
                      const cvector_t<data_t> &_mat)->void {
      const auto k0 = inds[0];
      const auto k1 = inds[1];
      data_[k0] *= _mat[0];
      data_[k1] *= _mat[1];
    };
    apply_lambda(lambda, areg_t<1>({{qubit}}), convert(diag));
  }
}

/*******************************************************************************
 *
 * NORMS
 *
 ******************************************************************************/
template <typename data_t>
double QubitVector<data_t>::norm() const {
  // Lambda function for norm
  auto lambda = [&](int_t k, double &val_re, double &val_im)->void {
    (void)val_im; // unused
    val_re += std::real(data_[k] * std::conj(data_[k]));
  };
  return std::real(apply_reduction_lambda(lambda));
}

template <typename data_t>
double QubitVector<data_t>::norm(const reg_t &qubits, const cvector_t<double> &mat) const {

  const uint_t N = qubits.size();

  // Error checking
  #ifdef DEBUG
  check_vector(mat, 2 * N);
  #endif

  // Static array optimized lambda functions
  switch (N) {
    case 1:
      return norm(qubits[0], mat);
    case 2: {
      // Lambda function for 2-qubit matrix norm
      auto lambda = [&](const areg_t<4> &inds, const cvector_t<data_t> &_mat,
                        double &val_re, double &val_im)->void {
        (void)val_im; // unused
        for (size_t i = 0; i < 4; i++) {
          std::complex<data_t> vi = 0;
          for (size_t j = 0; j < 4; j++)
            vi += _mat[i + 4 * j] * data_[inds[j]];
          val_re += std::real(vi * std::conj(vi));
        }
      };
      areg_t<2> qubits_arr = {{qubits[0], qubits[1]}};
      return std::real(apply_reduction_lambda(lambda, qubits_arr, convert(mat)));
    }
    case 3: {
      // Lambda function for 3-qubit matrix norm
      auto lambda = [&](const areg_t<8> &inds, const cvector_t<data_t> &_mat,
                        double &val_re, double &val_im)->void {
        (void)val_im; // unused
        for (size_t i = 0; i < 8; i++) {
          std::complex<data_t> vi = 0;
          for (size_t j = 0; j < 8; j++)
            vi += _mat[i + 8 * j] * data_[inds[j]];
          val_re += std::real(vi * std::conj(vi));
        }
      };
      areg_t<3> qubits_arr = {{qubits[0], qubits[1], qubits[2]}};
      return std::real(apply_reduction_lambda(lambda, qubits_arr, convert(mat)));
    }
    case 4: {
      // Lambda function for 4-qubit matrix norm
      auto lambda = [&](const areg_t<16> &inds, const cvector_t<data_t> &_mat,
                        double &val_re, double &val_im)->void {
        (void)val_im; // unused
        for (size_t i = 0; i < 16; i++) {
          std::complex<data_t> vi = 0;
          for (size_t j = 0; j < 16; j++)
            vi += _mat[i + 16 * j] * data_[inds[j]];
          val_re += std::real(vi * std::conj(vi));
        }
      };
      areg_t<4> qubits_arr = {{qubits[0], qubits[1], qubits[2], qubits[3]}};
      return std::real(apply_reduction_lambda(lambda, qubits_arr, convert(mat)));
    }
    default: {
      // Lambda function for N-qubit matrix norm
      const uint_t DIM = BITS[N];
      auto lambda = [&](const indexes_t &inds, const cvector_t<data_t> &_mat,
                        double &val_re, double &val_im)->void {
        (void)val_im; // unused
        for (size_t i = 0; i < DIM; i++) {
          std::complex<data_t> vi = 0;
          for (size_t j = 0; j < DIM; j++)
            vi += _mat[i + DIM * j] * data_[inds[j]];
          val_re += std::real(vi * std::conj(vi));
        }
      };
      // Use the lambda function
      return std::real(apply_reduction_lambda(lambda, qubits, convert(mat)));
    }
  } // end switch
}

template <typename data_t>
double QubitVector<data_t>::norm_diagonal(const reg_t &qubits, const cvector_t<double> &mat) const {

  const uint_t N = qubits.size();

  // Error checking
  #ifdef DEBUG
  check_vector(mat, N);
  #endif

  // Static array optimized lambda functions
  switch (N) {
    case 1:
      return norm_diagonal(qubits[0], mat);
    case 2: {
      // Lambda function for 2-qubit matrix norm
      auto lambda = [&](const areg_t<4> &inds, const cvector_t<data_t> &_mat,
                        double &val_re, double &val_im)->void {
        (void)val_im; // unused
        for (size_t i = 0; i < 4; i++) {
          const auto vi = _mat[i] * data_[inds[i]];
          val_re += std::real(vi * std::conj(vi));
        }
      };
      areg_t<2> qubits_arr = {{qubits[0], qubits[1]}};
      return std::real(apply_reduction_lambda(lambda, qubits_arr, convert(mat)));
    }
    case 3: {
      // Lambda function for 3-qubit matrix norm
      auto lambda = [&](const areg_t<8> &inds, const cvector_t<data_t> &_mat,
                        double &val_re, double &val_im)->void {
        (void)val_im; // unused
        for (size_t i = 0; i < 8; i++) {
          const auto vi = _mat[i] * data_[inds[i]];
          val_re += std::real(vi * std::conj(vi));
        }
      };
      areg_t<3> qubits_arr = {{qubits[0], qubits[1], qubits[2]}};
      return std::real(apply_reduction_lambda(lambda, qubits_arr, convert(mat)));
    }
    case 4: {
      // Lambda function for 4-qubit matrix norm
      auto lambda = [&](const areg_t<16> &inds, const cvector_t<data_t> &_mat,
                        double &val_re, double &val_im)->void {
        (void)val_im; // unused
        for (size_t i = 0; i < 16; i++) {
          const auto vi = _mat[i] * data_[inds[i]];
          val_re += std::real(vi * std::conj(vi));
        }
      };
      areg_t<4> qubits_arr = {{qubits[0], qubits[1], qubits[2], qubits[3]}};
      return std::real(apply_reduction_lambda(lambda, qubits_arr, convert(mat)));
    }
    default: {
      // Lambda function for N-qubit matrix norm
      const uint_t DIM = BITS[N];
      auto lambda = [&](const indexes_t &inds, const cvector_t<data_t> &_mat,
                        double &val_re, double &val_im)->void {
        (void)val_im; // unused
        for (size_t i = 0; i < DIM; i++) {
          const auto vi = _mat[i] * data_[inds[i]];
          val_re += std::real(vi * std::conj(vi));
        }
      };
      // Use the lambda function
      return std::real(apply_reduction_lambda(lambda, qubits, convert(mat)));
    }
  } // end switch
}

//------------------------------------------------------------------------------
// Single-qubit specialization
//------------------------------------------------------------------------------
template <typename data_t>
double QubitVector<data_t>::norm(const uint_t qubit, const cvector_t<double> &mat) const {
  // Error handling
  #ifdef DEBUG
  check_vector(mat, 2);
  #endif

  // Check if input matrix is diagonal, and if so use diagonal function.
  if (mat[1] == 0.0 && mat[2] == 0.0) {
    const cvector_t<double> diag = {{mat[0], mat[3]}};
    return norm_diagonal(qubit, diag);
  }

  // Lambda function for norm reduction to real value.
  auto lambda = [&](const areg_t<2> &inds,
                    const cvector_t<data_t> &_mat,
                    double &val_re,
                    double &val_im)->void {
    (void)val_im; // unused
    const auto v0 = _mat[0] * data_[inds[0]] + _mat[2] * data_[inds[1]];
    const auto v1 = _mat[1] * data_[inds[0]] + _mat[3] * data_[inds[1]];
    val_re += std::real(v0 * std::conj(v0)) + std::real(v1 * std::conj(v1));
  };
  return std::real(apply_reduction_lambda(lambda, areg_t<1>({{qubit}}), convert(mat)));
}

template <typename data_t>
double QubitVector<data_t>::norm_diagonal(const uint_t qubit, const cvector_t<double> &mat) const {
  // Error handling
  #ifdef DEBUG
  check_vector(mat, 1);
  #endif
  // Lambda function for norm reduction to real value.
  auto lambda = [&](const areg_t<2> &inds,
                    const cvector_t<data_t> &_mat,
                    double &val_re,
                    double &val_im)->void {
    (void)val_im; // unused
    const auto v0 = _mat[0] * data_[inds[0]];
    const auto v1 = _mat[1] * data_[inds[1]];
    val_re += std::real(v0 * std::conj(v0)) + std::real(v1 * std::conj(v1));
  };
  return std::real(apply_reduction_lambda(lambda, areg_t<1>({{qubit}}), convert(mat)));
}


/*******************************************************************************
 *
 * Probabilities
 *
 ******************************************************************************/
template <typename data_t>
double QubitVector<data_t>::probability(const uint_t outcome) const {
  return std::real(data_[outcome] * std::conj(data_[outcome]));
}

template <typename data_t>
std::vector<double> QubitVector<data_t>::probabilities() const {
  const int_t END = 1LL << num_qubits();
  std::vector<double> probs(END, 0.);
#pragma omp parallel for if (num_qubits_ > omp_threshold_ && omp_threads_ > 1) num_threads(omp_threads_)
  for (int_t j=0; j < END; j++) {
    probs[j] = probability(j);
  }
  return probs;
}

template <typename data_t>
std::vector<double> QubitVector<data_t>::probabilities(const reg_t &qubits) const {

  const size_t N = qubits.size();
  const int_t DIM = BITS[N];
  const int_t END = BITS[num_qubits() - N];

  // Error checking
  #ifdef DEBUG
  for (const auto &qubit : qubits)
    check_qubit(qubit);
  #endif

  auto qubits_sorted = qubits;
  std::sort(qubits_sorted.begin(), qubits_sorted.end());
  if ((N == num_qubits_) && (qubits == qubits_sorted))
    return probabilities();

  std::vector<double> probs(DIM, 0.);
  #pragma omp parallel if (num_qubits_ > omp_threshold_ && omp_threads_ > 1) num_threads(omp_threads_)
  {
    std::vector<data_t> probs_private(DIM, 0.);
    #pragma omp for
      for (int_t k = 0; k < END; k++) {
        auto idx = indexes(qubits, qubits_sorted, k);
        for (int_t m = 0; m < DIM; ++m) {
          probs_private[m] += probability(idx[m]);
        }
      }
    #pragma omp critical
    for (int_t m = 0; m < DIM; ++m) {
      probs[m] += probs_private[m];
    }
  }
  
  return probs;
}

//------------------------------------------------------------------------------
// Sample measure outcomes
//------------------------------------------------------------------------------
template <typename data_t>
reg_t QubitVector<data_t>::sample_measure(const std::vector<double> &rnds) const {

  const int_t END = 1LL << num_qubits();
  const int_t SHOTS = rnds.size();
  reg_t samples;
  samples.assign(SHOTS, 0);

  const int INDEX_SIZE = sample_measure_index_size_;
  const int_t INDEX_END = BITS[INDEX_SIZE];
  // Qubit number is below index size, loop over shots
  if (END < INDEX_END) {
    #pragma omp parallel if (num_qubits_ > omp_threshold_ && omp_threads_ > 1) num_threads(omp_threads_)
    {
      #pragma omp for
      for (int_t i = 0; i < SHOTS; ++i) {
        double rnd = rnds[i];
        double p = .0;
        int_t sample;
        for (sample = 0; sample < END - 1; ++sample) {
          p += probability(sample);
          if (rnd < p)
            break;
        }
        samples[i] = sample;
      }
    } // end omp parallel
  }
  // Qubit number is above index size, loop over index blocks
  else {
    // Initialize indexes
    std::vector<double> idxs;
    idxs.assign(INDEX_END, 0.0);
    uint_t loop = (END >> INDEX_SIZE);
    #pragma omp parallel if (num_qubits_ > omp_threshold_ && omp_threads_ > 1) num_threads(omp_threads_)
    {
      #pragma omp for
      for (int_t i = 0; i < INDEX_END; ++i) {
        uint_t base = loop * i;
        double total = .0;
        double p = .0;
        for (uint_t j = 0; j < loop; ++j) {
          uint_t k = base | j;
          p = probability(k);
          total += p;
        }
        idxs[i] = total;
      }
    } // end omp parallel

    #pragma omp parallel if (num_qubits_ > omp_threshold_ && omp_threads_ > 1) num_threads(omp_threads_)
    {
      #pragma omp for
      for (int_t i = 0; i < SHOTS; ++i) {
        double rnd = rnds[i];
        double p = .0;
        int_t sample = 0;
        for (uint_t j = 0; j < idxs.size(); ++j) {
          if (rnd < (p + idxs[j])) {
            break;
          }
          p += idxs[j];
          sample += loop;
        }

        for (; sample < END - 1; ++sample) {
          p += probability(sample);
          if (rnd < p){
            break;
          }
        }
        samples[i] = sample;
      }
    } // end omp parallel
  }
  return samples;
}

#ifdef SIMD_PPC64LE
//------------------------------------------------------------------------------
// Single-qubit matrices
//------------------------------------------------------------------------------

const __vector unsigned char vpat = {
                0x08, 0x09, 0x0A, 0x0B,
                0x0C, 0x0D, 0x0E, 0x0F,
                0x00, 0x01, 0x02, 0x03,
                0x04, 0x05, 0x06, 0x07};

inline __vector double complex_mul(__vector double vec, __vector double imag, __vector double real){
        __vector double perm = vec_perm(vec, vec, vpat);
        __vector double temp1 = vec_mul(real, vec);
        __vector double temp2 = vec_mul(imag, perm);
        return vec_sub(temp1, temp2);
}

template <typename data_t>
void QubitVector<data_t>::apply_matrix_ppc64le(const uint_t qubit,
                                       const cvector_t<double>& mat) {
  // Check if matrix is diagonal and if so use optimized lambda
  if (mat[1] == 0.0 && mat[2] == 0.0) {
    const cvector_t<double> diag = {{mat[0], mat[3]}};
    apply_diagonal_matrix(qubit, diag);
    return;
  }
  
  // Convert qubit to array register for lambda functions
  areg_t<1> qubits = {{qubit}};

  // Check if anti-diagonal matrix and if so use optimized lambda
  if(mat[0] == 0.0 && mat[3] == 0.0) {
    if (mat[1] == 1.0 && mat[2] == 1.0) {
      // X-matrix
      auto lambda = [&](const areg_t<2> &inds)->void {
        std::swap(data_[inds[0]], data_[inds[1]]);
      };
      apply_lambda(lambda, qubits);
      return;
    }
    if (mat[2] == 0.0) {
      // Non-unitary projector
      // possibly used in measure/reset/kraus update
      auto lambda = [&](const areg_t<2> &inds,
                        const cvector_t<data_t> &_mat)->void {
        data_[inds[1]] = _mat[1] * data_[inds[0]];
        data_[inds[0]] = 0.0;
      };
      apply_lambda(lambda, qubits, convert(mat));
      return;
    }
    if (mat[1] == 0.0) {
      // Non-unitary projector
      // possibly used in measure/reset/kraus update
      auto lambda = [&](const areg_t<2> &inds,
                        const cvector_t<data_t> &_mat)->void {
        data_[inds[0]] = _mat[2] * data_[inds[1]];
        data_[inds[1]] = 0.0;
      };
      apply_lambda(lambda, qubits, convert(mat));
      return;
    }
    // else we have a general anti-diagonal matrix
    auto lambda = [&](const areg_t<2> &inds,
                      const cvector_t<data_t> &_mat)->void {
      const std::complex<data_t> cache = data_[inds[0]];
      data_[inds[0]] = _mat[2] * data_[inds[1]];
      data_[inds[1]] = _mat[1] * cache;
    };
    apply_lambda(lambda, qubits, convert(mat));
    return;
  }
  if(type_double) {
    __vector double real[4];
    __vector double imag[4];

    int mat_vec_index = 0;
    for(int i = 0; i < 2; i++){
      for(int j = 0; j < 2; j++){ 
        real[mat_vec_index] = (__vector double){mat[i+j*2].real(), mat[i+j*2].real()}; 
	imag[mat_vec_index] = (__vector double){mat[i+j*2].imag(), mat[i+j*2].imag() * -1}; 
	mat_vec_index = mat_vec_index+1;
      }
    }

    auto lambda = [&](const areg_t<2> &inds, const cvector_t<data_t> &_mat)->void {
      __vector double vec[2];
      vec[0] = (__vector double){data_[inds[0]].real(), data_[inds[0]].imag()};
      vec[1] = (__vector double){data_[inds[1]].real(), data_[inds[1]].imag()};

      __vector double result = vec_add(complex_mul(vec[0],imag[0],real[0]), complex_mul(vec[1],imag[1],real[1]));
      data_[inds[0]] = std::complex<double>((double)result[0], (double)result[1]);
      result = vec_add(complex_mul(vec[0],imag[2],real[2]), complex_mul(vec[1],imag[3],real[3]));
      data_[inds[1]] = std::complex<double>((double)result[0], (double)result[1]);
    };
    apply_lambda(lambda, qubits, convert(mat));
  } else {
    // Otherwise general single-qubit matrix multiplication
    auto lambda = [&](const areg_t<2> &inds, const cvector_t<data_t> &_mat)->void {
      const auto cache = data_[inds[0]];
      data_[inds[0]] = _mat[0] * cache + _mat[2] * data_[inds[1]];
      data_[inds[1]] = _mat[1] * cache + _mat[3] * data_[inds[1]];
    };
    apply_lambda(lambda, qubits, convert(mat));
 }
}

/*******************************************************************************
 *
 * MATRIX MULTIPLICATION
 *
 ******************************************************************************/
template <typename data_t>
void QubitVector<data_t>::apply_matrix_ppc64le(const reg_t &qubits,
                                       const cvector_t<double> &mat) {

  const size_t N = qubits.size();
  // Error checking
  #ifdef DEBUG
  check_vector(mat, 2 * N);
  #endif

  // Static array optimized lambda functions
  switch (N) {
    case 1:
      apply_matrix(qubits[0], mat);
      return;
    case 2: {
      if(type_double) {
        __vector double real[16];
        __vector double imag[16]; 
	int mat_vec_index = 0;

        for(int i = 0; i < 4; i++){
          for(int j = 0; j < 4; j++){
	    real[mat_vec_index] = (__vector double){mat[i+j*4].real(), mat[i+j*4].real()};
	    imag[mat_vec_index] = (__vector double){mat[i+j*4].imag(), mat[i+j*4].imag() * -1};
            mat_vec_index = mat_vec_index+1;
	  }
        }

        auto lambda = [&](const areg_t<4> &inds, const cvector_t<data_t> &_mat)->void {
          __vector double vec[4];

          vec[0] = (__vector double){data_[inds[0]].real(), data_[inds[0]].imag()};
          vec[1] = (__vector double){data_[inds[1]].real(), data_[inds[1]].imag()};
          vec[2] = (__vector double){data_[inds[2]].real(), data_[inds[2]].imag()};
          vec[3] = (__vector double){data_[inds[3]].real(), data_[inds[3]].imag()};

          __vector double result = vec_add(vec_add(complex_mul(vec[0],imag[0],real[0]), complex_mul(vec[1],imag[1],real[1])),vec_add(complex_mul(vec[2],imag[2],real[2]), complex_mul(vec[3],imag[3],real[3])));
          data_[inds[0]] = std::complex<double>((double)result[0], (double)result[1]);
          result = vec_add(vec_add(complex_mul(vec[0],imag[4],real[4]), complex_mul(vec[1],imag[5],real[5])),vec_add(complex_mul(vec[2],imag[6],real[6]), complex_mul(vec[3],imag[7],real[7])));
          data_[inds[1]] = std::complex<double>((double)result[0], (double)result[1]);
          result = vec_add(vec_add(complex_mul(vec[0],imag[8],real[8]), complex_mul(vec[1],imag[9],real[9])),vec_add(complex_mul(vec[2],imag[10],real[10]), complex_mul(vec[3],imag[11],real[11])));
          data_[inds[2]] = std::complex<double>((double)result[0], (double)result[1]);
          result = vec_add(vec_add(complex_mul(vec[0],imag[12],real[12]), complex_mul(vec[1],imag[13],real[13])),vec_add(complex_mul(vec[2],imag[14],real[14]), complex_mul(vec[3],imag[15],real[15])));
          data_[inds[3]] = std::complex<double>((double)result[0], (double)result[1]);
        };
        apply_lambda(lambda, areg_t<2>({{qubits[0], qubits[1]}}), convert(mat));
	return;
      } else {
      // Lambda function for 2-qubit matrix multiplication
        auto lambda = [&](const areg_t<4> &inds, const cvector_t<data_t> &_mat)->void {
          std::array<std::complex<data_t>, 4> cache;
          for (size_t i = 0; i < 4; i++) {
            const auto ii = inds[i];
            cache[i] = data_[ii];
            data_[ii] = 0.;
          }
        // update state vector
          for (size_t i = 0; i < 4; i++)
            for (size_t j = 0; j < 4; j++)
              data_[inds[i]] += _mat[i + 4 * j] * cache[j];
        };
        apply_lambda(lambda, areg_t<2>({{qubits[0], qubits[1]}}), convert(mat));
      return;
      }
    }
    case 3: {
      if(type_double) {
        __vector double real[64];
        __vector double imag[64]; 
	int mat_vec_index = 0;
        for(int i = 0; i < 8; i++){
          for(int j = 0; j < 8; j++){
	    real[mat_vec_index] = (__vector double){mat[i+j*8].real(), mat[i+j*8].real()};
	    imag[mat_vec_index] = (__vector double){mat[i+j*8].imag(), mat[i+j*8].imag() * -1};
            mat_vec_index = mat_vec_index+1;
	  }
        }
        // Lambda function for 3-qubit matrix multiplication
        auto lambda = [&](const areg_t<8> &inds, const cvector_t<data_t> &_mat)->void {
          __vector double vec[4];
          __vector double temp[8];
          __vector double result;

          vec[0] = (__vector double){data_[inds[0]].real(), data_[inds[0]].imag()};
          vec[1] = (__vector double){data_[inds[1]].real(), data_[inds[1]].imag()};
          vec[2] = (__vector double){data_[inds[2]].real(), data_[inds[2]].imag()};
          vec[3] = (__vector double){data_[inds[3]].real(), data_[inds[3]].imag()};

          temp[0] = vec_add(vec_add(complex_mul(vec[0],imag[0],real[0]), complex_mul(vec[1],imag[1],real[1])),vec_add(complex_mul(vec[2],imag[2],real[2]), complex_mul(vec[3],imag[3],real[3])));
          temp[1] = vec_add(vec_add(complex_mul(vec[0],imag[8],real[8]), complex_mul(vec[1],imag[9],real[9])),vec_add(complex_mul(vec[2],imag[10],real[10]), complex_mul(vec[3],imag[11],real[11])));
          temp[2] = vec_add(vec_add(complex_mul(vec[0],imag[16],real[16]), complex_mul(vec[1],imag[17],real[17])),vec_add(complex_mul(vec[2],imag[18],real[18]), complex_mul(vec[3],imag[19],real[19])));
          temp[3] = vec_add(vec_add(complex_mul(vec[0],imag[24],real[24]), complex_mul(vec[1],imag[25],real[25])),vec_add(complex_mul(vec[2],imag[26],real[26]), complex_mul(vec[3],imag[27],real[27])));
          temp[4] = vec_add(vec_add(complex_mul(vec[0],imag[32],real[32]), complex_mul(vec[1],imag[33],real[33])),vec_add(complex_mul(vec[2],imag[34],real[34]), complex_mul(vec[3],imag[35],real[35])));
          temp[5] = vec_add(vec_add(complex_mul(vec[0],imag[40],real[40]), complex_mul(vec[1],imag[41],real[41])),vec_add(complex_mul(vec[2],imag[42],real[42]), complex_mul(vec[3],imag[43],real[43])));
          temp[6] = vec_add(vec_add(complex_mul(vec[0],imag[48],real[48]), complex_mul(vec[1],imag[49],real[49])),vec_add(complex_mul(vec[2],imag[50],real[50]), complex_mul(vec[3],imag[51],real[51])));
          temp[7] = vec_add(vec_add(complex_mul(vec[0],imag[56],real[56]), complex_mul(vec[1],imag[57],real[57])),vec_add(complex_mul(vec[2],imag[58],real[58]), complex_mul(vec[3],imag[59],real[59])));

          vec[0] = (__vector double){data_[inds[4]].real(), data_[inds[4]].imag()};
          vec[1] = (__vector double){data_[inds[5]].real(), data_[inds[5]].imag()};
          vec[2] = (__vector double){data_[inds[6]].real(), data_[inds[6]].imag()};
          vec[3] = (__vector double){data_[inds[7]].real(), data_[inds[7]].imag()};

          result = vec_add(temp[0], vec_add(vec_add(complex_mul(vec[0],imag[4],real[4]), complex_mul(vec[1],imag[5],real[5])),vec_add(complex_mul(vec[2],imag[6],real[6]), complex_mul(vec[3],imag[7],real[7]))));
          data_[inds[0]] = std::complex<double>((double)result[0], (double)result[1]);
          result = vec_add(temp[1], vec_add(vec_add(complex_mul(vec[0],imag[12],real[12]), complex_mul(vec[1],imag[13],real[13])),vec_add(complex_mul(vec[2],imag[14],real[14]), complex_mul(vec[3],imag[15],real[15]))));
          data_[inds[1]] = std::complex<double>((double)result[0], (double)result[1]);
          result = vec_add(temp[2], vec_add(vec_add(complex_mul(vec[0],imag[20],real[20]), complex_mul(vec[1],imag[21],real[21])),vec_add(complex_mul(vec[2],imag[22],real[22]), complex_mul(vec[3],imag[23],real[23]))));
          data_[inds[2]] = std::complex<double>((double)result[0], (double)result[1]);
          result = vec_add(temp[3], vec_add(vec_add(complex_mul(vec[0],imag[28],real[28]), complex_mul(vec[1],imag[29],real[29])),vec_add(complex_mul(vec[2],imag[30],real[30]), complex_mul(vec[3],imag[31],real[31]))));
          data_[inds[3]] = std::complex<double>((double)result[0], (double)result[1]);
          result = vec_add(temp[4], vec_add(vec_add(complex_mul(vec[0],imag[36],real[36]), complex_mul(vec[1],imag[37],real[37])),vec_add(complex_mul(vec[2],imag[38],real[38]), complex_mul(vec[3],imag[39],real[39]))));
          data_[inds[4]] = std::complex<double>((double)result[0], (double)result[1]);
          result = vec_add(temp[5], vec_add(vec_add(complex_mul(vec[0],imag[44],real[44]), complex_mul(vec[1],imag[45],real[45])),vec_add(complex_mul(vec[2],imag[46],real[46]), complex_mul(vec[3],imag[47],real[47]))));
          data_[inds[5]] = std::complex<double>((double)result[0], (double)result[1]);
          result = vec_add(temp[6], vec_add(vec_add(complex_mul(vec[0],imag[52],real[52]), complex_mul(vec[1],imag[53],real[53])),vec_add(complex_mul(vec[2],imag[54],real[54]), complex_mul(vec[3],imag[55],real[55]))));
          data_[inds[6]] = std::complex<double>((double)result[0], (double)result[1]);
          result = vec_add(temp[7], vec_add(vec_add(complex_mul(vec[0],imag[60],real[60]), complex_mul(vec[1],imag[61],real[61])),vec_add(complex_mul(vec[2],imag[62],real[62]), complex_mul(vec[3],imag[63],real[63]))));
          data_[inds[7]] = std::complex<double>((double)result[0], (double)result[1]);
        };
        apply_lambda(lambda, areg_t<3>({{qubits[0], qubits[1], qubits[2]}}), convert(mat));
        return;
      } else {
        // Lambda function for 3-qubit matrix multiplication
        auto lambda = [&](const areg_t<8> &inds, const cvector_t<data_t> &_mat)->void {
          std::array<std::complex<data_t>, 8> cache;
          for (size_t i = 0; i < 8; i++) {
            const auto ii = inds[i];
            cache[i] = data_[ii];
            data_[ii] = 0.;
          }
          // update state vector
          for (size_t i = 0; i < 8; i++)
            for (size_t j = 0; j < 8; j++)
              data_[inds[i]] += _mat[i + 8 * j] * cache[j];
        };
        apply_lambda(lambda, areg_t<3>({{qubits[0], qubits[1], qubits[2]}}), convert(mat));
        return;
      }
    }
    case 4: {
      if(type_double) {
        __vector double real[256];
        __vector double imag[256]; 
        int mat_vec_index = 0;
        for(int i = 0; i < 16; i++){
          for(int j = 0; j < 16; j++){
	    real[mat_vec_index] = (__vector double){mat[i+j*16].real(), mat[i+j*16].real()};
	    imag[mat_vec_index] = (__vector double){mat[i+j*16].imag(), mat[i+j*16].imag() * -1};
            mat_vec_index = mat_vec_index+1;
	  }
        }
	
      // Lambda function for 4-qubit matrix multiplication
      auto lambda = [&](const areg_t<16> &inds, const cvector_t<data_t> &_mat)->void {
        __vector double vec[4];
        __vector double temp[16];
        __vector double result;


        vec[0] = (__vector double){data_[inds[0]].real(), data_[inds[0]].imag()};
        vec[1] = (__vector double){data_[inds[1]].real(), data_[inds[1]].imag()};
        vec[2] = (__vector double){data_[inds[2]].real(), data_[inds[2]].imag()};
        vec[3] = (__vector double){data_[inds[3]].real(), data_[inds[3]].imag()};

        temp[0] = vec_add(vec_add(complex_mul(vec[0],imag[0],real[0]), complex_mul(vec[1],imag[1],real[1])),vec_add(complex_mul(vec[2],imag[2],real[2]), complex_mul(vec[3],imag[3],real[3])));
        temp[1] = vec_add(vec_add(complex_mul(vec[0],imag[16],real[16]), complex_mul(vec[1],imag[17],real[17])),vec_add(complex_mul(vec[2],imag[18],real[18]), complex_mul(vec[3],imag[19],real[19])));
        temp[2] = vec_add(vec_add(complex_mul(vec[0],imag[32],real[32]), complex_mul(vec[1],imag[33],real[33])),vec_add(complex_mul(vec[2],imag[34],real[34]), complex_mul(vec[3],imag[35],real[35])));
        temp[3] = vec_add(vec_add(complex_mul(vec[0],imag[48],real[48]), complex_mul(vec[1],imag[49],real[49])),vec_add(complex_mul(vec[2],imag[50],real[50]), complex_mul(vec[3],imag[51],real[51])));
        temp[4] = vec_add(vec_add(complex_mul(vec[0],imag[64],real[64]), complex_mul(vec[1],imag[65],real[65])),vec_add(complex_mul(vec[2],imag[66],real[66]), complex_mul(vec[3],imag[67],real[67])));
        temp[5] = vec_add(vec_add(complex_mul(vec[0],imag[80],real[80]), complex_mul(vec[1],imag[81],real[81])),vec_add(complex_mul(vec[2],imag[82],real[82]), complex_mul(vec[3],imag[83],real[83])));
        temp[6] = vec_add(vec_add(complex_mul(vec[0],imag[96],real[96]), complex_mul(vec[1],imag[97],real[97])),vec_add(complex_mul(vec[2],imag[98],real[98]), complex_mul(vec[3],imag[99],real[99])));
        temp[7] = vec_add(vec_add(complex_mul(vec[0],imag[112],real[112]), complex_mul(vec[1],imag[113],real[113])),vec_add(complex_mul(vec[2],imag[114],real[114]), complex_mul(vec[3],imag[115],real[115])));
        temp[8] = vec_add(vec_add(complex_mul(vec[0],imag[128],real[128]), complex_mul(vec[1],imag[129],real[129])),vec_add(complex_mul(vec[2],imag[130],real[130]), complex_mul(vec[3],imag[131],real[131])));
        temp[9] = vec_add(vec_add(complex_mul(vec[0],imag[144],real[144]), complex_mul(vec[1],imag[145],real[145])),vec_add(complex_mul(vec[2],imag[146],real[146]), complex_mul(vec[3],imag[147],real[147])));
        temp[10] = vec_add(vec_add(complex_mul(vec[0],imag[160],real[160]), complex_mul(vec[1],imag[161],real[161])),vec_add(complex_mul(vec[2],imag[162],real[162]), complex_mul(vec[3],imag[163],real[163])));
        temp[11] = vec_add(vec_add(complex_mul(vec[0],imag[176],real[176]), complex_mul(vec[1],imag[177],real[177])),vec_add(complex_mul(vec[2],imag[178],real[178]), complex_mul(vec[3],imag[179],real[179])));
        temp[12] = vec_add(vec_add(complex_mul(vec[0],imag[192],real[192]), complex_mul(vec[1],imag[193],real[193])),vec_add(complex_mul(vec[2],imag[194],real[194]), complex_mul(vec[3],imag[195],real[195])));
        temp[13] = vec_add(vec_add(complex_mul(vec[0],imag[208],real[208]), complex_mul(vec[1],imag[209],real[209])),vec_add(complex_mul(vec[2],imag[210],real[210]), complex_mul(vec[3],imag[211],real[211])));
        temp[14] = vec_add(vec_add(complex_mul(vec[0],imag[224],real[224]), complex_mul(vec[1],imag[225],real[225])),vec_add(complex_mul(vec[2],imag[226],real[226]), complex_mul(vec[3],imag[227],real[227])));
        temp[15] = vec_add(vec_add(complex_mul(vec[0],imag[240],real[240]), complex_mul(vec[1],imag[241],real[241])),vec_add(complex_mul(vec[2],imag[242],real[242]), complex_mul(vec[3],imag[243],real[243])));

        vec[0] = (__vector double){data_[inds[4]].real(), data_[inds[4]].imag()};
        vec[1] = (__vector double){data_[inds[5]].real(), data_[inds[5]].imag()};
        vec[2] = (__vector double){data_[inds[6]].real(), data_[inds[6]].imag()};
        vec[3] = (__vector double){data_[inds[7]].real(), data_[inds[7]].imag()};
    
        temp[0] = vec_add(temp[0], vec_add(vec_add(complex_mul(vec[0],imag[4],real[4]), complex_mul(vec[1],imag[5],real[5])),vec_add(complex_mul(vec[2],imag[6],real[6]), complex_mul(vec[3],imag[7],real[7]))));
        temp[1] = vec_add(temp[1], vec_add(vec_add(complex_mul(vec[0],imag[20],real[20]), complex_mul(vec[1],imag[21],real[21])),vec_add(complex_mul(vec[2],imag[22],real[22]), complex_mul(vec[3],imag[23],real[23]))));
        temp[2] = vec_add(temp[2], vec_add(vec_add(complex_mul(vec[0],imag[36],real[36]), complex_mul(vec[1],imag[37],real[37])),vec_add(complex_mul(vec[2],imag[38],real[38]), complex_mul(vec[3],imag[39],real[39]))));
        temp[3] = vec_add(temp[3], vec_add(vec_add(complex_mul(vec[0],imag[52],real[52]), complex_mul(vec[1],imag[53],real[53])),vec_add(complex_mul(vec[2],imag[54],real[54]), complex_mul(vec[3],imag[55],real[55]))));
        temp[4] = vec_add(temp[4], vec_add(vec_add(complex_mul(vec[0],imag[68],real[68]), complex_mul(vec[1],imag[69],real[69])),vec_add(complex_mul(vec[2],imag[70],real[70]), complex_mul(vec[3],imag[71],real[71]))));
        temp[5] = vec_add(temp[5], vec_add(vec_add(complex_mul(vec[0],imag[84],real[84]), complex_mul(vec[1],imag[85],real[85])),vec_add(complex_mul(vec[2],imag[86],real[86]), complex_mul(vec[3],imag[87],real[87]))));
        temp[6] = vec_add(temp[6], vec_add(vec_add(complex_mul(vec[0],imag[100],real[100]), complex_mul(vec[1],imag[101],real[101])),vec_add(complex_mul(vec[2],imag[102],real[102]), complex_mul(vec[3],imag[103],real[103]))));
        temp[7] = vec_add(temp[7], vec_add(vec_add(complex_mul(vec[0],imag[116],real[116]), complex_mul(vec[1],imag[117],real[117])),vec_add(complex_mul(vec[2],imag[118],real[118]), complex_mul(vec[3],imag[119],real[119]))));
        temp[8] = vec_add(temp[8], vec_add(vec_add(complex_mul(vec[0],imag[132],real[132]), complex_mul(vec[1],imag[133],real[133])),vec_add(complex_mul(vec[2],imag[134],real[134]), complex_mul(vec[3],imag[135],real[135]))));
        temp[9] = vec_add(temp[9], vec_add(vec_add(complex_mul(vec[0],imag[148],real[148]), complex_mul(vec[1],imag[149],real[149])),vec_add(complex_mul(vec[2],imag[150],real[150]), complex_mul(vec[3],imag[151],real[151]))));
        temp[10] = vec_add(temp[10], vec_add(vec_add(complex_mul(vec[0],imag[164],real[164]), complex_mul(vec[1],imag[165],real[165])),vec_add(complex_mul(vec[2],imag[166],real[166]), complex_mul(vec[3],imag[167],real[167]))));
        temp[11] = vec_add(temp[11], vec_add(vec_add(complex_mul(vec[0],imag[180],real[180]), complex_mul(vec[1],imag[181],real[181])),vec_add(complex_mul(vec[2],imag[182],real[182]), complex_mul(vec[3],imag[183],real[183]))));
        temp[12] = vec_add(temp[12], vec_add(vec_add(complex_mul(vec[0],imag[196],real[196]), complex_mul(vec[1],imag[197],real[197])),vec_add(complex_mul(vec[2],imag[198],real[198]), complex_mul(vec[3],imag[199],real[199]))));
        temp[13] = vec_add(temp[13], vec_add(vec_add(complex_mul(vec[0],imag[212],real[212]), complex_mul(vec[1],imag[213],real[213])),vec_add(complex_mul(vec[2],imag[214],real[214]), complex_mul(vec[3],imag[215],real[215]))));
        temp[14] = vec_add(temp[14], vec_add(vec_add(complex_mul(vec[0],imag[228],real[228]), complex_mul(vec[1],imag[229],real[229])),vec_add(complex_mul(vec[2],imag[230],real[230]), complex_mul(vec[3],imag[231],real[231]))));
        temp[15] = vec_add(temp[15], vec_add(vec_add(complex_mul(vec[0],imag[244],real[244]), complex_mul(vec[1],imag[245],real[245])),vec_add(complex_mul(vec[2],imag[246],real[246]), complex_mul(vec[3],imag[247],real[247]))));
    
        vec[0] = (__vector double){data_[inds[8]].real(), data_[inds[8]].imag()};
        vec[1] = (__vector double){data_[inds[9]].real(), data_[inds[9]].imag()};
        vec[2] = (__vector double){data_[inds[10]].real(), data_[inds[10]].imag()};
        vec[3] = (__vector double){data_[inds[11]].real(), data_[inds[11]].imag()};
    
        temp[0] = vec_add(temp[0], vec_add(vec_add(complex_mul(vec[0],imag[8],real[8]), complex_mul(vec[1],imag[9],real[9])),vec_add(complex_mul(vec[2],imag[10],real[10]), complex_mul(vec[3],imag[11],real[11]))));
        temp[1] = vec_add(temp[1], vec_add(vec_add(complex_mul(vec[0],imag[24],real[24]), complex_mul(vec[1],imag[25],real[25])),vec_add(complex_mul(vec[2],imag[26],real[26]), complex_mul(vec[3],imag[27],real[27]))));
        temp[2] = vec_add(temp[2], vec_add(vec_add(complex_mul(vec[0],imag[40],real[40]), complex_mul(vec[1],imag[41],real[41])),vec_add(complex_mul(vec[2],imag[42],real[42]), complex_mul(vec[3],imag[43],real[43]))));
        temp[3] = vec_add(temp[3], vec_add(vec_add(complex_mul(vec[0],imag[56],real[56]), complex_mul(vec[1],imag[57],real[57])),vec_add(complex_mul(vec[2],imag[58],real[58]), complex_mul(vec[3],imag[59],real[59]))));
        temp[4] = vec_add(temp[4], vec_add(vec_add(complex_mul(vec[0],imag[72],real[72]), complex_mul(vec[1],imag[73],real[73])),vec_add(complex_mul(vec[2],imag[74],real[74]), complex_mul(vec[3],imag[75],real[75]))));
        temp[5] = vec_add(temp[5], vec_add(vec_add(complex_mul(vec[0],imag[88],real[88]), complex_mul(vec[1],imag[89],real[89])),vec_add(complex_mul(vec[2],imag[90],real[90]), complex_mul(vec[3],imag[91],real[91]))));
        temp[6] = vec_add(temp[6], vec_add(vec_add(complex_mul(vec[0],imag[104],real[104]), complex_mul(vec[1],imag[105],real[105])),vec_add(complex_mul(vec[2],imag[106],real[106]), complex_mul(vec[3],imag[107],real[107]))));
        temp[7] = vec_add(temp[7], vec_add(vec_add(complex_mul(vec[0],imag[120],real[120]), complex_mul(vec[1],imag[121],real[121])),vec_add(complex_mul(vec[2],imag[122],real[122]), complex_mul(vec[3],imag[123],real[123]))));
        temp[8] = vec_add(temp[8], vec_add(vec_add(complex_mul(vec[0],imag[136],real[136]), complex_mul(vec[1],imag[137],real[137])),vec_add(complex_mul(vec[2],imag[138],real[138]), complex_mul(vec[3],imag[139],real[139]))));
        temp[9] = vec_add(temp[9], vec_add(vec_add(complex_mul(vec[0],imag[152],real[152]), complex_mul(vec[1],imag[153],real[153])),vec_add(complex_mul(vec[2],imag[154],real[154]), complex_mul(vec[3],imag[155],real[155]))));
        temp[10] = vec_add(temp[10], vec_add(vec_add(complex_mul(vec[0],imag[168],real[168]), complex_mul(vec[1],imag[169],real[169])),vec_add(complex_mul(vec[2],imag[170],real[170]), complex_mul(vec[3],imag[171],real[171]))));
        temp[11] = vec_add(temp[11], vec_add(vec_add(complex_mul(vec[0],imag[184],real[184]), complex_mul(vec[1],imag[185],real[185])),vec_add(complex_mul(vec[2],imag[186],real[186]), complex_mul(vec[3],imag[187],real[187]))));
        temp[12] = vec_add(temp[12], vec_add(vec_add(complex_mul(vec[0],imag[200],real[200]), complex_mul(vec[1],imag[201],real[201])),vec_add(complex_mul(vec[2],imag[202],real[202]), complex_mul(vec[3],imag[203],real[203]))));
        temp[13] = vec_add(temp[13], vec_add(vec_add(complex_mul(vec[0],imag[216],real[216]), complex_mul(vec[1],imag[217],real[217])),vec_add(complex_mul(vec[2],imag[218],real[218]), complex_mul(vec[3],imag[219],real[219]))));
        temp[14] = vec_add(temp[14], vec_add(vec_add(complex_mul(vec[0],imag[232],real[232]), complex_mul(vec[1],imag[233],real[233])),vec_add(complex_mul(vec[2],imag[234],real[234]), complex_mul(vec[3],imag[235],real[235]))));
        temp[15] = vec_add(temp[15], vec_add(vec_add(complex_mul(vec[0],imag[248],real[248]), complex_mul(vec[1],imag[249],real[249])),vec_add(complex_mul(vec[2],imag[250],real[250]), complex_mul(vec[3],imag[251],real[251]))));

        vec[0] = (__vector double){data_[inds[12]].real(), data_[inds[12]].imag()};
        vec[1] = (__vector double){data_[inds[13]].real(), data_[inds[13]].imag()};
        vec[2] = (__vector double){data_[inds[14]].real(), data_[inds[14]].imag()};
        vec[3] = (__vector double){data_[inds[15]].real(), data_[inds[15]].imag()};

        result = vec_add(temp[0], vec_add(vec_add(complex_mul(vec[0],imag[12],real[12]), complex_mul(vec[1],imag[13],real[13])),vec_add(complex_mul(vec[2],imag[14],real[14]), complex_mul(vec[3],imag[15],real[15]))));
    data_[inds[0]] = std::complex<double>((double)result[0], (double)result[1]);
        result = vec_add(temp[1], vec_add(vec_add(complex_mul(vec[0],imag[28],real[28]), complex_mul(vec[1],imag[29],real[29])),vec_add(complex_mul(vec[2],imag[30],real[30]), complex_mul(vec[3],imag[31],real[31]))));
        data_[inds[1]] = std::complex<double>((double)result[0], (double)result[1]);
        result = vec_add(temp[2], vec_add(vec_add(complex_mul(vec[0],imag[44],real[44]), complex_mul(vec[1],imag[45],real[45])),vec_add(complex_mul(vec[2],imag[46],real[46]), complex_mul(vec[3],imag[47],real[47]))));
        data_[inds[2]] = std::complex<double>((double)result[0], (double)result[1]);
        result = vec_add(temp[3], vec_add(vec_add(complex_mul(vec[0],imag[60],real[60]), complex_mul(vec[1],imag[61],real[61])),vec_add(complex_mul(vec[2],imag[62],real[62]), complex_mul(vec[3],imag[63],real[63]))));
        data_[inds[3]] = std::complex<double>((double)result[0], (double)result[1]);
        result = vec_add(temp[4], vec_add(vec_add(complex_mul(vec[0],imag[76],real[76]), complex_mul(vec[1],imag[77],real[77])),vec_add(complex_mul(vec[2],imag[78],real[78]), complex_mul(vec[3],imag[79],real[79]))));
        data_[inds[4]] = std::complex<double>((double)result[0], (double)result[1]);
        result = vec_add(temp[5], vec_add(vec_add(complex_mul(vec[0],imag[92],real[92]), complex_mul(vec[1],imag[93],real[93])),vec_add(complex_mul(vec[2],imag[94],real[94]), complex_mul(vec[3],imag[95],real[95]))));
        data_[inds[5]] = std::complex<double>((double)result[0], (double)result[1]);
        result = vec_add(temp[6], vec_add(vec_add(complex_mul(vec[0],imag[108],real[108]), complex_mul(vec[1],imag[109],real[109])),vec_add(complex_mul(vec[2],imag[110],real[110]), complex_mul(vec[3],imag[111],real[111]))));
        data_[inds[6]] = std::complex<double>((double)result[0], (double)result[1]);
        result = vec_add(temp[7], vec_add(vec_add(complex_mul(vec[0],imag[124],real[124]), complex_mul(vec[1],imag[125],real[125])),vec_add(complex_mul(vec[2],imag[126],real[126]), complex_mul(vec[3],imag[127],real[127]))));
        data_[inds[7]] = std::complex<double>((double)result[0], (double)result[1]);
        result = vec_add(temp[8], vec_add(vec_add(complex_mul(vec[0],imag[140],real[140]), complex_mul(vec[1],imag[141],real[141])),vec_add(complex_mul(vec[2],imag[142],real[142]), complex_mul(vec[3],imag[143],real[143]))));
        data_[inds[8]] = std::complex<double>((double)result[0], (double)result[1]);
        result = vec_add(temp[9], vec_add(vec_add(complex_mul(vec[0],imag[156],real[156]), complex_mul(vec[1],imag[157],real[157])),vec_add(complex_mul(vec[2],imag[158],real[158]), complex_mul(vec[3],imag[159],real[159]))));
        data_[inds[9]] = std::complex<double>((double)result[0], (double)result[1]);
        result = vec_add(temp[10], vec_add(vec_add(complex_mul(vec[0],imag[172],real[172]), complex_mul(vec[1],imag[173],real[173])),vec_add(complex_mul(vec[2],imag[174],real[174]), complex_mul(vec[3],imag[175],real[175]))));
        data_[inds[10]] = std::complex<double>((double)result[0], (double)result[1]);
        result = vec_add(temp[11], vec_add(vec_add(complex_mul(vec[0],imag[188],real[188]), complex_mul(vec[1],imag[189],real[189])),vec_add(complex_mul(vec[2],imag[190],real[190]), complex_mul(vec[3],imag[191],real[191]))));
        data_[inds[11]] = std::complex<double>((double)result[0], (double)result[1]);
        result = vec_add(temp[12], vec_add(vec_add(complex_mul(vec[0],imag[204],real[204]), complex_mul(vec[1],imag[205],real[205])),vec_add(complex_mul(vec[2],imag[206],real[206]), complex_mul(vec[3],imag[207],real[207]))));
        data_[inds[12]] = std::complex<double>((double)result[0], (double)result[1]);
        result = vec_add(temp[13], vec_add(vec_add(complex_mul(vec[0],imag[220],real[220]), complex_mul(vec[1],imag[221],real[221])),vec_add(complex_mul(vec[2],imag[222],real[222]), complex_mul(vec[3],imag[223],real[223]))));
        data_[inds[13]] = std::complex<double>((double)result[0], (double)result[1]);
        result = vec_add(temp[14], vec_add(vec_add(complex_mul(vec[0],imag[236],real[236]), complex_mul(vec[1],imag[237],real[237])),vec_add(complex_mul(vec[2],imag[238],real[238]), complex_mul(vec[3],imag[239],real[239]))));
        data_[inds[14]] = std::complex<double>((double)result[0], (double)result[1]);
        result = vec_add(temp[15], vec_add(vec_add(complex_mul(vec[0],imag[252],real[252]), complex_mul(vec[1],imag[253],real[253])),vec_add(complex_mul(vec[2],imag[254],real[254]), complex_mul(vec[3],imag[255],real[255]))));
        data_[inds[15]] = std::complex<double>((double)result[0], (double)result[1]);
      };
      apply_lambda(lambda, areg_t<4>({{qubits[0], qubits[1], qubits[2], qubits[3]}}), convert(mat));
      return;
    } else {
      // Lambda function for 4-qubit matrix multiplication
      auto lambda = [&](const areg_t<16> &inds, const cvector_t<data_t> &_mat)->void {
        std::array<std::complex<data_t>, 16> cache;
        for (size_t i = 0; i < 16; i++) {
          const auto ii = inds[i];
          cache[i] = data_[ii];
          data_[ii] = 0.;
        }
        // update state vector
        for (size_t i = 0; i < 16; i++)
          for (size_t j = 0; j < 16; j++)
            data_[inds[i]] += _mat[i + 16 * j] * cache[j];
      };
      apply_lambda(lambda, areg_t<4>({{qubits[0], qubits[1], qubits[2], qubits[3]}}), convert(mat));
      return;
      }
    }
    default: {
      const uint_t DIM = BITS[N];
      // Lambda function for N-qubit matrix multiplication
      auto lambda = [&](const indexes_t &inds, const cvector_t<data_t> &_mat)->void {
        auto cache = std::make_unique<std::complex<data_t>[]>(DIM);
        for (size_t i = 0; i < DIM; i++) {
          const auto ii = inds[i];
          cache[i] = data_[ii];
          data_[ii] = 0.;
        }
        // update state vector
        for (size_t i = 0; i < DIM; i++)
          for (size_t j = 0; j < DIM; j++)
            data_[inds[i]] += _mat[i + DIM * j] * cache[j];
      };
      apply_lambda(lambda, qubits, convert(mat));
    }
  } // end switch
}

template <typename data_t>
void QubitVector<data_t>::apply_diagonal_matrix_ppc64le(const reg_t &qubits,
                                                const cvector_t<double> &diag) {
  const int_t N = qubits.size();
  // Error checking
  #ifdef DEBUG
  check_vector(diag, N);
  #endif

  if (N == 1) {
    apply_diagonal_matrix(qubits[0], diag);
    return;
  }

  if(type_double) {
    auto lambda = [&](const areg_t<2> &inds, const cvector_t<data_t> &_diag)->void {
      for (int_t i = 0; i < 2; ++i) {
        const int_t k = inds[i];
        int_t iv = 0;
        for (int_t j = 0; j < N; j++)
          if ((k & (1ULL << qubits[j])) != 0)
            iv += (1 << j);
        if (_diag[iv] != (data_t) 1.0){
	   std::complex<double> test;
	    __vector double real = (__vector double){_diag[iv].real(), _diag[iv].real()};
	    __vector double imag = (__vector double){_diag[iv].imag(), _diag[iv].imag() * -1};
            __vector double vec = (__vector double){data_[k].real(), data_[k].imag()};

            __vector double result = complex_mul(vec,imag,real);
            data_[k] = std::complex<double>((double)result[0], (double)result[1]);
#if 0
            test = std::complex<double>((double)real[0], (double)real[1]);
	    std::cout << "real " << test << std::endl;
            test = std::complex<double>((double)imag[0], (double)imag[1]);
	    std::cout << "imag " << test << std::endl;
	    std::cout << "diag " << _diag[iv]<< std::endl;

            test = std::complex<double>((double)vec[0], (double)vec[1]);
	    std::cout << "vec " << test << std::endl;
	    std::cout << "data_k " << data_[k]<< std::endl;

            __vector double result = complex_mul(vec,imag,real);
            test = std::complex<double>((double)result[0], (double)result[1]);
	    std::cout << test << std::endl;

            data_[k] *= _diag[iv];
	    std::cout << "k " << data_[k]<< std::endl;
#endif
	}
      }
    };
    apply_lambda(lambda, areg_t<1>({{qubits[0]}}), convert(diag));
  } else {
    auto lambda = [&](const areg_t<2> &inds, const cvector_t<data_t> &_diag)->void {
      for (int_t i = 0; i < 2; ++i) {
        const int_t k = inds[i];
        int_t iv = 0;
        for (int_t j = 0; j < N; j++)
          if ((k & (1ULL << qubits[j])) != 0)
            iv += (1 << j);
        if (_diag[iv] != (data_t) 1.0)
          data_[k] *= _diag[iv];
      }
    };
    apply_lambda(lambda, areg_t<1>({{qubits[0]}}), convert(diag));
  }
}


#endif


//------------------------------------------------------------------------------
} // end namespace QV
//------------------------------------------------------------------------------

// ostream overload for templated qubitvector
template <typename data_t>
inline std::ostream &operator<<(std::ostream &out, const QV::QubitVector<data_t>&qv) {

  out << "[";
  size_t last = qv.size() - 1;
  for (size_t i = 0; i < qv.size(); ++i) {
    out << qv[i];
    if (i != last)
      out << ", ";
  }
  out << "]";
  return out;
}

//------------------------------------------------------------------------------
#endif // end module
