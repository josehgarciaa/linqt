#ifndef SPARSE_MATRIX
#define SPARSE_MATRIX

#include <assert.h> /* assert */
#include <iostream>
#include <vector>
#include <complex>
#include <fstream>
using namespace std;

namespace Sparse
{
bool OPERATOR_FromCSRFile(const std::string& input, int &dim, 
                          std::vector<int> &columns, std::vector<int> &rowIndex, 
                          std::vector<std::complex<double>> &values);
      // Reads a sparse matrix from a file in Compressed Sparse Row (CSR) format.
      // Parameters:
      // - input: The file path to the CSR formatted file.
      // - dim: Reference to an integer where the dimension of the square matrix will be stored.
      // - columns: Reference to a vector of integers to store the column indices of non-zero elements.
      // - rowIndex: Reference to a vector of integers to store the index in 'values' where each row starts.
      // - values: Reference to a vector of complex<double> to store the values of non-zero elements.
      // Returns:
      // - true if the file is successfully read and the data is loaded into the parameters.
      // - false if there's an error opening the file or reading its contents. 
};


class SparseMatrixBase {
public:
    typedef std::complex<double> ValueType;
    typedef std::vector<ValueType> VectorType;
    
    // Constructor to optionally set dimensions during initialization
    SparseMatrixBase(int numRows = 0, int numCols = 0, const std::string& id = "")
        : numRows_(numRows), numCols_(numCols), id_(id) {}

    // Virtual destructor to ensure proper cleanup in derived classes
    virtual ~SparseMatrixBase() {}

    // Accessors
    virtual int numRows() const { return numRows_; }
    virtual int numCols() const { return numCols_; }
    virtual int rank() const {
        return (numRows_ > numCols_) ? numCols_ : numRows_;
    }
    
    // Mutators
    virtual void setDimensions(int numRows, int numCols) {
        if (numRows < 0 || numCols < 0) {
            std::cerr << "Error: Dimensions must be non-negative." << std::endl;
            return;
        }
        numRows_ = numRows;
        numCols_ = numCols;
    }

    virtual void setId(const std::string& id) { id_ = id; }
    virtual std::string id() const { return id_; }

    // Utility functions
    virtual bool isIdentity(){ return (bool)( ID()=="1"); };

    // Add virtual functions for mathematical operations here
    // Example:
    // virtual SparseMatrixBase add(const SparseMatrixBase& other) const = 0;
    // virtual SparseMatrixBase multiply(const SparseMatrixBase& other) const = 0;

protected:
    int numRows_, numCols_;
    std::string id_;

    // Add common data structures for sparse representation here
    // For example:
    // VectorType values_;
    // std::vector<int> columns_, rowIndex_;
};

class SparseMatrixType_BASE
{
	
	
public:

	
	  int numRows() { return numRows_; };
	  int numCols() { return numCols_; };
	  int rank() { return ((this->numRows() > this->numCols()) ? this->numCols() : this->numRows()); };
	  void setDimensions(const int numRows, const int numCols)
	  {
		numRows_ = numRows;
		numCols_ = numCols;
	  };
		

};


// MKL LIBRARIES
#define MKL_Complex16 complex<double>
#include "mkl.h"
#include "mkl_spblas.h"
class SparseMatrixType  : public SparseMatrixType_BASE
{
public:
  SparseMatrixType()
  {
    descr.type = SPARSE_MATRIX_TYPE_HERMITIAN;
    descr.mode = SPARSE_FILL_MODE_UPPER;
    descr.diag = SPARSE_DIAG_NON_UNIT;
  }
  
  inline
  matrix_descr& mkl_descr()  { return  descr; }; 

  sparse_matrix_t& mkl_matrix()  { return Matrix; };
  
  string matrixType() const { return "CSR Matrix from MKL Library."; };
  void Multiply(const value_t a, const value_t *x, const value_t b, value_t *y);
  void Multiply(const value_t a, const vector_t& x, const value_t b, vector_t& y);
  void Rescale(const value_t a,const value_t b);
  inline 
  void Multiply(const value_t *x, value_t *y){ Multiply(value_t(1,0),x,value_t(0,0),y);};
  inline 
  void Multiply(const vector_t& x, vector_t& y){ Multiply(value_t(1,0),x,value_t(0,0),y);};


  void BatchMultiply(const int batchSize, const value_t a, const value_t *x, const value_t b, value_t *y);

  void ConvertFromCOO(vector<int> &rows, vector<int> &cols, vector<complex<double> > &vals);
  void ConvertFromCSR(vector<int> &rowIndex, vector<int> &cols, vector<complex<double> > &vals);

private:
  struct matrix_descr descr;
  sparse_matrix_t Matrix;
  vector<int> rows_;
  vector<int> cols_;
  vector<complex<double> > vals_;
};



class SparseMatrixBuilder
{
public:
  void setSparseMatrix(SparseMatrixType *b)
  {
    _matrix_type = b;
  };

public:
  void BuildOPFromCSRFile(const std::string input)
  {
    vector<int> columns, rowIndex;
    vector<complex<double> > values;
    int dim;
    Sparse::OPERATOR_FromCSRFile(input, dim, columns, rowIndex, values);
    _matrix_type->setDimensions(dim, dim);
    _matrix_type->ConvertFromCSR(rowIndex, columns, values);
    std::cout << "OPERATOR SUCCESSFULLY BUILD" << std::endl;
  }
  SparseMatrixType *_matrix_type;
};

#endif
