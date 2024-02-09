#include "sparse_matrix.hpp"






bool Sparse::OPERATOR_FromCSRFile(const std::string& input, int &dim, 
                                  std::vector<int> &columns, std::vector<int> &rowIndex, 
                                  std::vector<std::complex<double>> &values) {
  // Inform the user about the start of the file reading process.
  std::cout << "\nReading the CSR file located at: " << input << std::endl;


  // Attempt to open the matrix file.
  std::ifstream matrix_file(input);
  if (!matrix_file) {
    // If the file cannot be opened, log an error message and return false.
    std::cerr << "ERROR: Cannot open the matrix file" << std::endl;
    return false;
  }

  // Read the dimension of the matrix and the number of non-zero values (nnz).
  // The file format expects the first line to contain these two integers.
  int nnz;
  matrix_file >> dim >> nnz;
  if (!matrix_file) {
    // If reading the dimensions fails, log an error and exit.
    std::cerr << "ERROR: Failed to read dimensions from the matrix file" << std::endl;
    return false;
  }

  // Resize the vectors to accommodate the data based on the read dimensions.
  values.resize(nnz);
  columns.resize(nnz); 
  rowIndex.resize(dim + 1);

  // Read the complex values into the 'values' vector.
  double rev, imv;
  for (int i = 0; i < nnz && matrix_file >> rev >> imv; ++i) {
    values[i] = std::complex<double>(rev, imv);
  }

  // Read the column indices into the 'columns' vector.
  for (int i = 0; i < nnz && matrix_file >> columns[i]; ++i);

  // Read the row start indices into the 'rowIndex' vector.
  for (int i = 0; i <= dim && matrix_file >> rowIndex[i]; ++i);

  // Check for any errors that might have occurred during file reading.
  if (matrix_file.bad()) {
    std::cerr << "ERROR: An error occurred while reading the matrix file" << std::endl;
    return false;
  }

  // Successfully read the matrix file.
  std::cout << "Finish reading. Status: SUCCEEDED" << std::endl;
  return true;

};
