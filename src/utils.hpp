#ifndef __GmpEigenMatrix_Utils_HPP__
#define __GmpEigenMatrix_Utils_HPP__
#include "mex.h"

#include <iostream>
#include <vector>
#include <iomanip>
#include <Eigen/Sparse>
#include <Eigen/MPRealSupport>
#include <Eigen/LU>

/*
  This file contains several utilities which are used by the GmpEigenMatrix
  library.
*/

//using namespace std;
//using namespace mpfr;
//using namespace Eigen;


// We define a type for all matrix indices
//typedef long int IndexType;
typedef Eigen::Matrix<mpfr::mpreal,Eigen::Dynamic,Eigen::Dynamic>::Index IndexType;


// WARNING : The mpreal library uses a different default precision for all threads (!!!)
//   By default, this means that all threads can work with their own precision
//   level (even though with openmp we typicallz don't control which part of which
//   process goes to which thread). The following precision setting function
//   should thus be used instead of the mpreal one in order to make sure no
//   thread is left out upon setting of a new default precision level!
inline void setDefaultPrecForAllThreads(const mp_prec_t& precision)
{
    #pragma omp parallel
    {
        mpfr::mpreal::set_default_prec(precision);
  }
}


// This function is to deal with some old compilers who don't have std::to_string.
// To be removed eventually...
template <typename T> inline std::string ToString(const T& val)
{
    std::stringstream stream;
    stream << val;
    return stream.str();
}

// This function extracts exactly 15 digits from a double
inline std::string toString15(const double& val)
{
    // Note : the change of precision here doesn't affect other streams.
    std::stringstream stream;
    stream << std::setprecision(15) << val;
    return stream.str();
}

// This function produces a string of n spaces
inline std::string nSpaces(const int& n)
{
    std::string result("");
    for (int i(0); i < n; ++i)
        result = result + " ";
    return result;
}

// This function produces a string of n zeros
inline std::string nZeros(const int& n)
{
    std::string result("");
    for (int i(0); i < n; ++i)
        result = result + "0";
    return result;
}

// Here is a function which attempts to extract the n first digits of an mpreal
// number (with some rounding). If the precision of the number is smaller, it
// pads the rest with spaces. The function also takes a flag (modified by
// reference) that signals whether the rounding procedure created one more
// leading digit (e.g. it 99.999 is rounded to 100.0 at n=4). When this is the
// case, we know that if we had asked one more digit, it would have been a 0
// (i.e. the next digit in 100.0 as an approximation of 99.999 is 0, thus giving
// 100.00. We only get a non-zero digit is we go to one more digit.)
//
// The output of this function is really pure digits : it contains neither sign,
// nor dot nor exponent.
//
// NOTE : This function follows some IEEE rounding mode instead of the usual
// mathematical rounding, also known as round half away from zero. Therefore,
// numbers like 2.5 can be rounded down to 2 instead of being rounded up to 3.
//
// We could obtain correct rounding here by first computing the digits with
// a floor-like rounding. If the last digit is equal to 5, then we know that the
// opposite rounding is the good one, otherwise we already have the right rounding...
std::string firstDigitsToString(const mpfr::mpreal& x, const int& n, bool& oneMoreDigitFlag);



// Here is a function which produces a formatted output. It produces a string
// of exactly width characters that best represent the number under some
// constraints. The exponential shift can be used to inform the procedure that a
// global power of 10 will be added to whatever is printed. The dotPosition
// informs where the dot should be located with respect to the beginning of the
// string. Only digits within the precision of the number are returned. Un-
// defined terminating digits are left blank. Zeros before the first significant
// digit are also replaced by blanks. The last digit of the number is always
// rounded.
// Warning : the behavior is not guaranteed when width is too small (e.g. < 3)
std::string mprealToString(const mpfr::mpreal& x, int width, const int& exponentShift, int dotPosition, const bool& includeSign = true);



// Here is another function which produces a formatted output. It produces a
// string of exactly width characters that best represent the number. The first
// character is reserved for a sign, all the rest is free. The procedure
// automatically switches to scientific notation if needed, so as to best
// represent the number.
// Warning : the behavior is not guaranteed when width is too small (e.g. < 3)
// WARNING : The behavior was not tested in presence of rounding ! (the exponent
//           could have changed... -> we need to check that)
std::string mprealToString2(const mpfr::mpreal& x, const int& width, const bool& padOnTheRight = 0);




/* Now we define some mpreal matrix constructors */

// The zero matrix
inline Eigen::Matrix<mpfr::mpreal,Eigen::Dynamic,Eigen::Dynamic> zeroMatrix(const IndexType& rows, const IndexType& cols, const double& precision = 1.0)
{
    Eigen::Matrix<mpfr::mpreal,Eigen::Dynamic,Eigen::Dynamic> matrix;

    mp_prec_t precisionInBits(mpfr::digits2bits(precision));
    //mpfr::mpreal::set_default_prec(mpfr::digits2bits(precision));
    matrix.resize(rows, cols);
    for (IndexType j(0); j < cols; ++j)
        for (IndexType i(0); i < rows; ++i)
            matrix(i,j) = mpfr::mpreal("0", precisionInBits);

    return matrix;
}

// The identity matrix
inline Eigen::Matrix<mpfr::mpreal,Eigen::Dynamic,Eigen::Dynamic> identityMatrix(const IndexType& dim, const double& precision = 1.0)
{
    Eigen::Matrix<mpfr::mpreal,Eigen::Dynamic,Eigen::Dynamic> matrix;

    mp_prec_t precisionInBits(mpfr::digits2bits(precision));
    //mpfr::mpreal::set_default_prec(mpfr::digits2bits(precision));
    matrix.resize(dim, dim);
    for (IndexType j(0); j < dim; ++j)
        for (IndexType i(0); i < dim; ++i)
        {
            if (i == j)
                matrix(i,j) = mpfr::mpreal("1", precisionInBits);
            else
                matrix(i,j) = mpfr::mpreal("0", precisionInBits);
        }

    return matrix;
}

// The constant matrix
inline Eigen::Matrix<mpfr::mpreal,Eigen::Dynamic,Eigen::Dynamic> constantMatrix(const IndexType& rows, const IndexType& cols, const mpfr::mpreal& value)
{
    Eigen::Matrix<mpfr::mpreal,Eigen::Dynamic,Eigen::Dynamic> matrix;

    matrix.resize(rows, cols);
    for (IndexType j(0); j < cols; ++j)
        for (IndexType i(0); i < rows; ++i)
            matrix(i,j) = value;

    return matrix;
}

// The constant matrix in "sparse" format (This is actually a very bad format to
// use with value != 0, but matlab does produce such objects, as the result of
// the element-wise exponent on a sparse matrix)
inline Eigen::SparseMatrix<mpfr::mpreal> sparseConstantMatrix(const IndexType& rows, const IndexType& cols, const mpfr::mpreal& value)
{
    Eigen::SparseMatrix<mpfr::mpreal> matrix;

    matrix.resize(rows, cols);
    if (value == 0)
        return matrix;

    for (IndexType j(0); j < cols; ++j)
        for (IndexType i(0); i < rows; ++i)
            matrix.insert(i,j) = value;

    return matrix;
}



/* Here are some procedures to send and receive tables to matlab */

// This function sends a vector of type T to matlab (after recasting the types to
// doubles)
template <typename T> inline mxArray* vectorToMatlabDoubles(const std::vector<T>& vec)
{
    /* Get the size of the vector */
    mwSize m(vec.size());
    mwSize n(1);

    mxArray* plhs;
    plhs = mxCreateNumericMatrix(m, n, mxDOUBLE_CLASS, mxREAL);

    double* pointer(mxGetPr(plhs));
    for (mwIndex j = 0; j < vec.size(); ++j) {
        pointer[j] = (double)vec[j];
    }

    return plhs;
}

// This function sends a matrix of type T to matlab (after recasting the types to
// doubles)
template <typename T> inline mxArray* matrixToMatlabDoubles(const std::vector< std::vector <T> >& mat)
{
    /* Get the size of the matrix */
    mwSize m(mat.size());
    mwSize n(mat[0].size());

    mxArray* plhs;
    plhs = mxCreateNumericMatrix(m, n, mxDOUBLE_CLASS, mxREAL);

    double* pointer(mxGetPr(plhs));
    for (mwIndex j = 0; j < n; ++j) {
        // We first iterate on the columns
        for (mwIndex i = 0; i < m; ++i) {
            pointer[i] = (double)mat[i][j];
        }
        pointer += m;
    }

    return plhs;
}

// This function sends a matrix of type T to matlab (after recasting the types to
// doubles)
template <typename T> inline mxArray* matrixToMatlabDoubles(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mat)
{
    /* Get the size of the matrix */
    mwSize m(mat.rows());
    mwSize n(mat.cols());

    mxArray* plhs;
    plhs = mxCreateNumericMatrix(m, n, mxDOUBLE_CLASS, mxREAL);

    double* pointer(mxGetPr(plhs));
    for (mwIndex j = 0; j < n; ++j) {
        // We first iterate on the columns
        for (mwIndex i = 0; i < m; ++i) {
            pointer[i] = (double)mat(i,j);
        }
        pointer += m;
    }

    return plhs;
}

// This function sends a sparse matrix of type T to matlab (after recasting the
// types to doubles)
template <typename T> inline mxArray* matrixToMatlabDoubles(const Eigen::SparseMatrix<T>& mat)
{
    /* Get the size of the matrix */
    mwSize m(mat.rows());
    mwSize n(mat.cols());

    mxArray* plhs;
    plhs = mxCreateSparse(m, n, mat.nonZeros(), mxREAL);

    // We get the pointers to where we need to input data
    double *sr;
    mwIndex *irs,*jcs;
    sr  = mxGetPr(plhs);
    irs = mxGetIr(plhs);
    jcs = mxGetJc(plhs);

    // We copy the data into the matlab matrix
    jcs[0] = 0;
    for (IndexType k = 0; k < mat.outerSize(); ++k) {
        IndexType index(0);
        for (typename Eigen::SparseMatrix<T>::InnerIterator it(mat,k); it; ++it) {
            irs[index] = it.row();
            sr[index] = double(it.value());
            ++index;
        }
        jcs[k+1] = jcs[k] + index;
        irs += index;
        sr += index;
    }

    return plhs;
}


// This function creates a vector of given type from a matlab table of doubles
// If works for column and line vectors equally, and if given a matrix, sees it
// as one big vector.
template <typename T> inline std::vector<T> matlabDoublesToVector(const mxArray* prhs)
{
    /* Get the size and pointers to input data */
    mwSize m(mxGetM(prhs));
    mwSize n(mxGetN(prhs));
    double* pr(mxGetPr(prhs)); // We only care about the real data (imaginary part ignored if any)

    std::vector<T> vec;

    // Now we copy the data over
    for (mwIndex j = 0; j < n; ++j) {
        // We first iterate on the rows
        for (mwIndex i = 0; i < m; ++i) {
            vec.push_back((T)pr[i]);
        }
        pr += m;
    }

    return vec;
}

// This function creates a matrix of given type from a matlab table of doubles
// NOTE : Actually, the compiler doesn't seem to be able to distinguish this
//    function from the next one... so we changed its name...
template <typename T> inline std::vector< std::vector<T> > matlabDoublesToMatrixVV(const mxArray* prhs)
{
    /* Get the size and pointers to input data */
    mwSize m(mxGetM(prhs));
    mwSize n(mxGetN(prhs));
    double* pr(mxGetPr(prhs)); // We only care about the real data (imaginary part ignored if any)

    // We initialize the table
    std::vector< std::vector<T> > mat(m, std::vector<T> (n, (T)0));

    // Now we copy the data over
    for (mwIndex j = 0; j < n; ++j) {
        // We first iterate on the rows
        for (mwIndex i = 0; i < m; ++i) {
            mat[i][j] = (T)pr[i];
        }
        pr += m;
    }

    return mat;
}

// This function sends a matrix of type T to matlab (after recasting the types to
// doubles)
template <typename T> inline Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> matlabDoublesToMatrix(const mxArray* prhs)
{
    /* Get the size and pointers to input data */
    mwSize m(mxGetM(prhs));
    mwSize n(mxGetN(prhs));
    double* pr(mxGetPr(prhs)); // We only care about the real data (imaginary part ignored if any)

    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> mat;

    // Now we copy the data over
    for (mwIndex j = 0; j < n; ++j) {
        // We first iterate on the rows
        for (mwIndex i = 0; i < m; ++i) {
            mat(i,j) = (T)pr[i];
        }
        pr += m;
    }

    return mat;
}


#endif // __GmpEigenMatrix_Utils_HPP__
