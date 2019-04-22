#ifndef __SparseGmpEigenMatrix_HPP__
#define __SparseGmpEigenMatrix_HPP__
#include "mex.h"

#include <iostream>
#include <vector>
#include <Eigen/Sparse>
#include <Eigen/MPRealSupport>
#include <Eigen/OrderingMethods>
#include <Eigen/SparseLU>
#include "utils.hpp"
#include <Spectra/SymEigsSolver.h>
#include <Spectra/GenEigsSolver.h>
#include <Spectra/MatOp/SparseSymMatProd.h>
#include <Spectra/MatOp/SparseGenMatProd.h>
#include <Spectra/SymEigsShiftSolver.h>
#include <Spectra/GenEigsRealShiftSolver.h>
#include <Spectra/MatOp/SparseSymShiftSolve.h>
#include <Spectra/MatOp/SparseGenRealShiftSolve.h>

/*
  This file contains the description of our c++ class, including all its
  required procedures. It also includes a function which allows for the creation
  of a class instance from numerical data received from matlab, as well as a few
  output functions tailored for matlab.
*/

/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


class GmpEigenMatrix;

/* The class that we are interfacing to */
class SparseGmpEigenMatrix
{
    friend class GmpEigenMatrix;
public:
    /* The main constructor */
    SparseGmpEigenMatrix() : isComplex(false) {};
    /* Construction from one mpreal number */
    SparseGmpEigenMatrix(const mpfr::mpreal& value) : matrixR(1,1), isComplex(false) {matrixR.coeffRef(0,0) = value;};
    /* Construction from two mpreal numbers */
    SparseGmpEigenMatrix(const mpfr::mpreal& valueR, const mpfr::mpreal& valueI) : matrixR(1,1), matrixI(1,1), isComplex(valueI != 0) {
        if (valueR != 0)
            matrixR.coeffRef(0,0) = valueR;
        if (valueI != 0)
            matrixI.coeffRef(0,0) = valueI;
    };
    /* The copy constructor */
    SparseGmpEigenMatrix(const SparseGmpEigenMatrix& a) : matrixR(a.matrixR), matrixI(a.matrixI), isComplex(a.isComplex) {};
    /* Construction from a dense GmpEigenMatrix */
    SparseGmpEigenMatrix(const GmpEigenMatrix& a, const mpfr::mpreal& threshold = 0);
    /* Construction from a matlab struct containing all the class information. */
    SparseGmpEigenMatrix(const mxArray* prhs);
    /* Construction from a matlab table.*/
    SparseGmpEigenMatrix(const mxArray* rows, const mxArray* cols, const mxArray* values, const IndexType& m, const IndexType& n, const int& precision);
    /* Construction from a gem table.*/
    SparseGmpEigenMatrix(const mxArray* rows, const mxArray* cols, const GmpEigenMatrix& values, const IndexType& m, const IndexType& n, const int& precision);

    /* Destructor */
    virtual ~SparseGmpEigenMatrix() {};



    /* ------------------------------
       | Some procedures for matlab |
       ------------------------------ */

    /* Saving procedure for matlab */
    mxArray* saveobj() const;

    /* Display function, in which each number in the matrix is printed
       independently of the other ones. Here, the width include the sign, digits, dot,
       and exponent if they apply. Each number is processed independently. */
    void displayIndividual(int width) const;

    /* Extracting the raw data: the following function extracts the internal
       data from the class, and puts it in a basic c-type table to be accessed
       by matlab. The result is a sparse array in matlab. */
    mxArray* toDouble() const;

    /* This function returns a cell array with strings describing each matrix
       element */
    mxArray* toStrings(const int& precision) const;

    /* This function returns the number of significant digits of the real part in
       matlab format */
    //Eigen::Matrix <double, Eigen::Dynamic, Eigen::Dynamic> precision() const;
    mxArray* precision() const;

    /* This function removes all components in the matrix whose magnitude is smaller
       than the provided tolerance */
    SparseGmpEigenMatrix clean(const GmpEigenMatrix& tolerance) const;
    SparseGmpEigenMatrix& clean_new(const GmpEigenMatrix& tolerance) const;

    /* This function extracts the nonzero elements from the matrix */
    GmpEigenMatrix find(std::vector < double >& rows, std::vector < double >& cols) const;
    GmpEigenMatrix& find_new(std::vector < double >& rows, std::vector < double >& cols) const;

    /* This function returns the size in matlab format */
    mxArray* size() const;

    /* This function changes the size of the matrix */
    inline void resize(const IndexType& m2, const IndexType& n2)
    {
        matrixR.conservativeResize(m2,n2);
        if (isComplex)
            matrixI.conservativeResize(m2,n2);
    }

    /* Reshaping a matrix */
    SparseGmpEigenMatrix reshape(const IndexType& m2, const IndexType& n2) const;
    SparseGmpEigenMatrix& reshape_new(const IndexType& m2, const IndexType& n2) const;

    /* Concatenates two matrices */
    SparseGmpEigenMatrix vertcat(const SparseGmpEigenMatrix& b) const;
    SparseGmpEigenMatrix& vertcat_new(const SparseGmpEigenMatrix& b) const;
    SparseGmpEigenMatrix horzcat(const SparseGmpEigenMatrix& b) const;
    SparseGmpEigenMatrix& horzcat_new(const SparseGmpEigenMatrix& b) const;

    /* Extracting a sub-matrix */
    SparseGmpEigenMatrix block(const IndexType& i, const IndexType& j, const IndexType& rows, const IndexType& cols) const;
    SparseGmpEigenMatrix subsref(const std::vector<std::vector<IndexType> >& indices) const;
    SparseGmpEigenMatrix& subsref_new(const std::vector<std::vector<IndexType> >& indices) const;
    SparseGmpEigenMatrix subsref(const std::vector<IndexType>& indicesA, const std::vector<IndexType>& indicesB) const;
    SparseGmpEigenMatrix& subsref_new(const std::vector<IndexType>& indicesA, const std::vector<IndexType>& indicesB) const;

    /* Assigning value to a sub-matrix */
    void subsasgn(const std::vector<IndexType>& indices, const SparseGmpEigenMatrix& b);
    void subsasgn(const std::vector<IndexType>& indicesA, const std::vector<IndexType>& indicesB, const SparseGmpEigenMatrix& b);

    /* Diagonal matrix and elements */
    SparseGmpEigenMatrix diagCreate(const IndexType& k) const;
    SparseGmpEigenMatrix& diagCreate_new(const IndexType& k) const;
    SparseGmpEigenMatrix diagExtract(const IndexType& k) const;
    SparseGmpEigenMatrix& diagExtract_new(const IndexType& k) const;


    /* ------------------------
       | Some class operators |
       ------------------------ */

    /* Operators returning an object are defined for c++ calls. When the result
       of an operation needs to be returned to matlab, it must be allocated
       dynamically. This is what the functions with the "_new" ending do : they
       put their result in a dynamically allocated piece of memory which is not
       freed upon retunring to matlab.
           Note that for a proper memory management, the reference created by
       these functions should be passed to matlab after making a call to the
       "createMatlabIdFromObj" function. This allows to make sure that dynamically
       allocated memory is freed properly.

       NOTE : Some operations here don't preserve the zeros. Therefore, they
         are not well suited for sparse matrices : their result may not be
         sparse anymore (e.g. the exponental function tunrs all zeros to ones).
         We leave it to the matlab interface to choose whether or not to encode
         the result of these functions into sparse matrices or not. The results
         of the functions here are sparse matrices. */

    // Unary substraction operators b = -a
    SparseGmpEigenMatrix operator-() const;
    SparseGmpEigenMatrix& uminus_new() const;

    // Transpositions b = a.' and b = a'
    SparseGmpEigenMatrix transpose() const;
    SparseGmpEigenMatrix& transpose_new() const;
    SparseGmpEigenMatrix ctranspose() const;
    SparseGmpEigenMatrix& ctranspose_new() const;

    // Element-wise complex conjugation
    SparseGmpEigenMatrix conj() const;
    SparseGmpEigenMatrix& conj_new() const;

    // Real and imaginary parts
    SparseGmpEigenMatrix real() const;
    SparseGmpEigenMatrix& real_new() const;
    SparseGmpEigenMatrix imag() const;
    SparseGmpEigenMatrix& imag_new() const;

    // Various integer parts
    SparseGmpEigenMatrix round() const;
    SparseGmpEigenMatrix& round_new() const;
    SparseGmpEigenMatrix floor() const;
    SparseGmpEigenMatrix& floor_new() const;
    SparseGmpEigenMatrix ceil() const;
    SparseGmpEigenMatrix& ceil_new() const;
    SparseGmpEigenMatrix trunc() const;
    SparseGmpEigenMatrix& trunc_new() const;


    // Addition a += b and c = a+b
    SparseGmpEigenMatrix& operator+=(const SparseGmpEigenMatrix& b);
    SparseGmpEigenMatrix operator+(const SparseGmpEigenMatrix& b) const;
    SparseGmpEigenMatrix& plus_new(const SparseGmpEigenMatrix& b) const;

    // Substraction a -= b and c = a-b
    SparseGmpEigenMatrix& operator-=(const SparseGmpEigenMatrix& b);
    SparseGmpEigenMatrix operator-(const SparseGmpEigenMatrix& b) const;
    SparseGmpEigenMatrix& minus_new(const SparseGmpEigenMatrix& b) const;

    // Element-wise multiplication c = a.*b
    SparseGmpEigenMatrix times(const SparseGmpEigenMatrix& b) const;
    SparseGmpEigenMatrix& times_new(const SparseGmpEigenMatrix& b) const;
    // Element-wise multiplication with a full matrix
    SparseGmpEigenMatrix times_sf(const GmpEigenMatrix& b) const;
    SparseGmpEigenMatrix& times_sf_new(const GmpEigenMatrix& b) const;

    // Element-wise division on the right c = a./b
    SparseGmpEigenMatrix rdivide(const SparseGmpEigenMatrix& b) const;
    SparseGmpEigenMatrix& rdivide_new(const SparseGmpEigenMatrix& b) const;

    // Absolute value (or complex magnitude) c = abs(a)
    SparseGmpEigenMatrix abs() const;
    SparseGmpEigenMatrix& abs_new() const;

    // Complex norm b = angle(a)
    SparseGmpEigenMatrix angle() const;
    SparseGmpEigenMatrix& angle_new() const;

    // Element-wise exponential b = exp(a)
    SparseGmpEigenMatrix exp() const;
    SparseGmpEigenMatrix& exp_new() const;

    // Element-wise natural logarithms b = log(a)
    SparseGmpEigenMatrix log() const;
    SparseGmpEigenMatrix& log_new() const;

    // Element-wise power c = a.^b
    SparseGmpEigenMatrix power(const GmpEigenMatrix& b) const;
    SparseGmpEigenMatrix& power_new(const GmpEigenMatrix& b) const;

    // square roots
    SparseGmpEigenMatrix sqrt() const;
    SparseGmpEigenMatrix& sqrt_new() const;

    // matrix multiplication a *= b and c = a*b
    SparseGmpEigenMatrix& operator*=(const SparseGmpEigenMatrix& b);
    SparseGmpEigenMatrix operator*(const SparseGmpEigenMatrix& b) const;
    SparseGmpEigenMatrix& mtimes_new(const SparseGmpEigenMatrix& b) const;
    // multiplication with a full matrix
    GmpEigenMatrix mtimes_sf(const GmpEigenMatrix& b) const;
    GmpEigenMatrix& mtimes_sf_new(const GmpEigenMatrix& b) const;

    // tensor product c = kron(a,b)
    SparseGmpEigenMatrix kron(const SparseGmpEigenMatrix& b) const;
    SparseGmpEigenMatrix& kron_new(const SparseGmpEigenMatrix& b) const;
    // tensor product with a full matrix
    SparseGmpEigenMatrix kron_sf(const GmpEigenMatrix& b) const;
    SparseGmpEigenMatrix& kron_sf_new(const GmpEigenMatrix& b) const;


    /* ----------------------------------------------------------
       | Some mathematical functions defined for real arguments |
       |       (the real test is made on the matlab side.       |
       |       Here we simply ignore the imaginary parts)       |
       ---------------------------------------------------------- */

    // cubic roots
    SparseGmpEigenMatrix cbrt() const;
    SparseGmpEigenMatrix& cbrt_new() const;

    // trigonometric functions
    SparseGmpEigenMatrix sin() const;
    SparseGmpEigenMatrix& sin_new() const;
    SparseGmpEigenMatrix cos_nz() const;
    SparseGmpEigenMatrix tan() const;
    SparseGmpEigenMatrix& tan_new() const;
    SparseGmpEigenMatrix asin() const;
    SparseGmpEigenMatrix& asin_new() const;
    SparseGmpEigenMatrix atan() const;
    SparseGmpEigenMatrix& atan_new() const;



    /* ---------------------------------------
       |   Some linear algebra operations    |
       --------------------------------------- */

    IndexType rank() const;
    GmpEigenMatrix mldivide_sf(const GmpEigenMatrix& b) const;
    GmpEigenMatrix& mldivide_sf_new(const GmpEigenMatrix& b) const;
    SparseGmpEigenMatrix mldivide(const SparseGmpEigenMatrix& b) const;
    SparseGmpEigenMatrix& mldivide_new(const SparseGmpEigenMatrix& b) const;
    SparseGmpEigenMatrix inv() const;
    SparseGmpEigenMatrix& inv_new() const;
    GmpEigenMatrix eigs(const long int& nbEigenvalues, GmpEigenMatrix& V, const long int& type, const GmpEigenMatrix& sigma) const;
    GmpEigenMatrix& eigs_new(const long int& nbEigenvalues, GmpEigenMatrix& V, const long int& type, const GmpEigenMatrix& sigma) const;



    /* ---------------------------------------
       | Some tests and comparison operators |
       --------------------------------------- */

    // <, <=, >, >=, ==, !=, isequal tests
    Eigen::SparseMatrix <bool> operator<(const SparseGmpEigenMatrix& b) const;
    Eigen::SparseMatrix <bool> operator<=(const SparseGmpEigenMatrix& b) const;
    Eigen::SparseMatrix <bool> operator>(const SparseGmpEigenMatrix& b) const;
    Eigen::SparseMatrix <bool> operator>=(const SparseGmpEigenMatrix& b) const;
    Eigen::SparseMatrix <bool> eq(const SparseGmpEigenMatrix& b) const;
    Eigen::SparseMatrix <bool> ne(const SparseGmpEigenMatrix& b) const;
    Eigen::SparseMatrix <bool> isnan() const;
    Eigen::SparseMatrix <bool> isinf() const;
    bool identicalValues(const SparseGmpEigenMatrix& b) const;
    bool identicalValuesNaNok(const SparseGmpEigenMatrix& b) const;

    // isreal
    inline bool isreal() const { return (!isComplex); }

    // global minimum of the real part
    inline mpfr::mpreal globalMinReal() const
    {
        mpfr::mpreal minimum(mpfr::mpreal("Inf"));
        for (IndexType k = 0; k < matrixR.outerSize(); ++k)
            for (Eigen::SparseMatrix<mpfr::mpreal>::InnerIterator it(matrixR,k); it; ++it)
                if (it.value() < minimum)
                    minimum = it.value();
        return minimum;
    }
    // global maximum of the real part
    inline mpfr::mpreal globalMaxReal() const
    {
        mpfr::mpreal maximum(mpfr::mpreal("-Inf"));
        for (IndexType k = 0; k < matrixR.outerSize(); ++k)
            for (Eigen::SparseMatrix<mpfr::mpreal>::InnerIterator it(matrixR,k); it; ++it)
                if (it.value() > maximum)
                    maximum = it.value();
        return maximum;
    }

    // symmetry tests
    bool issymmetric() const;
    bool ishermitian() const;

    // column-wise minimum b = min(a)
    SparseGmpEigenMatrix colMin(std::vector<IndexType>& indices) const;
    inline SparseGmpEigenMatrix colMin() const
    {
        std::vector<IndexType> indices;
        return colMin(indices);
    }
    SparseGmpEigenMatrix& colMin_new(std::vector<IndexType>& indices) const;
    inline SparseGmpEigenMatrix& colMin_new() const
    {
        std::vector<IndexType> indices;
        return colMin_new(indices);
    }

    // line-wise minimum b = min(a,[],2)
    SparseGmpEigenMatrix rowMin(std::vector<IndexType>& indices) const;
    inline SparseGmpEigenMatrix rowMin() const
    {
        std::vector<IndexType> indices;
        return rowMin(indices);
    }
    SparseGmpEigenMatrix& rowMin_new(std::vector<IndexType>& indices) const;
    inline SparseGmpEigenMatrix& rowMin_new() const
    {
        std::vector<IndexType> indices;
        return rowMin_new(indices);
    }

    // element-wise minimum c = min(a, b)
    SparseGmpEigenMatrix ewMin(const SparseGmpEigenMatrix& b) const;
    SparseGmpEigenMatrix& ewMin_new(const SparseGmpEigenMatrix& b) const;

    // column-wise maximum b = max(a)
    SparseGmpEigenMatrix colMax(std::vector<IndexType>& indices) const;
    inline SparseGmpEigenMatrix colMax() const
    {
        std::vector<IndexType> indices;
        return colMax(indices);
    }
    SparseGmpEigenMatrix& colMax_new(std::vector<IndexType>& indices) const;
    inline SparseGmpEigenMatrix& colMax_new() const
    {
        std::vector<IndexType> indices;
        return colMax_new(indices);
    }

    // line-wise maximum b = max(a,[],2)
    SparseGmpEigenMatrix rowMax(std::vector<IndexType>& indices) const;
    inline SparseGmpEigenMatrix rowMax() const
    {
        std::vector<IndexType> indices;
        return rowMax(indices);
    }
    SparseGmpEigenMatrix& rowMax_new(std::vector<IndexType>& indices) const;
    inline SparseGmpEigenMatrix& rowMax_new() const
    {
        std::vector<IndexType> indices;
        return rowMax_new(indices);
    }

    // element-wise maximum c = max(a, b)
    SparseGmpEigenMatrix ewMax(const SparseGmpEigenMatrix& b) const;
    SparseGmpEigenMatrix& ewMax_new(const SparseGmpEigenMatrix& b) const;


    // sums and products
    // SparseGmpEigenMatrix colSum() const;
    // SparseGmpEigenMatrix& colSum_new() const;
    // SparseGmpEigenMatrix rowSum() const;
    // SparseGmpEigenMatrix& rowSum_new() const;
    SparseGmpEigenMatrix colProd() const;
    SparseGmpEigenMatrix& colProd_new() const;
    SparseGmpEigenMatrix rowProd() const;
    SparseGmpEigenMatrix& rowProd_new() const;

    // sorting functions
    SparseGmpEigenMatrix sort(const int& dim, const int& type, std::vector < std::vector < IndexType > >& index, std::vector < IndexType >& nbNegatives) const;
    SparseGmpEigenMatrix& sort_new(const int& dim, const int& type, std::vector < std::vector < IndexType > >& index, std::vector < IndexType >& nbNegatives) const;



    /* --------------------------------------------------
       | Useful functions to deal with complex matrices |
       -------------------------------------------------- */

    /* This function transforms complex matrices into twice as big real matrices */
    SparseGmpEigenMatrix complexIsometry() const;
    /* This function restores the complex matrix corresponding to a big real
       matrix created with the complexIsometry function. */
    SparseGmpEigenMatrix complexIsometryInverse() const;


    /* ------------------------
       | Some class operators |
       ------------------------ */


    /* This function tells whether the matrix is integer */
    inline bool isInt() const
    {
        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            for (Eigen::SparseMatrix<mpfr::mpreal>::InnerIterator it(matrixR,k); it; ++it) {
                if (!mpfr::isint(it.value()))
                    return false;
            }
        }
        if (isComplex) {
            for (IndexType k = 0; k < matrixI.outerSize(); ++k) {
                for (Eigen::SparseMatrix<mpfr::mpreal>::InnerIterator it(matrixI,k); it; ++it) {
                    if (!mpfr::isint(it.value()))
                        return false;
                }
            }
        }
        return true;
    }

    /* This function checks if the imaginary part of the matrix is zero.
       If it is zero, it erases the imaginary table and sets the isComplex
       property to false. If a nonzero imaginary component is detected, the
       imaginary part is kept and isComplex is set to true */
    inline void checkComplexity()
    {
        // First we check if the imaginary matrix is non-empty
        if (matrixI.cols()*matrixI.rows() == 0) {
            isComplex = false;
            return;
        }

        if (matrixI.nonZeros() != 0) {
            isComplex = true;
            return;
        }

        // Then the imaginary part is null so we erase it
        matrixI.resize(0,0);
        isComplex = false;

        return;
    }

    /* This function returns the number of matrix elements */
    inline IndexType numel() const
    {
        return matrixR.cols()*matrixR.rows();
    }

    /* This function returns the number of non-zero elements */
    inline IndexType nnz() const
    {
        if (!isComplex)
            return matrixR.nonZeros();

        IndexType co(0);
        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            Eigen::SparseMatrix<mpfr::mpreal>::InnerIterator itR(matrixR,k);
            Eigen::SparseMatrix<mpfr::mpreal>::InnerIterator itI(matrixI,k);

            while ((itR) || (itI)) {
                if ((itR) && (itI)) {
                    if (itR.row() < itI.row()) {
                        ++co;
                        ++itR;
                    } else if (itR.row() == itI.row()) {
                        ++co;
                        ++itR;
                        ++itI;
                    } else {
                        ++co;
                        ++itI;
                    }
                } else if (itR) {
                    ++co;
                    ++itR;
                } else {
                    ++co;
                    ++itI;
                }
            }
        }

        return co;
    }

public:
    Eigen::SparseMatrix<mpfr::mpreal> matrixR, matrixI;
    bool isComplex;
};


#endif // __SparseGmpEigenMatrix_HPP__
