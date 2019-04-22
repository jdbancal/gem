#ifndef __GmpEigenMatrix_HPP__
#define __GmpEigenMatrix_HPP__
#include "mex.h"

#include <iostream>
#include <vector>
#include <Eigen/MPRealSupport>
#include <Eigen/LU>
#include <Eigen/SVD>
#include "utils.hpp"
#include <Spectra/SymEigsSolver.h>
#include <Spectra/GenEigsSolver.h>
#include <Spectra/SymEigsShiftSolver.h>
#include <Spectra/GenEigsComplexShiftSolver.h>
#include <Spectra/MatOp/DenseSymShiftSolve.h>
#include <Spectra/MatOp/DenseGenComplexShiftSolve.h>

/*
  This file contains the description of our c++ class, including all its
  required procedures. It also includes a function which allows for the creation
  of a class instance from numerical data received from matlab, as well as a few
  output functions tailored for matlab.
*/

/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

class SparseGmpEigenMatrix;

/* The class that we are interfacing to */
class GmpEigenMatrix
{
    friend class SparseGmpEigenMatrix;

public:
    friend GmpEigenMatrix gemRand(const IndexType& m, const IndexType& n);
    friend GmpEigenMatrix& gemRand_new(const IndexType& m, const IndexType& n);

    /* The main constructor */
    GmpEigenMatrix() : isComplex(false) {};
    /* Construction from one mpreal number */
    GmpEigenMatrix(const mpfr::mpreal& value) : matrixR(1,1), isComplex(false) {matrixR(0,0) = value;};
    /* Construction from two mpreal numbers */
    GmpEigenMatrix(const mpfr::mpreal& valueR, const mpfr::mpreal& valueI) : matrixR(1,1), matrixI(1,1), isComplex(valueI != 0) {
        matrixR(0,0) = valueR;
        if (valueI != 0)
            matrixI(0,0) = valueI;
    };
    /* Construction from one double number */
    GmpEigenMatrix(const double& value) : matrixR(1,1), isComplex(false) {matrixR(0,0) = mpfr::mpreal(value);};
    /* Construction from two double numbers */
    GmpEigenMatrix(const double& valueR, const double& valueI) : matrixR(1,1), matrixI(1,1), isComplex(valueI != 0) {
        matrixR(0,0) = mpfr::mpreal(valueR);
        if (valueI != 0)
            matrixI(0,0) = mpfr::mpreal(valueI);
    };
    /* Construction from one int number */
    GmpEigenMatrix(const int& value) : matrixR(1,1), isComplex(false) {matrixR(0,0) = mpfr::mpreal(value);};
    /* Construction from two int numbers */
    GmpEigenMatrix(const int& valueR, const int& valueI) : matrixR(1,1), matrixI(1,1), isComplex(valueI != 0) {
        matrixR(0,0) = mpfr::mpreal(valueR);
        if (valueI != 0)
            matrixI(0,0) = mpfr::mpreal(valueI);
    };
    /* Construction from one mpreal matrix */
    GmpEigenMatrix(const Eigen::Matrix<mpfr::mpreal,Eigen::Dynamic,Eigen::Dynamic>& values) : matrixR(values), isComplex(false) {};
    /* Construction from two mpreal matrices */
    GmpEigenMatrix(const Eigen::Matrix<mpfr::mpreal,Eigen::Dynamic,Eigen::Dynamic>& valuesR, const Eigen::Matrix<mpfr::mpreal,Eigen::Dynamic,Eigen::Dynamic>& valuesI) : matrixR(valuesR), matrixI(valuesI), isComplex(true) {};
    /* The copy constructor */
    GmpEigenMatrix(const GmpEigenMatrix& a) : matrixR(a.matrixR), matrixI(a.matrixI), isComplex(a.isComplex) {};
    /* Construction from a sparse GmpEigenMatrix */
    GmpEigenMatrix(const SparseGmpEigenMatrix& a);
    /* Construction from a matlab struct containing all the class information. */
    GmpEigenMatrix(const mxArray* prhs);
    /* Construction from a matlab table or cell array.*/
    GmpEigenMatrix(const mxArray* prhs, const int& precision);

    /* Destructor */
    virtual ~GmpEigenMatrix() {};


    /* ------------------------------
       | Some procedures for matlab |
       ------------------------------ */

    /* Saving procedure for matlab */
    mxArray* saveobj() const;

    /* Printing of the object in matlab.
       Precision is the maximum number of significative digits that we wish to print.
       Due to alignment, all numbers will not be printed with the same number of
       significative digits. And if some numbers have fewer digits, the missing digits
       are left blank.

       Note : This function implements rounding as provided by the mpfr. Hence, 0.995
              is rounded up to 0.99 rather than 1.00, when cut after two digits (!)
              This is what they call 'correct rounding' or 'round half to even'.
       */
    void display(const int& precision) const;
    /* Another display function, in which each number in the matrix is printed
       independently of the other ones. Here, the width include the sign, digits, dot,
       and exponent if they apply. Each number is processed independently. */
    void displayIndividual(int width) const;

    /* Extracts the first real value of the matrix */
    inline mpfr::mpreal getFirstRealValue()
    {
        return matrixR(0,0);
    }

    /* Extracting the raw data: the following function extracts the internal
       data from the class, and puts it in a basic c-type table to be accessed
       by matlab */
    mxArray* toDouble() const;

    /* This function returns a cell array with strings describing each matrix
       element */
    mxArray* toStrings(const int& precision) const;

    /* This function returns the number of significant digits of the real part in
       matlab format */
    //Eigen::Matrix <double, Eigen::Dynamic, Eigen::Dynamic> precision() const;
    mxArray* precision() const;

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
    GmpEigenMatrix reshape(const IndexType& m2, const IndexType& n2) const;
    GmpEigenMatrix& reshape_new(const IndexType& m2, const IndexType& n2) const;

    /* Concatenates two matrices */
    GmpEigenMatrix vertcat(const GmpEigenMatrix& b) const;
    GmpEigenMatrix& vertcat_new(const GmpEigenMatrix& b) const;
    GmpEigenMatrix horzcat(const GmpEigenMatrix& b) const;
    GmpEigenMatrix& horzcat_new(const GmpEigenMatrix& b) const;

    /* Extracting a sub-matrix */
    GmpEigenMatrix block(const IndexType& i, const IndexType& j, const IndexType& rows, const IndexType& cols) const;
    GmpEigenMatrix subsref(const std::vector<std::vector<IndexType> >& indices) const;
    GmpEigenMatrix& subsref_new(const std::vector<std::vector<IndexType> >& indices) const;
    GmpEigenMatrix subsref(const std::vector<IndexType>& indicesA, const std::vector<IndexType>& indicesB) const;
    GmpEigenMatrix& subsref_new(const std::vector<IndexType>& indicesA, const std::vector<IndexType>& indicesB) const;

    /* Assigning value to a sub-matrix */
    void subsasgn(const std::vector<IndexType>& indices, const GmpEigenMatrix& b);
    void subsasgn(const std::vector<IndexType>& indicesA, const std::vector<IndexType>& indicesB, const GmpEigenMatrix& b);

    /* Diagonal matrix and elements */
    GmpEigenMatrix diagCreate(const IndexType& k) const;
    GmpEigenMatrix& diagCreate_new(const IndexType& k) const;
    GmpEigenMatrix diagExtract(const IndexType& k) const;
    GmpEigenMatrix& diagExtract_new(const IndexType& k) const;


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
      allocated memory is freed properly. */

    // Transpositions b = a.' and b = a'
    GmpEigenMatrix transpose() const;
    GmpEigenMatrix& transpose_new() const;
    GmpEigenMatrix ctranspose() const;
    GmpEigenMatrix& ctranspose_new() const;

    // Element-wise complex conjugation
    GmpEigenMatrix conj() const;
    GmpEigenMatrix& conj_new() const;

    // Real and imaginary parts
    GmpEigenMatrix real() const;
    GmpEigenMatrix& real_new() const;
    GmpEigenMatrix imag() const;
    GmpEigenMatrix& imag_new() const;

    // Various integer parts
    GmpEigenMatrix round() const;
    GmpEigenMatrix& round_new() const;
    GmpEigenMatrix floor() const;
    GmpEigenMatrix& floor_new() const;
    GmpEigenMatrix ceil() const;
    GmpEigenMatrix& ceil_new() const;
    GmpEigenMatrix trunc() const;
    GmpEigenMatrix& trunc_new() const;


    // Addition a += b and c = a+b
    GmpEigenMatrix& operator+=(const GmpEigenMatrix& b);
    GmpEigenMatrix operator+(const GmpEigenMatrix& b) const;
    GmpEigenMatrix& plus_new(const GmpEigenMatrix& b) const;

    // Substraction a -= b and c = a-b
    GmpEigenMatrix& operator-=(const GmpEigenMatrix& b);
    GmpEigenMatrix operator-(const GmpEigenMatrix& b) const;
    GmpEigenMatrix& minus_new(const GmpEigenMatrix& b) const;

    // Unary substraction operators b = -a
    GmpEigenMatrix operator-() const;
    GmpEigenMatrix& uminus_new() const;

    // Element-wise multiplication c = a.*b
    GmpEigenMatrix times(const GmpEigenMatrix& b) const;
    GmpEigenMatrix& times_new(const GmpEigenMatrix& b) const;
    // Element-wise multiplication with a sparse matrix
    SparseGmpEigenMatrix times_fs(const SparseGmpEigenMatrix& b) const;
    SparseGmpEigenMatrix& times_fs_new(const SparseGmpEigenMatrix& b) const;

    // Element-wise division on the right c = a./b
    GmpEigenMatrix rdivide(const GmpEigenMatrix& b) const;
    GmpEigenMatrix& rdivide_new(const GmpEigenMatrix& b) const;

    // Absolute value (or complex magnitude) c = abs(a)
    GmpEigenMatrix abs() const;
    GmpEigenMatrix& abs_new() const;

    // Complex norm b = angle(a)
    GmpEigenMatrix angle() const;
    GmpEigenMatrix& angle_new() const;

    // Element-wise exponential b = exp(a)
    GmpEigenMatrix exp() const;
    GmpEigenMatrix& exp_new() const;

    // Element-wise natural logarithms b = log(a)
    GmpEigenMatrix log() const;
    GmpEigenMatrix& log_new() const;

    // Element-wise power c = a.^b
    GmpEigenMatrix power(const GmpEigenMatrix& b) const;
    GmpEigenMatrix& power_new(const GmpEigenMatrix& b) const;

    // square roots
    GmpEigenMatrix sqrt() const;
    GmpEigenMatrix& sqrt_new() const;

    // matrix multiplication a *= b and c = a*b
    GmpEigenMatrix& operator*=(const GmpEigenMatrix& b);
    GmpEigenMatrix operator*(const GmpEigenMatrix& b) const;
    GmpEigenMatrix& mtimes_new(const GmpEigenMatrix& b) const;
    // multiplication with a sparse matrix c = a*b
    GmpEigenMatrix mtimes_fs(const SparseGmpEigenMatrix& b) const;
    GmpEigenMatrix& mtimes_fs_new(const SparseGmpEigenMatrix& b) const;

    // tensor product c = kron(a,b)
    GmpEigenMatrix kron(const GmpEigenMatrix& b) const;
    GmpEigenMatrix& kron_new(const GmpEigenMatrix& b) const;
    // tensor product with a sparse matrix c = kron(a,b)
    SparseGmpEigenMatrix kron_fs(const SparseGmpEigenMatrix& b) const;
    SparseGmpEigenMatrix& kron_fs_new(const SparseGmpEigenMatrix& b) const;


    /* ----------------------------------------------------------
       | Some mathematical functions defined for real arguments |
       |       (the real test is made on the matlab side.       |
       |       Here we simply ignore the imaginary parts)       |
       ---------------------------------------------------------- */

    // cubic roots
    GmpEigenMatrix cbrt() const;
    GmpEigenMatrix& cbrt_new() const;

    // trigonometric functions
    GmpEigenMatrix sin() const;
    GmpEigenMatrix& sin_new() const;
    GmpEigenMatrix cos() const;
    GmpEigenMatrix& cos_new() const;
    GmpEigenMatrix tan() const;
    GmpEigenMatrix& tan_new() const;

    GmpEigenMatrix sec() const;
    GmpEigenMatrix& sec_new() const;
    GmpEigenMatrix csc() const;
    GmpEigenMatrix& csc_new() const;
    GmpEigenMatrix cot() const;
    GmpEigenMatrix& cot_new() const;

    GmpEigenMatrix asin() const;
    GmpEigenMatrix& asin_new() const;
    GmpEigenMatrix acos() const;
    GmpEigenMatrix& acos_new() const;
    GmpEigenMatrix atan() const;
    GmpEigenMatrix& atan_new() const;

    // hypergeometric functions
    GmpEigenMatrix sinh() const;
    GmpEigenMatrix& sinh_new() const;
    GmpEigenMatrix cosh() const;
    GmpEigenMatrix& cosh_new() const;
    GmpEigenMatrix tanh() const;
    GmpEigenMatrix& tanh_new() const;



    /* ---------------------------------------
       |   Some linear algebra operations    |
       --------------------------------------- */

    IndexType rank() const;
    GmpEigenMatrix mldivide(const GmpEigenMatrix& b) const;
    GmpEigenMatrix& mldivide_new(const GmpEigenMatrix& b) const;
    GmpEigenMatrix inv() const;
    GmpEigenMatrix& inv_new() const;
    GmpEigenMatrix eig(GmpEigenMatrix& V) const;
    GmpEigenMatrix& eig_new(GmpEigenMatrix& V) const;
    GmpEigenMatrix eigs(const long int& nbEigenvalues, GmpEigenMatrix& V, const long int& type, const GmpEigenMatrix& sigma) const;
    GmpEigenMatrix& eigs_new(const long int& nbEigenvalues, GmpEigenMatrix& V, const long int& type, const GmpEigenMatrix& sigma) const;
    GmpEigenMatrix svd(GmpEigenMatrix& U, GmpEigenMatrix& V) const;
    GmpEigenMatrix& svd_new(GmpEigenMatrix& U, GmpEigenMatrix& V) const;




    /* ---------------------------------------
       | Some tests and comparison operators |
       --------------------------------------- */

    // <, <=, >, >=, ==, !=, isequal tests
    Eigen::Matrix <bool, Eigen::Dynamic, Eigen::Dynamic> operator<(const GmpEigenMatrix& b) const;
    Eigen::Matrix <bool, Eigen::Dynamic, Eigen::Dynamic> operator<=(const GmpEigenMatrix& b) const;
    Eigen::Matrix <bool, Eigen::Dynamic, Eigen::Dynamic> operator>(const GmpEigenMatrix& b) const;
    Eigen::Matrix <bool, Eigen::Dynamic, Eigen::Dynamic> operator>=(const GmpEigenMatrix& b) const;
    Eigen::Matrix <bool, Eigen::Dynamic, Eigen::Dynamic> eq(const GmpEigenMatrix& b) const;
    Eigen::Matrix <bool, Eigen::Dynamic, Eigen::Dynamic> ne(const GmpEigenMatrix& b) const;
    Eigen::Matrix <bool, Eigen::Dynamic, Eigen::Dynamic> isnan() const;
    Eigen::Matrix <bool, Eigen::Dynamic, Eigen::Dynamic> isinf() const;
    bool identicalValues(const GmpEigenMatrix& b) const;
    bool identicalValuesNaNok(const GmpEigenMatrix& b) const;

    // isreal
    inline bool isreal() const { return (!isComplex); }

    // symmetry tests
    bool issymmetric() const;
    bool ishermitian() const;

    // column-wise minimum b = min(a)
    GmpEigenMatrix colMin(std::vector<IndexType>& indices) const;
    inline GmpEigenMatrix colMin() const
    {
        std::vector<IndexType> indices;
        return colMin(indices);
    }
    GmpEigenMatrix& colMin_new(std::vector<IndexType>& indices) const;
    inline GmpEigenMatrix& colMin_new() const
    {
        std::vector<IndexType> indices;
        return colMin_new(indices);
    }

    // line-wise minimum b = min(a,[],2)
    GmpEigenMatrix rowMin(std::vector<IndexType>& indices) const;
    inline GmpEigenMatrix rowMin() const
    {
        std::vector<IndexType> indices;
        return rowMin(indices);
    }
    GmpEigenMatrix& rowMin_new(std::vector<IndexType>& indices) const;
    inline GmpEigenMatrix& rowMin_new() const
    {
        std::vector<IndexType> indices;
        return rowMin_new(indices);
    }

    // element-wise minimum c = min(a, b)
    GmpEigenMatrix ewMin(const GmpEigenMatrix& b) const;
    GmpEigenMatrix& ewMin_new(const GmpEigenMatrix& b) const;

    // column-wise maximum b = max(a)
    GmpEigenMatrix colMax(std::vector<IndexType>& indices) const;
    inline GmpEigenMatrix colMax() const
    {
        std::vector<IndexType> indices;
        return colMax(indices);
    }
    GmpEigenMatrix& colMax_new(std::vector<IndexType>& indices) const;
    inline GmpEigenMatrix& colMax_new() const
    {
        std::vector<IndexType> indices;
        return colMax_new(indices);
    }

    // line-wise maximum b = max(a,[],2)
    GmpEigenMatrix rowMax(std::vector<IndexType>& indices) const;
    inline GmpEigenMatrix rowMax() const
    {
        std::vector<IndexType> indices;
        return rowMax(indices);
    }
    GmpEigenMatrix& rowMax_new(std::vector<IndexType>& indices) const;
    inline GmpEigenMatrix& rowMax_new() const
    {
        std::vector<IndexType> indices;
        return rowMax_new(indices);
    }

    // element-wise maximum c = max(a, b)
    GmpEigenMatrix ewMax(const GmpEigenMatrix& b) const;
    GmpEigenMatrix& ewMax_new(const GmpEigenMatrix& b) const;


    // sums and products
    GmpEigenMatrix colSum() const;
    GmpEigenMatrix& colSum_new() const;
    GmpEigenMatrix rowSum() const;
    GmpEigenMatrix& rowSum_new() const;
    GmpEigenMatrix colProd() const;
    GmpEigenMatrix& colProd_new() const;
    GmpEigenMatrix rowProd() const;
    GmpEigenMatrix& rowProd_new() const;

    // sorting functions
    GmpEigenMatrix sort(const int& dim, const int& type, std::vector < std::vector < IndexType > >& index) const;
    GmpEigenMatrix& sort_new(const int& dim, const int& type, std::vector < std::vector < IndexType > >& index) const;
    std::vector < IndexType > sortrowsc(const std::vector < int > ascending) const;

    /* -------------------------------
       | Some mathematical constants |
       ------------------------------- */

    friend GmpEigenMatrix constLog2();
    friend GmpEigenMatrix& constLog2_new();
    friend GmpEigenMatrix constPi();
    friend GmpEigenMatrix& constPi_new();
    friend GmpEigenMatrix constEuler();
    friend GmpEigenMatrix& constEuler_new();
    friend GmpEigenMatrix constCatalan();
    friend GmpEigenMatrix& constCatalan_new();






    /* --------------------------------------------------
       | Useful functions to deal with complex matrices |
       -------------------------------------------------- */

    /* This function transforms complex matrices into twice as big real matrices */
    GmpEigenMatrix complexIsometry() const;
    /* This function restores the complex matrix corresponding to a big real
       matrix created with the complexIsometry function. */
    GmpEigenMatrix complexIsometryInverse() const;







    /* ------------------------
       | Some class operators |
       ------------------------ */


    /* This function tells whether the matrix is integer */
    inline bool isInt() const
    {
        for (IndexType j = 0; j < matrixR.cols(); ++j) {
            for (IndexType i = 0; i < matrixR.rows(); ++i) {
                if ((!mpfr::isint(matrixR(i,j))) || ((isComplex) && (!mpfr::isint(matrixI(i,j)))))
                    return false;
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

        for (IndexType j = 0; j < matrixI.cols(); ++j) {
            for (IndexType i = 0; i < matrixI.rows(); ++i) {
                if (!mpfr::iszero(matrixI(i,j))) {
                    //Then the imaginary part is not null so we keep it and adjust the complex property accordingly
                    isComplex = true;
                    return;
                }
            }
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

public: // Actually, friend classes like sgem.hpp should have access to these variables... so we make them public for now.
    Eigen::Matrix<mpfr::mpreal,Eigen::Dynamic,Eigen::Dynamic> matrixR, matrixI;
    bool isComplex;
};


/* ---------------------------------
   | Now we include the definition |
   |    of the friend functions    |
   --------------------------------- */


/* -------------------------------
   | Some mathematical constants |
   ------------------------------- */
// The imaginary unit constant 1i
inline GmpEigenMatrix constI()
{
    GmpEigenMatrix result(mpfr::mpreal(0));

    result.isComplex = true;
    result.matrixI.resize(1,1);
    result.matrixI(0,0) = 1;

    return result;
}

// The imaginary unit constant 1i
inline GmpEigenMatrix& constI_new()
{
    GmpEigenMatrix& result(*(new GmpEigenMatrix));

    result.isComplex = true;
    result.matrixR.resize(1,1);
    result.matrixR(0,0) = 0;
    result.matrixI.resize(1,1);
    result.matrixI(0,0) = 1;

    return result;
}

// The log(2) constant
inline GmpEigenMatrix constLog2()
{
   return GmpEigenMatrix(mpfr::const_log2());
}

// The log(2) constant
inline GmpEigenMatrix& constLog2_new()
{
   return *(new GmpEigenMatrix(mpfr::const_log2()));
}

// The pi constant
inline GmpEigenMatrix constPi()
{
    return GmpEigenMatrix(mpfr::const_pi());
}

// The pi constant
inline GmpEigenMatrix& constPi_new()
{
    return *(new GmpEigenMatrix(mpfr::const_pi()));
}

// The Euler constant
inline GmpEigenMatrix constEuler()
{
    return GmpEigenMatrix(mpfr::const_euler());
}

// The Euler constant
inline GmpEigenMatrix& constEuler_new()
{
    return *(new GmpEigenMatrix(mpfr::const_euler()));
}

// The Catalan constant
inline GmpEigenMatrix constCatalan()
{
    return GmpEigenMatrix(mpfr::const_catalan());
}

// The Catalan constant
inline GmpEigenMatrix& constCatalan_new()
{
    return *(new GmpEigenMatrix(mpfr::const_catalan()));
}

// Random matrix
GmpEigenMatrix gemRand(const IndexType& m, const IndexType& n);
GmpEigenMatrix& gemRand_new(const IndexType& m, const IndexType& n);





#endif // __GmpEigenMatrix_HPP__
