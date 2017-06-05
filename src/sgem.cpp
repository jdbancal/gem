#include "gem.hpp"
#include "sgem.hpp"
#include <unsupported/Eigen/MatrixFunctions>
#include <unsupported/Eigen/KroneckerProduct>

/*
  This file contains the implementation of our c++ class, including all its
  required procedures. It also includes a function which allows for the creation
  of a class instance from numerical data received from matlab, as well as a few
  output functions tailored for matlab.
*/

/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


using namespace std;
using namespace mpfr;
using namespace Eigen;


/* Construction from a dense GmpEigenMatrix */
SparseGmpEigenMatrix::SparseGmpEigenMatrix(const GmpEigenMatrix& a, const mpfr::mpreal& threshold) : matrixR(a.matrixR.sparseView(threshold, mpreal(1))), matrixI(a.matrixI.sparseView(threshold, mpreal(1))), isComplex(a.isComplex) {}


/* Construcion from a matlab struct containing all the class information. */
SparseGmpEigenMatrix::SparseGmpEigenMatrix(const mxArray* prhs) {
    if ((!mxIsStruct(prhs)) || ((mxGetNumberOfElements(prhs) != 1)))
        mexErrMsgTxt("The argument is not a struct array with one element.");

    // The structure should at least contain the following fields :
    // precisionR and matrixR
    const mxArray* sizeField = mxGetField(prhs, 0, "size");
    if (sizeField == NULL)
        mexErrMsgTxt("The field 'size' was not found in the provided structure.");
    if (!mxIsDouble(sizeField))
        mexErrMsgTxt("The size field should be a double array.");
    const mxArray* precisionField = mxGetField(prhs, 0, "precisionR");
    if (precisionField == NULL)
        mexErrMsgTxt("The field 'precisionR' was not found in the provided structure.");
    if (!mxIsDouble(precisionField))
        mexErrMsgTxt("The precisionR field should be a double array.");
    const mxArray* field = mxGetField(prhs, 0, "matrixR");
    if (field == NULL)
        mexErrMsgTxt("The field 'matrixR' was not found in the provided structure.");
    if (!mxIsCell(field))
        mexErrMsgTxt("The matrixR field should be a cell array.");

    if ((mxGetM(sizeField) != 1) && (mxGetN(sizeField) != 2))
        mexErrMsgTxt("The size field should have size 1x2.");
    double* sizePtr(mxGetPr(sizeField));
    mwSize m(sizePtr[0]);
    mwSize n(sizePtr[1]);
    mwSize nbElementsR(mxGetM(field));
    if (mxGetN(field) != 3)
        mexErrMsgTxt("The matrixR cell array should have size nx3.");
    if (mxGetN(precisionField) != 1)
        mexErrMsgTxt("PrecisionR should be a vertical vector.");
    if (mxGetM(precisionField) != nbElementsR)
        mexErrMsgTxt("Different fields have incompatible sizes.");

    // We set the size
    matrixR.resize(m,n);
//    matrixR.reserve(nbElementsR);

    // We first copy the real data to a triplet
    vector< Triplet<mpreal> > tripletList;
    tripletList.reserve(nbElementsR);
    double* precisionPtr(mxGetPr(precisionField));
    for (IndexType index = 0; index < nbElementsR; ++index) {
        const mxArray* iPtr(mxGetCell(field, index));
        const mxArray* jPtr(mxGetCell(field, nbElementsR+index));
        const mxArray* cellElementPtr(mxGetCell(field, 2*nbElementsR+index));
        if ((iPtr == NULL) || (jPtr == NULL) || (cellElementPtr == NULL))
            mexErrMsgTxt("Cell element is empty.");
        else {
            if (!(mxGetClassID(iPtr) == mxDOUBLE_CLASS))
                mexErrMsgTxt("Cell element is not an mxDOUBLE integer.");
            if (!(mxGetClassID(jPtr) == mxDOUBLE_CLASS))
                mexErrMsgTxt("Cell element is not an mxDOUBLE integer.");
            if (!mxIsChar(cellElementPtr))
                mexErrMsgTxt("Cell element is not a string.");

            // We get the numerical value of the i and j indices
            double* iPtr2(mxGetPr(iPtr));
            double* jPtr2(mxGetPr(jPtr));

            mwSize stringLength = mxGetNumberOfElements(cellElementPtr);
            if (stringLength == 0)
                mexErrMsgTxt("String is of length 0.");
            else {
                char c_string[stringLength+1];
                if (mxGetString(cellElementPtr, c_string, stringLength+1) != 0)
                    mexErrMsgTxt("Error in extracting a string from a cell array");
                tripletList.push_back(Triplet<mpreal>((*iPtr2)-1, (*jPtr2)-1, mpreal(string(c_string), precisionPtr[index])));
            }
        }
    }
    // Now we can assign the data to the sparse matrix
    matrixR.setFromTriplets(tripletList.begin(), tripletList.end());

    // We compress the output if possible
    matrixR.prune(0,0);
    matrixR.makeCompressed();


    // We check whether the matrix also has an imaginary part
    precisionField = mxGetField(prhs, 0, "precisionI");
    field = mxGetField(prhs, 0, "matrixI");
    if ((precisionField != NULL) && (field != NULL))
    {
        isComplex = true;

        if (!mxIsDouble(precisionField))
            mexErrMsgTxt("The precisionI field should be a double array.");
        if (!mxIsCell(field))
            mexErrMsgTxt("The matrixI field should be a cell array.");

        mwSize nbElementsI(mxGetM(field));
        if (mxGetN(field) != 3)
            mexErrMsgTxt("The matrixI cell array should have size nx3.");
        if (mxGetN(precisionField) != 1)
            mexErrMsgTxt("PrecisionI should be a vertical vector.");
        if (mxGetM(precisionField) != nbElementsI)
            mexErrMsgTxt("Different fields have incompatible sizes.");

        // We set the size
        matrixI.resize(m,n);

        // We first copy the real data to a triplet
        vector< Triplet<mpreal> > tripletList;
        tripletList.reserve(nbElementsI);
        double* precisionPtr(mxGetPr(precisionField));
        for (IndexType index = 0; index < nbElementsI; ++index) {
            const mxArray* iPtr(mxGetCell(field, index));
            const mxArray* jPtr(mxGetCell(field, nbElementsI+index));
            const mxArray* cellElementPtr(mxGetCell(field, 2*nbElementsI+index));
            if ((iPtr == NULL) || (jPtr == NULL) || (cellElementPtr == NULL))
                mexErrMsgTxt("Cell element is empty.");
            else {
                if (!(mxGetClassID(iPtr) == mxDOUBLE_CLASS))
                    mexErrMsgTxt("Cell element is not an mxDOUBLE integer.");
                if (!(mxGetClassID(jPtr) == mxDOUBLE_CLASS))
                    mexErrMsgTxt("Cell element is not an mxDOUBLE integer.");
                if (!mxIsChar(cellElementPtr))
                    mexErrMsgTxt("Cell element is not a string.");

                // We get the numerical value of the i and j indices
                double* iPtr2(mxGetPr(iPtr));
                double* jPtr2(mxGetPr(jPtr));

                mwSize stringLength = mxGetNumberOfElements(cellElementPtr);
                if (stringLength == 0)
                    mexErrMsgTxt("String is of length 0.");
                else {
                    char c_string[stringLength+1];
                    if (mxGetString(cellElementPtr, c_string, stringLength+1) != 0)
                        mexErrMsgTxt("Error in extracting a string from a cell array");
                    tripletList.push_back(Triplet<mpreal>((*iPtr2)-1, (*jPtr2)-1, mpreal(string(c_string), precisionPtr[index])));
                }
            }
        }
        // Now we can assign the data to the sparse matrix
        matrixI.setFromTriplets(tripletList.begin(), tripletList.end());

        // We compress the output if possible
        matrixI.prune(0,0);
        matrixI.makeCompressed();
        checkComplexity();
    } else
        isComplex = false;
}


/* Construction from a matlab sparse table. For double element, this
   constructor takes exactly 15 digits from them and sets all other digits
   to 0. */
SparseGmpEigenMatrix::SparseGmpEigenMatrix(const mxArray* rows, const mxArray* cols, const mxArray* values, const IndexType& m, const IndexType& n, const int& precision) {

    // We set the required precision
    mp_prec_t precisionInBits(mpfr::digits2bits(precision));

    /* Get the size and pointers to input data */
    double* prows(mxGetPr(rows));
    double* pcols(mxGetPr(cols));
    double* pr(mxGetPr(values));
    double* pi(mxGetPi(values));
    isComplex = (pi==NULL ? 0 : 1);

    // If the input is real, we set the imaginary data to nothing
    if (isComplex == 0)
        matrixI = SparseMatrix<mpfr::mpreal>();

    // We set the size
    matrixR.resize(m,n);
    if (isComplex)
        matrixI.resize(m,n);

    // The number of nonzero elements in the matrix
    mwSize nbElementsR(mxGetM(rows));

    // Now we copy the data over. First to triplet
    vector< Triplet<mpreal> > tripletList;
    tripletList.reserve(nbElementsR);
    for (IndexType index = 0; index < nbElementsR; ++index) {
        if ((!isComplex) || (pr[index] != 0)){
            tripletList.push_back(Triplet<mpreal>(prows[index]-1, pcols[index]-1, mpreal(toString15(pr[index]), precisionInBits)));
        }
    }

    // now into the sparse matrix
    matrixR.setFromTriplets(tripletList.begin(), tripletList.end());

    // This should not be needed if the , but just in case the struct array was
    // provided by a user, we need to make sure that there are no zeros in the
    // matrix. So we compress the output if possible
    matrixR.prune(0,0);
    matrixR.makeCompressed();

    if (isComplex) {
        // Now we copy the data over. First to triplet
        vector< Triplet<mpreal> > tripletList;
        tripletList.reserve(nbElementsR);
        for (IndexType index = 0; index < nbElementsR; ++index) {
            if (pi[index] != 0)
                tripletList.push_back(Triplet<mpreal>(prows[index]-1, pcols[index]-1, mpreal(toString15(pi[index]), precisionInBits)));
        }
        // now into the sparse matrix
        matrixI.setFromTriplets(tripletList.begin(), tripletList.end());

        // We compress the output if possible
        matrixI.prune(0,0);
        matrixI.makeCompressed();
        checkComplexity();
    }
}


/* Construction from a list of values in gem format */
SparseGmpEigenMatrix::SparseGmpEigenMatrix(const mxArray* rows, const mxArray* cols, const GmpEigenMatrix& values, const IndexType& m, const IndexType& n, const int& precision) {

    // We set the required precision
    mp_prec_t precisionInBits(mpfr::digits2bits(precision));

    /* Get the size and pointers to input data */
    double* prows(mxGetPr(rows));
    double* pcols(mxGetPr(cols));
    isComplex = values.isComplex;

    // If the input is real, we set the imaginary data to nothing
    if (isComplex == 0)
        matrixI = SparseMatrix<mpfr::mpreal>();

    // We set the size
    matrixR.resize(m,n);
    if (isComplex)
        matrixI.resize(m,n);

    // The number of nonzero elements in the matrix
    mwSize nbElementsR(mxGetM(rows));

    // Now we copy the data over. First to triplet
    vector< Triplet<mpreal> > tripletList;
    tripletList.reserve(nbElementsR);
    for (IndexType index = 0; index < nbElementsR; ++index) {
        if ((!isComplex) || (values.matrixR(index,0) != 0)){
            tripletList.push_back(Triplet<mpreal>(prows[index]-1, pcols[index]-1, values.matrixR(index,0)));
        }
    }

    // now into the sparse matrix
    matrixR.setFromTriplets(tripletList.begin(), tripletList.end());

    // This should not be needed if the , but just in case the struct array was
    // provided by a user, we need to make sure that there are no zeros in the
    // matrix. So we compress the output if possible
    matrixR.prune(0,0);
    matrixR.makeCompressed();

    if (isComplex) {
        // Now we copy the data over. First to triplet
        vector< Triplet<mpreal> > tripletList;
        tripletList.reserve(nbElementsR);
        for (IndexType index = 0; index < nbElementsR; ++index) {
            if (values.matrixI(index,0) != 0)
                tripletList.push_back(Triplet<mpreal>(prows[index]-1, pcols[index]-1, values.matrixI(index,0)));
        }
        // now into the sparse matrix
        matrixI.setFromTriplets(tripletList.begin(), tripletList.end());

        // We compress the output if possible
        matrixI.prune(0,0);
        matrixI.makeCompressed();
        checkComplexity();
    }
}


/* Saving procedure
   This function returns a structure to matlab which contains all the data
   that define the SparseGemEigenMatrix instance
   NOTE : This function returns indices in matlab format */
mxArray* SparseGmpEigenMatrix::saveobj() const
{
    // All possible fields of the structure we will return
    const char *field_names[] = {"size", "precisionR", "matrixR", "precisionI", "matrixI"};

    // We create the matlab structure to be populated
    mxArray* ptr;
    if (isComplex)
        ptr = mxCreateStructMatrix(1, 1, 5, field_names);
    else
        ptr = mxCreateStructMatrix(1, 1, 3, field_names);

    // Now, we add each field
    // We start with the matrix size
    mwSize m(matrixR.rows());
    mwSize n(matrixR.cols());
    mxArray* sizeField(mxCreateNumericMatrix(1, 2, mxDOUBLE_CLASS, mxREAL));
    double* pointerSize(mxGetPr(sizeField));
    pointerSize[0] = double(m);
    pointerSize[1] = double(n);

    // We save the size
    mxSetField(ptr, 0, "size", sizeField);

    // The number of elements we need to save
    IndexType nnzR(matrixR.nonZeros());

    // Now we save the real part
    mxArray* precisionField(mxCreateNumericMatrix(nnzR, 1, mxDOUBLE_CLASS, mxREAL));
    mxArray* field(mxCreateCellMatrix(nnzR, 3));

    // Now we iterate on all the elements of the cell array
    mwIndex index(0);
    double* pointerR(mxGetPr(precisionField));
    for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
        for (SparseMatrix<mpreal>::InnerIterator it(matrixR,k); it; ++it) {
            pointerR[index] = double(it.value().get_prec());
            mxSetCell(field, index, mxCreateDoubleScalar(it.row()+1));
            mxSetCell(field, nnzR + index, mxCreateDoubleScalar(it.col()+1));
            mxSetCell(field, 2*nnzR + index, mxCreateString(it.value().toString().c_str()));
            ++index;
        }
    }

    // We save these objects into the matlab structure
    mxSetField(ptr, 0, "precisionR", precisionField);
    mxSetField(ptr, 0, "matrixR", field);


    // Now we also add the imaginary part is there is one
    if (isComplex)
    {
        // The number of elements we need to save
        IndexType nnzI(matrixI.nonZeros());

        precisionField = mxCreateNumericMatrix(nnzI, 1, mxDOUBLE_CLASS, mxREAL);
        field = mxCreateCellMatrix(nnzI, 3);

        // Now we iterate on all the elements of the cell array
        mwIndex index(0);
        double* pointerI(mxGetPr(precisionField));
        for (IndexType k = 0; k < matrixI.outerSize(); ++k) {
            for (SparseMatrix<mpreal>::InnerIterator it(matrixI,k); it; ++it) {
                pointerI[index] = double(it.value().get_prec());
                mxSetCell(field, index, mxCreateDoubleScalar(it.row()+1));
                mxSetCell(field, nnzI + index, mxCreateDoubleScalar(it.col()+1));
                mxSetCell(field, 2*nnzI + index, mxCreateString(it.value().toString().c_str()));
                ++index;
            }
        }

        // We save these objects into the matlab structure
        mxSetField(ptr, 0, "precisionI", precisionField);
        mxSetField(ptr, 0, "matrixI", field);
    }

    return ptr;
}




/*
 Here is another display function, in which each number in the matrix is printed
 independently of the other ones. Here, the width include the sign, digits, dot,
 and exponent if they apply. Each number is processed independently.
*/
void SparseGmpEigenMatrix::displayIndividual(int width) const
{
    // If the matrix is empty, we print nothing
    if (matrixR.rows()*matrixR.cols() == 0)
    {
        mexPrintf("  []\n");
        return;
    }

    // If we have a matrix of zeros
    if ((matrixR.nonZeros() == 0) && (matrixI.nonZeros() == 0))
    {
        mexPrintf("  All zero sparse %ix%i matrix\n", matrixR.rows(), matrixR.cols());
        return;
    }


    if (!isComplex) {
        // This case is simpler
        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            for (SparseMatrix<mpreal>::InnerIterator it(matrixR,k); it; ++it) {
                mexPrintf("  (%i,%i)  %s\n", it.row()+1, it.col()+1, it.value().toString(std::min(bits2digits(it.value().get_prec()),width)).c_str());
            }
        }
    } else {
        // For each column, we should merge the lists of lines with nonzero elements...
        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
            SparseMatrix<mpreal>::InnerIterator itI(matrixI,k);

            while ((itR) || (itI)) {
                if ((itR) && (itI)) {
                    if (itR.row() < itI.row()) {
                        mexPrintf("  (%i,%i)  %s\n", itR.row()+1, itR.col()+1, itR.value().toString(std::min(bits2digits(itR.value().get_prec()),width)).c_str());
                        ++itR;
                    } else if (itR.row() == itI.row()) {
                        mexPrintf("  (%i,%i)  %s", itR.row()+1, itR.col()+1, itR.value().toString(std::min(bits2digits(itR.value().get_prec()),width)).c_str());
                        string chain(itI.value().toString(std::min(bits2digits(itI.value().get_prec()),width))+"i");
                        if (itI.value() < 0)
                            chain = " - " + chain.substr(1,chain.length()-1);
                        else
                            chain = " + " + chain;
                        mexPrintf("%s\n", chain.c_str());
                        ++itR;
                        ++itI;
                    } else {
                        mexPrintf("  (%i,%i)  %si\n", itI.row()+1, itI.col()+1, itI.value().toString(std::min(bits2digits(itI.value().get_prec()),width)).c_str());
                        ++itI;
                    }
                } else if (itR) {
                    mexPrintf("  (%i,%i)  %s\n", itR.row()+1, itR.col()+1, itR.value().toString(std::min(bits2digits(itR.value().get_prec()),width)).c_str());
                    ++itR;
                } else {
                    mexPrintf("  (%i,%i)  %si\n", itI.row()+1, itI.col()+1, itI.value().toString(std::min(bits2digits(itI.value().get_prec()),width)).c_str());
                    ++itI;
                }
            }
        }
    }
}




/* Extracting the raw data: the following function extracts the internal
   data from the class, and puts it in a basic c-type table to be accessed
   by matlab. The result is a sparse array in matlab. */
mxArray* SparseGmpEigenMatrix::toDouble() const
{
    // Let's allocate memory for the matrix
    mxArray* plhs;
    if (isComplex)
        plhs = mxCreateSparse(matrixR.rows(), matrixR.cols(), matrixR.nonZeros()+matrixI.nonZeros(), mxCOMPLEX);
    else
        plhs = mxCreateSparse(matrixR.rows(), matrixR.cols(), matrixR.nonZeros(), mxREAL);

    // We get the pointers to where we need to input data
    double *si,*sr;
    mwIndex *irs,*jcs;
    sr  = mxGetPr(plhs);
    si  = mxGetPi(plhs);
    irs = mxGetIr(plhs);
    jcs = mxGetJc(plhs);

    // We copy over the data
    if (!isComplex) {
        // This case is simpler
        jcs[0] = 0;
        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            IndexType index(0);
            for (SparseMatrix<mpreal>::InnerIterator it(matrixR,k); it; ++it) {
                irs[index] = it.row();
                sr[index] = it.value().toDouble();
                ++index;
            }
            jcs[k+1] = jcs[k] + index;
            irs += index;
            sr += index;
        }
    } else {
        // For each column, we should merge the lists of lines with nonzero elements...
        jcs[0] = 0;
        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
            SparseMatrix<mpreal>::InnerIterator itI(matrixI,k);

            IndexType index(0);
            while ((itR) || (itI)) {
                if ((itR) && (itI)) {
                    if (itR.row() < itI.row()) {
                        irs[index] = itR.row();
                        sr[index] = itR.value().toDouble();
                        si[index] = 0;
                        ++index;
                        ++itR;
                    } else if (itR.row() == itI.row()) {
                        irs[index] = itR.row();
                        sr[index] = itR.value().toDouble();
                        si[index] = itI.value().toDouble();
                        ++index;
                        ++itR;
                        ++itI;
                    } else {
                        irs[index] = itI.row();
                        sr[index] = 0;
                        si[index] = itI.value().toDouble();
                        ++index;
                        ++itI;
                    }
                } else if (itR) {
                    irs[index] = itR.row();
                    sr[index] = itR.value().toDouble();
                    si[index] = 0;
                    ++index;
                    ++itR;
                } else {
                    irs[index] = itI.row();
                    sr[index] = 0;
                    si[index] = itI.value().toDouble();
                    ++index;
                    ++itI;
                }
            }
            jcs[k+1] = jcs[k] + index;
            irs += index;
            sr += index;
            si += index;
        }
    }

    return plhs;
}


/* This function returns a cell array with strings describing each matrix
   element

   NOTE : Exceptionally, this function returns indices in matlab format
     (so that the matlab interface doesn't need to do it) */
mxArray* SparseGmpEigenMatrix::toStrings(const int& precision) const
{
    // We start by filling dynamically allocated tables (we don't know exactly
    // how many elements we'll end up with)
    vector<double> rows, cols;
    vector<string> values;
    rows.reserve(matrixR.nonZeros());
    cols.reserve(matrixR.nonZeros());
    values.reserve(matrixR.nonZeros());

    if (!isComplex) {
        // This case is simpler
        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            for (SparseMatrix<mpreal>::InnerIterator it(matrixR,k); it; ++it) {
                rows.push_back(it.row()+1);
                cols.push_back(it.col()+1);
                values.push_back(it.value().toString(precision));
            }
        }
    } else {
        // For each column, we should merge the lists of lines with nonzero elements...
        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
            SparseMatrix<mpreal>::InnerIterator itI(matrixI,k);

            while ((itR) || (itI)) {
                if ((itR) && (itI)) {
                    if (itR.row() < itI.row()) {
                        rows.push_back(itR.row()+1);
                        cols.push_back(itR.col()+1);
                        values.push_back(itR.value().toString(precision));
                        ++itR;
                    } else if (itR.row() == itI.row()) {
                        rows.push_back(itR.row()+1);
                        cols.push_back(itR.col()+1);
                        string chain(itI.value().toString(precision)+"i");
                        if (itI.value() > 0)
                            chain = "+" + chain;
                        values.push_back(itR.value().toString(precision) + chain);
                        ++itR;
                        ++itI;
                    } else {
                        rows.push_back(itI.row()+1);
                        cols.push_back(itI.col()+1);
                        values.push_back(itI.value().toString(precision)+"i");
                        ++itI;
                    }
                } else if (itR) {
                    rows.push_back(itR.row()+1);
                    cols.push_back(itR.col()+1);
                    values.push_back(itR.value().toString(precision));
                    ++itR;
                } else {
                    rows.push_back(itI.row()+1);
                    cols.push_back(itI.col()+1);
                    values.push_back(itI.value().toString(precision)+"i");
                    ++itI;
                }
            }
        }
    }

    // We send these strings to matlab
    mxArray* cellArray(mxCreateCellMatrix(values.size(), 3));
    for (IndexType i(0); i < values.size(); ++i) {
        mxSetCell(cellArray, i, mxCreateDoubleScalar(rows[i]));
        mxSetCell(cellArray, values.size()+i, mxCreateDoubleScalar(cols[i]));
        mxSetCell(cellArray, 2*values.size()+i, mxCreateString(values[i].c_str()));
    }

    return cellArray;
}


/* This function returns the number of significant digits in a table */
/*Matrix <double, Dynamic, Dynamic> GmpEigenMatrix::precision() const
{
    // Get the size and pointers to input data
    IndexType m(matrixR.rows());
    IndexType n(matrixR.cols());

    Matrix<double, Dynamic, Dynamic> table;
    table.resize(m,n);

    // Now we copy the precision to the table
    for (IndexType j = 0; j < n; ++j)
        for (IndexType i = 0; i < m; ++i)
            table(i,j) = mpfr::bits2digits(matrixR(i,j).get_prec());

    return table;
}*/

/* This function returns the number of significant digits of the real part in
   matlab format */
mxArray* SparseGmpEigenMatrix::precision() const
{
    // Create a sparse vector for all elements to be returned
    mxArray* plhs;
    plhs = mxCreateSparse(matrixR.rows(), matrixR.cols(), matrixR.nonZeros(), mxREAL);

    // We get the pointers to where we need to input data
    double *sr;
    mwIndex *irs,*jcs;
    sr  = mxGetPr(plhs);
    irs = mxGetIr(plhs);
    jcs = mxGetJc(plhs);

    // Now we copy the precision to the matlab table
    jcs[0] = 0;
    for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
        IndexType index(0);
        for (SparseMatrix<mpreal>::InnerIterator it(matrixR,k); it; ++it) {
            irs[index] = it.row();
            sr[index] = mpfr::bits2digits(it.value().get_prec());
            ++index;
        }
        jcs[k+1] = jcs[k] + index;
        irs += index;
        sr += index;
    }

    return plhs;
}



/* This function removes all components in the matrix whose magnitude is smaller
   than the provided tolerance */
SparseGmpEigenMatrix SparseGmpEigenMatrix::clean(const GmpEigenMatrix& tolerance) const
{
    SparseGmpEigenMatrix result(*this);

    // We extract the tolerance
    mpreal tol(tolerance.matrixR(0,0));

    for (IndexType k = 0; k < result.matrixR.outerSize(); ++k)
        for (SparseMatrix<mpreal>::InnerIterator it(result.matrixR,k); it; ++it)
            if (mpfr::abs(it.value()) < tol)
                it.valueRef() = 0;

    // We may have set some coefficients to zero...
    result.matrixR.prune(0, 0);
    result.matrixR.makeCompressed();

    if (isComplex) {
        for (IndexType k = 0; k < result.matrixI.outerSize(); ++k)
            for (SparseMatrix<mpreal>::InnerIterator it(result.matrixI,k); it; ++it)
                if (mpfr::abs(it.value()) < tol)
                    it.valueRef() = 0;

        // We may have set some coefficients to zero...
        result.matrixI.prune(0, 0);
        result.matrixI.makeCompressed();
        result.checkComplexity();
    }

    return result;
}

/* This function removes all components in the matrix whose magnitude is smaller
   than the provided tolerance */
SparseGmpEigenMatrix& SparseGmpEigenMatrix::clean_new(const GmpEigenMatrix& tolerance) const
{
    SparseGmpEigenMatrix& result(*(new SparseGmpEigenMatrix(*this)));

    // We extract the tolerance
    mpreal tol(tolerance.matrixR(0,0));

    for (IndexType k = 0; k < result.matrixR.outerSize(); ++k)
        for (SparseMatrix<mpreal>::InnerIterator it(result.matrixR,k); it; ++it)
            if (mpfr::abs(it.value()) < tol)
                it.valueRef() = 0;

    // We may have set some coefficients to zero...
    result.matrixR.prune(0, 0);
    result.matrixR.makeCompressed();

    if (isComplex) {
        for (IndexType k = 0; k < result.matrixI.outerSize(); ++k)
            for (SparseMatrix<mpreal>::InnerIterator it(result.matrixI,k); it; ++it)
                if (mpfr::abs(it.value()) < tol)
                    it.valueRef() = 0;

        // We may have set some coefficients to zero...
        result.matrixI.prune(0, 0);
        result.matrixI.makeCompressed();
        result.checkComplexity();
    }

    return result;
}



/* This function extracts the nonzero elements from the matrix */
GmpEigenMatrix SparseGmpEigenMatrix::find(vector < double >& rows, vector < double >& cols) const
{
    // We start by filling dynamically allocated tables (we don't know exactly
    // how many elements we'll end up with)
    GmpEigenMatrix values;
    rows.reserve(matrixR.nonZeros());
    cols.reserve(matrixR.nonZeros());
    values.matrixR.resize(matrixR.nonZeros() + matrixI.nonZeros(), 1);
    if (isComplex)
        values.matrixI.resize(matrixR.nonZeros() + matrixI.nonZeros(), 1);
    values.isComplex = isComplex;

    IndexType index(0);
    if (!isComplex) {
        // This case is simpler
        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            for (SparseMatrix<mpreal>::InnerIterator it(matrixR,k); it; ++it) {
                rows.push_back(it.row());
                cols.push_back(it.col());
                values.matrixR(index,0) = it.value();
                ++index;
            }
        }
    } else {
        // For each column, we should merge the lists of lines with nonzero elements...
        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
            SparseMatrix<mpreal>::InnerIterator itI(matrixI,k);

            while ((itR) || (itI)) {
                if ((itR) && (itI)) {
                    if (itR.row() < itI.row()) {
                        rows.push_back(itR.row());
                        cols.push_back(itR.col());
                        values.matrixR(index,0) = itR.value();
                        values.matrixI(index,0) = 0;
                        ++itR;
                        ++index;
                    } else if (itR.row() == itI.row()) {
                        rows.push_back(itR.row());
                        cols.push_back(itR.col());
                        values.matrixR(index,0) = itR.value();
                        values.matrixI(index,0) = itI.value();
                        ++itR;
                        ++itI;
                        ++index;
                    } else {
                        rows.push_back(itI.row());
                        cols.push_back(itI.col());
                        values.matrixR(index,0) = 0;
                        values.matrixI(index,0) = itI.value();
                        ++itI;
                        ++index;
                    }
                } else if (itR) {
                    rows.push_back(itR.row());
                    cols.push_back(itR.col());
                    values.matrixR(index,0) = itR.value();
                    values.matrixI(index,0) = 0;
                    ++itR;
                    ++index;
                } else {
                    rows.push_back(itI.row());
                    cols.push_back(itI.col());
                    values.matrixR(index,0) = 0;
                    values.matrixI(index,0) = itI.value();
                    ++itI;
                    ++index;
                }
            }
        }
    }

    // We readjust the size if needed;
    values.matrixR.conservativeResize(index,1);
    values.matrixI.conservativeResize(index,1);

    return values;
}

/* This function extracts the nonzero elements from the matrix */
GmpEigenMatrix& SparseGmpEigenMatrix::find_new(vector < double >& rows, vector < double >& cols) const
{
    // We start by filling dynamically allocated tables (we don't know exactly
    // how many elements we'll end up with)
    GmpEigenMatrix& values(*(new GmpEigenMatrix));
    rows.reserve(matrixR.nonZeros());
    cols.reserve(matrixR.nonZeros());
    values.matrixR.resize(matrixR.nonZeros() + matrixI.nonZeros(), 1);
    if (isComplex)
        values.matrixI.resize(matrixR.nonZeros() + matrixI.nonZeros(), 1);
    values.isComplex = isComplex;

    IndexType index(0);
    if (!isComplex) {
        // This case is simpler
        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            for (SparseMatrix<mpreal>::InnerIterator it(matrixR,k); it; ++it) {
                rows.push_back(it.row());
                cols.push_back(it.col());
                values.matrixR(index,0) = it.value();
                ++index;
            }
        }
    } else {
        // For each column, we should merge the lists of lines with nonzero elements...
        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
            SparseMatrix<mpreal>::InnerIterator itI(matrixI,k);

            while ((itR) || (itI)) {
                if ((itR) && (itI)) {
                    if (itR.row() < itI.row()) {
                        rows.push_back(itR.row());
                        cols.push_back(itR.col());
                        values.matrixR(index,0) = itR.value();
                        values.matrixI(index,0) = 0;
                        ++itR;
                        ++index;
                    } else if (itR.row() == itI.row()) {
                        rows.push_back(itR.row());
                        cols.push_back(itR.col());
                        values.matrixR(index,0) = itR.value();
                        values.matrixI(index,0) = itI.value();
                        ++itR;
                        ++itI;
                        ++index;
                    } else {
                        rows.push_back(itI.row());
                        cols.push_back(itI.col());
                        values.matrixR(index,0) = 0;
                        values.matrixI(index,0) = itI.value();
                        ++itI;
                        ++index;
                    }
                } else if (itR) {
                    rows.push_back(itR.row());
                    cols.push_back(itR.col());
                    values.matrixR(index,0) = itR.value();
                    values.matrixI(index,0) = 0;
                    ++itR;
                    ++index;
                } else {
                    rows.push_back(itI.row());
                    cols.push_back(itI.col());
                    values.matrixR(index,0) = 0;
                    values.matrixI(index,0) = itI.value();
                    ++itI;
                    ++index;
                }
            }
        }
    }

    // We readjust the size if needed;
    values.matrixR.conservativeResize(index,1);
    values.matrixI.conservativeResize(index,1);

    return values;
}


/* This function returns the size in matlab format */
mxArray* SparseGmpEigenMatrix::size() const
{
    /* Get the size and pointers to input data */
    /*mwSize m(matrixR.rows());
    mwSize n(matrixR.cols());*/

    mxArray* plhs(mxCreateNumericMatrix(1, 2, mxDOUBLE_CLASS, mxREAL));

    double* pointer(mxGetPr(plhs));

    pointer[0] = matrixR.rows();
    pointer[1] = matrixR.cols();

    return plhs;
}



/* Reshaping a matrix b = reshape(a, [m2, n2]) */
// NOTE 1 : Here we assume that m*n = m2*n2
// NOTE 2 : These functions (like probably a few more) could certainly be much
//    improved by relying more on the Eigen framework...
SparseGmpEigenMatrix SparseGmpEigenMatrix::reshape(const IndexType& m2, const IndexType& n2) const
{
    SparseGmpEigenMatrix result;

    result.matrixR.resize(m2, n2);
    result.matrixR.reserve(matrixR.nonZeros());
    for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
        for (SparseMatrix<mpreal>::InnerIterator it(matrixR,k); it; ++it) {
            IndexType index(matrixR.rows()*it.col()+it.row());
            result.matrixR.insert(index % m2, index/m2) = it.value();
        }
    }

    result.isComplex = isComplex;
    if (isComplex) {
        result.matrixI.resize(m2, n2);
        result.matrixI.reserve(matrixI.nonZeros());
        for (IndexType k = 0; k < matrixI.outerSize(); ++k) {
            for (SparseMatrix<mpreal>::InnerIterator it(matrixI,k); it; ++it) {
                IndexType index(matrixR.rows()*it.col()+it.row());
                result.matrixI.insert(index % m2, index/m2) = it.value();
            }
        }
    }

    return result;
}

/* Reshaping a matrix b = reshape(a, [m2, n2]) */
SparseGmpEigenMatrix& SparseGmpEigenMatrix::reshape_new(const IndexType& m2, const IndexType& n2) const
{
    SparseGmpEigenMatrix& result(*(new SparseGmpEigenMatrix));

    result.matrixR.resize(m2, n2);
    result.matrixR.reserve(matrixR.nonZeros());
    for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
        for (SparseMatrix<mpreal>::InnerIterator it(matrixR,k); it; ++it) {
            IndexType index(matrixR.rows()*it.col()+it.row());
            result.matrixR.insert(index % m2, index/m2) = it.value();
        }
    }

    result.isComplex = isComplex;
    if (isComplex) {
        result.matrixI.resize(m2, n2);
        result.matrixI.reserve(matrixI.nonZeros());
        for (IndexType k = 0; k < matrixI.outerSize(); ++k) {
            for (SparseMatrix<mpreal>::InnerIterator it(matrixI,k); it; ++it) {
                IndexType index(matrixR.rows()*it.col()+it.row());
                result.matrixI.insert(index % m2, index/m2) = it.value();
            }
        }
    }

    return result;
}


/* Vertical concatenation of two matrices c = [a; b]*/
SparseGmpEigenMatrix SparseGmpEigenMatrix::vertcat(const SparseGmpEigenMatrix& b) const
{
    SparseGmpEigenMatrix result(*this);

    result.matrixR.conservativeResize(matrixR.rows()+b.matrixR.rows(), matrixR.cols());
    // It seems we cannot just call result.matrixR.block(matrixR.rows(), 0, b.matrixR.rows(), b.matrixR.cols()) = b.matrixR;
    // So let's do this by hand...
    result.matrixR.reserve(b.matrixR.nonZeros());
    for (IndexType k = 0; k < b.matrixR.outerSize(); ++k)
        for (SparseMatrix<mpreal>::InnerIterator it(b.matrixR,k); it; ++it)
            result.matrixR.insert(matrixR.rows()+it.row(), it.col()) = it.value();

    if ((isComplex) || (b.isComplex)) {
        result.isComplex = true;
        result.matrixI.conservativeResize(matrixR.rows()+b.matrixR.rows(), matrixR.cols());
        if (b.isComplex) {
            result.matrixI.reserve(b.matrixI.nonZeros());
            for (IndexType k = 0; k < b.matrixI.outerSize(); ++k)
                for (SparseMatrix<mpreal>::InnerIterator it(b.matrixI,k); it; ++it)
                    result.matrixI.insert(matrixR.rows()+it.row(), it.col()) = it.value();
        }
    }

    return result;
}

/* Vertical concatenation of two matrices c = [a; b]*/
SparseGmpEigenMatrix& SparseGmpEigenMatrix::vertcat_new(const SparseGmpEigenMatrix& b) const
{
    SparseGmpEigenMatrix& result(*(new SparseGmpEigenMatrix(*this)));

    result.matrixR.conservativeResize(matrixR.rows()+b.matrixR.rows(), matrixR.cols());
    // It seems we cannot just call result.matrixR.block(matrixR.rows(), 0, b.matrixR.rows(), b.matrixR.cols()) = b.matrixR;
    // So let's do this by hand...
    result.matrixR.reserve(b.matrixR.nonZeros());
    for (IndexType k = 0; k < b.matrixR.outerSize(); ++k)
        for (SparseMatrix<mpreal>::InnerIterator it(b.matrixR,k); it; ++it)
            result.matrixR.insert(matrixR.rows()+it.row(), it.col()) = it.value();

    if ((isComplex) || (b.isComplex)) {
        result.isComplex = true;
        result.matrixI.conservativeResize(matrixR.rows()+b.matrixR.rows(), matrixR.cols());
        if (b.isComplex) {
            result.matrixI.reserve(b.matrixI.nonZeros());
            for (IndexType k = 0; k < b.matrixI.outerSize(); ++k)
                for (SparseMatrix<mpreal>::InnerIterator it(b.matrixI,k); it; ++it)
                    result.matrixI.insert(matrixR.rows()+it.row(), it.col()) = it.value();
        }
    }

    return result;
}

/* Horizontal concatenation of two matrices c = [a, b]*/
SparseGmpEigenMatrix SparseGmpEigenMatrix::horzcat(const SparseGmpEigenMatrix& b) const
{
    SparseGmpEigenMatrix result(*this);

    result.matrixR.conservativeResize(matrixR.rows(), matrixR.cols()+b.matrixR.cols());
    // It seems we cannot just call result.matrixR.block(0, matrixR.cols(), b.matrixR.rows(), b.matrixR.cols()) = b.matrixR;
    // So let's do this by hand...

    result.matrixR.reserve(b.matrixR.nonZeros());
    for (IndexType k = 0; k < b.matrixR.outerSize(); ++k)
        for (SparseMatrix<mpreal>::InnerIterator it(b.matrixR,k); it; ++it)
            result.matrixR.insert(it.row(), matrixR.cols()+it.col()) = it.value();

    if ((isComplex) || (b.isComplex)) {
        result.isComplex = true;
        result.matrixI.conservativeResize(matrixR.rows(), matrixR.cols()+b.matrixR.cols());
        if (b.isComplex) {
            result.matrixI.reserve(b.matrixI.nonZeros());
            for (IndexType k = 0; k < b.matrixI.outerSize(); ++k)
                for (SparseMatrix<mpreal>::InnerIterator it(b.matrixI,k); it; ++it)
                    result.matrixI.insert(it.row(), matrixR.cols()+it.col()) = it.value();
        }
    }

    return result;
}

/* Horizontal concatenation of two matrices c = [a, b]*/
SparseGmpEigenMatrix& SparseGmpEigenMatrix::horzcat_new(const SparseGmpEigenMatrix& b) const
{
    SparseGmpEigenMatrix& result(*(new SparseGmpEigenMatrix(*this)));

    result.matrixR.conservativeResize(matrixR.rows(), matrixR.cols()+b.matrixR.cols());
    // It seems we cannot just call result.matrixR.block(0, matrixR.cols(), b.matrixR.rows(), b.matrixR.cols()) = b.matrixR;
    // So let's do this by hand...
    result.matrixR.reserve(b.matrixR.nonZeros());
    for (IndexType k = 0; k < b.matrixR.outerSize(); ++k)
        for (SparseMatrix<mpreal>::InnerIterator it(b.matrixR,k); it; ++it)
            result.matrixR.insert(it.row(), matrixR.cols()+it.col()) = it.value();

    if ((isComplex) || (b.isComplex)) {
        result.isComplex = true;
        result.matrixI.conservativeResize(matrixR.rows(), matrixR.cols()+b.matrixR.cols());
        if (b.isComplex) {
            result.matrixI.reserve(b.matrixI.nonZeros());
            for (IndexType k = 0; k < b.matrixI.outerSize(); ++k)
                for (SparseMatrix<mpreal>::InnerIterator it(b.matrixI,k); it; ++it)
                    result.matrixI.insert(it.row(), matrixR.cols()+it.col()) = it.value();
        }
    }

    return result;
}




/* This function extracts a submatrix or size row x cols fomr the position (i,j) */
SparseGmpEigenMatrix SparseGmpEigenMatrix::block(const IndexType& i, const IndexType& j, const IndexType& rows, const IndexType& cols) const
{
    SparseGmpEigenMatrix result;

    result.matrixR = matrixR.block(i,j,rows,cols);
    if (isComplex)
        result.matrixI = matrixI.block(i,j,rows,cols);

    result.checkComplexity();

    return result;
}

/* This function extracts a sub-matrix with single set of index b = a([1,2,3]) */
SparseGmpEigenMatrix SparseGmpEigenMatrix::subsref(const vector<vector<IndexType> >& indices) const
{
    SparseGmpEigenMatrix result;

    result.matrixR.resize(indices.size(), indices[0].size());
    result.matrixR.reserve(std::min((IndexType) matrixR.nonZeros(), (IndexType) indices.size()));
    if (isComplex) {
        result.matrixI.resize(indices.size(), indices[0].size());
        result.matrixI.reserve(std::min((IndexType) matrixI.nonZeros(), (IndexType) indices.size()));
    }

    for (IndexType x(0); x < indices.size(); ++x) {
        for (IndexType y(0); y < indices[x].size(); ++y) {
            // We separate the number indices into i and j
            IndexType i, j;
            j = indices[x][y] / matrixR.rows();
            i = indices[x][y] % matrixR.rows();

            if (matrixR.coeff(i,j) != 0)
                result.matrixR.insert(x,y) = matrixR.coeff(i,j);
            if ((isComplex) && (matrixI.coeff(i,j) != 0))
                result.matrixI.insert(x,y) = matrixI.coeff(i,j);
        }
    }

    result.matrixR.makeCompressed();
    result.matrixI.makeCompressed();
    result.checkComplexity();

    return result;
}

/* This function extracts a sub-matrix with single set of index b = a([1,2,3]) */
SparseGmpEigenMatrix& SparseGmpEigenMatrix::subsref_new(const vector<vector<IndexType> >& indices) const
{
    SparseGmpEigenMatrix& result(*(new SparseGmpEigenMatrix));

    result.matrixR.resize(indices.size(), indices[0].size());
    result.matrixR.reserve(std::min((IndexType) matrixR.nonZeros(), (IndexType) indices.size()));
    if (isComplex) {
        result.matrixI.resize(indices.size(), indices[0].size());
        result.matrixI.reserve(std::min((IndexType) matrixI.nonZeros(), (IndexType) indices.size()));
    }

    for (IndexType x(0); x < indices.size(); ++x) {
        for (IndexType y(0); y < indices[x].size(); ++y) {
            // We separate the number indices into i and j
            IndexType i, j;
            j = indices[x][y] / matrixR.rows();
            i = indices[x][y] % matrixR.rows();

            if (matrixR.coeff(i,j) != 0)
                result.matrixR.insert(x,y) = matrixR.coeff(i,j);
            if ((isComplex) && (matrixI.coeff(i,j) != 0))
                result.matrixI.insert(x,y) = matrixI.coeff(i,j);
        }
    }

    result.matrixR.makeCompressed();
    result.matrixI.makeCompressed();
    result.checkComplexity();

    return result;
}

/* This function extracts a sub-matrix with two sets of indices b = a(1:2,2) */
SparseGmpEigenMatrix SparseGmpEigenMatrix::subsref(const vector<IndexType>& indicesA, const vector<IndexType>& indicesB) const
{
    SparseGmpEigenMatrix result;

    result.matrixR.resize(indicesA.size(), indicesB.size());
    result.matrixR.reserve(std::min((IndexType) matrixR.nonZeros(), (IndexType) (indicesA.size()*indicesB.size())));
    if (isComplex) {
        result.matrixI.resize(indicesA.size(), indicesB.size());
        result.matrixI.reserve(std::min((IndexType) matrixI.nonZeros(), (IndexType) (indicesA.size()*indicesB.size())));
    }

    for (IndexType i(0); i < indicesA.size(); ++i)
        for (IndexType j(0); j < indicesB.size(); ++j) {
            if (matrixR.coeff(indicesA[i], indicesB[j]) != 0)
                result.matrixR.insert(i,j) = matrixR.coeff(indicesA[i], indicesB[j]);
            if ((isComplex) && (matrixI.coeff(indicesA[i], indicesB[j]) != 0))
                result.matrixI.insert(i,j) = matrixI.coeff(indicesA[i], indicesB[j]);
        }

    result.matrixR.makeCompressed();
    result.matrixI.makeCompressed();
    result.checkComplexity();

    return result;
}

/* This function extracts a sub-matrix with two sets of indices b = a(1:2,2) */
SparseGmpEigenMatrix& SparseGmpEigenMatrix::subsref_new(const vector<IndexType>& indicesA, const vector<IndexType>& indicesB) const
{
    SparseGmpEigenMatrix& result(*(new SparseGmpEigenMatrix));

    result.matrixR.resize(indicesA.size(), indicesB.size());
    result.matrixR.reserve(std::min((IndexType) matrixR.nonZeros(), (IndexType) (indicesA.size()*indicesB.size())));
    if (isComplex) {
        result.matrixI.resize(indicesA.size(), indicesB.size());
        result.matrixI.reserve(std::min((IndexType) matrixI.nonZeros(), (IndexType) (indicesA.size()*indicesB.size())));
    }

    for (IndexType i(0); i < indicesA.size(); ++i)
        for (IndexType j(0); j < indicesB.size(); ++j) {
            if (matrixR.coeff(indicesA[i], indicesB[j]) != 0)
                result.matrixR.insert(i,j) = matrixR.coeff(indicesA[i], indicesB[j]);
            if ((isComplex) && (matrixI.coeff(indicesA[i], indicesB[j]) != 0))
                result.matrixI.insert(i,j) = matrixI.coeff(indicesA[i], indicesB[j]);
        }

    result.matrixR.makeCompressed();
    result.matrixI.makeCompressed();
    result.checkComplexity();

    return result;
}

/* Assigning value to a sub-matrix a([1 2 3]) = b */
void SparseGmpEigenMatrix::subsasgn(const std::vector<IndexType>& indices, const SparseGmpEigenMatrix& b)
{
    if ((!isComplex) && (b.isComplex)) {
        isComplex = true;
        matrixI.resize(matrixR.rows(), matrixR.cols());
        matrixI.reserve(b.matrixI.nonZeros());
    }

    for (IndexType x(0); x < indices.size(); ++x) {
        // We separate the number indices into i and j
        IndexType i, j;
        j = indices[x] / matrixR.rows();
        i = indices[x] % matrixR.rows();

        if ((matrixR.coeff(i,j) != 0) || (b.matrixR.coeff(x,0) != 0))
            matrixR.coeffRef(i,j) = b.matrixR.coeff(x,0);

        if ((b.isComplex) && (b.matrixI.coeff(x,0) != 0))
            matrixI.coeffRef(i,j) = b.matrixI.coeff(x,0);
        else if ((isComplex) && (matrixI.coeff(i,j) != 0))
            matrixI.coeffRef(i,j) = 0;
    }

    // We remove the elements which have been set to zero
    matrixR.prune(0, 0);
    matrixI.prune(0, 0);

    matrixR.makeCompressed();
    matrixI.makeCompressed();
    if (isComplex)
        checkComplexity();
}

/* Assigning value to a sub-matrix a(1:2,:) = b */
void SparseGmpEigenMatrix::subsasgn(const std::vector<IndexType>& indicesA, const std::vector<IndexType>& indicesB, const SparseGmpEigenMatrix& b)
{
    if ((!isComplex) && (b.isComplex)) {
        isComplex = true;
        matrixI.resize(matrixR.rows(), matrixR.cols());
        matrixI.reserve(b.matrixI.nonZeros());
    }

    for (IndexType i(0); i < indicesA.size(); ++i)
        for (IndexType j(0); j < indicesB.size(); ++j) {
            if ((matrixR.coeff(indicesA[i],indicesB[j]) != 0) || (b.matrixR.coeff(i,j) != 0))
                matrixR.coeffRef(indicesA[i], indicesB[j]) = b.matrixR.coeff(i,j);

            if ((b.isComplex) && (b.matrixI.coeff(i,j) != 0))
                matrixI.coeffRef(indicesA[i], indicesB[j]) = b.matrixI.coeff(i,j);
            else if ((isComplex) && (matrixI.coeff(indicesA[i], indicesB[j]) != 0))
                matrixI.coeffRef(indicesA[i], indicesB[j]) = 0;
        }

    // We remove the elements which have been set to zero
    matrixR.prune(0, 0);
    matrixI.prune(0, 0);

    matrixR.makeCompressed();
    matrixI.makeCompressed();
    if (isComplex)
        checkComplexity();
}

/* Creation of a diagonal matrix b = diag(a, k) */
SparseGmpEigenMatrix SparseGmpEigenMatrix::diagCreate(const IndexType& k) const
{
    SparseGmpEigenMatrix result;

    result.isComplex = isComplex;
    result.matrixR.resize(matrixR.rows()+std::abs(k), matrixR.rows()+std::abs(k));
    result.matrixR.reserve(matrixR.nonZeros());
    if (isComplex) {
        result.matrixI.resize(matrixR.rows()+std::abs(k), matrixR.rows()+std::abs(k));
        result.matrixI.reserve(matrixI.nonZeros());
    }

    // We copy the real data (the current matrix is assumed to be of size nx1)
    for (SparseMatrix<mpreal>::InnerIterator it(matrixR,0); it; ++it)
        result.matrixR.insert(max((IndexType)0,-k)+it.row(), max((IndexType)0,k)+it.row()) = it.value();

    // We copy the imaginary data
    if (isComplex) {
        for (SparseMatrix<mpreal>::InnerIterator it(matrixI,0); it; ++it)
            result.matrixI.insert(max((IndexType)0,-k)+it.row(), max((IndexType)0,k)+it.row()) = it.value();
    }

    return result;
}

/* Creation of a diagonal matrix b = diag(a, k) */
SparseGmpEigenMatrix& SparseGmpEigenMatrix::diagCreate_new(const IndexType& k) const
{
    SparseGmpEigenMatrix& result(*(new SparseGmpEigenMatrix));

    result.isComplex = isComplex;
    result.matrixR.resize(matrixR.rows()+std::abs(k), matrixR.rows()+std::abs(k));
    result.matrixR.reserve(matrixR.nonZeros());
    if (isComplex) {
        result.matrixI.resize(matrixR.rows()+std::abs(k), matrixR.rows()+std::abs(k));
        result.matrixI.reserve(matrixI.nonZeros());
    }

    // We copy the real data (the current matrix is assumed to be of size nx1)
    for (SparseMatrix<mpreal>::InnerIterator it(matrixR,0); it; ++it)
        result.matrixR.insert(max((IndexType)0,-k)+it.row(), max((IndexType)0,k)+it.row()) = it.value();

    // We copy the imaginary data
    if (isComplex) {
        for (SparseMatrix<mpreal>::InnerIterator it(matrixI,0); it; ++it)
            result.matrixI.insert(max((IndexType)0,-k)+it.row(), max((IndexType)0,k)+it.row()) = it.value();
    }

    return result;
}

/* Extraction of diagonal elements b = diag(a, k) */
SparseGmpEigenMatrix SparseGmpEigenMatrix::diagExtract(const IndexType& k) const
{
    SparseGmpEigenMatrix result, tmp;

    // First, we select the right block
    tmp.matrixR = matrixR.block(max((IndexType)0,-k), max((IndexType)0,k), matrixR.rows()-max((IndexType)0,-k), matrixR.cols()-max((IndexType)0,k));
    // Now we extract the diagonal
    result.matrixR.resize(min(tmp.matrixR.rows(), tmp.matrixR.cols()), 1);
    result.matrixR.reserve(min(min(tmp.matrixR.rows(), tmp.matrixR.cols()), tmp.matrixR.nonZeros()));
    for (IndexType k0 = 0; k0 < tmp.matrixR.outerSize(); ++k0) {
        for (SparseMatrix<mpreal>::InnerIterator it(tmp.matrixR,k0); it; ++it) {
            if (it.row() > it.col())
                break;
            else if (it.row() == it.col())
                result.matrixR.insert(it.row(), 0) = it.value();
        }
    }
    result.matrixR.makeCompressed();
    result.isComplex = isComplex;

    if (isComplex) {
        // First, we select the right block
        tmp.matrixI = matrixI.block(max((IndexType)0,-k), max((IndexType)0,k), matrixR.rows()-max((IndexType)0,-k), matrixR.cols()-max((IndexType)0,k));
        // Now we extract the diagonal
        result.matrixI.resize(min(tmp.matrixI.rows(), tmp.matrixI.cols()), 1);
        result.matrixI.reserve(min(min(tmp.matrixI.rows(), tmp.matrixI.cols()), tmp.matrixI.nonZeros()));
        for (IndexType k0 = 0; k0 < tmp.matrixI.outerSize(); ++k0) {
            for (SparseMatrix<mpreal>::InnerIterator it(tmp.matrixI,k0); it; ++it) {
                if (it.row() > it.col())
                    break;
                else if (it.row() == it.col())
                    result.matrixI.insert(it.row(), 0) = it.value();
            }
        }
        result.matrixI.makeCompressed();
        result.checkComplexity();
    }

    return result;
}

/* Extraction of diagonal elements b = diag(a, k) */
SparseGmpEigenMatrix& SparseGmpEigenMatrix::diagExtract_new(const IndexType& k) const
{
    SparseGmpEigenMatrix& result(*(new SparseGmpEigenMatrix));
    SparseGmpEigenMatrix tmp;

    // First, we select the right block
    tmp.matrixR = matrixR.block(max((IndexType)0,-k), max((IndexType)0,k), matrixR.rows()-max((IndexType)0,-k), matrixR.cols()-max((IndexType)0,k));
    // Now we extract the diagonal
    result.matrixR.resize(min(tmp.matrixR.rows(), tmp.matrixR.cols()), 1);
    result.matrixR.reserve(min(min(tmp.matrixR.rows(), tmp.matrixR.cols()), tmp.matrixR.nonZeros()));
    for (IndexType k0 = 0; k0 < tmp.matrixR.outerSize(); ++k0) {
        for (SparseMatrix<mpreal>::InnerIterator it(tmp.matrixR,k0); it; ++it) {
            if (it.row() > it.col())
                break;
            else if (it.row() == it.col())
                result.matrixR.insert(it.row(), 0) = it.value();
        }
    }
    result.matrixR.makeCompressed();
    result.isComplex = isComplex;

    if (isComplex) {
        // First, we select the right block
        tmp.matrixI = matrixI.block(max((IndexType)0,-k), max((IndexType)0,k), matrixR.rows()-max((IndexType)0,-k), matrixR.cols()-max((IndexType)0,k));
        // Now we extract the diagonal
        result.matrixI.resize(min(tmp.matrixI.rows(), tmp.matrixI.cols()), 1);
        result.matrixI.reserve(min(min(tmp.matrixI.rows(), tmp.matrixI.cols()), tmp.matrixI.nonZeros()));
        for (IndexType k0 = 0; k0 < tmp.matrixI.outerSize(); ++k0) {
            for (SparseMatrix<mpreal>::InnerIterator it(tmp.matrixI,k0); it; ++it) {
                if (it.row() > it.col())
                    break;
                else if (it.row() == it.col())
                    result.matrixI.insert(it.row(), 0) = it.value();
            }
        }
        result.matrixI.makeCompressed();
        result.checkComplexity();
    }

    return result;
}




/* Some class operators */



/* The unary substraction operators b = -a */
SparseGmpEigenMatrix SparseGmpEigenMatrix::operator-() const
{
    SparseGmpEigenMatrix result;

    result.matrixR = -matrixR;
    result.isComplex = isComplex;
    if (isComplex)
        result.matrixI = -matrixI;

    return result;
}

/* The unary substraction b = -a */
SparseGmpEigenMatrix& SparseGmpEigenMatrix::uminus_new() const
{
    SparseGmpEigenMatrix& result(*(new SparseGmpEigenMatrix));

    result.matrixR = -matrixR;
    result.isComplex = isComplex;
    if (isComplex)
        result.matrixI = -matrixI;

    return result;
}

/* The tansposition b = a' */
SparseGmpEigenMatrix SparseGmpEigenMatrix::transpose() const
{
    SparseGmpEigenMatrix result;

    result.isComplex = isComplex;
    result.matrixR = matrixR.transpose();
    if (isComplex)
        result.matrixI = matrixI.transpose();

    return result;
}

/* The tansposition b = a' */
SparseGmpEigenMatrix& SparseGmpEigenMatrix::transpose_new() const
{
    SparseGmpEigenMatrix& result(*(new SparseGmpEigenMatrix));

    result.isComplex = isComplex;
    result.matrixR = matrixR.transpose();
    if (isComplex)
        result.matrixI = matrixI.transpose();

    return result;
}

/* The conjugate transposition b = a' */
SparseGmpEigenMatrix SparseGmpEigenMatrix::ctranspose() const
{
    SparseGmpEigenMatrix result;

    result.isComplex = isComplex;
    result.matrixR = matrixR.transpose();
    if (isComplex)
        result.matrixI = -matrixI.transpose();

    return result;
}

/* The conjugate transposition b = a' */
SparseGmpEigenMatrix& SparseGmpEigenMatrix::ctranspose_new() const
{
    SparseGmpEigenMatrix& result(*(new SparseGmpEigenMatrix));

    result.isComplex = isComplex;
    result.matrixR = matrixR.transpose();
    if (isComplex)
        result.matrixI = -matrixI.transpose();

    return result;
}


/* Element-wise complex conjugation b = conj(a) */
SparseGmpEigenMatrix SparseGmpEigenMatrix::conj() const
{
    SparseGmpEigenMatrix result;

    result.isComplex = isComplex;
    result.matrixR = matrixR;
    if (isComplex)
        result.matrixI = -matrixI;

    return result;
}

/* Element-wise complex conjugation b = conj(a) */
SparseGmpEigenMatrix& SparseGmpEigenMatrix::conj_new() const
{
    SparseGmpEigenMatrix& result(*(new SparseGmpEigenMatrix));

    result.isComplex = isComplex;
    result.matrixR = matrixR;
    if (isComplex)
        result.matrixI = -matrixI;

    return result;
}

/* Real part b = real(a) */
SparseGmpEigenMatrix SparseGmpEigenMatrix::real() const
{
    SparseGmpEigenMatrix result;

    result.isComplex = false;
    result.matrixR = matrixR;

    return result;
}

/* Real part b = real(a) */
SparseGmpEigenMatrix& SparseGmpEigenMatrix::real_new() const
{
    SparseGmpEigenMatrix& result(*(new SparseGmpEigenMatrix));

    result.isComplex = false;
    result.matrixR = matrixR;

    return result;
}

/* Imaginary part b = imag(a) */
SparseGmpEigenMatrix SparseGmpEigenMatrix::imag() const
{
    SparseGmpEigenMatrix result;

    result.isComplex = false;
    if (isComplex)
        result.matrixR = matrixI;
    else
        result.matrixR.resize(matrixR.rows(), matrixR.cols());

    return result;
}

/* Imaginary part b = imag(a) */
SparseGmpEigenMatrix& SparseGmpEigenMatrix::imag_new() const
{
    SparseGmpEigenMatrix& result(*(new SparseGmpEigenMatrix));

    result.isComplex = false;
    if (isComplex)
        result.matrixR = matrixI;
    else
        result.matrixR.resize(matrixR.rows(), matrixR.cols());

    return result;
}


// Various integer parts
/* Rounds to nearest integer b = round(a) */
SparseGmpEigenMatrix SparseGmpEigenMatrix::round() const
{
    SparseGmpEigenMatrix result(*this);

    for (IndexType k = 0; k < matrixR.outerSize(); ++k)
        for (SparseMatrix<mpreal>::InnerIterator it(matrixR,k); it; ++it)
            it.valueRef() = mpfr::round(it.value());

    // We may have set some coefficients to zero...
    result.matrixR.prune(0, 0);
    result.matrixR.makeCompressed();

    if (isComplex) {
        for (IndexType k = 0; k < matrixI.outerSize(); ++k)
            for (SparseMatrix<mpreal>::InnerIterator it(matrixI,k); it; ++it)
                it.valueRef() = mpfr::round(it.value());

        // We may have set some coefficients to zero...
        result.matrixI.prune(0, 0);
        result.matrixI.makeCompressed();
        result.checkComplexity();
    }

    return result;
}

/* Rounds to nearest integer b = round(a) */
SparseGmpEigenMatrix& SparseGmpEigenMatrix::round_new() const
{
    SparseGmpEigenMatrix& result(*(new SparseGmpEigenMatrix(*this)));

    for (IndexType k = 0; k < result.matrixR.outerSize(); ++k)
        for (SparseMatrix<mpreal>::InnerIterator it(result.matrixR,k); it; ++it)
            it.valueRef() = mpfr::round(it.value());

    // We may have set some coefficients to zero...
    result.matrixR.prune(0, 0);
    result.matrixR.makeCompressed();

    if (isComplex) {
        for (IndexType k = 0; k < result.matrixI.outerSize(); ++k)
            for (SparseMatrix<mpreal>::InnerIterator it(result.matrixI,k); it; ++it)
                it.valueRef() = mpfr::round(it.value());

        // We may have set some coefficients to zero...
        result.matrixI.prune(0, 0);
        result.matrixI.makeCompressed();
        result.checkComplexity();
    }

    return result;
}

/* Next integer towards -Inf b = floor(a) */
SparseGmpEigenMatrix SparseGmpEigenMatrix::floor() const
{
    SparseGmpEigenMatrix result(*this);

    for (IndexType k = 0; k < result.matrixR.outerSize(); ++k)
        for (SparseMatrix<mpreal>::InnerIterator it(result.matrixR,k); it; ++it)
            it.valueRef() = mpfr::floor(it.value());

    // We may have set some coefficients to zero...
    result.matrixR.prune(0, 0);
    result.matrixR.makeCompressed();

    if (isComplex) {
        for (IndexType k = 0; k < result.matrixI.outerSize(); ++k)
            for (SparseMatrix<mpreal>::InnerIterator it(result.matrixI,k); it; ++it)
                it.valueRef() = mpfr::floor(it.value());

        // We may have set some coefficients to zero...
        result.matrixI.prune(0, 0);
        result.matrixI.makeCompressed();
        result.checkComplexity();
    }

    return result;
}

/* Next integer towards -Inf b = floor(a) */
SparseGmpEigenMatrix& SparseGmpEigenMatrix::floor_new() const
{
    SparseGmpEigenMatrix& result(*(new SparseGmpEigenMatrix(*this)));

    for (IndexType k = 0; k < result.matrixR.outerSize(); ++k)
        for (SparseMatrix<mpreal>::InnerIterator it(result.matrixR,k); it; ++it)
            it.valueRef() = mpfr::floor(it.value());

    // We may have set some coefficients to zero...
    result.matrixR.prune(0, 0);
    result.matrixR.makeCompressed();

    if (isComplex) {
        for (IndexType k = 0; k < result.matrixI.outerSize(); ++k)
            for (SparseMatrix<mpreal>::InnerIterator it(result.matrixI,k); it; ++it)
                it.valueRef() = mpfr::floor(it.value());

        // We may have set some coefficients to zero...
        result.matrixI.prune(0, 0);
        result.matrixI.makeCompressed();
        result.checkComplexity();
    }

    return result;
}

/* Next integer towards +Inf b = ceil(a) */
SparseGmpEigenMatrix SparseGmpEigenMatrix::ceil() const
{
    SparseGmpEigenMatrix result(*this);

    for (IndexType k = 0; k < result.matrixR.outerSize(); ++k)
        for (SparseMatrix<mpreal>::InnerIterator it(result.matrixR,k); it; ++it)
            it.valueRef() = mpfr::ceil(it.value());

    // We may have set some coefficients to zero...
    result.matrixR.prune(0, 0);
    result.matrixR.makeCompressed();

    if (isComplex) {
        for (IndexType k = 0; k < result.matrixI.outerSize(); ++k)
            for (SparseMatrix<mpreal>::InnerIterator it(result.matrixI,k); it; ++it)
                it.valueRef() = mpfr::ceil(it.value());

        // We may have set some coefficients to zero...
        result.matrixI.prune(0, 0);
        result.matrixI.makeCompressed();
        result.checkComplexity();
    }

    return result;
}

/* Next integer towards +Inf b = ceil(a) */
SparseGmpEigenMatrix& SparseGmpEigenMatrix::ceil_new() const
{
    SparseGmpEigenMatrix& result(*(new SparseGmpEigenMatrix(*this)));

    for (IndexType k = 0; k < result.matrixR.outerSize(); ++k)
        for (SparseMatrix<mpreal>::InnerIterator it(result.matrixR,k); it; ++it)
            it.valueRef() = mpfr::ceil(it.value());

    // We may have set some coefficients to zero...
    result.matrixR.prune(0, 0);
    result.matrixR.makeCompressed();

    if (isComplex) {
        for (IndexType k = 0; k < result.matrixI.outerSize(); ++k)
            for (SparseMatrix<mpreal>::InnerIterator it(result.matrixI,k); it; ++it)
                it.valueRef() = mpfr::ceil(it.value());

        // We may have set some coefficients to zero...
        result.matrixI.prune(0, 0);
        result.matrixI.makeCompressed();
        result.checkComplexity();
    }

    return result;
}

/* Next integer towards 0 b = trunc(a) */
SparseGmpEigenMatrix SparseGmpEigenMatrix::trunc() const
{
    SparseGmpEigenMatrix result(*this);

    for (IndexType k = 0; k < result.matrixR.outerSize(); ++k)
        for (SparseMatrix<mpreal>::InnerIterator it(result.matrixR,k); it; ++it)
            it.valueRef() = mpfr::trunc(it.value());

    // We may have set some coefficients to zero...
    result.matrixR.prune(0, 0);
    result.matrixR.makeCompressed();

    if (isComplex) {
        for (IndexType k = 0; k < result.matrixI.outerSize(); ++k)
            for (SparseMatrix<mpreal>::InnerIterator it(result.matrixI,k); it; ++it)
                it.valueRef() = mpfr::trunc(it.value());

        // We may have set some coefficients to zero...
        result.matrixI.prune(0, 0);
        result.matrixI.makeCompressed();
        result.checkComplexity();
    }

    return result;
}

/* Next integer towards 0 b = trunc(a) */
SparseGmpEigenMatrix& SparseGmpEigenMatrix::trunc_new() const
{
    SparseGmpEigenMatrix& result(*(new SparseGmpEigenMatrix(*this)));

    for (IndexType k = 0; k < result.matrixR.outerSize(); ++k)
        for (SparseMatrix<mpreal>::InnerIterator it(result.matrixR,k); it; ++it)
            it.valueRef() = mpfr::trunc(it.value());

    // We may have set some coefficients to zero...
    result.matrixR.prune(0, 0);
    result.matrixR.makeCompressed();

    if (isComplex) {
        for (IndexType k = 0; k < result.matrixI.outerSize(); ++k)
            for (SparseMatrix<mpreal>::InnerIterator it(result.matrixI,k); it; ++it)
                it.valueRef() = mpfr::trunc(it.value());

        // We may have set some coefficients to zero...
        result.matrixI.prune(0, 0);
        result.matrixI.makeCompressed();
        result.checkComplexity();
    }

    return result;
}


// Concerning the addition operators:
// Since addition of a scalar to a sparse matrix is not a sparse object anymore,
// We only support addition between two sparse matrices of the same size
SparseGmpEigenMatrix& SparseGmpEigenMatrix::operator+=(const SparseGmpEigenMatrix& b)
{
    matrixR += b.matrixR;
    if (b.isComplex) {
        if (isComplex) {
            matrixI += b.matrixI;
        } else {
            isComplex = true;
            matrixI = b.matrixI;
        }
    }

    // We may have set some coefficients to zero...
    matrixR.prune(0, 0);
    matrixR.makeCompressed();
    if (isComplex) {
        matrixI.prune(0, 0);
        matrixI.makeCompressed();
        checkComplexity();
    }

    return *this;
}

/* The addition operator c = a + b */
SparseGmpEigenMatrix SparseGmpEigenMatrix::operator+(const SparseGmpEigenMatrix& b) const
{
    SparseGmpEigenMatrix result;

    result.isComplex = isComplex;
    result.matrixR = matrixR + b.matrixR;
    if (b.isComplex) {
        if (isComplex) {
            result.matrixI = matrixI + b.matrixI;
        } else {
            result.isComplex = true;
            result.matrixI = b.matrixI;
        }
    } else if (isComplex) {
        result.matrixI = matrixI;
    }

    // We may have set some coefficients to zero...
    result.matrixR.prune(0, 0);
    result.matrixR.makeCompressed();
    if (result.isComplex) {
        result.matrixI.prune(0, 0);
        result.matrixI.makeCompressed();
        result.checkComplexity();
    }

    return result;
}


/* The following function is a special implementation of the + operator,
   which is tailored for the matlab interface: it computes the sum of two
   numbers, like operator+, but puts the result in a new instance of the
   class that is allocated dynamically with the 'new' function. This allows
   the result of this operation to remain in memory between two hands out of
   the interface to matlab.

   Here, we keep track of the created object through a reference for optimal
   clarity. A version with pointers is also given below.

   Note that any number produced with this function should be
   either manually managed, or encapsulated into a matlab object (in which
   case matlab will take care of clearing it when it needs it no more).
   Failure to do so will result in non-freed memory (just like alocation of
   memory with pointers needs to be carefully taken care of).*/
SparseGmpEigenMatrix& SparseGmpEigenMatrix::plus_new(const SparseGmpEigenMatrix& b) const
{
    SparseGmpEigenMatrix& result(*(new SparseGmpEigenMatrix));

    result.isComplex = isComplex;
    result.matrixR = matrixR + b.matrixR;
    if (b.isComplex) {
        if (isComplex) {
            result.matrixI = matrixI + b.matrixI;
        } else {
            result.isComplex = true;
            result.matrixI = b.matrixI;
        }
    } else if (isComplex) {
        result.matrixI = matrixI;
    }

    // We may have set some coefficients to zero...
    result.matrixR.prune(0, 0);
    result.matrixR.makeCompressed();
    if (result.isComplex) {
        result.matrixI.prune(0, 0);
        result.matrixI.makeCompressed();
        result.checkComplexity();
    }

    return result;
}

/* The substraction operator a -= b */
SparseGmpEigenMatrix& SparseGmpEigenMatrix::operator-=(const SparseGmpEigenMatrix& b)
{
    matrixR -= b.matrixR;
    if (b.isComplex) {
        if (isComplex) {
            matrixI -= b.matrixI;
        } else {
            isComplex = true;
            matrixI = -b.matrixI;
        }
    }

    // We may have set some coefficients to zero...
    matrixR.prune(0, 0);
    matrixR.makeCompressed();
    if (isComplex) {
        matrixI.prune(0, 0);
        matrixI.makeCompressed();
        checkComplexity();
    }

    return *this;
}

/* The substraction operator c = a-b */
SparseGmpEigenMatrix SparseGmpEigenMatrix::operator-(const SparseGmpEigenMatrix& b) const
{
    SparseGmpEigenMatrix result;

    result.isComplex = isComplex;
    result.matrixR = matrixR - b.matrixR;
    if (b.isComplex) {
        if (isComplex) {
            result.matrixI = matrixI - b.matrixI;
        } else {
            result.isComplex = true;
            result.matrixI = -b.matrixI;
        }
    } else if (isComplex) {
        result.matrixI = matrixI;
    }

    // We may have set some coefficients to zero...
    result.matrixR.prune(0, 0);
    result.matrixR.makeCompressed();
    if ((isComplex) && (b.isComplex)) {
        result.matrixI.prune(0, 0);
        result.matrixI.makeCompressed();
        result.checkComplexity();
    }

    return result;
}

/* The substraction c = a-b */
SparseGmpEigenMatrix& SparseGmpEigenMatrix::minus_new(const SparseGmpEigenMatrix& b) const
{
    SparseGmpEigenMatrix& result(*(new SparseGmpEigenMatrix));

    result.isComplex = isComplex;
    result.matrixR = matrixR - b.matrixR;
    if (b.isComplex) {
        if (isComplex) {
            result.matrixI = matrixI - b.matrixI;
        } else {
            result.isComplex = true;
            result.matrixI = -b.matrixI;
        }
    } else if (isComplex) {
        result.matrixI = matrixI;
    }

    // We may have set some coefficients to zero...
    result.matrixR.prune(0, 0);
    result.matrixR.makeCompressed();
    if ((isComplex) && (b.isComplex)) {
        result.matrixI.prune(0, 0);
        result.matrixI.makeCompressed();
        result.checkComplexity();
    }

    return result;
}


/* The element-wise multiplication c = a.*b */
SparseGmpEigenMatrix SparseGmpEigenMatrix::times(const SparseGmpEigenMatrix& b) const
{
    SparseGmpEigenMatrix result;

    // We want to support nxm matrix multiplied by 1x1 matrix (which we see as a scalar)
    if ((numel() != 1) && (b.numel() == 1)) {
        if (b.isComplex) {
            if (isComplex) {
                result.matrixR = matrixR*b.matrixR.coeff(0,0) - matrixI*b.matrixI.coeff(0,0);
                result.matrixI = matrixR*b.matrixI.coeff(0,0) + matrixI*b.matrixR.coeff(0,0);
            } else {
                result.matrixR = matrixR*b.matrixR.coeff(0,0);
                result.matrixI = matrixR*b.matrixI.coeff(0,0);
            }
        } else if (isComplex) {
            result.matrixR = matrixR*b.matrixR.coeff(0,0);
            result.matrixI = matrixI*b.matrixR.coeff(0,0);
        } else {
            result.matrixR = matrixR*b.matrixR.coeff(0,0);
        }

        // We may have set some coefficients to zero...
        result.matrixR.prune(0, 0);
        result.matrixR.makeCompressed();
        if ((isComplex) && (b.isComplex)) {
            result.matrixI.prune(0, 0);
            result.matrixI.makeCompressed();
        }
        result.checkComplexity();
    } else if ((numel() == 1) && (b.numel() != 1)) {
        if (b.isComplex) {
            if (isComplex) {
                result.matrixR = matrixR.coeff(0,0)*b.matrixR - matrixI.coeff(0,0)*b.matrixI;
                result.matrixI = matrixR.coeff(0,0)*b.matrixI + matrixI.coeff(0,0)*b.matrixR;
            } else {
                result.matrixR = matrixR.coeff(0,0)*b.matrixR;
                result.matrixI = matrixR.coeff(0,0)*b.matrixI;
            }
        } else if (isComplex) {
            result.matrixR = matrixR.coeff(0,0)*b.matrixR;
            result.matrixI = matrixI.coeff(0,0)*b.matrixR;
        } else {
            result.matrixR = matrixR.coeff(0,0)*b.matrixR;
        }

        // We may have set some coefficients to zero...
        result.matrixR.prune(0, 0);
        result.matrixR.makeCompressed();
        if ((isComplex) && (b.isComplex)) {
            result.matrixI.prune(0, 0);
            result.matrixI.makeCompressed();
        }
        result.checkComplexity();
    } else {
        // Eigen doesn't support this case, so we do it by hand...

        // setting the output size
        result.matrixR.resize(matrixR.rows(), matrixR.cols());
        if ((isComplex) || (b.isComplex))
            result.matrixI.resize(matrixR.rows(), matrixR.cols());

        // We will first copy the data to triplets
        vector< Triplet<mpreal> > tripletListR, tripletListI;
        tripletListR.reserve(std::min(matrixR.nonZeros(), b.matrixR.nonZeros()) + std::min(matrixI.nonZeros(), b.matrixI.nonZeros()));
        if ((isComplex) || (b.isComplex))
            tripletListI.reserve(std::min(matrixR.nonZeros(), b.matrixI.nonZeros()) + std::min(matrixI.nonZeros(), b.matrixR.nonZeros()));

        // Now for each column, we merge the lists of lines with nonzero elements of
        // both matrices. The cases are slightly different depending on the complexity
        // of both matrices
        if (b.isComplex) {
            if (isComplex) {
                for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
                    SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
                    SparseMatrix<mpreal>::InnerIterator itI(matrixI,k);
                    SparseMatrix<mpreal>::InnerIterator itRb(b.matrixR,k);
                    SparseMatrix<mpreal>::InnerIterator itIb(b.matrixI,k);

                    while (((itR) || (itI)) && ((itRb) || (itIb))) {
                        IndexType row, rowb;
                        if ((itR) && (itI))
                            row = min(itR.row(), itI.row());
                        else if (itR)
                            row = itR.row();
                        else
                            row = itI.row();
                        if ((itRb) && (itIb))
                            rowb = min(itRb.row(), itIb.row());
                        else if (itRb)
                            rowb = itRb.row();
                        else
                            rowb = itIb.row();

                        if (row < rowb) {
                            if ((!itI) || ((itR) && (itR.row() < itI.row()))) {
                                ++itR;
                            } else if ((itR) && (itI) && (itR.row() == itI.row())) {
                                ++itR;
                                ++itI;
                            } else {
                                ++itI;
                            }
                        } else if (row == rowb) {
                            if ((!itI) || ((itR) && (itR.row() < itI.row()))) {
                                if ((!itIb) || ((itRb) && (itRb.row() < itIb.row()))) {
                                    tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()*itRb.value()));
                                    ++itRb;
                                } else if ((itRb) && (itIb) && (itRb.row() == itIb.row())) {
                                    tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()*itRb.value()));
                                    tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itR.value()*itIb.value()));
                                    ++itRb;
                                    ++itIb;
                                } else {
                                    tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itR.value()*itIb.value()));
                                    ++itIb;
                                }
                                ++itR;
                            } else if ((itR) && (itI) && (itR.row() == itI.row())) {
                                if ((!itIb) || ((itRb) && (itRb.row() < itIb.row()))) {
                                    tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()*itRb.value()));
                                    tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itI.value()*itRb.value()));
                                    ++itRb;
                                } else if ((itRb) && (itIb) && (itRb.row() == itIb.row())) {
                                    tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()*itRb.value() - itI.value()*itIb.value()));
                                    tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itR.value()*itIb.value() + itI.value()*itRb.value()));
                                    ++itRb;
                                    ++itIb;
                                } else {
                                    tripletListR.push_back(Triplet<mpreal>(itR.row(), k, -itI.value()*itIb.value()));
                                    tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itR.value()*itIb.value()));
                                    ++itIb;
                                }
                                ++itR;
                                ++itI;
                            } else {
                                if ((!itIb) || ((itRb) && (itRb.row() < itIb.row()))) {
                                    tripletListI.push_back(Triplet<mpreal>(itI.row(), k, itI.value()*itRb.value()));
                                    ++itRb;
                                } else if ((itRb) && (itIb) && (itRb.row() == itIb.row())) {
                                    tripletListR.push_back(Triplet<mpreal>(itI.row(), k, -itI.value()*itIb.value()));
                                    tripletListI.push_back(Triplet<mpreal>(itI.row(), k, itI.value()*itRb.value()));
                                    ++itRb;
                                    ++itIb;
                                } else {
                                    tripletListR.push_back(Triplet<mpreal>(itI.row(), k, -itI.value()*itIb.value()));
                                    ++itIb;
                                }
                                ++itI;
                            }
                        } else {
                            if ((!itIb) || ((itRb) && (itRb.row() < itIb.row()))) {
                                ++itRb;
                            } else if ((itRb) && (itIb) && (itRb.row() == itIb.row())) {
                                ++itRb;
                                ++itIb;
                            } else {
                                ++itIb;
                            }
                        }
                    }
                }
            } else {
                for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
                    SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
                    SparseMatrix<mpreal>::InnerIterator itRb(b.matrixR,k);
                    SparseMatrix<mpreal>::InnerIterator itIb(b.matrixI,k);

                    while ((itR) && ((itRb) || (itIb))) {
                        IndexType rowb;
                        if ((itRb) && (itIb))
                            rowb = min(itRb.row(), itIb.row());
                        else if (itRb)
                            rowb = itRb.row();
                        else
                            rowb = itIb.row();

                        if (itR.row() < rowb) {
                            ++itR;
                        } else if (itR.row() == rowb) {
                            if ((!itIb) || ((itRb) && (itRb.row() < itIb.row()))) {
                                tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()*itRb.value()));
                                ++itRb;
                            } else if ((itRb) && (itIb) && (itRb.row() == itIb.row())) {
                                tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()*itRb.value()));
                                tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itR.value()*itIb.value()));
                                ++itRb;
                                ++itIb;
                            } else {
                                tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itR.value()*itIb.value()));
                                ++itIb;
                            }
                            ++itR;
                        } else {
                            if ((!itIb) || ((itRb) && (itRb.row() < itIb.row()))) {
                                ++itRb;
                            } else if ((itRb) && (itIb) && (itRb.row() == itIb.row())) {
                                ++itRb;
                                ++itIb;
                            } else {
                                ++itIb;
                            }
                        }
                    }
                }
            }
        } else if (isComplex) {
            for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
                SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
                SparseMatrix<mpreal>::InnerIterator itI(matrixI,k);
                SparseMatrix<mpreal>::InnerIterator itRb(b.matrixR,k);

                while (((itR) || (itI)) && (itRb)) {
                    IndexType row;
                    if ((itR) && (itI))
                        row = min(itR.row(), itI.row());
                    else if (itR)
                        row = itR.row();
                    else
                        row = itI.row();

                    if (row < itRb.row()) {
                        if ((!itI) || ((itR) && (itR.row() < itI.row()))) {
                            ++itR;
                        } else if ((itR) && (itI) && (itR.row() == itI.row())) {
                            ++itR;
                            ++itI;
                        } else {
                            ++itI;
                        }
                    } else if (row == itRb.row()) {
                        if ((!itI) || ((itR) && (itR.row() < itI.row()))) {
                            tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()*itRb.value()));
                            ++itRb;
                            ++itR;
                        } else if ((itR) && (itI) && (itR.row() == itI.row())) {
                            tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()*itRb.value()));
                            tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itI.value()*itRb.value()));
                            ++itRb;
                            ++itR;
                            ++itI;
                        } else {
                            tripletListI.push_back(Triplet<mpreal>(itI.row(), k, itI.value()*itRb.value()));
                            ++itRb;
                            ++itI;
                        }
                    } else {
                        ++itRb;
                    }
                }
            }
        } else {
            for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
                SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
                SparseMatrix<mpreal>::InnerIterator itRb(b.matrixR,k);

                while ((itR) && (itRb)) {
                    if (itR.row() < itRb.row()) {
                        ++itR;
                    } else if (itR.row() == itRb.row()) {
                        tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()*itRb.value()));
                        ++itRb;
                        ++itR;
                    } else {
                        ++itRb;
                    }
                }
            }
        }

        // Now we can assign the data to the sparse matrix
        result.matrixR.setFromTriplets(tripletListR.begin(), tripletListR.end());
        if ((isComplex) || (b.isComplex))
            result.matrixI.setFromTriplets(tripletListI.begin(), tripletListI.end());

        // We may have set some coefficients to zero...
        result.matrixR.prune(0, 0);
        result.matrixR.makeCompressed();
        if ((isComplex) && (b.isComplex)) {
            result.matrixI.prune(0, 0);
            result.matrixI.makeCompressed();
        }
        result.checkComplexity();
    }

    return result;
}

/* The element-wise multiplication c = a.*b */
SparseGmpEigenMatrix& SparseGmpEigenMatrix::times_new(const SparseGmpEigenMatrix& b) const
{
    SparseGmpEigenMatrix& result(*(new SparseGmpEigenMatrix));

    // We want to support nxm matrix multiplied by 1x1 matrix (which we see as a scalar)
    if ((numel() != 1) && (b.numel() == 1)) {
        if (b.isComplex) {
            if (isComplex) {
                result.matrixR = matrixR*b.matrixR.coeff(0,0) - matrixI*b.matrixI.coeff(0,0);
                result.matrixI = matrixR*b.matrixI.coeff(0,0) + matrixI*b.matrixR.coeff(0,0);
            } else {
                result.matrixR = matrixR*b.matrixR.coeff(0,0);
                result.matrixI = matrixR*b.matrixI.coeff(0,0);
            }
        } else if (isComplex) {
            result.matrixR = matrixR*b.matrixR.coeff(0,0);
            result.matrixI = matrixI*b.matrixR.coeff(0,0);
        } else {
            result.matrixR = matrixR*b.matrixR.coeff(0,0);
        }

        // We may have set some coefficients to zero...
        result.matrixR.prune(0, 0);
        result.matrixR.makeCompressed();
        if ((isComplex) && (b.isComplex)) {
            result.matrixI.prune(0, 0);
            result.matrixI.makeCompressed();
        }
        result.checkComplexity();
    } else if ((numel() == 1) && (b.numel() != 1)) {
        if (b.isComplex) {
            if (isComplex) {
                result.matrixR = matrixR.coeff(0,0)*b.matrixR - matrixI.coeff(0,0)*b.matrixI;
                result.matrixI = matrixR.coeff(0,0)*b.matrixI + matrixI.coeff(0,0)*b.matrixR;
            } else {
                result.matrixR = matrixR.coeff(0,0)*b.matrixR;
                result.matrixI = matrixR.coeff(0,0)*b.matrixI;
            }
        } else if (isComplex) {
            result.matrixR = matrixR.coeff(0,0)*b.matrixR;
            result.matrixI = matrixI.coeff(0,0)*b.matrixR;
        } else {
            result.matrixR = matrixR.coeff(0,0)*b.matrixR;
        }

        // We may have set some coefficients to zero...
        result.matrixR.prune(0, 0);
        result.matrixR.makeCompressed();
        if ((isComplex) && (b.isComplex)) {
            result.matrixI.prune(0, 0);
            result.matrixI.makeCompressed();
        }
        result.checkComplexity();
    } else {
        // Eigen doesn't support this case, so we do it by hand...

        // setting the output size
        result.matrixR.resize(matrixR.rows(), matrixR.cols());
        if ((isComplex) || (b.isComplex))
            result.matrixI.resize(matrixR.rows(), matrixR.cols());

        // We will first copy the data to triplets
        vector< Triplet<mpreal> > tripletListR, tripletListI;
        tripletListR.reserve(std::min(matrixR.nonZeros(), b.matrixR.nonZeros()) + std::min(matrixI.nonZeros(), b.matrixI.nonZeros()));
        if ((isComplex) || (b.isComplex))
            tripletListI.reserve(std::min(matrixR.nonZeros(), b.matrixI.nonZeros()) + std::min(matrixI.nonZeros(), b.matrixR.nonZeros()));

        // Now for each column, we merge the lists of lines with nonzero elements of
        // both matrices. The cases are slightly different depending on the complexity
        // of both matrices
        if (b.isComplex) {
            if (isComplex) {
                for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
                    SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
                    SparseMatrix<mpreal>::InnerIterator itI(matrixI,k);
                    SparseMatrix<mpreal>::InnerIterator itRb(b.matrixR,k);
                    SparseMatrix<mpreal>::InnerIterator itIb(b.matrixI,k);

                    while (((itR) || (itI)) && ((itRb) || (itIb))) {
                        IndexType row, rowb;
                        if ((itR) && (itI))
                            row = min(itR.row(), itI.row());
                        else if (itR)
                            row = itR.row();
                        else
                            row = itI.row();
                        if ((itRb) && (itIb))
                            rowb = min(itRb.row(), itIb.row());
                        else if (itRb)
                            rowb = itRb.row();
                        else
                            rowb = itIb.row();

                        if (row < rowb) {
                            if ((!itI) || ((itR) && (itR.row() < itI.row()))) {
                                ++itR;
                            } else if ((itR) && (itI) && (itR.row() == itI.row())) {
                                ++itR;
                                ++itI;
                            } else {
                                ++itI;
                            }
                        } else if (row == rowb) {
                            if ((!itI) || ((itR) && (itR.row() < itI.row()))) {
                                if ((!itIb) || ((itRb) && (itRb.row() < itIb.row()))) {
                                    tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()*itRb.value()));
                                    ++itRb;
                                } else if ((itRb) && (itIb) && (itRb.row() == itIb.row())) {
                                    tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()*itRb.value()));
                                    tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itR.value()*itIb.value()));
                                    ++itRb;
                                    ++itIb;
                                } else {
                                    tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itR.value()*itIb.value()));
                                    ++itIb;
                                }
                                ++itR;
                            } else if ((itR) && (itI) && (itR.row() == itI.row())) {
                                if ((!itIb) || ((itRb) && (itRb.row() < itIb.row()))) {
                                    tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()*itRb.value()));
                                    tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itI.value()*itRb.value()));
                                    ++itRb;
                                } else if ((itRb) && (itIb) && (itRb.row() == itIb.row())) {
                                    tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()*itRb.value() - itI.value()*itIb.value()));
                                    tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itR.value()*itIb.value() + itI.value()*itRb.value()));
                                    ++itRb;
                                    ++itIb;
                                } else {
                                    tripletListR.push_back(Triplet<mpreal>(itR.row(), k, -itI.value()*itIb.value()));
                                    tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itR.value()*itIb.value()));
                                    ++itIb;
                                }
                                ++itR;
                                ++itI;
                            } else {
                                if ((!itIb) || ((itRb) && (itRb.row() < itIb.row()))) {
                                    tripletListI.push_back(Triplet<mpreal>(itI.row(), k, itI.value()*itRb.value()));
                                    ++itRb;
                                } else if ((itRb) && (itIb) && (itRb.row() == itIb.row())) {
                                    tripletListR.push_back(Triplet<mpreal>(itI.row(), k, -itI.value()*itIb.value()));
                                    tripletListI.push_back(Triplet<mpreal>(itI.row(), k, itI.value()*itRb.value()));
                                    ++itRb;
                                    ++itIb;
                                } else {
                                    tripletListR.push_back(Triplet<mpreal>(itI.row(), k, -itI.value()*itIb.value()));
                                    ++itIb;
                                }
                                ++itI;
                            }
                        } else {
                            if ((!itIb) || ((itRb) && (itRb.row() < itIb.row()))) {
                                ++itRb;
                            } else if ((itRb) && (itIb) && (itRb.row() == itIb.row())) {
                                ++itRb;
                                ++itIb;
                            } else {
                                ++itIb;
                            }
                        }
                    }
                }
            } else {
                for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
                    SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
                    SparseMatrix<mpreal>::InnerIterator itRb(b.matrixR,k);
                    SparseMatrix<mpreal>::InnerIterator itIb(b.matrixI,k);

                    while ((itR) && ((itRb) || (itIb))) {
                        IndexType rowb;
                        if ((itRb) && (itIb))
                            rowb = min(itRb.row(), itIb.row());
                        else if (itRb)
                            rowb = itRb.row();
                        else
                            rowb = itIb.row();

                        if (itR.row() < rowb) {
                            ++itR;
                        } else if (itR.row() == rowb) {
                            if ((!itIb) || ((itRb) && (itRb.row() < itIb.row()))) {
                                tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()*itRb.value()));
                                ++itRb;
                            } else if ((itRb) && (itIb) && (itRb.row() == itIb.row())) {
                                tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()*itRb.value()));
                                tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itR.value()*itIb.value()));
                                ++itRb;
                                ++itIb;
                            } else {
                                tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itR.value()*itIb.value()));
                                ++itIb;
                            }
                            ++itR;
                        } else {
                            if ((!itIb) || ((itRb) && (itRb.row() < itIb.row()))) {
                                ++itRb;
                            } else if ((itRb) && (itIb) && (itRb.row() == itIb.row())) {
                                ++itRb;
                                ++itIb;
                            } else {
                                ++itIb;
                            }
                        }
                    }
                }
            }
        } else if (isComplex) {
            for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
                SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
                SparseMatrix<mpreal>::InnerIterator itI(matrixI,k);
                SparseMatrix<mpreal>::InnerIterator itRb(b.matrixR,k);

                while (((itR) || (itI)) && (itRb)) {
                    IndexType row;
                    if ((itR) && (itI))
                        row = min(itR.row(), itI.row());
                    else if (itR)
                        row = itR.row();
                    else
                        row = itI.row();

                    if (row < itRb.row()) {
                        if ((!itI) || ((itR) && (itR.row() < itI.row()))) {
                            ++itR;
                        } else if ((itR) && (itI) && (itR.row() == itI.row())) {
                            ++itR;
                            ++itI;
                        } else {
                            ++itI;
                        }
                    } else if (row == itRb.row()) {
                        if ((!itI) || ((itR) && (itR.row() < itI.row()))) {
                            tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()*itRb.value()));
                            ++itRb;
                            ++itR;
                        } else if ((itR) && (itI) && (itR.row() == itI.row())) {
                            tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()*itRb.value()));
                            tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itI.value()*itRb.value()));
                            ++itRb;
                            ++itR;
                            ++itI;
                        } else {
                            tripletListI.push_back(Triplet<mpreal>(itI.row(), k, itI.value()*itRb.value()));
                            ++itRb;
                            ++itI;
                        }
                    } else {
                        ++itRb;
                    }
                }
            }
        } else {
            for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
                SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
                SparseMatrix<mpreal>::InnerIterator itRb(b.matrixR,k);

                while ((itR) && (itRb)) {
                    if (itR.row() < itRb.row()) {
                        ++itR;
                    } else if (itR.row() == itRb.row()) {
                        tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()*itRb.value()));
                        ++itRb;
                        ++itR;
                    } else {
                        ++itRb;
                    }
                }
            }
        }

        // Now we can assign the data to the sparse matrix
        result.matrixR.setFromTriplets(tripletListR.begin(), tripletListR.end());
        if ((isComplex) || (b.isComplex))
            result.matrixI.setFromTriplets(tripletListI.begin(), tripletListI.end());

        // We may have set some coefficients to zero...
        result.matrixR.prune(0, 0);
        result.matrixR.makeCompressed();
        if ((isComplex) && (b.isComplex)) {
            result.matrixI.prune(0, 0);
            result.matrixI.makeCompressed();
        }
        result.checkComplexity();
    }

    return result;
}

/* Element-wise multiplication with a full matrix c = a.*b */
SparseGmpEigenMatrix SparseGmpEigenMatrix::times_sf(const GmpEigenMatrix& b) const
{
    SparseGmpEigenMatrix result;

    // setting the output size
    result.matrixR.resize(matrixR.rows(), matrixR.cols());
    if ((isComplex) || (b.isComplex))
        result.matrixI.resize(matrixR.rows(), matrixR.cols());

    // We will first copy the data to triplets
    vector< Triplet<mpreal> > tripletListR, tripletListI;
    tripletListR.reserve(matrixR.nonZeros() + matrixI.nonZeros());
    if ((isComplex) || (b.isComplex))
        tripletListI.reserve(matrixR.nonZeros() + matrixI.nonZeros());

    // Now for each column, we merge the lists of lines with nonzero elements of
    // both matrices. The cases are slightly different depending on the complexity
    // of both matrices
    if (b.isComplex) {
        if (isComplex) {
            for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
                SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
                SparseMatrix<mpreal>::InnerIterator itI(matrixI,k);
                while ((itR) || (itI)) {
                    if ((!itI) || ((itR) && (itR.row() < itI.row()))) {
                        tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()*b.matrixR(itR.row(),k)));
                        tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itR.value()*b.matrixI(itR.row(),k)));
                        ++itR;
                    } else if ((itR) && (itI) && (itR.row() == itI.row())) {
                        tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()*b.matrixR(itR.row(),k) - itI.value()*b.matrixI(itR.row(),k)));
                        tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itR.value()*b.matrixI(itR.row(),k) + itI.value()*b.matrixR(itR.row(),k)));
                        ++itR;
                        ++itI;
                    } else {
                        tripletListR.push_back(Triplet<mpreal>(itI.row(), k, -itI.value()*b.matrixI(itI.row(),k)));
                        tripletListI.push_back(Triplet<mpreal>(itI.row(), k, itI.value()*b.matrixR(itI.row(),k)));
                        ++itI;
                    }
                }
            }
        } else {
            for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
                SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
                while (itR) {
                    tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()*b.matrixR(itR.row(),k)));
                    tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itR.value()*b.matrixI(itR.row(),k)));
                    ++itR;
                }
            }
        }
    } else if (isComplex) {
        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
            SparseMatrix<mpreal>::InnerIterator itI(matrixI,k);
            while ((itR) || (itI)) {
                if ((!itI) || ((itR) && (itR.row() < itI.row()))) {
                    tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()*b.matrixR(itR.row(),k)));
                    ++itR;
                } else if ((itR) && (itI) && (itR.row() == itI.row())) {
                    tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()*b.matrixR(itR.row(),k)));
                    tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itI.value()*b.matrixR(itR.row(),k)));
                    ++itR;
                    ++itI;
                } else {
                    tripletListI.push_back(Triplet<mpreal>(itI.row(), k, itI.value()*b.matrixR(itI.row(),k)));
                    ++itI;
                }
            }
        }
    } else {
        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
            while (itR) {
                tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()*b.matrixR(itR.row(),k)));
                ++itR;
            }
        }
    }

    // Now we can assign the data to the sparse matrix
    result.matrixR.setFromTriplets(tripletListR.begin(), tripletListR.end());
    if ((isComplex) || (b.isComplex))
        result.matrixI.setFromTriplets(tripletListI.begin(), tripletListI.end());

    // We may have set some coefficients to zero...
    result.matrixR.prune(0, 0);
    result.matrixR.makeCompressed();
    if ((isComplex) && (b.isComplex)) {
        result.matrixI.prune(0, 0);
        result.matrixI.makeCompressed();
    }
    result.checkComplexity();

    return result;
}

SparseGmpEigenMatrix& SparseGmpEigenMatrix::times_sf_new(const GmpEigenMatrix& b) const
{
    SparseGmpEigenMatrix& result(*(new SparseGmpEigenMatrix));

    // setting the output size
    result.matrixR.resize(matrixR.rows(), matrixR.cols());
    if ((isComplex) || (b.isComplex))
        result.matrixI.resize(matrixR.rows(), matrixR.cols());

    // We will first copy the data to triplets
    vector< Triplet<mpreal> > tripletListR, tripletListI;
    tripletListR.reserve(matrixR.nonZeros() + matrixI.nonZeros());
    if ((isComplex) || (b.isComplex))
        tripletListI.reserve(matrixR.nonZeros() + matrixI.nonZeros());

    // Now for each column, we merge the lists of lines with nonzero elements of
    // both matrices. The cases are slightly different depending on the complexity
    // of both matrices
    if (b.isComplex) {
        if (isComplex) {
            for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
                SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
                SparseMatrix<mpreal>::InnerIterator itI(matrixI,k);
                while ((itR) || (itI)) {
                    if ((!itI) || ((itR) && (itR.row() < itI.row()))) {
                        tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()*b.matrixR(itR.row(),k)));
                        tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itR.value()*b.matrixI(itR.row(),k)));
                        ++itR;
                    } else if ((itR) && (itI) && (itR.row() == itI.row())) {
                        tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()*b.matrixR(itR.row(),k) - itI.value()*b.matrixI(itR.row(),k)));
                        tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itR.value()*b.matrixI(itR.row(),k) + itI.value()*b.matrixR(itR.row(),k)));
                        ++itR;
                        ++itI;
                    } else {
                        tripletListR.push_back(Triplet<mpreal>(itI.row(), k, -itI.value()*b.matrixI(itI.row(),k)));
                        tripletListI.push_back(Triplet<mpreal>(itI.row(), k, itI.value()*b.matrixR(itI.row(),k)));
                        ++itI;
                    }
                }
            }
        } else {
            for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
                SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
                while (itR) {
                    tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()*b.matrixR(itR.row(),k)));
                    tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itR.value()*b.matrixI(itR.row(),k)));
                    ++itR;
                }
            }
        }
    } else if (isComplex) {
        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
            SparseMatrix<mpreal>::InnerIterator itI(matrixI,k);
            while ((itR) || (itI)) {
                if ((!itI) || ((itR) && (itR.row() < itI.row()))) {
                    tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()*b.matrixR(itR.row(),k)));
                    ++itR;
                } else if ((itR) && (itI) && (itR.row() == itI.row())) {
                    tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()*b.matrixR(itR.row(),k)));
                    tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itI.value()*b.matrixR(itR.row(),k)));
                    ++itR;
                    ++itI;
                } else {
                    tripletListI.push_back(Triplet<mpreal>(itI.row(), k, itI.value()*b.matrixR(itI.row(),k)));
                    ++itI;
                }
            }
        }
    } else {
        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
            while (itR) {
                tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()*b.matrixR(itR.row(),k)));
                ++itR;
            }
        }
    }

    // Now we can assign the data to the sparse matrix
    result.matrixR.setFromTriplets(tripletListR.begin(), tripletListR.end());
    if ((isComplex) || (b.isComplex))
        result.matrixI.setFromTriplets(tripletListI.begin(), tripletListI.end());

    // We may have set some coefficients to zero...
    result.matrixR.prune(0, 0);
    result.matrixR.makeCompressed();
    if ((isComplex) && (b.isComplex)) {
        result.matrixI.prune(0, 0);
        result.matrixI.makeCompressed();
    }
    result.checkComplexity();

    return result;
}


/* The element-wise division on the right c = a./b */
SparseGmpEigenMatrix SparseGmpEigenMatrix::rdivide(const SparseGmpEigenMatrix& b) const
{
    SparseGmpEigenMatrix result;

    // We only support nxm matrix divided by 1x1 matrix (which we see as a scalar)
    if (b.isComplex) {
        if (isComplex) {
            //result.matrixR = (matrixR*b.matrixR.coeff(0,0) + matrixI*b.matrixI.coeff(0,0))/(pow(b.matrixR.coeff(0,0),2) + pow(b.matrixI.coeff(0,0),2)); // Not supported anymore on Eigen 3.3
            //result.matrixI = (-matrixR*b.matrixI.coeff(0,0) + matrixI*b.matrixR.coeff(0,0))/(pow(b.matrixR.coeff(0,0),2) + pow(b.matrixI.coeff(0,0),2));
            result.matrixR = (matrixR*b.matrixR.coeff(0,0) + matrixI*b.matrixI.coeff(0,0));
            result.matrixR = result.matrixR/(pow(b.matrixR.coeff(0,0),2) + pow(b.matrixI.coeff(0,0),2));
            result.matrixI = (-matrixR*b.matrixI.coeff(0,0) + matrixI*b.matrixR.coeff(0,0));
            result.matrixI = result.matrixI/(pow(b.matrixR.coeff(0,0),2) + pow(b.matrixI.coeff(0,0),2));
        } else {
            //result.matrixR = matrixR*b.matrixR.coeff(0,0)/(pow(b.matrixR.coeff(0,0),2) + pow(b.matrixI.coeff(0,0),2)); // Not supported anymore on Eigen 3.3
            //result.matrixI = -matrixR*b.matrixI.coeff(0,0)/(pow(b.matrixR.coeff(0,0),2) + pow(b.matrixI.coeff(0,0),2));
            result.matrixR = matrixR*b.matrixR.coeff(0,0);
            result.matrixR = result.matrixR/(pow(b.matrixR.coeff(0,0),2) + pow(b.matrixI.coeff(0,0),2));
            result.matrixI = -matrixR*b.matrixI.coeff(0,0);
            result.matrixI = result.matrixI/(pow(b.matrixR.coeff(0,0),2) + pow(b.matrixI.coeff(0,0),2));
        }
    } else if (isComplex) {
        result.matrixR = matrixR/b.matrixR.coeff(0,0);
        result.matrixI = matrixI/b.matrixR.coeff(0,0);
    } else {
        result.matrixR = matrixR/b.matrixR.coeff(0,0);
    }

    // We may have set some coefficients to zero...
    result.matrixR.prune(0, 0);
    result.matrixR.makeCompressed();
    if ((isComplex) || (b.isComplex)) {
        result.matrixI.prune(0, 0);
        result.matrixI.makeCompressed();
    }
    result.checkComplexity();

    return result;
}

/* The element-wise division on the right c = a./b */
SparseGmpEigenMatrix& SparseGmpEigenMatrix::rdivide_new(const SparseGmpEigenMatrix& b) const
{
    SparseGmpEigenMatrix& result(*(new SparseGmpEigenMatrix));

    // We only support nxm matrix divided by 1x1 matrix (which we see as a scalar)
    if (b.isComplex) {
        if (isComplex) {
            //result.matrixR = (matrixR*b.matrixR.coeff(0,0) + matrixI*b.matrixI.coeff(0,0))/(pow(b.matrixR.coeff(0,0),2) + pow(b.matrixI.coeff(0,0),2)); // Not supported anymore on Eigen 3.3
            //result.matrixI = (-matrixR*b.matrixI.coeff(0,0) + matrixI*b.matrixR.coeff(0,0))/(pow(b.matrixR.coeff(0,0),2) + pow(b.matrixI.coeff(0,0),2));
            result.matrixR = (matrixR*b.matrixR.coeff(0,0) + matrixI*b.matrixI.coeff(0,0));
            result.matrixR = result.matrixR/(pow(b.matrixR.coeff(0,0),2) + pow(b.matrixI.coeff(0,0),2));
            result.matrixI = (-matrixR*b.matrixI.coeff(0,0) + matrixI*b.matrixR.coeff(0,0));
            result.matrixI = result.matrixI/(pow(b.matrixR.coeff(0,0),2) + pow(b.matrixI.coeff(0,0),2));
        } else {
            //result.matrixR = matrixR*b.matrixR.coeff(0,0)/(pow(b.matrixR.coeff(0,0),2) + pow(b.matrixI.coeff(0,0),2)); // Not supported anymore on Eigen 3.3
            //result.matrixI = -matrixR*b.matrixI.coeff(0,0)/(pow(b.matrixR.coeff(0,0),2) + pow(b.matrixI.coeff(0,0),2));
            result.matrixR = matrixR*b.matrixR.coeff(0,0);
            result.matrixR = result.matrixR/(pow(b.matrixR.coeff(0,0),2) + pow(b.matrixI.coeff(0,0),2));
            result.matrixI = -matrixR*b.matrixI.coeff(0,0);
            result.matrixI = result.matrixI/(pow(b.matrixR.coeff(0,0),2) + pow(b.matrixI.coeff(0,0),2));
        }
    } else if (isComplex) {
        result.matrixR = matrixR/b.matrixR.coeff(0,0);
        result.matrixI = matrixI/b.matrixR.coeff(0,0);
    } else {
        result.matrixR = matrixR/b.matrixR.coeff(0,0);
    }

    // We may have set some coefficients to zero...
    result.matrixR.prune(0, 0);
    result.matrixR.makeCompressed();
    if ((isComplex) || (b.isComplex)) {
        result.matrixI.prune(0, 0);
        result.matrixI.makeCompressed();
    }
    result.checkComplexity();

    return result;
}


/* The absolute value (or complex magnitude) b = abs(a) */
// Note that this function is needed to define the operation power, so we
// don't call that function here.
SparseGmpEigenMatrix SparseGmpEigenMatrix::abs() const
{
    SparseGmpEigenMatrix result;

    result.isComplex = false;
    result.matrixR.resize(matrixR.rows(), matrixR.cols());
    if (!isComplex) {
        // This case is simpler
        result.matrixR.reserve(matrixR.nonZeros());
        for (IndexType k = 0; k < matrixR.outerSize(); ++k)
            for (SparseMatrix<mpreal>::InnerIterator it(matrixR,k); it; ++it)
                result.matrixR.insert(it.row(), it.col()) = mpfr::abs(it.value());
    } else {
        // For each column, we should merge the lists of lines with nonzero elements...
        result.matrixR.reserve(min(matrixR.rows()*matrixR.cols(), matrixR.nonZeros() + matrixI.nonZeros()));
        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
            SparseMatrix<mpreal>::InnerIterator itI(matrixI,k);

            while ((itR) || (itI)) {
                if ((itR) && (itI)) {
                    if (itR.row() < itI.row()) {
                        result.matrixR.insert(itR.row(), itR.col()) = mpfr::abs(itR.value());
                        ++itR;
                    } else if (itR.row() == itI.row()) {
                        result.matrixR.insert(itR.row(), itR.col()) = mpfr::sqrt((pow(itR.value(), 2) + pow(itI.value(), 2)));
                        ++itR;
                        ++itI;
                    } else {
                        result.matrixR.insert(itI.row(), itI.col()) = mpfr::abs(itI.value());
                        ++itI;
                    }
                } else if (itR) {
                    result.matrixR.insert(itR.row(), itR.col()) = mpfr::abs(itR.value());
                    ++itR;
                } else {
                    result.matrixR.insert(itI.row(), itI.col()) = mpfr::abs(itI.value());
                    ++itI;
                }
            }
        }
        result.matrixR.makeCompressed();
    }

    return result;
}

/* The absolute value (or complex magnitude) b = abs(a) */
SparseGmpEigenMatrix& SparseGmpEigenMatrix::abs_new() const
{
    SparseGmpEigenMatrix& result(*(new SparseGmpEigenMatrix));

    result.isComplex = false;
    result.matrixR.resize(matrixR.rows(), matrixR.cols());
    if (!isComplex) {
        // This case is simpler
        result.matrixR.reserve(matrixR.nonZeros());
        for (IndexType k = 0; k < matrixR.outerSize(); ++k)
            for (SparseMatrix<mpreal>::InnerIterator it(matrixR,k); it; ++it)
                result.matrixR.insert(it.row(), it.col()) = mpfr::abs(it.value());
    } else {
        // For each column, we should merge the lists of lines with nonzero elements...
        result.matrixR.reserve(min(matrixR.rows()*matrixR.cols(), matrixR.nonZeros() + matrixI.nonZeros()));
        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
            SparseMatrix<mpreal>::InnerIterator itI(matrixI,k);

            while ((itR) || (itI)) {
                if ((itR) && (itI)) {
                    if (itR.row() < itI.row()) {
                        result.matrixR.insert(itR.row(), itR.col()) = mpfr::abs(itR.value());
                        ++itR;
                    } else if (itR.row() == itI.row()) {
                        result.matrixR.insert(itR.row(), itR.col()) = mpfr::sqrt((pow(itR.value(), 2) + pow(itI.value(), 2)));
                        ++itR;
                        ++itI;
                    } else {
                        result.matrixR.insert(itI.row(), itI.col()) = mpfr::abs(itI.value());
                        ++itI;
                    }
                } else if (itR) {
                    result.matrixR.insert(itR.row(), itR.col()) = mpfr::abs(itR.value());
                    ++itR;
                } else {
                    result.matrixR.insert(itI.row(), itI.col()) = mpfr::abs(itI.value());
                    ++itI;
                }
            }
        }
        result.matrixR.makeCompressed();
    }

    return result;
}


/* The phase angle (i.e. complex argument) b = angle(a) */
SparseGmpEigenMatrix SparseGmpEigenMatrix::angle() const
{
    SparseGmpEigenMatrix result;

    result.isComplex = false;
    result.matrixR.resize(matrixR.rows(), matrixR.cols());
    if (!isComplex) {
        // This case is simpler
        result.matrixR.reserve(matrixR.nonZeros());
        for (IndexType k = 0; k < matrixR.outerSize(); ++k)
            for (SparseMatrix<mpreal>::InnerIterator it(matrixR,k); it; ++it)
                if ((it.value() < 0) || (mpfr::isnan(it.value())))
                    result.matrixR.insert(it.row(), it.col()) = atan2(mpreal(0), it.value());
    } else {
        // For each column, we should merge the lists of lines with nonzero elements...
        result.matrixR.reserve(min(matrixR.rows()*matrixR.cols(), matrixR.nonZeros() + matrixI.nonZeros()));
        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
            SparseMatrix<mpreal>::InnerIterator itI(matrixI,k);

            while ((itR) || (itI)) {
                if ((itR) && (itI)) {
                    if (itR.row() < itI.row()) {
                        if ((itR.value() < 0) || (mpfr::isnan(itR.value())))
                            result.matrixR.insert(itR.row(), itR.col()) = atan2(mpreal(0), itR.value());
                        ++itR;
                    } else if (itR.row() == itI.row()) {
                        result.matrixR.insert(itR.row(), itR.col()) = atan2(itI.value(), itR.value());
                        ++itR;
                        ++itI;
                    } else {
                        result.matrixR.insert(itI.row(), itI.col()) = atan2(itI.value(), mpreal(0));
                        ++itI;
                    }
                } else if (itR) {
                    if ((itR.value() < 0) || (mpfr::isnan(itR.value())))
                        result.matrixR.insert(itR.row(), itR.col()) = atan2(mpreal(0), itR.value());
                    ++itR;
                } else {
                    result.matrixR.insert(itI.row(), itI.col()) = atan2(itI.value(), mpreal(0));
                    ++itI;
                }
            }
        }
    }
    result.matrixR.makeCompressed();

    return result;
}

/* The phase angle (i.e. complex argument) b = angle(a) */
SparseGmpEigenMatrix& SparseGmpEigenMatrix::angle_new() const
{
    SparseGmpEigenMatrix& result(*(new SparseGmpEigenMatrix));

    result.isComplex = false;
    result.matrixR.resize(matrixR.rows(), matrixR.cols());
    if (!isComplex) {
        // This case is simpler
        result.matrixR.reserve(matrixR.nonZeros());
        for (IndexType k = 0; k < matrixR.outerSize(); ++k)
            for (SparseMatrix<mpreal>::InnerIterator it(matrixR,k); it; ++it)
                if ((it.value() < 0) || (mpfr::isnan(it.value())))
                    result.matrixR.insert(it.row(), it.col()) = atan2(mpreal(0), it.value());
    } else {
        // For each column, we should merge the lists of lines with nonzero elements...
        result.matrixR.reserve(min(matrixR.rows()*matrixR.cols(), matrixR.nonZeros() + matrixI.nonZeros()));
        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
            SparseMatrix<mpreal>::InnerIterator itI(matrixI,k);

            while ((itR) || (itI)) {
                if ((itR) && (itI)) {
                    if (itR.row() < itI.row()) {
                        if ((itR.value() < 0) || (mpfr::isnan(itR.value())))
                            result.matrixR.insert(itR.row(), itR.col()) = atan2(mpreal(0), itR.value());
                        ++itR;
                    } else if (itR.row() == itI.row()) {
                        result.matrixR.insert(itR.row(), itR.col()) = atan2(itI.value(), itR.value());
                        ++itR;
                        ++itI;
                    } else {
                        result.matrixR.insert(itI.row(), itI.col()) = atan2(itI.value(), mpreal(0));
                        ++itI;
                    }
                } else if (itR) {
                    if ((itR.value() < 0) || (mpfr::isnan(itR.value())))
                        result.matrixR.insert(itR.row(), itR.col()) = atan2(mpreal(0), itR.value());
                    ++itR;
                } else {
                    result.matrixR.insert(itI.row(), itI.col()) = atan2(itI.value(), mpreal(0));
                    ++itI;
                }
            }
        }
    }
    result.matrixR.makeCompressed();

    return result;
}



/* The element-wise exponential b = exp(a) */
SparseGmpEigenMatrix SparseGmpEigenMatrix::exp() const
{
    SparseGmpEigenMatrix result;

    if (!isComplex) {
        // This case is simpler
        result.isComplex = false;
        // exp(0) is :
        result.matrixR = sparseConstantMatrix(matrixR.rows(), matrixR.cols(), 1);
        for (IndexType k = 0; k < matrixR.outerSize(); ++k)
            for (SparseMatrix<mpreal>::InnerIterator it(matrixR,k); it; ++it)
                if (it.value() == mpreal("-Inf"))
                    result.matrixR.coeffRef(it.row(), it.col()) = 0;
                else
                    result.matrixR.coeffRef(it.row(), it.col()) = mpfr::exp(it.value());
        result.matrixR.prune(0, 0);
        result.matrixR.makeCompressed();
    } else {
        // For each column, we should merge the lists of lines with nonzero elements...
        // exp(0) is :
        result.matrixR = sparseConstantMatrix(matrixR.rows(), matrixR.cols(), 1);
        result.matrixI.resize(matrixR.rows(), matrixR.cols());
        result.matrixI.reserve(min(matrixR.rows()*matrixR.cols(), matrixR.nonZeros() + matrixI.nonZeros()));
        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
            SparseMatrix<mpreal>::InnerIterator itI(matrixI,k);

            while ((itR) || (itI)) {
                if ((itR) && (itI)) {
                    if (itR.row() < itI.row()) {
                        if (itR.value() == mpreal("-Inf"))
                            result.matrixR.coeffRef(itR.row(), itR.col()) = 0;
                        else
                            result.matrixR.coeffRef(itR.row(), itR.col()) = mpfr::exp(itR.value());
                        ++itR;
                    } else if (itR.row() == itI.row()) {
                        result.matrixR.coeffRef(itR.row(), itR.col()) = mpfr::exp(itR.value())*mpfr::cos(itI.value());
                        result.matrixI.insert(itR.row(), itR.col()) = mpfr::exp(itR.value())*mpfr::sin(itI.value());
                        ++itR;
                        ++itI;
                    } else {
                        result.matrixR.coeffRef(itI.row(), itI.col()) = mpfr::cos(itI.value());
                        result.matrixI.insert(itI.row(), itI.col()) = mpfr::sin(itI.value());
                        ++itI;
                    }
                } else if (itR) {
                    if (itR.value() == mpreal("-Inf"))
                        result.matrixR.coeffRef(itR.row(), itR.col()) = 0;
                    else
                        result.matrixR.coeffRef(itR.row(), itR.col()) = mpfr::exp(itR.value());
                    ++itR;
                } else {
                    result.matrixR.coeffRef(itI.row(), itI.col()) = mpfr::cos(itI.value());
                    result.matrixI.insert(itI.row(), itI.col()) = mpfr::sin(itI.value());
                    ++itI;
                }
            }
        }
        // We may have created some new zero coefficients...
        result.matrixR.prune(0, 0);
        result.matrixI.prune(0, 0);

        result.matrixR.makeCompressed();
        result.matrixI.makeCompressed();
        result.checkComplexity();
    }

    return result;
}

/* The element-wise exponential b = exp(a) */
SparseGmpEigenMatrix& SparseGmpEigenMatrix::exp_new() const
{
    SparseGmpEigenMatrix& result(*(new SparseGmpEigenMatrix));

    if (!isComplex) {
        // This case is simpler
        result.isComplex = false;
        // exp(0) is :
        result.matrixR = sparseConstantMatrix(matrixR.rows(), matrixR.cols(), 1);
        for (IndexType k = 0; k < matrixR.outerSize(); ++k)
            for (SparseMatrix<mpreal>::InnerIterator it(matrixR,k); it; ++it)
                if (it.value() == mpreal("-Inf"))
                    result.matrixR.coeffRef(it.row(), it.col()) = 0;
                else
                    result.matrixR.coeffRef(it.row(), it.col()) = mpfr::exp(it.value());
        result.matrixR.prune(0, 0);
        result.matrixR.makeCompressed();
    } else {
        // For each column, we should merge the lists of lines with nonzero elements...
        // exp(0) is :
        result.matrixR = sparseConstantMatrix(matrixR.rows(), matrixR.cols(), 1);
        result.matrixI.resize(matrixR.rows(), matrixR.cols());
        result.matrixI.reserve(min(matrixR.rows()*matrixR.cols(), matrixR.nonZeros() + matrixI.nonZeros()));
        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
            SparseMatrix<mpreal>::InnerIterator itI(matrixI,k);

            while ((itR) || (itI)) {
                if ((itR) && (itI)) {
                    if (itR.row() < itI.row()) {
                        if (itR.value() == mpreal("-Inf"))
                            result.matrixR.coeffRef(itR.row(), itR.col()) = 0;
                        else
                            result.matrixR.coeffRef(itR.row(), itR.col()) = mpfr::exp(itR.value());
                        ++itR;
                    } else if (itR.row() == itI.row()) {
                        result.matrixR.coeffRef(itR.row(), itR.col()) = mpfr::exp(itR.value())*mpfr::cos(itI.value());
                        result.matrixI.insert(itR.row(), itR.col()) = mpfr::exp(itR.value())*mpfr::sin(itI.value());
                        ++itR;
                        ++itI;
                    } else {
                        result.matrixR.coeffRef(itI.row(), itI.col()) = mpfr::cos(itI.value());
                        result.matrixI.insert(itI.row(), itI.col()) = mpfr::sin(itI.value());
                        ++itI;
                    }
                } else if (itR) {
                    if (itR.value() == mpreal("-Inf"))
                        result.matrixR.coeffRef(itR.row(), itR.col()) = 0;
                    else
                        result.matrixR.coeffRef(itR.row(), itR.col()) = mpfr::exp(itR.value());
                    ++itR;
                } else {
                    result.matrixR.coeffRef(itI.row(), itI.col()) = mpfr::cos(itI.value());
                    result.matrixI.insert(itI.row(), itI.col()) = mpfr::sin(itI.value());
                    ++itI;
                }
            }
        }
        // We may have created some new zero coefficients...
        result.matrixR.prune(0, 0);
        result.matrixI.prune(0, 0);

        result.matrixR.makeCompressed();
        result.matrixI.makeCompressed();
        result.checkComplexity();
    }

    return result;
}

/* The element-wise natural logarithm b = log(a) */
SparseGmpEigenMatrix SparseGmpEigenMatrix::log() const
{
    SparseGmpEigenMatrix result;

    SparseGmpEigenMatrix tmp(SparseGmpEigenMatrix::abs());
    // log(0) is -Inf
    result.matrixR = sparseConstantMatrix(matrixR.rows(), matrixR.cols(), mpreal("-Inf"));
    for (IndexType k = 0; k < tmp.matrixR.outerSize(); ++k)
        for (SparseMatrix<mpreal>::InnerIterator it(tmp.matrixR,k); it; ++it)
            if (it.value() == mpreal(1))
                result.matrixR.coeffRef(it.row(), it.col()) = 0;
            else
                result.matrixR.coeffRef(it.row(), it.col()) = mpfr::log(it.value());
    result.matrixR.prune(0, 0);
    result.matrixR.makeCompressed();

    tmp = SparseGmpEigenMatrix::angle();
    result.matrixI = tmp.matrixR;

    result.checkComplexity();

    return result;
}

/* The element-wise natural logarithm b = log(a) */
SparseGmpEigenMatrix& SparseGmpEigenMatrix::log_new() const
{
    SparseGmpEigenMatrix& result(*(new SparseGmpEigenMatrix));

    SparseGmpEigenMatrix tmp(SparseGmpEigenMatrix::abs());
    // log(0) is -Inf
    result.matrixR = sparseConstantMatrix(matrixR.rows(), matrixR.cols(), mpreal("-Inf"));
    for (IndexType k = 0; k < tmp.matrixR.outerSize(); ++k)
        for (SparseMatrix<mpreal>::InnerIterator it(tmp.matrixR,k); it; ++it)
            if (it.value() == mpreal(1))
                result.matrixR.coeffRef(it.row(), it.col()) = 0;
            else
                result.matrixR.coeffRef(it.row(), it.col()) = mpfr::log(it.value());
    result.matrixR.prune(0, 0);
    result.matrixR.makeCompressed();

    tmp = SparseGmpEigenMatrix::angle();
    result.matrixI = tmp.matrixR;

    result.checkComplexity();

    return result;
}


/* The element-wise power c = a.^b */
// Here, b is real ans b>0
SparseGmpEigenMatrix SparseGmpEigenMatrix::power(const GmpEigenMatrix& b) const
{
    SparseGmpEigenMatrix result(*this);

    if (isComplex) {
        if (b.numel() == 1) {
            for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
                SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
                SparseMatrix<mpreal>::InnerIterator itI(matrixI,k);

                IndexType row, col;
                while ((itR) || (itI)) {
                    if ((itR) && (itI)) {
                        if (itR.row() < itI.row()) {
                            row = itR.row();
                            col = itR.col();
                            ++itR;
                        } else if (itR.row() == itI.row()) {
                            row = itR.row();
                            col = itR.col();
                            ++itR;
                            ++itI;
                        } else {
                            row = itI.row();
                            col = itI.col();
                            ++itI;
                        }
                    } else if (itR) {
                        row = itR.row();
                        col = itR.col();
                        ++itR;
                    } else {
                        row = itI.row();
                        col = itI.col();
                        ++itI;
                    }
                    // We extract the element of interest and compute its power
                    SparseGmpEigenMatrix tmp(block(row, col, 1, 1));
                    tmp = tmp.log().times_sf(b).exp();
                    result.matrixR.coeffRef(row,col) = tmp.matrixI.coeff(0,0);
                    if (tmp.isComplex)
                        result.matrixI.coeffRef(row,col) = tmp.matrixR.coeff(0,0);
                    else
                        result.matrixI.coeffRef(row,col) = 0;
                }
            }
        } else {
            for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
                SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
                SparseMatrix<mpreal>::InnerIterator itI(matrixI,k);

                IndexType row, col;
                while ((itR) || (itI)) {
                    if ((itR) && (itI)) {
                        if (itR.row() < itI.row()) {
                            row = itR.row();
                            col = itR.col();
                            ++itR;
                        } else if (itR.row() == itI.row()) {
                            row = itR.row();
                            col = itR.col();
                            ++itR;
                            ++itI;
                        } else {
                            row = itI.row();
                            col = itI.col();
                            ++itI;
                        }
                    } else if (itR) {
                        row = itR.row();
                        col = itR.col();
                        ++itR;
                    } else {
                        row = itI.row();
                        col = itI.col();
                        ++itI;
                    }
                    // We extract the element of interest and compute its power
                    SparseGmpEigenMatrix tmp(block(row, col, 1, 1));
                    tmp = tmp.log().times_sf(b.block(row,col,1,1)).exp();
                    result.matrixR.coeffRef(row,col) = tmp.matrixR.coeff(0,0);
                    if (tmp.isComplex)
                        result.matrixI.coeffRef(row,col) = tmp.matrixI.coeff(0,0);
                    else
                        result.matrixI.coeffRef(row,col) = 0;
                }
            }
        }
        result.matrixR.prune(0, 0);
        result.matrixR.makeCompressed();
        result.matrixI.prune(0, 0);
        result.matrixI.makeCompressed();
        result.checkComplexity();
    } else {
        result.isComplex = false;
        if (b.numel() == 1) {
            for (IndexType k = 0; k < matrixR.outerSize(); ++k)
                for (SparseMatrix<mpreal>::InnerIterator it(matrixR,k); it; ++it)
                    result.matrixR.coeffRef(it.row(),it.col()) = pow(it.value(), b.matrixR(0,0));
        } else {
            for (IndexType k = 0; k < matrixR.outerSize(); ++k)
                for (SparseMatrix<mpreal>::InnerIterator it(matrixR,k); it; ++it)
                    result.matrixR.coeffRef(it.row(),it.col()) = pow(it.value(), b.matrixR(it.row(),it.col()));
        }
        result.matrixR.prune(0, 0);
        result.matrixR.makeCompressed();
    }

    return result;
}

/* The element-wise power c = a.^b */
SparseGmpEigenMatrix& SparseGmpEigenMatrix::power_new(const GmpEigenMatrix& b) const
{
    SparseGmpEigenMatrix& result(*(new SparseGmpEigenMatrix(*this)));

    if (isComplex) {
        if (b.numel() == 1) {
            for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
                SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
                SparseMatrix<mpreal>::InnerIterator itI(matrixI,k);

                IndexType row, col;
                while ((itR) || (itI)) {
                    if ((itR) && (itI)) {
                        if (itR.row() < itI.row()) {
                            row = itR.row();
                            col = itR.col();
                            ++itR;
                        } else if (itR.row() == itI.row()) {
                            row = itR.row();
                            col = itR.col();
                            ++itR;
                            ++itI;
                        } else {
                            row = itI.row();
                            col = itI.col();
                            ++itI;
                        }
                    } else if (itR) {
                        row = itR.row();
                        col = itR.col();
                        ++itR;
                    } else {
                        row = itI.row();
                        col = itI.col();
                        ++itI;
                    }
                    // We extract the element of interest and compute its power
                    SparseGmpEigenMatrix tmp(block(row, col, 1, 1));
                    tmp = tmp.log().times_sf(b).exp();
                    result.matrixR.coeffRef(row,col) = tmp.matrixR.coeff(0,0);
                    if (tmp.isComplex)
                        result.matrixI.coeffRef(row,col) = tmp.matrixI.coeff(0,0);
                    else
                        result.matrixI.coeffRef(row,col) = 0;
                }
            }
        } else {
            for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
                SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
                SparseMatrix<mpreal>::InnerIterator itI(matrixI,k);

                IndexType row, col;
                while ((itR) || (itI)) {
                    if ((itR) && (itI)) {
                        if (itR.row() < itI.row()) {
                            row = itR.row();
                            col = itR.col();
                            ++itR;
                        } else if (itR.row() == itI.row()) {
                            row = itR.row();
                            col = itR.col();
                            ++itR;
                            ++itI;
                        } else {
                            row = itI.row();
                            col = itI.col();
                            ++itI;
                        }
                    } else if (itR) {
                        row = itR.row();
                        col = itR.col();
                        ++itR;
                    } else {
                        row = itI.row();
                        col = itI.col();
                        ++itI;
                    }
                    // We extract the element of interest and compute its power
                    SparseGmpEigenMatrix tmp(block(row, col, 1, 1));
                    tmp = tmp.log().times_sf(b.block(row,col,1,1)).exp();
                    result.matrixR.coeffRef(row,col) = tmp.matrixR.coeff(0,0);
                    if (tmp.isComplex)
                        result.matrixI.coeffRef(row,col) = tmp.matrixI.coeff(0,0);
                    else
                        result.matrixI.coeffRef(row,col) = 0;
                }
            }
        }
        result.matrixR.prune(0, 0);
        result.matrixR.makeCompressed();
        result.matrixI.prune(0, 0);
        result.matrixI.makeCompressed();
        result.checkComplexity();
    } else {
        if (b.numel() == 1) {
            for (IndexType k = 0; k < matrixR.outerSize(); ++k)
                for (SparseMatrix<mpreal>::InnerIterator it(matrixR,k); it; ++it)
                    result.matrixR.coeffRef(it.row(),it.col()) = pow(it.value(), b.matrixR(0,0));
        } else {
            for (IndexType k = 0; k < matrixR.outerSize(); ++k)
                for (SparseMatrix<mpreal>::InnerIterator it(matrixR,k); it; ++it)
                    result.matrixR.coeffRef(it.row(),it.col()) = pow(it.value(), b.matrixR(it.row(),it.col()));
        }
        result.matrixR.prune(0, 0);
        result.matrixR.makeCompressed();
    }

    return result;
}


/* square root b = sqrt(a) */
SparseGmpEigenMatrix SparseGmpEigenMatrix::sqrt() const
{
    SparseGmpEigenMatrix result;

    mpreal globalMin(0);
    if (!isComplex) {
        // In this case we need to compute the smallest element in the matrix
        globalMin = globalMinReal();
    }

    if ((!isComplex) && (globalMin >= 0)) {
        // This case is simpler
        result.isComplex = false;
        result.matrixR.resize(matrixR.rows(), matrixR.cols());
        result.matrixR.reserve(matrixR.nonZeros());
        for (IndexType k = 0; k < matrixR.outerSize(); ++k)
            for (SparseMatrix<mpreal>::InnerIterator it(matrixR,k); it; ++it)
                result.matrixR.insert(it.row(), it.col()) = mpfr::sqrt(it.value());
    } else {
        // For each column, we should merge the lists of lines with nonzero elements...
        SparseGmpEigenMatrix norm(abs());
        result.matrixR.resize(matrixR.rows(), matrixR.cols());
        result.matrixR.reserve(norm.matrixR.nonZeros());
        result.matrixI.resize(matrixR.rows(), matrixR.cols());
        result.matrixI.reserve(norm.matrixR.nonZeros());
        for (IndexType k = 0; k < norm.matrixR.outerSize(); ++k)
            for (SparseMatrix<mpreal>::InnerIterator it(norm.matrixR,k); it; ++it) {
                // Re(sqrt(-Inf+c*1i)) = 0
                if ((matrixR.coeff(it.row(), it.col()) == mpreal("-Inf")) && (isfinite(matrixI.coeff(it.row(), it.col()))))
                    ;
                else
                    result.matrixR.insert(it.row(), it.col()) = mpfr::sqrt((it.value() + matrixR.coeff(it.row(), it.col()))/2);
                // Im(sqrt(NaN)) = 0 (by convention)
                // and
                // Im(sqrt(Inf+c*1i)) = 0
                if ((mpfr::isnan(matrixR.coeff(it.row(), it.col()))) && (matrixI.coeff(it.row(), it.col()) == 0))
                    ;
                else if ((matrixR.coeff(it.row(), it.col()) == mpreal("Inf")) && (isfinite(matrixI.coeff(it.row(), it.col()))))
                    ;
                else
                    result.matrixI.insert(it.row(), it.col()) = ((isComplex) && (matrixI.coeff(it.row(), it.col()) < 0) ? -mpfr::sqrt((it.value() - matrixR.coeff(it.row(), it.col()))/2) : mpfr::sqrt((it.value() - matrixR.coeff(it.row(), it.col()))/2));
            }

        // We may have created some new zero coefficients...
        result.matrixR.prune(0, 0);
        result.matrixI.prune(0, 0);

        result.matrixR.makeCompressed();
        result.matrixI.makeCompressed();
        result.checkComplexity();
    }

    return result;
}

/* square root b = sqrt(a) */
SparseGmpEigenMatrix& SparseGmpEigenMatrix::sqrt_new() const
{
    SparseGmpEigenMatrix& result(*(new SparseGmpEigenMatrix));

    mpreal globalMin(0);
    if (!isComplex) {
        // In this case we need to compute the smallest element in the matrix
        globalMin = globalMinReal();
    }

    if ((!isComplex) && (globalMin >= 0)) {
        // This case is simpler
        result.isComplex = false;
        result.matrixR.resize(matrixR.rows(), matrixR.cols());
        result.matrixR.reserve(matrixR.nonZeros());
        for (IndexType k = 0; k < matrixR.outerSize(); ++k)
            for (SparseMatrix<mpreal>::InnerIterator it(matrixR,k); it; ++it)
                result.matrixR.insert(it.row(), it.col()) = mpfr::sqrt(it.value());
    } else {
        // For each column, we should merge the lists of lines with nonzero elements...
        SparseGmpEigenMatrix norm(abs());
        result.matrixR.resize(matrixR.rows(), matrixR.cols());
        result.matrixR.reserve(norm.matrixR.nonZeros());
        result.matrixI.resize(matrixR.rows(), matrixR.cols());
        result.matrixI.reserve(norm.matrixR.nonZeros());
        for (IndexType k = 0; k < norm.matrixR.outerSize(); ++k)
            for (SparseMatrix<mpreal>::InnerIterator it(norm.matrixR,k); it; ++it) {
                // Re(sqrt(-Inf+c*1i)) = 0
                if ((matrixR.coeff(it.row(), it.col()) == mpreal("-Inf")) && (isfinite(matrixI.coeff(it.row(), it.col()))))
                    ;
                else
                    result.matrixR.insert(it.row(), it.col()) = mpfr::sqrt((it.value() + matrixR.coeff(it.row(), it.col()))/2);
                // Im(sqrt(NaN)) = 0 (by convention)
                // and
                // Im(sqrt(Inf+c*1i)) = 0
                if ((mpfr::isnan(matrixR.coeff(it.row(), it.col()))) && (matrixI.coeff(it.row(), it.col()) == 0))
                    ;
                else if ((matrixR.coeff(it.row(), it.col()) == mpreal("Inf")) && (isfinite(matrixI.coeff(it.row(), it.col()))))
                    ;
                else
                    result.matrixI.insert(it.row(), it.col()) = ((isComplex) && (matrixI.coeff(it.row(), it.col()) < 0) ? -mpfr::sqrt((it.value() - matrixR.coeff(it.row(), it.col()))/2) : mpfr::sqrt((it.value() - matrixR.coeff(it.row(), it.col()))/2));
            }

        // We may have created some new zero coefficients...
        result.matrixR.prune(0, 0);
        result.matrixI.prune(0, 0);

        result.matrixR.makeCompressed();
        result.matrixI.makeCompressed();
        result.checkComplexity();
    }

    return result;
}


/* Matrix multiplication a *= b */
SparseGmpEigenMatrix& SparseGmpEigenMatrix::operator*=(const SparseGmpEigenMatrix& b)
{
    SparseGmpEigenMatrix a((*this));

    (*this) = a*b;

    return *this;
}

/* Matrix multiplication c = a*b */
// By default, Eigen assumes aliasing for matrix multiplication, so we need to
// tell him whenever the matrix on the lhs is not on the rhs
SparseGmpEigenMatrix SparseGmpEigenMatrix::operator*(const SparseGmpEigenMatrix& b) const
{
    SparseGmpEigenMatrix result;

    // This function could be called to multiply two SparseGmpEigenMatrix objects
    // within this library. It should thus support multiplication by a scalar.
    if ((numel() == 1) || (b.numel() == 1)) {
        result = (*this).times(b);
        return result;
    }

    if (b.isComplex) {
        if (isComplex) {
            result.matrixR = matrixR*b.matrixR - matrixI*b.matrixI;
            result.matrixI = matrixR*b.matrixI + matrixI*b.matrixR;
            result.matrixR.prune(0,0);
            result.matrixI.prune(0,0);
            result.matrixR.makeCompressed();
            result.matrixI.makeCompressed();
            result.checkComplexity();
        } else {
            result.matrixR = (matrixR*b.matrixR).pruned();
            result.matrixI = (matrixR*b.matrixI).pruned();
            result.matrixR.prune(0,0);
            result.matrixI.prune(0,0);
            result.matrixR.makeCompressed();
            result.matrixI.makeCompressed();
            result.checkComplexity();
        }
    } else if (isComplex) {
        result.matrixR = (matrixR*b.matrixR).pruned();
        result.matrixI = (matrixI*b.matrixR).pruned();
        result.matrixR.prune(0,0);
        result.matrixI.prune(0,0);
        result.matrixR.makeCompressed();
        result.matrixI.makeCompressed();
        result.checkComplexity();
    } else {
        result.matrixR = (matrixR*b.matrixR).pruned();
        result.matrixR.prune(0,0);
        result.matrixR.makeCompressed();
        result.isComplex = false;
    }

    return result;
}

/* Matrix multiplication c = a*b */
SparseGmpEigenMatrix& SparseGmpEigenMatrix::mtimes_new(const SparseGmpEigenMatrix& b) const
{
    SparseGmpEigenMatrix& result(*(new SparseGmpEigenMatrix));

    if (b.isComplex) {
        if (isComplex) {
            result.matrixR = matrixR*b.matrixR - matrixI*b.matrixI;
            result.matrixI = matrixR*b.matrixI + matrixI*b.matrixR;
            result.matrixR.prune(0,0);
            result.matrixI.prune(0,0);
            result.matrixR.makeCompressed();
            result.matrixI.makeCompressed();
            result.checkComplexity();
        } else {
            result.matrixR = matrixR*b.matrixR;//).pruned(); // Using the pruned multiplication seems not to always remove all zeros, e.g. in m=[1 2; 3 5]; inv(sgem(m))*sgem(m). Moreover, it now gives segmentation faults... (!)
            result.matrixI = matrixR*b.matrixI;//).pruned();
            result.matrixR.prune(0,0);
            result.matrixI.prune(0,0);
            result.matrixR.makeCompressed();
            result.matrixI.makeCompressed();
            result.checkComplexity();
        }
    } else if (isComplex) {
        result.matrixR = matrixR*b.matrixR;//).pruned();
        result.matrixI = matrixI*b.matrixR;//).pruned();
        result.matrixR.prune(0,0);
        result.matrixI.prune(0,0);
        result.matrixR.makeCompressed();
        result.matrixI.makeCompressed();
        result.checkComplexity();
    } else {
        result.matrixR = matrixR*b.matrixR;//).pruned();
        result.matrixR.prune(0,0);
        result.matrixR.makeCompressed();
        result.isComplex = false;
    }

    return result;
}

/* multiplication with a full matrix c = a*b */
GmpEigenMatrix SparseGmpEigenMatrix::mtimes_sf(const GmpEigenMatrix& b) const
{
    GmpEigenMatrix result;

    if (b.isComplex) {
        if (isComplex) {
            result.matrixR.noalias() = matrixR*b.matrixR - matrixI*b.matrixI;
            result.matrixI.noalias() = matrixR*b.matrixI + matrixI*b.matrixR;
            result.checkComplexity();
        } else {
            result.matrixR.noalias() = matrixR*b.matrixR;
            result.matrixI.noalias() = matrixR*b.matrixI;
            result.checkComplexity();
        }
    } else if (isComplex) {
        result.matrixR.noalias() = matrixR*b.matrixR;
        result.matrixI.noalias() = matrixI*b.matrixR;
        result.checkComplexity();
    } else {
        result.matrixR.noalias() = matrixR*b.matrixR;
        result.isComplex = false;
    }

    return result;
}

/* multiplication with a full matrix c = a*b */
GmpEigenMatrix& SparseGmpEigenMatrix::mtimes_sf_new(const GmpEigenMatrix& b) const
{
    GmpEigenMatrix& result(*(new GmpEigenMatrix));

    if (b.isComplex) {
        if (isComplex) {
            result.matrixR.noalias() = matrixR*b.matrixR - matrixI*b.matrixI;
            result.matrixI.noalias() = matrixR*b.matrixI + matrixI*b.matrixR;
            result.checkComplexity();
        } else {
            result.matrixR.noalias() = matrixR*b.matrixR;
            result.matrixI.noalias() = matrixR*b.matrixI;
            result.checkComplexity();
        }
    } else if (isComplex) {
        result.matrixR.noalias() = matrixR*b.matrixR;
        result.matrixI.noalias() = matrixI*b.matrixR;
        result.checkComplexity();
    } else {
        result.matrixR.noalias() = matrixR*b.matrixR;
        result.isComplex = false;
    }

    return result;
}

/* Tensor product c = kron(a,b) */
SparseGmpEigenMatrix SparseGmpEigenMatrix::kron(const SparseGmpEigenMatrix& b) const
{
    SparseGmpEigenMatrix result;

    if (b.isComplex) {
        if (isComplex) {
            result.matrixR = kroneckerProduct(matrixR, b.matrixR);
            //result.matrixR -= kroneckerProduct(matrixI, b.matrixI); // Doesn't wotk on Eigen 3.3 anymore
            result.matrixR = result.matrixR - kroneckerProduct(matrixI, b.matrixI);
            result.matrixI = kroneckerProduct(matrixR, b.matrixI);
            //result.matrixI += kroneckerProduct(matrixI, b.matrixR);
            result.matrixI = result.matrixI + kroneckerProduct(matrixI, b.matrixR);
            result.matrixR.prune(0,0);
            result.matrixI.prune(0,0);
            result.matrixR.makeCompressed();
            result.matrixI.makeCompressed();
            result.checkComplexity();
        } else {
            result.matrixR = kroneckerProduct(matrixR, b.matrixR);
            result.matrixI = kroneckerProduct(matrixR, b.matrixI);
            result.checkComplexity();
        }
    } else if (isComplex) {
        result.matrixR = kroneckerProduct(matrixR, b.matrixR);
        result.matrixI = kroneckerProduct(matrixI, b.matrixR);
        result.checkComplexity();
    } else {
        result.matrixR = kroneckerProduct(matrixR, b.matrixR);
        result.isComplex = false;
    }

    return result;
}

/* Tensor product c = kron(a,b) */
SparseGmpEigenMatrix& SparseGmpEigenMatrix::kron_new(const SparseGmpEigenMatrix& b) const
{
    SparseGmpEigenMatrix& result(*(new SparseGmpEigenMatrix));

    if (b.isComplex) {
        if (isComplex) {
            result.matrixR = kroneckerProduct(matrixR, b.matrixR);
            //result.matrixR -= kroneckerProduct(matrixI, b.matrixI); // Doesn't wotk on Eigen 3.3 anymore
            result.matrixR = result.matrixR - kroneckerProduct(matrixI, b.matrixI);
            result.matrixI = kroneckerProduct(matrixR, b.matrixI);
            //result.matrixI += kroneckerProduct(matrixI, b.matrixR);
            result.matrixI = result.matrixI + kroneckerProduct(matrixI, b.matrixR);
            result.matrixR.prune(0,0);
            result.matrixI.prune(0,0);
            result.matrixR.makeCompressed();
            result.matrixI.makeCompressed();
            result.checkComplexity();
        } else {
            result.matrixR = kroneckerProduct(matrixR, b.matrixR);
            result.matrixI = kroneckerProduct(matrixR, b.matrixI);
            result.checkComplexity();
        }
    } else if (isComplex) {
        result.matrixR = kroneckerProduct(matrixR, b.matrixR);
        result.matrixI = kroneckerProduct(matrixI, b.matrixR);
        result.checkComplexity();
    } else {
        result.matrixR = kroneckerProduct(matrixR, b.matrixR);
        result.isComplex = false;
    }

    return result;
}

/* Tensor product with a full matrix c = kron(a,b) */
SparseGmpEigenMatrix SparseGmpEigenMatrix::kron_sf(const GmpEigenMatrix& b) const
{
    SparseGmpEigenMatrix result;

    if (b.isComplex) {
        if (isComplex) {
            result.matrixR = kroneckerProduct(matrixR, b.matrixR);
            //result.matrixR -= kroneckerProduct(matrixI, b.matrixI); // Doesn't wotk on Eigen 3.3 anymore
            result.matrixR = result.matrixR - kroneckerProduct(matrixI, b.matrixI);
            result.matrixI = kroneckerProduct(matrixR, b.matrixI);
            //result.matrixI += kroneckerProduct(matrixI, b.matrixR);
            result.matrixI = result.matrixI + kroneckerProduct(matrixI, b.matrixR);
            result.matrixR.prune(0,0);
            result.matrixI.prune(0,0);
            result.matrixR.makeCompressed();
            result.matrixI.makeCompressed();
            result.checkComplexity();
        } else {
            result.matrixR = kroneckerProduct(matrixR, b.matrixR);
            result.matrixI = kroneckerProduct(matrixR, b.matrixI);
            result.matrixR.prune(0,0);
            result.matrixI.prune(0,0);
            result.matrixR.makeCompressed();
            result.matrixI.makeCompressed();
            result.checkComplexity();
        }
    } else if (isComplex) {
        result.matrixR = kroneckerProduct(matrixR, b.matrixR);
        result.matrixI = kroneckerProduct(matrixI, b.matrixR);
        result.matrixR.prune(0,0);
        result.matrixI.prune(0,0);
        result.matrixR.makeCompressed();
        result.matrixI.makeCompressed();
        result.checkComplexity();
    } else {
        result.matrixR = kroneckerProduct(matrixR, b.matrixR);
        result.matrixR.prune(0,0);
        result.matrixR.makeCompressed();
        result.isComplex = false;
    }

    return result;
}

/* Tensor product with a full matrix c = kron(a,b) */
SparseGmpEigenMatrix& SparseGmpEigenMatrix::kron_sf_new(const GmpEigenMatrix& b) const
{
    SparseGmpEigenMatrix& result(*(new SparseGmpEigenMatrix));

    if (b.isComplex) {
        if (isComplex) {
            result.matrixR = kroneckerProduct(matrixR, b.matrixR);
            //result.matrixR -= kroneckerProduct(matrixI, b.matrixI); // Doesn't wotk on Eigen 3.3 anymore
            result.matrixR = result.matrixR - kroneckerProduct(matrixI, b.matrixI);
            result.matrixI = kroneckerProduct(matrixR, b.matrixI);
            //result.matrixI += kroneckerProduct(matrixI, b.matrixR);
            result.matrixI = result.matrixI + kroneckerProduct(matrixI, b.matrixR);
            result.matrixR.prune(0,0);
            result.matrixI.prune(0,0);
            result.matrixR.makeCompressed();
            result.matrixI.makeCompressed();
            result.checkComplexity();
        } else {
            result.matrixR = kroneckerProduct(matrixR, b.matrixR);
            result.matrixI = kroneckerProduct(matrixR, b.matrixI);
            result.matrixR.prune(0,0);
            result.matrixI.prune(0,0);
            result.matrixR.makeCompressed();
            result.matrixI.makeCompressed();
            result.checkComplexity();
        }
    } else if (isComplex) {
        result.matrixR = kroneckerProduct(matrixR, b.matrixR);
        result.matrixI = kroneckerProduct(matrixI, b.matrixR);
        result.matrixR.prune(0,0);
        result.matrixI.prune(0,0);
        result.matrixR.makeCompressed();
        result.matrixI.makeCompressed();
        result.checkComplexity();
    } else {
        result.matrixR = kroneckerProduct(matrixR, b.matrixR);
        result.matrixR.prune(0,0);
        result.matrixR.makeCompressed();
        result.isComplex = false;
    }

    return result;
}



/* ----------------------------------------------------------
   | Some mathematical functions defined for real arguments |
   |       (the real test is made on the matlab side.       |
   |       Here we simply ignore the imaginary parts)       |
   ---------------------------------------------------------- */


/* cubic root b = cqrt(a) -- only works if the real part is positive and the imaginary part zero (!) */
SparseGmpEigenMatrix SparseGmpEigenMatrix::cbrt() const
{
   SparseGmpEigenMatrix result;

   result.isComplex = false;
   result.matrixR.resize(matrixR.rows(), matrixR.cols());
   result.matrixR.reserve(matrixR.nonZeros());
   for (IndexType k = 0; k < matrixR.outerSize(); ++k)
       for (SparseMatrix<mpreal>::InnerIterator it(matrixR,k); it; ++it)
           result.matrixR.insert(it.row(), it.col()) = mpfr::cbrt(it.value());

   return result;
}

/* cubic root b = cqrt(a) -- only works if the real part is positive and the imaginary part zero (!) */
SparseGmpEigenMatrix& SparseGmpEigenMatrix::cbrt_new() const
{
   SparseGmpEigenMatrix& result(*(new SparseGmpEigenMatrix));

   result.isComplex = false;
   result.matrixR.resize(matrixR.rows(), matrixR.cols());
   result.matrixR.reserve(matrixR.nonZeros());
   for (IndexType k = 0; k < matrixR.outerSize(); ++k)
       for (SparseMatrix<mpreal>::InnerIterator it(matrixR,k); it; ++it)
           result.matrixR.insert(it.row(), it.col()) = mpfr::cbrt(it.value());

   return result;
}

/* sine function b = sin(a) */
SparseGmpEigenMatrix SparseGmpEigenMatrix::sin() const
{
   SparseGmpEigenMatrix result;

   result.isComplex = isComplex;
   if (isComplex) {
       result.matrixR.resize(matrixR.rows(),matrixR.cols());
       result.matrixI.resize(matrixR.rows(),matrixR.cols());

       // We will first copy the data to triplets
       vector< Triplet<mpreal> > tripletListR, tripletListI;
       tripletListR.reserve(matrixR.nonZeros());
       tripletListI.reserve(matrixI.nonZeros());

       for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
           SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
           SparseMatrix<mpreal>::InnerIterator itI(matrixI,k);

           while ((itR) || (itI)) {
               if ((itR) && (itI)) {
                   if (itR.row() < itI.row()) {
                       tripletListR.push_back(Triplet<mpreal>(itR.row(), k, mpfr::sin(itR.value())));
                       ++itR;
                   } else if (itR.row() == itI.row()) {
                       tripletListR.push_back(Triplet<mpreal>(itR.row(), k, mpfr::sin(itR.value())*mpfr::cosh(itI.value())));
                       tripletListI.push_back(Triplet<mpreal>(itR.row(), k, mpfr::cos(itR.value())*mpfr::sinh(itI.value())));
                       ++itR;
                       ++itI;
                   } else {
                       tripletListI.push_back(Triplet<mpreal>(itR.row(), k, mpfr::sinh(itR.value())));
                       ++itI;
                   }
               } else if (itR) {
                   tripletListR.push_back(Triplet<mpreal>(itR.row(), k, mpfr::sin(itR.value())));
                   ++itR;
               } else {
                   tripletListI.push_back(Triplet<mpreal>(itR.row(), k, mpfr::sinh(itR.value())));
                   ++itI;
               }
           }
       }
       result.matrixR.setFromTriplets(tripletListR.begin(), tripletListR.end());
       result.matrixI.setFromTriplets(tripletListI.begin(), tripletListI.end());

       // We compress the output if possible
       result.matrixR.prune(0,0);
       result.matrixI.prune(0,0);
       result.matrixR.makeCompressed();
       result.matrixI.makeCompressed();
       result.checkComplexity();
   } else {
       result.matrixR.resize(matrixR.rows(),matrixR.cols());
       result.matrixR.reserve(matrixR.nonZeros());
       for (IndexType k = 0; k < matrixR.outerSize(); ++k)
           for (SparseMatrix<mpreal>::InnerIterator it(matrixR,k); it; ++it)
               result.matrixR.insert(it.row(), it.col()) = mpfr::sin(it.value());
       result.matrixR.prune(0,0);
       result.matrixR.makeCompressed();
   }

   return result;
}

/* sine function b = sin(a) */
SparseGmpEigenMatrix& SparseGmpEigenMatrix::sin_new() const
{
   SparseGmpEigenMatrix& result(*(new SparseGmpEigenMatrix));

   result.isComplex = isComplex;
   if (isComplex) {
       result.matrixR.resize(matrixR.rows(),matrixR.cols());
       result.matrixI.resize(matrixR.rows(),matrixR.cols());

       // We will first copy the data to triplets
       vector< Triplet<mpreal> > tripletListR, tripletListI;
       tripletListR.reserve(matrixR.nonZeros());
       tripletListI.reserve(matrixI.nonZeros());

       for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
           SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
           SparseMatrix<mpreal>::InnerIterator itI(matrixI,k);

           while ((itR) || (itI)) {
               if ((itR) && (itI)) {
                   if (itR.row() < itI.row()) {
                       tripletListR.push_back(Triplet<mpreal>(itR.row(), k, mpfr::sin(itR.value())));
                       ++itR;
                   } else if (itR.row() == itI.row()) {
                       tripletListR.push_back(Triplet<mpreal>(itR.row(), k, mpfr::sin(itR.value())*mpfr::cosh(itI.value())));
                       tripletListI.push_back(Triplet<mpreal>(itR.row(), k, mpfr::cos(itR.value())*mpfr::sinh(itI.value())));
                       ++itR;
                       ++itI;
                   } else {
                       tripletListI.push_back(Triplet<mpreal>(itR.row(), k, mpfr::sinh(itR.value())));
                       ++itI;
                   }
               } else if (itR) {
                   tripletListR.push_back(Triplet<mpreal>(itR.row(), k, mpfr::sin(itR.value())));
                   ++itR;
               } else {
                   tripletListI.push_back(Triplet<mpreal>(itR.row(), k, mpfr::sinh(itR.value())));
                   ++itI;
               }
           }
       }
       result.matrixR.setFromTriplets(tripletListR.begin(), tripletListR.end());
       result.matrixI.setFromTriplets(tripletListI.begin(), tripletListI.end());

       // We compress the output if possible
       result.matrixR.prune(0,0);
       result.matrixI.prune(0,0);
       result.matrixR.makeCompressed();
       result.matrixI.makeCompressed();
       result.checkComplexity();
   } else {
       result.matrixR.resize(matrixR.rows(),matrixR.cols());
       result.matrixR.reserve(matrixR.nonZeros());
       for (IndexType k = 0; k < matrixR.outerSize(); ++k)
           for (SparseMatrix<mpreal>::InnerIterator it(matrixR,k); it; ++it)
               result.matrixR.insert(it.row(), it.col()) = mpfr::sin(it.value());
       result.matrixR.prune(0,0);
       result.matrixR.makeCompressed();
   }

   return result;
}

/* cosine function b = cos(a) on non-zero elements
   WARNING : this is an internal function, which only computes the cosine of
   the non-zero matrix elements */
SparseGmpEigenMatrix SparseGmpEigenMatrix::cos_nz() const
{
   SparseGmpEigenMatrix result;

   result.isComplex = isComplex;
   if (isComplex) {
       result.matrixR.resize(matrixR.rows(),matrixR.cols());
       result.matrixI.resize(matrixR.rows(),matrixR.cols());

       // We will first copy the data to triplets
       vector< Triplet<mpreal> > tripletListR, tripletListI;
       tripletListR.reserve(matrixR.nonZeros());
       tripletListI.reserve(matrixI.nonZeros());

       for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
           SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
           SparseMatrix<mpreal>::InnerIterator itI(matrixI,k);

           while ((itR) || (itI)) {
               if ((itR) && (itI)) {
                   if (itR.row() < itI.row()) {
                       tripletListR.push_back(Triplet<mpreal>(itR.row(), k, mpfr::cos(itR.value())));
                       ++itR;
                   } else if (itR.row() == itI.row()) {
                       tripletListR.push_back(Triplet<mpreal>(itR.row(), k, mpfr::cos(itR.value())*mpfr::cosh(itI.value())));
                       tripletListI.push_back(Triplet<mpreal>(itR.row(), k, -mpfr::sin(itR.value())*mpfr::sinh(itI.value())));
                       ++itR;
                       ++itI;
                   } else {
                       tripletListI.push_back(Triplet<mpreal>(itR.row(), k, mpfr::cosh(itR.value())));
                       ++itI;
                   }
               } else if (itR) {
                   tripletListR.push_back(Triplet<mpreal>(itR.row(), k, mpfr::cos(itR.value())));
                   ++itR;
               } else {
                   tripletListI.push_back(Triplet<mpreal>(itR.row(), k, mpfr::cosh(itR.value())));
                   ++itI;
               }
           }
       }
       result.matrixR.setFromTriplets(tripletListR.begin(), tripletListR.end());
       result.matrixI.setFromTriplets(tripletListI.begin(), tripletListI.end());

       // We compress the output if possible
       result.matrixR.prune(0,0);
       result.matrixI.prune(0,0);
       result.matrixR.makeCompressed();
       result.matrixI.makeCompressed();
       result.checkComplexity();
   } else {
       result.matrixR.resize(matrixR.rows(),matrixR.cols());
       result.matrixR.reserve(matrixR.nonZeros());
       for (IndexType k = 0; k < matrixR.outerSize(); ++k)
           for (SparseMatrix<mpreal>::InnerIterator it(matrixR,k); it; ++it)
               result.matrixR.insert(it.row(), it.col()) = mpfr::cos(it.value());
       result.matrixR.prune(0,0);
       result.matrixR.makeCompressed();
       result.checkComplexity();
   }

   return result;
}

/* tangent function b = tan(a) */
SparseGmpEigenMatrix SparseGmpEigenMatrix::tan() const
{
   SparseGmpEigenMatrix result;

   result.isComplex = isComplex;
   if (isComplex) {
       result.matrixR.resize(matrixR.rows(),matrixR.cols());
       result.matrixI.resize(matrixR.rows(),matrixR.cols());

       // We will first copy the data to triplets
       vector< Triplet<mpreal> > tripletListR, tripletListI;
       tripletListR.reserve(matrixR.nonZeros());
       tripletListI.reserve(matrixI.nonZeros());

       for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
           SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
           SparseMatrix<mpreal>::InnerIterator itI(matrixI,k);

           while ((itR) || (itI)) {
               if ((itR) && (itI)) {
                   if (itR.row() < itI.row()) {
                       tripletListR.push_back(Triplet<mpreal>(itR.row(), k, mpfr::tan(itR.value())));
                       ++itR;
                   } else if (itR.row() == itI.row()) {
                       tripletListR.push_back(Triplet<mpreal>(itR.row(), k, mpfr::sin(2*itR.value())/(mpfr::cos(2*itR.value()) + mpfr::cosh(2*itI.value()))));
                       tripletListI.push_back(Triplet<mpreal>(itR.row(), k, mpfr::sinh(2*itI.value())/(mpfr::cos(2*itR.value()) + mpfr::cosh(2*itI.value()))));
                       ++itR;
                       ++itI;
                   } else {
                       tripletListI.push_back(Triplet<mpreal>(itR.row(), k, mpfr::tanh(itR.value())));
                       ++itI;
                   }
               } else if (itR) {
                   tripletListR.push_back(Triplet<mpreal>(itR.row(), k, mpfr::tan(itR.value())));
                   ++itR;
               } else {
                   tripletListI.push_back(Triplet<mpreal>(itR.row(), k, mpfr::tanh(itR.value())));
                   ++itI;
               }
           }
       }

       result.matrixR.setFromTriplets(tripletListR.begin(), tripletListR.end());
       result.matrixI.setFromTriplets(tripletListI.begin(), tripletListI.end());

       // We compress the output if possible
       result.matrixR.prune(0,0);
       result.matrixI.prune(0,0);
       result.matrixR.makeCompressed();
       result.matrixI.makeCompressed();
       result.checkComplexity();
   } else {
       result.matrixR.resize(matrixR.rows(),matrixR.cols());
       result.matrixR.reserve(matrixR.nonZeros());
       for (IndexType k = 0; k < matrixR.outerSize(); ++k)
           for (SparseMatrix<mpreal>::InnerIterator it(matrixR,k); it; ++it)
               result.matrixR.insert(it.row(), it.col()) = mpfr::tan(it.value());
       result.matrixR.prune(0,0);
       result.matrixR.makeCompressed();
   }

   return result;
}

/* tangent function b = tan(a) */
SparseGmpEigenMatrix& SparseGmpEigenMatrix::tan_new() const
{
   SparseGmpEigenMatrix& result(*(new SparseGmpEigenMatrix));

   result.isComplex = isComplex;
   if (isComplex) {
       result.matrixR.resize(matrixR.rows(),matrixR.cols());
       result.matrixI.resize(matrixR.rows(),matrixR.cols());

       // We will first copy the data to triplets
       vector< Triplet<mpreal> > tripletListR, tripletListI;
       tripletListR.reserve(matrixR.nonZeros());
       tripletListI.reserve(matrixI.nonZeros());

       for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
           SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
           SparseMatrix<mpreal>::InnerIterator itI(matrixI,k);

           while ((itR) || (itI)) {
               if ((itR) && (itI)) {
                   if (itR.row() < itI.row()) {
                       tripletListR.push_back(Triplet<mpreal>(itR.row(), k, mpfr::tan(itR.value())));
                       ++itR;
                   } else if (itR.row() == itI.row()) {
                       tripletListR.push_back(Triplet<mpreal>(itR.row(), k, mpfr::sin(2*itR.value())/(mpfr::cos(2*itR.value()) + mpfr::cosh(2*itI.value()))));
                       tripletListI.push_back(Triplet<mpreal>(itR.row(), k, mpfr::sinh(2*itI.value())/(mpfr::cos(2*itR.value()) + mpfr::cosh(2*itI.value()))));
                       ++itR;
                       ++itI;
                   } else {
                       tripletListI.push_back(Triplet<mpreal>(itR.row(), k, mpfr::tanh(itR.value())));
                       ++itI;
                   }
               } else if (itR) {
                   tripletListR.push_back(Triplet<mpreal>(itR.row(), k, mpfr::tan(itR.value())));
                   ++itR;
               } else {
                   tripletListI.push_back(Triplet<mpreal>(itR.row(), k, mpfr::tanh(itR.value())));
                   ++itI;
               }
           }
       }

       result.matrixR.setFromTriplets(tripletListR.begin(), tripletListR.end());
       result.matrixI.setFromTriplets(tripletListI.begin(), tripletListI.end());

       // We compress the output if possible
       result.matrixR.prune(0,0);
       result.matrixI.prune(0,0);
       result.matrixR.makeCompressed();
       result.matrixI.makeCompressed();
       result.checkComplexity();
   } else {
       result.matrixR.resize(matrixR.rows(),matrixR.cols());
       result.matrixR.reserve(matrixR.nonZeros());
       for (IndexType k = 0; k < matrixR.outerSize(); ++k)
           for (SparseMatrix<mpreal>::InnerIterator it(matrixR,k); it; ++it)
               result.matrixR.insert(it.row(), it.col()) = mpfr::tan(it.value());
       result.matrixR.prune(0,0);
       result.matrixR.makeCompressed();
   }

   return result;
}

/* arc sine function b = asin(a) */
SparseGmpEigenMatrix SparseGmpEigenMatrix::asin() const
{
    SparseGmpEigenMatrix result;

    result.isComplex = true;
    result.matrixR.resize(matrixR.rows(),matrixR.cols());
    result.matrixI.resize(matrixR.rows(),matrixR.cols());

    // We will first copy the data to triplets
    vector< Triplet<mpreal> > tripletListR, tripletListI;
    tripletListR.reserve(matrixR.nonZeros());
    tripletListI.reserve(matrixI.nonZeros());

    if (isComplex) {
        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
            SparseMatrix<mpreal>::InnerIterator itI(matrixI,k);

            IndexType row;
            GmpEigenMatrix value;
            while ((itR) || (itI)) {
                if ((itR) && (itI)) {
                    if (itR.row() < itI.row()) {
                        value = GmpEigenMatrix(itR.value());
                        row = itR.row();
                        ++itR;
                    } else if (itR.row() == itI.row()) {
                        value = GmpEigenMatrix(itR.value(), itI.value());
                        row = itR.row();
                        ++itR;
                        ++itI;
                    } else {
                        value = GmpEigenMatrix(0, itI.value());
                        row = itI.row();
                        ++itI;
                    }
               } else if (itR) {
                    value = GmpEigenMatrix(itR.value());
                    row = itR.row();
                    ++itR;
               } else {
                    value = GmpEigenMatrix(0, itI.value());
                    row = itI.row();
                    ++itI;
               }

               if ((value.matrixI(0,0) == 0) && (value.matrixR(0,0) >= -1) && (value.matrixR(0,0) <= 1)) {
                    tripletListR.push_back(Triplet<mpreal>(row, k, mpfr::asin(value.matrixR(0,0))));
               } else {
                    GmpEigenMatrix image;
                    image = -constI()*(constI()*value + (GmpEigenMatrix(1) - value.power(GmpEigenMatrix(2))).sqrt()).log();
                    if (image.matrixR(0,0) != 0)
                        tripletListR.push_back(Triplet<mpreal>(row, k, image.matrixR(0,0)));
                    if (image.matrixI(0,0) != 0)
                        tripletListI.push_back(Triplet<mpreal>(row, k, image.matrixI(0,0)));
               }
            }
        }
    } else {
        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);

            IndexType row;
            GmpEigenMatrix value;
            while (itR) {
                if ((itR.value() >= -1) && (itR.value() <= 1)) {
                    tripletListR.push_back(Triplet<mpreal>(itR.row(), k, mpfr::asin(itR.value())));
                } else {
                    GmpEigenMatrix image;
                    image = -constI()*(constI()*itR.value() + (GmpEigenMatrix(1) - pow(itR.value(),2)).sqrt()).log();
                    if (image.matrixR(0,0) != 0)
                        tripletListR.push_back(Triplet<mpreal>(itR.row(), k, image.matrixR(0,0)));
                    if (image.matrixI(0,0) != 0)
                        tripletListI.push_back(Triplet<mpreal>(itR.row(), k, image.matrixI(0,0)));
                }
                ++itR;
            }
        }
    }

    result.matrixR.setFromTriplets(tripletListR.begin(), tripletListR.end());
    result.matrixI.setFromTriplets(tripletListI.begin(), tripletListI.end());

    // We compress the output if possible
    result.matrixR.prune(0,0);
    result.matrixI.prune(0,0);
    result.checkComplexity();

    return result;
}

/* arc sine function b = asin(a) */
SparseGmpEigenMatrix& SparseGmpEigenMatrix::asin_new() const
{
    SparseGmpEigenMatrix& result(*(new SparseGmpEigenMatrix));

    result.isComplex = true;
    result.matrixR.resize(matrixR.rows(),matrixR.cols());
    result.matrixI.resize(matrixR.rows(),matrixR.cols());

    // We will first copy the data to triplets
    vector< Triplet<mpreal> > tripletListR, tripletListI;
    tripletListR.reserve(matrixR.nonZeros());
    tripletListI.reserve(matrixI.nonZeros());

    if (isComplex) {
        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
            SparseMatrix<mpreal>::InnerIterator itI(matrixI,k);

            IndexType row;
            GmpEigenMatrix value;
            while ((itR) || (itI)) {
                if ((itR) && (itI)) {
                    if (itR.row() < itI.row()) {
                        value = GmpEigenMatrix(itR.value());
                        row = itR.row();
                        ++itR;
                    } else if (itR.row() == itI.row()) {
                        value = GmpEigenMatrix(itR.value(), itI.value());
                        row = itR.row();
                        ++itR;
                        ++itI;
                    } else {
                        value = GmpEigenMatrix(0, itI.value());
                        row = itI.row();
                        ++itI;
                    }
               } else if (itR) {
                    value = GmpEigenMatrix(itR.value());
                    row = itR.row();
                    ++itR;
               } else {
                    value = GmpEigenMatrix(0, itI.value());
                    row = itI.row();
                    ++itI;
               }

               if ((value.matrixI(0,0) == 0) && (value.matrixR(0,0) >= -1) && (value.matrixR(0,0) <= 1)) {
                    tripletListR.push_back(Triplet<mpreal>(row, k, mpfr::asin(value.matrixR(0,0))));
               } else {
                    GmpEigenMatrix image;
                    image = -constI()*(constI()*value + (GmpEigenMatrix(1) - value.power(GmpEigenMatrix(2))).sqrt()).log();
                    if (image.matrixR(0,0) != 0)
                        tripletListR.push_back(Triplet<mpreal>(row, k, image.matrixR(0,0)));
                    if (image.matrixI(0,0) != 0)
                        tripletListI.push_back(Triplet<mpreal>(row, k, image.matrixI(0,0)));
               }
            }
        }
    } else {
        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);

            IndexType row;
            GmpEigenMatrix value;
            while (itR) {
                if ((itR.value() >= -1) && (itR.value() <= 1)) {
                    tripletListR.push_back(Triplet<mpreal>(itR.row(), k, mpfr::asin(itR.value())));
                } else {
                    GmpEigenMatrix image;
                    image = -constI()*(constI()*itR.value() + (GmpEigenMatrix(1) - pow(itR.value(),2)).sqrt()).log();
                    if (image.matrixR(0,0) != 0)
                        tripletListR.push_back(Triplet<mpreal>(itR.row(), k, image.matrixR(0,0)));
                    if (image.matrixI(0,0) != 0)
                        tripletListI.push_back(Triplet<mpreal>(itR.row(), k, image.matrixI(0,0)));
                }
                ++itR;
            }
        }
    }

    result.matrixR.setFromTriplets(tripletListR.begin(), tripletListR.end());
    result.matrixI.setFromTriplets(tripletListI.begin(), tripletListI.end());

    // We compress the output if possible
    result.matrixR.prune(0,0);
    result.matrixI.prune(0,0);
    result.checkComplexity();

    return result;
}

/* arc tangent function b = atan(a) */
SparseGmpEigenMatrix SparseGmpEigenMatrix::atan() const
{
    SparseGmpEigenMatrix result;

    if (isComplex) {
        result.isComplex = true;
        result.matrixR.resize(matrixR.rows(),matrixR.cols());
        result.matrixI.resize(matrixR.rows(),matrixR.cols());

        // We will first copy the data to triplets
        vector< Triplet<mpreal> > tripletListR, tripletListI;
        tripletListR.reserve(matrixR.nonZeros());
        tripletListI.reserve(matrixI.nonZeros());

        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
            SparseMatrix<mpreal>::InnerIterator itI(matrixI,k);

            IndexType row;
            GmpEigenMatrix value;
            while ((itR) || (itI)) {
                if ((itR) && (itI)) {
                    if (itR.row() < itI.row()) {
                        value = GmpEigenMatrix(itR.value());
                        row = itR.row();
                        ++itR;
                    } else if (itR.row() == itI.row()) {
                        value = GmpEigenMatrix(itR.value(), itI.value());
                        row = itR.row();
                        ++itR;
                        ++itI;
                    } else {
                        value = GmpEigenMatrix(0, itI.value());
                        row = itI.row();
                        ++itI;
                    }
               } else if (itR) {
                    value = GmpEigenMatrix(itR.value());
                    row = itR.row();
                    ++itR;
               } else {
                    value = GmpEigenMatrix(0, itI.value());
                    row = itI.row();
                    ++itI;
               }

               if (value.matrixI(0,0) == 0) {
                    tripletListR.push_back(Triplet<mpreal>(row, k, mpfr::atan(value.matrixR(0,0))));
               } else {
                    GmpEigenMatrix image;
                    image = GmpEigenMatrix(0,mpreal("0.5"))*((GmpEigenMatrix(1)-constI()*value).log() - (GmpEigenMatrix(1)+constI()*value).log());
                    if (image.matrixR(0,0) != 0)
                        tripletListR.push_back(Triplet<mpreal>(row, k, image.matrixR(0,0)));
                    if (image.matrixI(0,0) != 0)
                        tripletListI.push_back(Triplet<mpreal>(row, k, image.matrixI(0,0)));
               }
            }
        }

        result.matrixR.setFromTriplets(tripletListR.begin(), tripletListR.end());
        result.matrixI.setFromTriplets(tripletListI.begin(), tripletListI.end());

        // We compress the output if possible
        result.matrixR.prune(0,0);
        result.matrixI.prune(0,0);
        result.checkComplexity();
    } else {
        result.isComplex = false;
        result.matrixR.resize(matrixR.rows(),matrixR.cols());

        // We will first copy the data to triplets
        vector< Triplet<mpreal> > tripletListR;
        tripletListR.reserve(matrixR.nonZeros());

        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);

            IndexType row;
            GmpEigenMatrix value;
            while (itR) {
                tripletListR.push_back(Triplet<mpreal>(itR.row(), k, mpfr::atan(itR.value())));
                ++itR;
            }
        }
        result.matrixR.setFromTriplets(tripletListR.begin(), tripletListR.end());

        // We compress the output if possible
        result.matrixR.prune(0,0);
    }

    return result;
}

/* arc tangent function b = atan(a) */
SparseGmpEigenMatrix& SparseGmpEigenMatrix::atan_new() const
{
    SparseGmpEigenMatrix& result(*(new SparseGmpEigenMatrix));

    if (isComplex) {
        result.isComplex = true;
        result.matrixR.resize(matrixR.rows(),matrixR.cols());
        result.matrixI.resize(matrixR.rows(),matrixR.cols());

        // We will first copy the data to triplets
        vector< Triplet<mpreal> > tripletListR, tripletListI;
        tripletListR.reserve(matrixR.nonZeros());
        tripletListI.reserve(matrixI.nonZeros());

        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
            SparseMatrix<mpreal>::InnerIterator itI(matrixI,k);

            IndexType row;
            GmpEigenMatrix value;
            while ((itR) || (itI)) {
                if ((itR) && (itI)) {
                    if (itR.row() < itI.row()) {
                        value = GmpEigenMatrix(itR.value());
                        row = itR.row();
                        ++itR;
                    } else if (itR.row() == itI.row()) {
                        value = GmpEigenMatrix(itR.value(), itI.value());
                        row = itR.row();
                        ++itR;
                        ++itI;
                    } else {
                        value = GmpEigenMatrix(0, itI.value());
                        row = itI.row();
                        ++itI;
                    }
               } else if (itR) {
                    value = GmpEigenMatrix(itR.value());
                    row = itR.row();
                    ++itR;
               } else {
                    value = GmpEigenMatrix(0, itI.value());
                    row = itI.row();
                    ++itI;
               }

               if (value.matrixI(0,0) == 0) {
                    tripletListR.push_back(Triplet<mpreal>(row, k, mpfr::atan(value.matrixR(0,0))));
               } else {
                    GmpEigenMatrix image;
                    image = GmpEigenMatrix(0,mpreal("0.5"))*((GmpEigenMatrix(1)-constI()*value).log() - (GmpEigenMatrix(1)+constI()*value).log());
                    if (image.matrixR(0,0) != 0)
                        tripletListR.push_back(Triplet<mpreal>(row, k, image.matrixR(0,0)));
                    if (image.matrixI(0,0) != 0)
                        tripletListI.push_back(Triplet<mpreal>(row, k, image.matrixI(0,0)));
               }
            }
        }

        result.matrixR.setFromTriplets(tripletListR.begin(), tripletListR.end());
        result.matrixI.setFromTriplets(tripletListI.begin(), tripletListI.end());

        // We compress the output if possible
        result.matrixR.prune(0,0);
        result.matrixI.prune(0,0);
        result.checkComplexity();
    } else {
        result.isComplex = false;
        result.matrixR.resize(matrixR.rows(),matrixR.cols());

        // We will first copy the data to triplets
        vector< Triplet<mpreal> > tripletListR;
        tripletListR.reserve(matrixR.nonZeros());

        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);

            IndexType row;
            GmpEigenMatrix value;
            while (itR) {
                tripletListR.push_back(Triplet<mpreal>(itR.row(), k, mpfr::atan(itR.value())));
                ++itR;
            }
        }
        result.matrixR.setFromTriplets(tripletListR.begin(), tripletListR.end());

        // We compress the output if possible
        result.matrixR.prune(0,0);
    }

    return result;
}








/* ---------------------------------------
   |   Some linear algebra operations    |
   --------------------------------------- */

/* Matrix rank b = rank(a) */
IndexType SparseGmpEigenMatrix::rank() const
{
    IndexType result;

    if (isComplex) {
        result = complexIsometry().rank()/2;
    } else {
        // We use a QR solver to compute the matrix rank

        //matrixR.makeCompressed(); // The matrix is supposed to be alread compressed at this stage...
        SparseQR< SparseMatrix<mpreal, ColMajor>, COLAMDOrdering<int> > solver(matrixR);

        result = solver.rank();
    }

    return result;
}

/* Linear system solving c = a\b, with dense b*/
GmpEigenMatrix SparseGmpEigenMatrix::mldivide_sf(const GmpEigenMatrix& b) const
{
    GmpEigenMatrix result;

    if ((isComplex) || (b.isComplex)) {
        result = complexIsometry().mldivide_sf(b.complexIsometry()).complexIsometryInverse();
    } else {
        SparseQR< SparseMatrix<mpreal, ColMajor>, COLAMDOrdering<int> > solver(matrixR);

        result.isComplex = false;
        result.matrixR = solver.solve(b.matrixR);
        result.matrixI.resize(0,0);
    }

    return result;
}

/* Linear system solving c = a\b, with dense b*/
GmpEigenMatrix& SparseGmpEigenMatrix::mldivide_sf_new(const GmpEigenMatrix& b) const
{
    GmpEigenMatrix& result(*(new GmpEigenMatrix));

    if ((isComplex) || (b.isComplex)) {
        result = complexIsometry().mldivide_sf(b.complexIsometry()).complexIsometryInverse();
    } else {
        SparseQR< SparseMatrix<mpreal, ColMajor>, COLAMDOrdering<int> > solver(matrixR);

        result.isComplex = false;
        result.matrixR = solver.solve(b.matrixR);
        result.matrixI.resize(0,0);
    }

    return result;
}

/* Linear system solving c = a\b, with sparse b*/
SparseGmpEigenMatrix SparseGmpEigenMatrix::mldivide(const SparseGmpEigenMatrix& b) const
{
    SparseGmpEigenMatrix result;

    if ((isComplex) || (b.isComplex)) {
        result = complexIsometry().mldivide(b.complexIsometry()).complexIsometryInverse();
    } else {
        SparseQR< SparseMatrix<mpreal, ColMajor>, COLAMDOrdering<int> > solver;

        solver.analyzePattern(matrixR);
        solver.factorize(matrixR);

        result.isComplex = false;
        result.matrixI.resize(0,0);

        // We solve for each column iteratively
        result.matrixR.resize(b.matrixR.rows(), b.matrixR.cols());
        result.matrixR.reserve(b.matrixR.nonZeros()*2);
        for (IndexType i(0); i < b.matrixR.cols(); ++i) {
            Matrix<mpreal, Dynamic, 1> colb(b.matrixR.col(i));
            SparseMatrix<mpreal> solb(solver.solve(colb).sparseView(mpreal(0),mpreal(1)));

            SparseMatrix<mpreal>::InnerIterator itR(solb,0);
            while (itR) {
                result.matrixR.insert(itR.row(), i) = itR.value();
                ++itR;
            }
        }
        result.matrixR.makeCompressed();
    }

    return result;
}

/* Linear system solving c = a\b, with sparse b*/
SparseGmpEigenMatrix& SparseGmpEigenMatrix::mldivide_new(const SparseGmpEigenMatrix& b) const
{
    SparseGmpEigenMatrix& result(*(new SparseGmpEigenMatrix));

    if ((isComplex) || (b.isComplex)) {
        result = complexIsometry().mldivide(b.complexIsometry()).complexIsometryInverse();
    } else {
        SparseQR< SparseMatrix<mpreal, ColMajor>, COLAMDOrdering<int> > solver;

        solver.analyzePattern(matrixR);
        solver.factorize(matrixR);

        result.isComplex = false;
        result.matrixI.resize(0,0);

        // We solve for each column iteratively
        result.matrixR.resize(b.matrixR.rows(), b.matrixR.cols());
        result.matrixR.reserve(b.matrixR.nonZeros()*2);
        for (IndexType i(0); i < b.matrixR.cols(); ++i) {
            Matrix<mpreal, Dynamic, 1> colb(b.matrixR.col(i));
            SparseMatrix<mpreal> solb(solver.solve(colb).sparseView(mpreal(0),mpreal(1)));

            SparseMatrix<mpreal>::InnerIterator itR(solb,0);
            while (itR) {
                result.matrixR.insert(itR.row(), i) = itR.value();
                ++itR;
            }
        }
        result.matrixR.makeCompressed();
    }

    return result;
}

/* Matrix inverse b = inv(a) */
SparseGmpEigenMatrix SparseGmpEigenMatrix::inv() const
{
    SparseGmpEigenMatrix result;

    if (isComplex) {
        result = complexIsometry().inv().complexIsometryInverse();
    } else {
        // Since Eigen does not provide an inverse function for sparse matrices,
        // we use an LU solver to compute the full form of each column of the
        // inverse matrix one by one (this avoids storing big matrices...).

        //matrixR.makeCompressed(); // The matrix is supposed to be alread compressed at this stage...
        SparseLU< SparseMatrix<mpreal, ColMajor>, COLAMDOrdering<int> > solver;

        solver.analyzePattern(matrixR);
        solver.factorize(matrixR);

        Matrix<mpreal, Dynamic, 1> I(matrixR.rows(), 1);
        GmpEigenMatrix X;

        // We start with the first column
        I(0,0) = 1;
        X.matrixR = solver.solve(I);
        result = SparseGmpEigenMatrix(X);

        // The next columns
        for (IndexType i(1); i < matrixR.rows(); ++i) {
            I(i-1,0) = 0;
            I(i,0) = 1;
            X.matrixR = solver.solve(I);
            result = result.horzcat(SparseGmpEigenMatrix(X));
        }
    }

    return result;
}

/* Matrix inverse b = inv(a) */
SparseGmpEigenMatrix& SparseGmpEigenMatrix::inv_new() const
{
    SparseGmpEigenMatrix& result(*(new SparseGmpEigenMatrix));

    if (isComplex) {
        result = complexIsometry().inv().complexIsometryInverse();
    } else {
        // Since Eigen does not provide an inverse function for sparse matrices,
        // we use an LU solver to compute the full form of each column of the
        // inverse matrix one by one (this avoids storing big matrices...).

        //matrixR.makeCompressed(); // The matrix is supposed to be alread compressed at this stage...
        SparseLU< SparseMatrix<mpreal, ColMajor>, COLAMDOrdering<int> > solver;

        solver.analyzePattern(matrixR);
        solver.factorize(matrixR);

        Matrix<mpreal, Dynamic, 1> I(matrixR.rows(), 1);
        GmpEigenMatrix X;

        // We start with the first column
        I(0,0) = 1;
        X.matrixR = solver.solve(I);
        result = SparseGmpEigenMatrix(X);

        // The next columns
        for (IndexType i(1); i < matrixR.rows(); ++i) {
            I(i-1,0) = 0;
            I(i,0) = 1;
            X.matrixR = solver.solve(I);
            result = result.horzcat(SparseGmpEigenMatrix(X));
        }
    }

    return result;
}

GmpEigenMatrix SparseGmpEigenMatrix::eigs(const long int& nbEigenvalues, GmpEigenMatrix& V, const long int& type, const GmpEigenMatrix& sigma) const
{
    GmpEigenMatrix result;

    // NOTE : At the moment, it seems that spectra only garantees of precision
    // of eps^(2/3). Hence only ~34 digits of precision with 50-digits mpreals

    // NOTE : In presence of degenerate eigenvalues, it seems that the second
    // largest eigenvalue can sometimes be found before the second copy of the
    // largest eigenvalue...

    if (ishermitian()) {
        if (isComplex) {
            IndexType size(matrixR.rows());

            // First, we solve a real version of the problem
            GmpEigenMatrix Dreal, Vreal; // These matrices will hold an isometry of the eigenvalues and eigenvectors
            Dreal = complexIsometry().eigs(2*nbEigenvalues, Vreal, type, sigma);

            // Since the eigenvalues are sorted and are supposed to come by pair,
            // we can easily remove duplicates. So we select the eigenvalues and
            // form the corresponding eigenvectors
            result.isComplex = false;
            result.matrixR.setZero(nbEigenvalues, nbEigenvalues);
            result.matrixI.resize(0,0);
            V.matrixR.setZero(size,nbEigenvalues);
            V.matrixI.setZero(size,nbEigenvalues);
            for (IndexType i(0); i < nbEigenvalues; ++i) {
                result.matrixR(i,i) = Dreal.matrixR(2*i,2*i);
                // This way of extracting the eigenvector always works ;-)
                V.matrixR.block(0,i,size,1) = Vreal.matrixR.block(0,2*i,size,1);
                V.matrixI.block(0,i,size,1) = -Vreal.matrixR.block(size,2*i,size,1);
            }
            V.checkComplexity();
        } else {
            //long int ncv(min(max(1+nbEigenvalues,2*nbEigenvalues),matrixR.rows())); // the tightest bounds... don't always converge
            long int ncv(min(3+max(1+nbEigenvalues,2*nbEigenvalues),matrixR.rows()));
            switch (type) {
                case 1: {
                    // Construct matrix operation object using the wrapper class
                    Spectra::SparseSymMatProd<mpreal> op(matrixR);

                    // Construct eigen solver object, requesting desired eigenvalues
                    Spectra::SymEigsSolver< mpreal, Spectra::LARGEST_MAGN, Spectra::SparseSymMatProd<mpreal> > eigs(&op, nbEigenvalues, ncv);

                    // Initialize and compute
                    eigs.init();
                    int maxIter(1000);
                    mpreal tolerance(pow(10,-mpfr::bits2digits(mpfr::mpreal::get_default_prec())));
                    int nconv = eigs.compute(maxIter, tolerance, Spectra::LARGEST_MAGN);

                    // Check for error
                    if(eigs.info() != Spectra::SUCCESSFUL)
                        mexErrMsgTxt("Eigenvalue decomposition failed.");

                    // Retrieve results
                    result.isComplex = false;
                    result.matrixR = eigs.eigenvalues();
                    result.matrixI.resize(0,0);
                    result = result.diagCreate(0);

                    V.isComplex = false;
                    V.matrixR = eigs.eigenvectors();
                    V.matrixI.resize(0,0);

                    break;
                }
                case 2: {
                    // Construct matrix operation object using the wrapper class
                    Spectra::SparseSymShiftSolve<mpreal> op(matrixR);

                    // Construct eigen solver object, requesting desired eigenvalues
                    Spectra::SymEigsShiftSolver< mpreal, Spectra::LARGEST_MAGN, Spectra::SparseSymShiftSolve<mpreal> > eigs(&op, nbEigenvalues, ncv, sigma.matrixR(0,0));

                    // Initialize and compute
                    eigs.init();
                    int maxIter(1000);
                    mpreal tolerance(pow(10,-mpfr::bits2digits(mpfr::mpreal::get_default_prec())));
                    int nconv = eigs.compute(maxIter, tolerance, Spectra::SMALLEST_MAGN);

                    // Check for error
                    if(eigs.info() != Spectra::SUCCESSFUL)
                        mexErrMsgTxt("Eigenvalue decomposition failed.");

                    // Retrieve results
                    result.isComplex = false;
                    result.matrixR = eigs.eigenvalues();
                    result.matrixI.resize(0,0);
                    result = result.diagCreate(0);

                    V.isComplex = false;
                    V.matrixR = eigs.eigenvectors();
                    V.matrixI.resize(0,0);

                    break;
                }
            }
        }
    } else {
        if (isComplex) {
            IndexType size(matrixR.rows());

            // First, we solve a real version of the problem
            GmpEigenMatrix Dreal, Vreal; // These matrices will hold an isometry of the eigenvalues and eigenvectors
            Dreal = complexIsometry().eigs(2*nbEigenvalues, Vreal, type, sigma);

            // Since the eigenvalues are sorted and are supposed to come by pair,
            // we now remove duplicates. Just like in the eig function, we need
            // to take special care for complex eigenvalues : they come in
            // conjugate form, but only one of them is correct. For pairs of
            // real eigenvalues, any of the two eigenvalues can be chosen
            // equivalently.

            // So we select the eigenvalues and
            // form the corresponding eigenvectors
            result.isComplex = false;
            result.matrixR.setZero(nbEigenvalues, nbEigenvalues);
            result.matrixI.setZero(nbEigenvalues, nbEigenvalues);
            V.matrixR.setZero(size,nbEigenvalues);
            V.matrixI.setZero(size,nbEigenvalues);
            for (IndexType i(0); i < nbEigenvalues; ++i) {
                // This way of extracting the eigenvector always works ;-)
                V.matrixR.block(0,i,size,1) = Vreal.matrixR.block(0,2*i,size,1);
                V.matrixI.block(0,i,size,1) = -Vreal.matrixR.block(size,2*i,size,1);
                if (Vreal.isComplex) {
                    if ( (Vreal.matrixR.block(0,2*i,size,1)+Vreal.matrixI.block(size,2*i,size,1)).array().abs().maxCoeff() >= (Vreal.matrixI.block(0,2*i,size,1)+Vreal.matrixR.block(size,2*i,size,1)).array().abs().maxCoeff() ) {
                        // The eigenvalue is correct
                        result.matrixR(i,i) = Dreal.matrixR(2*i,2*i);
                        if (Dreal.isComplex)
                            result.matrixI(i,i) = Dreal.matrixI(2*i,2*i);
                    } else {
                        // The eigenvalue needs to be conjugated
                        // (it should also appear in the next column, but to be
                        // safe we compute it from this column)
                        result.matrixR(i,i) = Dreal.matrixR(2*i,2*i);
                        if (Dreal.isComplex)
                            result.matrixI(i,i) = -Dreal.matrixI(2*i,2*i);
                    }
                } else {
                    // We just pick the eigenvalue on the diagonal
                    result.matrixR(i,i) = Dreal.matrixR(2*i,2*i);
                    if (Dreal.isComplex)
                        result.matrixI(i,i) = Dreal.matrixI(2*i,2*i);
                }
            }
            result.checkComplexity();
            V.checkComplexity();
        } else {
            //long int ncv(min(max(2+nbEigenvalues,2*nbEigenvalues),matrixR.rows())); // the tightest bounds... don't always converge
            long int ncv(min(3+max(2+nbEigenvalues,2*nbEigenvalues),matrixR.rows()));
            switch (type) {
                case 1: {
                    // Construct matrix operation object using the wrapper class
                    Spectra::SparseGenMatProd<mpreal> op(matrixR);

                    // Construct eigen solver object, requesting desired eigenvalues
                    Spectra::GenEigsSolver< mpreal, Spectra::LARGEST_MAGN, Spectra::SparseGenMatProd<mpreal> > eigs(&op, nbEigenvalues, ncv);

                    // Initialize and compute
                    eigs.init();
                    int maxIter(1000);
                    mpreal tolerance(pow(10,-mpfr::bits2digits(mpfr::mpreal::get_default_prec())));
                    int nconv = eigs.compute(maxIter, tolerance, Spectra::LARGEST_MAGN);

                    // Check for error
                    if(eigs.info() != Spectra::SUCCESSFUL)
                        mexErrMsgTxt("Eigenvalue decomposition failed.");

                    // Retrieve results
                    Matrix< complex<mpreal>, Dynamic, 1> evalues(eigs.eigenvalues());
                    Matrix< complex<mpreal>, Dynamic, Dynamic> evectors(eigs.eigenvectors());

                    result.matrixR = evalues.real();
                    result.matrixI = evalues.imag();
                    result.checkComplexity();
                    result = result.diagCreate(0);

                    V.matrixR = evectors.real();
                    V.matrixI = evectors.imag();
                    V.checkComplexity();

                    break;
                }
                case 2: {
                    // Construct matrix operation object using the wrapper class
                    Spectra::SparseGenRealShiftSolve<mpreal> op(matrixR);

                    // Construct eigen solver object, requesting desired eigenvalues
                    Spectra::GenEigsRealShiftSolver< mpreal, Spectra::LARGEST_MAGN, Spectra::SparseGenRealShiftSolve<mpreal> > eigs(&op, nbEigenvalues, ncv, sigma.matrixR(0,0));

                    // Initialize and compute
                    eigs.init();
                    int maxIter(1000);
                    mpreal tolerance(pow(10,-mpfr::bits2digits(mpfr::mpreal::get_default_prec())));
                    int nconv = eigs.compute(maxIter, tolerance, Spectra::SMALLEST_MAGN);

                    // Check for error
                    if(eigs.info() != Spectra::SUCCESSFUL)
                        mexErrMsgTxt("Eigenvalue decomposition failed.");

                    // Retrieve results
                    Matrix< complex<mpreal>, Dynamic, 1> evalues(eigs.eigenvalues());
                    Matrix< complex<mpreal>, Dynamic, Dynamic> evectors(eigs.eigenvectors());

                    result.matrixR = evalues.real();
                    result.matrixI = evalues.imag();
                    result.checkComplexity();
                    result = result.diagCreate(0);

                    V.matrixR = evectors.real();
                    V.matrixI = evectors.imag();
                    V.checkComplexity();

                    break;
                }
            }
        }
    }

    return result;
}

GmpEigenMatrix& SparseGmpEigenMatrix::eigs_new(const long int& nbEigenvalues, GmpEigenMatrix& V, const long int& type, const GmpEigenMatrix& sigma) const
{
    GmpEigenMatrix& result(*(new GmpEigenMatrix));

    // NOTE : At the moment, it seems that spectra only garantees of precision
    // of eps^(2/3). Hence only ~34 digits of precision with 50-digits mpreals

    // NOTE : In presence of degenerate eigenvalues, it seems that the second
    // largest eigenvalue can sometimes be found before the second copy of the
    // largest eigenvalue...

    if (ishermitian()) {
        if (isComplex) {
            IndexType size(matrixR.rows());

            // First, we solve a real version of the problem
            GmpEigenMatrix Dreal, Vreal; // These matrices will hold an isometry of the eigenvalues and eigenvectors
            Dreal = complexIsometry().eigs(2*nbEigenvalues, Vreal, type, sigma);

            // Since the eigenvalues are sorted and are supposed to come by pair,
            // we can easily remove duplicates. So we select the eigenvalues and
            // form the corresponding eigenvectors
            result.isComplex = false;
            result.matrixR.setZero(nbEigenvalues, nbEigenvalues);
            result.matrixI.resize(0,0);
            V.matrixR.setZero(size,nbEigenvalues);
            V.matrixI.setZero(size,nbEigenvalues);
            for (IndexType i(0); i < nbEigenvalues; ++i) {
                result.matrixR(i,i) = Dreal.matrixR(2*i,2*i);
                // This way of extracting the eigenvector always works ;-)
                V.matrixR.block(0,i,size,1) = Vreal.matrixR.block(0,2*i,size,1);
                V.matrixI.block(0,i,size,1) = -Vreal.matrixR.block(size,2*i,size,1);
            }
            V.checkComplexity();
        } else {
            //long int ncv(min(max(1+nbEigenvalues,2*nbEigenvalues),matrixR.rows())); // the tightest bounds... don't always converge
            long int ncv(min(3+max(1+nbEigenvalues,2*nbEigenvalues),matrixR.rows()));
            switch (type) {
                case 1: {
                    // Construct matrix operation object using the wrapper class
                    Spectra::SparseSymMatProd<mpreal> op(matrixR);

                    // Construct eigen solver object, requesting desired eigenvalues
                    Spectra::SymEigsSolver< mpreal, Spectra::LARGEST_MAGN, Spectra::SparseSymMatProd<mpreal> > eigs(&op, nbEigenvalues, ncv);

                    // Initialize and compute
                    eigs.init();
                    int maxIter(1000);
                    mpreal tolerance(pow(10,-mpfr::bits2digits(mpfr::mpreal::get_default_prec())));
                    int nconv = eigs.compute(maxIter, tolerance, Spectra::LARGEST_MAGN);

                    // Check for error
                    if(eigs.info() != Spectra::SUCCESSFUL)
                        mexErrMsgTxt("Eigenvalue decomposition failed.");

                    // Retrieve results
                    result.isComplex = false;
                    result.matrixR = eigs.eigenvalues();
                    result.matrixI.resize(0,0);
                    result = result.diagCreate(0);

                    V.isComplex = false;
                    V.matrixR = eigs.eigenvectors();
                    V.matrixI.resize(0,0);

                    break;
                }
                case 2: {
                    // Construct matrix operation object using the wrapper class
                    Spectra::SparseSymShiftSolve<mpreal> op(matrixR);

                    // Construct eigen solver object, requesting desired eigenvalues
                    Spectra::SymEigsShiftSolver< mpreal, Spectra::LARGEST_MAGN, Spectra::SparseSymShiftSolve<mpreal> > eigs(&op, nbEigenvalues, ncv, sigma.matrixR(0,0));

                    // Initialize and compute
                    eigs.init();
                    int maxIter(1000);
                    mpreal tolerance(pow(10,-mpfr::bits2digits(mpfr::mpreal::get_default_prec())));
                    int nconv = eigs.compute(maxIter, tolerance, Spectra::SMALLEST_MAGN);

                    // Check for error
                    if(eigs.info() != Spectra::SUCCESSFUL)
                        mexErrMsgTxt("Eigenvalue decomposition failed.");

                    // Retrieve results
                    result.isComplex = false;
                    result.matrixR = eigs.eigenvalues();
                    result.matrixI.resize(0,0);
                    result = result.diagCreate(0);

                    V.isComplex = false;
                    V.matrixR = eigs.eigenvectors();
                    V.matrixI.resize(0,0);

                    break;
                }
            }
        }
    } else {
        if (isComplex) {
            IndexType size(matrixR.rows());

            // First, we solve a real version of the problem
            GmpEigenMatrix Dreal, Vreal; // These matrices will hold an isometry of the eigenvalues and eigenvectors
            Dreal = complexIsometry().eigs(2*nbEigenvalues, Vreal, type, sigma);

            // Since the eigenvalues are sorted and are supposed to come by pair,
            // we now remove duplicates. Just like in the eig function, we need
            // to take special care for complex eigenvalues : they come in
            // conjugate form, but only one of them is correct. For pairs of
            // real eigenvalues, any of the two eigenvalues can be chosen
            // equivalently.

            // So we select the eigenvalues and
            // form the corresponding eigenvectors
            result.isComplex = false;
            result.matrixR.setZero(nbEigenvalues, nbEigenvalues);
            result.matrixI.setZero(nbEigenvalues, nbEigenvalues);
            V.matrixR.setZero(size,nbEigenvalues);
            V.matrixI.setZero(size,nbEigenvalues);
            for (IndexType i(0); i < nbEigenvalues; ++i) {
                // This way of extracting the eigenvector always works ;-)
                V.matrixR.block(0,i,size,1) = Vreal.matrixR.block(0,2*i,size,1);
                V.matrixI.block(0,i,size,1) = -Vreal.matrixR.block(size,2*i,size,1);
                if (Vreal.isComplex) {
                    if ( (Vreal.matrixR.block(0,2*i,size,1)+Vreal.matrixI.block(size,2*i,size,1)).array().abs().maxCoeff() >= (Vreal.matrixI.block(0,2*i,size,1)+Vreal.matrixR.block(size,2*i,size,1)).array().abs().maxCoeff() ) {
                        // The eigenvalue is correct
                        result.matrixR(i,i) = Dreal.matrixR(2*i,2*i);
                        if (Dreal.isComplex)
                            result.matrixI(i,i) = Dreal.matrixI(2*i,2*i);
                    } else {
                        // The eigenvalue needs to be conjugated
                        // (it should also appear in the next column, but to be
                        // safe we compute it from this column)
                        result.matrixR(i,i) = Dreal.matrixR(2*i,2*i);
                        if (Dreal.isComplex)
                            result.matrixI(i,i) = -Dreal.matrixI(2*i,2*i);
                    }
                } else {
                    // We just pick the eigenvalue on the diagonal
                    result.matrixR(i,i) = Dreal.matrixR(2*i,2*i);
                    if (Dreal.isComplex)
                        result.matrixI(i,i) = Dreal.matrixI(2*i,2*i);
                }
            }
            result.checkComplexity();
            V.checkComplexity();
        } else {
            //long int ncv(min(max(2+nbEigenvalues,2*nbEigenvalues),matrixR.rows())); // the tightest bounds... don't always converge
            long int ncv(min(3+max(2+nbEigenvalues,2*nbEigenvalues),matrixR.rows()));
            switch (type) {
                case 1: {
                    // Construct matrix operation object using the wrapper class
                    Spectra::SparseGenMatProd<mpreal> op(matrixR);

                    // Construct eigen solver object, requesting desired eigenvalues
                    Spectra::GenEigsSolver< mpreal, Spectra::LARGEST_MAGN, Spectra::SparseGenMatProd<mpreal> > eigs(&op, nbEigenvalues, ncv);

                    // Initialize and compute
                    eigs.init();
                    int maxIter(1000);
                    mpreal tolerance(pow(10,-mpfr::bits2digits(mpfr::mpreal::get_default_prec())));
                    int nconv = eigs.compute(maxIter, tolerance, Spectra::LARGEST_MAGN);

                    // Check for error
                    if(eigs.info() != Spectra::SUCCESSFUL)
                        mexErrMsgTxt("Eigenvalue decomposition failed.");

                    // Retrieve results
                    Matrix< complex<mpreal>, Dynamic, 1> evalues(eigs.eigenvalues());
                    Matrix< complex<mpreal>, Dynamic, Dynamic> evectors(eigs.eigenvectors());

                    result.matrixR = evalues.real();
                    result.matrixI = evalues.imag();
                    result.checkComplexity();
                    result = result.diagCreate(0);

                    V.matrixR = evectors.real();
                    V.matrixI = evectors.imag();
                    V.checkComplexity();

                    break;
                }
                case 2: {
                    // Construct matrix operation object using the wrapper class
                    Spectra::SparseGenRealShiftSolve<mpreal> op(matrixR);

                    // Construct eigen solver object, requesting desired eigenvalues
                    Spectra::GenEigsRealShiftSolver< mpreal, Spectra::LARGEST_MAGN, Spectra::SparseGenRealShiftSolve<mpreal> > eigs(&op, nbEigenvalues, ncv, sigma.matrixR(0,0));

                    // Initialize and compute
                    eigs.init();
                    int maxIter(1000);
                    mpreal tolerance(pow(10,-mpfr::bits2digits(mpfr::mpreal::get_default_prec())));
                    int nconv = eigs.compute(maxIter, tolerance, Spectra::SMALLEST_MAGN);

                    // Check for error
                    if(eigs.info() != Spectra::SUCCESSFUL)
                        mexErrMsgTxt("Eigenvalue decomposition failed.");

                    // Retrieve results
                    Matrix< complex<mpreal>, Dynamic, 1> evalues(eigs.eigenvalues());
                    Matrix< complex<mpreal>, Dynamic, Dynamic> evectors(eigs.eigenvectors());

                    result.matrixR = evalues.real();
                    result.matrixI = evalues.imag();
                    result.checkComplexity();
                    result = result.diagCreate(0);

                    V.matrixR = evectors.real();
                    V.matrixI = evectors.imag();
                    V.checkComplexity();

                    break;
                }
            }
        }
    }

    return result;
}







/* -----------------------------
   | Some comparison operators |
   ----------------------------- */

// The comparison operators only care about the real part of the matrices
SparseMatrix <bool> SparseGmpEigenMatrix::operator<(const SparseGmpEigenMatrix& b) const
{
    SparseMatrix <bool> result;

    if ((numel() != 1) && (b.numel() == 1)) {
        // setting the output size
        result.resize(matrixR.rows(), matrixR.cols());

        // We will first copy the data to triplets
        vector< Triplet<bool> > tripletList;
        tripletList.reserve(matrixR.nonZeros()); // WARNING : Here, we assume that the real part of the scalar b is smaller or equal to 0 (otherwise it will make the sparse matrix full...)

        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);

            while (itR) {
                if (itR.value() < b.matrixR.coeff(0,0))
                    tripletList.push_back(Triplet<bool>(itR.row(), k, true));
                ++itR;
            }
        }

        // Now we can assign the data to the sparse matrix
        result.setFromTriplets(tripletList.begin(), tripletList.end());
    } else if ((numel() == 1) && (b.numel() != 1)) {
        result = b.operator>(*this);
    } else {
        // setting the output size
        result.resize(matrixR.rows(), matrixR.cols());

        // We will first copy the data to triplets
        vector< Triplet<bool> > tripletList;
        tripletList.reserve(matrixR.nonZeros() + b.matrixR.nonZeros());

        // Now we perform the comparison between the two matrices
        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
            SparseMatrix<mpreal>::InnerIterator itRb(b.matrixR,k);

            while ((itR) || (itRb)) {
                if ((!itRb) || ((itR) && (itR.row() < itRb.row()))) {
                    if (itR.value() < 0)
                        tripletList.push_back(Triplet<bool>(itR.row(), k, true));
                    ++itR;
                } else if ((itR) && (itRb) && (itR.row() == itRb.row())) {
                    if (itR.value() < itRb.value())
                        tripletList.push_back(Triplet<bool>(itR.row(), k, true));
                    ++itRb;
                    ++itR;
                } else {
                    if (0 < itRb.value())
                        tripletList.push_back(Triplet<bool>(itRb.row(), k, true));
                    ++itRb;
                }
            }
        }

        // Now we can assign the data to the sparse matrix
        result.setFromTriplets(tripletList.begin(), tripletList.end());
    }

    return result;
}

SparseMatrix <bool> SparseGmpEigenMatrix::operator<=(const SparseGmpEigenMatrix& b) const
{
    SparseMatrix <bool> result;

    if ((numel() != 1) && (b.numel() == 1)) {
        // setting the output size
        result.resize(matrixR.rows(), matrixR.cols());

        // We will first copy the data to triplets
        vector< Triplet<bool> > tripletList;
        tripletList.reserve(matrixR.nonZeros()); // WARNING : Here, we assume that the real part of the scalar b is smaller than 0 (otherwise it will make the sparse matrix full...)

        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);

            while (itR) {
                if (itR.value() <= b.matrixR.coeff(0,0))
                    tripletList.push_back(Triplet<bool>(itR.row(), k, true));
                ++itR;
            }
        }

        // Now we can assign the data to the sparse matrix
        result.setFromTriplets(tripletList.begin(), tripletList.end());
    } else if ((numel() == 1) && (b.numel() != 1)) {
        result = b.operator>=(*this);
    } else {
        // setting the output size
        result.resize(matrixR.rows(), matrixR.cols());

        // We will first copy the data to triplets
        vector< Triplet<bool> > tripletList;
        tripletList.reserve(matrixR.nonZeros() + b.matrixR.nonZeros());

        // Now we perform the comparison between the two matrices
        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
            SparseMatrix<mpreal>::InnerIterator itRb(b.matrixR,k);

            while ((itR) || (itRb)) {
                if ((!itRb) || ((itR) && (itR.row() < itRb.row()))) {
                    if (itR.value() <= 0)
                        tripletList.push_back(Triplet<bool>(itR.row(), k, true));
                    ++itR;
                } else if ((itR) && (itRb) && (itR.row() == itRb.row())) {
                    if (itR.value() <= itRb.value())
                        tripletList.push_back(Triplet<bool>(itR.row(), k, true));
                    ++itRb;
                    ++itR;
                } else {
                    if (0 <= itRb.value())
                        tripletList.push_back(Triplet<bool>(itRb.row(), k, true));
                    ++itRb;
                }
            }
        }

        // Now we can assign the data to the sparse matrix
        result.setFromTriplets(tripletList.begin(), tripletList.end());
    }

    return result;
}

SparseMatrix <bool> SparseGmpEigenMatrix::operator>(const SparseGmpEigenMatrix& b) const
{
    SparseMatrix <bool> result;

    if ((numel() != 1) && (b.numel() == 1)) {
        // setting the output size
        result.resize(matrixR.rows(), matrixR.cols());

        // We will first copy the data to triplets
        vector< Triplet<bool> > tripletList;
        tripletList.reserve(matrixR.nonZeros()); // WARNING : Here, we assume that the real part of the scalar b is larger or equal to 0 (otherwise it will make the sparse matrix full...)

        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);

            while (itR) {
                if (itR.value() > b.matrixR.coeff(0,0))
                    tripletList.push_back(Triplet<bool>(itR.row(), k, true));
                ++itR;
            }
        }

        // Now we can assign the data to the sparse matrix
        result.setFromTriplets(tripletList.begin(), tripletList.end());
    } else if ((numel() == 1) && (b.numel() != 1)) {
        result = b.operator<(*this);
    } else {
        // setting the output size
        result.resize(matrixR.rows(), matrixR.cols());

        // We will first copy the data to triplets
        vector< Triplet<bool> > tripletList;
        tripletList.reserve(matrixR.nonZeros() + b.matrixR.nonZeros());

        // Now we perform the comparison between the two matrices
        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
            SparseMatrix<mpreal>::InnerIterator itRb(b.matrixR,k);

            while ((itR) || (itRb)) {
                if ((!itRb) || ((itR) && (itR.row() < itRb.row()))) {
                    if (itR.value() > 0)
                        tripletList.push_back(Triplet<bool>(itR.row(), k, true));
                    ++itR;
                } else if ((itR) && (itRb) && (itR.row() == itRb.row())) {
                    if (itR.value() > itRb.value())
                        tripletList.push_back(Triplet<bool>(itR.row(), k, true));
                    ++itRb;
                    ++itR;
                } else {
                    if (0 > itRb.value())
                        tripletList.push_back(Triplet<bool>(itRb.row(), k, true));
                    ++itRb;
                }
            }
        }

        // Now we can assign the data to the sparse matrix
        result.setFromTriplets(tripletList.begin(), tripletList.end());
    }

    return result;
}

SparseMatrix <bool> SparseGmpEigenMatrix::operator>=(const SparseGmpEigenMatrix& b) const
{
    SparseMatrix <bool> result;

    if ((numel() != 1) && (b.numel() == 1)) {
        // setting the output size
        result.resize(matrixR.rows(), matrixR.cols());

        // We will first copy the data to triplets
        vector< Triplet<bool> > tripletList;
        tripletList.reserve(matrixR.nonZeros()); // WARNING : Here, we assume that the real part of the scalar b is larger than 0 (otherwise it will make the sparse matrix full...)

        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);

            while (itR) {
                if (itR.value() >= b.matrixR.coeff(0,0))
                    tripletList.push_back(Triplet<bool>(itR.row(), k, true));
                ++itR;
            }
        }

        // Now we can assign the data to the sparse matrix
        result.setFromTriplets(tripletList.begin(), tripletList.end());
    } else if ((numel() == 1) && (b.numel() != 1)) {
        result = b.operator<=(*this);
    } else {
        // setting the output size
        result.resize(matrixR.rows(), matrixR.cols());

        // We will first copy the data to triplets
        vector< Triplet<bool> > tripletList;
        tripletList.reserve(matrixR.nonZeros() + b.matrixR.nonZeros());

        // Now we perform the comparison between the two matrices
        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
            SparseMatrix<mpreal>::InnerIterator itRb(b.matrixR,k);

            while ((itR) || (itRb)) {
                if ((!itRb) || ((itR) && (itR.row() < itRb.row()))) {
                    if (itR.value() >= 0)
                        tripletList.push_back(Triplet<bool>(itR.row(), k, true));
                    ++itR;
                } else if ((itR) && (itRb) && (itR.row() == itRb.row())) {
                    if (itR.value() >= itRb.value())
                        tripletList.push_back(Triplet<bool>(itR.row(), k, true));
                    ++itRb;
                    ++itR;
                } else {
                    if (0 >= itRb.value())
                        tripletList.push_back(Triplet<bool>(itRb.row(), k, true));
                    ++itRb;
                }
            }
        }

        // Now we can assign the data to the sparse matrix
        result.setFromTriplets(tripletList.begin(), tripletList.end());
    }

    return result;
}

/* Test of equality  c = (a == b)
   NOTE : Here we make some explicit use of the assumption that any element
   that is indexed in a sparse matrix is nonzero. We could do much more with
   this assumption... */
SparseMatrix <bool> SparseGmpEigenMatrix::eq(const SparseGmpEigenMatrix& b) const
{
    SparseMatrix <bool> result;

    // We want to support comparison between an nxm matrix and a 1x1 matrix
    if ((numel() != 1) && (b.numel() == 1)) {
        // setting the output size
        result.resize(matrixR.rows(), matrixR.cols());

        // We will first copy the data to triplets
        vector< Triplet<bool> > tripletList;

        if (b.isComplex) {
            if (isComplex) {
                tripletList.reserve(matrixR.nonZeros() + matrixI.nonZeros()); // WARNING : Here, we assume that either the real part or the imaginary part of the scalar b is not 0 (otherwise it will make the sparse matrix full...)

                for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
                    SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
                    SparseMatrix<mpreal>::InnerIterator itI(matrixI,k);

                    while ((itR) || (itI)) {
                        if ((!itI) || ((itR) && (itR.row() < itI.row()))) {
                            if ((itR.value() == b.matrixR.coeff(0,0)) && (0 == b.matrixI.coeff(0,0)))
                                tripletList.push_back(Triplet<bool>(itR.row(), k, true));
                            ++itR;
                        } else if ((itR) && (itI) && (itR.row() == itI.row())) {
                            if ((itR.value() == b.matrixR.coeff(0,0)) && (itI.value() == b.matrixI.coeff(0,0)))
                                tripletList.push_back(Triplet<bool>(itR.row(), k, true));
                            ++itR;
                            ++itI;
                        } else {
                            if ((0 == b.matrixR.coeff(0,0)) && (itI.value() == b.matrixI.coeff(0,0)))
                                tripletList.push_back(Triplet<bool>(itI.row(), k, true));
                            ++itI;
                        }
                    }
                }
            } else {
                // The matrix is real, but the element to be compared to is complex, so no element in the matrix equals the scalar : the comparison matrix is empty.
                if (b.matrixI.coeff(0,0) == 0)
                    // We give an error if the above assessement is not correct.
                    mexErrMsgTxt("Error in SparseGmpEigenMatrix::eq while copmaring to a complex scalar: the complex scalar has null imaginary part.");
            }
        } else {
            if (isComplex) {
                tripletList.reserve(matrixR.nonZeros()); // WARNING : Here, we assume that the real part of the scalar b is not 0 (otherwise it will make the sparse matrix full...)

                for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
                    SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
                    SparseMatrix<mpreal>::InnerIterator itI(matrixI,k);

                    while ((itR) || (itI)) {
                        if ((!itI) || ((itR) && (itR.row() < itI.row()))) {
                            if (itR.value() == b.matrixR.coeff(0,0))
                                tripletList.push_back(Triplet<bool>(itR.row(), k, true));
                            ++itR;
                        } else if ((itR) && (itI) && (itR.row() == itI.row())) {
                            if ((itR.value() == b.matrixR.coeff(0,0)) && (itI.value() == 0))
                                tripletList.push_back(Triplet<bool>(itR.row(), k, true));
                            ++itR;
                            ++itI;
                        } else {
                            if ((0 == b.matrixR.coeff(0,0)) && (itI.value() == 0))
                                tripletList.push_back(Triplet<bool>(itI.row(), k, true));
                            ++itI;
                        }
                    }
                }
            } else {
                tripletList.reserve(matrixR.nonZeros()); // WARNING : Here, we assume that the real part of the scalar b is not 0 (otherwise it will make the sparse matrix full...)

                for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
                    SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);

                    while (itR) {
                        if (itR.value() == b.matrixR.coeff(0,0))
                            tripletList.push_back(Triplet<bool>(itR.row(), k, true));
                        ++itR;
                    }
                }
            }
        }

        // Now we can assign the data to the sparse matrix
        result.setFromTriplets(tripletList.begin(), tripletList.end());
    } else if ((numel() == 1) && (b.numel() != 1)) {
        result = b.eq(*this);
    } else {
        // setting the output size
        result.resize(matrixR.rows(), matrixR.cols());

        // We will first copy the data to triplets
        vector< Triplet<bool> > tripletList;
        tripletList.reserve(numel()); // It is hard to estimate the number of elements that will be created here : any two zeros will give a one...
        //tripletList.reserve(matrixR.nonZeros() + b.matrixR.nonZeros() + max(0,numel()-matrixR.nonZeros()-b.matrixR.nonZeros())); // WARNING : This may be a bad estimate if most of both matrices are full of non-zeros... (in such a case, however, the user should prefer converting the matrices to full ones first...)

        if (b.isComplex) {
            if (isComplex) {
                for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
                    SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
                    SparseMatrix<mpreal>::InnerIterator itI(matrixI,k);
                    SparseMatrix<mpreal>::InnerIterator itRb(b.matrixR,k);
                    SparseMatrix<mpreal>::InnerIterator itIb(b.matrixI,k);
                    IndexType previousNonZero(-1);

                    while ((itR) || (itI) || (itRb) || (itIb)) {
                        IndexType row, rowb;
                        if ((itR) || (itI)) {
                            if ((itR) && (itI))
                                row = min(itR.row(), itI.row());
                            else if (itR)
                                row = itR.row();
                            else
                                row = itI.row();
                        }
                        if ((itRb) || (itIb)) {
                            if ((itRb) && (itIb))
                                rowb = min(itRb.row(), itIb.row());
                            else if (itRb)
                                rowb = itRb.row();
                            else
                                rowb = itIb.row();
                        }

                        if (((((itR) || (itI)) && ((itRb) || (itIb))) && (row < rowb)) || ((!itRb) && (!itIb))) {
                            if ((!itI) || ((itR) && (itR.row() < itI.row()))) {
                                for (IndexType i(previousNonZero+1); i < itR.row(); ++i)
                                    tripletList.push_back(Triplet<bool>(i, k, true));
                                previousNonZero = itR.row();
                                if (itR.value() == 0)
                                    tripletList.push_back(Triplet<bool>(itR.row(), k, true));
                                ++itR;
                            } else if ((itR) && (itI) && (itR.row() == itI.row())) {
                                for (IndexType i(previousNonZero+1); i < itR.row(); ++i)
                                    tripletList.push_back(Triplet<bool>(i, k, true));
                                previousNonZero = itR.row();
                                if ((itR.value() == 0) && (itI.value() == 0))
                                    tripletList.push_back(Triplet<bool>(itR.row(), k, true));
                                ++itR;
                                ++itI;
                            } else {
                                for (IndexType i(previousNonZero+1); i < itI.row(); ++i)
                                    tripletList.push_back(Triplet<bool>(i, k, true));
                                previousNonZero = itI.row();
                                if (itI.value() == 0)
                                    tripletList.push_back(Triplet<bool>(itI.row(), k, true));
                                ++itI;
                            }
                        } else if ((((itR) || (itI)) && ((itRb) || (itIb))) && (row == rowb)) {
                            if ((!itI) || ((itR) && (itR.row() < itI.row()))) {
                                if ((!itIb) || ((itRb) && (itRb.row() < itIb.row()))) {
                                    for (IndexType i(previousNonZero+1); i < itR.row(); ++i)
                                        tripletList.push_back(Triplet<bool>(i, k, true));
                                    previousNonZero = itR.row();
                                    if (itR.value() == itRb.value())
                                        tripletList.push_back(Triplet<bool>(itR.row(), k, true));
                                    ++itRb;
                                } else if ((itRb) && (itIb) && (itRb.row() == itIb.row())) {
                                    for (IndexType i(previousNonZero+1); i < itR.row(); ++i)
                                        tripletList.push_back(Triplet<bool>(i, k, true));
                                    previousNonZero = itR.row();
                                    if ((itR.value() == itRb.value()) && (0 == itIb.value()))
                                        tripletList.push_back(Triplet<bool>(itR.row(), k, true));
                                    ++itRb;
                                    ++itIb;
                                } else {
                                    for (IndexType i(previousNonZero+1); i < itR.row(); ++i)
                                        tripletList.push_back(Triplet<bool>(i, k, true));
                                    previousNonZero = itR.row();
                                    if (0 == itIb.value())
                                        tripletList.push_back(Triplet<bool>(itR.row(), k, true));
                                    ++itIb;
                                }
                                ++itR;
                            } else if ((itR) && (itI) && (itR.row() == itI.row())) {
                                if ((!itIb) || ((itRb) && (itRb.row() < itIb.row()))) {
                                    for (IndexType i(previousNonZero+1); i < itR.row(); ++i)
                                        tripletList.push_back(Triplet<bool>(i, k, true));
                                    previousNonZero = itR.row();
                                    if ((itR.value() == itRb.value()) && (itI.value() == 0))
                                        tripletList.push_back(Triplet<bool>(itR.row(), k, true));
                                    ++itRb;
                                } else if ((itRb) && (itIb) && (itRb.row() == itIb.row())) {
                                    for (IndexType i(previousNonZero+1); i < itR.row(); ++i)
                                        tripletList.push_back(Triplet<bool>(i, k, true));
                                    previousNonZero = itR.row();
                                    if ((itR.value() == itRb.value()) && (itI.value() == itIb.value()))
                                        tripletList.push_back(Triplet<bool>(itR.row(), k, true));
                                    ++itRb;
                                    ++itIb;
                                } else {
                                    for (IndexType i(previousNonZero+1); i < itR.row(); ++i)
                                        tripletList.push_back(Triplet<bool>(i, k, true));
                                    previousNonZero = itR.row();
                                    if ((itR.value() == 0) && (itI.value() == itIb.value()))
                                        tripletList.push_back(Triplet<bool>(itR.row(), k, true));
                                    ++itIb;
                                }
                                ++itR;
                                ++itI;
                            } else {
                                if ((!itIb) || ((itRb) && (itRb.row() < itIb.row()))) {
                                    for (IndexType i(previousNonZero+1); i < itI.row(); ++i)
                                        tripletList.push_back(Triplet<bool>(i, k, true));
                                    previousNonZero = itI.row();
                                    if ((0 == itRb.value()) && (itI.value() == 0))
                                        tripletList.push_back(Triplet<bool>(itI.row(), k, true));
                                    ++itRb;
                                } else if ((itRb) && (itIb) && (itRb.row() == itIb.row())) {
                                    for (IndexType i(previousNonZero+1); i < itI.row(); ++i)
                                        tripletList.push_back(Triplet<bool>(i, k, true));
                                    previousNonZero = itI.row();
                                    if ((0 == itRb.value()) && (itI.value() == itIb.value()))
                                        tripletList.push_back(Triplet<bool>(itI.row(), k, true));
                                    ++itRb;
                                    ++itIb;
                                } else {
                                    for (IndexType i(previousNonZero+1); i < itI.row(); ++i)
                                        tripletList.push_back(Triplet<bool>(i, k, true));
                                    previousNonZero = itI.row();
                                    if (itI.value() == itIb.value())
                                        tripletList.push_back(Triplet<bool>(itI.row(), k, true));
                                    ++itIb;
                                }
                                ++itI;
                            }
                        } else {
                            if ((!itIb) || ((itRb) && (itRb.row() < itIb.row()))) {
                                for (IndexType i(previousNonZero+1); i < itRb.row(); ++i)
                                    tripletList.push_back(Triplet<bool>(i, k, true));
                                previousNonZero = itRb.row();
                                if (0 == itRb.value())
                                    tripletList.push_back(Triplet<bool>(itRb.row(), k, true));
                                ++itRb;
                            } else if ((itRb) && (itIb) && (itRb.row() == itIb.row())) {
                                for (IndexType i(previousNonZero+1); i < itRb.row(); ++i)
                                    tripletList.push_back(Triplet<bool>(i, k, true));
                                previousNonZero = itRb.row();
                                if ((0 == itRb.value()) && (0 == itIb.value()))
                                    tripletList.push_back(Triplet<bool>(itRb.row(), k, true));
                                ++itRb;
                                ++itIb;
                            } else {
                                for (IndexType i(previousNonZero+1); i < itIb.row(); ++i)
                                    tripletList.push_back(Triplet<bool>(i, k, true));
                                previousNonZero = itIb.row();
                                if (0 == itIb.value())
                                    tripletList.push_back(Triplet<bool>(itIb.row(), k, true));
                                ++itIb;
                            }
                        }
                    }

                    // If there are still zeros
                    for (IndexType i(previousNonZero+1); i < matrixR.rows(); ++i)
                        tripletList.push_back(Triplet<bool>(i, k, true));
                }
            } else {
                for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
                    SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
                    SparseMatrix<mpreal>::InnerIterator itRb(b.matrixR,k);
                    SparseMatrix<mpreal>::InnerIterator itIb(b.matrixI,k);
                    IndexType previousNonZero(-1);

                    while ((itR) || (itRb) || (itIb)) {
                        IndexType row, rowb;
                        if (itR)
                            row = itR.row();
                        if ((itRb) || (itIb)) {
                            if ((itRb) && (itIb))
                                rowb = min(itRb.row(), itIb.row());
                            else if (itRb)
                                rowb = itRb.row();
                            else
                                rowb = itIb.row();
                        }

                        if ((((itR) && ((itRb) || (itIb))) && (row < rowb)) || ((!itRb) && (!itIb))) {
                            for (IndexType i(previousNonZero+1); i < itR.row(); ++i)
                                tripletList.push_back(Triplet<bool>(i, k, true));
                            previousNonZero = itR.row();
                            if (itR.value() == 0)
                                tripletList.push_back(Triplet<bool>(itR.row(), k, true));
                            ++itR;
                        } else if (((itR) && ((itRb) || (itIb))) && (row == rowb)) {
                            if ((!itIb) || ((itRb) && (itRb.row() < itIb.row()))) {
                                for (IndexType i(previousNonZero+1); i < itR.row(); ++i)
                                    tripletList.push_back(Triplet<bool>(i, k, true));
                                previousNonZero = itR.row();
                                if (itR.value() == itRb.value())
                                    tripletList.push_back(Triplet<bool>(itR.row(), k, true));
                                ++itRb;
                            } else if ((itRb) && (itIb) && (itRb.row() == itIb.row())) {
                                for (IndexType i(previousNonZero+1); i < itR.row(); ++i)
                                    tripletList.push_back(Triplet<bool>(i, k, true));
                                previousNonZero = itR.row();
                                if ((itR.value() == itRb.value()) && (0 == itIb.value()))
                                    tripletList.push_back(Triplet<bool>(itR.row(), k, true));
                                ++itRb;
                                ++itIb;
                            } else {
                                for (IndexType i(previousNonZero+1); i < itR.row(); ++i)
                                    tripletList.push_back(Triplet<bool>(i, k, true));
                                previousNonZero = itR.row();
                                if (0 == itIb.value())
                                    tripletList.push_back(Triplet<bool>(itR.row(), k, true));
                                ++itIb;
                            }
                            ++itR;
                        } else {
                            if ((!itIb) || ((itRb) && (itRb.row() < itIb.row()))) {
                                for (IndexType i(previousNonZero+1); i < itRb.row(); ++i)
                                    tripletList.push_back(Triplet<bool>(i, k, true));
                                previousNonZero = itRb.row();
                                if (0 == itRb.value())
                                    tripletList.push_back(Triplet<bool>(itRb.row(), k, true));
                                ++itRb;
                            } else if ((itRb) && (itIb) && (itRb.row() == itIb.row())) {
                                for (IndexType i(previousNonZero+1); i < itRb.row(); ++i)
                                    tripletList.push_back(Triplet<bool>(i, k, true));
                                previousNonZero = itRb.row();
                                if ((0 == itRb.value()) && (0 == itIb.value()))
                                    tripletList.push_back(Triplet<bool>(itRb.row(), k, true));
                                ++itRb;
                                ++itIb;
                            } else {
                                for (IndexType i(previousNonZero+1); i < itIb.row(); ++i)
                                    tripletList.push_back(Triplet<bool>(i, k, true));
                                previousNonZero = itIb.row();
                                if (0 == itIb.value())
                                    tripletList.push_back(Triplet<bool>(itIb.row(), k, true));
                                ++itIb;
                            }
                        }
                    }

                    // If there are still zeros
                    for (IndexType i(previousNonZero+1); i < matrixR.rows(); ++i)
                        tripletList.push_back(Triplet<bool>(i, k, true));
                }
            }
        } else if (isComplex) {
            for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
                SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
                SparseMatrix<mpreal>::InnerIterator itI(matrixI,k);
                SparseMatrix<mpreal>::InnerIterator itRb(b.matrixR,k);
                IndexType previousNonZero(-1);

                while ((itR) || (itI) || (itRb)) {
                    IndexType row, rowb;
                    if ((itR) || (itI)) {
                        if ((itR) && (itI))
                            row = min(itR.row(), itI.row());
                        else if (itR)
                            row = itR.row();
                        else
                            row = itI.row();
                    }
                    if (itRb)
                        rowb = itRb.row();

                    if (((((itR) || (itI)) && (itRb)) && (row < rowb)) || (!itRb)) {
                        if ((!itI) || ((itR) && (itR.row() < itI.row()))) {
                            for (IndexType i(previousNonZero+1); i < itR.row(); ++i)
                                tripletList.push_back(Triplet<bool>(i, k, true));
                            previousNonZero = itR.row();
                            if (itR.value() == 0)
                                tripletList.push_back(Triplet<bool>(itR.row(), k, true));
                            ++itR;
                        } else if ((itR) && (itI) && (itR.row() == itI.row())) {
                            for (IndexType i(previousNonZero+1); i < itR.row(); ++i)
                                tripletList.push_back(Triplet<bool>(i, k, true));
                            previousNonZero = itR.row();
                            if ((itR.value() == 0) && (itI.value() == 0))
                                tripletList.push_back(Triplet<bool>(itR.row(), k, true));
                            ++itR;
                            ++itI;
                        } else {
                            for (IndexType i(previousNonZero+1); i < itI.row(); ++i)
                                tripletList.push_back(Triplet<bool>(i, k, true));
                            previousNonZero = itI.row();
                            if (itI.value() == 0)
                                tripletList.push_back(Triplet<bool>(itI.row(), k, true));
                            ++itI;
                        }
                    } else if ((((itR) || (itI)) && (itRb)) && (row == rowb)) {
                        if ((!itI) || ((itR) && (itR.row() < itI.row()))) {
                            for (IndexType i(previousNonZero+1); i < itR.row(); ++i)
                                tripletList.push_back(Triplet<bool>(i, k, true));
                            previousNonZero = itR.row();
                            if (itR.value() == itRb.value())
                                tripletList.push_back(Triplet<bool>(itR.row(), k, true));
                            ++itRb;
                        } else if ((itR) && (itI) && (itR.row() == itI.row())) {
                            for (IndexType i(previousNonZero+1); i < itR.row(); ++i)
                                tripletList.push_back(Triplet<bool>(i, k, true));
                            previousNonZero = itR.row();
                            if ((itR.value() == itRb.value()) && (itI.value() == 0))
                                tripletList.push_back(Triplet<bool>(itR.row(), k, true));
                            ++itRb;
                            ++itR;
                            ++itI;
                        } else {
                            for (IndexType i(previousNonZero+1); i < itI.row(); ++i)
                                tripletList.push_back(Triplet<bool>(i, k, true));
                            previousNonZero = itI.row();
                            if ((0 == itRb.value()) && (itI.value() == 0))
                                tripletList.push_back(Triplet<bool>(itI.row(), k, true));
                            ++itRb;
                            ++itI;
                        }
                    } else {
                        for (IndexType i(previousNonZero+1); i < itRb.row(); ++i)
                            tripletList.push_back(Triplet<bool>(i, k, true));
                        previousNonZero = itRb.row();
                        if (0 == itRb.value())
                            tripletList.push_back(Triplet<bool>(itRb.row(), k, true));
                        ++itRb;
                    }
                }

                // If there are still zeros
                for (IndexType i(previousNonZero+1); i < matrixR.rows(); ++i)
                    tripletList.push_back(Triplet<bool>(i, k, true));
            }
        } else {
            for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
                SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
                SparseMatrix<mpreal>::InnerIterator itRb(b.matrixR,k);
                IndexType previousNonZero(-1);

                while ((itR) || (itRb)) {
                    if ((((itR) && (itRb)) && (itR.row() < itRb.row())) || (!itRb)) {
                        for (IndexType i(previousNonZero+1); i < itR.row(); ++i)
                            tripletList.push_back(Triplet<bool>(i, k, true));
                        previousNonZero = itR.row();
                        if (itR.value() == 0)
                            tripletList.push_back(Triplet<bool>(itR.row(), k, true));
                        ++itR;
                    } else if (((itR) && (itRb)) && (itR.row() == itRb.row())) {
                        for (IndexType i(previousNonZero+1); i < itR.row(); ++i)
                            tripletList.push_back(Triplet<bool>(i, k, true));
                        previousNonZero = itR.row();
                        if (itR.value() == itRb.value())
                            tripletList.push_back(Triplet<bool>(itR.row(), k, true));
                        ++itRb;
                    } else {
                        for (IndexType i(previousNonZero+1); i < itRb.row(); ++i)
                            tripletList.push_back(Triplet<bool>(i, k, true));
                        previousNonZero = itRb.row();
                        if (0 == itRb.value())
                            tripletList.push_back(Triplet<bool>(itRb.row(), k, true));
                        ++itRb;
                    }
                }

                // If there are still zeros
                for (IndexType i(previousNonZero+1); i < matrixR.rows(); ++i)
                    tripletList.push_back(Triplet<bool>(i, k, true));
            }
        }

        // Now we can assign the data to the sparse matrix
        result.setFromTriplets(tripletList.begin(), tripletList.end());
    }

    return result;
}

/* Test of non-equality  c = (a != b) */
/* NOTE : Here, we assume a lot about b when it is a scalar. This could be used
     to simplify the code in this case. */
SparseMatrix <bool> SparseGmpEigenMatrix::ne(const SparseGmpEigenMatrix& b) const
{
    SparseMatrix <bool> result;

    // We want to support comparison between an nxm matrix and a 1x1 matrix
    if ((numel() != 1) && (b.numel() == 1)) {
        // setting the output size
        result.resize(matrixR.rows(), matrixR.cols());

        // We will first copy the data to triplets
        vector< Triplet<bool> > tripletList;

        if (b.isComplex) {
            if (isComplex) {
                tripletList.reserve(matrixR.nonZeros() + matrixI.nonZeros()); // Here, b is assumed to be zero

                for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
                    SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
                    SparseMatrix<mpreal>::InnerIterator itI(matrixI,k);

                    while ((itR) || (itI)) {
                        if ((!itI) || ((itR) && (itR.row() < itI.row()))) {
                            if ((itR.value() != b.matrixR.coeff(0,0)) || (0 != b.matrixI.coeff(0,0)))
                                tripletList.push_back(Triplet<bool>(itR.row(), k, true));
                            ++itR;
                        } else if ((itR) && (itI) && (itR.row() == itI.row())) {
                            if ((itR.value() != b.matrixR.coeff(0,0)) || (itI.value() != b.matrixI.coeff(0,0)))
                                tripletList.push_back(Triplet<bool>(itR.row(), k, true));
                            ++itR;
                            ++itI;
                        } else {
                            if ((0 != b.matrixR.coeff(0,0)) || (itI.value() != b.matrixI.coeff(0,0)))
                                tripletList.push_back(Triplet<bool>(itI.row(), k, true));
                            ++itI;
                        }
                    }
                }
            } else {
                tripletList.reserve(matrixR.nonZeros()); // Here, b is assumed to be zero

                for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
                    SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);

                    while (itR) {
                        if ((itR.value() != b.matrixR.coeff(0,0)) || (0 != b.matrixI.coeff(0,0)))
                            tripletList.push_back(Triplet<bool>(itR.row(), k, true));
                        ++itR;
                    }
                }
            }
        } else {
            if (isComplex) {
                tripletList.reserve(matrixR.nonZeros()); // WARNING : Here, we assume that the real part of the scalar b is 0 (otherwise it will make the sparse matrix full...)

                for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
                    SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
                    SparseMatrix<mpreal>::InnerIterator itI(matrixI,k);

                    while ((itR) || (itI)) {
                        if ((!itI) || ((itR) && (itR.row() < itI.row()))) {
                            if (itR.value() != b.matrixR.coeff(0,0))
                                tripletList.push_back(Triplet<bool>(itR.row(), k, true));
                            ++itR;
                        } else if ((itR) && (itI) && (itR.row() == itI.row())) {
                            if ((itR.value() != b.matrixR.coeff(0,0)) || (itI.value() != 0))
                                tripletList.push_back(Triplet<bool>(itR.row(), k, true));
                            ++itR;
                            ++itI;
                        } else {
                            if ((0 != b.matrixR.coeff(0,0)) || (itI.value() != 0))
                                tripletList.push_back(Triplet<bool>(itI.row(), k, true));
                            ++itI;
                        }
                    }
                }
            } else {
                tripletList.reserve(matrixR.nonZeros()); // WARNING : Here, we assume that the real part of the scalar b is 0 (otherwise it will make the sparse matrix full...)

                for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
                    SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);

                    while (itR) {
                        if (itR.value() != b.matrixR.coeff(0,0))
                            tripletList.push_back(Triplet<bool>(itR.row(), k, true));
                        ++itR;
                    }
                }
            }
        }

        // Now we can assign the data to the sparse matrix
        result.setFromTriplets(tripletList.begin(), tripletList.end());
    } else if ((numel() == 1) && (b.numel() != 1)) {
        result = b.ne(*this);
    } else {
        // setting the output size
        result.resize(matrixR.rows(), matrixR.cols());

        // We will first copy the data to triplets
        vector< Triplet<bool> > tripletList;

        if (b.isComplex) {
            if (isComplex) {
                tripletList.reserve(matrixR.nonZeros() + matrixI.nonZeros() + b.matrixR.nonZeros() + b.matrixI.nonZeros());

                for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
                    SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
                    SparseMatrix<mpreal>::InnerIterator itI(matrixI,k);
                    SparseMatrix<mpreal>::InnerIterator itRb(b.matrixR,k);
                    SparseMatrix<mpreal>::InnerIterator itIb(b.matrixI,k);

                    while ((itR) || (itI) || (itRb) || (itIb)) {
                        IndexType row, rowb;
                        if ((itR) || (itI)) {
                            if ((itR) && (itI))
                                row = min(itR.row(), itI.row());
                            else if (itR)
                                row = itR.row();
                            else
                                row = itI.row();
                        }
                        if ((itRb) || (itIb)) {
                            if ((itRb) && (itIb))
                                rowb = min(itRb.row(), itIb.row());
                            else if (itRb)
                                rowb = itRb.row();
                            else
                                rowb = itIb.row();
                        }

                        if (((((itR) || (itI)) && ((itRb) || (itIb))) && (row < rowb)) || ((!itRb) && (!itIb))) {
                            if ((!itI) || ((itR) && (itR.row() < itI.row()))) {
                                if (itR.value() != 0)
                                    tripletList.push_back(Triplet<bool>(itR.row(), k, true));
                                ++itR;
                            } else if ((itR) && (itI) && (itR.row() == itI.row())) {
                                if ((itR.value() != 0) || (itI.value() != 0))
                                    tripletList.push_back(Triplet<bool>(itR.row(), k, true));
                                ++itR;
                                ++itI;
                            } else {
                                if (itI.value() != 0)
                                    tripletList.push_back(Triplet<bool>(itI.row(), k, true));
                                ++itI;
                            }
                        } else if ((((itR) || (itI)) && ((itRb) || (itIb))) && (row == rowb)) {
                            if ((!itI) || ((itR) && (itR.row() < itI.row()))) {
                                if ((!itIb) || ((itRb) && (itRb.row() < itIb.row()))) {
                                    if (itR.value() != itRb.value())
                                        tripletList.push_back(Triplet<bool>(itR.row(), k, true));
                                    ++itRb;
                                } else if ((itRb) && (itIb) && (itRb.row() == itIb.row())) {
                                    if ((itR.value() != itRb.value()) || (0 != itIb.value()))
                                        tripletList.push_back(Triplet<bool>(itR.row(), k, true));
                                    ++itRb;
                                    ++itIb;
                                } else {
                                    if (0 != itIb.value())
                                        tripletList.push_back(Triplet<bool>(itR.row(), k, true));
                                    ++itIb;
                                }
                                ++itR;
                            } else if ((itR) && (itI) && (itR.row() == itI.row())) {
                                if ((!itIb) || ((itRb) && (itRb.row() < itIb.row()))) {
                                    if ((itR.value() != itRb.value()) || (itI.value() != 0))
                                        tripletList.push_back(Triplet<bool>(itR.row(), k, true));
                                    ++itRb;
                                } else if ((itRb) && (itIb) && (itRb.row() == itIb.row())) {
                                    if ((itR.value() != itRb.value()) || (itI.value() != itIb.value()))
                                        tripletList.push_back(Triplet<bool>(itR.row(), k, true));
                                    ++itRb;
                                    ++itIb;
                                } else {
                                    if ((itR.value() != 0) || (itI.value() != itIb.value()))
                                        tripletList.push_back(Triplet<bool>(itR.row(), k, true));
                                    ++itIb;
                                }
                                ++itR;
                                ++itI;
                            } else {
                                if ((!itIb) || ((itRb) && (itRb.row() < itIb.row()))) {
                                    if ((0 != itRb.value()) || (itI.value() != 0))
                                        tripletList.push_back(Triplet<bool>(itI.row(), k, true));
                                    ++itRb;
                                } else if ((itRb) && (itIb) && (itRb.row() == itIb.row())) {
                                    if ((0 != itRb.value()) || (itI.value() != itIb.value()))
                                        tripletList.push_back(Triplet<bool>(itI.row(), k, true));
                                    ++itRb;
                                    ++itIb;
                                } else {
                                    if (itI.value() != itIb.value())
                                        tripletList.push_back(Triplet<bool>(itI.row(), k, true));
                                    ++itIb;
                                }
                                ++itI;
                            }
                        } else {
                            if ((!itIb) || ((itRb) && (itRb.row() < itIb.row()))) {
                                if (0 != itRb.value())
                                    tripletList.push_back(Triplet<bool>(itRb.row(), k, true));
                                ++itRb;
                            } else if ((itRb) && (itIb) && (itRb.row() == itIb.row())) {
                                if ((0 != itRb.value()) || (0 != itIb.value()))
                                    tripletList.push_back(Triplet<bool>(itRb.row(), k, true));
                                ++itRb;
                                ++itIb;
                            } else {
                                if (0 != itIb.value())
                                    tripletList.push_back(Triplet<bool>(itIb.row(), k, true));
                                ++itIb;
                            }
                        }
                    }
                }
            } else {
                tripletList.reserve(matrixR.nonZeros() + b.matrixR.nonZeros() + b.matrixI.nonZeros());

                for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
                    SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
                    SparseMatrix<mpreal>::InnerIterator itRb(b.matrixR,k);
                    SparseMatrix<mpreal>::InnerIterator itIb(b.matrixI,k);

                    while ((itR) || (itRb) || (itIb)) {
                        IndexType row, rowb;
                        if (itR)
                            row = itR.row();
                        if ((itRb) || (itIb)) {
                            if ((itRb) && (itIb))
                                rowb = min(itRb.row(), itIb.row());
                            else if (itRb)
                                rowb = itRb.row();
                            else
                                rowb = itIb.row();
                        }

                        if ((((itR) && ((itRb) || (itIb))) && (row < rowb)) || ((!itRb) && (!itIb))) {
                            if (itR.value() != 0)
                                tripletList.push_back(Triplet<bool>(itR.row(), k, true));
                            ++itR;
                        } else if (((itR) && ((itRb) || (itIb))) && (row == rowb)) {
                            if ((!itIb) || ((itRb) && (itRb.row() < itIb.row()))) {
                                if (itR.value() != itRb.value())
                                    tripletList.push_back(Triplet<bool>(itR.row(), k, true));
                                ++itRb;
                            } else if ((itRb) && (itIb) && (itRb.row() == itIb.row())) {
                                if ((itR.value() != itRb.value()) || (0 != itIb.value()))
                                    tripletList.push_back(Triplet<bool>(itR.row(), k, true));
                                ++itRb;
                                ++itIb;
                            } else {
                                if (0 != itIb.value())
                                    tripletList.push_back(Triplet<bool>(itR.row(), k, true));
                                ++itIb;
                            }
                            ++itR;
                        } else {
                            if ((!itIb) || ((itRb) && (itRb.row() < itIb.row()))) {
                                if (0 != itRb.value())
                                    tripletList.push_back(Triplet<bool>(itRb.row(), k, true));
                                ++itRb;
                            } else if ((itRb) && (itIb) && (itRb.row() == itIb.row())) {
                                if ((0 != itRb.value()) || (0 != itIb.value()))
                                    tripletList.push_back(Triplet<bool>(itRb.row(), k, true));
                                ++itRb;
                                ++itIb;
                            } else {
                                if (0 != itIb.value())
                                    tripletList.push_back(Triplet<bool>(itIb.row(), k, true));
                                ++itIb;
                            }
                        }
                    }
                }
            }
        } else if (isComplex) {
            tripletList.reserve(matrixR.nonZeros() + matrixI.nonZeros() + b.matrixR.nonZeros());
            for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
                SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
                SparseMatrix<mpreal>::InnerIterator itI(matrixI,k);
                SparseMatrix<mpreal>::InnerIterator itRb(b.matrixR,k);

                while ((itR) || (itI) || (itRb)) {
                    IndexType row, rowb;
                    if ((itR) || (itI)) {
                        if ((itR) && (itI))
                            row = min(itR.row(), itI.row());
                        else if (itR)
                            row = itR.row();
                        else
                            row = itI.row();
                    }
                    if (itRb)
                        rowb = itRb.row();

                    if (((((itR) || (itI)) && (itRb)) && (row < rowb)) || (!itRb)) {
                        if ((!itI) || ((itR) && (itR.row() < itI.row()))) {
                            if (itR.value() != 0)
                                tripletList.push_back(Triplet<bool>(itR.row(), k, true));
                            ++itR;
                        } else if ((itR) && (itI) && (itR.row() == itI.row())) {
                            if ((itR.value() != 0) || (itI.value() != 0))
                                tripletList.push_back(Triplet<bool>(itR.row(), k, true));
                            ++itR;
                            ++itI;
                        } else {
                            if (itI.value() != 0)
                                tripletList.push_back(Triplet<bool>(itI.row(), k, true));
                            ++itI;
                        }
                    } else if ((((itR) || (itI)) && (itRb)) && (row == rowb)) {
                        if ((!itI) || ((itR) && (itR.row() < itI.row()))) {
                            if (itR.value() != itRb.value())
                                tripletList.push_back(Triplet<bool>(itR.row(), k, true));
                            ++itRb;
                        } else if ((itR) && (itI) && (itR.row() == itI.row())) {
                            if ((itR.value() != itRb.value()) || (itI.value() != 0))
                                tripletList.push_back(Triplet<bool>(itR.row(), k, true));
                            ++itRb;
                            ++itR;
                            ++itI;
                        } else {
                            if ((0 != itRb.value()) || (itI.value() != 0))
                                tripletList.push_back(Triplet<bool>(itI.row(), k, true));
                            ++itRb;
                            ++itI;
                        }
                    } else {
                        if (0 != itRb.value())
                            tripletList.push_back(Triplet<bool>(itRb.row(), k, true));
                        ++itRb;
                    }
                }
            }
        } else {
            tripletList.reserve(matrixR.nonZeros() + b.matrixR.nonZeros());
            for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
                SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
                SparseMatrix<mpreal>::InnerIterator itRb(b.matrixR,k);

                while ((itR) || (itRb)) {
                    if ((((itR) && (itRb)) && (itR.row() < itRb.row())) || (!itRb)) {
                        if (itR.value() != 0)
                            tripletList.push_back(Triplet<bool>(itR.row(), k, true));
                        ++itR;
                    } else if (((itR) && (itRb)) && (itR.row() == itRb.row())) {
                        if (itR.value() != itRb.value())
                            tripletList.push_back(Triplet<bool>(itR.row(), k, true));
                        ++itRb;
                    } else {
                        if (0 != itRb.value())
                            tripletList.push_back(Triplet<bool>(itRb.row(), k, true));
                        ++itRb;
                    }
                }
            }
        }

        // Now we can assign the data to the sparse matrix
        result.setFromTriplets(tripletList.begin(), tripletList.end());
    }

    return result;
}

/* Test whether the value is nan: b = isnan(a) */
SparseMatrix <bool> SparseGmpEigenMatrix::isnan() const
{
    SparseMatrix <bool> result(matrixR.rows(), matrixR.cols());

    // We will first copy the data to triplets
    vector< Triplet<bool> > tripletList;

    if (isComplex) {
        tripletList.reserve(matrixR.nonZeros() + matrixI.nonZeros());
        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
            SparseMatrix<mpreal>::InnerIterator itI(matrixI,k);

            while ((itR) || (itI)) {
                if ((itR) && (itI)) {
                    if (itR.row() < itI.row()) {
                        if (mpfr::isnan(itR.value()))
                            tripletList.push_back(Triplet<bool>(itR.row(), k, true));
                        ++itR;
                    } else if (itR.row() == itI.row()) {
                        if ((mpfr::isnan(itR.value())) || (mpfr::isnan(itI.value())))
                            tripletList.push_back(Triplet<bool>(itR.row(), k, true));
                        ++itR;
                        ++itI;
                    } else {
                        if (mpfr::isnan(itI.value()))
                            tripletList.push_back(Triplet<bool>(itI.row(), k, true));
                        ++itI;
                    }
                } else if (itR) {
                    if (mpfr::isnan(itR.value()))
                        tripletList.push_back(Triplet<bool>(itR.row(), k, true));
                    ++itR;
                } else {
                    if (mpfr::isnan(itI.value()))
                        tripletList.push_back(Triplet<bool>(itI.row(), k, true));
                    ++itI;
                }
            }                
        }
    } else {
        tripletList.reserve(matrixR.nonZeros());
        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);

            while (itR) {
                if (mpfr::isnan(itR.value()))
                    tripletList.push_back(Triplet<bool>(itR.row(), k, true));
                ++itR;
            }
        }
    }

    // Now we can assign the data to the sparse matrix
    result.setFromTriplets(tripletList.begin(), tripletList.end());

    return result;
}
/* Test whether the value is +inf or -inf: b = isinf(a) */
SparseMatrix <bool> SparseGmpEigenMatrix::isinf() const
{
    SparseMatrix <bool> result(matrixR.rows(), matrixR.cols());

    // We will first copy the data to triplets
    vector< Triplet<bool> > tripletList;

    if (isComplex) {
        tripletList.reserve(matrixR.nonZeros() + matrixI.nonZeros());
        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
            SparseMatrix<mpreal>::InnerIterator itI(matrixI,k);

            while ((itR) || (itI)) {
                if ((itR) && (itI)) {
                    if (itR.row() < itI.row()) {
                        if (mpfr::isinf(itR.value()))
                            tripletList.push_back(Triplet<bool>(itR.row(), k, true));
                        ++itR;
                    } else if (itR.row() == itI.row()) {
                        if ((mpfr::isinf(itR.value())) || (mpfr::isinf(itI.value())))
                            tripletList.push_back(Triplet<bool>(itR.row(), k, true));
                        ++itR;
                        ++itI;
                    } else {
                        if (mpfr::isinf(itI.value()))
                            tripletList.push_back(Triplet<bool>(itI.row(), k, true));
                        ++itI;
                    }
                } else if (itR) {
                    if (mpfr::isinf(itR.value()))
                        tripletList.push_back(Triplet<bool>(itR.row(), k, true));
                    ++itR;
                } else {
                    if (mpfr::isinf(itI.value()))
                        tripletList.push_back(Triplet<bool>(itI.row(), k, true));
                    ++itI;
                }
            }                
        }
    } else {
        tripletList.reserve(matrixR.nonZeros());
        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);

            while (itR) {
                if (mpfr::isinf(itR.value()))
                    tripletList.push_back(Triplet<bool>(itR.row(), k, true));
                ++itR;
            }
        }
    }

    // Now we can assign the data to the sparse matrix
    result.setFromTriplets(tripletList.begin(), tripletList.end());

    return result;
}

/* Test whether the numerical values match
   this function is called by c = isequal(a, b)
   at this stage, we already know that both tables have the same size, that both
   are either real or complex, and that all their precision match. */
bool SparseGmpEigenMatrix::identicalValues(const SparseGmpEigenMatrix& b) const
{
    if (isComplex) {
        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
            SparseMatrix<mpreal>::InnerIterator itI(matrixI,k);
            SparseMatrix<mpreal>::InnerIterator itRb(b.matrixR,k);
            SparseMatrix<mpreal>::InnerIterator itIb(b.matrixI,k);

            while ((itR) || (itI) || (itRb) || (itIb)) {
                IndexType row, rowb;
                if ((itR) || (itI)) {
                    if ((itR) && (itI))
                        row = min(itR.row(), itI.row());
                    else if (itR)
                        row = itR.row();
                    else
                        row = itI.row();
                }
                if ((itRb) || (itIb)) {
                    if ((itRb) && (itIb))
                        rowb = min(itRb.row(), itIb.row());
                    else if (itRb)
                        rowb = itRb.row();
                    else
                        rowb = itIb.row();
                }

                if (((((itR) || (itI)) && ((itRb) || (itIb))) && (row < rowb)) || ((!itRb) && (!itIb))) {
                    if ((!itI) || ((itR) && (itR.row() < itI.row()))) {
                        if (itR.value() != 0)
                            return false;
                        ++itR;
                    } else if ((itR) && (itI) && (itR.row() == itI.row())) {
                        if ((itR.value() != 0) || (itI.value() != 0))
                            return false;
                        ++itR;
                        ++itI;
                    } else {
                        if (itI.value() != 0)
                            return false;
                        ++itI;
                    }
                } else if ((((itR) || (itI)) && ((itRb) || (itIb))) && (row == rowb)) {
                    if ((!itI) || ((itR) && (itR.row() < itI.row()))) {
                        if ((!itIb) || ((itRb) && (itRb.row() < itIb.row()))) {
                            if (itR.value() != itRb.value())
                                return false;
                            ++itRb;
                        } else if ((itRb) && (itIb) && (itRb.row() == itIb.row())) {
                            if ((itR.value() != itRb.value()) || (0 != itIb.value()))
                                return false;
                            ++itRb;
                            ++itIb;
                        } else {
                            if (0 != itIb.value())
                                return false;
                            ++itIb;
                        }
                        ++itR;
                    } else if ((itR) && (itI) && (itR.row() == itI.row())) {
                        if ((!itIb) || ((itRb) && (itRb.row() < itIb.row()))) {
                            if ((itR.value() != itRb.value()) || (itI.value() != 0))
                                return false;
                            ++itRb;
                        } else if ((itRb) && (itIb) && (itRb.row() == itIb.row())) {
                            if ((itR.value() != itRb.value()) || (itI.value() != itIb.value()))
                                return false;
                            ++itRb;
                            ++itIb;
                        } else {
                            if ((itR.value() != 0) || (itI.value() != itIb.value()))
                                return false;
                            ++itIb;
                        }
                        ++itR;
                        ++itI;
                    } else {
                        if ((!itIb) || ((itRb) && (itRb.row() < itIb.row()))) {
                            if ((0 != itRb.value()) || (itI.value() != 0))
                                return false;
                            ++itRb;
                        } else if ((itRb) && (itIb) && (itRb.row() == itIb.row())) {
                            if ((0 != itRb.value()) || (itI.value() != itIb.value()))
                                return false;
                            ++itRb;
                            ++itIb;
                        } else {
                            if (itI.value() != itIb.value())
                                return false;
                            ++itIb;
                        }
                        ++itI;
                    }
                } else {
                    if ((!itIb) || ((itRb) && (itRb.row() < itIb.row()))) {
                        if (0 != itRb.value())
                            return false;
                        ++itRb;
                    } else if ((itRb) && (itIb) && (itRb.row() == itIb.row())) {
                        if ((0 != itRb.value()) || (0 != itIb.value()))
                            return false;
                        ++itRb;
                        ++itIb;
                    } else {
                        if (0 != itIb.value())
                            return false;
                        ++itIb;
                    }
                }
            }
        }
    } else {
        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
            SparseMatrix<mpreal>::InnerIterator itRb(b.matrixR,k);

            while ((itR) || (itRb)) {
                if ((((itR) && (itRb)) && (itR.row() < itRb.row())) || (!itRb)) {
                    if (itR.value() != 0)
                        return false;
                    ++itR;
                } else if (((itR) && (itRb)) && (itR.row() == itRb.row())) {
                    if (itR.value() != itRb.value())
                        return false;
                    ++itR;
                    ++itRb;
                } else {
                    if (0 != itRb.value())
                        return false;
                    ++itRb;
                }
            }
        }
    }

    return true;
}

/* Test whether the numerical values match
   this function is called by c = isequal(a, b)
   This function makes no difference between NaN numbers.
   at this stage, we already know that both tables have the same size, that both
   are either real or complex, and that all their precision match. */
bool SparseGmpEigenMatrix::identicalValuesNaNok(const SparseGmpEigenMatrix& b) const
{
    if (isComplex) {
        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
            SparseMatrix<mpreal>::InnerIterator itI(matrixI,k);
            SparseMatrix<mpreal>::InnerIterator itRb(b.matrixR,k);
            SparseMatrix<mpreal>::InnerIterator itIb(b.matrixI,k);

            while ((itR) || (itI) || (itRb) || (itIb)) {
                IndexType row, rowb;
                if ((itR) || (itI)) {
                    if ((itR) && (itI))
                        row = min(itR.row(), itI.row());
                    else if (itR)
                        row = itR.row();
                    else
                        row = itI.row();
                }
                if ((itRb) || (itIb)) {
                    if ((itRb) && (itIb))
                        rowb = min(itRb.row(), itIb.row());
                    else if (itRb)
                        rowb = itRb.row();
                    else
                        rowb = itIb.row();
                }

                if (((((itR) || (itI)) && ((itRb) || (itIb))) && (row < rowb)) || ((!itRb) && (!itIb))) {
                    if ((!itI) || ((itR) && (itR.row() < itI.row()))) {
                        if (itR.value() != 0)
                            return false;
                        ++itR;
                    } else if ((itR) && (itI) && (itR.row() == itI.row())) {
                        if ((itR.value() != 0) || (itI.value() != 0))
                            return false;
                        ++itR;
                        ++itI;
                    } else {
                        if (itI.value() != 0)
                            return false;
                        ++itI;
                    }
                } else if ((((itR) || (itI)) && ((itRb) || (itIb))) && (row == rowb)) {
                    if ((!itI) || ((itR) && (itR.row() < itI.row()))) {
                        if ((!itIb) || ((itRb) && (itRb.row() < itIb.row()))) {
                            if (!((itR.value() == itRb.value()) || (mpfr::isnan(itR.value()) && mpfr::isnan(itRb.value()))))
                                return false;
                            ++itRb;
                        } else if ((itRb) && (itIb) && (itRb.row() == itIb.row())) {
                            if ((!((itR.value() == itRb.value()) || (mpfr::isnan(itR.value()) && mpfr::isnan(itRb.value())))) || (0 != itIb.value()))
                                return false;
                            ++itRb;
                            ++itIb;
                        } else {
                            if (0 != itIb.value())
                                return false;
                            ++itIb;
                        }
                        ++itR;
                    } else if ((itR) && (itI) && (itR.row() == itI.row())) {
                        if ((!itIb) || ((itRb) && (itRb.row() < itIb.row()))) {
                            if ((!((itR.value() == itRb.value()) || (mpfr::isnan(itR.value()) && mpfr::isnan(itRb.value())))) || (itI.value() != 0))
                                return false;
                            ++itRb;
                        } else if ((itRb) && (itIb) && (itRb.row() == itIb.row())) {
                            if ((!((itR.value() == itRb.value()) || (mpfr::isnan(itR.value()) && mpfr::isnan(itRb.value())))) || (!((itI.value() == itIb.value()) || (mpfr::isnan(itI.value()) && mpfr::isnan(itIb.value())))))
                                return false;
                            ++itRb;
                            ++itIb;
                        } else {
                            if ((itR.value() != 0) || (!((itI.value() == itIb.value()) || (mpfr::isnan(itI.value()) && mpfr::isnan(itIb.value())))))
                                return false;
                            ++itIb;
                        }
                        ++itR;
                        ++itI;
                    } else {
                        if ((!itIb) || ((itRb) && (itRb.row() < itIb.row()))) {
                            if ((0 != itRb.value()) || (itI.value() != 0))
                                return false;
                            ++itRb;
                        } else if ((itRb) && (itIb) && (itRb.row() == itIb.row())) {
                            if ((0 != itRb.value()) || (!((itI.value() == itIb.value()) || (mpfr::isnan(itI.value()) && mpfr::isnan(itIb.value())))))
                                return false;
                            ++itRb;
                            ++itIb;
                        } else {
                            if (!((itI.value() == itIb.value()) || (mpfr::isnan(itI.value()) && mpfr::isnan(itIb.value()))))
                                return false;
                            ++itIb;
                        }
                        ++itI;
                    }
                } else {
                    if ((!itIb) || ((itRb) && (itRb.row() < itIb.row()))) {
                        if (0 != itRb.value())
                            return false;
                        ++itRb;
                    } else if ((itRb) && (itIb) && (itRb.row() == itIb.row())) {
                        if ((0 != itRb.value()) || (0 != itIb.value()))
                            return false;
                        ++itRb;
                        ++itIb;
                    } else {
                        if (0 != itIb.value())
                            return false;
                        ++itIb;
                    }
                }
            }
        }
    } else {
        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
            SparseMatrix<mpreal>::InnerIterator itRb(b.matrixR,k);

            while ((itR) || (itRb)) {
                if ((((itR) && (itRb)) && (itR.row() < itRb.row())) || (!itRb)) {
                    if (itR.value() != 0)
                        return false;
                    ++itR;
                } else if (((itR) && (itRb)) && (itR.row() == itRb.row())) {
                    if (!((itR.value() == itRb.value()) || (mpfr::isnan(itR.value()) && mpfr::isnan(itRb.value()))))
                        return false;
                    ++itR;
                    ++itRb;
                } else {
                    if (0 != itRb.value())
                        return false;
                    ++itRb;
                }
            }
        }
    }

    return true;
}





// symmetry tests
bool SparseGmpEigenMatrix::issymmetric() const
{
    if (isComplex) {
        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
            SparseMatrix<mpreal>::InnerIterator itI(matrixI,k);

            while ((itR) || (itI)) {
                if ((!itI) || ((itR) && (itR.row() < itI.row()))) {
                    if (itR.value() != matrixR.coeff(itR.col(), itR.row()))
                        return false;
                    ++itR;
                } else if ((itR) && (itI) && (itR.row() == itI.row())) {
                    if ((itR.value() != matrixR.coeff(itR.col(), itR.row())) || (itI.value() != matrixI.coeff(itR.col(), itR.row())))
                        return false;
                    ++itR;
                    ++itI;
                } else {
                    if (itI.value() != matrixI.coeff(itI.col(), itI.row()))
                        return false;
                    ++itI;
                }
            }
        }
    } else {
        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);

            while (itR) {
                if (itR.value() != matrixR.coeff(itR.col(), itR.row()))
                    return false;
                ++itR;
            }
        }
    }

    return true;
}

bool SparseGmpEigenMatrix::ishermitian() const
{
    if (isComplex) {
        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
            SparseMatrix<mpreal>::InnerIterator itI(matrixI,k);

            while ((itR) || (itI)) {
                if ((!itI) || ((itR) && (itR.row() < itI.row()))) {
                    if (itR.value() != matrixR.coeff(itR.col(), itR.row()))
                        return false;
                    ++itR;
                } else if ((itR) && (itI) && (itR.row() == itI.row())) {
                    if ((itR.value() != matrixR.coeff(itR.col(), itR.row())) || (itI.value() != -matrixI.coeff(itR.col(), itR.row())))
                        return false;
                    ++itR;
                    ++itI;
                } else {
                    if (itI.value() != -matrixI.coeff(itI.col(), itI.row()))
                        return false;
                    ++itI;
                }
            }
        }
    } else {
        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);

            while (itR) {
                if (itR.value() != matrixR.coeff(itR.col(), itR.row()))
                    return false;
                ++itR;
            }
        }
    }

    return true;
}




// Note : Eigen's minCoeff function behaves differently than matlab's min
// function (e.g. if gives no result in presence of NaN instead of ignoring them,
// so we write our own implementation)

// column-wise minimum b = min(a)
SparseGmpEigenMatrix SparseGmpEigenMatrix::colMin(vector<IndexType>& indices) const
{
    SparseGmpEigenMatrix result;

    // We initialize the index vector
    indices.clear();
    indices.resize(matrixR.cols(), 0);

    if (isComplex) {
        // We initialize the result matrix
        result.isComplex = true;
        result.matrixR.resize(1,matrixR.cols());
        result.matrixI.resize(1,matrixI.cols());
        result.matrixR.reserve(min(matrixR.rows(),matrixR.nonZeros()));
        result.matrixI.reserve(min(matrixI.rows(),matrixI.nonZeros()));

        // We find the minimum
        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            mpreal minValue;
            mpreal minAngle;
            minValue.setInf(+1);
            minAngle.setInf(+1);
            IndexType firstZero(0); // This variable is incremented until there is a gap in column enumeration

            SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
            SparseMatrix<mpreal>::InnerIterator itI(matrixI,k);
            while ((itR) || (itI)) {
                if ((itR) && (itI)) {
                    if (itR.row() < itI.row()) {
                        if (firstZero == itR.row())
                            ++firstZero;
                        if (!mpfr::isnan(itR.value())) {
                            if ((pow(itR.value(),2) < minValue) || ((pow(itR.value(),2) == minValue) && (atan2(0, itR.value()) < minAngle))) {
                                indices[k] = itR.row();
                                minValue = pow(itR.value(),2);
                                minAngle = atan2(0, itR.value());
                            }
                        }
                        ++itR;
                    } else if (itR.row() == itI.row()) {
                        if (firstZero == itR.row())
                            ++firstZero;
                        if ((!mpfr::isnan(itR.value())) && (!mpfr::isnan(itI.value()))) {
                            if (((pow(itR.value(),2) + pow(itI.value(),2)) < minValue) || (((pow(itR.value(),2) + pow(itI.value(),2)) == minValue) && (atan2(itI.value(), itR.value()) < minAngle))) {
                                indices[k] = itR.row();
                                minValue = pow(itR.value(),2) + pow(itI.value(),2);
                                minAngle = atan2(itI.value(), itR.value());
                            }
                        }
                        ++itR;
                        ++itI;
                    } else {
                        if (firstZero == itI.row())
                            ++firstZero;
                        if (!mpfr::isnan(itI.value())) {
                            if ((pow(itI.value(),2) < minValue) || ((pow(itI.value(),2) == minValue) && (atan2(itI.value(), 0) < minAngle))) {
                                indices[k] = itI.row();
                                minValue = pow(itI.value(),2);
                                minAngle = atan2(itI.value(), 0);
                            }
                        }
                        ++itI;
                    }
                } else if (itR) {
                    if (firstZero == itR.row())
                        ++firstZero;
                    if (!mpfr::isnan(itR.value())) {
                        if ((pow(itR.value(),2) < minValue) || ((pow(itR.value(),2) == minValue) && (atan2(0, itR.value()) < minAngle))) {
                            indices[k] = itR.row();
                            minValue = pow(itR.value(),2);
                            minAngle = atan2(0, itR.value());
                        }
                    }
                    ++itR;
                } else {
                    if (firstZero == itI.row())
                        ++firstZero;
                    if (!mpfr::isnan(itI.value())) {
                        if ((pow(itI.value(),2) < minValue) || ((pow(itI.value(),2) == minValue) && (atan2(itI.value(), 0) < minAngle))) {
                            indices[k] = itI.row();
                            minValue = pow(itI.value(),2);
                            minAngle = atan2(itI.value(), 0);
                        }
                    }
                    ++itI;
                }
            }

            if ((firstZero < matrixR.rows()) && (0 < minValue)) {// If the column contains zeros
                indices[k] = firstZero;
                minValue = 0;
            }
            if (minValue != 0) {
                if (matrixR.coeff(indices[k],k) != 0)
                    result.matrixR.insert(0,k) = matrixR.coeff(indices[k],k);
                if (matrixI.coeff(indices[k],k) != 0)
                    result.matrixI.insert(0,k) = matrixI.coeff(indices[k],k);
            }
        }
        result.matrixR.makeCompressed();
        result.matrixI.makeCompressed();

        result.checkComplexity();
    } else {
        // We initialize the result matrix
        result.isComplex = false;
        result.matrixR.resize(1,matrixR.cols());
        result.matrixR.reserve(min(matrixR.rows(), matrixR.nonZeros()));

        // We find the minimum
        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            mpreal minValue;
            minValue.setInf(+1);
            IndexType firstZero(0); // This variable is incremented until there is a gap in column enumeration
            for (SparseMatrix<mpreal>::InnerIterator it(matrixR,k); it; ++it) {
                if (firstZero == it.row())
                    ++firstZero;
                if ((!mpfr::isnan(it.value())) && (it.value() < minValue)) {
                    indices[k] = it.row();
                    minValue = it.value();
                }
            }

            if ((firstZero < matrixR.rows()) && (0 < minValue)) {// If the column contains zeros
                indices[k] = firstZero;
                minValue = 0;
            }
            if (minValue != 0)
                result.matrixR.insert(0,k) = minValue;
        }
        result.matrixR.makeCompressed();
    }

    return result;
}

// column-wise minimum b = min(a)
SparseGmpEigenMatrix& SparseGmpEigenMatrix::colMin_new(vector<IndexType>& indices) const
{
    SparseGmpEigenMatrix& result(*(new SparseGmpEigenMatrix));

    // We initialize the index vector
    indices.clear();
    indices.resize(matrixR.cols(), 0);

    if (isComplex) {
        // We initialize the result matrix
        result.isComplex = true;
        result.matrixR.resize(1,matrixR.cols());
        result.matrixI.resize(1,matrixI.cols());
        result.matrixR.reserve(min(matrixR.rows(),matrixR.nonZeros()));
        result.matrixI.reserve(min(matrixI.rows(),matrixI.nonZeros()));

        // We find the minimum
        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            mpreal minValue;
            mpreal minAngle;
            minValue.setInf(+1);
            minAngle.setInf(+1);
            IndexType firstZero(0); // This variable is incremented until there is a gap in column enumeration

            SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
            SparseMatrix<mpreal>::InnerIterator itI(matrixI,k);
            while ((itR) || (itI)) {
                if ((itR) && (itI)) {
                    if (itR.row() < itI.row()) {
                        if (firstZero == itR.row())
                            ++firstZero;
                        if (!mpfr::isnan(itR.value())) {
                            if ((pow(itR.value(),2) < minValue) || ((pow(itR.value(),2) == minValue) && (atan2(0, itR.value()) < minAngle))) {
                                indices[k] = itR.row();
                                minValue = pow(itR.value(),2);
                                minAngle = atan2(0, itR.value());
                            }
                        }
                        ++itR;
                    } else if (itR.row() == itI.row()) {
                        if (firstZero == itR.row())
                            ++firstZero;
                        if ((!mpfr::isnan(itR.value())) && (!mpfr::isnan(itI.value()))) {
                            if (((pow(itR.value(),2) + pow(itI.value(),2)) < minValue) || (((pow(itR.value(),2) + pow(itI.value(),2)) == minValue) && (atan2(itI.value(), itR.value()) < minAngle))) {
                                indices[k] = itR.row();
                                minValue = pow(itR.value(),2) + pow(itI.value(),2);
                                minAngle = atan2(itI.value(), itR.value());
                            }
                        }
                        ++itR;
                        ++itI;
                    } else {
                        if (firstZero == itI.row())
                            ++firstZero;
                        if (!mpfr::isnan(itI.value())) {
                            if ((pow(itI.value(),2) < minValue) || ((pow(itI.value(),2) == minValue) && (atan2(itI.value(), 0) < minAngle))) {
                                indices[k] = itI.row();
                                minValue = pow(itI.value(),2);
                                minAngle = atan2(itI.value(), 0);
                            }
                        }
                        ++itI;
                    }
                } else if (itR) {
                    if (firstZero == itR.row())
                        ++firstZero;
                    if (!mpfr::isnan(itR.value())) {
                        if ((pow(itR.value(),2) < minValue) || ((pow(itR.value(),2) == minValue) && (atan2(0, itR.value()) < minAngle))) {
                            indices[k] = itR.row();
                            minValue = pow(itR.value(),2);
                            minAngle = atan2(0, itR.value());
                        }
                    }
                    ++itR;
                } else {
                    if (firstZero == itI.row())
                        ++firstZero;
                    if (!mpfr::isnan(itI.value())) {
                        if ((pow(itI.value(),2) < minValue) || ((pow(itI.value(),2) == minValue) && (atan2(itI.value(), 0) < minAngle))) {
                            indices[k] = itI.row();
                            minValue = pow(itI.value(),2);
                            minAngle = atan2(itI.value(), 0);
                        }
                    }
                    ++itI;
                }
            }

            if ((firstZero < matrixR.rows()) && (0 < minValue)) {// If the column contains zeros
                indices[k] = firstZero;
                minValue = 0;
            }
            if (minValue != 0) {
                if (matrixR.coeff(indices[k],k) != 0)
                    result.matrixR.insert(0,k) = matrixR.coeff(indices[k],k);
                if (matrixI.coeff(indices[k],k) != 0)
                    result.matrixI.insert(0,k) = matrixI.coeff(indices[k],k);
            }
        }
        result.matrixR.makeCompressed();
        result.matrixI.makeCompressed();

        result.checkComplexity();
    } else {
        // We initialize the result matrix
        result.isComplex = false;
        result.matrixR.resize(1,matrixR.cols());
        result.matrixR.reserve(min(matrixR.rows(), matrixR.nonZeros()));

        // We find the minimum
        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            mpreal minValue;
            minValue.setInf(+1);
            IndexType firstZero(0); // This variable is incremented until there is a gap in column enumeration
            for (SparseMatrix<mpreal>::InnerIterator it(matrixR,k); it; ++it) {
                if (firstZero == it.row())
                    ++firstZero;
                if ((!mpfr::isnan(it.value())) && (it.value() < minValue)) {
                    indices[k] = it.row();
                    minValue = it.value();
                }
            }

            if ((firstZero < matrixR.rows()) && (0 < minValue)) {// If the column contains zeros
                indices[k] = firstZero;
                minValue = 0;
            }
            if (minValue != 0)
                result.matrixR.insert(0,k) = minValue;
        }
        result.matrixR.makeCompressed();
    }

    return result;
}

// line-wise minimum b = min(a,[],2)
SparseGmpEigenMatrix SparseGmpEigenMatrix::rowMin(vector<IndexType>& indices) const
{
    SparseGmpEigenMatrix result;

    // We initialize the index vector
    indices.clear();
    indices.resize(matrixR.rows(), 0);

    if (isComplex) {
        // We initialize the result matrix
        result.isComplex = true;
        result.matrixR.resize(matrixR.rows(),1);
        result.matrixI.resize(matrixI.rows(),1);
        result.matrixR.reserve(min(matrixR.cols(),matrixR.nonZeros()));
        result.matrixI.reserve(min(matrixI.cols(),matrixI.nonZeros()));
        vector < mpreal > minValue(matrixR.rows(),const_infinity(+1));
        vector < mpreal > minAngle(matrixR.rows(),const_infinity(+1));

        // We find the minimums
        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            IndexType previousNonZero(-1);

            SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
            SparseMatrix<mpreal>::InnerIterator itI(matrixI,k);

            // If there is nothing in this column
            if ((!itR) && (!itI)) {
                for (IndexType i(0); i < matrixR.rows(); ++i) {
                    if (0 < minValue[i]) {
                        indices[i] = k;
                        minValue[i] = 0;
                    }
                }
            }

            while ((itR) || (itI)) {
                if ((itR) && (itI)) {
                    if (itR.row() < itI.row()) {
                        for (IndexType i(previousNonZero+1); i < itR.row(); ++i) {
                            if (0 < minValue[i]) {
                                indices[i] = k;
                                minValue[i] = 0;
                            }
                        }
                        previousNonZero = itR.row();
                        if (!mpfr::isnan(itR.value())) {
                            if ((pow(itR.value(),2) < minValue[itR.row()]) || ((pow(itR.value(),2) == minValue[itR.row()]) && (atan2(0, itR.value()) < minAngle[itR.row()]))) {
                                indices[itR.row()] = k;
                                minValue[itR.row()] = pow(itR.value(),2);
                                minAngle[itR.row()] = atan2(0, itR.value());
                            }
                        }
                        ++itR;
                    } else if (itR.row() == itI.row()) {
                        for (IndexType i(previousNonZero+1); i < itR.row(); ++i) {
                            if (0 < minValue[i]) {
                                indices[i] = k;
                                minValue[i] = 0;
                            }
                        }
                        previousNonZero = itR.row();
                        if ((!mpfr::isnan(itR.value())) && (!mpfr::isnan(itI.value()))) {
                            if (((pow(itR.value(),2) + pow(itI.value(),2)) < minValue[itR.row()]) || (((pow(itR.value(),2) + pow(itI.value(),2)) == minValue[itR.row()]) && (atan2(itI.value(), itR.value()) < minAngle[itR.row()]))) {
                                indices[itR.row()] = k;
                                minValue[itR.row()] = pow(itR.value(),2) + pow(itI.value(),2);
                                minAngle[itR.row()] = atan2(itI.value(), itR.value());
                            }
                        }
                        ++itR;
                        ++itI;
                    } else {
                        for (IndexType i(previousNonZero+1); i < itI.row(); ++i) {
                            if (0 < minValue[i]) {
                                indices[i] = k;
                                minValue[i] = 0;
                            }
                        }
                        previousNonZero = itR.row();
                        if (!mpfr::isnan(itI.value())) {
                            if ((pow(itI.value(),2) < minValue[itI.row()]) || ((pow(itI.value(),2) == minValue[itI.row()]) && (atan2(itI.value(), 0) < minAngle[itI.row()]))) {
                                indices[itI.row()] = k;
                                minValue[itI.row()] = pow(itI.value(),2);
                                minAngle[itI.row()] = atan2(itI.value(), 0);
                            }
                        }
                        ++itI;
                    }
                } else if (itR) {
                    for (IndexType i(previousNonZero+1); i < itR.row(); ++i) {
                        if (0 < minValue[i]) {
                            indices[i] = k;
                            minValue[i] = 0;
                        }
                    }
                    previousNonZero = itR.row();
                    if (!mpfr::isnan(itR.value())) {
                        if ((pow(itR.value(),2) < minValue[itR.row()]) || ((pow(itR.value(),2) == minValue[itR.row()]) && (atan2(0, itR.value()) < minAngle[itR.row()]))) {
                            indices[itR.row()] = k;
                            minValue[itR.row()] = pow(itR.value(),2);
                            minAngle[itR.row()] = atan2(0, itR.value());
                        }
                    }
                    ++itR;
                } else {
                    for (IndexType i(previousNonZero+1); i < itI.row(); ++i) {
                        if (0 < minValue[i]) {
                            indices[i] = k;
                            minValue[i] = 0;
                        }
                    }
                    previousNonZero = itI.row();
                    if (!mpfr::isnan(itI.value())) {
                        if ((pow(itI.value(),2) < minValue[itI.row()]) || ((pow(itI.value(),2) == minValue[itI.row()]) && (atan2(itI.value(), 0) < minAngle[itI.row()]))) {
                            indices[itI.row()] = k;
                            minValue[itI.row()] = pow(itI.value(),2);
                            minAngle[itI.row()] = atan2(itI.value(), 0);
                        }
                    }
                    ++itI;
                }
            }

            // If there are still remaining zeros
            for (IndexType i(previousNonZero+1); i < matrixR.rows(); ++i) {
                if (0 < minValue[i]) {
                    indices[i] = k;
                    minValue[i] = 0;
                }
            }
        }

        // We copy the minimums
        for (IndexType i = 0; i < matrixR.rows(); ++i) {
            if (minValue[i] != 0) {
                if (matrixR.coeff(i,indices[i]) != 0)
                    result.matrixR.insert(i,0) = matrixR.coeff(i,indices[i]);
                if (matrixI.coeff(i,indices[i]) != 0)
                    result.matrixI.insert(i,0) = matrixI.coeff(i,indices[i]);
            }
        }
        result.matrixR.makeCompressed();
        result.matrixI.makeCompressed();

        result.checkComplexity();
    } else {
        // We initialize the result matrix
        result.isComplex = false;
        result.matrixR.resize(matrixR.rows(),1);
        result.matrixR.reserve(min(matrixR.cols(),matrixR.nonZeros()));
        vector < mpreal > minValue(matrixR.rows(),const_infinity(+1));

        // We find the minimums
        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            IndexType previousNonZero(-1);
            SparseMatrix<mpreal>::InnerIterator it(matrixR,k);

            // If there is nothing in this column
            if (!it) {
                for (IndexType i(0); i < matrixR.rows(); ++i) {
                    if (0 < minValue[i]) {
                        indices[i] = k;
                        minValue[i] = 0;
                    }
                }
            }

            for (; it; ++it) {
                for (IndexType i(previousNonZero+1); i < it.row(); ++i) {
                    if (0 < minValue[i]) {
                        indices[i] = k;
                        minValue[i] = 0;
                    }
                }
                previousNonZero = it.row();
                if ((!mpfr::isnan(it.value())) && (it.value() < minValue[it.row()])) {
                    indices[it.row()] = k;
                    minValue[it.row()] = it.value();
                }
            }

            // If there are still zeros
            for (IndexType i(previousNonZero+1); i < matrixR.rows(); ++i) {
                if (0 < minValue[i]) {
                    indices[i] = k;
                    minValue[i] = 0;
                }
            }
        }

        // We copy the minimums
        for (IndexType i = 0; i < matrixR.rows(); ++i) {
            if (minValue[i] != 0) {
                result.matrixR.insert(i,0) = minValue[i];
            }
        }
        result.matrixR.makeCompressed();
    }

    return result;
}

// line-wise minimum b = min(a,[],2)
SparseGmpEigenMatrix& SparseGmpEigenMatrix::rowMin_new(vector<IndexType>& indices) const
{
    SparseGmpEigenMatrix& result(*(new SparseGmpEigenMatrix));

    // We initialize the index vector
    indices.clear();
    indices.resize(matrixR.rows(), 0);

    if (isComplex) {
        // We initialize the result matrix
        result.isComplex = true;
        result.matrixR.resize(matrixR.rows(),1);
        result.matrixI.resize(matrixI.rows(),1);
        result.matrixR.reserve(min(matrixR.cols(),matrixR.nonZeros()));
        result.matrixI.reserve(min(matrixI.cols(),matrixI.nonZeros()));
        vector < mpreal > minValue(matrixR.rows(),const_infinity(+1));
        vector < mpreal > minAngle(matrixR.rows(),const_infinity(+1));

        // We find the minimums
        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            IndexType previousNonZero(-1);

            SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
            SparseMatrix<mpreal>::InnerIterator itI(matrixI,k);

            // If there is nothing in this column
            if ((!itR) && (!itI)) {
                for (IndexType i(0); i < matrixR.rows(); ++i) {
                    if (0 < minValue[i]) {
                        indices[i] = k;
                        minValue[i] = 0;
                    }
                }
            }

            while ((itR) || (itI)) {
                if ((itR) && (itI)) {
                    if (itR.row() < itI.row()) {
                        for (IndexType i(previousNonZero+1); i < itR.row(); ++i) {
                            if (0 < minValue[i]) {
                                indices[i] = k;
                                minValue[i] = 0;
                            }
                        }
                        previousNonZero = itR.row();
                        if (!mpfr::isnan(itR.value())) {
                            if ((pow(itR.value(),2) < minValue[itR.row()]) || ((pow(itR.value(),2) == minValue[itR.row()]) && (atan2(0, itR.value()) < minAngle[itR.row()]))) {
                                indices[itR.row()] = k;
                                minValue[itR.row()] = pow(itR.value(),2);
                                minAngle[itR.row()] = atan2(0, itR.value());
                            }
                        }
                        ++itR;
                    } else if (itR.row() == itI.row()) {
                        for (IndexType i(previousNonZero+1); i < itR.row(); ++i) {
                            if (0 < minValue[i]) {
                                indices[i] = k;
                                minValue[i] = 0;
                            }
                        }
                        previousNonZero = itR.row();
                        if ((!mpfr::isnan(itR.value())) && (!mpfr::isnan(itI.value()))) {
                            if (((pow(itR.value(),2) + pow(itI.value(),2)) < minValue[itR.row()]) || (((pow(itR.value(),2) + pow(itI.value(),2)) == minValue[itR.row()]) && (atan2(itI.value(), itR.value()) < minAngle[itR.row()]))) {
                                indices[itR.row()] = k;
                                minValue[itR.row()] = pow(itR.value(),2) + pow(itI.value(),2);
                                minAngle[itR.row()] = atan2(itI.value(), itR.value());
                            }
                        }
                        ++itR;
                        ++itI;
                    } else {
                        for (IndexType i(previousNonZero+1); i < itI.row(); ++i) {
                            if (0 < minValue[i]) {
                                indices[i] = k;
                                minValue[i] = 0;
                            }
                        }
                        previousNonZero = itR.row();
                        if (!mpfr::isnan(itI.value())) {
                            if ((pow(itI.value(),2) < minValue[itI.row()]) || ((pow(itI.value(),2) == minValue[itI.row()]) && (atan2(itI.value(), 0) < minAngle[itI.row()]))) {
                                indices[itI.row()] = k;
                                minValue[itI.row()] = pow(itI.value(),2);
                                minAngle[itI.row()] = atan2(itI.value(), 0);
                            }
                        }
                        ++itI;
                    }
                } else if (itR) {
                    for (IndexType i(previousNonZero+1); i < itR.row(); ++i) {
                        if (0 < minValue[i]) {
                            indices[i] = k;
                            minValue[i] = 0;
                        }
                    }
                    previousNonZero = itR.row();
                    if (!mpfr::isnan(itR.value())) {
                        if ((pow(itR.value(),2) < minValue[itR.row()]) || ((pow(itR.value(),2) == minValue[itR.row()]) && (atan2(0, itR.value()) < minAngle[itR.row()]))) {
                            indices[itR.row()] = k;
                            minValue[itR.row()] = pow(itR.value(),2);
                            minAngle[itR.row()] = atan2(0, itR.value());
                        }
                    }
                    ++itR;
                } else {
                    for (IndexType i(previousNonZero+1); i < itI.row(); ++i) {
                        if (0 < minValue[i]) {
                            indices[i] = k;
                            minValue[i] = 0;
                        }
                    }
                    previousNonZero = itI.row();
                    if (!mpfr::isnan(itI.value())) {
                        if ((pow(itI.value(),2) < minValue[itI.row()]) || ((pow(itI.value(),2) == minValue[itI.row()]) && (atan2(itI.value(), 0) < minAngle[itI.row()]))) {
                            indices[itI.row()] = k;
                            minValue[itI.row()] = pow(itI.value(),2);
                            minAngle[itI.row()] = atan2(itI.value(), 0);
                        }
                    }
                    ++itI;
                }
            }

            // If there are still remaining zeros
            for (IndexType i(previousNonZero+1); i < matrixR.rows(); ++i) {
                if (0 < minValue[i]) {
                    indices[i] = k;
                    minValue[i] = 0;
                }
            }
        }

        // We copy the minimums
        for (IndexType i = 0; i < matrixR.rows(); ++i) {
            if (minValue[i] != 0) {
                if (matrixR.coeff(i,indices[i]) != 0)
                    result.matrixR.insert(i,0) = matrixR.coeff(i,indices[i]);
                if (matrixI.coeff(i,indices[i]) != 0)
                    result.matrixI.insert(i,0) = matrixI.coeff(i,indices[i]);
            }
        }
        result.matrixR.makeCompressed();
        result.matrixI.makeCompressed();

        result.checkComplexity();
    } else {
        // We initialize the result matrix
        result.isComplex = false;
        result.matrixR.resize(matrixR.rows(),1);
        result.matrixR.reserve(min(matrixR.cols(),matrixR.nonZeros()));
        vector < mpreal > minValue(matrixR.rows(),const_infinity(+1));

        // We find the minimums
        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            IndexType previousNonZero(-1);
            SparseMatrix<mpreal>::InnerIterator it(matrixR,k);

            // If there is nothing in this column
            if (!it) {
                for (IndexType i(0); i < matrixR.rows(); ++i) {
                    if (0 < minValue[i]) {
                        indices[i] = k;
                        minValue[i] = 0;
                    }
                }
            }

            for (; it; ++it) {
                for (IndexType i(previousNonZero+1); i < it.row(); ++i) {
                    if (0 < minValue[i]) {
                        indices[i] = k;
                        minValue[i] = 0;
                    }
                }
                previousNonZero = it.row();
                if ((!mpfr::isnan(it.value())) && (it.value() < minValue[it.row()])) {
                    indices[it.row()] = k;
                    minValue[it.row()] = it.value();
                }
            }

            // If there are still zeros
            for (IndexType i(previousNonZero+1); i < matrixR.rows(); ++i) {
                if (0 < minValue[i]) {
                    indices[i] = k;
                    minValue[i] = 0;
                }
            }
        }

        // We copy the minimums
        for (IndexType i = 0; i < matrixR.rows(); ++i) {
            if (minValue[i] != 0) {
                result.matrixR.insert(i,0) = minValue[i];
            }
        }
        result.matrixR.makeCompressed();
    }

    return result;
}

// element-wise minimum c = min(a, b)
SparseGmpEigenMatrix SparseGmpEigenMatrix::ewMin(const SparseGmpEigenMatrix& b) const
{
    SparseGmpEigenMatrix result;

    if ((numel() != 1) && (b.numel() == 1)) {
        // setting the output size
        result.matrixR.resize(matrixR.rows(), matrixR.cols());
        if ((isComplex) || (b.isComplex))
            result.matrixI.resize(matrixR.rows(), matrixR.cols());

        // We will first copy the data to triplets
        vector< Triplet<mpreal> > tripletListR, tripletListI;
        tripletListR.reserve(matrixR.nonZeros()); // WARNING : Here, we assume that the scalar b is not smaller than 0 (otherwise it will make the sparse matrix full...)
        if ((isComplex) || (b.isComplex))
            tripletListI.reserve(matrixI.nonZeros());

        if (b.isComplex) {
            if (isComplex) {
                for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
                    SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
                    SparseMatrix<mpreal>::InnerIterator itI(matrixI,k);

                    while ((itR) || (itI)) {
                        IndexType row;
                        if ((itR) && (itI))
                            row = min(itR.row(), itI.row());
                        else if (itR)
                            row = itR.row();
                        else
                            row = itI.row();

                        if ((!itI) || ((itR) && (itR.row() < itI.row()))) {
                            if ((mpfr::isnan(pow(b.matrixR.coeff(0,0), 2) + pow(b.matrixI.coeff(0,0), 2))) || ((pow(itR.value(), 2) < pow(b.matrixR.coeff(0,0), 2) + pow(b.matrixI.coeff(0,0), 2)) || ((pow(itR.value(), 2) == pow(b.matrixR.coeff(0,0), 2) + pow(b.matrixI.coeff(0,0), 2)) && (atan2(0, itR.value()) <= atan2(b.matrixI.coeff(0,0), b.matrixR.coeff(0,0)))))) {
                                tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                            } else {
                                tripletListR.push_back(Triplet<mpreal>(itR.row(), k, b.matrixR.coeff(0,0)));
                                tripletListI.push_back(Triplet<mpreal>(itR.row(), k, b.matrixI.coeff(0,0)));
                            }
                            ++itR;
                        } else if ((itR) && (itI) && (itR.row() == itI.row())) {
                            if ((mpfr::isnan(pow(b.matrixR.coeff(0,0), 2) + pow(b.matrixI.coeff(0,0), 2))) || ((pow(itR.value(), 2) + pow(itI.value(), 2) < pow(b.matrixR.coeff(0,0), 2) + pow(b.matrixI.coeff(0,0), 2)) || ((pow(itR.value(), 2) + pow(itI.value(), 2) == pow(b.matrixR.coeff(0,0), 2) + pow(b.matrixI.coeff(0,0), 2)) && (atan2(itI.value(), itR.value()) <= atan2(b.matrixI.coeff(0,0), b.matrixR.coeff(0,0)))))) {
                                tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                                tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itI.value()));
                            } else {
                                tripletListR.push_back(Triplet<mpreal>(itR.row(), k, b.matrixR.coeff(0,0)));
                                tripletListI.push_back(Triplet<mpreal>(itR.row(), k, b.matrixI.coeff(0,0)));
                            }
                            ++itR;
                            ++itI;
                        } else {
                            if ((mpfr::isnan(pow(b.matrixR.coeff(0,0), 2) + pow(b.matrixI.coeff(0,0), 2))) || ((pow(itI.value(), 2) < pow(b.matrixR.coeff(0,0), 2) + pow(b.matrixI.coeff(0,0), 2)) || ((pow(itI.value(), 2) == pow(b.matrixR.coeff(0,0), 2) + pow(b.matrixI.coeff(0,0), 2)) && (atan2(itI.value(), 0) <= atan2(b.matrixI.coeff(0,0), b.matrixR.coeff(0,0)))))) {
                                tripletListI.push_back(Triplet<mpreal>(itI.row(), k, itI.value()));
                            } else {
                                tripletListR.push_back(Triplet<mpreal>(itI.row(), k, b.matrixR.coeff(0,0)));
                                tripletListI.push_back(Triplet<mpreal>(itI.row(), k, b.matrixI.coeff(0,0)));
                            }
                            ++itI;
                        }
                    }
                }
            } else {
                for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
                    SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);

                    while (itR) {
                        if ((mpfr::isnan(pow(b.matrixR.coeff(0,0), 2) + pow(b.matrixI.coeff(0,0), 2))) || ((pow(itR.value(), 2) < pow(b.matrixR.coeff(0,0), 2) + pow(b.matrixI.coeff(0,0), 2)) || ((pow(itR.value(), 2) == pow(b.matrixR.coeff(0,0), 2) + pow(b.matrixI.coeff(0,0), 2)) && (atan2(0, itR.value()) <= atan2(b.matrixI.coeff(0,0), b.matrixR.coeff(0,0)))))) {
                            tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                        } else {
                            tripletListR.push_back(Triplet<mpreal>(itR.row(), k, b.matrixR.coeff(0,0)));
                            tripletListI.push_back(Triplet<mpreal>(itR.row(), k, b.matrixI.coeff(0,0)));
                        }
                        ++itR;
                    }
                }
            }
        } else if (isComplex) {
            for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
                SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
                SparseMatrix<mpreal>::InnerIterator itI(matrixI,k);

                while ((itR) || (itI)) {
                    IndexType row;
                    if ((itR) && (itI))
                        row = min(itR.row(), itI.row());
                    else if (itR)
                        row = itR.row();
                    else
                        row = itI.row();

                    if ((!itI) || ((itR) && (itR.row() < itI.row()))) {
                        if ((mpfr::isnan(pow(b.matrixR.coeff(0,0), 2))) || ((pow(itR.value(), 2) < pow(b.matrixR.coeff(0,0), 2)) || ((pow(itR.value(), 2) == pow(b.matrixR.coeff(0,0), 2)) && (atan2(0, itR.value()) <= atan2(0, b.matrixR.coeff(0,0)))))) {
                            tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                        } else {
                            tripletListR.push_back(Triplet<mpreal>(itR.row(), k, b.matrixR.coeff(0,0)));
                        }
                        ++itR;
                    } else if ((itR) && (itI) && (itR.row() == itI.row())) {
                        if ((mpfr::isnan(pow(b.matrixR.coeff(0,0), 2))) || ((pow(itR.value(), 2) + pow(itI.value(), 2) < pow(b.matrixR.coeff(0,0), 2)) || ((pow(itR.value(), 2) + pow(itI.value(), 2) == pow(b.matrixR.coeff(0,0), 2)) && (atan2(itI.value(), itR.value()) <= atan2(0, b.matrixR.coeff(0,0)))))) {
                            tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                            tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itI.value()));
                        } else {
                            tripletListR.push_back(Triplet<mpreal>(itR.row(), k, b.matrixR.coeff(0,0)));
                        }
                        ++itR;
                        ++itI;
                    } else {
                        if ((mpfr::isnan(pow(b.matrixR.coeff(0,0), 2))) || ((pow(itI.value(), 2) < pow(b.matrixR.coeff(0,0), 2)) || ((pow(itI.value(), 2) == pow(b.matrixR.coeff(0,0), 2)) && (atan2(itI.value(), 0) <= atan2(0, b.matrixR.coeff(0,0)))))) {
                            tripletListI.push_back(Triplet<mpreal>(itI.row(), k, itI.value()));
                        } else {
                            tripletListR.push_back(Triplet<mpreal>(itI.row(), k, b.matrixR.coeff(0,0)));
                        }
                        ++itI;
                    }
                }
            }
        } else {
            for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
                SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);

                while (itR) {
                    if ((mpfr::isnan(b.matrixR.coeff(0,0))) || (itR.value() <= b.matrixR.coeff(0,0)))
                        tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                    else
                        tripletListR.push_back(Triplet<mpreal>(itR.row(), k, b.matrixR.coeff(0,0)));
                    ++itR;
                }
            }
        }

        // Now we can assign the data to the sparse matrix
        result.matrixR.setFromTriplets(tripletListR.begin(), tripletListR.end());
        if ((isComplex) || (b.isComplex))
            result.matrixI.setFromTriplets(tripletListI.begin(), tripletListI.end());

        result.checkComplexity();
    } else if ((numel() == 1) && (b.numel() != 1)) {
        result = b.ewMin(*this);
    } else {
        // setting the output size
        result.matrixR.resize(matrixR.rows(), matrixR.cols());
        if ((isComplex) || (b.isComplex))
            result.matrixI.resize(matrixR.rows(), matrixR.cols());

        // We will first copy the data to triplets
        vector< Triplet<mpreal> > tripletListR, tripletListI;
        tripletListR.reserve(matrixR.nonZeros() + b.matrixR.nonZeros());
        if ((isComplex) || (b.isComplex))
            tripletListI.reserve(matrixI.nonZeros() + b.matrixI.nonZeros());

        // Now for each column, we merge the lists of lines with nonzero elements of
        // both matrices. The cases are slightly different depending on the complexity
        // of both matrices
        if (b.isComplex) {
            if (isComplex) {
                for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
                    SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
                    SparseMatrix<mpreal>::InnerIterator itI(matrixI,k);
                    SparseMatrix<mpreal>::InnerIterator itRb(b.matrixR,k);
                    SparseMatrix<mpreal>::InnerIterator itIb(b.matrixI,k);

                    while ((itR) || (itI) || (itRb) || (itIb)) {
                        IndexType row, rowb;
                        if ((itR) || (itI)) {
                            if ((itR) && (itI))
                                row = min(itR.row(), itI.row());
                            else if (itR)
                                row = itR.row();
                            else
                                row = itI.row();
                        }
                        if ((itRb) || (itIb)) {
                            if ((itRb) && (itIb))
                                rowb = min(itRb.row(), itIb.row());
                            else if (itRb)
                                rowb = itRb.row();
                            else
                                rowb = itIb.row();
                        }

                        if (((((itR) || (itI)) && ((itRb) || (itIb))) && (row < rowb)) || ((!itRb) && (!itIb))) {
                            if ((!itI) || ((itR) && (itR.row() < itI.row()))) {
                                if (mpfr::abs(itR.value()) < 0)
                                    tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                                ++itR;
                            } else if ((itR) && (itI) && (itR.row() == itI.row())) {
                                if (pow(itR.value(), 2) + pow(itI.value(), 2) < 0) {
                                    tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                                    tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itI.value()));
                                }
                                ++itR;
                                ++itI;
                            } else {
                                if (mpfr::abs(itI.value()) < 0)
                                    tripletListI.push_back(Triplet<mpreal>(itI.row(), k, itI.value()));
                                ++itI;
                            }
                        } else if ((((itR) || (itI)) && ((itRb) || (itIb))) && (row == rowb)) {
                            if ((!itI) || ((itR) && (itR.row() < itI.row()))) {
                                if ((!itIb) || ((itRb) && (itRb.row() < itIb.row()))) {
                                    if ((mpfr::isnan(itRb.value())) || (mpfr::abs(itR.value()) <= mpfr::abs(itRb.value())))
                                        tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                                    else
                                        tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itRb.value()));
                                    ++itRb;
                                } else if ((itRb) && (itIb) && (itRb.row() == itIb.row())) {
                                    if ((mpfr::isnan(pow(itRb.value(), 2) + pow(itIb.value(), 2))) || ((pow(itR.value(), 2) < pow(itRb.value(), 2) + pow(itIb.value(), 2)) || ((pow(itR.value(), 2) == pow(itRb.value(), 2) + pow(itIb.value(), 2)) && (atan2(0, itR.value()) <= atan2(itIb.value(), itRb.value()))))) {
                                        tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                                    } else {
                                        tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itRb.value()));
                                        tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itIb.value()));
                                    }
                                    ++itRb;
                                    ++itIb;
                                } else {
                                    if ((mpfr::isnan(itIb.value())) || ((pow(itR.value(), 2) < pow(itIb.value(), 2)) || ((pow(itR.value(), 2) == pow(itIb.value(), 2)) && (atan2(0, itR.value()) <= atan2(itIb.value(), 0)))))
                                        tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                                    else
                                        tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itIb.value()));
                                    ++itIb;
                                }
                                ++itR;
                            } else if ((itR) && (itI) && (itR.row() == itI.row())) {
                                if ((!itIb) || ((itRb) && (itRb.row() < itIb.row()))) {
                                    if ((mpfr::isnan(itRb.value())) || ((pow(itR.value(), 2) + pow(itI.value(), 2) < pow(itRb.value(), 2)) || ((pow(itR.value(), 2) + pow(itI.value(), 2) == pow(itRb.value(), 2)) && (atan2(itI.value(), itR.value()) <= atan2(0, itRb.value()))))) {
                                        tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                                        tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itI.value()));
                                    } else {
                                        tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itRb.value()));
                                    }
                                    ++itRb;
                                } else if ((itRb) && (itIb) && (itRb.row() == itIb.row())) {
                                    if ((mpfr::isnan(pow(itRb.value(), 2) + pow(itIb.value(), 2))) || ((pow(itR.value(), 2) + pow(itI.value(), 2) < pow(itRb.value(), 2) + pow(itIb.value(), 2)) || ((pow(itR.value(), 2) + pow(itI.value(), 2) == pow(itRb.value(), 2) + pow(itIb.value(), 2)) && (atan2(itI.value(), itR.value()) <= atan2(itIb.value(), itRb.value()))))) {
                                        tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                                        tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itI.value()));
                                    } else {
                                        tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itRb.value()));
                                        tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itIb.value()));
                                    }
                                    ++itRb;
                                    ++itIb;
                                } else {
                                    if ((mpfr::isnan(itIb.value())) || ((pow(itR.value(), 2) + pow(itI.value(), 2) < pow(itIb.value(), 2)) || ((pow(itR.value(), 2) + pow(itI.value(), 2) == pow(itIb.value(), 2)) && (atan2(itI.value(), itR.value()) <= atan2(itIb.value(), 0))))) {
                                        tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                                        tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itI.value()));
                                    } else {
                                        tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itIb.value()));
                                    }
                                    ++itIb;
                                }
                                ++itR;
                                ++itI;
                            } else {
                                if ((!itIb) || ((itRb) && (itRb.row() < itIb.row()))) {
                                    if ((mpfr::isnan(itRb.value())) || ((pow(itI.value(), 2) < pow(itRb.value(), 2)) || ((pow(itI.value(), 2) == pow(itRb.value(), 2)) && (atan2(itI.value(), 0) <= atan2(0, itRb.value()))))) {
                                        tripletListI.push_back(Triplet<mpreal>(itI.row(), k, itI.value()));
                                    } else {
                                        tripletListR.push_back(Triplet<mpreal>(itI.row(), k, itRb.value()));
                                    }
                                    ++itRb;
                                } else if ((itRb) && (itIb) && (itRb.row() == itIb.row())) {
                                    if ((mpfr::isnan(pow(itRb.value(), 2) + pow(itIb.value(), 2))) || ((pow(itI.value(), 2) < pow(itRb.value(), 2) + pow(itIb.value(), 2)) || ((pow(itI.value(), 2) == pow(itRb.value(), 2) + pow(itIb.value(), 2)) && (atan2(itI.value(), 0) <= atan2(itIb.value(), itRb.value()))))) {
                                        tripletListI.push_back(Triplet<mpreal>(itI.row(), k, itI.value()));
                                    } else {
                                        tripletListR.push_back(Triplet<mpreal>(itI.row(), k, itRb.value()));
                                        tripletListI.push_back(Triplet<mpreal>(itI.row(), k, itIb.value()));
                                    }
                                    ++itRb;
                                    ++itIb;
                                } else {
                                    if ((mpfr::isnan(itIb.value())) || ((pow(itI.value(), 2) < pow(itIb.value(), 2)) || ((pow(itI.value(), 2) == pow(itIb.value(), 2)) && (atan2(itI.value(), 0) <= atan2(itIb.value(), 0))))) {
                                        tripletListI.push_back(Triplet<mpreal>(itI.row(), k, itI.value()));
                                    } else {
                                        tripletListI.push_back(Triplet<mpreal>(itI.row(), k, itIb.value()));
                                    }
                                    ++itIb;
                                }
                                ++itI;
                            }
                        } else {
                            if ((!itIb) || ((itRb) && (itRb.row() < itIb.row()))) {
                                if (mpfr::abs(itRb.value()) < 0)
                                    tripletListR.push_back(Triplet<mpreal>(itRb.row(), k, itRb.value()));
                                ++itRb;
                            } else if ((itRb) && (itIb) && (itRb.row() == itIb.row())) {
                                if (pow(itRb.value(), 2) + pow(itIb.value(), 2) < 0) {
                                    tripletListR.push_back(Triplet<mpreal>(itRb.row(), k, itRb.value()));
                                    tripletListI.push_back(Triplet<mpreal>(itRb.row(), k, itIb.value()));
                                }
                                ++itRb;
                                ++itIb;
                            } else {
                                if (mpfr::abs(itIb.value()) < 0)
                                    tripletListI.push_back(Triplet<mpreal>(itIb.row(), k, itIb.value()));
                                ++itIb;
                            }
                        }
                    }
                }
            } else {
                for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
                    SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
                    SparseMatrix<mpreal>::InnerIterator itRb(b.matrixR,k);
                    SparseMatrix<mpreal>::InnerIterator itIb(b.matrixI,k);

                    while ((itR) || (itRb) || (itIb)) {
                        IndexType row, rowb;
                        if (itR)
                            row = itR.row();
                        if ((itRb) || (itIb)) {
                            if ((itRb) && (itIb))
                                rowb = min(itRb.row(), itIb.row());
                            else if (itRb)
                                rowb = itRb.row();
                            else
                                rowb = itIb.row();
                        }

                        if ((((itR) && ((itRb) || (itIb))) && (row < rowb)) || ((!itRb) && (!itIb))) {
                            if (mpfr::abs(itR.value()) < 0)
                                tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                            ++itR;
                        } else if (((itR) && ((itRb) || (itIb))) && (row == rowb)) {
                            if ((!itIb) || ((itRb) && (itRb.row() < itIb.row()))) {
                                if ((mpfr::isnan(itRb.value())) || (mpfr::abs(itR.value()) <= mpfr::abs(itRb.value())))
                                    tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                                else
                                    tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itRb.value()));
                                ++itRb;
                            } else if ((itRb) && (itIb) && (itRb.row() == itIb.row())) {
                                if ((mpfr::isnan(pow(itRb.value(), 2) + pow(itIb.value(), 2))) || ((pow(itR.value(), 2) < pow(itRb.value(), 2) + pow(itIb.value(), 2)) || ((pow(itR.value(), 2) == pow(itRb.value(), 2) + pow(itIb.value(), 2)) && (atan2(0, itR.value()) <= atan2(itIb.value(), itRb.value()))))) {
                                    tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                                } else {
                                    tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itRb.value()));
                                    tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itIb.value()));
                                }
                                ++itRb;
                                ++itIb;
                            } else {
                                if ((mpfr::isnan(itIb.value())) || ((pow(itR.value(), 2) < pow(itIb.value(), 2)) || ((pow(itR.value(), 2) == pow(itIb.value(), 2)) && (atan2(0, itR.value()) <= atan2(itIb.value(), 0)))))
                                    tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                                else
                                    tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itIb.value()));
                                ++itIb;
                            }
                            ++itR;
                        } else {
                            if ((!itIb) || ((itRb) && (itRb.row() < itIb.row()))) {
                                if (mpfr::abs(itRb.value()) < 0)
                                    tripletListR.push_back(Triplet<mpreal>(itRb.row(), k, itRb.value()));
                                ++itRb;
                            } else if ((itRb) && (itIb) && (itRb.row() == itIb.row())) {
                                if (pow(itRb.value(), 2) + pow(itIb.value(), 2) < 0) {
                                    tripletListR.push_back(Triplet<mpreal>(itRb.row(), k, itRb.value()));
                                    tripletListI.push_back(Triplet<mpreal>(itRb.row(), k, itIb.value()));
                                }
                                ++itRb;
                                ++itIb;
                            } else {
                                if (mpfr::abs(itIb.value()) < 0)
                                    tripletListI.push_back(Triplet<mpreal>(itIb.row(), k, itIb.value()));
                                ++itIb;
                            }
                        }
                    }
                }
            }
        } else if (isComplex) {
            for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
                SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
                SparseMatrix<mpreal>::InnerIterator itI(matrixI,k);
                SparseMatrix<mpreal>::InnerIterator itRb(b.matrixR,k);

                while ((itR) || (itI) || (itRb)) {
                    IndexType row, rowb;
                    if ((itR) || (itI)) {
                        if ((itR) && (itI))
                            row = min(itR.row(), itI.row());
                        else if (itR)
                            row = itR.row();
                        else
                            row = itI.row();
                    }
                    if (itRb)
                        rowb = itRb.row();

                    if (((((itR) || (itI)) && (itRb)) && (row < rowb)) || (!itRb)) {
                        if ((!itI) || ((itR) && (itR.row() < itI.row()))) {
                            if (mpfr::abs(itR.value()) < 0)
                                tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                            ++itR;
                        } else if ((itR) && (itI) && (itR.row() == itI.row())) {
                            if (pow(itR.value(), 2) + pow(itI.value(), 2) < 0) {
                                tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                                tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itI.value()));
                            }
                            ++itR;
                            ++itI;
                        } else {
                            if (mpfr::abs(itI.value()) < 0)
                                tripletListI.push_back(Triplet<mpreal>(itI.row(), k, itI.value()));
                            ++itI;
                        }
                    } else if ((((itR) || (itI)) && (itRb)) && (row == rowb)) {
                        if ((!itI) || ((itR) && (itR.row() < itI.row()))) {
                            if ((mpfr::isnan(itRb.value())) || (mpfr::abs(itR.value()) <= mpfr::abs(itRb.value())))
                                tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                            else
                                tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itRb.value()));
                            ++itRb;
                            ++itR;
                        } else if ((itR) && (itI) && (itR.row() == itI.row())) {
                            if ((mpfr::isnan(itRb.value())) || ((pow(itR.value(), 2) + pow(itI.value(), 2) < pow(itRb.value(), 2)) || ((pow(itR.value(), 2) + pow(itI.value(), 2) == pow(itRb.value(), 2)) && (atan2(itI.value(), itR.value()) <= atan2(0, itRb.value()))))) {
                                tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                                tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itI.value()));
                            } else {
                                tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itRb.value()));
                            }
                            ++itRb;
                            ++itR;
                            ++itI;
                        } else {
                            if ((mpfr::isnan(itRb.value())) || ((pow(itI.value(), 2) < pow(itRb.value(), 2)) || ((pow(itI.value(), 2) == pow(itRb.value(), 2)) && (atan2(itI.value(), 0) <= atan2(0, itRb.value()))))) {
                                tripletListI.push_back(Triplet<mpreal>(itI.row(), k, itI.value()));
                            } else {
                                tripletListR.push_back(Triplet<mpreal>(itI.row(), k, itRb.value()));
                            }
                            ++itRb;
                            ++itI;
                        }
                    } else {
                        if (mpfr::abs(itRb.value()) < 0)
                            tripletListR.push_back(Triplet<mpreal>(itRb.row(), k, itRb.value()));
                        ++itRb;
                    }
                }
            }
        } else {
            for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
                SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
                SparseMatrix<mpreal>::InnerIterator itRb(b.matrixR,k);

                while ((itR) || (itRb)) {
                    if ((!itRb) || ((itR) && (itR.row() < itRb.row()))) {
                        if (itR.value() < 0)
                            tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                        ++itR;
                    } else if ((itR) && (itRb) && (itR.row() == itRb.row())) {
                        if ((mpfr::isnan(itRb.value())) || (itR.value() <= itRb.value()))
                            tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                        else
                            tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itRb.value()));
                        ++itRb;
                        ++itR;
                    } else {
                        if (itRb.value() < 0)
                            tripletListR.push_back(Triplet<mpreal>(itRb.row(), k, itRb.value()));
                        ++itRb;
                    }
                }
            }
        }

        // Now we can assign the data to the sparse matrix
        result.matrixR.setFromTriplets(tripletListR.begin(), tripletListR.end());
        if ((isComplex) || (b.isComplex))
            result.matrixI.setFromTriplets(tripletListI.begin(), tripletListI.end());

        result.checkComplexity();
    }

    return result;
}

// element-wise minimum c = min(a, b)
SparseGmpEigenMatrix& SparseGmpEigenMatrix::ewMin_new(const SparseGmpEigenMatrix& b) const
{
    SparseGmpEigenMatrix& result(*(new SparseGmpEigenMatrix));

    if ((numel() != 1) && (b.numel() == 1)) {
        // setting the output size
        result.matrixR.resize(matrixR.rows(), matrixR.cols());
        if ((isComplex) || (b.isComplex))
            result.matrixI.resize(matrixR.rows(), matrixR.cols());

        // We will first copy the data to triplets
        vector< Triplet<mpreal> > tripletListR, tripletListI;
        tripletListR.reserve(matrixR.nonZeros()); // WARNING : Here, we assume that the scalar b is not smaller than 0 (otherwise it will make the sparse matrix full...)
        if ((isComplex) || (b.isComplex))
            tripletListI.reserve(matrixI.nonZeros());

        if (b.isComplex) {
            if (isComplex) {
                for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
                    SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
                    SparseMatrix<mpreal>::InnerIterator itI(matrixI,k);

                    while ((itR) || (itI)) {
                        IndexType row;
                        if ((itR) && (itI))
                            row = min(itR.row(), itI.row());
                        else if (itR)
                            row = itR.row();
                        else
                            row = itI.row();

                        if ((!itI) || ((itR) && (itR.row() < itI.row()))) {
                            if ((mpfr::isnan(pow(b.matrixR.coeff(0,0), 2) + pow(b.matrixI.coeff(0,0), 2))) || ((pow(itR.value(), 2) < pow(b.matrixR.coeff(0,0), 2) + pow(b.matrixI.coeff(0,0), 2)) || ((pow(itR.value(), 2) == pow(b.matrixR.coeff(0,0), 2) + pow(b.matrixI.coeff(0,0), 2)) && (atan2(0, itR.value()) <= atan2(b.matrixI.coeff(0,0), b.matrixR.coeff(0,0)))))) {
                                tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                            } else {
                                tripletListR.push_back(Triplet<mpreal>(itR.row(), k, b.matrixR.coeff(0,0)));
                                tripletListI.push_back(Triplet<mpreal>(itR.row(), k, b.matrixI.coeff(0,0)));
                            }
                            ++itR;
                        } else if ((itR) && (itI) && (itR.row() == itI.row())) {
                            if ((mpfr::isnan(pow(b.matrixR.coeff(0,0), 2) + pow(b.matrixI.coeff(0,0), 2))) || ((pow(itR.value(), 2) + pow(itI.value(), 2) < pow(b.matrixR.coeff(0,0), 2) + pow(b.matrixI.coeff(0,0), 2)) || ((pow(itR.value(), 2) + pow(itI.value(), 2) == pow(b.matrixR.coeff(0,0), 2) + pow(b.matrixI.coeff(0,0), 2)) && (atan2(itI.value(), itR.value()) <= atan2(b.matrixI.coeff(0,0), b.matrixR.coeff(0,0)))))) {
                                tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                                tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itI.value()));
                            } else {
                                tripletListR.push_back(Triplet<mpreal>(itR.row(), k, b.matrixR.coeff(0,0)));
                                tripletListI.push_back(Triplet<mpreal>(itR.row(), k, b.matrixI.coeff(0,0)));
                            }
                            ++itR;
                            ++itI;
                        } else {
                            if ((mpfr::isnan(pow(b.matrixR.coeff(0,0), 2) + pow(b.matrixI.coeff(0,0), 2))) || ((pow(itI.value(), 2) < pow(b.matrixR.coeff(0,0), 2) + pow(b.matrixI.coeff(0,0), 2)) || ((pow(itI.value(), 2) == pow(b.matrixR.coeff(0,0), 2) + pow(b.matrixI.coeff(0,0), 2)) && (atan2(itI.value(), 0) <= atan2(b.matrixI.coeff(0,0), b.matrixR.coeff(0,0)))))) {
                                tripletListI.push_back(Triplet<mpreal>(itI.row(), k, itI.value()));
                            } else {
                                tripletListR.push_back(Triplet<mpreal>(itI.row(), k, b.matrixR.coeff(0,0)));
                                tripletListI.push_back(Triplet<mpreal>(itI.row(), k, b.matrixI.coeff(0,0)));
                            }
                            ++itI;
                        }
                    }
                }
            } else {
                for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
                    SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);

                    while (itR) {
                        if ((mpfr::isnan(pow(b.matrixR.coeff(0,0), 2) + pow(b.matrixI.coeff(0,0), 2))) || ((pow(itR.value(), 2) < pow(b.matrixR.coeff(0,0), 2) + pow(b.matrixI.coeff(0,0), 2)) || ((pow(itR.value(), 2) == pow(b.matrixR.coeff(0,0), 2) + pow(b.matrixI.coeff(0,0), 2)) && (atan2(0, itR.value()) <= atan2(b.matrixI.coeff(0,0), b.matrixR.coeff(0,0)))))) {
                            tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                        } else {
                            tripletListR.push_back(Triplet<mpreal>(itR.row(), k, b.matrixR.coeff(0,0)));
                            tripletListI.push_back(Triplet<mpreal>(itR.row(), k, b.matrixI.coeff(0,0)));
                        }
                        ++itR;
                    }
                }
            }
        } else if (isComplex) {
            for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
                SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
                SparseMatrix<mpreal>::InnerIterator itI(matrixI,k);

                while ((itR) || (itI)) {
                    IndexType row;
                    if ((itR) && (itI))
                        row = min(itR.row(), itI.row());
                    else if (itR)
                        row = itR.row();
                    else
                        row = itI.row();

                    if ((!itI) || ((itR) && (itR.row() < itI.row()))) {
                        if ((mpfr::isnan(pow(b.matrixR.coeff(0,0), 2))) || ((pow(itR.value(), 2) < pow(b.matrixR.coeff(0,0), 2)) || ((pow(itR.value(), 2) == pow(b.matrixR.coeff(0,0), 2)) && (atan2(0, itR.value()) <= atan2(0, b.matrixR.coeff(0,0)))))) {
                            tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                        } else {
                            tripletListR.push_back(Triplet<mpreal>(itR.row(), k, b.matrixR.coeff(0,0)));
                        }
                        ++itR;
                    } else if ((itR) && (itI) && (itR.row() == itI.row())) {
                        if ((mpfr::isnan(pow(b.matrixR.coeff(0,0), 2))) || ((pow(itR.value(), 2) + pow(itI.value(), 2) < pow(b.matrixR.coeff(0,0), 2)) || ((pow(itR.value(), 2) + pow(itI.value(), 2) == pow(b.matrixR.coeff(0,0), 2)) && (atan2(itI.value(), itR.value()) <= atan2(0, b.matrixR.coeff(0,0)))))) {
                            tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                            tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itI.value()));
                        } else {
                            tripletListR.push_back(Triplet<mpreal>(itR.row(), k, b.matrixR.coeff(0,0)));
                        }
                        ++itR;
                        ++itI;
                    } else {
                        if ((mpfr::isnan(pow(b.matrixR.coeff(0,0), 2))) || ((pow(itI.value(), 2) < pow(b.matrixR.coeff(0,0), 2)) || ((pow(itI.value(), 2) == pow(b.matrixR.coeff(0,0), 2)) && (atan2(itI.value(), 0) <= atan2(0, b.matrixR.coeff(0,0)))))) {
                            tripletListI.push_back(Triplet<mpreal>(itI.row(), k, itI.value()));
                        } else {
                            tripletListR.push_back(Triplet<mpreal>(itI.row(), k, b.matrixR.coeff(0,0)));
                        }
                        ++itI;
                    }
                }
            }
        } else {
            for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
                SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);

                while (itR) {
                    if ((mpfr::isnan(b.matrixR.coeff(0,0))) || (itR.value() <= b.matrixR.coeff(0,0)))
                        tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                    else
                        tripletListR.push_back(Triplet<mpreal>(itR.row(), k, b.matrixR.coeff(0,0)));
                    ++itR;
                }
            }
        }

        // Now we can assign the data to the sparse matrix
        result.matrixR.setFromTriplets(tripletListR.begin(), tripletListR.end());
        if ((isComplex) || (b.isComplex))
            result.matrixI.setFromTriplets(tripletListI.begin(), tripletListI.end());

        result.checkComplexity();
    } else if ((numel() == 1) && (b.numel() != 1)) {
        result = b.ewMin(*this);
    } else {
        // setting the output size
        result.matrixR.resize(matrixR.rows(), matrixR.cols());
        if ((isComplex) || (b.isComplex))
            result.matrixI.resize(matrixR.rows(), matrixR.cols());

        // We will first copy the data to triplets
        vector< Triplet<mpreal> > tripletListR, tripletListI;
        tripletListR.reserve(matrixR.nonZeros() + b.matrixR.nonZeros());
        if ((isComplex) || (b.isComplex))
            tripletListI.reserve(matrixI.nonZeros() + b.matrixI.nonZeros());

        // Now for each column, we merge the lists of lines with nonzero elements of
        // both matrices. The cases are slightly different depending on the complexity
        // of both matrices
        if (b.isComplex) {
            if (isComplex) {
                for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
                    SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
                    SparseMatrix<mpreal>::InnerIterator itI(matrixI,k);
                    SparseMatrix<mpreal>::InnerIterator itRb(b.matrixR,k);
                    SparseMatrix<mpreal>::InnerIterator itIb(b.matrixI,k);

                    while ((itR) || (itI) || (itRb) || (itIb)) {
                        IndexType row, rowb;
                        if ((itR) || (itI)) {
                            if ((itR) && (itI))
                                row = min(itR.row(), itI.row());
                            else if (itR)
                                row = itR.row();
                            else
                                row = itI.row();
                        }
                        if ((itRb) || (itIb)) {
                            if ((itRb) && (itIb))
                                rowb = min(itRb.row(), itIb.row());
                            else if (itRb)
                                rowb = itRb.row();
                            else
                                rowb = itIb.row();
                        }

                        if (((((itR) || (itI)) && ((itRb) || (itIb))) && (row < rowb)) || ((!itRb) && (!itIb))) {
                            if ((!itI) || ((itR) && (itR.row() < itI.row()))) {
                                if (mpfr::abs(itR.value()) < 0)
                                    tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                                ++itR;
                            } else if ((itR) && (itI) && (itR.row() == itI.row())) {
                                if (pow(itR.value(), 2) + pow(itI.value(), 2) < 0) {
                                    tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                                    tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itI.value()));
                                }
                                ++itR;
                                ++itI;
                            } else {
                                if (mpfr::abs(itI.value()) < 0)
                                    tripletListI.push_back(Triplet<mpreal>(itI.row(), k, itI.value()));
                                ++itI;
                            }
                        } else if ((((itR) || (itI)) && ((itRb) || (itIb))) && (row == rowb)) {
                            if ((!itI) || ((itR) && (itR.row() < itI.row()))) {
                                if ((!itIb) || ((itRb) && (itRb.row() < itIb.row()))) {
                                    if ((mpfr::isnan(itRb.value())) || (mpfr::abs(itR.value()) <= mpfr::abs(itRb.value())))
                                        tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                                    else
                                        tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itRb.value()));
                                    ++itRb;
                                } else if ((itRb) && (itIb) && (itRb.row() == itIb.row())) {
                                    if ((mpfr::isnan(pow(itRb.value(), 2) + pow(itIb.value(), 2))) || ((pow(itR.value(), 2) < pow(itRb.value(), 2) + pow(itIb.value(), 2)) || ((pow(itR.value(), 2) == pow(itRb.value(), 2) + pow(itIb.value(), 2)) && (atan2(0, itR.value()) <= atan2(itIb.value(), itRb.value()))))) {
                                        tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                                    } else {
                                        tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itRb.value()));
                                        tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itIb.value()));
                                    }
                                    ++itRb;
                                    ++itIb;
                                } else {
                                    if ((mpfr::isnan(itIb.value())) || ((pow(itR.value(), 2) < pow(itIb.value(), 2)) || ((pow(itR.value(), 2) == pow(itIb.value(), 2)) && (atan2(0, itR.value()) <= atan2(itIb.value(), 0)))))
                                        tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                                    else
                                        tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itIb.value()));
                                    ++itIb;
                                }
                                ++itR;
                            } else if ((itR) && (itI) && (itR.row() == itI.row())) {
                                if ((!itIb) || ((itRb) && (itRb.row() < itIb.row()))) {
                                    if ((mpfr::isnan(itRb.value())) || ((pow(itR.value(), 2) + pow(itI.value(), 2) < pow(itRb.value(), 2)) || ((pow(itR.value(), 2) + pow(itI.value(), 2) == pow(itRb.value(), 2)) && (atan2(itI.value(), itR.value()) <= atan2(0, itRb.value()))))) {
                                        tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                                        tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itI.value()));
                                    } else {
                                        tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itRb.value()));
                                    }
                                    ++itRb;
                                } else if ((itRb) && (itIb) && (itRb.row() == itIb.row())) {
                                    if ((mpfr::isnan(pow(itRb.value(), 2) + pow(itIb.value(), 2))) || ((pow(itR.value(), 2) + pow(itI.value(), 2) < pow(itRb.value(), 2) + pow(itIb.value(), 2)) || ((pow(itR.value(), 2) + pow(itI.value(), 2) == pow(itRb.value(), 2) + pow(itIb.value(), 2)) && (atan2(itI.value(), itR.value()) <= atan2(itIb.value(), itRb.value()))))) {
                                        tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                                        tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itI.value()));
                                    } else {
                                        tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itRb.value()));
                                        tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itIb.value()));
                                    }
                                    ++itRb;
                                    ++itIb;
                                } else {
                                    if ((mpfr::isnan(itIb.value())) || ((pow(itR.value(), 2) + pow(itI.value(), 2) < pow(itIb.value(), 2)) || ((pow(itR.value(), 2) + pow(itI.value(), 2) == pow(itIb.value(), 2)) && (atan2(itI.value(), itR.value()) <= atan2(itIb.value(), 0))))) {
                                        tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                                        tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itI.value()));
                                    } else {
                                        tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itIb.value()));
                                    }
                                    ++itIb;
                                }
                                ++itR;
                                ++itI;
                            } else {
                                if ((!itIb) || ((itRb) && (itRb.row() < itIb.row()))) {
                                    if ((mpfr::isnan(itRb.value())) || ((pow(itI.value(), 2) < pow(itRb.value(), 2)) || ((pow(itI.value(), 2) == pow(itRb.value(), 2)) && (atan2(itI.value(), 0) <= atan2(0, itRb.value()))))) {
                                        tripletListI.push_back(Triplet<mpreal>(itI.row(), k, itI.value()));
                                    } else {
                                        tripletListR.push_back(Triplet<mpreal>(itI.row(), k, itRb.value()));
                                    }
                                    ++itRb;
                                } else if ((itRb) && (itIb) && (itRb.row() == itIb.row())) {
                                    if ((mpfr::isnan(pow(itRb.value(), 2) + pow(itIb.value(), 2))) || ((pow(itI.value(), 2) < pow(itRb.value(), 2) + pow(itIb.value(), 2)) || ((pow(itI.value(), 2) == pow(itRb.value(), 2) + pow(itIb.value(), 2)) && (atan2(itI.value(), 0) <= atan2(itIb.value(), itRb.value()))))) {
                                        tripletListI.push_back(Triplet<mpreal>(itI.row(), k, itI.value()));
                                    } else {
                                        tripletListR.push_back(Triplet<mpreal>(itI.row(), k, itRb.value()));
                                        tripletListI.push_back(Triplet<mpreal>(itI.row(), k, itIb.value()));
                                    }
                                    ++itRb;
                                    ++itIb;
                                } else {
                                    if ((mpfr::isnan(itIb.value())) || ((pow(itI.value(), 2) < pow(itIb.value(), 2)) || ((pow(itI.value(), 2) == pow(itIb.value(), 2)) && (atan2(itI.value(), 0) <= atan2(itIb.value(), 0))))) {
                                        tripletListI.push_back(Triplet<mpreal>(itI.row(), k, itI.value()));
                                    } else {
                                        tripletListI.push_back(Triplet<mpreal>(itI.row(), k, itIb.value()));
                                    }
                                    ++itIb;
                                }
                                ++itI;
                            }
                        } else {
                            if ((!itIb) || ((itRb) && (itRb.row() < itIb.row()))) {
                                if (mpfr::abs(itRb.value()) < 0)
                                    tripletListR.push_back(Triplet<mpreal>(itRb.row(), k, itRb.value()));
                                ++itRb;
                            } else if ((itRb) && (itIb) && (itRb.row() == itIb.row())) {
                                if (pow(itRb.value(), 2) + pow(itIb.value(), 2) < 0) {
                                    tripletListR.push_back(Triplet<mpreal>(itRb.row(), k, itRb.value()));
                                    tripletListI.push_back(Triplet<mpreal>(itRb.row(), k, itIb.value()));
                                }
                                ++itRb;
                                ++itIb;
                            } else {
                                if (mpfr::abs(itIb.value()) < 0)
                                    tripletListI.push_back(Triplet<mpreal>(itIb.row(), k, itIb.value()));
                                ++itIb;
                            }
                        }
                    }
                }
            } else {
                for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
                    SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
                    SparseMatrix<mpreal>::InnerIterator itRb(b.matrixR,k);
                    SparseMatrix<mpreal>::InnerIterator itIb(b.matrixI,k);

                    while ((itR) || (itRb) || (itIb)) {
                        IndexType row, rowb;
                        if (itR)
                            row = itR.row();
                        if ((itRb) || (itIb)) {
                            if ((itRb) && (itIb))
                                rowb = min(itRb.row(), itIb.row());
                            else if (itRb)
                                rowb = itRb.row();
                            else
                                rowb = itIb.row();
                        }

                        if ((((itR) && ((itRb) || (itIb))) && (row < rowb)) || ((!itRb) && (!itIb))) {
                            if (mpfr::abs(itR.value()) < 0)
                                tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                            ++itR;
                        } else if (((itR) && ((itRb) || (itIb))) && (row == rowb)) {
                            if ((!itIb) || ((itRb) && (itRb.row() < itIb.row()))) {
                                if ((mpfr::isnan(itRb.value())) || (mpfr::abs(itR.value()) <= mpfr::abs(itRb.value())))
                                    tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                                else
                                    tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itRb.value()));
                                ++itRb;
                            } else if ((itRb) && (itIb) && (itRb.row() == itIb.row())) {
                                if ((mpfr::isnan(pow(itRb.value(), 2) + pow(itIb.value(), 2))) || ((pow(itR.value(), 2) < pow(itRb.value(), 2) + pow(itIb.value(), 2)) || ((pow(itR.value(), 2) == pow(itRb.value(), 2) + pow(itIb.value(), 2)) && (atan2(0, itR.value()) <= atan2(itIb.value(), itRb.value()))))) {
                                    tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                                } else {
                                    tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itRb.value()));
                                    tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itIb.value()));
                                }
                                ++itRb;
                                ++itIb;
                            } else {
                                if ((mpfr::isnan(itIb.value())) || ((pow(itR.value(), 2) < pow(itIb.value(), 2)) || ((pow(itR.value(), 2) == pow(itIb.value(), 2)) && (atan2(0, itR.value()) <= atan2(itIb.value(), 0)))))
                                    tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                                else
                                    tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itIb.value()));
                                ++itIb;
                            }
                            ++itR;
                        } else {
                            if ((!itIb) || ((itRb) && (itRb.row() < itIb.row()))) {
                                if (mpfr::abs(itRb.value()) < 0)
                                    tripletListR.push_back(Triplet<mpreal>(itRb.row(), k, itRb.value()));
                                ++itRb;
                            } else if ((itRb) && (itIb) && (itRb.row() == itIb.row())) {
                                if (pow(itRb.value(), 2) + pow(itIb.value(), 2) < 0) {
                                    tripletListR.push_back(Triplet<mpreal>(itRb.row(), k, itRb.value()));
                                    tripletListI.push_back(Triplet<mpreal>(itRb.row(), k, itIb.value()));
                                }
                                ++itRb;
                                ++itIb;
                            } else {
                                if (mpfr::abs(itIb.value()) < 0)
                                    tripletListI.push_back(Triplet<mpreal>(itIb.row(), k, itIb.value()));
                                ++itIb;
                            }
                        }
                    }
                }
            }
        } else if (isComplex) {
            for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
                SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
                SparseMatrix<mpreal>::InnerIterator itI(matrixI,k);
                SparseMatrix<mpreal>::InnerIterator itRb(b.matrixR,k);

                while ((itR) || (itI) || (itRb)) {
                    IndexType row, rowb;
                    if ((itR) || (itI)) {
                        if ((itR) && (itI))
                            row = min(itR.row(), itI.row());
                        else if (itR)
                            row = itR.row();
                        else
                            row = itI.row();
                    }
                    if (itRb)
                        rowb = itRb.row();

                    if (((((itR) || (itI)) && (itRb)) && (row < rowb)) || (!itRb)) {
                        if ((!itI) || ((itR) && (itR.row() < itI.row()))) {
                            if (mpfr::abs(itR.value()) < 0)
                                tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                            ++itR;
                        } else if ((itR) && (itI) && (itR.row() == itI.row())) {
                            if (pow(itR.value(), 2) + pow(itI.value(), 2) < 0) {
                                tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                                tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itI.value()));
                            }
                            ++itR;
                            ++itI;
                        } else {
                            if (mpfr::abs(itI.value()) < 0)
                                tripletListI.push_back(Triplet<mpreal>(itI.row(), k, itI.value()));
                            ++itI;
                        }
                    } else if ((((itR) || (itI)) && (itRb)) && (row == rowb)) {
                        if ((!itI) || ((itR) && (itR.row() < itI.row()))) {
                            if ((mpfr::isnan(itRb.value())) || (mpfr::abs(itR.value()) <= mpfr::abs(itRb.value())))
                                tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                            else
                                tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itRb.value()));
                            ++itRb;
                            ++itR;
                        } else if ((itR) && (itI) && (itR.row() == itI.row())) {
                            if ((mpfr::isnan(itRb.value())) || ((pow(itR.value(), 2) + pow(itI.value(), 2) < pow(itRb.value(), 2)) || ((pow(itR.value(), 2) + pow(itI.value(), 2) == pow(itRb.value(), 2)) && (atan2(itI.value(), itR.value()) <= atan2(0, itRb.value()))))) {
                                tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                                tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itI.value()));
                            } else {
                                tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itRb.value()));
                            }
                            ++itRb;
                            ++itR;
                            ++itI;
                        } else {
                            if ((mpfr::isnan(itRb.value())) || ((pow(itI.value(), 2) < pow(itRb.value(), 2)) || ((pow(itI.value(), 2) == pow(itRb.value(), 2)) && (atan2(itI.value(), 0) <= atan2(0, itRb.value()))))) {
                                tripletListI.push_back(Triplet<mpreal>(itI.row(), k, itI.value()));
                            } else {
                                tripletListR.push_back(Triplet<mpreal>(itI.row(), k, itRb.value()));
                            }
                            ++itRb;
                            ++itI;
                        }
                    } else {
                        if (mpfr::abs(itRb.value()) < 0)
                            tripletListR.push_back(Triplet<mpreal>(itRb.row(), k, itRb.value()));
                        ++itRb;
                    }
                }
            }
        } else {
            for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
                SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
                SparseMatrix<mpreal>::InnerIterator itRb(b.matrixR,k);

                while ((itR) || (itRb)) {
                    if ((!itRb) || ((itR) && (itR.row() < itRb.row()))) {
                        if (itR.value() < 0)
                            tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                        ++itR;
                    } else if ((itR) && (itRb) && (itR.row() == itRb.row())) {
                        if ((mpfr::isnan(itRb.value())) || (itR.value() <= itRb.value()))
                            tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                        else
                            tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itRb.value()));
                        ++itRb;
                        ++itR;
                    } else {
                        if (itRb.value() < 0)
                            tripletListR.push_back(Triplet<mpreal>(itRb.row(), k, itRb.value()));
                        ++itRb;
                    }
                }
            }
        }

        // Now we can assign the data to the sparse matrix
        result.matrixR.setFromTriplets(tripletListR.begin(), tripletListR.end());
        if ((isComplex) || (b.isComplex))
            result.matrixI.setFromTriplets(tripletListI.begin(), tripletListI.end());

        result.checkComplexity();
    }

    return result;
}

// column-wise maximum b = max(a)
SparseGmpEigenMatrix SparseGmpEigenMatrix::colMax(vector<IndexType>& indices) const
{
    SparseGmpEigenMatrix result;

    // We initialize the index vector
    indices.clear();
    indices.resize(matrixR.cols(), 0);

    if (isComplex) {
        // We initialize the result matrix
        result.isComplex = true;
        result.matrixR.resize(1,matrixR.cols());
        result.matrixI.resize(1,matrixI.cols());
        result.matrixR.reserve(min(matrixR.rows(),matrixR.nonZeros()));
        result.matrixI.reserve(min(matrixI.rows(),matrixI.nonZeros()));

        // We find the maximum
        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            mpreal maxValue;
            mpreal maxAngle;
            maxValue.setInf(-1);
            maxAngle.setInf(-1);
            IndexType firstZero(0); // This variable is incremented until there is a gap in column enumeration

            SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
            SparseMatrix<mpreal>::InnerIterator itI(matrixI,k);
            while ((itR) || (itI)) {
                if ((itR) && (itI)) {
                    if (itR.row() < itI.row()) {
                        if (firstZero == itR.row())
                            ++firstZero;
                        if (!mpfr::isnan(itR.value())) {
                            if ((pow(itR.value(),2) > maxValue) || ((pow(itR.value(),2) == maxValue) && (atan2(0, itR.value()) > maxAngle))) {
                                indices[k] = itR.row();
                                maxValue = pow(itR.value(),2);
                                maxAngle = atan2(0, itR.value());
                            }
                        }
                        ++itR;
                    } else if (itR.row() == itI.row()) {
                        if (firstZero == itR.row())
                            ++firstZero;
                        if ((!mpfr::isnan(itR.value())) && (!mpfr::isnan(itI.value()))) {
                            if (((pow(itR.value(),2) + pow(itI.value(),2)) > maxValue) || (((pow(itR.value(),2) + pow(itI.value(),2)) == maxValue) && (atan2(itI.value(), itR.value()) > maxAngle))) {
                                indices[k] = itR.row();
                                maxValue = pow(itR.value(),2) + pow(itI.value(),2);
                                maxAngle = atan2(itI.value(), itR.value());
                            }
                        }
                        ++itR;
                        ++itI;
                    } else {
                        if (firstZero == itI.row())
                            ++firstZero;
                        if (!mpfr::isnan(itI.value())) {
                            if ((pow(itI.value(),2) > maxValue) || ((pow(itI.value(),2) == maxValue) && (atan2(itI.value(), 0) > maxAngle))) {
                                indices[k] = itI.row();
                                maxValue = pow(itI.value(),2);
                                maxAngle = atan2(itI.value(), 0);
                            }
                        }
                        ++itI;
                    }
                } else if (itR) {
                    if (firstZero == itR.row())
                        ++firstZero;
                    if (!mpfr::isnan(itR.value())) {
                        if ((pow(itR.value(),2) > maxValue) || ((pow(itR.value(),2) == maxValue) && (atan2(0, itR.value()) > maxAngle))) {
                            indices[k] = itR.row();
                            maxValue = pow(itR.value(),2);
                            maxAngle = atan2(0, itR.value());
                        }
                    }
                    ++itR;
                } else {
                    if (firstZero == itI.row())
                        ++firstZero;
                    if (!mpfr::isnan(itI.value())) {
                        if ((pow(itI.value(),2) > maxValue) || ((pow(itI.value(),2) == maxValue) && (atan2(itI.value(), 0) > maxAngle))) {
                            indices[k] = itI.row();
                            maxValue = pow(itI.value(),2);
                            maxAngle = atan2(itI.value(), 0);
                        }
                    }
                    ++itI;
                }
            }

            if ((firstZero < matrixR.rows()) && (0 > maxValue)) {// If the column contains zeros
                indices[k] = firstZero;
                maxValue = 0;
            }
            if (maxValue != 0) {
                if (matrixR.coeff(indices[k],k) != 0)
                    result.matrixR.insert(0,k) = matrixR.coeff(indices[k],k);
                if (matrixI.coeff(indices[k],k) != 0)
                    result.matrixI.insert(0,k) = matrixI.coeff(indices[k],k);
            }
        }
        result.matrixR.makeCompressed();
        result.matrixI.makeCompressed();

        result.checkComplexity();
    } else {
        // We initialize the result matrix
        result.isComplex = false;
        result.matrixR.resize(1,matrixR.cols());
        result.matrixR.reserve(min(matrixR.rows(), matrixR.nonZeros()));

        // We find the maximum
        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            mpreal maxValue;
            maxValue.setInf(-1);
            IndexType firstZero(0); // This variable is incremented until there is a gap in column enumeration
            for (SparseMatrix<mpreal>::InnerIterator it(matrixR,k); it; ++it) {
                if (firstZero == it.row())
                    ++firstZero;
                if ((!mpfr::isnan(it.value())) && (it.value() > maxValue)) {
                    indices[k] = it.row();
                    maxValue = it.value();
                }
            }

            if ((firstZero < matrixR.rows()) && (0 > maxValue)) {// If the column contains zeros
                indices[k] = firstZero;
                maxValue = 0;
            }
            if (maxValue != 0)
                result.matrixR.insert(0,k) = maxValue;
        }
        result.matrixR.makeCompressed();
    }

    return result;
}

// column-wise maximum b = max(a)
SparseGmpEigenMatrix& SparseGmpEigenMatrix::colMax_new(vector<IndexType>& indices) const
{
    SparseGmpEigenMatrix& result(*(new SparseGmpEigenMatrix));

    // We initialize the index vector
    indices.clear();
    indices.resize(matrixR.cols(), 0);

    if (isComplex) {
        // We initialize the result matrix
        result.isComplex = true;
        result.matrixR.resize(1,matrixR.cols());
        result.matrixI.resize(1,matrixI.cols());
        result.matrixR.reserve(min(matrixR.rows(),matrixR.nonZeros()));
        result.matrixI.reserve(min(matrixI.rows(),matrixI.nonZeros()));

        // We find the maximum
        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            mpreal maxValue;
            mpreal maxAngle;
            maxValue.setInf(-1);
            maxAngle.setInf(-1);
            IndexType firstZero(0); // This variable is incremented until there is a gap in column enumeration

            SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
            SparseMatrix<mpreal>::InnerIterator itI(matrixI,k);
            while ((itR) || (itI)) {
                if ((itR) && (itI)) {
                    if (itR.row() < itI.row()) {
                        if (firstZero == itR.row())
                            ++firstZero;
                        if (!mpfr::isnan(itR.value())) {
                            if ((pow(itR.value(),2) > maxValue) || ((pow(itR.value(),2) == maxValue) && (atan2(0, itR.value()) > maxAngle))) {
                                indices[k] = itR.row();
                                maxValue = pow(itR.value(),2);
                                maxAngle = atan2(0, itR.value());
                            }
                        }
                        ++itR;
                    } else if (itR.row() == itI.row()) {
                        if (firstZero == itR.row())
                            ++firstZero;
                        if ((!mpfr::isnan(itR.value())) && (!mpfr::isnan(itI.value()))) {
                            if (((pow(itR.value(),2) + pow(itI.value(),2)) > maxValue) || (((pow(itR.value(),2) + pow(itI.value(),2)) == maxValue) && (atan2(itI.value(), itR.value()) > maxAngle))) {
                                indices[k] = itR.row();
                                maxValue = pow(itR.value(),2) + pow(itI.value(),2);
                                maxAngle = atan2(itI.value(), itR.value());
                            }
                        }
                        ++itR;
                        ++itI;
                    } else {
                        if (firstZero == itI.row())
                            ++firstZero;
                        if (!mpfr::isnan(itI.value())) {
                            if ((pow(itI.value(),2) > maxValue) || ((pow(itI.value(),2) == maxValue) && (atan2(itI.value(), 0) > maxAngle))) {
                                indices[k] = itI.row();
                                maxValue = pow(itI.value(),2);
                                maxAngle = atan2(itI.value(), 0);
                            }
                        }
                        ++itI;
                    }
                } else if (itR) {
                    if (firstZero == itR.row())
                        ++firstZero;
                    if (!mpfr::isnan(itR.value())) {
                        if ((pow(itR.value(),2) > maxValue) || ((pow(itR.value(),2) == maxValue) && (atan2(0, itR.value()) > maxAngle))) {
                            indices[k] = itR.row();
                            maxValue = pow(itR.value(),2);
                            maxAngle = atan2(0, itR.value());
                        }
                    }
                    ++itR;
                } else {
                    if (firstZero == itI.row())
                        ++firstZero;
                    if (!mpfr::isnan(itI.value())) {
                        if ((pow(itI.value(),2) > maxValue) || ((pow(itI.value(),2) == maxValue) && (atan2(itI.value(), 0) > maxAngle))) {
                            indices[k] = itI.row();
                            maxValue = pow(itI.value(),2);
                            maxAngle = atan2(itI.value(), 0);
                        }
                    }
                    ++itI;
                }
            }

            if ((firstZero < matrixR.rows()) && (0 > maxValue)) {// If the column contains zeros
                indices[k] = firstZero;
                maxValue = 0;
            }
            if (maxValue != 0) {
                if (matrixR.coeff(indices[k],k) != 0)
                    result.matrixR.insert(0,k) = matrixR.coeff(indices[k],k);
                if (matrixI.coeff(indices[k],k) != 0)
                    result.matrixI.insert(0,k) = matrixI.coeff(indices[k],k);
            }
        }
        result.matrixR.makeCompressed();
        result.matrixI.makeCompressed();

        result.checkComplexity();
    } else {
        // We initialize the result matrix
        result.isComplex = false;
        result.matrixR.resize(1,matrixR.cols());
        result.matrixR.reserve(min(matrixR.rows(), matrixR.nonZeros()));

        // We find the maximum
        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            mpreal maxValue;
            maxValue.setInf(-1);
            IndexType firstZero(0); // This variable is incremented until there is a gap in column enumeration
            for (SparseMatrix<mpreal>::InnerIterator it(matrixR,k); it; ++it) {
                if (firstZero == it.row())
                    ++firstZero;
                if ((!mpfr::isnan(it.value())) && (it.value() > maxValue)) {
                    indices[k] = it.row();
                    maxValue = it.value();
                }
            }

            if ((firstZero < matrixR.rows()) && (0 > maxValue)) {// If the column contains zeros
                indices[k] = firstZero;
                maxValue = 0;
            }
            if (maxValue != 0)
                result.matrixR.insert(0,k) = maxValue;
        }
        result.matrixR.makeCompressed();
    }

    return result;
}

// line-wise maximum b = max(a,[],2)
SparseGmpEigenMatrix SparseGmpEigenMatrix::rowMax(vector<IndexType>& indices) const
{
    SparseGmpEigenMatrix result;

    // We initialize the index vector
    indices.clear();
    indices.resize(matrixR.rows(), 0);

    if (isComplex) {
        // We initialize the result matrix
        result.isComplex = true;
        result.matrixR.resize(matrixR.rows(),1);
        result.matrixI.resize(matrixI.rows(),1);
        result.matrixR.reserve(min(matrixR.cols(),matrixR.nonZeros()));
        result.matrixI.reserve(min(matrixI.cols(),matrixI.nonZeros()));
        vector < mpreal > maxValue(matrixR.rows(),const_infinity(-1));
        vector < mpreal > maxAngle(matrixR.rows(),const_infinity(-1));

        // We find the minimums
        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            IndexType previousNonZero(-1);

            SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
            SparseMatrix<mpreal>::InnerIterator itI(matrixI,k);

            // If there is nothing in this column
            if ((!itR) && (!itI)) {
                for (IndexType i(0); i < matrixR.rows(); ++i) {
                    if (0 > maxValue[i]) {
                        indices[i] = k;
                        maxValue[i] = 0;
                    }
                }
            }

            while ((itR) || (itI)) {
                if ((itR) && (itI)) {
                    if (itR.row() < itI.row()) {
                        for (IndexType i(previousNonZero+1); i < itR.row(); ++i) {
                            if (0 > maxValue[i]) {
                                indices[i] = k;
                                maxValue[i] = 0;
                            }
                        }
                        previousNonZero = itR.row();
                        if (!mpfr::isnan(itR.value())) {
                            if ((pow(itR.value(),2) > maxValue[itR.row()]) || ((pow(itR.value(),2) == maxValue[itR.row()]) && (atan2(0, itR.value()) > maxAngle[itR.row()]))) {
                                indices[itR.row()] = k;
                                maxValue[itR.row()] = pow(itR.value(),2);
                                maxAngle[itR.row()] = atan2(0, itR.value());
                            }
                        }
                        ++itR;
                    } else if (itR.row() == itI.row()) {
                        for (IndexType i(previousNonZero+1); i < itR.row(); ++i) {
                            if (0 > maxValue[i]) {
                                indices[i] = k;
                                maxValue[i] = 0;
                            }
                        }
                        previousNonZero = itR.row();
                        if ((!mpfr::isnan(itR.value())) && (!mpfr::isnan(itI.value()))) {
                            if (((pow(itR.value(),2) + pow(itI.value(),2)) > maxValue[itR.row()]) || (((pow(itR.value(),2) + pow(itI.value(),2)) == maxValue[itR.row()]) && (atan2(itI.value(), itR.value()) > maxAngle[itR.row()]))) {
                                indices[itR.row()] = k;
                                maxValue[itR.row()] = pow(itR.value(),2) + pow(itI.value(),2);
                                maxAngle[itR.row()] = atan2(itI.value(), itR.value());
                            }
                        }
                        ++itR;
                        ++itI;
                    } else {
                        for (IndexType i(previousNonZero+1); i < itI.row(); ++i) {
                            if (0 > maxValue[i]) {
                                indices[i] = k;
                                maxValue[i] = 0;
                            }
                        }
                        previousNonZero = itR.row();
                        if (!mpfr::isnan(itI.value())) {
                            if ((pow(itI.value(),2) > maxValue[itI.row()]) || ((pow(itI.value(),2) == maxValue[itI.row()]) && (atan2(itI.value(), 0) > maxAngle[itI.row()]))) {
                                indices[itI.row()] = k;
                                maxValue[itI.row()] = pow(itI.value(),2);
                                maxAngle[itI.row()] = atan2(itI.value(), 0);
                            }
                        }
                        ++itI;
                    }
                } else if (itR) {
                    for (IndexType i(previousNonZero+1); i < itR.row(); ++i) {
                        if (0 > maxValue[i]) {
                            indices[i] = k;
                            maxValue[i] = 0;
                        }
                    }
                    previousNonZero = itR.row();
                    if (!mpfr::isnan(itR.value())) {
                        if ((pow(itR.value(),2) > maxValue[itR.row()]) || ((pow(itR.value(),2) == maxValue[itR.row()]) && (atan2(0, itR.value()) > maxAngle[itR.row()]))) {
                            indices[itR.row()] = k;
                            maxValue[itR.row()] = pow(itR.value(),2);
                            maxAngle[itR.row()] = atan2(0, itR.value());
                        }
                    }
                    ++itR;
                } else {
                    for (IndexType i(previousNonZero+1); i < itI.row(); ++i) {
                        if (0 > maxValue[i]) {
                            indices[i] = k;
                            maxValue[i] = 0;
                        }
                    }
                    previousNonZero = itI.row();
                    if (!mpfr::isnan(itI.value())) {
                        if ((pow(itI.value(),2) > maxValue[itI.row()]) || ((pow(itI.value(),2) == maxValue[itI.row()]) && (atan2(itI.value(), 0) > maxAngle[itI.row()]))) {
                            indices[itI.row()] = k;
                            maxValue[itI.row()] = pow(itI.value(),2);
                            maxAngle[itI.row()] = atan2(itI.value(), 0);
                        }
                    }
                    ++itI;
                }
            }

            // If there are still remaining zeros
            for (IndexType i(previousNonZero+1); i < matrixR.rows(); ++i) {
                if (0 > maxValue[i]) {
                    indices[i] = k;
                    maxValue[i] = 0;
                }
            }
        }

        // We copy the minimums
        for (IndexType i = 0; i < matrixR.rows(); ++i) {
            if (maxValue[i] != 0) {
                if (matrixR.coeff(i,indices[i]) != 0)
                    result.matrixR.insert(i,0) = matrixR.coeff(i,indices[i]);
                if (matrixI.coeff(i,indices[i]) != 0)
                    result.matrixI.insert(i,0) = matrixI.coeff(i,indices[i]);
            }
        }
        result.matrixR.makeCompressed();
        result.matrixI.makeCompressed();

        result.checkComplexity();
    } else {
        // We initialize the result matrix
        result.isComplex = false;
        result.matrixR.resize(matrixR.rows(),1);
        result.matrixR.reserve(min(matrixR.cols(),matrixR.nonZeros()));
        vector < mpreal > maxValue(matrixR.rows(),const_infinity(-1));

        // We find the minimums
        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            IndexType previousNonZero(-1);
            SparseMatrix<mpreal>::InnerIterator it(matrixR,k);

            // If there is nothing in this column
            if (!it) {
                for (IndexType i(0); i < matrixR.rows(); ++i) {
                    if (0 > maxValue[i]) {
                        indices[i] = k;
                        maxValue[i] = 0;
                    }
                }
            }

            for (; it; ++it) {
                for (IndexType i(previousNonZero+1); i < it.row(); ++i) {
                    if (0 > maxValue[i]) {
                        indices[i] = k;
                        maxValue[i] = 0;
                    }
                }
                previousNonZero = it.row();
                if ((!mpfr::isnan(it.value())) && (it.value() > maxValue[it.row()])) {
                    indices[it.row()] = k;
                    maxValue[it.row()] = it.value();
                }
            }

            // If there are still zeros
            for (IndexType i(previousNonZero+1); i < matrixR.rows(); ++i) {
                if (0 > maxValue[i]) {
                    indices[i] = k;
                    maxValue[i] = 0;
                }
            }
        }

        // We copy the minimums
        for (IndexType i = 0; i < matrixR.rows(); ++i) {
            if (maxValue[i] != 0) {
                result.matrixR.insert(i,0) = maxValue[i];
            }
        }
        result.matrixR.makeCompressed();
    }

    return result;
}

// line-wise maximum b = max(a,[],2)
SparseGmpEigenMatrix& SparseGmpEigenMatrix::rowMax_new(vector<IndexType>& indices) const
{
    SparseGmpEigenMatrix& result(*(new SparseGmpEigenMatrix));

    // We initialize the index vector
    indices.clear();
    indices.resize(matrixR.rows(), 0);

    if (isComplex) {
        // We initialize the result matrix
        result.isComplex = true;
        result.matrixR.resize(matrixR.rows(),1);
        result.matrixI.resize(matrixI.rows(),1);
        result.matrixR.reserve(min(matrixR.cols(),matrixR.nonZeros()));
        result.matrixI.reserve(min(matrixI.cols(),matrixI.nonZeros()));
        vector < mpreal > maxValue(matrixR.rows(),const_infinity(-1));
        vector < mpreal > maxAngle(matrixR.rows(),const_infinity(-1));

        // We find the minimums
        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            IndexType previousNonZero(-1);

            SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
            SparseMatrix<mpreal>::InnerIterator itI(matrixI,k);

            // If there is nothing in this column
            if ((!itR) && (!itI)) {
                for (IndexType i(0); i < matrixR.rows(); ++i) {
                    if (0 > maxValue[i]) {
                        indices[i] = k;
                        maxValue[i] = 0;
                    }
                }
            }

            while ((itR) || (itI)) {
                if ((itR) && (itI)) {
                    if (itR.row() < itI.row()) {
                        for (IndexType i(previousNonZero+1); i < itR.row(); ++i) {
                            if (0 > maxValue[i]) {
                                indices[i] = k;
                                maxValue[i] = 0;
                            }
                        }
                        previousNonZero = itR.row();
                        if (!mpfr::isnan(itR.value())) {
                            if ((pow(itR.value(),2) > maxValue[itR.row()]) || ((pow(itR.value(),2) == maxValue[itR.row()]) && (atan2(0, itR.value()) > maxAngle[itR.row()]))) {
                                indices[itR.row()] = k;
                                maxValue[itR.row()] = pow(itR.value(),2);
                                maxAngle[itR.row()] = atan2(0, itR.value());
                            }
                        }
                        ++itR;
                    } else if (itR.row() == itI.row()) {
                        for (IndexType i(previousNonZero+1); i < itR.row(); ++i) {
                            if (0 > maxValue[i]) {
                                indices[i] = k;
                                maxValue[i] = 0;
                            }
                        }
                        previousNonZero = itR.row();
                        if ((!mpfr::isnan(itR.value())) && (!mpfr::isnan(itI.value()))) {
                            if (((pow(itR.value(),2) + pow(itI.value(),2)) > maxValue[itR.row()]) || (((pow(itR.value(),2) + pow(itI.value(),2)) == maxValue[itR.row()]) && (atan2(itI.value(), itR.value()) > maxAngle[itR.row()]))) {
                                indices[itR.row()] = k;
                                maxValue[itR.row()] = pow(itR.value(),2) + pow(itI.value(),2);
                                maxAngle[itR.row()] = atan2(itI.value(), itR.value());
                            }
                        }
                        ++itR;
                        ++itI;
                    } else {
                        for (IndexType i(previousNonZero+1); i < itI.row(); ++i) {
                            if (0 > maxValue[i]) {
                                indices[i] = k;
                                maxValue[i] = 0;
                            }
                        }
                        previousNonZero = itR.row();
                        if (!mpfr::isnan(itI.value())) {
                            if ((pow(itI.value(),2) > maxValue[itI.row()]) || ((pow(itI.value(),2) == maxValue[itI.row()]) && (atan2(itI.value(), 0) > maxAngle[itI.row()]))) {
                                indices[itI.row()] = k;
                                maxValue[itI.row()] = pow(itI.value(),2);
                                maxAngle[itI.row()] = atan2(itI.value(), 0);
                            }
                        }
                        ++itI;
                    }
                } else if (itR) {
                    for (IndexType i(previousNonZero+1); i < itR.row(); ++i) {
                        if (0 > maxValue[i]) {
                            indices[i] = k;
                            maxValue[i] = 0;
                        }
                    }
                    previousNonZero = itR.row();
                    if (!mpfr::isnan(itR.value())) {
                        if ((pow(itR.value(),2) > maxValue[itR.row()]) || ((pow(itR.value(),2) == maxValue[itR.row()]) && (atan2(0, itR.value()) > maxAngle[itR.row()]))) {
                            indices[itR.row()] = k;
                            maxValue[itR.row()] = pow(itR.value(),2);
                            maxAngle[itR.row()] = atan2(0, itR.value());
                        }
                    }
                    ++itR;
                } else {
                    for (IndexType i(previousNonZero+1); i < itI.row(); ++i) {
                        if (0 > maxValue[i]) {
                            indices[i] = k;
                            maxValue[i] = 0;
                        }
                    }
                    previousNonZero = itI.row();
                    if (!mpfr::isnan(itI.value())) {
                        if ((pow(itI.value(),2) > maxValue[itI.row()]) || ((pow(itI.value(),2) == maxValue[itI.row()]) && (atan2(itI.value(), 0) > maxAngle[itI.row()]))) {
                            indices[itI.row()] = k;
                            maxValue[itI.row()] = pow(itI.value(),2);
                            maxAngle[itI.row()] = atan2(itI.value(), 0);
                        }
                    }
                    ++itI;
                }
            }

            // If there are still remaining zeros
            for (IndexType i(previousNonZero+1); i < matrixR.rows(); ++i) {
                if (0 > maxValue[i]) {
                    indices[i] = k;
                    maxValue[i] = 0;
                }
            }
        }

        // We copy the minimums
        for (IndexType i = 0; i < matrixR.rows(); ++i) {
            if (maxValue[i] != 0) {
                if (matrixR.coeff(i,indices[i]) != 0)
                    result.matrixR.insert(i,0) = matrixR.coeff(i,indices[i]);
                if (matrixI.coeff(i,indices[i]) != 0)
                    result.matrixI.insert(i,0) = matrixI.coeff(i,indices[i]);
            }
        }
        result.matrixR.makeCompressed();
        result.matrixI.makeCompressed();

        result.checkComplexity();
    } else {
        // We initialize the result matrix
        result.isComplex = false;
        result.matrixR.resize(matrixR.rows(),1);
        result.matrixR.reserve(min(matrixR.cols(),matrixR.nonZeros()));
        vector < mpreal > maxValue(matrixR.rows(),const_infinity(-1));

        // We find the minimums
        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            IndexType previousNonZero(-1);
            SparseMatrix<mpreal>::InnerIterator it(matrixR,k);

            // If there is nothing in this column
            if (!it) {
                for (IndexType i(0); i < matrixR.rows(); ++i) {
                    if (0 > maxValue[i]) {
                        indices[i] = k;
                        maxValue[i] = 0;
                    }
                }
            }

            for (; it; ++it) {
                for (IndexType i(previousNonZero+1); i < it.row(); ++i) {
                    if (0 > maxValue[i]) {
                        indices[i] = k;
                        maxValue[i] = 0;
                    }
                }
                previousNonZero = it.row();
                if ((!mpfr::isnan(it.value())) && (it.value() > maxValue[it.row()])) {
                    indices[it.row()] = k;
                    maxValue[it.row()] = it.value();
                }
            }

            // If there are still zeros
            for (IndexType i(previousNonZero+1); i < matrixR.rows(); ++i) {
                if (0 > maxValue[i]) {
                    indices[i] = k;
                    maxValue[i] = 0;
                }
            }
        }

        // We copy the minimums
        for (IndexType i = 0; i < matrixR.rows(); ++i) {
            if (maxValue[i] != 0) {
                result.matrixR.insert(i,0) = maxValue[i];
            }
        }
        result.matrixR.makeCompressed();
    }

    return result;
}

// element-wise maximum c = max(a, b)
SparseGmpEigenMatrix SparseGmpEigenMatrix::ewMax(const SparseGmpEigenMatrix& b) const
{
    SparseGmpEigenMatrix result;

    if ((numel() != 1) && (b.numel() == 1)) {
        // WARNING : Here, we assume that the scalar b is not larger than 0, and that no matrix is complex (otherwise it will make the sparse matrix full...)

        // setting the output size
        result.matrixR.resize(matrixR.rows(), matrixR.cols());
        result.isComplex = false;
        if ((isComplex) || (b.isComplex))
            cout << "Error: unexpected situation. Expect unexpected result." << endl;

        // We will first copy the data to triplets
        vector< Triplet<mpreal> > tripletListR;
        tripletListR.reserve(matrixR.nonZeros());

        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);

            while (itR) {
                if ((mpfr::isnan(b.matrixR.coeff(0,0))) || (itR.value() >= b.matrixR.coeff(0,0)))
                    tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                else
                    tripletListR.push_back(Triplet<mpreal>(itR.row(), k, b.matrixR.coeff(0,0)));
                ++itR;
            }
        }

        // Now we can assign the data to the sparse matrix
        result.matrixR.setFromTriplets(tripletListR.begin(), tripletListR.end());
    } else if ((numel() == 1) && (b.numel() != 1)) {
        result = b.ewMax(*this);
    } else {
        // setting the output size
        result.matrixR.resize(matrixR.rows(), matrixR.cols());
        if ((isComplex) || (b.isComplex))
            result.matrixI.resize(matrixR.rows(), matrixR.cols());

        // We will first copy the data to triplets
        vector< Triplet<mpreal> > tripletListR, tripletListI;
        tripletListR.reserve(matrixR.nonZeros() + b.matrixR.nonZeros());
        if ((isComplex) || (b.isComplex))
            tripletListI.reserve(matrixI.nonZeros() + b.matrixI.nonZeros());

        // Now for each column, we merge the lists of lines with nonzero elements of
        // both matrices. The cases are slightly different depending on the complexity
        // of both matrices
        if (b.isComplex) {
            if (isComplex) {
                for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
                    SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
                    SparseMatrix<mpreal>::InnerIterator itI(matrixI,k);
                    SparseMatrix<mpreal>::InnerIterator itRb(b.matrixR,k);
                    SparseMatrix<mpreal>::InnerIterator itIb(b.matrixI,k);

                    while ((itR) || (itI) || (itRb) || (itIb)) {
                        IndexType row, rowb;
                        if ((itR) || (itI)) {
                            if ((itR) && (itI))
                                row = min(itR.row(), itI.row());
                            else if (itR)
                                row = itR.row();
                            else
                                row = itI.row();
                        }
                        if ((itRb) || (itIb)) {
                            if ((itRb) && (itIb))
                                rowb = min(itRb.row(), itIb.row());
                            else if (itRb)
                                rowb = itRb.row();
                            else
                                rowb = itIb.row();
                        }

                        if (((((itR) || (itI)) && ((itRb) || (itIb))) && (row < rowb)) || ((!itRb) && (!itIb))) {
                            if ((!itI) || ((itR) && (itR.row() < itI.row()))) {
                                if (mpfr::abs(itR.value()) > 0) // Note : this will be wrong if itR.value() is NaN
                                    tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                                ++itR;
                            } else if ((itR) && (itI) && (itR.row() == itI.row())) {
                                if (pow(itR.value(), 2) + pow(itI.value(), 2) > 0) {
                                    tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                                    tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itI.value()));
                                }
                                ++itR;
                                ++itI;
                            } else {
                                if (mpfr::abs(itI.value()) > 0)
                                    tripletListI.push_back(Triplet<mpreal>(itI.row(), k, itI.value()));
                                ++itI;
                            }
                        } else if ((((itR) || (itI)) && ((itRb) || (itIb))) && (row == rowb)) {
                            if ((!itI) || ((itR) && (itR.row() < itI.row()))) {
                                if ((!itIb) || ((itRb) && (itRb.row() < itIb.row()))) {
                                    if ((mpfr::isnan(itRb.value())) || (mpfr::abs(itR.value()) >= mpfr::abs(itRb.value())))
                                        tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                                    else
                                        tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itRb.value()));
                                    ++itRb;
                                } else if ((itRb) && (itIb) && (itRb.row() == itIb.row())) {
                                    if ((mpfr::isnan(pow(itRb.value(), 2) + pow(itIb.value(), 2))) || ((pow(itR.value(), 2) > pow(itRb.value(), 2) + pow(itIb.value(), 2)) || ((pow(itR.value(), 2) == pow(itRb.value(), 2) + pow(itIb.value(), 2)) && (atan2(0, itR.value()) >= atan2(itIb.value(), itRb.value()))))) {
                                        tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                                    } else {
                                        tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itRb.value()));
                                        tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itIb.value()));
                                    }
                                    ++itRb;
                                    ++itIb;
                                } else {
                                    if ((mpfr::isnan(itIb.value())) || ((pow(itR.value(), 2) > pow(itIb.value(), 2)) || ((pow(itR.value(), 2) == pow(itIb.value(), 2)) && (atan2(0, itR.value()) >= atan2(itIb.value(), 0)))))
                                        tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                                    else
                                        tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itIb.value()));
                                    ++itIb;
                                }
                                ++itR;
                            } else if ((itR) && (itI) && (itR.row() == itI.row())) {
                                if ((!itIb) || ((itRb) && (itRb.row() < itIb.row()))) {
                                    if ((mpfr::isnan(itRb.value())) || ((pow(itR.value(), 2) + pow(itI.value(), 2) > pow(itRb.value(), 2)) || ((pow(itR.value(), 2) + pow(itI.value(), 2) == pow(itRb.value(), 2)) && (atan2(itI.value(), itR.value()) >= atan2(0, itRb.value()))))) {
                                        tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                                        tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itI.value()));
                                    } else {
                                        tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itRb.value()));
                                    }
                                    ++itRb;
                                } else if ((itRb) && (itIb) && (itRb.row() == itIb.row())) {
                                    if ((mpfr::isnan(pow(itRb.value(), 2) + pow(itIb.value(), 2))) || ((pow(itR.value(), 2) + pow(itI.value(), 2) > pow(itRb.value(), 2) + pow(itIb.value(), 2)) || ((pow(itR.value(), 2) + pow(itI.value(), 2) == pow(itRb.value(), 2) + pow(itIb.value(), 2)) && (atan2(itI.value(), itR.value()) >= atan2(itIb.value(), itRb.value()))))) {
                                        tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                                        tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itI.value()));
                                    } else {
                                        tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itRb.value()));
                                        tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itIb.value()));
                                    }
                                    ++itRb;
                                    ++itIb;
                                } else {
                                    if ((mpfr::isnan(itIb.value())) || ((pow(itR.value(), 2) + pow(itI.value(), 2) > pow(itIb.value(), 2)) || ((pow(itR.value(), 2) + pow(itI.value(), 2) == pow(itIb.value(), 2)) && (atan2(itI.value(), itR.value()) >= atan2(itIb.value(), 0))))) {
                                        tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                                        tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itI.value()));
                                    } else {
                                        tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itIb.value()));
                                    }
                                    ++itIb;
                                }
                                ++itR;
                                ++itI;
                            } else {
                                if ((!itIb) || ((itRb) && (itRb.row() < itIb.row()))) {
                                    if ((mpfr::isnan(itRb.value())) || ((pow(itI.value(), 2) > pow(itRb.value(), 2)) || ((pow(itI.value(), 2) == pow(itRb.value(), 2)) && (atan2(itI.value(), 0) >= atan2(0, itRb.value()))))) {
                                        tripletListI.push_back(Triplet<mpreal>(itI.row(), k, itI.value()));
                                    } else {
                                        tripletListR.push_back(Triplet<mpreal>(itI.row(), k, itRb.value()));
                                    }
                                    ++itRb;
                                } else if ((itRb) && (itIb) && (itRb.row() == itIb.row())) {
                                    if ((mpfr::isnan(pow(itRb.value(), 2) + pow(itIb.value(), 2))) || ((pow(itI.value(), 2) > pow(itRb.value(), 2) + pow(itIb.value(), 2)) || ((pow(itI.value(), 2) == pow(itRb.value(), 2) + pow(itIb.value(), 2)) && (atan2(itI.value(), 0) >= atan2(itIb.value(), itRb.value()))))) {
                                        tripletListI.push_back(Triplet<mpreal>(itI.row(), k, itI.value()));
                                    } else {
                                        tripletListR.push_back(Triplet<mpreal>(itI.row(), k, itRb.value()));
                                        tripletListI.push_back(Triplet<mpreal>(itI.row(), k, itIb.value()));
                                    }
                                    ++itRb;
                                    ++itIb;
                                } else {
                                    if ((mpfr::isnan(itIb.value())) || ((pow(itI.value(), 2) > pow(itIb.value(), 2)) || ((pow(itI.value(), 2) == pow(itIb.value(), 2)) && (atan2(itI.value(), 0) >= atan2(itIb.value(), 0))))) {
                                        tripletListI.push_back(Triplet<mpreal>(itI.row(), k, itI.value()));
                                    } else {
                                        tripletListI.push_back(Triplet<mpreal>(itI.row(), k, itIb.value()));
                                    }
                                    ++itIb;
                                }
                                ++itI;
                            }
                        } else {
                            if ((!itIb) || ((itRb) && (itRb.row() < itIb.row()))) {
                                if (mpfr::abs(itRb.value()) > 0)
                                    tripletListR.push_back(Triplet<mpreal>(itRb.row(), k, itRb.value()));
                                ++itRb;
                            } else if ((itRb) && (itIb) && (itRb.row() == itIb.row())) {
                                if (pow(itRb.value(), 2) + pow(itIb.value(), 2) > 0) {
                                    tripletListR.push_back(Triplet<mpreal>(itRb.row(), k, itRb.value()));
                                    tripletListI.push_back(Triplet<mpreal>(itRb.row(), k, itIb.value()));
                                }
                                ++itRb;
                                ++itIb;
                            } else {
                                if (mpfr::abs(itIb.value()) > 0)
                                    tripletListI.push_back(Triplet<mpreal>(itIb.row(), k, itIb.value()));
                                ++itIb;
                            }
                        }
                    }
                }
            } else {
                for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
                    SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
                    SparseMatrix<mpreal>::InnerIterator itRb(b.matrixR,k);
                    SparseMatrix<mpreal>::InnerIterator itIb(b.matrixI,k);

                    while ((itR) || (itRb) || (itIb)) {
                        IndexType row, rowb;
                        if (itR)
                            row = itR.row();
                        if ((itRb) || (itIb)) {
                            if ((itRb) && (itIb))
                                rowb = min(itRb.row(), itIb.row());
                            else if (itRb)
                                rowb = itRb.row();
                            else
                                rowb = itIb.row();
                        }

                        if ((((itR) && ((itRb) || (itIb))) && (row < rowb)) || ((!itRb) && (!itIb))) {
                            if (mpfr::abs(itR.value()) > 0)
                                tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                            ++itR;
                        } else if (((itR) && ((itRb) || (itIb))) && (row == rowb)) {
                            if ((!itIb) || ((itRb) && (itRb.row() < itIb.row()))) {
                                if ((mpfr::isnan(itRb.value())) || (mpfr::abs(itR.value()) >= mpfr::abs(itRb.value())))
                                    tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                                else
                                    tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itRb.value()));
                                ++itRb;
                            } else if ((itRb) && (itIb) && (itRb.row() == itIb.row())) {
                                if ((mpfr::isnan(pow(itRb.value(), 2) + pow(itIb.value(), 2))) || ((pow(itR.value(), 2) > pow(itRb.value(), 2) + pow(itIb.value(), 2)) || ((pow(itR.value(), 2) == pow(itRb.value(), 2) + pow(itIb.value(), 2)) && (atan2(0, itR.value()) >= atan2(itIb.value(), itRb.value()))))) {
                                    tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                                } else {
                                    tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itRb.value()));
                                    tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itIb.value()));
                                }
                                ++itRb;
                                ++itIb;
                            } else {
                                if ((mpfr::isnan(itIb.value())) || ((pow(itR.value(), 2) > pow(itIb.value(), 2)) || ((pow(itR.value(), 2) == pow(itIb.value(), 2)) && (atan2(0, itR.value()) >= atan2(itIb.value(), 0)))))
                                    tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                                else
                                    tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itIb.value()));
                                ++itIb;
                            }
                            ++itR;
                        } else {
                            if ((!itIb) || ((itRb) && (itRb.row() < itIb.row()))) {
                                if (mpfr::abs(itRb.value()) > 0)
                                    tripletListR.push_back(Triplet<mpreal>(itRb.row(), k, itRb.value()));
                                ++itRb;
                            } else if ((itRb) && (itIb) && (itRb.row() == itIb.row())) {
                                if (pow(itRb.value(), 2) + pow(itIb.value(), 2) > 0) {
                                    tripletListR.push_back(Triplet<mpreal>(itRb.row(), k, itRb.value()));
                                    tripletListI.push_back(Triplet<mpreal>(itRb.row(), k, itIb.value()));
                                }
                                ++itRb;
                                ++itIb;
                            } else {
                                if (mpfr::abs(itIb.value()) > 0)
                                    tripletListI.push_back(Triplet<mpreal>(itIb.row(), k, itIb.value()));
                                ++itIb;
                            }
                        }
                    }
                }
            }
        } else if (isComplex) {
            for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
                SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
                SparseMatrix<mpreal>::InnerIterator itI(matrixI,k);
                SparseMatrix<mpreal>::InnerIterator itRb(b.matrixR,k);

                while ((itR) || (itI) || (itRb)) {
                    IndexType row, rowb;
                    if ((itR) || (itI)) {
                        if ((itR) && (itI))
                            row = min(itR.row(), itI.row());
                        else if (itR)
                            row = itR.row();
                        else
                            row = itI.row();
                    }
                    if (itRb)
                        rowb = itRb.row();

                    if (((((itR) || (itI)) && (itRb)) && (row < rowb)) || (!itRb)) {
                        if ((!itI) || ((itR) && (itR.row() < itI.row()))) {
                            if (mpfr::abs(itR.value()) > 0)
                                tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                            ++itR;
                        } else if ((itR) && (itI) && (itR.row() == itI.row())) {
                            if (pow(itR.value(), 2) + pow(itI.value(), 2) > 0) {
                                tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                                tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itI.value()));
                            }
                            ++itR;
                            ++itI;
                        } else {
                            if (mpfr::abs(itI.value()) > 0)
                                tripletListI.push_back(Triplet<mpreal>(itI.row(), k, itI.value()));
                            ++itI;
                        }
                    } else if ((((itR) || (itI)) && (itRb)) && (row == rowb)) {
                        if ((!itI) || ((itR) && (itR.row() < itI.row()))) {
                            if ((mpfr::isnan(itRb.value())) || (mpfr::abs(itR.value()) >= mpfr::abs(itRb.value())))
                                tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                            else
                                tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itRb.value()));
                            ++itRb;
                            ++itR;
                        } else if ((itR) && (itI) && (itR.row() == itI.row())) {
                            if ((mpfr::isnan(itRb.value())) || ((pow(itR.value(), 2) + pow(itI.value(), 2) > pow(itRb.value(), 2)) || ((pow(itR.value(), 2) + pow(itI.value(), 2) == pow(itRb.value(), 2)) && (atan2(itI.value(), itR.value()) >= atan2(0, itRb.value()))))) {
                                tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                                tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itI.value()));
                            } else {
                                tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itRb.value()));
                            }
                            ++itRb;
                            ++itR;
                            ++itI;
                        } else {
                            if ((mpfr::isnan(itRb.value())) || ((pow(itI.value(), 2) > pow(itRb.value(), 2)) || ((pow(itI.value(), 2) == pow(itRb.value(), 2)) && (atan2(itI.value(), 0) >= atan2(0, itRb.value()))))) {
                                tripletListI.push_back(Triplet<mpreal>(itI.row(), k, itI.value()));
                            } else {
                                tripletListR.push_back(Triplet<mpreal>(itI.row(), k, itRb.value()));
                            }
                            ++itRb;
                            ++itI;
                        }
                    } else {
                        if (mpfr::abs(itRb.value()) > 0)
                            tripletListR.push_back(Triplet<mpreal>(itRb.row(), k, itRb.value()));
                        ++itRb;
                    }
                }
            }
        } else {
            for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
                SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
                SparseMatrix<mpreal>::InnerIterator itRb(b.matrixR,k);

                while ((itR) || (itRb)) {
                    if ((!itRb) || ((itR) && (itR.row() < itRb.row()))) {
                        if (itR.value() > 0)
                            tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                        ++itR;
                    } else if ((itR) && (itRb) && (itR.row() == itRb.row())) {
                        if ((mpfr::isnan(itRb.value())) || (itR.value() >= itRb.value()))
                            tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                        else
                            tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itRb.value()));
                        ++itRb;
                        ++itR;
                    } else {
                        if (itRb.value() > 0)
                            tripletListR.push_back(Triplet<mpreal>(itRb.row(), k, itRb.value()));
                        ++itRb;
                    }
                }
            }
        }

        // Now we can assign the data to the sparse matrix
        result.matrixR.setFromTriplets(tripletListR.begin(), tripletListR.end());
        if ((isComplex) || (b.isComplex))
            result.matrixI.setFromTriplets(tripletListI.begin(), tripletListI.end());

        result.checkComplexity();
    }

    return result;
}

// element-wise maximum c = max(a, b)
SparseGmpEigenMatrix& SparseGmpEigenMatrix::ewMax_new(const SparseGmpEigenMatrix& b) const
{
    SparseGmpEigenMatrix& result(*(new SparseGmpEigenMatrix));

    if ((numel() != 1) && (b.numel() == 1)) {
        // WARNING : Here, we assume that the scalar b is not larger than 0, and that no matrix is complex (otherwise it will make the sparse matrix full...)

        // setting the output size
        result.matrixR.resize(matrixR.rows(), matrixR.cols());
        result.isComplex = false;
        if ((isComplex) || (b.isComplex))
            cout << "Error: unexpected situation. Expect unexpected result." << endl;

        // We will first copy the data to triplets
        vector< Triplet<mpreal> > tripletListR;
        tripletListR.reserve(matrixR.nonZeros());

        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);

            while (itR) {
                if ((mpfr::isnan(b.matrixR.coeff(0,0))) || (itR.value() >= b.matrixR.coeff(0,0)))
                    tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                else
                    tripletListR.push_back(Triplet<mpreal>(itR.row(), k, b.matrixR.coeff(0,0)));
                ++itR;
            }
        }

        // Now we can assign the data to the sparse matrix
        result.matrixR.setFromTriplets(tripletListR.begin(), tripletListR.end());
    } else if ((numel() == 1) && (b.numel() != 1)) {
        result = b.ewMax(*this);
    } else {
        // setting the output size
        result.matrixR.resize(matrixR.rows(), matrixR.cols());
        if ((isComplex) || (b.isComplex))
            result.matrixI.resize(matrixR.rows(), matrixR.cols());

        // We will first copy the data to triplets
        vector< Triplet<mpreal> > tripletListR, tripletListI;
        tripletListR.reserve(matrixR.nonZeros() + b.matrixR.nonZeros());
        if ((isComplex) || (b.isComplex))
            tripletListI.reserve(matrixI.nonZeros() + b.matrixI.nonZeros());

        // Now for each column, we merge the lists of lines with nonzero elements of
        // both matrices. The cases are slightly different depending on the complexity
        // of both matrices
        if (b.isComplex) {
            if (isComplex) {
                for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
                    SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
                    SparseMatrix<mpreal>::InnerIterator itI(matrixI,k);
                    SparseMatrix<mpreal>::InnerIterator itRb(b.matrixR,k);
                    SparseMatrix<mpreal>::InnerIterator itIb(b.matrixI,k);

                    while ((itR) || (itI) || (itRb) || (itIb)) {
                        IndexType row, rowb;
                        if ((itR) || (itI)) {
                            if ((itR) && (itI))
                                row = min(itR.row(), itI.row());
                            else if (itR)
                                row = itR.row();
                            else
                                row = itI.row();
                        }
                        if ((itRb) || (itIb)) {
                            if ((itRb) && (itIb))
                                rowb = min(itRb.row(), itIb.row());
                            else if (itRb)
                                rowb = itRb.row();
                            else
                                rowb = itIb.row();
                        }

                        if (((((itR) || (itI)) && ((itRb) || (itIb))) && (row < rowb)) || ((!itRb) && (!itIb))) {
                            if ((!itI) || ((itR) && (itR.row() < itI.row()))) {
                                if (mpfr::abs(itR.value()) > 0) // Note : this will be wrong if itR.value() is NaN
                                    tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                                ++itR;
                            } else if ((itR) && (itI) && (itR.row() == itI.row())) {
                                if (pow(itR.value(), 2) + pow(itI.value(), 2) > 0) {
                                    tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                                    tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itI.value()));
                                }
                                ++itR;
                                ++itI;
                            } else {
                                if (mpfr::abs(itI.value()) > 0)
                                    tripletListI.push_back(Triplet<mpreal>(itI.row(), k, itI.value()));
                                ++itI;
                            }
                        } else if ((((itR) || (itI)) && ((itRb) || (itIb))) && (row == rowb)) {
                            if ((!itI) || ((itR) && (itR.row() < itI.row()))) {
                                if ((!itIb) || ((itRb) && (itRb.row() < itIb.row()))) {
                                    if ((mpfr::isnan(itRb.value())) || (mpfr::abs(itR.value()) >= mpfr::abs(itRb.value())))
                                        tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                                    else
                                        tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itRb.value()));
                                    ++itRb;
                                } else if ((itRb) && (itIb) && (itRb.row() == itIb.row())) {
                                    if ((mpfr::isnan(pow(itRb.value(), 2) + pow(itIb.value(), 2))) || ((pow(itR.value(), 2) > pow(itRb.value(), 2) + pow(itIb.value(), 2)) || ((pow(itR.value(), 2) == pow(itRb.value(), 2) + pow(itIb.value(), 2)) && (atan2(0, itR.value()) >= atan2(itIb.value(), itRb.value()))))) {
                                        tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                                    } else {
                                        tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itRb.value()));
                                        tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itIb.value()));
                                    }
                                    ++itRb;
                                    ++itIb;
                                } else {
                                    if ((mpfr::isnan(itIb.value())) || ((pow(itR.value(), 2) > pow(itIb.value(), 2)) || ((pow(itR.value(), 2) == pow(itIb.value(), 2)) && (atan2(0, itR.value()) >= atan2(itIb.value(), 0)))))
                                        tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                                    else
                                        tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itIb.value()));
                                    ++itIb;
                                }
                                ++itR;
                            } else if ((itR) && (itI) && (itR.row() == itI.row())) {
                                if ((!itIb) || ((itRb) && (itRb.row() < itIb.row()))) {
                                    if ((mpfr::isnan(itRb.value())) || ((pow(itR.value(), 2) + pow(itI.value(), 2) > pow(itRb.value(), 2)) || ((pow(itR.value(), 2) + pow(itI.value(), 2) == pow(itRb.value(), 2)) && (atan2(itI.value(), itR.value()) >= atan2(0, itRb.value()))))) {
                                        tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                                        tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itI.value()));
                                    } else {
                                        tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itRb.value()));
                                    }
                                    ++itRb;
                                } else if ((itRb) && (itIb) && (itRb.row() == itIb.row())) {
                                    if ((mpfr::isnan(pow(itRb.value(), 2) + pow(itIb.value(), 2))) || ((pow(itR.value(), 2) + pow(itI.value(), 2) > pow(itRb.value(), 2) + pow(itIb.value(), 2)) || ((pow(itR.value(), 2) + pow(itI.value(), 2) == pow(itRb.value(), 2) + pow(itIb.value(), 2)) && (atan2(itI.value(), itR.value()) >= atan2(itIb.value(), itRb.value()))))) {
                                        tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                                        tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itI.value()));
                                    } else {
                                        tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itRb.value()));
                                        tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itIb.value()));
                                    }
                                    ++itRb;
                                    ++itIb;
                                } else {
                                    if ((mpfr::isnan(itIb.value())) || ((pow(itR.value(), 2) + pow(itI.value(), 2) > pow(itIb.value(), 2)) || ((pow(itR.value(), 2) + pow(itI.value(), 2) == pow(itIb.value(), 2)) && (atan2(itI.value(), itR.value()) >= atan2(itIb.value(), 0))))) {
                                        tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                                        tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itI.value()));
                                    } else {
                                        tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itIb.value()));
                                    }
                                    ++itIb;
                                }
                                ++itR;
                                ++itI;
                            } else {
                                if ((!itIb) || ((itRb) && (itRb.row() < itIb.row()))) {
                                    if ((mpfr::isnan(itRb.value())) || ((pow(itI.value(), 2) > pow(itRb.value(), 2)) || ((pow(itI.value(), 2) == pow(itRb.value(), 2)) && (atan2(itI.value(), 0) >= atan2(0, itRb.value()))))) {
                                        tripletListI.push_back(Triplet<mpreal>(itI.row(), k, itI.value()));
                                    } else {
                                        tripletListR.push_back(Triplet<mpreal>(itI.row(), k, itRb.value()));
                                    }
                                    ++itRb;
                                } else if ((itRb) && (itIb) && (itRb.row() == itIb.row())) {
                                    if ((mpfr::isnan(pow(itRb.value(), 2) + pow(itIb.value(), 2))) || ((pow(itI.value(), 2) > pow(itRb.value(), 2) + pow(itIb.value(), 2)) || ((pow(itI.value(), 2) == pow(itRb.value(), 2) + pow(itIb.value(), 2)) && (atan2(itI.value(), 0) >= atan2(itIb.value(), itRb.value()))))) {
                                        tripletListI.push_back(Triplet<mpreal>(itI.row(), k, itI.value()));
                                    } else {
                                        tripletListR.push_back(Triplet<mpreal>(itI.row(), k, itRb.value()));
                                        tripletListI.push_back(Triplet<mpreal>(itI.row(), k, itIb.value()));
                                    }
                                    ++itRb;
                                    ++itIb;
                                } else {
                                    if ((mpfr::isnan(itIb.value())) || ((pow(itI.value(), 2) > pow(itIb.value(), 2)) || ((pow(itI.value(), 2) == pow(itIb.value(), 2)) && (atan2(itI.value(), 0) >= atan2(itIb.value(), 0))))) {
                                        tripletListI.push_back(Triplet<mpreal>(itI.row(), k, itI.value()));
                                    } else {
                                        tripletListI.push_back(Triplet<mpreal>(itI.row(), k, itIb.value()));
                                    }
                                    ++itIb;
                                }
                                ++itI;
                            }
                        } else {
                            if ((!itIb) || ((itRb) && (itRb.row() < itIb.row()))) {
                                if (mpfr::abs(itRb.value()) > 0)
                                    tripletListR.push_back(Triplet<mpreal>(itRb.row(), k, itRb.value()));
                                ++itRb;
                            } else if ((itRb) && (itIb) && (itRb.row() == itIb.row())) {
                                if (pow(itRb.value(), 2) + pow(itIb.value(), 2) > 0) {
                                    tripletListR.push_back(Triplet<mpreal>(itRb.row(), k, itRb.value()));
                                    tripletListI.push_back(Triplet<mpreal>(itRb.row(), k, itIb.value()));
                                }
                                ++itRb;
                                ++itIb;
                            } else {
                                if (mpfr::abs(itIb.value()) > 0)
                                    tripletListI.push_back(Triplet<mpreal>(itIb.row(), k, itIb.value()));
                                ++itIb;
                            }
                        }
                    }
                }
            } else {
                for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
                    SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
                    SparseMatrix<mpreal>::InnerIterator itRb(b.matrixR,k);
                    SparseMatrix<mpreal>::InnerIterator itIb(b.matrixI,k);

                    while ((itR) || (itRb) || (itIb)) {
                        IndexType row, rowb;
                        if (itR)
                            row = itR.row();
                        if ((itRb) || (itIb)) {
                            if ((itRb) && (itIb))
                                rowb = min(itRb.row(), itIb.row());
                            else if (itRb)
                                rowb = itRb.row();
                            else
                                rowb = itIb.row();
                        }

                        if ((((itR) && ((itRb) || (itIb))) && (row < rowb)) || ((!itRb) && (!itIb))) {
                            if (mpfr::abs(itR.value()) > 0)
                                tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                            ++itR;
                        } else if (((itR) && ((itRb) || (itIb))) && (row == rowb)) {
                            if ((!itIb) || ((itRb) && (itRb.row() < itIb.row()))) {
                                if ((mpfr::isnan(itRb.value())) || (mpfr::abs(itR.value()) >= mpfr::abs(itRb.value())))
                                    tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                                else
                                    tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itRb.value()));
                                ++itRb;
                            } else if ((itRb) && (itIb) && (itRb.row() == itIb.row())) {
                                if ((mpfr::isnan(pow(itRb.value(), 2) + pow(itIb.value(), 2))) || ((pow(itR.value(), 2) > pow(itRb.value(), 2) + pow(itIb.value(), 2)) || ((pow(itR.value(), 2) == pow(itRb.value(), 2) + pow(itIb.value(), 2)) && (atan2(0, itR.value()) >= atan2(itIb.value(), itRb.value()))))) {
                                    tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                                } else {
                                    tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itRb.value()));
                                    tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itIb.value()));
                                }
                                ++itRb;
                                ++itIb;
                            } else {
                                if ((mpfr::isnan(itIb.value())) || ((pow(itR.value(), 2) > pow(itIb.value(), 2)) || ((pow(itR.value(), 2) == pow(itIb.value(), 2)) && (atan2(0, itR.value()) >= atan2(itIb.value(), 0)))))
                                    tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                                else
                                    tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itIb.value()));
                                ++itIb;
                            }
                            ++itR;
                        } else {
                            if ((!itIb) || ((itRb) && (itRb.row() < itIb.row()))) {
                                if (mpfr::abs(itRb.value()) > 0)
                                    tripletListR.push_back(Triplet<mpreal>(itRb.row(), k, itRb.value()));
                                ++itRb;
                            } else if ((itRb) && (itIb) && (itRb.row() == itIb.row())) {
                                if (pow(itRb.value(), 2) + pow(itIb.value(), 2) > 0) {
                                    tripletListR.push_back(Triplet<mpreal>(itRb.row(), k, itRb.value()));
                                    tripletListI.push_back(Triplet<mpreal>(itRb.row(), k, itIb.value()));
                                }
                                ++itRb;
                                ++itIb;
                            } else {
                                if (mpfr::abs(itIb.value()) > 0)
                                    tripletListI.push_back(Triplet<mpreal>(itIb.row(), k, itIb.value()));
                                ++itIb;
                            }
                        }
                    }
                }
            }
        } else if (isComplex) {
            for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
                SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
                SparseMatrix<mpreal>::InnerIterator itI(matrixI,k);
                SparseMatrix<mpreal>::InnerIterator itRb(b.matrixR,k);

                while ((itR) || (itI) || (itRb)) {
                    IndexType row, rowb;
                    if ((itR) || (itI)) {
                        if ((itR) && (itI))
                            row = min(itR.row(), itI.row());
                        else if (itR)
                            row = itR.row();
                        else
                            row = itI.row();
                    }
                    if (itRb)
                        rowb = itRb.row();

                    if (((((itR) || (itI)) && (itRb)) && (row < rowb)) || (!itRb)) {
                        if ((!itI) || ((itR) && (itR.row() < itI.row()))) {
                            if (mpfr::abs(itR.value()) > 0)
                                tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                            ++itR;
                        } else if ((itR) && (itI) && (itR.row() == itI.row())) {
                            if (pow(itR.value(), 2) + pow(itI.value(), 2) > 0) {
                                tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                                tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itI.value()));
                            }
                            ++itR;
                            ++itI;
                        } else {
                            if (mpfr::abs(itI.value()) > 0)
                                tripletListI.push_back(Triplet<mpreal>(itI.row(), k, itI.value()));
                            ++itI;
                        }
                    } else if ((((itR) || (itI)) && (itRb)) && (row == rowb)) {
                        if ((!itI) || ((itR) && (itR.row() < itI.row()))) {
                            if ((mpfr::isnan(itRb.value())) || (mpfr::abs(itR.value()) >= mpfr::abs(itRb.value())))
                                tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                            else
                                tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itRb.value()));
                            ++itRb;
                            ++itR;
                        } else if ((itR) && (itI) && (itR.row() == itI.row())) {
                            if ((mpfr::isnan(itRb.value())) || ((pow(itR.value(), 2) + pow(itI.value(), 2) > pow(itRb.value(), 2)) || ((pow(itR.value(), 2) + pow(itI.value(), 2) == pow(itRb.value(), 2)) && (atan2(itI.value(), itR.value()) >= atan2(0, itRb.value()))))) {
                                tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                                tripletListI.push_back(Triplet<mpreal>(itR.row(), k, itI.value()));
                            } else {
                                tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itRb.value()));
                            }
                            ++itRb;
                            ++itR;
                            ++itI;
                        } else {
                            if ((mpfr::isnan(itRb.value())) || ((pow(itI.value(), 2) > pow(itRb.value(), 2)) || ((pow(itI.value(), 2) == pow(itRb.value(), 2)) && (atan2(itI.value(), 0) >= atan2(0, itRb.value()))))) {
                                tripletListI.push_back(Triplet<mpreal>(itI.row(), k, itI.value()));
                            } else {
                                tripletListR.push_back(Triplet<mpreal>(itI.row(), k, itRb.value()));
                            }
                            ++itRb;
                            ++itI;
                        }
                    } else {
                        if (mpfr::abs(itRb.value()) > 0)
                            tripletListR.push_back(Triplet<mpreal>(itRb.row(), k, itRb.value()));
                        ++itRb;
                    }
                }
            }
        } else {
            for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
                SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
                SparseMatrix<mpreal>::InnerIterator itRb(b.matrixR,k);

                while ((itR) || (itRb)) {
                    if ((!itRb) || ((itR) && (itR.row() < itRb.row()))) {
                        if (itR.value() > 0)
                            tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                        ++itR;
                    } else if ((itR) && (itRb) && (itR.row() == itRb.row())) {
                        if ((mpfr::isnan(itRb.value())) || (itR.value() >= itRb.value()))
                            tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itR.value()));
                        else
                            tripletListR.push_back(Triplet<mpreal>(itR.row(), k, itRb.value()));
                        ++itRb;
                        ++itR;
                    } else {
                        if (itRb.value() > 0)
                            tripletListR.push_back(Triplet<mpreal>(itRb.row(), k, itRb.value()));
                        ++itRb;
                    }
                }
            }
        }

        // Now we can assign the data to the sparse matrix
        result.matrixR.setFromTriplets(tripletListR.begin(), tripletListR.end());
        if ((isComplex) || (b.isComplex))
            result.matrixI.setFromTriplets(tripletListI.begin(), tripletListI.end());

        result.checkComplexity();
    }

    return result;
}


// column-wise product b = prod(a)
SparseGmpEigenMatrix SparseGmpEigenMatrix::colProd() const
{
    SparseGmpEigenMatrix result;

    if (isComplex) {
        // We initialize the result matrix
        result.matrixR.resize(1,matrixR.cols());
        result.matrixI.resize(1,matrixI.cols());
        result.matrixR.reserve(min(matrixR.rows(),matrixR.nonZeros()+matrixI.nonZeros()));
        result.matrixI.reserve(min(matrixI.rows(),matrixR.nonZeros()+matrixI.nonZeros()));

        // We compute the products
        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            GmpEigenMatrix product(mpreal(1));
            IndexType previousIndex(0); // This variable is incremented until there is a gap in column enumeration

            SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
            SparseMatrix<mpreal>::InnerIterator itI(matrixI,k);
            while ((itR) || (itI)) {
                if ((itR) && (itI)) {
                    if (itR.row() < itI.row()) {
                        if (previousIndex != itR.row()) {
                            product *= GmpEigenMatrix(mpreal(0));
                            break;
                        }
                        ++previousIndex;
                        product *= GmpEigenMatrix(itR.value());
                        ++itR;
                    } else if (itR.row() == itI.row()) {
                        if (previousIndex != itR.row()) {
                            product *= GmpEigenMatrix(mpreal(0));
                            break;
                        }
                        ++previousIndex;
                        product *= GmpEigenMatrix(itR.value(), itI.value());
                        ++itR;
                        ++itI;
                    } else {
                        if (previousIndex != itI.row()) {
                            product *= GmpEigenMatrix(mpreal(0));
                            break;
                        }
                        ++previousIndex;
                        product *= GmpEigenMatrix(0, itI.value());
                        ++itI;
                    }
                } else if (itR) {
                    if (previousIndex != itR.row()) {
                        product *= GmpEigenMatrix(mpreal(0));
                        break;
                    }
                    ++previousIndex;
                    product *= GmpEigenMatrix(itR.value());
                    ++itR;
                } else {
                    if (previousIndex != itI.row()) {
                        product *= GmpEigenMatrix(mpreal(0));
                        break;
                    }
                    ++previousIndex;
                    product *= GmpEigenMatrix(0, itI.value());
                    ++itI;
                }
            }

            if (previousIndex < matrixR.rows())
                product *= GmpEigenMatrix(mpreal(0));
            if (product.matrixR(0,0) != 0)
                result.matrixR.insert(0,k) = product.matrixR(0,0);
            if (product.matrixI(0,0) != 0)
                result.matrixI.insert(0,k) = product.matrixI(0,0);
        }
        result.matrixR.makeCompressed();
        result.matrixI.makeCompressed();

        result.checkComplexity();
    } else {
        // We initialize the result matrix
        result.isComplex = false;
        result.matrixR.resize(1,matrixR.cols());
        result.matrixR.reserve(min(matrixR.rows(), matrixR.nonZeros()));

        // We compute the products
        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            mpreal product(1);
            IndexType previousIndex(0); // This variable is incremented until there is a gap in column enumeration
            for (SparseMatrix<mpreal>::InnerIterator it(matrixR,k); it; ++it) {
                if (previousIndex != it.row()) {
                    product *= 0;
                    break;
                }
                ++previousIndex;
                product *= it.value();
            }

            if (previousIndex < matrixR.rows())
                product *= 0;
            if (product != 0)
                result.matrixR.insert(0,k) = product;
        }
        result.matrixR.makeCompressed();
    }

    return result;
}

// column-wise product b = prod(a)
SparseGmpEigenMatrix& SparseGmpEigenMatrix::colProd_new() const
{
    SparseGmpEigenMatrix& result(*(new SparseGmpEigenMatrix));

    if (isComplex) {
        // We initialize the result matrix
        result.matrixR.resize(1,matrixR.cols());
        result.matrixI.resize(1,matrixI.cols());
        result.matrixR.reserve(min(matrixR.rows(),matrixR.nonZeros()+matrixI.nonZeros()));
        result.matrixI.reserve(min(matrixI.rows(),matrixR.nonZeros()+matrixI.nonZeros()));

        // We compute the products
        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            GmpEigenMatrix product(mpreal(1));
            IndexType previousIndex(0); // This variable is incremented until there is a gap in column enumeration

            SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
            SparseMatrix<mpreal>::InnerIterator itI(matrixI,k);
            while ((itR) || (itI)) {
                if ((itR) && (itI)) {
                    if (itR.row() < itI.row()) {
                        if (previousIndex != itR.row()) {
                            product *= GmpEigenMatrix(mpreal(0));
                            break;
                        }
                        ++previousIndex;
                        product *= GmpEigenMatrix(itR.value());
                        ++itR;
                    } else if (itR.row() == itI.row()) {
                        if (previousIndex != itR.row()) {
                            product *= GmpEigenMatrix(mpreal(0));
                            break;
                        }
                        ++previousIndex;
                        product *= GmpEigenMatrix(itR.value(), itI.value());
                        ++itR;
                        ++itI;
                    } else {
                        if (previousIndex != itI.row()) {
                            product *= GmpEigenMatrix(mpreal(0));
                            break;
                        }
                        ++previousIndex;
                        product *= GmpEigenMatrix(0, itI.value());
                        ++itI;
                    }
                } else if (itR) {
                    if (previousIndex != itR.row()) {
                        product *= GmpEigenMatrix(mpreal(0));
                        break;
                    }
                    ++previousIndex;
                    product *= GmpEigenMatrix(itR.value());
                    ++itR;
                } else {
                    if (previousIndex != itI.row()) {
                        product *= GmpEigenMatrix(mpreal(0));
                        break;
                    }
                    ++previousIndex;
                    product *= GmpEigenMatrix(0, itI.value());
                    ++itI;
                }
            }

            if (previousIndex < matrixR.rows())
                product *= GmpEigenMatrix(mpreal(0));
            if (product.matrixR(0,0) != 0)
                result.matrixR.insert(0,k) = product.matrixR(0,0);
            if (product.matrixI(0,0) != 0)
                result.matrixI.insert(0,k) = product.matrixI(0,0);
        }
        result.matrixR.makeCompressed();
        result.matrixI.makeCompressed();

        result.checkComplexity();
    } else {
        // We initialize the result matrix
        result.isComplex = false;
        result.matrixR.resize(1,matrixR.cols());
        result.matrixR.reserve(min(matrixR.rows(), matrixR.nonZeros()));

        // We compute the products
        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            mpreal product(1);
            IndexType previousIndex(0); // This variable is incremented until there is a gap in column enumeration
            for (SparseMatrix<mpreal>::InnerIterator it(matrixR,k); it; ++it) {
                if (previousIndex != it.row()) {
                    product *= 0;
                    break;
                }
                ++previousIndex;
                product *= it.value();
            }

            if (previousIndex < matrixR.rows())
                product *= 0;
            if (product != 0)
                result.matrixR.insert(0,k) = product;
        }
        result.matrixR.makeCompressed();
    }

    return result;
}

// column-wise product b = prod(a,2)
SparseGmpEigenMatrix SparseGmpEigenMatrix::rowProd() const
{
    SparseGmpEigenMatrix result;

    if (isComplex) {
        // We initialize the result matrix
        GmpEigenMatrix products(constantMatrix(matrixR.rows(),1,1));
        // We add an explicit empty complex part
        products.isComplex = true;
        products.matrixI = constantMatrix(matrixR.rows(),1,0);

        // We compute the products
        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            IndexType previousNonZero(-1);

            SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
            SparseMatrix<mpreal>::InnerIterator itI(matrixI,k);

            while ((itR) || (itI)) {
                if ((itR) && (itI)) {
                    if (itR.row() < itI.row()) {
                        for (IndexType i(previousNonZero+1); i < itR.row(); ++i)
                            products.subsasgn(vector<IndexType>(1,i), products.block(i,0,1,1)*0);
                        previousNonZero = itR.row();
                        GmpEigenMatrix tmp(products.block(itR.row(),0,1,1));
                        tmp *= GmpEigenMatrix(itR.value());
                        products.subsasgn(vector<IndexType>(1,itR.row()), tmp);
                        ++itR;
                    } else if (itR.row() == itI.row()) {
                        for (IndexType i(previousNonZero+1); i < itR.row(); ++i)
                            products.subsasgn(vector<IndexType>(1,i), products.block(i,0,1,1)*0);
                        previousNonZero = itR.row();
                        GmpEigenMatrix tmp(products.block(itR.row(),0,1,1));
                        tmp *= GmpEigenMatrix(itR.value(), itI.value());
                        products.subsasgn(vector<IndexType>(1,itR.row()), tmp);
                        ++itR;
                        ++itI;
                    } else {
                        for (IndexType i(previousNonZero+1); i < itI.row(); ++i)
                            products.subsasgn(vector<IndexType>(1,i), products.block(i,0,1,1)*0);
                        previousNonZero = itI.row();
                        GmpEigenMatrix tmp(products.block(itI.row(),0,1,1));
                        tmp *= GmpEigenMatrix(0, itI.value());
                        products.subsasgn(vector<IndexType>(1,itI.row()), tmp);
                        ++itI;
                    }
                } else if (itR) {
                    for (IndexType i(previousNonZero+1); i < itR.row(); ++i)
                        products.subsasgn(vector<IndexType>(1,i), products.block(i,0,1,1)*0);
                    previousNonZero = itR.row();
                    GmpEigenMatrix tmp(products.block(itR.row(),0,1,1));
                    tmp *= GmpEigenMatrix(itR.value());
                    products.subsasgn(vector<IndexType>(1,itR.row()), tmp);
                    ++itR;
                } else {
                    for (IndexType i(previousNonZero+1); i < itI.row(); ++i)
                        products.subsasgn(vector<IndexType>(1,i), products.block(i,0,1,1)*0);
                    previousNonZero = itI.row();
                    GmpEigenMatrix tmp(products.block(itI.row(),0,1,1));
                    tmp *= GmpEigenMatrix(0, itI.value());
                    products.subsasgn(vector<IndexType>(1,itI.row()), tmp);
                    ++itI;
                }
            }

            // If there are still remaining zeros
            for (IndexType i(previousNonZero+1); i < matrixR.rows(); ++i) {
                products.subsasgn(vector<IndexType>(1,i), products.block(i,0,1,1)*0);
            }
        }

        // We copy the products into the output matrix
        result.matrixR = products.matrixR.sparseView(mpreal(0),mpreal(1));
        result.matrixI = products.matrixI.sparseView(mpreal(0),mpreal(1));

        result.checkComplexity();
    } else {
        // We initialize the result matrix
        GmpEigenMatrix products(constantMatrix(matrixR.rows(),1,1));

        // We compute the products
        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            IndexType previousNonZero(-1);

            SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);

            while (itR) {
                for (IndexType i(previousNonZero+1); i < itR.row(); ++i)
                    products.matrixR(i,0) *= 0;
                previousNonZero = itR.row();
                products.matrixR(itR.row(),0) *= itR.value();
                ++itR;
            }

            // If there are still remaining zeros
            for (IndexType i(previousNonZero+1); i < matrixR.rows(); ++i)
                products.matrixR(i,0) *= 0;
        }

        // We copy the products into the output matrix
        result.matrixR = products.matrixR.sparseView(mpreal(0),mpreal(1));
        result.isComplex = false;
    }

    return result;
}

// column-wise product b = prod(a,2)
SparseGmpEigenMatrix& SparseGmpEigenMatrix::rowProd_new() const
{
    SparseGmpEigenMatrix& result(*(new SparseGmpEigenMatrix));

    if (isComplex) {
        // We initialize the result matrix
        GmpEigenMatrix products(constantMatrix(matrixR.rows(),1,1));
        // We add an explicit empty complex part
        products.isComplex = true;
        products.matrixI = constantMatrix(matrixR.rows(),1,0);

        // We compute the products
        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            IndexType previousNonZero(-1);

            SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
            SparseMatrix<mpreal>::InnerIterator itI(matrixI,k);

            while ((itR) || (itI)) {
                if ((itR) && (itI)) {
                    if (itR.row() < itI.row()) {
                        for (IndexType i(previousNonZero+1); i < itR.row(); ++i)
                            products.subsasgn(vector<IndexType>(1,i), products.block(i,0,1,1)*0);
                        previousNonZero = itR.row();
                        GmpEigenMatrix tmp(products.block(itR.row(),0,1,1));
                        tmp *= GmpEigenMatrix(itR.value());
                        products.subsasgn(vector<IndexType>(1,itR.row()), tmp);
                        ++itR;
                    } else if (itR.row() == itI.row()) {
                        for (IndexType i(previousNonZero+1); i < itR.row(); ++i)
                            products.subsasgn(vector<IndexType>(1,i), products.block(i,0,1,1)*0);
                        previousNonZero = itR.row();
                        GmpEigenMatrix tmp(products.block(itR.row(),0,1,1));
                        tmp *= GmpEigenMatrix(itR.value(), itI.value());
                        products.subsasgn(vector<IndexType>(1,itR.row()), tmp);
                        ++itR;
                        ++itI;
                    } else {
                        for (IndexType i(previousNonZero+1); i < itI.row(); ++i)
                            products.subsasgn(vector<IndexType>(1,i), products.block(i,0,1,1)*0);
                        previousNonZero = itI.row();
                        GmpEigenMatrix tmp(products.block(itI.row(),0,1,1));
                        tmp *= GmpEigenMatrix(0, itI.value());
                        products.subsasgn(vector<IndexType>(1,itI.row()), tmp);
                        ++itI;
                    }
                } else if (itR) {
                    for (IndexType i(previousNonZero+1); i < itR.row(); ++i)
                        products.subsasgn(vector<IndexType>(1,i), products.block(i,0,1,1)*0);
                    previousNonZero = itR.row();
                    GmpEigenMatrix tmp(products.block(itR.row(),0,1,1));
                    tmp *= GmpEigenMatrix(itR.value());
                    products.subsasgn(vector<IndexType>(1,itR.row()), tmp);
                    ++itR;
                } else {
                    for (IndexType i(previousNonZero+1); i < itI.row(); ++i)
                        products.subsasgn(vector<IndexType>(1,i), products.block(i,0,1,1)*0);
                    previousNonZero = itI.row();
                    GmpEigenMatrix tmp(products.block(itI.row(),0,1,1));
                    tmp *= GmpEigenMatrix(0, itI.value());
                    products.subsasgn(vector<IndexType>(1,itI.row()), tmp);
                    ++itI;
                }
            }

            // If there are still remaining zeros
            for (IndexType i(previousNonZero+1); i < matrixR.rows(); ++i) {
                products.subsasgn(vector<IndexType>(1,i), products.block(i,0,1,1)*0);
            }
        }

        // We copy the products into the output matrix
        result.matrixR = products.matrixR.sparseView(mpreal(0),mpreal(1));
        result.matrixI = products.matrixI.sparseView(mpreal(0),mpreal(1));

        result.checkComplexity();
    } else {
        // We initialize the result matrix
        GmpEigenMatrix products(constantMatrix(matrixR.rows(),1,1));

        // We compute the products
        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            IndexType previousNonZero(-1);

            SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);

            while (itR) {
                for (IndexType i(previousNonZero+1); i < itR.row(); ++i)
                    products.matrixR(i,0) *= 0;
                previousNonZero = itR.row();
                products.matrixR(itR.row(),0) *= itR.value();
                ++itR;
            }

            // If there are still remaining zeros
            for (IndexType i(previousNonZero+1); i < matrixR.rows(); ++i)
                products.matrixR(i,0) *= 0;
        }

        // We copy the products into the output matrix
        result.matrixR = products.matrixR.sparseView(mpreal(0),mpreal(1));
        result.isComplex = false;
    }

    return result;
}


/* This function sorts out elements in a matrix */
SparseGmpEigenMatrix SparseGmpEigenMatrix::sort(const int& dim, const int& type, vector < vector < IndexType > >& index, vector < IndexType >& nbNegatives) const
{
    SparseGmpEigenMatrix result(*this);

    // We prepare the matrices
    result.matrixR.resize(matrixR.rows(),matrixR.cols());
    result.matrixR.reserve(matrixR.nonZeros());
    if (isComplex) {
        result.isComplex = true;
        result.matrixI.resize(matrixR.rows(), matrixR.cols());
        result.matrixI.reserve(matrixI.nonZeros());
    }

    // We assume that index and nbNegatives are already empty
    //index.clear();
    //nbNegatives.clear();

    if (dim == 0) {
        for (IndexType k(0); k < matrixR.outerSize(); ++k) {
            // First, we find the non-zero elements which are in this column
            // as well as how many elements are smaller than zero
            vector < IndexType > indexMap(0, 0);
            nbNegatives.push_back(0);
            if (isComplex) {
                SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
                SparseMatrix<mpreal>::InnerIterator itI(matrixI,k);
                while ((itR) || (itI)) {
                    if ((itR) && (itI)) {
                        if (itR.row() < itI.row()) {
                            indexMap.push_back(itR.row());
                            ++itR;
                        } else if (itR.row() == itI.row()) {
                            indexMap.push_back(itR.row());
                            ++itR;
                            ++itI;
                        } else {
                            indexMap.push_back(itI.row());
                            ++itI;
                        }
                    } else if (itR) {
                        indexMap.push_back(itR.row());
                        ++itR;
                    } else {
                        indexMap.push_back(itI.row());
                        ++itI;
                    }
                }
            } else {
                SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
                while (itR) {
                    indexMap.push_back(itR.row());
                    if (itR.value() < 0)
                        ++nbNegatives[k];
                    ++itR;
                }
            }

            // Now we sort them
            if (isComplex)
                std::sort(indexMap.begin(), indexMap.end(), compareThroughIndicesComplex< Matrix< mpreal, Dynamic, 1 >, Matrix< mpreal, Dynamic, 1 > >(matrixR.col(k), matrixI.col(k)));
            else
                std::sort(indexMap.begin(), indexMap.end(), compareThroughIndices< Matrix< mpreal, Dynamic, 1 > >(matrixR.col(k)));

            if (type == 1)
                std::reverse(indexMap.begin(), indexMap.end());

            index.push_back(indexMap);

            IndexType startBlock, endBlock;
            if (type == 0) {
                startBlock = nbNegatives[k];
                endBlock = indexMap.size() - nbNegatives[k];
            } else {
                startBlock = indexMap.size() - nbNegatives[k];
                endBlock = nbNegatives[k];
            }

            // We construct the sorted column
            for (IndexType i(0); i < indexMap.size(); ++i) {
                IndexType x(i);
                if (i >= startBlock)
                    x = matrixR.innerSize()-indexMap.size()+i;

                if (isComplex) {
                    if (matrixR.coeff(indexMap[i],k) != 0)
                        result.matrixR.insert(x,k) = matrixR.coeff(indexMap[i],k);
                    if (matrixI.coeff(indexMap[i],k) != 0)
                        result.matrixI.insert(x,k) = matrixI.coeff(indexMap[i],k);
                } else {
                    result.matrixR.insert(x,k) = matrixR.coeff(indexMap[i],k);
                }
            }
        }
    } else {
        // We start by transposing the matrix
        SparseGmpEigenMatrix transpo(transpose());

        for (IndexType k(0); k < transpo.matrixR.outerSize(); ++k) {
            // First, we find the non-zero elements which are in this column
            // as well as how many elements are smaller than zero
            vector < IndexType > indexMap(0, 0);
            nbNegatives.push_back(0);
            if (isComplex) {
                SparseMatrix<mpreal>::InnerIterator itR(transpo.matrixR,k);
                SparseMatrix<mpreal>::InnerIterator itI(transpo.matrixI,k);
                while ((itR) || (itI)) {
                    if ((itR) && (itI)) {
                        if (itR.row() < itI.row()) {
                            indexMap.push_back(itR.row());
                            ++itR;
                        } else if (itR.row() == itI.row()) {
                            indexMap.push_back(itR.row());
                            ++itR;
                            ++itI;
                        } else {
                            indexMap.push_back(itI.row());
                            ++itI;
                        }
                    } else if (itR) {
                        indexMap.push_back(itR.row());
                        ++itR;
                    } else {
                        indexMap.push_back(itI.row());
                        ++itI;
                    }
                }
            } else {
                SparseMatrix<mpreal>::InnerIterator itR(transpo.matrixR,k);
                while (itR) {
                    indexMap.push_back(itR.row());
                    if (itR.value() < 0)
                        ++nbNegatives[k];
                    ++itR;
                }
            }

            // Now we sort them
            if (isComplex)
                std::sort(indexMap.begin(), indexMap.end(), compareThroughIndicesComplex< Matrix< mpreal, Dynamic, 1 >, Matrix< mpreal, Dynamic, 1 > >(transpo.matrixR.col(k), transpo.matrixI.col(k)));
            else
                std::sort(indexMap.begin(), indexMap.end(), compareThroughIndices< Matrix< mpreal, Dynamic, 1 > >(transpo.matrixR.col(k)));

            if (type == 1)
                std::reverse(indexMap.begin(), indexMap.end());

            index.push_back(indexMap);

            IndexType startBlock, endBlock;
            if (type == 0) {
                startBlock = nbNegatives[k];
                endBlock = indexMap.size() - nbNegatives[k];
            } else {
                startBlock = indexMap.size() - nbNegatives[k];
                endBlock = nbNegatives[k];
            }

            // We construct the sorted column
            for (IndexType i(0); i < indexMap.size(); ++i) {
                IndexType x(i);
                if (i >= startBlock)
                    x = matrixR.innerSize()-indexMap.size()+i;

                if (isComplex) {
                    if (transpo.matrixR.coeff(indexMap[i],k) != 0)
                        result.matrixR.insert(k,x) = transpo.matrixR.coeff(indexMap[i],k);
                    if (transpo.matrixI.coeff(indexMap[i],k) != 0)
                        result.matrixI.insert(k,x) = transpo.matrixI.coeff(indexMap[i],k);
                } else {
                    result.matrixR.insert(k,x) = transpo.matrixR.coeff(indexMap[i],k);
                }
            }
        }
    }

    result.matrixR.makeCompressed();
    if (isComplex)
        result.matrixI.makeCompressed();

    return result;
}

/* This function sorts out elements in the columns of a matrix */
SparseGmpEigenMatrix& SparseGmpEigenMatrix::sort_new(const int& dim, const int& type, vector < vector < IndexType > >& index, vector < IndexType >& nbNegatives) const
{
    SparseGmpEigenMatrix& result(*(new SparseGmpEigenMatrix));

    // We prepare the matrices
    result.matrixR.resize(matrixR.rows(),matrixR.cols());
    result.matrixR.reserve(matrixR.nonZeros());
    if (isComplex) {
        result.isComplex = true;
        result.matrixI.resize(matrixR.rows(), matrixR.cols());
        result.matrixI.reserve(matrixI.nonZeros());
    }

    // We assume that index and nbNegatives are already empty
    //index.clear();
    //nbNegatives.clear();

    if (dim == 0) {
        for (IndexType k(0); k < matrixR.outerSize(); ++k) {
            // First, we find the non-zero elements which are in this column
            // as well as how many elements are smaller than zero
            vector < IndexType > indexMap(0, 0);
            nbNegatives.push_back(0);
            if (isComplex) {
                SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
                SparseMatrix<mpreal>::InnerIterator itI(matrixI,k);
                while ((itR) || (itI)) {
                    if ((itR) && (itI)) {
                        if (itR.row() < itI.row()) {
                            indexMap.push_back(itR.row());
                            ++itR;
                        } else if (itR.row() == itI.row()) {
                            indexMap.push_back(itR.row());
                            ++itR;
                            ++itI;
                        } else {
                            indexMap.push_back(itI.row());
                            ++itI;
                        }
                    } else if (itR) {
                        indexMap.push_back(itR.row());
                        ++itR;
                    } else {
                        indexMap.push_back(itI.row());
                        ++itI;
                    }
                }
            } else {
                SparseMatrix<mpreal>::InnerIterator itR(matrixR,k);
                while (itR) {
                    indexMap.push_back(itR.row());
                    if (itR.value() < 0)
                        ++nbNegatives[k];
                    ++itR;
                }
            }

            // Now we sort them
            if (isComplex)
                std::sort(indexMap.begin(), indexMap.end(), compareThroughIndicesComplex< Matrix< mpreal, Dynamic, 1 >, Matrix< mpreal, Dynamic, 1 > >(matrixR.col(k), matrixI.col(k)));
            else
                std::sort(indexMap.begin(), indexMap.end(), compareThroughIndices< Matrix< mpreal, Dynamic, 1 > >(matrixR.col(k)));

            if (type == 1)
                std::reverse(indexMap.begin(), indexMap.end());

            index.push_back(indexMap);

            IndexType startBlock, endBlock;
            if (type == 0) {
                startBlock = nbNegatives[k];
                endBlock = indexMap.size() - nbNegatives[k];
            } else {
                startBlock = indexMap.size() - nbNegatives[k];
                endBlock = nbNegatives[k];
            }

            // We construct the sorted column
            for (IndexType i(0); i < indexMap.size(); ++i) {
                IndexType x(i);
                if (i >= startBlock)
                    x = matrixR.innerSize()-indexMap.size()+i;

                if (isComplex) {
                    if (matrixR.coeff(indexMap[i],k) != 0)
                        result.matrixR.insert(x,k) = matrixR.coeff(indexMap[i],k);
                    if (matrixI.coeff(indexMap[i],k) != 0)
                        result.matrixI.insert(x,k) = matrixI.coeff(indexMap[i],k);
                } else {
                    result.matrixR.insert(x,k) = matrixR.coeff(indexMap[i],k);
                }
            }
        }
    } else {
        // We start by transposing the matrix
        SparseGmpEigenMatrix transpo(transpose());

        for (IndexType k(0); k < transpo.matrixR.outerSize(); ++k) {
            // First, we find the non-zero elements which are in this column
            // as well as how many elements are smaller than zero
            vector < IndexType > indexMap(0, 0);
            nbNegatives.push_back(0);
            if (isComplex) {
                SparseMatrix<mpreal>::InnerIterator itR(transpo.matrixR,k);
                SparseMatrix<mpreal>::InnerIterator itI(transpo.matrixI,k);
                while ((itR) || (itI)) {
                    if ((itR) && (itI)) {
                        if (itR.row() < itI.row()) {
                            indexMap.push_back(itR.row());
                            ++itR;
                        } else if (itR.row() == itI.row()) {
                            indexMap.push_back(itR.row());
                            ++itR;
                            ++itI;
                        } else {
                            indexMap.push_back(itI.row());
                            ++itI;
                        }
                    } else if (itR) {
                        indexMap.push_back(itR.row());
                        ++itR;
                    } else {
                        indexMap.push_back(itI.row());
                        ++itI;
                    }
                }
            } else {
                SparseMatrix<mpreal>::InnerIterator itR(transpo.matrixR,k);
                while (itR) {
                    indexMap.push_back(itR.row());
                    if (itR.value() < 0)
                        ++nbNegatives[k];
                    ++itR;
                }
            }

            // Now we sort them
            if (isComplex)
                std::sort(indexMap.begin(), indexMap.end(), compareThroughIndicesComplex< Matrix< mpreal, Dynamic, 1 >, Matrix< mpreal, Dynamic, 1 > >(transpo.matrixR.col(k), transpo.matrixI.col(k)));
            else
                std::sort(indexMap.begin(), indexMap.end(), compareThroughIndices< Matrix< mpreal, Dynamic, 1 > >(transpo.matrixR.col(k)));

            if (type == 1)
                std::reverse(indexMap.begin(), indexMap.end());

            index.push_back(indexMap);

            IndexType startBlock, endBlock;
            if (type == 0) {
                startBlock = nbNegatives[k];
                endBlock = indexMap.size() - nbNegatives[k];
            } else {
                startBlock = indexMap.size() - nbNegatives[k];
                endBlock = nbNegatives[k];
            }

            // We construct the sorted column
            for (IndexType i(0); i < indexMap.size(); ++i) {
                IndexType x(i);
                if (i >= startBlock)
                    x = matrixR.innerSize()-indexMap.size()+i;

                if (isComplex) {
                    if (transpo.matrixR.coeff(indexMap[i],k) != 0)
                        result.matrixR.insert(k,x) = transpo.matrixR.coeff(indexMap[i],k);
                    if (transpo.matrixI.coeff(indexMap[i],k) != 0)
                        result.matrixI.insert(k,x) = transpo.matrixI.coeff(indexMap[i],k);
                } else {
                    result.matrixR.insert(k,x) = transpo.matrixR.coeff(indexMap[i],k);
                }
            }
        }
    }

    result.matrixR.makeCompressed();
    if (isComplex)
        result.matrixI.makeCompressed();

    return result;
}



/* This function transforms complex matrices into twice as big real matrices */
SparseGmpEigenMatrix SparseGmpEigenMatrix::complexIsometry() const
{
    SparseGmpEigenMatrix result;
    SparseGmpEigenMatrix tmpR;
    SparseGmpEigenMatrix tmpI;

    tmpR.matrixR = matrixR;
    if (isComplex)
        tmpI.matrixR = matrixI;
    else
        tmpI.matrixR.resize(matrixR.rows(),matrixR.cols());

    result = (tmpR.horzcat(tmpI)).vertcat((-tmpI).horzcat(tmpR));

    return result;
}

/* This function restores the complex matrix corresponding to a big real
   matrix created with the complexIsometry function. */
SparseGmpEigenMatrix SparseGmpEigenMatrix::complexIsometryInverse() const
{
    SparseGmpEigenMatrix result;

    if (isComplex) {
        mexErrMsgTxt("Error in complexIsometryInverse : the provided matrix is not real.");
    }

    result.matrixR = matrixR.block(0, 0, matrixR.rows()/2, matrixR.cols()/2);
    result.matrixI = matrixR.block(0, matrixR.cols()/2, matrixR.rows()/2, matrixR.cols()/2);
    result.checkComplexity();

    return result;
}
