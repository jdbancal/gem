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



/* Construction from a sparse GmpEigenMatrix */
GmpEigenMatrix::GmpEigenMatrix(const SparseGmpEigenMatrix& a) : matrixR(a.matrixR), matrixI(a.matrixI), isComplex(a.isComplex) {}


/* Construcion from a matlab struct containing all the class information. */
GmpEigenMatrix::GmpEigenMatrix(const mxArray* prhs) {
    if ((!mxIsStruct(prhs)) || ((mxGetNumberOfElements(prhs) != 1)))
        mexErrMsgTxt("The argument is not a struct array with one element.");

    // The structure should at least contain the following fields :
    // precisionR and matrixR
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

    mwSize m(mxGetM(field));
    mwSize n(mxGetN(field));
    if ((mxGetM(precisionField) != m) || (mxGetN(precisionField) != n))
        mexErrMsgTxt("Different fields have incompatible sizes.");

    // We set the size
    matrixR.resize(m,n);

    // Now we copy the real data
    double* precisionPtr(mxGetPr(precisionField));
    mwIndex index(0);
    for (mwIndex j = 0; j < n; ++j) {
        // We first iterate on the rows
        for (mwIndex i = 0; i < m; ++i) {
            const mxArray* cellElementPtr(mxGetCell(field, index));
            if (cellElementPtr == NULL)
                mexErrMsgTxt("Cell element is empty.");
            else {
                if (!mxIsChar(cellElementPtr))
                    mexErrMsgTxt("Cell element is not a string.");
                else {
                    mwSize stringLength = mxGetNumberOfElements(cellElementPtr);
                    if (stringLength == 0)
                        mexErrMsgTxt("String is of length 0.");
                    else {
                        char c_string[stringLength+1];
                        if (mxGetString(cellElementPtr, c_string, stringLength+1) != 0)
                            mexErrMsgTxt("Error in extracting a string from a cell array");
                        matrixR(i,j) = mpreal(string(c_string), precisionPtr[i]);
                    }
                }
            }
            index += 1;
        }
        precisionPtr += m;
    }

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

        if ((mxGetM(precisionField) != m) || (mxGetN(precisionField) != n) || (mxGetM(field) != m) || (mxGetN(field) != n))
            mexErrMsgTxt("Different fields have incompatible sizes.");

        // We set the size
        matrixI.resize(m,n);

        // Now we copy the imaginary data
        double* precisionPtr(mxGetPr(precisionField));
        mwIndex index(0);
        for (mwIndex j = 0; j < n; ++j) {
            // We first iterate on the rows
            for (mwIndex i = 0; i < m; ++i) {
                const mxArray* cellElementPtr(mxGetCell(field, index));
                if (cellElementPtr == NULL)
                    mexErrMsgTxt("Cell element is empty.");
                else {
                    if (!mxIsChar(cellElementPtr))
                        mexErrMsgTxt("Cell element is not a string.");
                    else {
                        mwSize stringLength = mxGetNumberOfElements(cellElementPtr);
                        if (stringLength == 0)
                            mexErrMsgTxt("String is of length 0.");
                        else {
                            char c_string[stringLength+1];
                            if (mxGetString(cellElementPtr, c_string, stringLength+1) != 0)
                                mexErrMsgTxt("Error in extracting a string from a cell array");
                            matrixI(i,j) = mpreal(string(c_string), precisionPtr[i]);
                        }
                    }
                }
                index += 1;
            }
            precisionPtr += m;
        }
    } else
        isComplex = false;
}

/* Construction from a matlab table or cell array. For double element, this
   constructor takes exactly 15 digits from them and sets all other digits
   to 0. For strings, it takes all the provided digits (up to precision). */
GmpEigenMatrix::GmpEigenMatrix(const mxArray* prhs, const int& precision) {

    // We set the required precision
    mp_prec_t precisionInBits(mpfr::digits2bits(precision));

    // Now we check whether the matlab object is an array or a cell array
    if (mxIsDouble(prhs))
    { // The matlab object is an array
        /* Get the size and pointers to input data */
        mwSize m(mxGetM(prhs));
        mwSize n(mxGetN(prhs));
        double* pr(mxGetPr(prhs));
        double* pi(mxGetPi(prhs));
        isComplex = (pi==NULL ? 0 : 1);

        // If the input is real, we set the imaginary data to nothing
        if (isComplex == 0)
            matrixI = Matrix<mpreal,Dynamic,Dynamic>();

        // We set the size
        matrixR.resize(m,n);
        if (isComplex)
            matrixI.resize(m,n);

        // Now we copy the data over
        for (mwIndex j = 0; j < n; ++j) {
            // We first iterate on the rows
            for (mwIndex i = 0; i < m; ++i) {
                matrixR(i,j) = mpreal(toString15(pr[i]), precisionInBits);
                if (isComplex) {
                    matrixI(i,j) = mpreal(toString15(pi[i]), precisionInBits);
                }
            }
            pr += m;
            pi += m;
        }
    }
    else if (mxIsCell(prhs))
    { // The matlab object is a cell array
        mwSize m(mxGetM(prhs));
        mwSize n(mxGetN(prhs));

        // We set the size
        matrixR.resize(m,n);
        if (isComplex)
            matrixI.resize(m,n);

        // We don't support construction from strings with an imaginary part
        isComplex = false;

        // Now we iterate on all the elements of the cell array
        mwIndex index(0);
        for (mwIndex j = 0; j < n; ++j) {
            // We first iterate on the rows
            for (mwIndex i = 0; i < m; ++i) {
                const mxArray* cellElementPtr(mxGetCell(prhs, index));
                if (cellElementPtr == NULL)
                {
                    // The cell element is empty
                    matrixR(i,j) = mpreal("0", precisionInBits);
                } else {
                    // We check whether the element is a number or a string
                    if (mxIsDouble(cellElementPtr)) {
                        double* pr(mxGetPr(cellElementPtr));
                        matrixR(i,j) = mpreal(toString15(pr[0]), precisionInBits);
                    } else if(mxIsChar(cellElementPtr)) {
                        mwSize stringLength = mxGetNumberOfElements(cellElementPtr);
                        if (stringLength == 0)
                            matrixR(i,j) = mpreal("0", precisionInBits);
                        else {
                            char c_string[stringLength+1];
                            if (mxGetString(cellElementPtr, c_string, stringLength+1) != 0)
                                mexErrMsgTxt("Error in extracting a string from a cell array");
                            matrixR(i,j) = mpreal(string(c_string), precisionInBits);
                        }
                    } else
                        mexErrMsgTxt("Type of cell array element not supported.");
                }
                index += 1;
            }
        }
    }
}


/* Saving procedure
   This function returns a structure to matlab which contains all the data
   that define the GemEigenMatrix instance */
mxArray* GmpEigenMatrix::saveobj() const
{
    // All possible fields of the structure we will return
    const char *field_names[] = {"precisionR", "matrixR", "precisionI", "matrixI"};

    // We create the matlab structure to be populated
    mxArray* ptr;
    if (isComplex)
        ptr = mxCreateStructMatrix(1, 1, 4, field_names);
    else
        ptr = mxCreateStructMatrix(1, 1, 2, field_names);

    // Now, we add each field
    // We start with the real part
    mwSize m(matrixR.rows());
    mwSize n(matrixR.cols());
    mxArray* precisionField(mxCreateNumericMatrix(m, n, mxDOUBLE_CLASS, mxREAL));
    mxArray* field(mxCreateCellMatrix(m, n));

    // Now we iterate on all the elements of the cell array
    mwIndex index(0);
    double* pointerR(mxGetPr(precisionField));
    for (mwIndex j = 0; j < n; ++j) {
        // We first iterate on the rows
        for (mwIndex i = 0; i < m; ++i) {
            pointerR[i] = double(matrixR(i,j).get_prec());
            mxSetCell(field, index, mxCreateString(matrixR(i,j).toString().c_str()));
            ++index;
        }
        pointerR += m;
    }

    // We save these objects into the matlab structure
    mxSetField(ptr, 0, "precisionR", precisionField);
    mxSetField(ptr, 0, "matrixR", field);

    // Now we also add the imaginary part is there is one
    if (isComplex)
    {
        precisionField = mxCreateNumericMatrix(m, n, mxDOUBLE_CLASS, mxREAL);
        field = mxCreateCellMatrix(m, n);

        // Now we iterate on all the elements of the cell array
        mwIndex index(0);
        double* pointerI(mxGetPr(precisionField));
        for (mwIndex j = 0; j < n; ++j) {
            // We first iterate on the rows
            for (mwIndex i = 0; i < m; ++i) {
                pointerI[i] = double(matrixI(i,j).get_prec());
                mxSetCell(field, index, mxCreateString(matrixI(i,j).toString().c_str()));
                ++index;
            }
            pointerI += m;
        }

        // We save these objects into the matlab structure
        mxSetField(ptr, 0, "precisionI", precisionField);
        mxSetField(ptr, 0, "matrixI", field);
    }

    return ptr;
}


/* Printing of the object in matlab.
   Precision is the maximum number of significative digits that we wish to print.
   Due to alignment, all numbers will not be printed with the same number of
   significative digits. And if some numbers have fewer digits, the missing digits
   are left blank.

   Note : This function implements rounding as provided by the mpfr. Hence, 0.995
          is rounded up to 0.99 rather than 1.00, when cut after two digits (!)
          This is what they call 'correct rounding' or 'round half to even'.
   */
void GmpEigenMatrix::display(const int& precision) const
{
    // If the matrix is empty, we print nothing
    if (matrixR.rows()*matrixR.cols() == 0)
    {
        mexPrintf("  []\n");
        return;
    }

    // We start by checking whether the contains only integers. Then we compute the
    // power of 10 corresponding to the largest number in the matrix, after being
    // rounded to just the required precision (e.g. 9.95 is rounded to 10.0 if we
    // can only use two digits of precision: the first digit needed is given by
    //   floor(log10(x/(1-5*10^(-1-precision)))))
    // where precision is the number of significant digits at which we aim to cut
    // the description of x).
    // Special care needs to be taken when computing the number of digits required
    // to represent zero, nan or inf.
    //
    // Note : here we assume rounding is as mathematical rounding would be. This
    // always leaves enough room for mpfr's default rounding...
    //
    // Note : we construct all mpfr numbers from strings, otherwise they
    //        are not well converted (they may not finish with zeros)
    bool integerMatrix(isInt());
//        cout << "integerMatrix = " << integerMatrix << endl;
    int largestExponent((int)((iszero(matrixR(0,0)) || mpfr::isnan(matrixR(0,0)) || mpfr::isinf(matrixR(0,0))) ? INT_MIN : mpfr::floor(mpfr::log10(mpfr::abs((matrixR(0,0)/(mpreal("1.0") - mpreal("5.0")*mpfr::pow(10,-1-precision)))))).toLong()));
    for (mwIndex i = 0; i < matrixR.rows(); ++i) {
        for (mwIndex j = 0; j < matrixR.cols(); ++j) {
            largestExponent = max(largestExponent, (int)((iszero(matrixR(i,j)) || mpfr::isnan(matrixR(i,j)) || mpfr::isinf(matrixR(i,j))) ? INT_MIN : mpfr::floor(mpfr::log10(mpfr::abs((matrixR(i,j)/(mpreal("1.0") - mpreal("5.0")*mpfr::pow(10,-1-precision)))))).toLong()));
            if (isComplex)
                largestExponent = max(largestExponent, (int)((iszero(matrixI(i,j)) || mpfr::isnan(matrixI(i,j)) || mpfr::isinf(matrixI(i,j))) ? INT_MIN : mpfr::floor(mpfr::log10(mpfr::abs((matrixI(i,j)/(mpreal("1.0") - mpreal("5.0")*mpfr::pow(10,-1-precision)))))).toLong()));
        }
    }
    // In case the matrix only contains zeros, NaNs and Infs, the global exponent is 0
    if (largestExponent == INT_MIN)
        largestExponent = 0;
//    cout << "largestExponent = " << largestExponent << endl;


    // Now we print the table according to the following guidelines:
    //  - If the whole matrix is integer and the number of digits needed
    //    to represent the largest integer part is not larger than precision, we
    //    just us the minimal number of digits needed (and don't put any dot)
    //  - If the number of digits needed to print the integer parts of some
    //    numbers in the matrix is larger than the width, we use the
    //    scientific notation
    //  - If some numbers are not integer, we use the full width to print
    //    the numbers, starting with the digit corresponding the largest
    //    number in the table.
    //    The numbers are then completed with zeros within their precision,
    //    with spaces behind.
    // Note that the width that really matters here is width-1, because we
    // keep one character for the sign.

    // Let's determine whether we need to use a global exponent, and fix the
    // display width (which includes the number of digits to print plus a dot (if
    // applicable))
    int exponentShift(0);
    int width(precision);
    if (integerMatrix)
    {
        if (1+largestExponent > width)
        {
            mexPrintf("  1e%i *\n\n",largestExponent);
            exponentShift = largestExponent;
            width = width + 1; // We need one more character to print the dot
        } else {
            width = std::min(width, 1+largestExponent); // We choose the smallest possible width
        }
//        cout << "Width is " << width << endl;
    } else {
        if (1+largestExponent > width-2) // Width minus 2 digits for '.0' digit
        {
            mexPrintf("  1e%i *\n\n",largestExponent);
            exponentShift = largestExponent;
        }
        if (largestExponent < -4) // We allow for the largest number to be just a few orders of magnitude below 0 before switching to the scientific notation
        {
            mexPrintf("  1e%i *\n\n",largestExponent);
            exponentShift = largestExponent;
        }
        width = width + 1; // In all these cases we need one more character to print the dot (which doesn't count as a digit)
    }
//        cout << "width that will be used is = " << width << endl;
//        cout << "exponentShift = " << exponentShift << endl;
//        cout << "isComplex = " << isComplex << endl;

    // We never use less than 3 digits, because we need to be able to print at least
    // NaN or Inf
    width = max(3, width);

    // Now, we decide where to put the dot, knowing that :
    // We will use the scientific notation iff exponentShift is not zero,
    // and no decimals iff (integerMatrix && (exponentShift == 0))
    int dotPosition(width);
    if (integerMatrix && (exponentShift == 0))
    {
        dotPosition = width;
    } else {
        if (exponentShift != 0)
        {
            dotPosition = 1;// i.e. after the sign and first digit (the sign is not counted in the width)
        } else {
            if (largestExponent > 0)
            {
                dotPosition = largestExponent + 1;// i.e. where the dot should be according to the largest number in the table
            } else {
                dotPosition = 1;// Then we will start all numbers by 0.(...)
            }
        }
    }

//        cout << "dotPosition = " << dotPosition << endl;

    // Now we write each component of the matrix according to these rules
    for (mwIndex i = 0; i < matrixR.rows(); ++i) {
        int lineWidth(0);
        for (mwIndex j = 0; j < matrixR.cols(); ++j) {
            mexPrintf("  ");
            if (iszero(matrixR(i,j)))
                mexPrintf("%s", (nSpaces(width)+"0").c_str());
            else if (mpfr::isnan(matrixR(i,j)))
                mexPrintf("%s", (nSpaces(width-2)+"NaN").c_str());
            else if (mpfr::isinf(matrixR(i,j))) {
                if (sgn(matrixR(i,j)) < 0)
                    mexPrintf("%s", (nSpaces(width-3)+"-Inf").c_str());
                else
                    mexPrintf("%s", (nSpaces(width-2)+"Inf").c_str());
            } else
                mexPrintf("%s", mprealToString(matrixR(i,j), width+1, exponentShift, dotPosition+1).c_str());
//                cout << "width + 1 = " << width+1 << endl;
//                cout << "exponentShift = " << exponentShift << endl;
//                cout << "dotPosition + 1 = " << dotPosition+1 << endl;
//                cout << "mprealToString(matrixR(i,j), width+1, exponentShift, dotPosition+1) = " << mprealToString(matrixR(i,j), width+1, exponentShift, dotPosition+1) << endl;
            lineWidth += 2+width+1;
            if (isComplex) {
                if (iszero(matrixI(i,j)))
                    mexPrintf("%s", nSpaces(3+width+1).c_str());
                else if (mpfr::isnan(matrixI(i,j)))
                    mexPrintf("%s", (" + "+nSpaces(width-3)+"NaNi").c_str());
                else if (mpfr::isinf(matrixI(i,j))) {
                    if (sgn(matrixI(i,j)) < 0)
                        mexPrintf("%s", (" - "+nSpaces(width-3)+"Infi").c_str());
                    else
                        mexPrintf("%s", (" + "+nSpaces(width-3)+"Infi").c_str());
                } else {
                    // We get the string for the imaginary part (we don't reserve one character for the "i")
                    string imagText(mprealToString(matrixI(i,j), width, exponentShift, dotPosition, 0));
                    // And now we introduce the i, by replacing the first blank
                    // in the last blanks block by it, or by adding it at the end
                    size_t lastDigit = imagText.find_last_not_of(" ");
                    if ((lastDigit != string::npos) && (lastDigit < width-1)) {
                        imagText.replace(1+lastDigit,1,"i");
                        imagText = imagText + " ";
                    } else {
                        imagText = imagText + "i";
                    }
                    if (matrixI(i,j) < 0)
                        imagText = " - " + imagText;
                    else
                        imagText = " + " + imagText;
                    mexPrintf("%s", imagText.c_str());
                }
                lineWidth += 3+width+1;
            }
            if (lineWidth >= 25000)
            {
                mexPrintf("... Output truncated.  Text exceeds maximum line length of 25,000 characters for Command Window display.");
                break;
            }
        }
        mexPrintf("\n");
    }
}


/*
 Here is another display function, in which each number in the matrix is printed
 independently of the other ones. Here, the width include the sign, digits, dot,
 and exponent if they apply. Each number is processed independently.
*/
void GmpEigenMatrix::displayIndividual(int width) const
{
    // If the matrix is empty, we print nothing
    if (matrixR.rows()*matrixR.cols() == 0)
    {
        mexPrintf("  []\n");
        return;
    }

    if (width < 0) {
        // In this case, we give all possible digits for all numbers

        // First we check which will be the longest number...
        for (mwIndex i = 0; i < matrixR.rows(); ++i) {
            for (mwIndex j = 0; j < matrixR.cols(); ++j) {
                if (iszero(matrixR(i,j))) {}
                else if (mpfr::isnan(matrixR(i,j)))
                    width = max(width, 3);
                else if (mpfr::isinf(matrixR(i,j))) {
                    if (sgn(matrixR(i,j)) < 0)
                        width = max(width, 4);
                    else
                        width = max(width, 3);
                } else
                    width = max(width, (int) matrixR(i,j).toString().length());
                if (isComplex) {
                    if (iszero(matrixI(i,j))) {}
                    else if (mpfr::isnan(matrixI(i,j)))
                        width = max(width, 3);
                    else if (mpfr::isinf(matrixI(i,j))) {
                        if (sgn(matrixI(i,j)) < 0)
                            width = max(width, 3);
                        else
                            width = max(width, 3);
                    } else {
                        width = max(width, (int) matrixI(i,j).toString().length());
                    }
                }
            }
        }

        // Now that we fixed the width, we go ahead and print the matrix
        for (mwIndex i = 0; i < matrixR.rows(); ++i) {
            int lineWidth(0);
            for (mwIndex j = 0; j < matrixR.cols(); ++j) {
                mexPrintf("  ");
                if (iszero(matrixR(i,j)))
                    mexPrintf("%s", (nSpaces(width-1)+"0").c_str());
                else if (mpfr::isnan(matrixR(i,j)))
                    mexPrintf("%s", (nSpaces(width-3)+"NaN").c_str());
                else if (mpfr::isinf(matrixR(i,j))) {
                    if (sgn(matrixR(i,j)) < 0)
                        mexPrintf("%s", (nSpaces(width-4)+"-Inf").c_str());
                    else
                        mexPrintf("%s", (nSpaces(width-3)+"Inf").c_str());
                } else {
                    string chain(matrixR(i,j).toString());
                    mexPrintf("%s", (nSpaces(width-chain.length())+chain).c_str());
                }
                lineWidth += 2+width;
                if (isComplex) {
                    if (iszero(matrixI(i,j)))
                        mexPrintf("%s", nSpaces(3+width+1).c_str());
                    else if (mpfr::isnan(matrixI(i,j)))
                        mexPrintf("%s", (" + "+nSpaces(width-3)+"NaNi").c_str());
                    else if (mpfr::isinf(matrixI(i,j))) {
                        if (sgn(matrixI(i,j)) < 0)
                            mexPrintf("%s", (" - "+nSpaces(width-3)+"Infi").c_str());
                        else
                            mexPrintf("%s", (" + "+nSpaces(width-3)+"Infi").c_str());
                    } else {
                        string chain(matrixI(i,j).toString()+"i");
                        if (matrixI(i,j) < 0)
                            chain = " - " + chain.substr(1,chain.length()-1);
                        else
                            chain = " + " + chain;
                        mexPrintf("%s", (chain+nSpaces(4+width-chain.length())).c_str());
                    }
                    lineWidth += 3+width+1;
                }
                if (lineWidth >= 25000)
                {
                    mexPrintf("... Output truncated.  Text exceeds maximum line length of 25,000 characters for Command Window display.");
                    break;
                }
            }
            mexPrintf("\n");
        }

        return;
    }


    // WARNING : This part of the function uses a function which may not work
    //           properly in presence of rounding(!!!)

    // We make sure the width is not too small to avoid problems
    if (width < 5)
        width = 5;

    // We write each component of the matrix independently
    for (mwIndex i = 0; i < matrixR.rows(); ++i) {
        for (mwIndex j = 0; j < matrixR.cols(); ++j) {
            mexPrintf("  ");
            if (iszero(matrixR(i,j)))
                mexPrintf("%s", (nSpaces(width)+"0").c_str());
            else if (mpfr::isnan(matrixR(i,j)))
                mexPrintf("%s", (nSpaces(width-2)+"NaN").c_str());
            else if (mpfr::isinf(matrixR(i,j))) {
                if (sgn(matrixR(i,j)) < 0)
                    mexPrintf("%s", (nSpaces(width-3)+"-Inf").c_str());
                else
                    mexPrintf("%s", (nSpaces(width-2)+"Inf").c_str());
            } else
                mexPrintf("%s", mprealToString2(matrixR(i,j), width, false).c_str());
            if (isComplex) {
                if (iszero(matrixI(i,j)))
                    mexPrintf("%s", nSpaces(2+width+1).c_str());
                else if (mpfr::isnan(matrixI(i,j)))
                    mexPrintf("%s", (" + "+nSpaces(width-4)+"NaNi").c_str());
                else if (mpfr::isinf(matrixI(i,j))) {
                    if (sgn(matrixI(i,j)) < 0)
                        mexPrintf("%s", (" - "+nSpaces(width-4)+"Infi").c_str());
                    else
                        mexPrintf("%s", (nSpaces(width-1)+"Infi").c_str());
                } else {
                    // We get the string for the imaginary part (we don't reserve one character for the "i")
                    string imagText(mprealToString2(matrixI(i,j), width, true));
                    // And now we introduce the i, by replacing the first blank
                    // in the last blanks block by it, or by adding it at the end
                    size_t lastDigit = imagText.find_last_not_of(" ");
                    if ((lastDigit != string::npos) && (lastDigit < width-1)) {
                        imagText.replace(1+lastDigit,1,"i");
                        imagText = imagText + " ";
                    } else {
                        imagText = imagText + "i";
                    }
                    if (matrixI(i,j) < 0)
                        imagText = " - " + imagText.substr(1,imagText.length()-1);
                    else
                        imagText = " + " + imagText.substr(1,imagText.length()-1);
                    mexPrintf("%s", imagText.c_str());
                }
            }
        }
        mexPrintf("\n");
    }
}





/* Extracting the raw data: the following function extracts the internal
   data from the class, and puts it in a basic c-type table to be accessed
   by matlab */
mxArray* GmpEigenMatrix::toDouble() const
{
    /* Get the size and pointers to input data */
    mwSize m(matrixR.rows());
    mwSize n(matrixR.cols());

    mxArray* plhs;
    if (isComplex)
        plhs = mxCreateNumericMatrix(m, n, mxDOUBLE_CLASS, mxCOMPLEX);
    else
        plhs = mxCreateNumericMatrix(m, n, mxDOUBLE_CLASS, mxREAL);

    // Now we copy the data over (actually, only one item should be copied
    // here, but we keep the loop for when we'll like to copy tables)
    double* pointerR(mxGetPr(plhs));
    double* pointerI(mxGetPi(plhs));
    for (mwIndex j = 0; j < n; ++j) {
        // We first iterate on the columns
        for (mwIndex i = 0; i < m; ++i) {
            pointerR[i] = matrixR(i,j).toDouble();
            if (isComplex) {
                pointerI[i] = matrixI(i,j).toDouble();
            }
        }
        pointerR += m;
        pointerI += m;
    }
    return plhs;
}


/* This function returns a cell array with strings describing each matrix
   element */
mxArray* GmpEigenMatrix::toStrings(const int& precision) const
{
    // Now, we add each field
    // We start with the real part
    mwSize m(matrixR.rows());
    mwSize n(matrixR.cols());
    mxArray* cellArray(mxCreateCellMatrix(m, n));

    // Now we iterate on all the elements of the cell array
    mwIndex index(0);
    for (mwIndex j = 0; j < n; ++j) {
        // We first iterate on the columns
        for (mwIndex i = 0; i < m; ++i) {
            string text(matrixR(i,j).toString(precision));
            if ((isComplex) && (matrixI(i,j) != 0))
            {
                if (matrixR(i,j)==0)
                    text = "";
                else if (sgn(matrixI(i,j)) >= 0)
                    text += "+";
                text += matrixI(i,j).toString(precision) + "i";
            }
            mxSetCell(cellArray, index, mxCreateString(text.c_str()));
            ++index;
        }
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

/* This function returns the number of significant digits in matlab format */
mxArray* GmpEigenMatrix::precision() const
{
    /* Get the size and pointers to input data */
    mwSize m(matrixR.rows());
    mwSize n(matrixR.cols());

    mxArray* plhs;
    plhs = mxCreateNumericMatrix(m, n, mxDOUBLE_CLASS, mxREAL);

    // Now we copy the precision to the matlab table
    double* pointerR(mxGetPr(plhs));
    for (mwIndex j = 0; j < n; ++j) {
        // We first iterate on the columns
        for (mwIndex i = 0; i < m; ++i) {
            pointerR[i] = mpfr::bits2digits(matrixR(i,j).get_prec());
        }
        pointerR += m;
    }

    return plhs;
}


/* This function extracts the nonzero elements from the matrix */
GmpEigenMatrix GmpEigenMatrix::find(vector < double >& rows, vector < double >& cols) const
{
    // We start by filling dynamically allocated tables (we don't know exactly
    // how many elements we'll end up with)
    vector < vector < IndexType > > items;
    for (IndexType j = 0; j < matrixR.cols(); ++j) {
        for (IndexType i = 0; i < matrixR.rows(); ++i) {
            if ((matrixR(i,j) != 0) || (isComplex && (matrixI(i,j) != 0))) {
                rows.push_back(i);
                cols.push_back(j);
                items.push_back(std::vector< IndexType > (1, matrixR.rows()*j + i));
            }
        }
    }

    // Now we extract the corresponding non-zero values
    return subsref(items);
}


/* This function extracts the nonzero elements from the matrix */
GmpEigenMatrix& GmpEigenMatrix::find_new(vector < double >& rows, vector < double >& cols) const
{
    // We start by filling dynamically allocated tables (we don't know exactly
    // how many elements we'll end up with)
    vector < vector < IndexType > > items;
    for (IndexType j = 0; j < matrixR.cols(); ++j) {
        for (IndexType i = 0; i < matrixR.rows(); ++i) {
            if ((matrixR(i,j) != 0) || (isComplex && (matrixI(i,j) != 0))) {
                rows.push_back(i);
                cols.push_back(j);
                items.push_back(std::vector< IndexType > (1, matrixR.rows()*j + i));
            }
        }
    }

    // Now we extract the corresponding non-zero values
    return subsref_new(items);
}


/* This function returns the size in matlab format */
mxArray* GmpEigenMatrix::size() const
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
GmpEigenMatrix GmpEigenMatrix::reshape(const IndexType& m2, const IndexType& n2) const
{
    GmpEigenMatrix result;

    result.isComplex = isComplex;
    result.matrixR.resize(m2,n2);
    if (isComplex)
        result.matrixI.resize(m2,n2);

    IndexType m(matrixR.rows());
    IndexType n(matrixR.cols());
    IndexType i(0), j(0), i2(0), j2(0);
    for (IndexType x(0); x < m*n ; ++x) {
        result.matrixR(i2,j2) = matrixR(i,j);
        if (isComplex)
            result.matrixI(i2,j2) = matrixI(i,j);
        ++i;
        if (i == m) {i = 0; j += 1;}
        ++i2;
        if (i2 == m2) {i2 = 0; j2 += 1;}
    }

    return result;
}

/* Reshaping a matrix b = reshape(a, [m2, n2]) */
GmpEigenMatrix& GmpEigenMatrix::reshape_new(const IndexType& m2, const IndexType& n2) const
{
    GmpEigenMatrix& result(*(new GmpEigenMatrix));

    result.isComplex = isComplex;
    result.matrixR.resize(m2,n2);
    if (isComplex)
        result.matrixI.resize(m2,n2);

    IndexType m(matrixR.rows());
    IndexType n(matrixR.cols());
    IndexType i(0), j(0), i2(0), j2(0);
    for (IndexType x(0); x < m*n ; ++x) {
        result.matrixR(i2,j2) = matrixR(i,j);
        if (isComplex)
            result.matrixI(i2,j2) = matrixI(i,j);
        ++i;
        if (i == m) {i = 0; j += 1;}
        ++i2;
        if (i2 == m2) {i2 = 0; j2 += 1;}
    }

    return result;
}


/* Vertical concatenation of two matrices c = [a; b]*/
GmpEigenMatrix GmpEigenMatrix::vertcat(const GmpEigenMatrix& b) const
{
    GmpEigenMatrix result(*this);

    result.matrixR.conservativeResize(matrixR.rows()+b.matrixR.rows(), matrixR.cols());
    result.matrixR.block(matrixR.rows(), 0, b.matrixR.rows(), b.matrixR.cols()) = b.matrixR;
    if ((isComplex) || (b.isComplex)) {
        result.isComplex = true;
        result.matrixI.conservativeResize(matrixR.rows()+b.matrixR.rows(), matrixR.cols());
        if (b.isComplex)
            result.matrixI.block(matrixR.rows(), 0, b.matrixR.rows(), b.matrixR.cols()) = b.matrixI;
    }

    return result;
}

/* Vertical concatenation of two matrices c = [a; b]*/
GmpEigenMatrix& GmpEigenMatrix::vertcat_new(const GmpEigenMatrix& b) const
{
    GmpEigenMatrix& result(*(new GmpEigenMatrix(*this)));

    result.matrixR.conservativeResize(matrixR.rows()+b.matrixR.rows(), matrixR.cols());
    result.matrixR.block(matrixR.rows(), 0, b.matrixR.rows(), b.matrixR.cols()) = b.matrixR;
    if ((isComplex) || (b.isComplex)) {
        result.isComplex = true;
        result.matrixI.conservativeResize(matrixR.rows()+b.matrixR.rows(), matrixR.cols());
        if (b.isComplex)
            result.matrixI.block(matrixR.rows(), 0, b.matrixR.rows(), b.matrixR.cols()) = b.matrixI;
    }

    return result;
}

/* Horizontal concatenation of two matrices c = [a, b]*/
GmpEigenMatrix GmpEigenMatrix::horzcat(const GmpEigenMatrix& b) const
{
    GmpEigenMatrix result(*this);

    result.matrixR.conservativeResize(matrixR.rows(), matrixR.cols()+b.matrixR.cols());
    result.matrixR.block(0, matrixR.cols(), b.matrixR.rows(), b.matrixR.cols()) = b.matrixR;
    if ((isComplex) || (b.isComplex)) {
        result.isComplex = true;
        result.matrixI.conservativeResize(matrixR.rows(), matrixR.cols()+b.matrixR.cols());
        if (b.isComplex)
            result.matrixI.block(0, matrixR.cols(), b.matrixR.rows(), b.matrixR.cols()) = b.matrixI;
    }

    return result;
}

/* Horizontal concatenation of two matrices c = [a, b]*/
GmpEigenMatrix& GmpEigenMatrix::horzcat_new(const GmpEigenMatrix& b) const
{
    GmpEigenMatrix& result(*(new GmpEigenMatrix(*this)));

    result.matrixR.conservativeResize(matrixR.rows(), matrixR.cols()+b.matrixR.cols());
    result.matrixR.block(0, matrixR.cols(), b.matrixR.rows(), b.matrixR.cols()) = b.matrixR;
    if ((isComplex) || (b.isComplex)) {
        result.isComplex = true;
        result.matrixI.conservativeResize(matrixR.rows(), matrixR.cols()+b.matrixR.cols());
        if (b.isComplex)
            result.matrixI.block(0, matrixR.cols(), b.matrixR.rows(), b.matrixR.cols()) = b.matrixI;
    }

    return result;
}




/* This function extracts a submatrix or size row x cols fomr the position (i,j) */
GmpEigenMatrix GmpEigenMatrix::block(const IndexType& i, const IndexType& j, const IndexType& rows, const IndexType& cols) const
{
    GmpEigenMatrix result;

    result.matrixR = matrixR.block(i,j,rows,cols);
    if (isComplex)
        result.matrixI = matrixI.block(i,j,rows,cols);

    result.checkComplexity();

    return result;
}

/* This function extracts a sub-matrix with single set of index b = a([1,2,3])
   or b = a([1 2; 2 2]) */
GmpEigenMatrix GmpEigenMatrix::subsref(const vector<vector<IndexType> >& indices) const
{
    GmpEigenMatrix result;

    result.matrixR.resize(indices.size(), indices[0].size());
    if (isComplex)
        result.matrixI.resize(indices.size(), indices[0].size());

    for (IndexType x(0); x < indices.size(); ++x) {
        for (IndexType y(0); y < indices[x].size(); ++y) {
            // We separate the number indices into i and j
            IndexType i, j;
            j = indices[x][y] / matrixR.rows();
            i = indices[x][y] % matrixR.rows();

            result.matrixR(x,y) = matrixR(i,j);
            if (isComplex)
                result.matrixI(x,y) = matrixI(i,j);
        }
    }

    result.checkComplexity();

    return result;
}

/* This function extracts a sub-matrix with single set of index b = a([1,2,3])
   or b = a([1 2; 2 2]) */
GmpEigenMatrix& GmpEigenMatrix::subsref_new(const vector<vector<IndexType> >& indices) const
{
    GmpEigenMatrix& result(*(new GmpEigenMatrix));

    result.matrixR.resize(indices.size(), indices[0].size());
    if (isComplex)
        result.matrixI.resize(indices.size(), indices[0].size());

    for (IndexType x(0); x < indices.size(); ++x) {
        for (IndexType y(0); y < indices[x].size(); ++y) {
            // We separate the number indices into i and j
            IndexType i, j;
            j = indices[x][y] / matrixR.rows();
            i = indices[x][y] % matrixR.rows();

            result.matrixR(x,y) = matrixR(i,j);
            if (isComplex)
                result.matrixI(x,y) = matrixI(i,j);
        }
    }

    result.checkComplexity();

    return result;
}

/* This function extracts a sub-matrix with two sets of indices b = a(1:2,2) */
GmpEigenMatrix GmpEigenMatrix::subsref(const vector<IndexType>& indicesA, const vector<IndexType>& indicesB) const
{
    GmpEigenMatrix result;

    result.matrixR.resize(indicesA.size(), indicesB.size());
    if (isComplex)
        result.matrixI.resize(indicesA.size(), indicesB.size());

    for (IndexType i(0); i < indicesA.size(); ++i)
        for (IndexType j(0); j < indicesB.size(); ++j) {
            result.matrixR(i,j) = matrixR(indicesA[i], indicesB[j]);
            if (isComplex)
                result.matrixI(i,j) = matrixI(indicesA[i], indicesB[j]);
        }

    result.checkComplexity();

    return result;
}

/* This function extracts a sub-matrix with two sets of indices b = a(1:2,2) */
GmpEigenMatrix& GmpEigenMatrix::subsref_new(const vector<IndexType>& indicesA, const vector<IndexType>& indicesB) const
{
    GmpEigenMatrix& result(*(new GmpEigenMatrix));

    result.matrixR.resize(indicesA.size(), indicesB.size());
    if (isComplex)
        result.matrixI.resize(indicesA.size(), indicesB.size());

    for (IndexType i(0); i < indicesA.size(); ++i)
        for (IndexType j(0); j < indicesB.size(); ++j) {
            result.matrixR(i,j) = matrixR(indicesA[i], indicesB[j]);
            if (isComplex)
                result.matrixI(i,j) = matrixI(indicesA[i], indicesB[j]);
        }

    result.checkComplexity();

    return result;
}

/* Assigning value to a sub-matrix a([1 2 3]) = b */
void GmpEigenMatrix::subsasgn(const std::vector<IndexType>& indices, const GmpEigenMatrix& b)
{
    if ((!isComplex) && (b.isComplex)) {
        isComplex = true;
        matrixI = zeroMatrix(matrixR.rows(), matrixR.cols());
    }

    for (IndexType x(0); x < indices.size(); ++x) {
        // We separate the number indices into i and j
        IndexType i, j;
        j = indices[x] / matrixR.rows();
        i = indices[x] % matrixR.rows();

        matrixR(i,j) = b.matrixR(x,0);
        if (b.isComplex)
            matrixI(i,j) = b.matrixI(x,0);
        else if (isComplex)
            matrixI(i,j) = "0";
    }

    if (isComplex)
        checkComplexity();
}

/* Assigning value to a sub-matrix a(1:2,:) = b */
void GmpEigenMatrix::subsasgn(const std::vector<IndexType>& indicesA, const std::vector<IndexType>& indicesB, const GmpEigenMatrix& b)
{
    if ((!isComplex) && (b.isComplex)) {
        isComplex = true;
        matrixI = zeroMatrix(matrixR.rows(), matrixR.cols());
    }

    for (IndexType i(0); i < indicesA.size(); ++i)
        for (IndexType j(0); j < indicesB.size(); ++j) {
            matrixR(indicesA[i], indicesB[j]) = b.matrixR(i,j);
            if (b.isComplex)
                matrixI(indicesA[i], indicesB[j]) = b.matrixI(i,j);
            else if (isComplex)
                matrixI(indicesA[i], indicesB[j]) = "0";
        }

    if (isComplex)
        checkComplexity();
}

/* Creation of a diagonal matrix b = diag(a, k) */
GmpEigenMatrix GmpEigenMatrix::diagCreate(const IndexType& k) const
{
    GmpEigenMatrix result;

    result.isComplex = isComplex;
    result.matrixR.resize(matrixR.rows()+std::abs(k), matrixR.rows()+std::abs(k));
    if (isComplex)
        result.matrixI.resize(matrixR.rows()+std::abs(k), matrixR.rows()+std::abs(k));

    for (IndexType i(0); i < matrixR.rows(); ++i) {
        result.matrixR(max((IndexType)0,-k)+i, max((IndexType)0,k)+i) = matrixR(i, 0);
        if (isComplex)
            result.matrixI(max((IndexType)0,-k)+i, max((IndexType)0,k)+i) = matrixI(i, 0);
    }

    return result;
}

/* Creation of a diagonal matrix b = diag(a, k) */
GmpEigenMatrix& GmpEigenMatrix::diagCreate_new(const IndexType& k) const
{
    GmpEigenMatrix& result(*(new GmpEigenMatrix));

    result.isComplex = isComplex;
    result.matrixR.resize(matrixR.rows()+std::abs(k), matrixR.rows()+std::abs(k));
    if (isComplex)
        result.matrixI.resize(matrixR.rows()+std::abs(k), matrixR.rows()+std::abs(k));

    for (IndexType i(0); i < matrixR.rows(); ++i) {
        result.matrixR(max((IndexType)0,-k)+i, max((IndexType)0,k)+i) = matrixR(i, 0);
        if (isComplex)
            result.matrixI(max((IndexType)0,-k)+i, max((IndexType)0,k)+i) = matrixI(i, 0);
    }

    return result;
}

/* Extraction of diagonal elements b = diag(a, k) */
GmpEigenMatrix GmpEigenMatrix::diagExtract(const IndexType& k) const
{
    GmpEigenMatrix result;

    result.matrixR = matrixR.block(max((IndexType)0,-k), max((IndexType)0,k), matrixR.rows()-max((IndexType)0,-k), matrixR.cols()-max((IndexType)0,k)).diagonal();
    if (isComplex)
        result.matrixI = matrixI.block(max((IndexType)0,-k), max((IndexType)0,k), matrixR.rows()-max((IndexType)0,-k), matrixR.cols()-max((IndexType)0,k)).diagonal();
    result.checkComplexity();

    return result;
}

/* Extraction of diagonal elements b = diag(a, k) */
GmpEigenMatrix& GmpEigenMatrix::diagExtract_new(const IndexType& k) const
{
    GmpEigenMatrix& result(*(new GmpEigenMatrix));

    result.matrixR = matrixR.block(max((IndexType)0,-k), max((IndexType)0,k), matrixR.rows()-max((IndexType)0,-k), matrixR.cols()-max((IndexType)0,k)).diagonal();
    if (isComplex)
        result.matrixI = matrixI.block(max((IndexType)0,-k), max((IndexType)0,k), matrixR.rows()-max((IndexType)0,-k), matrixR.cols()-max((IndexType)0,k)).diagonal();
    result.checkComplexity();

    return result;
}







/* Some class operators */



/* The unary substraction operators b = -a */
GmpEigenMatrix GmpEigenMatrix::operator-() const
{
    GmpEigenMatrix result;

    result.matrixR = -matrixR;
    result.isComplex = isComplex;
    if (isComplex)
        result.matrixI = -matrixI;

    return result;
}

/* The unary substraction b = -a */
GmpEigenMatrix& GmpEigenMatrix::uminus_new() const
{
    GmpEigenMatrix& result(*(new GmpEigenMatrix));

    result.matrixR = -matrixR;
    result.isComplex = isComplex;
    if (isComplex)
        result.matrixI = -matrixI;

    return result;
}

/* The tansposition b = a' */
GmpEigenMatrix GmpEigenMatrix::transpose() const
{
    GmpEigenMatrix result;

    result.isComplex = isComplex;
    result.matrixR = matrixR.transpose();
    if (isComplex)
        result.matrixI = matrixI.transpose();

    return result;
}

/* The tansposition b = a' */
GmpEigenMatrix& GmpEigenMatrix::transpose_new() const
{
    GmpEigenMatrix& result(*(new GmpEigenMatrix));

    result.isComplex = isComplex;
    result.matrixR = matrixR.transpose();
    if (isComplex)
        result.matrixI = matrixI.transpose();

    return result;
}

/* The conjugate transposition b = a' */
GmpEigenMatrix GmpEigenMatrix::ctranspose() const
{
    GmpEigenMatrix result;

    result.isComplex = isComplex;
    result.matrixR = matrixR.transpose();
    if (isComplex)
        result.matrixI = -matrixI.transpose();

    return result;
}

/* The conjugate transposition b = a' */
GmpEigenMatrix& GmpEigenMatrix::ctranspose_new() const
{
    GmpEigenMatrix& result(*(new GmpEigenMatrix));

    result.isComplex = isComplex;
    result.matrixR = matrixR.transpose();
    if (isComplex)
        result.matrixI = -matrixI.transpose();

    return result;
}


/* Element-wise complex conjugation b = conj(a) */
GmpEigenMatrix GmpEigenMatrix::conj() const
{
    GmpEigenMatrix result;

    result.isComplex = isComplex;
    result.matrixR = matrixR;
    if (isComplex)
        result.matrixI = -matrixI;

    return result;
}

/* Element-wise complex conjugation b = conj(a) */
GmpEigenMatrix& GmpEigenMatrix::conj_new() const
{
    GmpEigenMatrix& result(*(new GmpEigenMatrix));

    result.isComplex = isComplex;
    result.matrixR = matrixR;
    if (isComplex)
        result.matrixI = -matrixI;

    return result;
}

/* Real part b = real(a) */
GmpEigenMatrix GmpEigenMatrix::real() const
{
    GmpEigenMatrix result;

    result.isComplex = false;
    result.matrixR = matrixR;

    return result;
}

/* Real part b = real(a) */
GmpEigenMatrix& GmpEigenMatrix::real_new() const
{
    GmpEigenMatrix& result(*(new GmpEigenMatrix));

    result.isComplex = false;
    result.matrixR = matrixR;

    return result;
}

/* Imaginary part b = imag(a) */
GmpEigenMatrix GmpEigenMatrix::imag() const
{
    GmpEigenMatrix result;

    result.isComplex = false;
    if (isComplex)
        result.matrixR = matrixI;
    else
        result.matrixR = zeroMatrix(matrixR.rows(), matrixR.cols(), mpfr::mpreal::get_default_prec());

    return result;
}

/* Imaginary part b = imag(a) */
GmpEigenMatrix& GmpEigenMatrix::imag_new() const
{
    GmpEigenMatrix& result(*(new GmpEigenMatrix));

    result.isComplex = false;
    if (isComplex)
        result.matrixR = matrixI;
    else
        result.matrixR = zeroMatrix(matrixR.rows(), matrixR.cols(), mpfr::mpreal::get_default_prec());

    return result;
}


// Various integer parts
/* Rounds to nearest integer b = round(a) */
GmpEigenMatrix GmpEigenMatrix::round() const
{
    GmpEigenMatrix result;

    result.matrixR.resize(matrixR.rows(), matrixR.cols());
    if (isComplex)
        result.matrixI.resize(matrixR.rows(), matrixR.cols());

    for (IndexType j = 0; j < matrixR.cols(); ++j) {
        for (IndexType i = 0; i < matrixR.rows(); ++i) {
            result.matrixR(i,j) = mpfr::round(matrixR(i,j));
            if (isComplex)
                result.matrixI(i,j) = mpfr::round(matrixI(i,j));
        }
    }

    result.checkComplexity();

    return result;
}

/* Rounds to nearest integer b = round(a) */
GmpEigenMatrix& GmpEigenMatrix::round_new() const
{
    GmpEigenMatrix& result(*(new GmpEigenMatrix));

    result.matrixR.resize(matrixR.rows(), matrixR.cols());
    if (isComplex)
        result.matrixI.resize(matrixR.rows(), matrixR.cols());

    for (IndexType j = 0; j < matrixR.cols(); ++j) {
        for (IndexType i = 0; i < matrixR.rows(); ++i) {
            result.matrixR(i,j) = mpfr::round(matrixR(i,j));
            if (isComplex)
                result.matrixI(i,j) = mpfr::round(matrixI(i,j));
        }
    }

    result.checkComplexity();

    return result;
}

/* Next integer towards -Inf b = floor(a) */
GmpEigenMatrix GmpEigenMatrix::floor() const
{
    GmpEigenMatrix result;

    result.matrixR.resize(matrixR.rows(), matrixR.cols());
    if (isComplex)
        result.matrixI.resize(matrixR.rows(), matrixR.cols());

    for (IndexType j = 0; j < matrixR.cols(); ++j) {
        for (IndexType i = 0; i < matrixR.rows(); ++i) {
            result.matrixR(i,j) = mpfr::floor(matrixR(i,j));
            if (isComplex)
                result.matrixI(i,j) = mpfr::floor(matrixI(i,j));
        }
    }

    result.checkComplexity();

    return result;
}

/* Next integer towards -Inf b = floor(a) */
GmpEigenMatrix& GmpEigenMatrix::floor_new() const
{
    GmpEigenMatrix& result(*(new GmpEigenMatrix));

    result.matrixR.resize(matrixR.rows(), matrixR.cols());
    if (isComplex)
        result.matrixI.resize(matrixR.rows(), matrixR.cols());

    for (IndexType j = 0; j < matrixR.cols(); ++j) {
        for (IndexType i = 0; i < matrixR.rows(); ++i) {
            result.matrixR(i,j) = mpfr::floor(matrixR(i,j));
            if (isComplex)
                result.matrixI(i,j) = mpfr::floor(matrixI(i,j));
        }
    }

    result.checkComplexity();

    return result;
}

/* Next integer towards +Inf b = ceil(a) */
GmpEigenMatrix GmpEigenMatrix::ceil() const
{
    GmpEigenMatrix result;

    result.matrixR.resize(matrixR.rows(), matrixR.cols());
    if (isComplex)
        result.matrixI.resize(matrixR.rows(), matrixR.cols());

    for (IndexType j = 0; j < matrixR.cols(); ++j) {
        for (IndexType i = 0; i < matrixR.rows(); ++i) {
            result.matrixR(i,j) = mpfr::ceil(matrixR(i,j));
            if (isComplex)
                result.matrixI(i,j) = mpfr::ceil(matrixI(i,j));
        }
    }

    result.checkComplexity();

    return result;
}

/* Next integer towards +Inf b = ceil(a) */
GmpEigenMatrix& GmpEigenMatrix::ceil_new() const
{
    GmpEigenMatrix& result(*(new GmpEigenMatrix));

    result.matrixR.resize(matrixR.rows(), matrixR.cols());
    if (isComplex)
        result.matrixI.resize(matrixR.rows(), matrixR.cols());

    for (IndexType j = 0; j < matrixR.cols(); ++j) {
        for (IndexType i = 0; i < matrixR.rows(); ++i) {
            result.matrixR(i,j) = mpfr::ceil(matrixR(i,j));
            if (isComplex)
                result.matrixI(i,j) = mpfr::ceil(matrixI(i,j));
        }
    }

    result.checkComplexity();

    return result;
}

/* Next integer towards 0 b = trunc(a) */
GmpEigenMatrix GmpEigenMatrix::trunc() const
{
    GmpEigenMatrix result;

    result.matrixR.resize(matrixR.rows(), matrixR.cols());
    if (isComplex)
        result.matrixI.resize(matrixR.rows(), matrixR.cols());

    for (IndexType j = 0; j < matrixR.cols(); ++j) {
        for (IndexType i = 0; i < matrixR.rows(); ++i) {
            result.matrixR(i,j) = mpfr::trunc(matrixR(i,j));
            if (isComplex)
                result.matrixI(i,j) = mpfr::trunc(matrixI(i,j));
        }
    }

    result.checkComplexity();

    return result;
}

/* Next integer towards 0 b = trunc(a) */
GmpEigenMatrix& GmpEigenMatrix::trunc_new() const
{
    GmpEigenMatrix& result(*(new GmpEigenMatrix));

    result.matrixR.resize(matrixR.rows(), matrixR.cols());
    if (isComplex)
        result.matrixI.resize(matrixR.rows(), matrixR.cols());

    for (IndexType j = 0; j < matrixR.cols(); ++j) {
        for (IndexType i = 0; i < matrixR.rows(); ++i) {
            result.matrixR(i,j) = mpfr::trunc(matrixR(i,j));
            if (isComplex)
                result.matrixI(i,j) = mpfr::trunc(matrixI(i,j));
        }
    }

    result.checkComplexity();

    return result;
}


// Concerning the addition operators:
// We want to understand addition by a 1x1 matrix as scalar addition,
// which is not natively supported by Eigen, so we adjust the addition
// procedure by hand.
GmpEigenMatrix& GmpEigenMatrix::operator+=(const GmpEigenMatrix& b)
{
    if ((numel() != 1) && (b.numel() == 1)) {
        matrixR = matrixR.array() + b.matrixR(0,0);
        if (b.isComplex) {
            if (isComplex) {
                matrixI = matrixI.array() + b.matrixI(0,0);
                checkComplexity();
            } else {
                isComplex = true;
                matrixI = constantMatrix(matrixR.rows(),matrixR.cols(), b.matrixI(0,0));
            }
        }
    } else if ((numel() == 1) && (b.numel() != 1)) {
        matrixR = matrixR(0,0) + b.matrixR.array();
        if (b.isComplex) {
            if (isComplex) {
                matrixI = matrixI(0,0) + b.matrixI.array();
                checkComplexity();
            } else {
                isComplex = true;
                matrixI = b.matrixI;
            }
        } else if (isComplex) {
            matrixI = constantMatrix(matrixR.rows(),matrixR.cols(), matrixI(0,0));
        }
    } else {
        matrixR += b.matrixR;
        if (b.isComplex) {
            if (isComplex) {
                matrixI += b.matrixI;
                checkComplexity();
            } else {
                isComplex = true;
                matrixI = b.matrixI;
            }
        }
    }

    return *this;
}

/* The addition operator c = a + b */
GmpEigenMatrix GmpEigenMatrix::operator+(const GmpEigenMatrix& b) const
{
    GmpEigenMatrix result;

    result.isComplex = isComplex;
    if ((numel() != 1) && (b.numel() == 1)) {
        result.matrixR = matrixR.array() + b.matrixR(0,0);
        if (b.isComplex) {
            if (isComplex) {
                result.matrixI = matrixI.array() + b.matrixI(0,0);
                result.checkComplexity();
            } else {
                result.isComplex = true;
                result.matrixI = constantMatrix(matrixR.rows(),matrixR.cols(), b.matrixI(0,0));
            }
        } else if (isComplex) {
            result.matrixI = matrixI;
        }
    } else if ((numel() == 1) && (b.numel() != 1)) {
        result.matrixR = matrixR(0,0) + b.matrixR.array();
        if (b.isComplex) {
            if (isComplex) {
                result.matrixI = matrixI(0,0) + b.matrixI.array();
                result.checkComplexity();
            } else {
                result.isComplex = true;
                result.matrixI = b.matrixI;
            }
        } else if (isComplex) {
            result.matrixI = constantMatrix(matrixR.rows(),matrixR.cols(), matrixI(0,0));
        }
    } else {
        result.matrixR = matrixR + b.matrixR;
        if (b.isComplex) {
            if (isComplex) {
                result.matrixI = matrixI + b.matrixI;
                result.checkComplexity();
            } else {
                result.isComplex = true;
                result.matrixI = b.matrixI;
            }
        } else if (isComplex) {
            result.matrixI = matrixI;
        }
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
GmpEigenMatrix& GmpEigenMatrix::plus_new(const GmpEigenMatrix& b) const
{
    GmpEigenMatrix& result(*(new GmpEigenMatrix));

    result.isComplex = isComplex;
    if ((numel() != 1) && (b.numel() == 1)) {
        result.matrixR = matrixR.array() + b.matrixR(0,0);
        if (b.isComplex) {
            if (isComplex) {
                result.matrixI = matrixI.array() + b.matrixI(0,0);
                result.checkComplexity();
            } else {
                result.isComplex = true;
                result.matrixI = constantMatrix(matrixR.rows(),matrixR.cols(), b.matrixI(0,0));
            }
        } else if (isComplex) {
            result.matrixI = matrixI;
        }
    } else if ((numel() == 1) && (b.numel() != 1)) {
        result.matrixR = matrixR(0,0) + b.matrixR.array();
        if (b.isComplex) {
            if (isComplex) {
                result.matrixI = matrixI(0,0) + b.matrixI.array();
                result.checkComplexity();
            } else {
                result.isComplex = true;
                result.matrixI = b.matrixI;
            }
        } else if (isComplex) {
            result.matrixI = constantMatrix(matrixR.rows(),matrixR.cols(), matrixI(0,0));
        }
    } else {
        result.matrixR = matrixR + b.matrixR;
        if (b.isComplex) {
            if (isComplex) {
                result.matrixI = matrixI + b.matrixI;
                result.checkComplexity();
            } else {
                result.isComplex = true;
                result.matrixI = b.matrixI;
            }
        } else if (isComplex) {
            result.matrixI = matrixI;
        }
    }

    return result;
}

/* The substraction operator a -= b */
GmpEigenMatrix& GmpEigenMatrix::operator-=(const GmpEigenMatrix& b)
{
    if ((numel() != 1) && (b.numel() == 1)) {
        matrixR = matrixR.array() - b.matrixR(0,0);
        if (b.isComplex) {
            if (isComplex) {
                matrixI = matrixI.array() - b.matrixI(0,0);
                checkComplexity();
            } else {
                isComplex = true;
                matrixI = constantMatrix(matrixR.rows(),matrixR.cols(), -b.matrixI(0,0));
            }
        }
    } else if ((numel() == 1) && (b.numel() != 1)) {
        matrixR = matrixR(0,0) - b.matrixR.array();
        if (b.isComplex) {
            if (isComplex) {
                matrixI = matrixI(0,0) - b.matrixI.array();
                checkComplexity();
            } else {
                isComplex = true;
                matrixI = b.matrixI;
            }
        } else if (isComplex) {
            matrixI = constantMatrix(matrixR.rows(),matrixR.cols(), matrixI(0,0));
        }
    } else {
        matrixR -= b.matrixR;
        if (b.isComplex) {
            if (isComplex) {
                matrixI -= b.matrixI;
                checkComplexity();
            } else {
                isComplex = true;
                matrixI = -b.matrixI;
            }
        }
    }

    return *this;
}

/* The substraction operator c = a-b */
GmpEigenMatrix GmpEigenMatrix::operator-(const GmpEigenMatrix& b) const
{
    GmpEigenMatrix result;

    result.isComplex = isComplex;
    if ((numel() != 1) && (b.numel() == 1)) {
        result.matrixR = matrixR.array() - b.matrixR(0,0);
        if (b.isComplex) {
            if (isComplex) {
                result.matrixI = matrixI.array() - b.matrixI(0,0);
                result.checkComplexity();
            } else {
                result.isComplex = true;
                result.matrixI = constantMatrix(matrixR.rows(),matrixR.cols(), -b.matrixI(0,0));
            }
        } else if (isComplex) {
            result.matrixI = matrixI;
        }
    } else if ((numel() == 1) && (b.numel() != 1)) {
        result.matrixR = matrixR(0,0) - b.matrixR.array();
        if (b.isComplex) {
            if (isComplex) {
                result.matrixI = matrixI(0,0) - b.matrixI.array();
                result.checkComplexity();
            } else {
                result.isComplex = true;
                result.matrixI = -b.matrixI;
            }
        } else if (isComplex) {
            result.matrixI = constantMatrix(matrixR.rows(),matrixR.cols(), matrixI(0,0));
        }
    } else {
        result.matrixR = matrixR - b.matrixR;
        if (b.isComplex) {
            if (isComplex) {
                result.matrixI = matrixI - b.matrixI;
                result.checkComplexity();
            } else {
                result.isComplex = true;
                result.matrixI = -b.matrixI;
            }
        } else if (isComplex) {
            result.matrixI = matrixI;
        }
    }

    return result;
}

/* The substraction c = a-b */
GmpEigenMatrix& GmpEigenMatrix::minus_new(const GmpEigenMatrix& b) const
{
    GmpEigenMatrix& result(*(new GmpEigenMatrix));

    result.isComplex = isComplex;
    if ((numel() != 1) && (b.numel() == 1)) {
        result.matrixR = matrixR.array() - b.matrixR(0,0);
        if (b.isComplex) {
            if (isComplex) {
                result.matrixI = matrixI.array() - b.matrixI(0,0);
                result.checkComplexity();
            } else {
                result.isComplex = true;
                result.matrixI = constantMatrix(matrixR.rows(),matrixR.cols(), -b.matrixI(0,0));
            }
        } else if (isComplex) {
            result.matrixI = matrixI;
        }
    } else if ((numel() == 1) && (b.numel() != 1)) {
        result.matrixR = matrixR(0,0) - b.matrixR.array();
        if (b.isComplex) {
            if (isComplex) {
                result.matrixI = matrixI(0,0) - b.matrixI.array();
                result.checkComplexity();
            } else {
                result.isComplex = true;
                result.matrixI = -b.matrixI;
            }
        } else if (isComplex) {
            result.matrixI = constantMatrix(matrixR.rows(),matrixR.cols(), matrixI(0,0));
        }
    } else {
        result.matrixR = matrixR - b.matrixR;
        if (b.isComplex) {
            if (isComplex) {
                result.matrixI = matrixI - b.matrixI;
                result.checkComplexity();
            } else {
                result.isComplex = true;
                result.matrixI = -b.matrixI;
            }
        } else if (isComplex) {
            result.matrixI = matrixI;
        }
    }

    return result;
}


/* The element-wise multiplication c = a.*b */
GmpEigenMatrix GmpEigenMatrix::times(const GmpEigenMatrix& b) const
{
    GmpEigenMatrix result;

    // We want to support nxm matrix multiplied by 1x1 matrix (which we see as a scalar)
    if ((numel() != 1) && (b.numel() == 1)) {
        if (b.isComplex) {
            if (isComplex) {
                result.matrixR = matrixR*b.matrixR(0,0) - matrixI*b.matrixI(0,0);
                result.matrixI = matrixR*b.matrixI(0,0) + matrixI*b.matrixR(0,0);
                result.checkComplexity();
            } else {
                result.matrixR = matrixR*b.matrixR(0,0);
                result.matrixI = matrixR*b.matrixI(0,0);
                result.isComplex = true;
            }
        } else if (isComplex) {
            result.matrixR = matrixR*b.matrixR(0,0);
            result.matrixI = matrixI*b.matrixR(0,0);
            result.isComplex = true;
        } else {
            result.matrixR = matrixR*b.matrixR(0,0);
            result.isComplex = false;
        }
    } else if ((numel() == 1) && (b.numel() != 1)) {
        if (b.isComplex) {
            if (isComplex) {
                result.matrixR = matrixR(0,0)*b.matrixR - matrixI(0,0)*b.matrixI;
                result.matrixI = matrixR(0,0)*b.matrixI + matrixI(0,0)*b.matrixR;
                result.checkComplexity();
            } else {
                result.matrixR = matrixR(0,0)*b.matrixR;
                result.matrixI = matrixR(0,0)*b.matrixI;
                result.isComplex = true;
            }
        } else if (isComplex) {
            result.matrixR = matrixR(0,0)*b.matrixR;
            result.matrixI = matrixI(0,0)*b.matrixR;
            result.isComplex = true;
        } else {
            result.matrixR = matrixR(0,0)*b.matrixR;
            result.isComplex = false;
        }
    } else {
        if (b.isComplex) {
            if (isComplex) {
                result.matrixR = matrixR.array()*b.matrixR.array() - matrixI.array()*b.matrixI.array();
                result.matrixI = matrixR.array()*b.matrixI.array() + matrixI.array()*b.matrixR.array();
                result.checkComplexity();
            } else {
                result.matrixR = matrixR.array()*b.matrixR.array();
                result.matrixI = matrixR.array()*b.matrixI.array();
                result.checkComplexity();
            }
        } else if (isComplex) {
            result.matrixR = matrixR.array()*b.matrixR.array();
            result.matrixI = matrixI.array()*b.matrixR.array();
            result.checkComplexity();
        } else {
            result.matrixR = matrixR.array()*b.matrixR.array();
            result.isComplex = false;
        }
    }

    return result;
}

/* The element-wise multiplication c = a.*b */
GmpEigenMatrix& GmpEigenMatrix::times_new(const GmpEigenMatrix& b) const
{
    GmpEigenMatrix& result(*(new GmpEigenMatrix));

    // We want to support nxm matrix multiplied by 1x1 matrix (which we see as a scalar)
    if ((numel() != 1) && (b.numel() == 1)) {
        if (b.isComplex) {
            if (isComplex) {
                result.matrixR = matrixR*b.matrixR(0,0) - matrixI*b.matrixI(0,0);
                result.matrixI = matrixR*b.matrixI(0,0) + matrixI*b.matrixR(0,0);
                result.checkComplexity();
            } else {
                result.matrixR = matrixR*b.matrixR(0,0);
                result.matrixI = matrixR*b.matrixI(0,0);
                result.isComplex = true;
            }
        } else if (isComplex) {
            result.matrixR = matrixR*b.matrixR(0,0);
            result.matrixI = matrixI*b.matrixR(0,0);
            result.isComplex = true;
        } else {
            result.matrixR = matrixR*b.matrixR(0,0);
            result.isComplex = false;
        }
    } else if ((numel() == 1) && (b.numel() != 1)) {
        if (b.isComplex) {
            if (isComplex) {
                result.matrixR = matrixR(0,0)*b.matrixR - matrixI(0,0)*b.matrixI;
                result.matrixI = matrixR(0,0)*b.matrixI + matrixI(0,0)*b.matrixR;
                result.checkComplexity();
            } else {
                result.matrixR = matrixR(0,0)*b.matrixR;
                result.matrixI = matrixR(0,0)*b.matrixI;
                result.isComplex = true;
            }
        } else if (isComplex) {
            result.matrixR = matrixR(0,0)*b.matrixR;
            result.matrixI = matrixI(0,0)*b.matrixR;
            result.isComplex = true;
        } else {
            result.matrixR = matrixR(0,0)*b.matrixR;
            result.isComplex = false;
        }
    } else {
        if (b.isComplex) {
            if (isComplex) {
                result.matrixR = matrixR.array()*b.matrixR.array() - matrixI.array()*b.matrixI.array();
                result.matrixI = matrixR.array()*b.matrixI.array() + matrixI.array()*b.matrixR.array();
                result.checkComplexity();
            } else {
                result.matrixR = matrixR.array()*b.matrixR.array();
                result.matrixI = matrixR.array()*b.matrixI.array();
                result.checkComplexity();
            }
        } else if (isComplex) {
            result.matrixR = matrixR.array()*b.matrixR.array();
            result.matrixI = matrixI.array()*b.matrixR.array();
            result.checkComplexity();
        } else {
            result.matrixR = matrixR.array()*b.matrixR.array();
            result.isComplex = false;
        }
    }

    return result;
}


/* Element-wise multiplication with a sparse matrix c = a.*b */
SparseGmpEigenMatrix GmpEigenMatrix::times_fs(const SparseGmpEigenMatrix& b) const
{
    SparseGmpEigenMatrix result;

    // setting the output size
    result.matrixR.resize(matrixR.rows(), matrixR.cols());
    if ((isComplex) || (b.isComplex))
        result.matrixI.resize(matrixR.rows(), matrixR.cols());

    // We will first copy the data to triplets
    vector< Triplet<mpreal> > tripletListR, tripletListI;
    tripletListR.reserve(b.matrixR.nonZeros() + b.matrixI.nonZeros());
    if ((isComplex) || (b.isComplex))
        tripletListI.reserve(b.matrixR.nonZeros() + b.matrixR.nonZeros());

    // Now for each column, we merge the lists of lines with nonzero elements of
    // both matrices. The cases are slightly different depending on the complexity
    // of both matrices
    if (b.isComplex) {
        if (isComplex) {
            for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
                SparseMatrix<mpreal>::InnerIterator itRb(b.matrixR,k);
                SparseMatrix<mpreal>::InnerIterator itIb(b.matrixI,k);
                while ((itRb) || (itIb)) {
                    if ((!itIb) || ((itRb) && (itRb.row() < itIb.row()))) {
                        tripletListR.push_back(Triplet<mpreal>(itRb.row(), k, matrixR(itRb.row(), k)*itRb.value()));
                        tripletListI.push_back(Triplet<mpreal>(itRb.row(), k, matrixI(itRb.row(), k)*itRb.value()));
                        ++itRb;
                    } else if ((itRb) && (itIb) && (itRb.row() == itIb.row())) {
                        tripletListR.push_back(Triplet<mpreal>(itRb.row(), k, matrixR(itRb.row(), k)*itRb.value() - matrixI(itRb.row(), k)*itIb.value()));
                        tripletListI.push_back(Triplet<mpreal>(itRb.row(), k, matrixR(itRb.row(), k)*itIb.value() + matrixI(itRb.row(), k)*itRb.value()));
                        ++itRb;
                        ++itIb;
                    } else {
                        tripletListR.push_back(Triplet<mpreal>(itIb.row(), k, -matrixI(itIb.row(), k)*itIb.value()));
                        tripletListI.push_back(Triplet<mpreal>(itIb.row(), k, matrixR(itIb.row(), k)*itIb.value()));
                        ++itIb;
                    }
                }
            }
        } else {
            for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
                SparseMatrix<mpreal>::InnerIterator itRb(b.matrixR,k);
                SparseMatrix<mpreal>::InnerIterator itIb(b.matrixI,k);
                while ((itRb) || (itIb)) {
                    if ((!itIb) || ((itRb) && (itRb.row() < itIb.row()))) {
                        tripletListR.push_back(Triplet<mpreal>(itRb.row(), k, matrixR(itRb.row(), k)*itRb.value()));
                        ++itRb;
                    } else if ((itRb) && (itIb) && (itRb.row() == itIb.row())) {
                        tripletListR.push_back(Triplet<mpreal>(itRb.row(), k, matrixR(itRb.row(), k)*itRb.value()));
                        tripletListI.push_back(Triplet<mpreal>(itRb.row(), k, matrixR(itRb.row(), k)*itIb.value()));
                        ++itRb;
                        ++itIb;
                    } else {
                        tripletListI.push_back(Triplet<mpreal>(itIb.row(), k, matrixR(itIb.row(), k)*itIb.value()));
                        ++itIb;
                    }
                }
            }
        }
    } else if (isComplex) {
        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            SparseMatrix<mpreal>::InnerIterator itRb(b.matrixR,k);
            while (itRb) {
                tripletListR.push_back(Triplet<mpreal>(itRb.row(), k, matrixR(itRb.row(), k)*itRb.value()));
                tripletListI.push_back(Triplet<mpreal>(itRb.row(), k, matrixI(itRb.row(), k)*itRb.value()));
                ++itRb;
            }
        }
    } else {
        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            SparseMatrix<mpreal>::InnerIterator itRb(b.matrixR,k);
            while (itRb) {
                tripletListR.push_back(Triplet<mpreal>(itRb.row(), k, matrixR(itRb.row(), k)*itRb.value()));
                ++itRb;
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

/* Element-wise multiplication with a sparse matrix c = a.*b */
SparseGmpEigenMatrix& GmpEigenMatrix::times_fs_new(const SparseGmpEigenMatrix& b) const
{
    SparseGmpEigenMatrix& result(*(new SparseGmpEigenMatrix));

    // setting the output size
    result.matrixR.resize(matrixR.rows(), matrixR.cols());
    if ((isComplex) || (b.isComplex))
        result.matrixI.resize(matrixR.rows(), matrixR.cols());

    // We will first copy the data to triplets
    vector< Triplet<mpreal> > tripletListR, tripletListI;
    tripletListR.reserve(b.matrixR.nonZeros() + b.matrixI.nonZeros());
    if ((isComplex) || (b.isComplex))
        tripletListI.reserve(b.matrixR.nonZeros() + b.matrixR.nonZeros());

    // Now for each column, we merge the lists of lines with nonzero elements of
    // both matrices. The cases are slightly different depending on the complexity
    // of both matrices
    if (b.isComplex) {
        if (isComplex) {
            for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
                SparseMatrix<mpreal>::InnerIterator itRb(b.matrixR,k);
                SparseMatrix<mpreal>::InnerIterator itIb(b.matrixI,k);
                while ((itRb) || (itIb)) {
                    if ((!itIb) || ((itRb) && (itRb.row() < itIb.row()))) {
                        tripletListR.push_back(Triplet<mpreal>(itRb.row(), k, matrixR(itRb.row(), k)*itRb.value()));
                        tripletListI.push_back(Triplet<mpreal>(itRb.row(), k, matrixI(itRb.row(), k)*itRb.value()));
                        ++itRb;
                    } else if ((itRb) && (itIb) && (itRb.row() == itIb.row())) {
                        tripletListR.push_back(Triplet<mpreal>(itRb.row(), k, matrixR(itRb.row(), k)*itRb.value() - matrixI(itRb.row(), k)*itIb.value()));
                        tripletListI.push_back(Triplet<mpreal>(itRb.row(), k, matrixR(itRb.row(), k)*itIb.value() + matrixI(itRb.row(), k)*itRb.value()));
                        ++itRb;
                        ++itIb;
                    } else {
                        tripletListR.push_back(Triplet<mpreal>(itIb.row(), k, -matrixI(itIb.row(), k)*itIb.value()));
                        tripletListI.push_back(Triplet<mpreal>(itIb.row(), k, matrixR(itIb.row(), k)*itIb.value()));
                        ++itIb;
                    }
                }
            }
        } else {
            for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
                SparseMatrix<mpreal>::InnerIterator itRb(b.matrixR,k);
                SparseMatrix<mpreal>::InnerIterator itIb(b.matrixI,k);
                while ((itRb) || (itIb)) {
                    if ((!itIb) || ((itRb) && (itRb.row() < itIb.row()))) {
                        tripletListR.push_back(Triplet<mpreal>(itRb.row(), k, matrixR(itRb.row(), k)*itRb.value()));
                        ++itRb;
                    } else if ((itRb) && (itIb) && (itRb.row() == itIb.row())) {
                        tripletListR.push_back(Triplet<mpreal>(itRb.row(), k, matrixR(itRb.row(), k)*itRb.value()));
                        tripletListI.push_back(Triplet<mpreal>(itRb.row(), k, matrixR(itRb.row(), k)*itIb.value()));
                        ++itRb;
                        ++itIb;
                    } else {
                        tripletListI.push_back(Triplet<mpreal>(itIb.row(), k, matrixR(itIb.row(), k)*itIb.value()));
                        ++itIb;
                    }
                }
            }
        }
    } else if (isComplex) {
        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            SparseMatrix<mpreal>::InnerIterator itRb(b.matrixR,k);
            while (itRb) {
                tripletListR.push_back(Triplet<mpreal>(itRb.row(), k, matrixR(itRb.row(), k)*itRb.value()));
                tripletListI.push_back(Triplet<mpreal>(itRb.row(), k, matrixI(itRb.row(), k)*itRb.value()));
                ++itRb;
            }
        }
    } else {
        for (IndexType k = 0; k < matrixR.outerSize(); ++k) {
            SparseMatrix<mpreal>::InnerIterator itRb(b.matrixR,k);
            while (itRb) {
                tripletListR.push_back(Triplet<mpreal>(itRb.row(), k, matrixR(itRb.row(), k)*itRb.value()));
                ++itRb;
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
GmpEigenMatrix GmpEigenMatrix::rdivide(const GmpEigenMatrix& b) const
{
    GmpEigenMatrix result;

    // We want to support nxm matrix divided by 1x1 matrix (which we see as a scalar)
    if ((numel() != 1) && (b.numel() == 1)) {
        if (b.isComplex) {
            if (isComplex) {
                result.matrixR = (matrixR*b.matrixR(0,0) + matrixI*b.matrixI(0,0))/(mpfr::pow(b.matrixR(0,0),2) + mpfr::pow(b.matrixI(0,0),2));
                result.matrixI = (-matrixR*b.matrixI(0,0) + matrixI*b.matrixR(0,0))/(mpfr::pow(b.matrixR(0,0),2) + mpfr::pow(b.matrixI(0,0),2));
                result.checkComplexity();
            } else {
                result.matrixR = (matrixR*b.matrixR(0,0))/(mpfr::pow(b.matrixR(0,0),2) + mpfr::pow(b.matrixI(0,0),2));
                result.matrixI = (-matrixR*b.matrixI(0,0))/(mpfr::pow(b.matrixR(0,0),2) + mpfr::pow(b.matrixI(0,0),2));
                result.isComplex = true;
            }
        } else if (isComplex) {
            result.matrixR = matrixR/b.matrixR(0,0);
            result.matrixI = matrixI/b.matrixR(0,0);
            result.isComplex = true;
        } else {
            result.matrixR = matrixR/b.matrixR(0,0);
            result.isComplex = false;
        }
    } else if ((numel() == 1) && (b.numel() != 1)) {
        if (b.isComplex) {
            if (isComplex) {
                result.matrixR = (matrixR(0,0)*b.matrixR + matrixI(0,0)*b.matrixI).array()/(b.matrixR.array().pow(2) + b.matrixI.array().pow(2));
                result.matrixI = (-matrixR(0,0)*b.matrixI + matrixI(0,0)*b.matrixR).array()/(b.matrixR.array().pow(2) + b.matrixI.array().pow(2));
                result.checkComplexity();
            } else {
                result.matrixR = (matrixR(0,0)*b.matrixR).array()/(b.matrixR.array().pow(2) + b.matrixI.array().pow(2));
                result.matrixI = (-matrixR(0,0)*b.matrixI).array()/(b.matrixR.array().pow(2) + b.matrixI.array().pow(2));
                result.isComplex = true;
            }
        } else if (isComplex) {
            result.matrixR = constantMatrix(b.matrixR.rows(),b.matrixR.cols(),matrixR(0,0)).array()/b.matrixR.array();
            result.matrixI = constantMatrix(b.matrixR.rows(),b.matrixR.cols(),matrixI(0,0)).array()/b.matrixR.array();
            result.isComplex = true;
        } else {
            result.matrixR = constantMatrix(b.matrixR.rows(),b.matrixR.cols(),matrixR(0,0)).array()/b.matrixR.array();
            result.isComplex = false;
        }
    } else {
        if (b.isComplex) {
            if (isComplex) {
                result.matrixR = (matrixR.array()*b.matrixR.array() + matrixI.array()*b.matrixI.array())/(b.matrixR.array().pow(2) + b.matrixI.array().pow(2));
                result.matrixI = (-matrixR.array()*b.matrixI.array() + matrixI.array()*b.matrixR.array())/(b.matrixR.array().pow(2) + b.matrixI.array().pow(2));
                result.checkComplexity();
            } else {
                result.matrixR = (matrixR.array()*b.matrixR.array())/(b.matrixR.array().pow(2) + b.matrixI.array().pow(2));
                result.matrixI = (-matrixR.array()*b.matrixI.array())/(b.matrixR.array().pow(2) + b.matrixI.array().pow(2));
                result.checkComplexity();
            }
        } else if (isComplex) {
            result.matrixR = matrixR.array()/b.matrixR.array();
            result.matrixI = matrixI.array()/b.matrixR.array();
            result.checkComplexity();
        } else {
            result.matrixR = matrixR.array()/b.matrixR.array();
            result.isComplex = false;
        }
    }

    return result;
}

/* The element-wise division on the right c = a./b */
GmpEigenMatrix& GmpEigenMatrix::rdivide_new(const GmpEigenMatrix& b) const
{
    GmpEigenMatrix& result(*(new GmpEigenMatrix));

    // We want to support nxm matrix divided by 1x1 matrix (which we see as a scalar)
    if ((numel() != 1) && (b.numel() == 1)) {
        if (b.isComplex) {
            if (isComplex) {
                result.matrixR = (matrixR*b.matrixR(0,0) + matrixI*b.matrixI(0,0))/(mpfr::pow(b.matrixR(0,0),2) + mpfr::pow(b.matrixI(0,0),2));
                result.matrixI = (-matrixR*b.matrixI(0,0) + matrixI*b.matrixR(0,0))/(mpfr::pow(b.matrixR(0,0),2) + mpfr::pow(b.matrixI(0,0),2));
                result.checkComplexity();
            } else {
                result.matrixR = (matrixR*b.matrixR(0,0))/(mpfr::pow(b.matrixR(0,0),2) + mpfr::pow(b.matrixI(0,0),2));
                result.matrixI = (-matrixR*b.matrixI(0,0))/(mpfr::pow(b.matrixR(0,0),2) + mpfr::pow(b.matrixI(0,0),2));
                result.isComplex = true;
            }
        } else if (isComplex) {
            result.matrixR = matrixR/b.matrixR(0,0);
            result.matrixI = matrixI/b.matrixR(0,0);
            result.isComplex = true;
        } else {
            result.matrixR = matrixR/b.matrixR(0,0);
            result.isComplex = false;
        }
    } else if ((numel() == 1) && (b.numel() != 1)) {
        if (b.isComplex) {
            if (isComplex) {
                result.matrixR = (matrixR(0,0)*b.matrixR + matrixI(0,0)*b.matrixI).array()/(b.matrixR.array().pow(2) + b.matrixI.array().pow(2));
                result.matrixI = (-matrixR(0,0)*b.matrixI + matrixI(0,0)*b.matrixR).array()/(b.matrixR.array().pow(2) + b.matrixI.array().pow(2));
                result.checkComplexity();
            } else {
                result.matrixR = (matrixR(0,0)*b.matrixR).array()/(b.matrixR.array().pow(2) + b.matrixI.array().pow(2));
                result.matrixI = (-matrixR(0,0)*b.matrixI).array()/(b.matrixR.array().pow(2) + b.matrixI.array().pow(2));
                result.isComplex = true;
            }
        } else if (isComplex) {
            result.matrixR = constantMatrix(b.matrixR.rows(),b.matrixR.cols(),matrixR(0,0)).array()/b.matrixR.array();
            result.matrixI = constantMatrix(b.matrixR.rows(),b.matrixR.cols(),matrixI(0,0)).array()/b.matrixR.array();
            result.isComplex = true;
        } else {
            result.matrixR = constantMatrix(b.matrixR.rows(),b.matrixR.cols(),matrixR(0,0)).array()/b.matrixR.array();
            result.isComplex = false;
        }
    } else {
        if (b.isComplex) {
            if (isComplex) {
                result.matrixR = (matrixR.array()*b.matrixR.array() + matrixI.array()*b.matrixI.array())/(b.matrixR.array().pow(2) + b.matrixI.array().pow(2));
                result.matrixI = (-matrixR.array()*b.matrixI.array() + matrixI.array()*b.matrixR.array())/(b.matrixR.array().pow(2) + b.matrixI.array().pow(2));
                result.checkComplexity();
            } else {
                result.matrixR = (matrixR.array()*b.matrixR.array())/(b.matrixR.array().pow(2) + b.matrixI.array().pow(2));
                result.matrixI = (-matrixR.array()*b.matrixI.array())/(b.matrixR.array().pow(2) + b.matrixI.array().pow(2));
                result.checkComplexity();
            }
        } else if (isComplex) {
            result.matrixR = matrixR.array()/b.matrixR.array();
            result.matrixI = matrixI.array()/b.matrixR.array();
            result.checkComplexity();
        } else {
            result.matrixR = matrixR.array()/b.matrixR.array();
            result.isComplex = false;
        }
    }

    return result;
}


/* The absolute value (or complex magnitude) b = abs(a) */
// Note that this function is needed to define the operation power, so we
// don't call that function here.
GmpEigenMatrix GmpEigenMatrix::abs() const
{
    GmpEigenMatrix result;

    result.isComplex = false;
    if (isComplex)
        result.matrixR = (matrixR.array().pow(2) + matrixI.array().pow(2)).sqrt();
    else
        result.matrixR = matrixR.array().abs();

    return result;
}

/* The absolute value (or complex magnitude) b = abs(a) */
GmpEigenMatrix& GmpEigenMatrix::abs_new() const
{
    GmpEigenMatrix& result(*(new GmpEigenMatrix));

    result.isComplex = false;
    if (isComplex)
        result.matrixR = (matrixR.array().pow(2) + matrixI.array().pow(2)).sqrt();
    else
        result.matrixR = matrixR.array().abs();

    return result;
}


/* The phase angle (i.e. complex argument) b = angle(a) */
GmpEigenMatrix GmpEigenMatrix::angle() const
{
    GmpEigenMatrix result;

    // Eigen doesn't implement atan2, so we iterate over all elements by hand
    // NOTE : Maybe we could use Eigen::CwiseUnaryOp for that as well
    // or maybe mat1.unaryExpr(std::ptr_fun(foo)) , but std::ptr_fun is
    // depreciated since C++11, so we keep with just loops for now
    result.isComplex = false;
    result.matrixR.resize(matrixR.rows(),matrixR.cols());
    if (isComplex) {
        for (IndexType j(0); j < matrixR.cols(); ++j)
            for (IndexType i(0); i < matrixR.rows(); ++i)
                result.matrixR(i,j) = atan2(matrixI(i,j), matrixR(i,j));
    } else {
        for (IndexType j(0); j < matrixR.cols(); ++j)
            for (IndexType i(0); i < matrixR.rows(); ++i)
                result.matrixR(i,j) = atan2(mpreal("0"), matrixR(i,j));
    }

    return result;
}

/* The phase angle (i.e. complex argument) b = angle(a) */
GmpEigenMatrix& GmpEigenMatrix::angle_new() const
{
    GmpEigenMatrix& result(*(new GmpEigenMatrix));

    // Eigen doesn't implement atan2, so we iterate over all elements by hand
    result.isComplex = false;
    result.matrixR.resize(matrixR.rows(),matrixR.cols());
    if (isComplex) {
        for (IndexType j(0); j < matrixR.cols(); ++j)
            for (IndexType i(0); i < matrixR.rows(); ++i)
                result.matrixR(i,j) = atan2(matrixI(i,j), matrixR(i,j));
    } else {
        for (IndexType j(0); j < matrixR.cols(); ++j)
            for (IndexType i(0); i < matrixR.rows(); ++i)
                result.matrixR(i,j) = atan2(mpreal("0"), matrixR(i,j));
    }

    return result;
}



/* The element-wise exponential b = exp(a) */
GmpEigenMatrix GmpEigenMatrix::exp() const
{
    GmpEigenMatrix result;

    if (isComplex) {
        result.matrixR = matrixR.array().exp()*matrixI.array().cos();
        result.matrixI = matrixR.array().exp()*matrixI.array().sin();
        result.checkComplexity();
    } else {
        result.isComplex = false;
        result.matrixR = matrixR.array().exp();
    }

    return result;
}

/* The element-wise exponential b = exp(a) */
GmpEigenMatrix& GmpEigenMatrix::exp_new() const
{
    GmpEigenMatrix& result(*(new GmpEigenMatrix));

    if (isComplex) {
        result.matrixR = matrixR.array().exp()*matrixI.array().cos();
        result.matrixI = matrixR.array().exp()*matrixI.array().sin();
        result.checkComplexity();
    } else {
        result.isComplex = false;
        result.matrixR = matrixR.array().exp();
    }

    return result;
}

/* The element-wise natural logarithm b = log(a) */
GmpEigenMatrix GmpEigenMatrix::log() const
{
    GmpEigenMatrix result;

    GmpEigenMatrix tmp(GmpEigenMatrix::abs());
    result.matrixR = tmp.matrixR.array().log();

    tmp = GmpEigenMatrix::angle();
    result.matrixI = tmp.matrixR;

    result.checkComplexity();

    return result;
}

/* The element-wise natural logarithm b = log(a) */
GmpEigenMatrix& GmpEigenMatrix::log_new() const
{
    GmpEigenMatrix& result(*(new GmpEigenMatrix));

    GmpEigenMatrix tmp(GmpEigenMatrix::abs());
    result.matrixR = tmp.matrixR.array().log();

    tmp = GmpEigenMatrix::angle();
    result.matrixI = tmp.matrixR;

    result.checkComplexity();

    return result;
}


/* The element-wise power c = a.^b */
// knowing that a^b = exp(b*log(a))
GmpEigenMatrix GmpEigenMatrix::power(const GmpEigenMatrix& b) const
{
    GmpEigenMatrix result;

    if ((isComplex) || (b.isComplex) || (!b.isInt())) {
        result = b.times(log()).exp();

        // Since this formula give 0^x = NaN, we need to enforce by hand that 0^0 = 1, 0^(-2) = Inf, 0^1i = NaN, and 0^3 = 0, but also 0^Inf=0, 0^(-Inf)=0 ...
        if ((b.numel() == 1) && (!b.isComplex)) {
            for (IndexType j = 0; j < matrixR.cols(); ++j) {
                for (IndexType i = 0; i < matrixR.rows(); ++i) {
                    if ((matrixR(i,j) == 0) && ((!isComplex) || ((isComplex) && (matrixI(i,j) == 0)))) {
                        if (b.matrixR(0,0) == 0) {
                            result.matrixR(i,j) = 1;
                            if (result.isComplex)
                                result.matrixI(i,j) = 0;
                        } else if (b.matrixR(0,0) < 0) {
                            result.matrixR(i,j) = mpreal("Inf");
                            if (result.isComplex)
                                result.matrixI(i,j) = 0;
                        } else if (b.matrixR(0,0) > 0) {
                            result.matrixR(i,j) = 0;
                            if (result.isComplex)
                                result.matrixI(i,j) = 0;
                        }
                    }
                }
            }
        } else if (b.numel() != 1) {
            for (IndexType j = 0; j < matrixR.cols(); ++j) {
                for (IndexType i = 0; i < matrixR.rows(); ++i) {
                    if (((!b.isComplex) || ((b.isComplex) && (b.matrixI(i,j) == 0))) && (matrixR(i,j) == 0) && ((!isComplex) || ((isComplex) && (matrixI(i,j) == 0)))) {
                        if (b.matrixR(i,j) == 0) {
                            result.matrixR(i,j) = 1;
                            if (result.isComplex)
                                result.matrixI(i,j) = 0;
                        } else if (b.matrixR(i,j) < 0) {
                            result.matrixR(i,j) = mpreal("Inf");
                            if (result.isComplex)
                                result.matrixI(i,j) = 0;
                        } else if (b.matrixR(i,j) > 0) {
                            result.matrixR(i,j) = 0;
                            if (result.isComplex)
                                result.matrixI(i,j) = 0;
                        }
                    }
                }
            }
        }
    } else {
        result.isComplex = false;
        if ((numel() != 1) && (b.numel() == 1)) {
            result.matrixR = matrixR.array().pow(b.matrixR(0,0));
        } else if ((numel() == 1) && (b.numel() != 1)) {
            result.matrixR.resize(b.matrixR.rows(),b.matrixR.cols());
            for (IndexType j(0); j < result.matrixR.cols(); ++j)
                for (IndexType i(0); i < result.matrixR.rows(); ++i)
                    result.matrixR(i,j) = mpfr::pow(matrixR(0,0), b.matrixR(i,j));
        } else {
            result.matrixR.resize(matrixR.rows(),matrixR.cols());
            for (IndexType j(0); j < result.matrixR.cols(); ++j)
                for (IndexType i(0); i < result.matrixR.rows(); ++i)
                    result.matrixR(i,j) = mpfr::pow(matrixR(i,j), b.matrixR(i,j));
        }
    }

    return result;
}

/* The element-wise power c = a.^b */
GmpEigenMatrix& GmpEigenMatrix::power_new(const GmpEigenMatrix& b) const
{
    GmpEigenMatrix& result(*(new GmpEigenMatrix));

    if ((isComplex) || (b.isComplex) || (!b.isInt())) {
        result = b.times(log()).exp();

        // Since this formula give 0^x = NaN, we need to enforce by hand that 0^0 = 1, 0^(-2) = Inf, 0^1i = NaN, and 0^3 = 0
        if ((b.numel() == 1) && (!b.isComplex)) {
            for (IndexType j = 0; j < matrixR.cols(); ++j) {
                for (IndexType i = 0; i < matrixR.rows(); ++i) {
                    if ((matrixR(i,j) == 0) && ((!isComplex) || ((isComplex) && (matrixI(i,j) == 0)))) {
                        if (b.matrixR(0,0) == 0) {
                            result.matrixR(i,j) = 1;
                            if (result.isComplex)
                                result.matrixI(i,j) = 0;
                        } else if (b.matrixR(0,0) < 0) {
                            result.matrixR(i,j) = mpreal("Inf");
                            if (result.isComplex)
                                result.matrixI(i,j) = 0;
                        } else if (b.matrixR(0,0) > 0) {
                            result.matrixR(i,j) = 0;
                            if (result.isComplex)
                                result.matrixI(i,j) = 0;
                        }
                    }
                }
            }
        } else if (b.numel() != 1) {
            for (IndexType j = 0; j < matrixR.cols(); ++j) {
                for (IndexType i = 0; i < matrixR.rows(); ++i) {
                    if (((!b.isComplex) || ((b.isComplex) && (b.matrixI(i,j) == 0))) && (matrixR(i,j) == 0) && ((!isComplex) || ((isComplex) && (matrixI(i,j) == 0)))) {
                        if (b.matrixR(i,j) == 0) {
                            result.matrixR(i,j) = 1;
                            if (result.isComplex)
                                result.matrixI(i,j) = 0;
                        } else if (b.matrixR(i,j) < 0) {
                            result.matrixR(i,j) = mpreal("Inf");
                            if (result.isComplex)
                                result.matrixI(i,j) = 0;
                        } else if (b.matrixR(i,j) > 0) {
                            result.matrixR(i,j) = 0;
                            if (result.isComplex)
                                result.matrixI(i,j) = 0;
                        }
                    }
                }
            }
        }
    } else {
        result.isComplex = false;
        if ((numel() != 1) && (b.numel() == 1)) {
            result.matrixR = matrixR.array().pow(b.matrixR(0,0));
        } else if ((numel() == 1) && (b.numel() != 1)) {
            result.matrixR.resize(b.matrixR.rows(),b.matrixR.cols());
            for (IndexType j(0); j < result.matrixR.cols(); ++j)
                for (IndexType i(0); i < result.matrixR.rows(); ++i)
                    result.matrixR(i,j) = mpfr::pow(matrixR(0,0), b.matrixR(i,j));
        } else {
            result.matrixR.resize(matrixR.rows(),matrixR.cols());
            for (IndexType j(0); j < result.matrixR.cols(); ++j)
                for (IndexType i(0); i < result.matrixR.rows(); ++i)
                    result.matrixR(i,j) = mpfr::pow(matrixR(i,j), b.matrixR(i,j));
        }
    }

    return result;
}


/* square root b = sqrt(a) */
GmpEigenMatrix GmpEigenMatrix::sqrt() const
{
    GmpEigenMatrix result;

    mpreal globalMin(0);
    if (!isComplex) {
        // In this case we need to compute the smallest element in the matrix
        GmpEigenMatrix globalMinMatrix(rowMin().colMin());
        globalMin = globalMinMatrix.matrixR(0,0);
    }
    if ((isComplex) || (globalMin < 0)) {
        result.matrixR.resize(matrixR.rows(),matrixR.cols());
        result.matrixI.resize(matrixR.rows(),matrixR.cols());
        GmpEigenMatrix norm(abs());
        for (IndexType j(0); j < result.matrixR.cols(); ++j)
            for (IndexType i(0); i < result.matrixR.rows(); ++i)
            {
                // Re(sqrt(-Inf+c*1i)) = 0
                if ((matrixR(i,j) == mpreal("-Inf")) && (isfinite(matrixI(i,j))))
                    result.matrixR(i,j) = 0;
                else
                    result.matrixR(i,j) = mpfr::sqrt((norm.matrixR(i,j) + matrixR(i,j))/2);
                // Im(sqrt(NaN)) = 0 (by convention)
                // and
                // Im(sqrt(Inf+c*1i)) = 0
                if ((mpfr::isnan(matrixR(i,j))) && (matrixI(i,j) == 0))
                    result.matrixI(i,j) = 0;
                else if ((matrixR(i,j) == mpreal("Inf")) && (isfinite(matrixI(i,j))))
                    result.matrixI(i,j) = 0;
                else
                    result.matrixI(i,j) = ((isComplex) && (matrixI(i,j) < 0) ? -mpfr::sqrt((norm.matrixR(i,j) - matrixR(i,j))/2) : mpfr::sqrt((norm.matrixR(i,j) - matrixR(i,j))/2));
            }
    } else
        result.matrixR = matrixR.array().sqrt();

    result.checkComplexity();

    return result;
}

/* square root b = sqrt(a) */
GmpEigenMatrix& GmpEigenMatrix::sqrt_new() const
{
    GmpEigenMatrix& result(*(new GmpEigenMatrix));

    mpreal globalMin(0);
    if (!isComplex) {
        // In this case we need to compute the smallest element in the matrix
        GmpEigenMatrix globalMinMatrix(rowMin().colMin());
        globalMin = globalMinMatrix.matrixR(0,0);
    }
    if ((isComplex) || (globalMin < 0)) {
        result.matrixR.resize(matrixR.rows(),matrixR.cols());
        result.matrixI.resize(matrixR.rows(),matrixR.cols());
        GmpEigenMatrix norm(abs());
        for (IndexType j(0); j < result.matrixR.cols(); ++j)
            for (IndexType i(0); i < result.matrixR.rows(); ++i)
            {
                // Re(sqrt(-Inf+c*1i)) = 0
                if ((matrixR(i,j) == mpreal("-Inf")) && (isfinite(matrixI(i,j))))
                    result.matrixR(i,j) = 0;
                else
                    result.matrixR(i,j) = mpfr::sqrt((norm.matrixR(i,j) + matrixR(i,j))/2);
                // Im(sqrt(NaN)) = 0 (by convention)
                // and
                // Im(sqrt(Inf+c*1i)) = 0
                if ((mpfr::isnan(matrixR(i,j))) && (matrixI(i,j) == 0))
                    result.matrixI(i,j) = 0;
                else if ((matrixR(i,j) == mpreal("Inf")) && (isfinite(matrixI(i,j))))
                    result.matrixI(i,j) = 0;
                else
                    result.matrixI(i,j) = ((isComplex) && (matrixI(i,j) < 0) ? -mpfr::sqrt((norm.matrixR(i,j) - matrixR(i,j))/2) : mpfr::sqrt((norm.matrixR(i,j) - matrixR(i,j))/2));
            }
    } else
        result.matrixR = matrixR.array().sqrt();

    result.checkComplexity();

    return result;
}


/* Matrix multiplication a *= b */
GmpEigenMatrix& GmpEigenMatrix::operator*=(const GmpEigenMatrix& b)
{
    GmpEigenMatrix a((*this));

    (*this) = a*b;

    return *this;
}

/* Matrix multiplication c = a*b */
// By default, Eigen assumes aliasing for matrix multiplication, so we need to
// tell him whenever the matrix on the lhs is not on the rhs
GmpEigenMatrix GmpEigenMatrix::operator*(const GmpEigenMatrix& b) const
{
    GmpEigenMatrix result;

    // This function could be called to multiply two GmpEigenMatrix objects
    // within this library. It should thus support multiplication by a scalar.
    if ((numel() == 1) || (b.numel() == 1)) {
        result = (*this).times(b);
        return result;
    }

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

/* Matrix multiplication c = a*b */
GmpEigenMatrix& GmpEigenMatrix::mtimes_new(const GmpEigenMatrix& b) const
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


/* multiplication with a sparse matrix c = a*b */
GmpEigenMatrix GmpEigenMatrix::mtimes_fs(const SparseGmpEigenMatrix& b) const
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

/* multiplication with a sparse matrix c = a*b */
GmpEigenMatrix& GmpEigenMatrix::mtimes_fs_new(const SparseGmpEigenMatrix& b) const
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
GmpEigenMatrix GmpEigenMatrix::kron(const GmpEigenMatrix& b) const
{
    GmpEigenMatrix result;

    if (b.isComplex) {
        if (isComplex) {
            result.matrixR = kroneckerProduct(matrixR, b.matrixR) - kroneckerProduct(matrixI, b.matrixI);
            result.matrixI = kroneckerProduct(matrixR, b.matrixI) + kroneckerProduct(matrixI, b.matrixR);
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
GmpEigenMatrix& GmpEigenMatrix::kron_new(const GmpEigenMatrix& b) const
{
    GmpEigenMatrix& result(*(new GmpEigenMatrix));

    if (b.isComplex) {
        if (isComplex) {
            result.matrixR = kroneckerProduct(matrixR, b.matrixR) - kroneckerProduct(matrixI, b.matrixI);
            result.matrixI = kroneckerProduct(matrixR, b.matrixI) + kroneckerProduct(matrixI, b.matrixR);
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

/* Tensor product with a sparse matrix c = kron(a,b) */
SparseGmpEigenMatrix GmpEigenMatrix::kron_fs(const SparseGmpEigenMatrix& b) const
{
    SparseGmpEigenMatrix result;

    if (b.isComplex) {
        if (isComplex) {
            result.matrixR = kroneckerProduct(matrixR, b.matrixR);
            //result.matrixR -= kroneckerProduct(matrixI, b.matrixI); // This is not supported anymore on Eigen 3.3...
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

/* Tensor product with a sparse matrix c = kron(a,b) */
SparseGmpEigenMatrix& GmpEigenMatrix::kron_fs_new(const SparseGmpEigenMatrix& b) const
{
    SparseGmpEigenMatrix& result(*(new SparseGmpEigenMatrix));

    if (b.isComplex) {
        if (isComplex) {
            result.matrixR = kroneckerProduct(matrixR, b.matrixR);
            //result.matrixR -= kroneckerProduct(matrixI, b.matrixI); // This is not supported anymore on Eigen 3.3...
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
GmpEigenMatrix GmpEigenMatrix::cbrt() const
{
    GmpEigenMatrix result;

    result.matrixR.resize(matrixR.rows(),matrixR.cols());
    result.isComplex = false;
    for (IndexType j(0); j < result.matrixR.cols(); ++j)
        for (IndexType i(0); i < result.matrixR.rows(); ++i)
            result.matrixR(i,j) = mpfr::cbrt(matrixR(i,j));

    return result;
}

/* cubic root b = cqrt(a) -- only works if the real part is positive and the imaginary part zero (!) */
GmpEigenMatrix& GmpEigenMatrix::cbrt_new() const
{
    GmpEigenMatrix& result(*(new GmpEigenMatrix));

    result.matrixR.resize(matrixR.rows(),matrixR.cols());
    result.isComplex = false;
    for (IndexType j(0); j < result.matrixR.cols(); ++j)
        for (IndexType i(0); i < result.matrixR.rows(); ++i)
            result.matrixR(i,j) = mpfr::cbrt(matrixR(i,j));

    return result;
}

/* sine function b = sin(a) */
GmpEigenMatrix GmpEigenMatrix::sin() const
{
    GmpEigenMatrix result;

    result.matrixR.resize(matrixR.rows(),matrixR.cols());
    result.isComplex = isComplex;
    if (isComplex) {
        result.matrixI.resize(matrixR.rows(),matrixR.cols());
//        result = -constI()*GmpEigenMatrix(1/mpreal(2))*((constI()*(*this)).exp() - (-constI()*(*this)).exp());
        for (IndexType j(0); j < result.matrixR.cols(); ++j)
            for (IndexType i(0); i < result.matrixR.rows(); ++i) {
                result.matrixR(i,j) = mpfr::sin(matrixR(i,j))*mpfr::cosh(matrixI(i,j));
                result.matrixI(i,j) = mpfr::cos(matrixR(i,j))*mpfr::sinh(matrixI(i,j));
            }
        result.checkComplexity();
    } else {
        for (IndexType j(0); j < result.matrixR.cols(); ++j)
            for (IndexType i(0); i < result.matrixR.rows(); ++i)
                result.matrixR(i,j) = mpfr::sin(matrixR(i,j));
    }

    return result;
}

/* sine function b = sin(a) */
GmpEigenMatrix& GmpEigenMatrix::sin_new() const
{
    GmpEigenMatrix& result(*(new GmpEigenMatrix));

    result.matrixR.resize(matrixR.rows(),matrixR.cols());
    result.isComplex = isComplex;
    if (isComplex) {
        result.matrixI.resize(matrixR.rows(),matrixR.cols());
//        result = -constI()*GmpEigenMatrix(1/mpreal(2))*((constI()*(*this)).exp() - (-constI()*(*this)).exp());
        for (IndexType j(0); j < result.matrixR.cols(); ++j)
            for (IndexType i(0); i < result.matrixR.rows(); ++i) {
                result.matrixR(i,j) = mpfr::sin(matrixR(i,j))*mpfr::cosh(matrixI(i,j));
                result.matrixI(i,j) = mpfr::cos(matrixR(i,j))*mpfr::sinh(matrixI(i,j));
            }
        result.checkComplexity();
    } else {
        for (IndexType j(0); j < result.matrixR.cols(); ++j)
            for (IndexType i(0); i < result.matrixR.rows(); ++i)
                result.matrixR(i,j) = mpfr::sin(matrixR(i,j));
    }

    return result;
}

/* cosine function b = cos(a) */
GmpEigenMatrix GmpEigenMatrix::cos() const
{
    GmpEigenMatrix result;

    result.matrixR.resize(matrixR.rows(),matrixR.cols());
    result.isComplex = isComplex;
    if (isComplex) {
        result.matrixI.resize(matrixR.rows(),matrixR.cols());
        for (IndexType j(0); j < result.matrixR.cols(); ++j)
            for (IndexType i(0); i < result.matrixR.rows(); ++i) {
                result.matrixR(i,j) = mpfr::cos(matrixR(i,j))*mpfr::cosh(matrixI(i,j));
                result.matrixI(i,j) = -mpfr::sin(matrixR(i,j))*mpfr::sinh(matrixI(i,j));
            }
        result.checkComplexity();
    } else {
        for (IndexType j(0); j < result.matrixR.cols(); ++j)
            for (IndexType i(0); i < result.matrixR.rows(); ++i)
                result.matrixR(i,j) = mpfr::cos(matrixR(i,j));
    }

    return result;
}

/* cosine function b = cos(a) */
GmpEigenMatrix& GmpEigenMatrix::cos_new() const
{
    GmpEigenMatrix& result(*(new GmpEigenMatrix));

    result.matrixR.resize(matrixR.rows(),matrixR.cols());
    result.isComplex = isComplex;
    if (isComplex) {
        result.matrixI.resize(matrixR.rows(),matrixR.cols());
        for (IndexType j(0); j < result.matrixR.cols(); ++j)
            for (IndexType i(0); i < result.matrixR.rows(); ++i) {
                result.matrixR(i,j) = mpfr::cos(matrixR(i,j))*mpfr::cosh(matrixI(i,j));
                result.matrixI(i,j) = -mpfr::sin(matrixR(i,j))*mpfr::sinh(matrixI(i,j));
            }
        result.checkComplexity();
    } else {
        for (IndexType j(0); j < result.matrixR.cols(); ++j)
            for (IndexType i(0); i < result.matrixR.rows(); ++i)
                result.matrixR(i,j) = mpfr::cos(matrixR(i,j));
    }

    return result;
}

/* tangent function b = tan(a) */
GmpEigenMatrix GmpEigenMatrix::tan() const
{
    GmpEigenMatrix result;

    result.matrixR.resize(matrixR.rows(),matrixR.cols());
    result.isComplex = isComplex;
    if (isComplex) {
        result = ((*this).sin()).rdivide((*this).cos());
        result.checkComplexity();
    } else {
        for (IndexType j(0); j < result.matrixR.cols(); ++j)
            for (IndexType i(0); i < result.matrixR.rows(); ++i)
                result.matrixR(i,j) = mpfr::tan(matrixR(i,j));
    }

    return result;
}

/* tangent function b = tan(a) */
GmpEigenMatrix& GmpEigenMatrix::tan_new() const
{
    GmpEigenMatrix& result(*(new GmpEigenMatrix));

    result.matrixR.resize(matrixR.rows(),matrixR.cols());
    result.isComplex = isComplex;
    if (isComplex) {
        result = ((*this).sin()).rdivide((*this).cos());
        result.checkComplexity();
    } else {
        for (IndexType j(0); j < result.matrixR.cols(); ++j)
            for (IndexType i(0); i < result.matrixR.rows(); ++i)
                result.matrixR(i,j) = mpfr::tan(matrixR(i,j));
    }

    return result;
}

/* secant function b = sec(a) */
GmpEigenMatrix GmpEigenMatrix::sec() const
{
    GmpEigenMatrix result;

    result.isComplex = isComplex;
    if (isComplex) {
        result.matrixR.resize(matrixR.rows(),matrixR.cols());
        result.matrixI.resize(matrixR.rows(),matrixR.cols());
        for (IndexType j(0); j < result.matrixR.cols(); ++j)
            for (IndexType i(0); i < result.matrixR.rows(); ++i) {
                result.matrixR(i,j) = 2*(mpfr::cos(matrixR(i,j))*mpfr::cosh(matrixI(i,j)))/(mpfr::cos(2*matrixR(i,j)) + mpfr::cosh(2*matrixI(i,j)));
                result.matrixI(i,j) = 2*(mpfr::sin(matrixR(i,j))*mpfr::sinh(matrixI(i,j)))/(mpfr::cos(2*matrixR(i,j)) + mpfr::cosh(2*matrixI(i,j)));
            }
        result.checkComplexity();
    } else {
        result.matrixR.resize(matrixR.rows(),matrixR.cols());
        for (IndexType j(0); j < result.matrixR.cols(); ++j)
            for (IndexType i(0); i < result.matrixR.rows(); ++i)
                result.matrixR(i,j) = mpfr::sec(matrixR(i,j));
    }

    return result;
}

/* secant function b = sec(a) */
GmpEigenMatrix& GmpEigenMatrix::sec_new() const
{
    GmpEigenMatrix& result(*(new GmpEigenMatrix));

    result.isComplex = isComplex;
    if (isComplex) {
        result.matrixR.resize(matrixR.rows(),matrixR.cols());
        result.matrixI.resize(matrixR.rows(),matrixR.cols());
        for (IndexType j(0); j < result.matrixR.cols(); ++j)
            for (IndexType i(0); i < result.matrixR.rows(); ++i) {
                result.matrixR(i,j) = 2*(mpfr::cos(matrixR(i,j))*mpfr::cosh(matrixI(i,j)))/(mpfr::cos(2*matrixR(i,j)) + mpfr::cosh(2*matrixI(i,j)));
                result.matrixI(i,j) = 2*(mpfr::sin(matrixR(i,j))*mpfr::sinh(matrixI(i,j)))/(mpfr::cos(2*matrixR(i,j)) + mpfr::cosh(2*matrixI(i,j)));
            }
        result.checkComplexity();
    } else {
        result.matrixR.resize(matrixR.rows(),matrixR.cols());
        for (IndexType j(0); j < result.matrixR.cols(); ++j)
            for (IndexType i(0); i < result.matrixR.rows(); ++i)
                result.matrixR(i,j) = mpfr::sec(matrixR(i,j));
    }

    return result;
}

/* cosecant function b = csc(a) */
GmpEigenMatrix GmpEigenMatrix::csc() const
{
    GmpEigenMatrix result;

    result.isComplex = isComplex;
    if (isComplex) {
        result.matrixR.resize(matrixR.rows(),matrixR.cols());
        result.matrixI.resize(matrixR.rows(),matrixR.cols());
        for (IndexType j(0); j < result.matrixR.cols(); ++j)
            for (IndexType i(0); i < result.matrixR.rows(); ++i) {
                result.matrixR(i,j) = 2*(mpfr::sin(matrixR(i,j))*mpfr::cosh(matrixI(i,j)))/(-mpfr::cos(2*matrixR(i,j)) + mpfr::cosh(2*matrixI(i,j)));
                result.matrixI(i,j) = -2*(mpfr::cos(matrixR(i,j))*mpfr::sinh(matrixI(i,j)))/(-mpfr::cos(2*matrixR(i,j)) + mpfr::cosh(2*matrixI(i,j)));
            }
        result.checkComplexity();
    } else {
        result.matrixR.resize(matrixR.rows(),matrixR.cols());
        for (IndexType j(0); j < result.matrixR.cols(); ++j)
            for (IndexType i(0); i < result.matrixR.rows(); ++i)
                result.matrixR(i,j) = mpfr::csc(matrixR(i,j));
    }

    return result;
}

/* cosecant function b = csc(a) */
GmpEigenMatrix& GmpEigenMatrix::csc_new() const
{
    GmpEigenMatrix& result(*(new GmpEigenMatrix));

    result.isComplex = isComplex;
    if (isComplex) {
        result.matrixR.resize(matrixR.rows(),matrixR.cols());
        result.matrixI.resize(matrixR.rows(),matrixR.cols());
        for (IndexType j(0); j < result.matrixR.cols(); ++j)
            for (IndexType i(0); i < result.matrixR.rows(); ++i) {
                result.matrixR(i,j) = 2*(mpfr::sin(matrixR(i,j))*mpfr::cosh(matrixI(i,j)))/(-mpfr::cos(2*matrixR(i,j)) + mpfr::cosh(2*matrixI(i,j)));
                result.matrixI(i,j) = -2*(mpfr::cos(matrixR(i,j))*mpfr::sinh(matrixI(i,j)))/(-mpfr::cos(2*matrixR(i,j)) + mpfr::cosh(2*matrixI(i,j)));
            }
        result.checkComplexity();
    } else {
        result.matrixR.resize(matrixR.rows(),matrixR.cols());
        for (IndexType j(0); j < result.matrixR.cols(); ++j)
            for (IndexType i(0); i < result.matrixR.rows(); ++i)
                result.matrixR(i,j) = mpfr::csc(matrixR(i,j));
    }

    return result;
}

/* cotangent function b = cot(a) */
GmpEigenMatrix GmpEigenMatrix::cot() const
{
    GmpEigenMatrix result;

    result.isComplex = isComplex;
    if (isComplex) {
        result.matrixR.resize(matrixR.rows(),matrixR.cols());
        result.matrixI.resize(matrixR.rows(),matrixR.cols());
        for (IndexType j(0); j < result.matrixR.cols(); ++j)
            for (IndexType i(0); i < result.matrixR.rows(); ++i) {
                result.matrixR(i,j) = (mpfr::sin(2*matrixR(i,j)))/(-mpfr::cos(2*matrixR(i,j)) + mpfr::cosh(2*matrixI(i,j)));
                result.matrixI(i,j) = (-mpfr::sinh(2*matrixI(i,j)))/(-mpfr::cos(2*matrixR(i,j)) + mpfr::cosh(2*matrixI(i,j)));
            }
        result.checkComplexity();
    } else {
        result.matrixR.resize(matrixR.rows(),matrixR.cols());
        for (IndexType j(0); j < result.matrixR.cols(); ++j)
            for (IndexType i(0); i < result.matrixR.rows(); ++i)
                result.matrixR(i,j) = mpfr::cot(matrixR(i,j));
    }

    return result;
}

/* cotangent function b = cot(a) */
GmpEigenMatrix& GmpEigenMatrix::cot_new() const
{
    GmpEigenMatrix& result(*(new GmpEigenMatrix));

    result.isComplex = isComplex;
    if (isComplex) {
        result.matrixR.resize(matrixR.rows(),matrixR.cols());
        result.matrixI.resize(matrixR.rows(),matrixR.cols());
        for (IndexType j(0); j < result.matrixR.cols(); ++j)
            for (IndexType i(0); i < result.matrixR.rows(); ++i) {
                result.matrixR(i,j) = (mpfr::sin(2*matrixR(i,j)))/(-mpfr::cos(2*matrixR(i,j)) + mpfr::cosh(2*matrixI(i,j)));
                result.matrixI(i,j) = (-mpfr::sinh(2*matrixI(i,j)))/(-mpfr::cos(2*matrixR(i,j)) + mpfr::cosh(2*matrixI(i,j)));
            }
        result.checkComplexity();
    } else {
        result.matrixR.resize(matrixR.rows(),matrixR.cols());
        for (IndexType j(0); j < result.matrixR.cols(); ++j)
            for (IndexType i(0); i < result.matrixR.rows(); ++i)
                result.matrixR(i,j) = mpfr::cot(matrixR(i,j));
    }

    return result;
}

/* arc sine function b = asin(a) */
GmpEigenMatrix GmpEigenMatrix::asin() const
{
    GmpEigenMatrix result;

    // We compute the analytic extension (valid for all complex numbers
    // including real numbers outside [-1,1])
    result = -constI()*(constI()*(*this) + (GmpEigenMatrix(1) - (*this).power(GmpEigenMatrix(2))).sqrt()).log();

    // Now we could also compute the images of real numbers between -1 and 1 if we want
    for (IndexType j(0); j < result.matrixR.cols(); ++j)
        for (IndexType i(0); i < result.matrixR.rows(); ++i)
            if (((!isComplex) || ((isComplex) && (matrixI(i,j) == 0))) && (matrixR(i,j) >= -1) && (matrixR(i,j) <= 1) && (result.isComplex)) {
//                result.matrixR(i,j) = mpfr::asin(matrixR(i,j));
                result.matrixI(i,j) = 0; // At least, we remove the imaginary parts which should be zero
            }
    result.checkComplexity();

    return result;
}

/* arc sine function b = asin(a) */
GmpEigenMatrix& GmpEigenMatrix::asin_new() const
{
    GmpEigenMatrix& result(*(new GmpEigenMatrix));

    // We compute the analytic extension (valid for all complex numbers
    // including real numbers outside [-1,1])
    result = -constI()*(constI()*(*this) + (GmpEigenMatrix(1) - (*this).power(GmpEigenMatrix(2))).sqrt()).log();

    // Now we could also compute the images of real numbers between -1 and 1 if we want
    for (IndexType j(0); j < result.matrixR.cols(); ++j)
        for (IndexType i(0); i < result.matrixR.rows(); ++i)
            if (((!isComplex) || ((isComplex) && (matrixI(i,j) == 0))) && (matrixR(i,j) >= -1) && (matrixR(i,j) <= 1) && (result.isComplex)) {
//                result.matrixR(i,j) = mpfr::asin(matrixR(i,j));
                result.matrixI(i,j) = 0;
            }
    result.checkComplexity();

    return result;
}

/* arc cosine function b = acos(a) */
GmpEigenMatrix GmpEigenMatrix::acos() const
{
    // We rely on the formula acos(x) = 1/2*(pi-2*asin(x))
    return GmpEigenMatrix(mpreal("0.5"))*(constPi() - GmpEigenMatrix(2)*(*this).asin());
}

/* arc cosine function b = acos(a) */
GmpEigenMatrix& GmpEigenMatrix::acos_new() const
{
    GmpEigenMatrix& result(*(new GmpEigenMatrix));

    // We rely on the formula acos(x) = 1/2*(pi-2*asin(x))
    result = GmpEigenMatrix(mpreal("0.5"))*(constPi() - GmpEigenMatrix(2)*(*this).asin());

    return result;
}

/* arc tangent function b = atan(a) */
GmpEigenMatrix GmpEigenMatrix::atan() const
{
    GmpEigenMatrix result;

    result.isComplex = isComplex;
    if (isComplex) {
        // We compute the analytic extension (valid for all complex numbers)
        result = GmpEigenMatrix(0,mpreal("0.5"))*((GmpEigenMatrix(1)-constI()*(*this)).log() - (GmpEigenMatrix(1)+constI()*(*this)).log());

        result.checkComplexity();
    } else {
        result.matrixR.resize(matrixR.rows(),matrixR.cols());
        for (IndexType j(0); j < result.matrixR.cols(); ++j)
            for (IndexType i(0); i < result.matrixR.rows(); ++i)
                result.matrixR(i,j) = mpfr::atan(matrixR(i,j));
    }

    return result;
}

/* arc tangent function b = atan(a) */
GmpEigenMatrix& GmpEigenMatrix::atan_new() const
{
    GmpEigenMatrix& result(*(new GmpEigenMatrix));

    result.isComplex = isComplex;
    if (isComplex) {
        // We compute the analytic extension (valid for all complex numbers)
        result = GmpEigenMatrix(0,mpreal("0.5"))*((GmpEigenMatrix(1)-constI()*(*this)).log() - (GmpEigenMatrix(1)+constI()*(*this)).log());

        result.checkComplexity();
    } else {
        result.matrixR.resize(matrixR.rows(),matrixR.cols());
        for (IndexType j(0); j < result.matrixR.cols(); ++j)
            for (IndexType i(0); i < result.matrixR.rows(); ++i)
                result.matrixR(i,j) = mpfr::atan(matrixR(i,j));
    }

    return result;
}








/* hyperbolic sine function b = sinh(a) */
GmpEigenMatrix GmpEigenMatrix::sinh() const
{
    GmpEigenMatrix result;

    result.matrixR.resize(matrixR.rows(),matrixR.cols());
    result.isComplex = isComplex;
    if (isComplex) {
        result.matrixI.resize(matrixR.rows(),matrixR.cols());
//        result = -constI()*GmpEigenMatrix(1/mpreal(2))*((constI()*(*this)).exp() - (-constI()*(*this)).exp());
        for (IndexType j(0); j < result.matrixR.cols(); ++j)
            for (IndexType i(0); i < result.matrixR.rows(); ++i) {
                result.matrixR(i,j) = mpfr::sinh(matrixR(i,j))*mpfr::cos(matrixI(i,j));
                result.matrixI(i,j) = mpfr::cosh(matrixR(i,j))*mpfr::sin(matrixI(i,j));
            }
        result.checkComplexity();
    } else {
        for (IndexType j(0); j < result.matrixR.cols(); ++j)
            for (IndexType i(0); i < result.matrixR.rows(); ++i)
                result.matrixR(i,j) = mpfr::sinh(matrixR(i,j));
    }

    return result;
}

/* hyperbolic sine function b = sinh(a) */
GmpEigenMatrix& GmpEigenMatrix::sinh_new() const
{
    GmpEigenMatrix& result(*(new GmpEigenMatrix));

    result.matrixR.resize(matrixR.rows(),matrixR.cols());
    result.isComplex = isComplex;
    if (isComplex) {
        result.matrixI.resize(matrixR.rows(),matrixR.cols());
//        result = -constI()*GmpEigenMatrix(1/mpreal(2))*((constI()*(*this)).exp() - (-constI()*(*this)).exp());
        for (IndexType j(0); j < result.matrixR.cols(); ++j)
            for (IndexType i(0); i < result.matrixR.rows(); ++i) {
                result.matrixR(i,j) = mpfr::sinh(matrixR(i,j))*mpfr::cos(matrixI(i,j));
                result.matrixI(i,j) = mpfr::cosh(matrixR(i,j))*mpfr::sin(matrixI(i,j));
            }
        result.checkComplexity();
    } else {
        for (IndexType j(0); j < result.matrixR.cols(); ++j)
            for (IndexType i(0); i < result.matrixR.rows(); ++i)
                result.matrixR(i,j) = mpfr::sinh(matrixR(i,j));
    }

    return result;
}

/* hyperbolic cosine function b = cosh(a) */
GmpEigenMatrix GmpEigenMatrix::cosh() const
{
    GmpEigenMatrix result;

    result.matrixR.resize(matrixR.rows(),matrixR.cols());
    result.isComplex = isComplex;
    if (isComplex) {
        result.matrixI.resize(matrixR.rows(),matrixR.cols());
        for (IndexType j(0); j < result.matrixR.cols(); ++j)
            for (IndexType i(0); i < result.matrixR.rows(); ++i) {
                result.matrixR(i,j) = mpfr::cosh(matrixR(i,j))*mpfr::cos(matrixI(i,j));
                result.matrixI(i,j) = mpfr::sinh(matrixR(i,j))*mpfr::sin(matrixI(i,j));
            }
        result.checkComplexity();
    } else {
        for (IndexType j(0); j < result.matrixR.cols(); ++j)
            for (IndexType i(0); i < result.matrixR.rows(); ++i)
                result.matrixR(i,j) = mpfr::cosh(matrixR(i,j));
    }

    return result;
}

/* hyperbolic cosine function b = cosh(a) */
GmpEigenMatrix& GmpEigenMatrix::cosh_new() const
{
    GmpEigenMatrix& result(*(new GmpEigenMatrix));

    result.matrixR.resize(matrixR.rows(),matrixR.cols());
    result.isComplex = isComplex;
    if (isComplex) {
        result.matrixI.resize(matrixR.rows(),matrixR.cols());
        for (IndexType j(0); j < result.matrixR.cols(); ++j)
            for (IndexType i(0); i < result.matrixR.rows(); ++i) {
                result.matrixR(i,j) = mpfr::cosh(matrixR(i,j))*mpfr::cos(matrixI(i,j));
                result.matrixI(i,j) = mpfr::sinh(matrixR(i,j))*mpfr::sin(matrixI(i,j));
            }
        result.checkComplexity();
    } else {
        for (IndexType j(0); j < result.matrixR.cols(); ++j)
            for (IndexType i(0); i < result.matrixR.rows(); ++i)
                result.matrixR(i,j) = mpfr::cosh(matrixR(i,j));
    }

    return result;
}

/* For tanh, use the fact that tan(ix) = i*tahn(x), plus the expansion of tanh in terms of sinh and cosh. */






/* ---------------------------------------
   |   Some linear algebra operations    |
   --------------------------------------- */

/* Matrix rank b = rank(a) */
IndexType GmpEigenMatrix::rank() const
{
   IndexType result;

   if (isComplex) {
       result = complexIsometry().rank()/2;
   } else {
       // We use a QR solver to compute the matrix rank

       //matrixR.makeCompressed(); // The matrix is supposed to be alread compressed at this stage...
       ColPivHouseholderQR< Matrix<mpreal, Dynamic, Dynamic> > solver(matrixR);

       result = solver.rank();
   }

   return result;
}

/* Linear system solving c = a\b*/
GmpEigenMatrix GmpEigenMatrix::mldivide(const GmpEigenMatrix& b) const
{
    GmpEigenMatrix result;

    if ((isComplex) || (b.isComplex)) {
        result = complexIsometry().mldivide(b.complexIsometry()).complexIsometryInverse();
    } else {
        result.isComplex = false;
        result.matrixR = matrixR.colPivHouseholderQr().solve(b.matrixR);
        result.matrixI.resize(0,0);
    }

    return result;
}

/* Linear system solving c = a\b*/
GmpEigenMatrix& GmpEigenMatrix::mldivide_new(const GmpEigenMatrix& b) const
{
    GmpEigenMatrix& result(*(new GmpEigenMatrix));

    if ((isComplex) || (b.isComplex)) {
        result = complexIsometry().mldivide(b.complexIsometry()).complexIsometryInverse();
    } else {
        result.isComplex = false;
        result.matrixR = matrixR.colPivHouseholderQr().solve(b.matrixR);
        result.matrixI.resize(0,0);
    }

    return result;
}

/* Matrix inverse b = inv(a) */
GmpEigenMatrix GmpEigenMatrix::inv() const
{
    GmpEigenMatrix result;

    if (isComplex) {
        result = complexIsometry().inv().complexIsometryInverse();
    } else {
        result.matrixR = matrixR.inverse();
    }

    return result;
}

/* Matrix inverse b = inv(a) */
GmpEigenMatrix& GmpEigenMatrix::inv_new() const
{
    GmpEigenMatrix& result(*(new GmpEigenMatrix));

    if (isComplex) {
        result = complexIsometry().inv().complexIsometryInverse();
    } else {
        result.matrixR = matrixR.inverse();
    }

    return result;
}

/* Eigen decopmosition: similar to [V D] = eig(A) */
GmpEigenMatrix GmpEigenMatrix::eig(GmpEigenMatrix& V) const
{
    GmpEigenMatrix result;

    if (ishermitian()) {
        if (isComplex) {
            IndexType size(matrixR.rows());

            // First, we solve a real version of the problem
            GmpEigenMatrix Dreal, Vreal; // These matrices will hold an isometry of the eigenvalues and eigenvectors
            SelfAdjointEigenSolver< Matrix<mpreal,Dynamic,Dynamic> > es(complexIsometry().matrixR);

            Dreal.isComplex = false;
            Dreal.matrixR = es.eigenvalues();
            Dreal.matrixI.resize(0,0);

            Vreal.isComplex = false;
            Vreal.matrixR = es.eigenvectors();
            Vreal.matrixI.resize(0,0);

            // Since the eigenvalues are sorted in increasing order, we can easily remove duplicates
            // So we select the eigenvalues and form the corresponding eigenvectors
            result.isComplex = false;
            result.matrixR.setZero(size,size);
            result.matrixI.resize(0,0);
            V.matrixR.setZero(size,size);
            V.matrixI.setZero(size,size);
            for (IndexType i(0); i < size; ++i) {
                result.matrixR(i,i) = Dreal.matrixR(2*i,0);
                // This way of extracting the eigenvector always works ;-)
                V.matrixR.block(0,i,size,1) = Vreal.matrixR.block(0,2*i,size,1);
                V.matrixI.block(0,i,size,1) = -Vreal.matrixR.block(size,2*i,size,1);
            }
            V.checkComplexity();
        } else {
            SelfAdjointEigenSolver< Matrix<mpreal,Dynamic,Dynamic> > es(matrixR);

            result.isComplex = false;
            result.matrixR = es.eigenvalues();
            result = result.diagCreate(0);

            V.isComplex = false;
            V.matrixR = es.eigenvectors();
            V.matrixI.resize(0,0);
        }
    } else {
        if (isComplex) {
            // NOTE : According to http://www.holoborodko.com/pavel/mpfr/ , using
            // std::complex<mpfr::mpreal> can lead to dramatic losses of precision
            // in the eigenvalue computation, because Eigen's algorithm uses
            // std::abs(std::complex) on those numbers, which are then implicitely
            // transformed as doubles. (c.f. comment from April 20, 2016)
            // But we can avoid std::complex!

            IndexType size(matrixR.rows());

            // First, we solve a real version of the problem
            GmpEigenMatrix Dreal, Vreal; // These matrices will hold an isometry of the eigenvalues and eigenvectors
            EigenSolver< Matrix<mpreal,Dynamic,Dynamic> > es(complexIsometry().matrixR);

            Dreal.isComplex = false;
            Dreal.matrixR = es.pseudoEigenvalueMatrix();
            Dreal.matrixI.resize(0,0);

            Vreal.isComplex = false;
            Vreal.matrixR = es.pseudoEigenvectors();
            Vreal.matrixI.resize(0,0);

            // Now we analyse the result and extract the complex eigenvalues and eigenvectors.
            // Note that the solver does not provide the eigenvalues in an ordered
            // way. So we cannot trust the order of the columns of Vrel.
            // However, the solver always groups together in a 2x2 block any couple
            // of complex eigenvalues. The only order we need to check for them is
            // thus whether the first coefficient of the block is real and the
            // second imaginary, or the reverse, i.e. whether the eigenvectors
            // column pairs are of the form [r i; -i r] or [i r; r -i].
            // Concerning the real eigenvalues, however, random order of the columns
            // means that in order to remove the double real eigenvalues, we need to
            // sort them out.

            // So we start by identifying the complex 2x2 blocks and the real 1x1 blocks
            GmpEigenMatrix eigenvalues;
            eigenvalues.isComplex = true;
            eigenvalues.matrixR.resize(size,1);
            eigenvalues.matrixI.resize(size,1);
            V.isComplex = true;
            V.matrixR.resize(size, size);
            V.matrixI.resize(size, size);
            vector < pair < mpreal, IndexType > > realEigenvaluesAndColumns;
            IndexType i(0), nbEigsWritten(0);
            mpreal eps(10*mpfr::pow(10, -mpfr::bits2digits(mpreal::get_default_prec()))); // This is the critical point at which we judge whether an eigenvector is complex or not
            while (i < Vreal.matrixR.cols()) {
                if ((i == Vreal.matrixR.cols()-1) || (Dreal.matrixR(i,i+1) == 0)
                || ((Dreal.matrixR(i,i+1) < eps) && (((Vreal.matrixR.block(0,i,size,1).array().abs() - Vreal.matrixR.block(size,i+1,size,1).array().abs()).abs().maxCoeff() > eps) || ((Vreal.matrixR.block(0,i+1,size,1).array().abs() - Vreal.matrixR.block(size,i,size,1).array().abs()).abs().maxCoeff() > eps)))) {
                    // The next eigenvalue is real
                    realEigenvaluesAndColumns.push_back(make_pair(Dreal.matrixR(i,i), i));
                    ++i;
                } else {
                    // The next eigenvalue is complex
                    // We directly take care of this case
                    if ( (Vreal.matrixR.block(0,i,size,1)+Vreal.matrixR.block(size,i+1,size,1)).array().abs().maxCoeff() >= (Vreal.matrixR.block(0,i+1,size,1)+Vreal.matrixR.block(size,i,size,1)).array().abs().maxCoeff() ) {
                        // We are in the [r i; -i r] case
                        eigenvalues.matrixR(nbEigsWritten,0) = Dreal.matrixR(i,i);
                        eigenvalues.matrixI(nbEigsWritten,0) = Dreal.matrixR(i,i+1);
                        V.matrixR.block(0,nbEigsWritten,size,1) = Vreal.matrixR.block(0,i,size,1);
                        V.matrixI.block(0,nbEigsWritten,size,1) = Vreal.matrixR.block(0,i+1,size,1);
                    } else {
                        // We are in the [i r; r -i] case
                        eigenvalues.matrixR(nbEigsWritten,0) = Dreal.matrixR(i,i);
                        eigenvalues.matrixI(nbEigsWritten,0) = Dreal.matrixR(i+1,i); // this is the same as -result(i,i+1)
                        V.matrixR.block(0,nbEigsWritten,size,1) = Vreal.matrixR.block(0,i+1,size,1);
                        V.matrixI.block(0,nbEigsWritten,size,1) = Vreal.matrixR.block(0,i,size,1);
                    }
                    ++nbEigsWritten;
                    ++i; ++i; // iterate by two columns
                }
            }

            // Now we add the real eigenvalues and the corresponding eigenvectors
            // First, we sort them
            std::sort(realEigenvaluesAndColumns.begin(),realEigenvaluesAndColumns.end());
            // Now we can extract the values and eigenvectors one by one
            i = 0;
            while (i < realEigenvaluesAndColumns.size()) {
                eigenvalues.matrixR(nbEigsWritten,0) = realEigenvaluesAndColumns[i].first;
                // This way of extracting the eigenvector always works ;-)
                V.matrixR.block(0,nbEigsWritten,size,1) = Vreal.matrixR.block(0,realEigenvaluesAndColumns[i].second,size,1);
                V.matrixI.block(0,nbEigsWritten,size,1) = -Vreal.matrixR.block(size,realEigenvaluesAndColumns[i].second,size,1);
                ++nbEigsWritten;
                ++i; ++i; // We iterate eigenvalues by pairs
            }
            eigenvalues.checkComplexity();
            V.checkComplexity();

            // We can now put the eigenvalues in a matrix
            result = eigenvalues.diagCreate(0);
        } else {
            EigenSolver< Matrix<mpreal,Dynamic,Dynamic> > es(matrixR);

            result.isComplex = false;
            result.matrixR = es.pseudoEigenvalueMatrix();

            V.isComplex = false;
            V.matrixR = es.pseudoEigenvectors();
            V.matrixI.resize(0,0);

            // Now, we take care of complex eigenvectors, if there are any.
            // Their presence is signified by 2x2 blocks on the pseudo eigenvalue
            // matrix.
            vector < IndexType > indicesAll, indicesBlock, indicesBlockP1, indicesFullBlock;
            for (IndexType i(0); i < matrixR.rows(); ++i) {
                indicesAll.push_back(i);
                if ((i < matrixR.rows()-1) && (result.matrixR(i,i+1) != 0)) {
                    indicesBlock.push_back(i);
                    indicesBlockP1.push_back(i+1);
                    indicesFullBlock.push_back(i);
                    indicesFullBlock.push_back(i+1);
                }
            }
            V.subsasgn(indicesAll, indicesBlock, V.subsref(indicesAll, indicesBlock) + V.subsref(indicesAll, indicesBlockP1).times(constI()));
            V.subsasgn(indicesAll, indicesBlockP1, V.subsref(indicesAll, indicesBlock).conj());

            // Now we also recombine the eigenvalues
            if (indicesBlock.size() > 0) {
                result.isComplex = true;
                result.matrixI.resize(result.matrixR.rows(), result.matrixR.cols());
                for (IndexType i(0); i < indicesBlock.size(); ++i) {
                    IndexType j(indicesBlock[i]);
                    result.matrixI(j,j) = result.matrixR(j,j+1);
                    result.matrixI(j+1,j+1) = result.matrixR(j+1,j);
                    result.matrixR(j,j+1) = 0;
                    result.matrixR(j+1,j) = 0;
                }
            }
        }
    }

    return result;
}

/* Eigen decopmosition: similar to [V D] = eig(A) */
GmpEigenMatrix& GmpEigenMatrix::eig_new(GmpEigenMatrix& V) const
{
    GmpEigenMatrix& result(*(new GmpEigenMatrix));

    if (ishermitian()) {
        if (isComplex) {
            IndexType size(matrixR.rows());

            // First, we solve a real version of the problem
            GmpEigenMatrix Dreal, Vreal; // These matrices will hold an isometry of the eigenvalues and eigenvectors
            SelfAdjointEigenSolver< Matrix<mpreal,Dynamic,Dynamic> > es(complexIsometry().matrixR);

            Dreal.isComplex = false;
            Dreal.matrixR = es.eigenvalues();
            Dreal.matrixI.resize(0,0);

            Vreal.isComplex = false;
            Vreal.matrixR = es.eigenvectors();
            Vreal.matrixI.resize(0,0);

            // Since the eigenvalues are sorted in increasing order, we can easily remove duplicates
            // So we select the eigenvalues and form the corresponding eigenvectors
            result.isComplex = false;
            result.matrixR.setZero(size,size);
            result.matrixI.resize(0,0);
            V.matrixR.setZero(size,size);
            V.matrixI.setZero(size,size);
            for (IndexType i(0); i < size; ++i) {
                result.matrixR(i,i) = Dreal.matrixR(2*i,0);
                // This way of extracting the eigenvector always works ;-)
                V.matrixR.block(0,i,size,1) = Vreal.matrixR.block(0,2*i,size,1);
                V.matrixI.block(0,i,size,1) = -Vreal.matrixR.block(size,2*i,size,1);
            }
            V.checkComplexity();
        } else {
            SelfAdjointEigenSolver< Matrix<mpreal,Dynamic,Dynamic> > es(matrixR);

            result.isComplex = false;
            result.matrixR = es.eigenvalues();
            result = result.diagCreate(0);

            V.isComplex = false;
            V.matrixR = es.eigenvectors();
            V.matrixI.resize(0,0);
        }
    } else {
        if (isComplex) {
            // NOTE : According to http://www.holoborodko.com/pavel/mpfr/ , using
            // std::complex<mpfr::mpreal> can lead to dramatic losses of precision
            // in the eigenvalue computation, because Eigen's algorithm uses
            // std::abs(std::complex) on those numbers, which are then implicitely
            // transformed as doubles. (c.f. comment from April 20, 2016)
            // But we can avoid std::complex!

            IndexType size(matrixR.rows());

            // First, we solve a real version of the problem
            GmpEigenMatrix Dreal, Vreal; // These matrices will hold an isometry of the eigenvalues and eigenvectors
            EigenSolver< Matrix<mpreal,Dynamic,Dynamic> > es(complexIsometry().matrixR);

            Dreal.isComplex = false;
            Dreal.matrixR = es.pseudoEigenvalueMatrix();
            Dreal.matrixI.resize(0,0);

            Vreal.isComplex = false;
            Vreal.matrixR = es.pseudoEigenvectors();
            Vreal.matrixI.resize(0,0);

            // Now we analyse the result and extract the complex eigenvalues and eigenvectors.
            // Note that the solver does not provide the eigenvalues in an ordered
            // way. So we cannot trust the order of the columns of Vrel.
            // However, the solver always groups together in a 2x2 block any couple
            // of complex eigenvalues. The only order we need to check for them is
            // thus whether the first coefficient of the block is real and the
            // second imaginary, or the reverse, i.e. whether the eigenvectors
            // column pairs are of the form [r i; -i r] or [i r; r -i].
            // Concerning the real eigenvalues, however, random order of the columns
            // means that in order to remove the double real eigenvalues, we need to
            // sort them out.

            // So we start by identifying the complex 2x2 blocks and the real 1x1 blocks
            GmpEigenMatrix eigenvalues;
            eigenvalues.isComplex = true;
            eigenvalues.matrixR.resize(size,1);
            eigenvalues.matrixI.resize(size,1);
            V.isComplex = true;
            V.matrixR.resize(size, size);
            V.matrixI.resize(size, size);
            vector < pair < mpreal, IndexType > > realEigenvaluesAndColumns;
            IndexType i(0), nbEigsWritten(0);
            mpreal eps(10*mpfr::pow(10, -mpfr::bits2digits(mpreal::get_default_prec()))); // This is the critical point at which we judge whether an eigenvector is complex or not
            while (i < Vreal.matrixR.cols()) {
                if ((i == Vreal.matrixR.cols()-1) || (Dreal.matrixR(i,i+1) == 0)
                || ((Dreal.matrixR(i,i+1) < eps) && (((Vreal.matrixR.block(0,i,size,1).array().abs() - Vreal.matrixR.block(size,i+1,size,1).array().abs()).abs().maxCoeff() > eps) || ((Vreal.matrixR.block(0,i+1,size,1).array().abs() - Vreal.matrixR.block(size,i,size,1).array().abs()).abs().maxCoeff() > eps)))) {
                    // The next eigenvalue is real
                    realEigenvaluesAndColumns.push_back(make_pair(Dreal.matrixR(i,i), i));
                    ++i;
                } else {
                    // The next eigenvalue is complex
                    // We directly take care of this case
                    if ( (Vreal.matrixR.block(0,i,size,1)+Vreal.matrixR.block(size,i+1,size,1)).array().abs().maxCoeff() >= (Vreal.matrixR.block(0,i+1,size,1)+Vreal.matrixR.block(size,i,size,1)).array().abs().maxCoeff() ) {
                        // We are in the [r i; -i r] case
                        eigenvalues.matrixR(nbEigsWritten,0) = Dreal.matrixR(i,i);
                        eigenvalues.matrixI(nbEigsWritten,0) = Dreal.matrixR(i,i+1);
                        V.matrixR.block(0,nbEigsWritten,size,1) = Vreal.matrixR.block(0,i,size,1);
                        V.matrixI.block(0,nbEigsWritten,size,1) = Vreal.matrixR.block(0,i+1,size,1);
                    } else {
                        // We are in the [i r; r -i] case
                        eigenvalues.matrixR(nbEigsWritten,0) = Dreal.matrixR(i,i);
                        eigenvalues.matrixI(nbEigsWritten,0) = Dreal.matrixR(i+1,i); // this is the same as -result(i,i+1)
                        V.matrixR.block(0,nbEigsWritten,size,1) = Vreal.matrixR.block(0,i+1,size,1);
                        V.matrixI.block(0,nbEigsWritten,size,1) = Vreal.matrixR.block(0,i,size,1);
                    }
                    ++nbEigsWritten;
                    ++i; ++i; // iterate by two columns
                }
            }

            // Now we add the real eigenvalues and the corresponding eigenvectors
            // First, we sort them
            std::sort(realEigenvaluesAndColumns.begin(),realEigenvaluesAndColumns.end());
            // Now we can extract the values and eigenvectors one by one
            i = 0;
            while (i < realEigenvaluesAndColumns.size()) {
                eigenvalues.matrixR(nbEigsWritten,0) = realEigenvaluesAndColumns[i].first;
                // This way of extracting the eigenvector always works ;-)
                V.matrixR.block(0,nbEigsWritten,size,1) = Vreal.matrixR.block(0,realEigenvaluesAndColumns[i].second,size,1);
                V.matrixI.block(0,nbEigsWritten,size,1) = -Vreal.matrixR.block(size,realEigenvaluesAndColumns[i].second,size,1);
                ++nbEigsWritten;
                ++i; ++i; // We iterate eigenvalues by pairs
            }
            eigenvalues.checkComplexity();
            V.checkComplexity();

            // We can now put the eigenvalues in a matrix
            result = eigenvalues.diagCreate(0);
        } else {
            EigenSolver< Matrix<mpreal,Dynamic,Dynamic> > es(matrixR);

            result.isComplex = false;
            result.matrixR = es.pseudoEigenvalueMatrix();

            V.isComplex = false;
            V.matrixR = es.pseudoEigenvectors();
            V.matrixI.resize(0,0);

            // Now, we take care of complex eigenvectors, if there are any.
            // Their presence is signified by 2x2 blocks on the pseudo eigenvalue
            // matrix.
            vector < IndexType > indicesAll, indicesBlock, indicesBlockP1, indicesFullBlock;
            for (IndexType i(0); i < matrixR.rows(); ++i) {
                indicesAll.push_back(i);
                if ((i < matrixR.rows()-1) && (result.matrixR(i,i+1) != 0)) {
                    indicesBlock.push_back(i);
                    indicesBlockP1.push_back(i+1);
                    indicesFullBlock.push_back(i);
                    indicesFullBlock.push_back(i+1);
                }
            }
            V.subsasgn(indicesAll, indicesBlock, V.subsref(indicesAll, indicesBlock) + V.subsref(indicesAll, indicesBlockP1).times(constI()));
            V.subsasgn(indicesAll, indicesBlockP1, V.subsref(indicesAll, indicesBlock).conj());

            // Now we also recombine the eigenvalues
            if (indicesBlock.size() > 0) {
                result.isComplex = true;
                result.matrixI.resize(result.matrixR.rows(), result.matrixR.cols());
                for (IndexType i(0); i < indicesBlock.size(); ++i) {
                    IndexType j(indicesBlock[i]);
                    result.matrixI(j,j) = result.matrixR(j,j+1);
                    result.matrixI(j+1,j+1) = result.matrixR(j+1,j);
                    result.matrixR(j,j+1) = 0;
                    result.matrixR(j+1,j) = 0;
                }
            }
        }
    }

    return result;
}

/* Partial eigenvalue decopmosition: similar to [V D] = eigs(A, k) */
GmpEigenMatrix GmpEigenMatrix::eigs(const long int& nbEigenvalues, GmpEigenMatrix& V, const long int& type, const GmpEigenMatrix& sigma) const
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
            // Construct matrix operation object using the wrapper class
            Spectra::DenseSymMatProd<mpreal> op(matrixR);

            //long int ncv(min(max(1+nbEigenvalues,2*nbEigenvalues),matrixR.rows())); // the tightest bounds... don't always converge
            long int ncv(min(3+max(1+nbEigenvalues,2*nbEigenvalues),matrixR.rows()));
            switch (type) {
                case 1: {
                    // Construct eigen solver object, requesting desired eigenvalues
                    Spectra::SymEigsSolver< mpreal, Spectra::LARGEST_MAGN, Spectra::DenseSymMatProd<mpreal> > eigs(&op, nbEigenvalues, ncv);

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
                    Spectra::DenseSymShiftSolve<mpreal> op(matrixR);

                    // Construct eigen solver object, requesting desired eigenvalues
                    Spectra::SymEigsShiftSolver< mpreal, Spectra::LARGEST_MAGN, Spectra::DenseSymShiftSolve<mpreal> > eigs(&op, nbEigenvalues, ncv, sigma.matrixR(0,0));

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
            //long int ncv(min(max(2+nbEigenvalues,2*nbEigenvalues),matrixR.rows()));  // the tightest bounds... don't always converge
            long int ncv(min(3+max(1+nbEigenvalues,2*nbEigenvalues),matrixR.rows()));
            switch (type) {
                case 1: {
                    // Construct matrix operation object using the wrapper class
                    Spectra::DenseGenMatProd<mpreal> op(matrixR);

                    // Construct eigen solver object, requesting desired eigenvalues
                    Spectra::GenEigsSolver< mpreal, Spectra::LARGEST_MAGN, Spectra::DenseGenMatProd<mpreal> > eigs(&op, nbEigenvalues, ncv);

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
                    Spectra::DenseGenComplexShiftSolve<mpreal> op(matrixR);

                    // Construct eigen solver object, requesting desired eigenvalues
                    mpreal sigmaI(0);
                    if (sigma.isComplex)
                        sigmaI = sigma.matrixI(0,0);
                    Spectra::GenEigsComplexShiftSolver< mpreal, Spectra::LARGEST_MAGN, Spectra::DenseGenComplexShiftSolve<mpreal> > eigs(&op, nbEigenvalues, ncv, sigma.matrixR(0,0), sigmaI);

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

/* Partial eigenvalue decopmosition: similar to [V D] = eigs(A, k) */
GmpEigenMatrix& GmpEigenMatrix::eigs_new(const long int& nbEigenvalues, GmpEigenMatrix& V, const long int& type, const GmpEigenMatrix& sigma) const
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
            // Construct matrix operation object using the wrapper class
            Spectra::DenseSymMatProd<mpreal> op(matrixR);

            //long int ncv(min(max(1+nbEigenvalues,2*nbEigenvalues),matrixR.rows())); // the tightest bounds... don't always converge
            long int ncv(min(3+max(1+nbEigenvalues,2*nbEigenvalues),matrixR.rows()));
            switch (type) {
                case 1: {
                    // Construct eigen solver object, requesting desired eigenvalues
                    Spectra::SymEigsSolver< mpreal, Spectra::LARGEST_MAGN, Spectra::DenseSymMatProd<mpreal> > eigs(&op, nbEigenvalues, ncv);

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
                    Spectra::DenseSymShiftSolve<mpreal> op(matrixR);

                    // Construct eigen solver object, requesting desired eigenvalues
                    Spectra::SymEigsShiftSolver< mpreal, Spectra::LARGEST_MAGN, Spectra::DenseSymShiftSolve<mpreal> > eigs(&op, nbEigenvalues, ncv, sigma.matrixR(0,0));

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
            //long int ncv(min(max(2+nbEigenvalues,2*nbEigenvalues),matrixR.rows()));  // the tightest bounds... don't always converge
            long int ncv(min(3+max(1+nbEigenvalues,2*nbEigenvalues),matrixR.rows()));
            switch (type) {
                case 1: {
                    // Construct matrix operation object using the wrapper class
                    Spectra::DenseGenMatProd<mpreal> op(matrixR);

                    // Construct eigen solver object, requesting desired eigenvalues
                    Spectra::GenEigsSolver< mpreal, Spectra::LARGEST_MAGN, Spectra::DenseGenMatProd<mpreal> > eigs(&op, nbEigenvalues, ncv);

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
                    Spectra::DenseGenComplexShiftSolve<mpreal> op(matrixR);

                    // Construct eigen solver object, requesting desired eigenvalues
                    mpreal sigmaI(0);
                    if (sigma.isComplex)
                        sigmaI = sigma.matrixI(0,0);
                    Spectra::GenEigsComplexShiftSolver< mpreal, Spectra::LARGEST_MAGN, Spectra::DenseGenComplexShiftSolve<mpreal> > eigs(&op, nbEigenvalues, ncv, sigma.matrixR(0,0), sigmaI);

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

/* Singular value decopmosition: similar to [U S V] = svd(A) */
GmpEigenMatrix GmpEigenMatrix::svd(GmpEigenMatrix& U, GmpEigenMatrix& V) const
{
    GmpEigenMatrix result;

    if (isComplex) {
        // First, we solve a real version of the problem
        GmpEigenMatrix Sreal, Ureal, Vreal; // These matrices will hold an isometry of the eigenvalues and eigenvectors
        JacobiSVD< Matrix<mpreal,Dynamic,Dynamic> > svd(complexIsometry().matrixR, ComputeThinU | ComputeThinV);

        Sreal.isComplex = false;
        Sreal.matrixR = svd.singularValues();
        Sreal.matrixI.resize(0,0);

        Ureal.isComplex = false;
        Ureal.matrixR = svd.matrixU();
        Ureal.matrixI.resize(0,0);

        Vreal.isComplex = false;
        Vreal.matrixR = svd.matrixV();
        Vreal.matrixI.resize(0,0);

        // Since the singular values are sorted in decreasing order, we can easily remove duplicates
        // So we select the singular values and form the corresponding singular vectors
        IndexType sizeS(Sreal.matrixR.rows()), sizeU(Ureal.matrixR.rows()), sizeV(Vreal.matrixR.rows());
        result.isComplex = false;
        result.matrixR.setZero(sizeS/2,sizeS/2);
        result.matrixI.resize(0,0);
        U.matrixR.setZero(sizeU/2,sizeS/2);
        U.matrixI.setZero(sizeU/2,sizeS/2);
        V.matrixR.setZero(sizeV/2,sizeS/2);
        V.matrixI.setZero(sizeV/2,sizeS/2);
        for (IndexType i(0); i < sizeS/2; ++i) {
            result.matrixR(i,i) = Sreal.matrixR(2*i,0);
            // This way of extracting the singular vectors always works ;-)
            U.matrixR.block(0,i,sizeU/2,1) = Ureal.matrixR.block(0,2*i,sizeU/2,1);
            U.matrixI.block(0,i,sizeU/2,1) = -Ureal.matrixR.block(sizeU/2,2*i,sizeU/2,1);
            V.matrixR.block(0,i,sizeV/2,1) = Vreal.matrixR.block(0,2*i,sizeV/2,1);
            V.matrixI.block(0,i,sizeV/2,1) = -Vreal.matrixR.block(sizeV/2,2*i,sizeV/2,1);
        }
        U.checkComplexity();
        V.checkComplexity();
    } else {
        JacobiSVD< Matrix<mpreal,Dynamic,Dynamic> > svd(matrixR, ComputeThinU | ComputeThinV);

        result.isComplex = false;
        result.matrixR = svd.singularValues();
        result = result.diagCreate(0);

        U.isComplex = false;
        U.matrixR = svd.matrixU();
        U.matrixI.resize(0,0);

        V.isComplex = false;
        V.matrixR = svd.matrixV();
        V.matrixI.resize(0,0);
    }

    return result;
}

/* Singular value decopmosition: similar to [U S V] = svd(A) */
GmpEigenMatrix& GmpEigenMatrix::svd_new(GmpEigenMatrix& U, GmpEigenMatrix& V) const
{
    GmpEigenMatrix& result(*(new GmpEigenMatrix));

    if (isComplex) {
        // First, we solve a real version of the problem
        GmpEigenMatrix Sreal, Ureal, Vreal; // These matrices will hold an isometry of the eigenvalues and eigenvectors
        JacobiSVD< Matrix<mpreal,Dynamic,Dynamic> > svd(complexIsometry().matrixR, ComputeThinU | ComputeThinV);

        Sreal.isComplex = false;
        Sreal.matrixR = svd.singularValues();
        Sreal.matrixI.resize(0,0);

        Ureal.isComplex = false;
        Ureal.matrixR = svd.matrixU();
        Ureal.matrixI.resize(0,0);

        Vreal.isComplex = false;
        Vreal.matrixR = svd.matrixV();
        Vreal.matrixI.resize(0,0);

        // Since the singular values are sorted in decreasing order, we can easily remove duplicates
        // So we select the singular values and form the corresponding singular vectors
        IndexType sizeS(Sreal.matrixR.rows()), sizeU(Ureal.matrixR.rows()), sizeV(Vreal.matrixR.rows());
        result.isComplex = false;
        result.matrixR.setZero(sizeS/2,sizeS/2);
        result.matrixI.resize(0,0);
        U.matrixR.setZero(sizeU/2,sizeS/2);
        U.matrixI.setZero(sizeU/2,sizeS/2);
        V.matrixR.setZero(sizeV/2,sizeS/2);
        V.matrixI.setZero(sizeV/2,sizeS/2);
        for (IndexType i(0); i < sizeS/2; ++i) {
            result.matrixR(i,i) = Sreal.matrixR(2*i,0);
            // This way of extracting the singular vectors always works ;-)
            U.matrixR.block(0,i,sizeU/2,1) = Ureal.matrixR.block(0,2*i,sizeU/2,1);
            U.matrixI.block(0,i,sizeU/2,1) = -Ureal.matrixR.block(sizeU/2,2*i,sizeU/2,1);
            V.matrixR.block(0,i,sizeV/2,1) = Vreal.matrixR.block(0,2*i,sizeV/2,1);
            V.matrixI.block(0,i,sizeV/2,1) = -Vreal.matrixR.block(sizeV/2,2*i,sizeV/2,1);
        }
        U.checkComplexity();
        V.checkComplexity();
    } else {
        JacobiSVD< Matrix<mpreal,Dynamic,Dynamic> > svd(matrixR, ComputeThinU | ComputeThinV);

        result.isComplex = false;
        result.matrixR = svd.singularValues();
        result = result.diagCreate(0);

        U.isComplex = false;
        U.matrixR = svd.matrixU();
        U.matrixI.resize(0,0);

        V.isComplex = false;
        V.matrixR = svd.matrixV();
        V.matrixI.resize(0,0);
    }

    return result;
}







/* -----------------------------
   | Some comparison operators |
   ----------------------------- */


Matrix <bool, Dynamic, Dynamic> GmpEigenMatrix::operator<(const GmpEigenMatrix& b) const
{
    Matrix <bool, Dynamic, Dynamic> result;

    // We want to support comparison between an nxm matrix and a 1x1 matrix
    if ((numel() != 1) && (b.numel() == 1)) {
        result = (matrixR.array() < b.matrixR(0,0));
    } else if ((numel() == 1) && (b.numel() != 1)) {
        result = (matrixR(0,0) < b.matrixR.array());
    } else {
        result = (matrixR.array() < b.matrixR.array());
    }

    return result;
}

Matrix <bool, Dynamic, Dynamic> GmpEigenMatrix::operator<=(const GmpEigenMatrix& b) const
{
    Matrix <bool, Dynamic, Dynamic> result;

    // We want to support comparison between an nxm matrix and a 1x1 matrix
    if ((numel() != 1) && (b.numel() == 1)) {
        result = (matrixR.array() <= b.matrixR(0,0));
    } else if ((numel() == 1) && (b.numel() != 1)) {
        result = (matrixR(0,0) <= b.matrixR.array());
    } else {
        result = (matrixR.array() <= b.matrixR.array());
    }

    return result;
}

Matrix <bool, Dynamic, Dynamic> GmpEigenMatrix::operator>(const GmpEigenMatrix& b) const
{
    Matrix <bool, Dynamic, Dynamic> result;

    // We want to support comparison between an nxm matrix and a 1x1 matrix
    if ((numel() != 1) && (b.numel() == 1)) {
        result = (matrixR.array() > b.matrixR(0,0));
    } else if ((numel() == 1) && (b.numel() != 1)) {
        result = (matrixR(0,0) > b.matrixR.array());
    } else {
        result = (matrixR.array() > b.matrixR.array());
    }

    return result;
}

Matrix <bool, Dynamic, Dynamic> GmpEigenMatrix::operator>=(const GmpEigenMatrix& b) const
{
    Matrix <bool, Dynamic, Dynamic> result;

    // We want to support comparison between an nxm matrix and a 1x1 matrix
    if ((numel() != 1) && (b.numel() == 1)) {
        result = (matrixR.array() >= b.matrixR(0,0));
    } else if ((numel() == 1) && (b.numel() != 1)) {
        result = (matrixR(0,0) >= b.matrixR.array());
    } else {
        result = (matrixR.array() >= b.matrixR.array());
    }

    return result;
}

/* Test of equality  c = (a == b) */
Matrix <bool, Dynamic, Dynamic> GmpEigenMatrix::eq(const GmpEigenMatrix& b) const
{
    Matrix <bool, Dynamic, Dynamic> result;

    // We want to support comparison between an nxm matrix and a 1x1 matrix
    if ((numel() != 1) && (b.numel() == 1)) {
        result = (matrixR.array() == b.matrixR(0,0));
        if (b.isComplex) {
            if (isComplex) {
                result = result.array() * (matrixI.array() == b.matrixI(0,0));
            } else {
                result = result * iszero(b.matrixI(0,0)); // This should always give a table of false --- '&&' operator is not supported by Eigen here
            }
        } else if (isComplex) {
            for (IndexType j(0); j < result.cols(); ++j)
                for (IndexType i(0); i < result.rows(); ++i)
                    result(i,j) = result(i,j) && iszero(matrixI(i,j));
        }
    } else if ((numel() == 1) && (b.numel() != 1)) {
        result = (matrixR(0,0) == b.matrixR.array());
        if (b.isComplex) {
            if (isComplex) {
                result = result.array() * (matrixI(0,0) == b.matrixI.array());
            } else {
                for (IndexType j(0); j < result.cols(); ++j)
                    for (IndexType i(0); i < result.rows(); ++i)
                        result(i,j) = result(i,j) && iszero(b.matrixI(i,j));
            }
        } else if (isComplex) {
            result = result * iszero(matrixI(0,0)); // This should always give a table of false
        }
    } else {
        result = (matrixR.array() == b.matrixR.array());
        if (b.isComplex) {
            if (isComplex) {
                result = result.array() * (matrixI.array() == b.matrixI.array());
            } else {
                for (IndexType j(0); j < result.cols(); ++j)
                    for (IndexType i(0); i < result.rows(); ++i)
                        result(i,j) = result(i,j) && iszero(b.matrixI(i,j));
            }
        } else if (isComplex) {
            for (IndexType j(0); j < result.cols(); ++j)
                for (IndexType i(0); i < result.rows(); ++i)
                    result(i,j) = result(i,j) && iszero(matrixI(i,j));
        }
    }

    return result;
}

/* Test of non-equality  c = (a != b) */
Matrix <bool, Dynamic, Dynamic> GmpEigenMatrix::ne(const GmpEigenMatrix& b) const
{
    Matrix <bool, Dynamic, Dynamic> result;

    // We want to support comparison between an nxm matrix and a 1x1 matrix
    if ((numel() != 1) && (b.numel() == 1)) {
        result = (matrixR.array() != b.matrixR(0,0));
        if (b.isComplex) {
            if (isComplex) {
                result = result.array() || (matrixI.array() != b.matrixI(0,0));
            } else {
                //result = result.array() || (!iszero(b.matrixI(0,0))); // This should always give a table of true --- it doesn't work anymore on Eigen 3.3
                if (!iszero(b.matrixI(0,0)))
                    result = result.array() || (!result.array());
            }
        } else if (isComplex) {
            for (IndexType j(0); j < result.cols(); ++j)
                for (IndexType i(0); i < result.rows(); ++i)
                    result(i,j) = result(i,j) || (!iszero(matrixI(i,j)));
        }
    } else if ((numel() == 1) && (b.numel() != 1)) {
        result = (matrixR(0,0) != b.matrixR.array());
        if (b.isComplex) {
            if (isComplex) {
                result = result.array() || (matrixI(0,0) != b.matrixI.array());
            } else {
                for (IndexType j(0); j < result.cols(); ++j)
                    for (IndexType i(0); i < result.rows(); ++i)
                        result(i,j) = result(i,j) || (!iszero(b.matrixI(i,j)));
            }
        } else if (isComplex) {
            //result = result.array() || (!iszero(matrixI(0,0))); // This should always give a table of true --- it doesn't work anymore on Eigen 3.3
            if (!iszero(matrixI(0,0)))
                result = result.array() || (!result.array());
        }
    } else {
        result = (matrixR.array() != b.matrixR.array());
        if (b.isComplex) {
            if (isComplex) {
                result = result.array() || (matrixI.array() != b.matrixI.array());
            } else {
                for (IndexType j(0); j < result.cols(); ++j)
                    for (IndexType i(0); i < result.rows(); ++i)
                        result(i,j) = result(i,j) || (!iszero(b.matrixI(i,j)));
            }
        } else if (isComplex) {
            for (IndexType j(0); j < result.cols(); ++j)
                for (IndexType i(0); i < result.rows(); ++i)
                    result(i,j) = result(i,j) || (!iszero(matrixI(i,j)));
        }
    }

    return result;
}

/* Test whether the elements are nan  b = isnan(a) */
Matrix <bool, Dynamic, Dynamic> GmpEigenMatrix::isnan() const
{
    Matrix <bool, Dynamic, Dynamic> result(matrixR.rows(), matrixR.cols());

    for (IndexType j(0); j < matrixR.cols(); ++j)
        for (IndexType i(0); i < matrixR.rows(); ++i)
            if (isComplex)
                result(i,j) = mpfr::isnan(matrixR(i,j)) || mpfr::isnan(matrixI(i,j));
            else
                result(i,j) = mpfr::isnan(matrixR(i,j));

    return result;
}

/* Test whether the elements are +/-inf  b = isinf(a) */
Matrix <bool, Dynamic, Dynamic> GmpEigenMatrix::isinf() const
{
    Matrix <bool, Dynamic, Dynamic> result(matrixR.rows(), matrixR.cols());

    for (IndexType j(0); j < matrixR.cols(); ++j)
        for (IndexType i(0); i < matrixR.rows(); ++i)
            if (isComplex)
                result(i,j) = mpfr::isinf(matrixR(i,j)) || mpfr::isinf(matrixI(i,j));
            else
                result(i,j) = mpfr::isinf(matrixR(i,j));

    return result;
}

/* Test whether the numerical values match
   this function is called by c = isequal(a, b)
   at this stage, we already know that both tables have the same size, that both
   are either real or complex, and that all their precision match. */
bool GmpEigenMatrix::identicalValues(const GmpEigenMatrix& b) const
{
    for (IndexType j(0); j < matrixR.cols(); ++j)
        for (IndexType i(0); i < matrixR.rows(); ++i)
        {
            bool result(matrixR(i,j) == b.matrixR(i,j));
            if ((isComplex) && (b.isComplex))
                result = result && (matrixI(i,j) == b.matrixI(i,j));
            if (!result)
                return false;
        }

    return true;
}

/* Test whether the numerical values match
   this function is called by c = isequaln(a, b)
   This function makes no difference between NaN numbers.
   at this stage, we already know that both tables have the same size, that both
   are either real or complex, and that all their precision match. */
bool GmpEigenMatrix::identicalValuesNaNok(const GmpEigenMatrix& b) const
{
    for (IndexType j(0); j < matrixR.cols(); ++j)
        for (IndexType i(0); i < matrixR.rows(); ++i)
        {
            bool result((matrixR(i,j) == b.matrixR(i,j)) || (mpfr::isnan(matrixR(i,j)) && mpfr::isnan(b.matrixR(i,j))));
            if ((isComplex) && (b.isComplex))
                result = result && ((matrixI(i,j) == b.matrixI(i,j)) || (mpfr::isnan(matrixI(i,j)) && mpfr::isnan(b.matrixI(i,j))));
            if (!result)
                return false;
        }

    return true;
}


// symmetry tests
bool GmpEigenMatrix::issymmetric() const
{
    if (isComplex) {
        for (IndexType i(0); i < matrixR.rows(); ++i)
            for (IndexType j(i+1); j < matrixR.cols(); ++j)
                if ((matrixR(i,j) != matrixR(j,i)) || (matrixI(i,j) != matrixI(j,i)))
                    return false;
    } else {
        for (IndexType i(0); i < matrixR.rows(); ++i)
            for (IndexType j(i+1); j < matrixR.cols(); ++j)
                if (matrixR(i,j) != matrixR(j,i))
                    return false;
    }

    return true;
}

bool GmpEigenMatrix::ishermitian() const
{
    if (isComplex) {
        for (IndexType i(0); i < matrixR.rows(); ++i) {
            for (IndexType j(i); j < matrixR.cols(); ++j)
                if ((matrixR(i,j) != matrixR(j,i)) || (matrixI(i,j) != -matrixI(j,i)))
                    return false;
        }
    } else {
        for (IndexType i(0); i < matrixR.rows(); ++i)
            for (IndexType j(i+1); j < matrixR.cols(); ++j)
                if (matrixR(i,j) != matrixR(j,i))
                    return false;
    }

    return true;
}



// Note : Eigen's minCoeff function behaves differently than matlab's min
// function (e.g. if gives no result in presence of NaN instead of ignoring them,
// so we write our own implementation)

// column-wise minimum b = min(a)
GmpEigenMatrix GmpEigenMatrix::colMin(vector<IndexType>& indices) const
{
    GmpEigenMatrix result;

    // We initialize the index vector
    indices.clear();
    indices.resize(matrixR.cols(), 0);

    if (isComplex) {
        // We initialize the result matrix
        result.matrixR.resize(1,matrixR.cols());
        result.matrixI.resize(1,matrixI.cols());
        result.isComplex = true;

        // We find the minimum
        for (IndexType j(0); j < matrixR.cols(); ++j)
        {
            mpreal minValue;
            mpreal minAngle;
            minValue.setInf(+1);
            minAngle.setInf(+1);
            for (IndexType i(0); i < matrixR.rows(); ++i) {
                if ((!mpfr::isnan(matrixR(i,j))) && (!mpfr::isnan(matrixI(i,j)))) {
                    if (((mpfr::pow(matrixR(i,j),2) + mpfr::pow(matrixI(i,j),2)) < minValue) || (((mpfr::pow(matrixR(i,j),2) + mpfr::pow(matrixI(i,j),2)) == minValue) && (atan2(matrixI(i,j), matrixR(i,j)) < minAngle))) {
                        indices[j] = i;
                        minValue = mpfr::pow(matrixR(i,j),2) + mpfr::pow(matrixI(i,j),2);
                        minAngle = atan2(matrixI(i,j), matrixR(i,j));
                    }
                }
            }
            result.matrixR(j) = matrixR(indices[j],j);
            result.matrixI(j) = matrixI(indices[j],j);
        }

        result.checkComplexity();
    } else {
        // We initialize the result matrix
        result.matrixR.resize(1,matrixR.cols());
        result.isComplex = false;

        // We find the minimum
        for (IndexType j(0); j < matrixR.cols(); ++j)
        {
            mpreal minValue;
            minValue.setInf(+1);
            for (IndexType i(0); i < matrixR.rows(); ++i) {
                if ((!mpfr::isnan(matrixR(i,j))) && (matrixR(i,j) < minValue)) {
                    indices[j] = i;
                    minValue = matrixR(i,j);
                }
            }
            result.matrixR(j) = matrixR(indices[j],j);
        }
    }

    return result;
}

// column-wise minimum b = min(a)
GmpEigenMatrix& GmpEigenMatrix::colMin_new(vector<IndexType>& indices) const
{
    GmpEigenMatrix& result(*(new GmpEigenMatrix));

    // We initialize the index vector
    indices.clear();
    indices.resize(matrixR.cols(), 0);

    if (isComplex) {
        // We initialize the result matrix
        result.matrixR.resize(1,matrixR.cols());
        result.matrixI.resize(1,matrixI.cols());
        result.isComplex = true;

        // We find the minimum
        for (IndexType j(0); j < matrixR.cols(); ++j)
        {
            mpreal minValue;
            mpreal minAngle;
            minValue.setInf(+1);
            minAngle.setInf(+1);
            for (IndexType i(0); i < matrixR.rows(); ++i) {
                if ((!mpfr::isnan(matrixR(i,j))) && (!mpfr::isnan(matrixI(i,j)))) {
                    if (((mpfr::pow(matrixR(i,j),2) + mpfr::pow(matrixI(i,j),2)) < minValue) || (((mpfr::pow(matrixR(i,j),2) + mpfr::pow(matrixI(i,j),2)) == minValue) && (atan2(matrixI(i,j), matrixR(i,j)) < minAngle))) {
                        indices[j] = i;
                        minValue = mpfr::pow(matrixR(i,j),2) + mpfr::pow(matrixI(i,j),2);
                        minAngle = atan2(matrixI(i,j), matrixR(i,j));
                    }
                }
            }
            result.matrixR(j) = matrixR(indices[j],j);
            result.matrixI(j) = matrixI(indices[j],j);
        }

        result.checkComplexity();
    } else {
        // We initialize the result matrix
        result.matrixR.resize(1,matrixR.cols());
        result.isComplex = false;

        // We find the minimum
        for (IndexType j(0); j < matrixR.cols(); ++j)
        {
            mpreal minValue;
            minValue.setInf(+1);
            for (IndexType i(0); i < matrixR.rows(); ++i) {
                if ((!mpfr::isnan(matrixR(i,j))) && (matrixR(i,j) < minValue)) {
                    indices[j] = i;
                    minValue = matrixR(i,j);
                }
            }
            result.matrixR(j) = matrixR(indices[j],j);
        }
    }

    return result;
}

// line-wise minimum b = min(a,[],2)
GmpEigenMatrix GmpEigenMatrix::rowMin(vector<IndexType>& indices) const
{
    GmpEigenMatrix result;

    // We initialize the index vector
    indices.clear();
    indices.resize(matrixR.rows(), 0);

    if (isComplex) {
        // We initialize the result matrix
        result.matrixR.resize(matrixR.rows(),1);
        result.matrixI.resize(matrixI.rows(),1);
        result.isComplex = true;

        // We find the minimum
        for (IndexType i(0); i < matrixR.rows(); ++i)
        {
            mpreal minValue;
            mpreal minAngle;
            minValue.setInf(+1);
            minAngle.setInf(+1);
            for (IndexType j(0); j < matrixR.cols(); ++j) {
                if ((!mpfr::isnan(matrixR(i,j))) && (!mpfr::isnan(matrixI(i,j)))) {
                    if (((mpfr::pow(matrixR(i,j),2) + mpfr::pow(matrixI(i,j),2)) < minValue) || (((mpfr::pow(matrixR(i,j),2) + mpfr::pow(matrixI(i,j),2)) == minValue) && (atan2(matrixI(i,j), matrixR(i,j)) < minAngle))) {
                        indices[i] = j;
                        minValue = mpfr::pow(matrixR(i,j),2) + mpfr::pow(matrixI(i,j),2);
                        minAngle = atan2(matrixI(i,j), matrixR(i,j));
                    }
                }
            }
            result.matrixR(i) = matrixR(i,indices[i]);
            result.matrixI(i) = matrixI(i,indices[i]);
        }

        result.checkComplexity();
    } else {
        // We initialize the result matrix
        result.matrixR.resize(matrixR.rows(),1);
        result.isComplex = false;

        // We find the minimum
        for (IndexType i(0); i < matrixR.rows(); ++i)
        {
            mpreal minValue;
            minValue.setInf(+1);
            for (IndexType j(0); j < matrixR.cols(); ++j) {
                if ((!mpfr::isnan(matrixR(i,j))) && (matrixR(i,j) < minValue)) {
                    indices[i] = j;
                    minValue = matrixR(i,j);
                }
            }
            result.matrixR(i) = matrixR(i,indices[i]);
        }
    }

    return result;
}

// line-wise minimum b = min(a,[],2)
GmpEigenMatrix& GmpEigenMatrix::rowMin_new(vector<IndexType>& indices) const
{
    GmpEigenMatrix& result(*(new GmpEigenMatrix));

    // We initialize the index vector
    indices.clear();
    indices.resize(matrixR.rows(), 0);

    if (isComplex) {
        // We initialize the result matrix
        result.matrixR.resize(matrixR.rows(),1);
        result.matrixI.resize(matrixI.rows(),1);
        result.isComplex = true;

        // We find the minimum
        for (IndexType i(0); i < matrixR.rows(); ++i)
        {
            mpreal minValue;
            mpreal minAngle;
            minValue.setInf(+1);
            minAngle.setInf(+1);
            for (IndexType j(0); j < matrixR.cols(); ++j) {
                if ((!mpfr::isnan(matrixR(i,j))) && (!mpfr::isnan(matrixI(i,j)))) {
                    if (((mpfr::pow(matrixR(i,j),2) + mpfr::pow(matrixI(i,j),2)) < minValue) || (((mpfr::pow(matrixR(i,j),2) + mpfr::pow(matrixI(i,j),2)) == minValue) && (atan2(matrixI(i,j), matrixR(i,j)) < minAngle))) {
                        indices[i] = j;
                        minValue = mpfr::pow(matrixR(i,j),2) + mpfr::pow(matrixI(i,j),2);
                        minAngle = atan2(matrixI(i,j), matrixR(i,j));
                    }
                }
            }
            result.matrixR(i) = matrixR(i,indices[i]);
            result.matrixI(i) = matrixI(i,indices[i]);
        }

        result.checkComplexity();
    } else {
        // We initialize the result matrix
        result.matrixR.resize(matrixR.rows(),1);
        result.isComplex = false;

        // We find the minimum
        for (IndexType i(0); i < matrixR.rows(); ++i)
        {
            mpreal minValue;
            minValue.setInf(+1);
            for (IndexType j(0); j < matrixR.cols(); ++j) {
                if ((!mpfr::isnan(matrixR(i,j))) && (matrixR(i,j) < minValue)) {
                    indices[i] = j;
                    minValue = matrixR(i,j);
                }
            }
            result.matrixR(i) = matrixR(i,indices[i]);
        }
    }

    return result;
}

// element-wise minimum c = min(a, b)
GmpEigenMatrix GmpEigenMatrix::ewMin(const GmpEigenMatrix& b) const
{
    GmpEigenMatrix result;

    if ((numel() != 1) && (b.numel() == 1)) {
        if ((b.isComplex) || (isComplex)) {
            result.matrixR.resize(matrixR.rows(),matrixR.cols());
            result.matrixI.resize(matrixR.rows(),matrixR.cols());
            // We find the minimums
            GmpEigenMatrix aAbs(abs()), aAngle(angle()), bAbs(b.abs()), bAngle(b.angle());
            for (IndexType i(0); i < matrixR.rows(); ++i) {
                for (IndexType j(0); j < matrixR.cols(); ++j) {
                    if ((mpfr::isnan(aAbs.matrixR(i,j))) && (mpfr::isnan(bAbs.matrixR(0,0)))) {
                        result.matrixR(i,j) = aAbs.matrixR(i,j);
                        result.matrixI(i,j) = 0;
                    } else if ((mpfr::isnan(bAbs.matrixR(0,0))) || (aAbs.matrixR(i,j) < bAbs.matrixR(0,0)) || ((aAbs.matrixR(i,j) == bAbs.matrixR(0,0)) && (aAngle.matrixR(i,j) <= bAngle.matrixR(0,0)))) {
                        result.matrixR(i,j) = matrixR(i,j);
                        if (isComplex)
                            result.matrixI(i,j) = matrixI(i,j);
                        else
                            result.matrixI(i,j) = 0;
                    } else {
                        result.matrixR(i,j) = b.matrixR(0,0);
                        if (b.isComplex)
                            result.matrixI(i,j) = b.matrixI(0,0);
                        else
                            result.matrixI(i,j) = 0;
                    }
                }
            }
            result.checkComplexity();
        } else {
            result.matrixR.resize(matrixR.rows(),matrixR.cols());
            for (IndexType i(0); i < matrixR.rows(); ++i) {
                for (IndexType j(0); j < matrixR.cols(); ++j) {
                    if ((mpfr::isnan(matrixR(i,j))) && (mpfr::isnan(b.matrixR(0,0)))) {
                        result.matrixR(i,j) = matrixR(i,j);
                    } else if ((mpfr::isnan(b.matrixR(0,0))) || (matrixR(i,j) <= b.matrixR(0,0))) {
                        result.matrixR(i,j) = matrixR(i,j);
                    } else {
                        result.matrixR(i,j) = b.matrixR(0,0);
                    }
                }
            }
            result.isComplex = false;
        }
    } else if ((numel() == 1) && (b.numel() != 1)) {
        if ((b.isComplex) || (isComplex)) {
            result.matrixR.resize(b.matrixR.rows(),b.matrixR.cols());
            result.matrixI.resize(b.matrixR.rows(),b.matrixR.cols());
            // We find the minimums
            GmpEigenMatrix aAbs(abs()), aAngle(angle()), bAbs(b.abs()), bAngle(b.angle());
            for (IndexType i(0); i < b.matrixR.rows(); ++i) {
                for (IndexType j(0); j < b.matrixR.cols(); ++j) {
                    if ((mpfr::isnan(aAbs.matrixR(0,0))) && (mpfr::isnan(bAbs.matrixR(i,j)))) {
                        result.matrixR(i,j) = aAbs.matrixR(0,0);
                        result.matrixI(i,j) = 0;
                    } else if ((mpfr::isnan(bAbs.matrixR(i,j))) || (aAbs.matrixR(0,0) < bAbs.matrixR(i,j)) || ((aAbs.matrixR(0,0) == bAbs.matrixR(i,j)) && (aAngle.matrixR(0,0) <= bAngle.matrixR(i,j)))) {
                        result.matrixR(i,j) = matrixR(0,0);
                        if (isComplex)
                            result.matrixI(i,j) = matrixI(0,0);
                        else
                            result.matrixI(i,j) = 0;
                    } else {
                        result.matrixR(i,j) = b.matrixR(i,j);
                        if (b.isComplex)
                            result.matrixI(i,j) = b.matrixI(i,j);
                        else
                            result.matrixI(i,j) = 0;
                    }
                }
            }
            result.checkComplexity();
        } else {
            result.matrixR.resize(b.matrixR.rows(),b.matrixR.cols());
            for (IndexType i(0); i < b.matrixR.rows(); ++i) {
                for (IndexType j(0); j < b.matrixR.cols(); ++j) {
                    if ((mpfr::isnan(matrixR(0,0))) && (mpfr::isnan(b.matrixR(i,j)))) {
                        result.matrixR(i,j) = matrixR(0,0);
                    } else if ((mpfr::isnan(b.matrixR(i,j))) || (matrixR(0,0) <= b.matrixR(i,j))) {
                        result.matrixR(i,j) = matrixR(0,0);
                    } else {
                        result.matrixR(i,j) = b.matrixR(i,j);
                    }
                }
            }
            result.isComplex = false;
        }
    } else {
        if ((b.isComplex) || (isComplex)) {
            result.matrixR.resize(matrixR.rows(),matrixR.cols());
            result.matrixI.resize(matrixR.rows(),matrixR.cols());
            // We find the minimums
            GmpEigenMatrix aAbs(abs()), aAngle(angle()), bAbs(b.abs()), bAngle(b.angle());
            for (IndexType i(0); i < matrixR.rows(); ++i) {
                for (IndexType j(0); j < matrixR.cols(); ++j) {
                    if ((mpfr::isnan(aAbs.matrixR(i,j))) && (mpfr::isnan(bAbs.matrixR(i,j)))) {
                        result.matrixR(i,j) = aAbs.matrixR(i,j);
                        result.matrixI(i,j) = 0;
                    } else if ((mpfr::isnan(bAbs.matrixR(i,j))) || (aAbs.matrixR(i,j) < bAbs.matrixR(i,j)) || ((aAbs.matrixR(i,j) == bAbs.matrixR(i,j)) && (aAngle.matrixR(i,j) <= bAngle.matrixR(i,j)))) {
                        result.matrixR(i,j) = matrixR(i,j);
                        if (isComplex)
                            result.matrixI(i,j) = matrixI(i,j);
                        else
                            result.matrixI(i,j) = 0;
                    } else {
                        result.matrixR(i,j) = b.matrixR(i,j);
                        if (b.isComplex)
                            result.matrixI(i,j) = b.matrixI(i,j);
                        else
                            result.matrixI(i,j) = 0;
                    }
                }
            }
            result.checkComplexity();
        } else {
            result.matrixR.resize(matrixR.rows(),matrixR.cols());
            for (IndexType i(0); i < matrixR.rows(); ++i) {
                for (IndexType j(0); j < matrixR.cols(); ++j) {
                    if ((mpfr::isnan(matrixR(i,j))) && (mpfr::isnan(b.matrixR(i,j)))) {
                        result.matrixR(i,j) = matrixR(i,j);
                    } else if ((mpfr::isnan(b.matrixR(i,j))) || (matrixR(i,j) <= b.matrixR(i,j))) {
                        result.matrixR(i,j) = matrixR(i,j);
                    } else {
                        result.matrixR(i,j) = b.matrixR(i,j);
                    }
                }
            }
            result.isComplex = false;
        }
    }

    return result;
}

// element-wise minimum c = min(a, b)
GmpEigenMatrix& GmpEigenMatrix::ewMin_new(const GmpEigenMatrix& b) const
{
    GmpEigenMatrix& result(*(new GmpEigenMatrix));

    if ((numel() != 1) && (b.numel() == 1)) {
        if ((b.isComplex) || (isComplex)) {
            result.matrixR.resize(matrixR.rows(),matrixR.cols());
            result.matrixI.resize(matrixR.rows(),matrixR.cols());
            // We find the minimums
            GmpEigenMatrix aAbs(abs()), aAngle(angle()), bAbs(b.abs()), bAngle(b.angle());
            for (IndexType i(0); i < matrixR.rows(); ++i) {
                for (IndexType j(0); j < matrixR.cols(); ++j) {
                    if ((mpfr::isnan(aAbs.matrixR(i,j))) && (mpfr::isnan(bAbs.matrixR(0,0)))) {
                        result.matrixR(i,j) = aAbs.matrixR(i,j);
                        result.matrixI(i,j) = 0;
                    } else if ((mpfr::isnan(bAbs.matrixR(0,0))) || (aAbs.matrixR(i,j) < bAbs.matrixR(0,0)) || ((aAbs.matrixR(i,j) == bAbs.matrixR(0,0)) && (aAngle.matrixR(i,j) <= bAngle.matrixR(0,0)))) {
                        result.matrixR(i,j) = matrixR(i,j);
                        if (isComplex)
                            result.matrixI(i,j) = matrixI(i,j);
                        else
                            result.matrixI(i,j) = 0;
                    } else {
                        result.matrixR(i,j) = b.matrixR(0,0);
                        if (b.isComplex)
                            result.matrixI(i,j) = b.matrixI(0,0);
                        else
                            result.matrixI(i,j) = 0;
                    }
                }
            }
            result.checkComplexity();
        } else {
            result.matrixR.resize(matrixR.rows(),matrixR.cols());
            for (IndexType i(0); i < matrixR.rows(); ++i) {
                for (IndexType j(0); j < matrixR.cols(); ++j) {
                    if ((mpfr::isnan(matrixR(i,j))) && (mpfr::isnan(b.matrixR(0,0)))) {
                        result.matrixR(i,j) = matrixR(i,j);
                    } else if ((mpfr::isnan(b.matrixR(0,0))) || (matrixR(i,j) <= b.matrixR(0,0))) {
                        result.matrixR(i,j) = matrixR(i,j);
                    } else {
                        result.matrixR(i,j) = b.matrixR(0,0);
                    }
                }
            }
            result.isComplex = false;
        }
    } else if ((numel() == 1) && (b.numel() != 1)) {
        if ((b.isComplex) || (isComplex)) {
            result.matrixR.resize(b.matrixR.rows(),b.matrixR.cols());
            result.matrixI.resize(b.matrixR.rows(),b.matrixR.cols());
            // We find the minimums
            GmpEigenMatrix aAbs(abs()), aAngle(angle()), bAbs(b.abs()), bAngle(b.angle());
            for (IndexType i(0); i < b.matrixR.rows(); ++i) {
                for (IndexType j(0); j < b.matrixR.cols(); ++j) {
                    if ((mpfr::isnan(aAbs.matrixR(0,0))) && (mpfr::isnan(bAbs.matrixR(i,j)))) {
                        result.matrixR(i,j) = aAbs.matrixR(0,0);
                        result.matrixI(i,j) = 0;
                    } else if ((mpfr::isnan(bAbs.matrixR(i,j))) || (aAbs.matrixR(0,0) < bAbs.matrixR(i,j)) || ((aAbs.matrixR(0,0) == bAbs.matrixR(i,j)) && (aAngle.matrixR(0,0) <= bAngle.matrixR(i,j)))) {
                        result.matrixR(i,j) = matrixR(0,0);
                        if (isComplex)
                            result.matrixI(i,j) = matrixI(0,0);
                        else
                            result.matrixI(i,j) = 0;
                    } else {
                        result.matrixR(i,j) = b.matrixR(i,j);
                        if (b.isComplex)
                            result.matrixI(i,j) = b.matrixI(i,j);
                        else
                            result.matrixI(i,j) = 0;
                    }
                }
            }
            result.checkComplexity();
        } else {
            result.matrixR.resize(b.matrixR.rows(),b.matrixR.cols());
            for (IndexType i(0); i < b.matrixR.rows(); ++i) {
                for (IndexType j(0); j < b.matrixR.cols(); ++j) {
                    if ((mpfr::isnan(matrixR(0,0))) && (mpfr::isnan(b.matrixR(i,j)))) {
                        result.matrixR(i,j) = matrixR(0,0);
                    } else if ((mpfr::isnan(b.matrixR(i,j))) || (matrixR(0,0) <= b.matrixR(i,j))) {
                        result.matrixR(i,j) = matrixR(0,0);
                    } else {
                        result.matrixR(i,j) = b.matrixR(i,j);
                    }
                }
            }
            result.isComplex = false;
        }
    } else {
        if ((b.isComplex) || (isComplex)) {
            result.matrixR.resize(matrixR.rows(),matrixR.cols());
            result.matrixI.resize(matrixR.rows(),matrixR.cols());
            // We find the minimums
            GmpEigenMatrix aAbs(abs()), aAngle(angle()), bAbs(b.abs()), bAngle(b.angle());
            for (IndexType i(0); i < matrixR.rows(); ++i) {
                for (IndexType j(0); j < matrixR.cols(); ++j) {
                    if ((mpfr::isnan(aAbs.matrixR(i,j))) && (mpfr::isnan(bAbs.matrixR(i,j)))) {
                        result.matrixR(i,j) = aAbs.matrixR(i,j);
                        result.matrixI(i,j) = 0;
                    } else if ((mpfr::isnan(bAbs.matrixR(i,j))) || (aAbs.matrixR(i,j) < bAbs.matrixR(i,j)) || ((aAbs.matrixR(i,j) == bAbs.matrixR(i,j)) && (aAngle.matrixR(i,j) <= bAngle.matrixR(i,j)))) {
                        result.matrixR(i,j) = matrixR(i,j);
                        if (isComplex)
                            result.matrixI(i,j) = matrixI(i,j);
                        else
                            result.matrixI(i,j) = 0;
                    } else {
                        result.matrixR(i,j) = b.matrixR(i,j);
                        if (b.isComplex)
                            result.matrixI(i,j) = b.matrixI(i,j);
                        else
                            result.matrixI(i,j) = 0;
                    }
                }
            }
            result.checkComplexity();
        } else {
            result.matrixR.resize(matrixR.rows(),matrixR.cols());
            for (IndexType i(0); i < matrixR.rows(); ++i) {
                for (IndexType j(0); j < matrixR.cols(); ++j) {
                    if ((mpfr::isnan(matrixR(i,j))) && (mpfr::isnan(b.matrixR(i,j)))) {
                        result.matrixR(i,j) = matrixR(i,j);
                    } else if ((mpfr::isnan(b.matrixR(i,j))) || (matrixR(i,j) <= b.matrixR(i,j))) {
                        result.matrixR(i,j) = matrixR(i,j);
                    } else {
                        result.matrixR(i,j) = b.matrixR(i,j);
                    }
                }
            }
            result.isComplex = false;
        }
    }

    return result;
}

// column-wise maximum b = max(a)
GmpEigenMatrix GmpEigenMatrix::colMax(vector<IndexType>& indices) const
{
    GmpEigenMatrix result;

    // We initialize the index vector
    indices.clear();
    indices.resize(matrixR.cols(), 0);

    if (isComplex) {
        // We initialize the result matrix
        result.matrixR.resize(1,matrixR.cols());
        result.matrixI.resize(1,matrixI.cols());
        result.isComplex = true;

        // We find the minimum
        for (IndexType j(0); j < matrixR.cols(); ++j)
        {
            mpreal maxValue;
            mpreal maxAngle;
            maxValue.setInf(-1);
            maxAngle.setInf(-1);
            for (IndexType i(0); i < matrixR.rows(); ++i) {
                if ((!mpfr::isnan(matrixR(i,j))) && (!mpfr::isnan(matrixI(i,j)))) {
                    if (((mpfr::pow(matrixR(i,j),2) + mpfr::pow(matrixI(i,j),2)) > maxValue) || (((mpfr::pow(matrixR(i,j),2) + mpfr::pow(matrixI(i,j),2)) == maxValue) && (atan2(matrixI(i,j), matrixR(i,j)) > maxAngle))) {
                        indices[j] = i;
                        maxValue = mpfr::pow(matrixR(i,j),2) + mpfr::pow(matrixI(i,j),2);
                        maxAngle = atan2(matrixI(i,j), matrixR(i,j));
                    }
                }
            }
            result.matrixR(j) = matrixR(indices[j],j);
            result.matrixI(j) = matrixI(indices[j],j);
        }

        result.checkComplexity();
    } else {
        // We initialize the result matrix
        result.matrixR.resize(1,matrixR.cols());
        result.isComplex = false;

        // We find the minimum
        for (IndexType j(0); j < matrixR.cols(); ++j)
        {
            mpreal maxValue;
            maxValue.setInf(-1);
            for (IndexType i(0); i < matrixR.rows(); ++i) {
                if ((!mpfr::isnan(matrixR(i,j))) && (matrixR(i,j) > maxValue)) {
                    indices[j] = i;
                    maxValue = matrixR(i,j);
                }
            }
            result.matrixR(j) = matrixR(indices[j],j);
        }
    }

    return result;
}

// column-wise maximum b = max(a)
GmpEigenMatrix& GmpEigenMatrix::colMax_new(vector<IndexType>& indices) const
{
    GmpEigenMatrix& result(*(new GmpEigenMatrix));

    // We initialize the index vector
    indices.clear();
    indices.resize(matrixR.cols(), 0);

    if (isComplex) {
        // We initialize the result matrix
        result.matrixR.resize(1,matrixR.cols());
        result.matrixI.resize(1,matrixI.cols());
        result.isComplex = true;

        // We find the minimum
        for (IndexType j(0); j < matrixR.cols(); ++j)
        {
            mpreal maxValue;
            mpreal maxAngle;
            maxValue.setInf(-1);
            maxAngle.setInf(-1);
            for (IndexType i(0); i < matrixR.rows(); ++i) {
                if ((!mpfr::isnan(matrixR(i,j))) && (!mpfr::isnan(matrixI(i,j)))) {
                    if (((mpfr::pow(matrixR(i,j),2) + mpfr::pow(matrixI(i,j),2)) > maxValue) || (((mpfr::pow(matrixR(i,j),2) + mpfr::pow(matrixI(i,j),2)) == maxValue) && (atan2(matrixI(i,j), matrixR(i,j)) > maxAngle))) {
                        indices[j] = i;
                        maxValue = mpfr::pow(matrixR(i,j),2) + mpfr::pow(matrixI(i,j),2);
                        maxAngle = atan2(matrixI(i,j), matrixR(i,j));
                    }
                }
            }
            result.matrixR(j) = matrixR(indices[j],j);
            result.matrixI(j) = matrixI(indices[j],j);
        }

        result.checkComplexity();
    } else {
        // We initialize the result matrix
        result.matrixR.resize(1,matrixR.cols());
        result.isComplex = false;

        // We find the minimum
        for (IndexType j(0); j < matrixR.cols(); ++j)
        {
            mpreal maxValue;
            maxValue.setInf(-1);
            for (IndexType i(0); i < matrixR.rows(); ++i) {
                if ((!mpfr::isnan(matrixR(i,j))) && (matrixR(i,j) > maxValue)) {
                    indices[j] = i;
                    maxValue = matrixR(i,j);
                }
            }
            result.matrixR(j) = matrixR(indices[j],j);
        }
    }

    return result;
}

// line-wise maximum b = max(a,[],2)
GmpEigenMatrix GmpEigenMatrix::rowMax(vector<IndexType>& indices) const
{
    GmpEigenMatrix result;

    // We initialize the index vector
    indices.clear();
    indices.resize(matrixR.rows(), 0);

    if (isComplex) {
        // We initialize the result matrix
        result.matrixR.resize(matrixR.rows(),1);
        result.matrixI.resize(matrixI.rows(),1);
        result.isComplex = true;

        // We find the minimum
        for (IndexType i(0); i < matrixR.rows(); ++i)
        {
            mpreal maxValue;
            mpreal maxAngle;
            maxValue.setInf(-1);
            maxAngle.setInf(-1);
            for (IndexType j(0); j < matrixR.cols(); ++j) {
                if ((!mpfr::isnan(matrixR(i,j))) && (!mpfr::isnan(matrixI(i,j)))) {
                    if (((mpfr::pow(matrixR(i,j),2) + mpfr::pow(matrixI(i,j),2)) > maxValue) || (((mpfr::pow(matrixR(i,j),2) + mpfr::pow(matrixI(i,j),2)) == maxValue) && (atan2(matrixI(i,j), matrixR(i,j)) > maxAngle))) {
                        indices[i] = j;
                        maxValue = mpfr::pow(matrixR(i,j),2) + mpfr::pow(matrixI(i,j),2);
                        maxAngle = atan2(matrixI(i,j), matrixR(i,j));
                    }
                }
            }
            result.matrixR(i) = matrixR(i,indices[i]);
            result.matrixI(i) = matrixI(i,indices[i]);
        }

        result.checkComplexity();
    } else {
        // We initialize the result matrix
        result.matrixR.resize(matrixR.rows(),1);
        result.isComplex = false;

        // We find the minimum
        for (IndexType i(0); i < matrixR.rows(); ++i)
        {
            mpreal maxValue;
            maxValue.setInf(-1);
            for (IndexType j(0); j < matrixR.cols(); ++j) {
                if ((!mpfr::isnan(matrixR(i,j))) && (matrixR(i,j) > maxValue)) {
                    indices[i] = j;
                    maxValue = matrixR(i,j);
                }
            }
            result.matrixR(i) = matrixR(i,indices[i]);
        }
    }

    return result;
}

// line-wise maximum b = max(a,[],2)
GmpEigenMatrix& GmpEigenMatrix::rowMax_new(vector<IndexType>& indices) const
{
    GmpEigenMatrix& result(*(new GmpEigenMatrix));

    // We initialize the index vector
    indices.clear();
    indices.resize(matrixR.rows(), 0);

    if (isComplex) {
        // We initialize the result matrix
        result.matrixR.resize(matrixR.rows(),1);
        result.matrixI.resize(matrixI.rows(),1);
        result.isComplex = true;

        // We find the minimum
        for (IndexType i(0); i < matrixR.rows(); ++i)
        {
            mpreal maxValue;
            mpreal maxAngle;
            maxValue.setInf(-1);
            maxAngle.setInf(-1);
            for (IndexType j(0); j < matrixR.cols(); ++j) {
                if ((!mpfr::isnan(matrixR(i,j))) && (!mpfr::isnan(matrixI(i,j)))) {
                    if (((mpfr::pow(matrixR(i,j),2) + mpfr::pow(matrixI(i,j),2)) > maxValue) || (((mpfr::pow(matrixR(i,j),2) + mpfr::pow(matrixI(i,j),2)) == maxValue) && (atan2(matrixI(i,j), matrixR(i,j)) > maxAngle))) {
                        indices[i] = j;
                        maxValue = mpfr::pow(matrixR(i,j),2) + mpfr::pow(matrixI(i,j),2);
                        maxAngle = atan2(matrixI(i,j), matrixR(i,j));
                    }
                }
            }
            result.matrixR(i) = matrixR(i,indices[i]);
            result.matrixI(i) = matrixI(i,indices[i]);
        }

        result.checkComplexity();
    } else {
        // We initialize the result matrix
        result.matrixR.resize(matrixR.rows(),1);
        result.isComplex = false;

        // We find the minimum
        for (IndexType i(0); i < matrixR.rows(); ++i)
        {
            mpreal maxValue;
            maxValue.setInf(-1);
            for (IndexType j(0); j < matrixR.cols(); ++j) {
                if ((!mpfr::isnan(matrixR(i,j))) && (matrixR(i,j) > maxValue)) {
                    indices[i] = j;
                    maxValue = matrixR(i,j);
                }
            }
            result.matrixR(i) = matrixR(i,indices[i]);
        }
    }

    return result;
}

// element-wise maximum c = max(a, b)
GmpEigenMatrix GmpEigenMatrix::ewMax(const GmpEigenMatrix& b) const
{
    GmpEigenMatrix result;

    if ((numel() != 1) && (b.numel() == 1)) {
        if ((b.isComplex) || (isComplex)) {
            result.matrixR.resize(matrixR.rows(),matrixR.cols());
            result.matrixI.resize(matrixR.rows(),matrixR.cols());
            // We find the minimums
            GmpEigenMatrix aAbs(abs()), aAngle(angle()), bAbs(b.abs()), bAngle(b.angle());
            for (IndexType i(0); i < matrixR.rows(); ++i) {
                for (IndexType j(0); j < matrixR.cols(); ++j) {
                    if ((mpfr::isnan(aAbs.matrixR(i,j))) && (mpfr::isnan(bAbs.matrixR(0,0)))) {
                        result.matrixR(i,j) = aAbs.matrixR(i,j);
                        result.matrixI(i,j) = 0;
                    } else if ((mpfr::isnan(bAbs.matrixR(0,0))) || (aAbs.matrixR(i,j) > bAbs.matrixR(0,0)) || ((aAbs.matrixR(i,j) == bAbs.matrixR(0,0)) && (aAngle.matrixR(i,j) >= bAngle.matrixR(0,0)))) {
                        result.matrixR(i,j) = matrixR(i,j);
                        if (isComplex)
                            result.matrixI(i,j) = matrixI(i,j);
                        else
                            result.matrixI(i,j) = 0;
                    } else {
                        result.matrixR(i,j) = b.matrixR(0,0);
                        if (b.isComplex)
                            result.matrixI(i,j) = b.matrixI(0,0);
                        else
                            result.matrixI(i,j) = 0;
                    }
                }
            }
            result.checkComplexity();
        } else {
            result.matrixR.resize(matrixR.rows(),matrixR.cols());
            for (IndexType i(0); i < matrixR.rows(); ++i) {
                for (IndexType j(0); j < matrixR.cols(); ++j) {
                    if ((mpfr::isnan(matrixR(i,j))) && (mpfr::isnan(b.matrixR(0,0)))) {
                        result.matrixR(i,j) = matrixR(i,j);
                    } else if ((mpfr::isnan(b.matrixR(0,0))) || (matrixR(i,j) >= b.matrixR(0,0))) {
                        result.matrixR(i,j) = matrixR(i,j);
                    } else {
                        result.matrixR(i,j) = b.matrixR(0,0);
                    }
                }
            }
            result.isComplex = false;
        }
    } else if ((numel() == 1) && (b.numel() != 1)) {
        if ((b.isComplex) || (isComplex)) {
            result.matrixR.resize(b.matrixR.rows(),b.matrixR.cols());
            result.matrixI.resize(b.matrixR.rows(),b.matrixR.cols());
            // We find the minimums
            GmpEigenMatrix aAbs(abs()), aAngle(angle()), bAbs(b.abs()), bAngle(b.angle());
            for (IndexType i(0); i < b.matrixR.rows(); ++i) {
                for (IndexType j(0); j < b.matrixR.cols(); ++j) {
                    if ((mpfr::isnan(aAbs.matrixR(0,0))) && (mpfr::isnan(bAbs.matrixR(i,j)))) {
                        result.matrixR(i,j) = aAbs.matrixR(0,0);
                        result.matrixI(i,j) = 0;
                    } else if ((mpfr::isnan(bAbs.matrixR(i,j))) || (aAbs.matrixR(0,0) > bAbs.matrixR(i,j)) || ((aAbs.matrixR(0,0) == bAbs.matrixR(i,j)) && (aAngle.matrixR(0,0) >= bAngle.matrixR(i,j)))) {
                        result.matrixR(i,j) = matrixR(0,0);
                        if (isComplex)
                            result.matrixI(i,j) = matrixI(0,0);
                        else
                            result.matrixI(i,j) = 0;
                    } else {
                        result.matrixR(i,j) = b.matrixR(i,j);
                        if (b.isComplex)
                            result.matrixI(i,j) = b.matrixI(i,j);
                        else
                            result.matrixI(i,j) = 0;
                    }
                }
            }
            result.checkComplexity();
        } else {
            result.matrixR.resize(b.matrixR.rows(),b.matrixR.cols());
            for (IndexType i(0); i < b.matrixR.rows(); ++i) {
                for (IndexType j(0); j < b.matrixR.cols(); ++j) {
                    if ((mpfr::isnan(matrixR(0,0))) && (mpfr::isnan(b.matrixR(i,j)))) {
                        result.matrixR(i,j) = matrixR(0,0);
                    } else if ((mpfr::isnan(b.matrixR(i,j))) || (matrixR(0,0) >= b.matrixR(i,j))) {
                        result.matrixR(i,j) = matrixR(0,0);
                    } else {
                        result.matrixR(i,j) = b.matrixR(i,j);
                    }
                }
            }
            result.isComplex = false;
        }
    } else {
        if ((b.isComplex) || (isComplex)) {
            result.matrixR.resize(matrixR.rows(),matrixR.cols());
            result.matrixI.resize(matrixR.rows(),matrixR.cols());
            // We find the minimums
            GmpEigenMatrix aAbs(abs()), aAngle(angle()), bAbs(b.abs()), bAngle(b.angle());
            for (IndexType i(0); i < matrixR.rows(); ++i) {
                for (IndexType j(0); j < matrixR.cols(); ++j) {
                    if ((mpfr::isnan(aAbs.matrixR(i,j))) && (mpfr::isnan(bAbs.matrixR(i,j)))) {
                        result.matrixR(i,j) = aAbs.matrixR(i,j);
                        result.matrixI(i,j) = 0;
                    } else if ((mpfr::isnan(bAbs.matrixR(i,j))) || (aAbs.matrixR(i,j) > bAbs.matrixR(i,j)) || ((aAbs.matrixR(i,j) == bAbs.matrixR(i,j)) && (aAngle.matrixR(i,j) >= bAngle.matrixR(i,j)))) {
                        result.matrixR(i,j) = matrixR(i,j);
                        if (isComplex)
                            result.matrixI(i,j) = matrixI(i,j);
                        else
                            result.matrixI(i,j) = 0;
                    } else {
                        result.matrixR(i,j) = b.matrixR(i,j);
                        if (b.isComplex)
                            result.matrixI(i,j) = b.matrixI(i,j);
                        else
                            result.matrixI(i,j) = 0;
                    }
                }
            }
            result.checkComplexity();
        } else {
            result.matrixR.resize(matrixR.rows(),matrixR.cols());
            for (IndexType i(0); i < matrixR.rows(); ++i) {
                for (IndexType j(0); j < matrixR.cols(); ++j) {
                    if ((mpfr::isnan(matrixR(i,j))) && (mpfr::isnan(b.matrixR(i,j)))) {
                        result.matrixR(i,j) = matrixR(i,j);
                    } else if ((mpfr::isnan(b.matrixR(i,j))) || (matrixR(i,j) >= b.matrixR(i,j))) {
                        result.matrixR(i,j) = matrixR(i,j);
                    } else {
                        result.matrixR(i,j) = b.matrixR(i,j);
                    }
                }
            }
            result.isComplex = false;
        }
    }

    return result;
}

// element-wise maximum c = max(a, b)
GmpEigenMatrix& GmpEigenMatrix::ewMax_new(const GmpEigenMatrix& b) const
{
    GmpEigenMatrix& result(*(new GmpEigenMatrix));

    if ((numel() != 1) && (b.numel() == 1)) {
        if ((b.isComplex) || (isComplex)) {
            result.matrixR.resize(matrixR.rows(),matrixR.cols());
            result.matrixI.resize(matrixR.rows(),matrixR.cols());
            // We find the minimums
            GmpEigenMatrix aAbs(abs()), aAngle(angle()), bAbs(b.abs()), bAngle(b.angle());
            for (IndexType i(0); i < matrixR.rows(); ++i) {
                for (IndexType j(0); j < matrixR.cols(); ++j) {
                    if ((mpfr::isnan(aAbs.matrixR(i,j))) && (mpfr::isnan(bAbs.matrixR(0,0)))) {
                        result.matrixR(i,j) = aAbs.matrixR(i,j);
                        result.matrixI(i,j) = 0;
                    } else if ((mpfr::isnan(bAbs.matrixR(0,0))) || (aAbs.matrixR(i,j) > bAbs.matrixR(0,0)) || ((aAbs.matrixR(i,j) == bAbs.matrixR(0,0)) && (aAngle.matrixR(i,j) >= bAngle.matrixR(0,0)))) {
                        result.matrixR(i,j) = matrixR(i,j);
                        if (isComplex)
                            result.matrixI(i,j) = matrixI(i,j);
                        else
                            result.matrixI(i,j) = 0;
                    } else {
                        result.matrixR(i,j) = b.matrixR(0,0);
                        if (b.isComplex)
                            result.matrixI(i,j) = b.matrixI(0,0);
                        else
                            result.matrixI(i,j) = 0;
                    }
                }
            }
            result.checkComplexity();
        } else {
            result.matrixR.resize(matrixR.rows(),matrixR.cols());
            for (IndexType i(0); i < matrixR.rows(); ++i) {
                for (IndexType j(0); j < matrixR.cols(); ++j) {
                    if ((mpfr::isnan(matrixR(i,j))) && (mpfr::isnan(b.matrixR(0,0)))) {
                        result.matrixR(i,j) = matrixR(i,j);
                    } else if ((mpfr::isnan(b.matrixR(0,0))) || (matrixR(i,j) >= b.matrixR(0,0))) {
                        result.matrixR(i,j) = matrixR(i,j);
                    } else {
                        result.matrixR(i,j) = b.matrixR(0,0);
                    }
                }
            }
            result.isComplex = false;
        }
    } else if ((numel() == 1) && (b.numel() != 1)) {
        if ((b.isComplex) || (isComplex)) {
            result.matrixR.resize(b.matrixR.rows(),b.matrixR.cols());
            result.matrixI.resize(b.matrixR.rows(),b.matrixR.cols());
            // We find the minimums
            GmpEigenMatrix aAbs(abs()), aAngle(angle()), bAbs(b.abs()), bAngle(b.angle());
            for (IndexType i(0); i < b.matrixR.rows(); ++i) {
                for (IndexType j(0); j < b.matrixR.cols(); ++j) {
                    if ((mpfr::isnan(aAbs.matrixR(0,0))) && (mpfr::isnan(bAbs.matrixR(i,j)))) {
                        result.matrixR(i,j) = aAbs.matrixR(0,0);
                        result.matrixI(i,j) = 0;
                    } else if ((mpfr::isnan(bAbs.matrixR(i,j))) || (aAbs.matrixR(0,0) > bAbs.matrixR(i,j)) || ((aAbs.matrixR(0,0) == bAbs.matrixR(i,j)) && (aAngle.matrixR(0,0) >= bAngle.matrixR(i,j)))) {
                        result.matrixR(i,j) = matrixR(0,0);
                        if (isComplex)
                            result.matrixI(i,j) = matrixI(0,0);
                        else
                            result.matrixI(i,j) = 0;
                    } else {
                        result.matrixR(i,j) = b.matrixR(i,j);
                        if (b.isComplex)
                            result.matrixI(i,j) = b.matrixI(i,j);
                        else
                            result.matrixI(i,j) = 0;
                    }
                }
            }
            result.checkComplexity();
        } else {
            result.matrixR.resize(b.matrixR.rows(),b.matrixR.cols());
            for (IndexType i(0); i < b.matrixR.rows(); ++i) {
                for (IndexType j(0); j < b.matrixR.cols(); ++j) {
                    if ((mpfr::isnan(matrixR(0,0))) && (mpfr::isnan(b.matrixR(i,j)))) {
                        result.matrixR(i,j) = matrixR(0,0);
                    } else if ((mpfr::isnan(b.matrixR(i,j))) || (matrixR(0,0) >= b.matrixR(i,j))) {
                        result.matrixR(i,j) = matrixR(0,0);
                    } else {
                        result.matrixR(i,j) = b.matrixR(i,j);
                    }
                }
            }
            result.isComplex = false;
        }
    } else {
        if ((b.isComplex) || (isComplex)) {
            result.matrixR.resize(matrixR.rows(),matrixR.cols());
            result.matrixI.resize(matrixR.rows(),matrixR.cols());
            // We find the minimums
            GmpEigenMatrix aAbs(abs()), aAngle(angle()), bAbs(b.abs()), bAngle(b.angle());
            for (IndexType i(0); i < matrixR.rows(); ++i) {
                for (IndexType j(0); j < matrixR.cols(); ++j) {
                    if ((mpfr::isnan(aAbs.matrixR(i,j))) && (mpfr::isnan(bAbs.matrixR(i,j)))) {
                        result.matrixR(i,j) = aAbs.matrixR(i,j);
                        result.matrixI(i,j) = 0;
                    } else if ((mpfr::isnan(bAbs.matrixR(i,j))) || (aAbs.matrixR(i,j) > bAbs.matrixR(i,j)) || ((aAbs.matrixR(i,j) == bAbs.matrixR(i,j)) && (aAngle.matrixR(i,j) >= bAngle.matrixR(i,j)))) {
                        result.matrixR(i,j) = matrixR(i,j);
                        if (isComplex)
                            result.matrixI(i,j) = matrixI(i,j);
                        else
                            result.matrixI(i,j) = 0;
                    } else {
                        result.matrixR(i,j) = b.matrixR(i,j);
                        if (b.isComplex)
                            result.matrixI(i,j) = b.matrixI(i,j);
                        else
                            result.matrixI(i,j) = 0;
                    }
                }
            }
            result.checkComplexity();
        } else {
            result.matrixR.resize(matrixR.rows(),matrixR.cols());
            for (IndexType i(0); i < matrixR.rows(); ++i) {
                for (IndexType j(0); j < matrixR.cols(); ++j) {
                    if ((mpfr::isnan(matrixR(i,j))) && (mpfr::isnan(b.matrixR(i,j)))) {
                        result.matrixR(i,j) = matrixR(i,j);
                    } else if ((mpfr::isnan(b.matrixR(i,j))) || (matrixR(i,j) >= b.matrixR(i,j))) {
                        result.matrixR(i,j) = matrixR(i,j);
                    } else {
                        result.matrixR(i,j) = b.matrixR(i,j);
                    }
                }
            }
            result.isComplex = false;
        }
    }

    return result;
}

// column-wise sum b = sum(a)
GmpEigenMatrix GmpEigenMatrix::colSum() const
{
    GmpEigenMatrix result;

    result.matrixR = matrixR.colwise().sum();
    if (isComplex)
        result.matrixI = matrixI.colwise().sum();

    result.checkComplexity();

    return result;
}

// column-wise sum b = sum(a)
GmpEigenMatrix& GmpEigenMatrix::colSum_new() const
{
    GmpEigenMatrix& result(*(new GmpEigenMatrix));

    result.matrixR = matrixR.colwise().sum();
    if (isComplex)
        result.matrixI = matrixI.colwise().sum();

    result.checkComplexity();

    return result;
}

// row-wise sum b = sum(a,2)
GmpEigenMatrix GmpEigenMatrix::rowSum() const
{
    GmpEigenMatrix result;

    result.matrixR = matrixR.rowwise().sum();
    if (isComplex)
        result.matrixI = matrixI.rowwise().sum();

    result.checkComplexity();

    return result;
}

// row-wise sum b = sum(a,2)
GmpEigenMatrix& GmpEigenMatrix::rowSum_new() const
{
    GmpEigenMatrix& result(*(new GmpEigenMatrix));

    result.matrixR = matrixR.rowwise().sum();
    if (isComplex)
        result.matrixI = matrixI.rowwise().sum();

    result.checkComplexity();

    return result;
}

// column-wise product b = prod(a)
GmpEigenMatrix GmpEigenMatrix::colProd() const
{
    GmpEigenMatrix result(block(0,0,1,matrixR.cols()));

    for (IndexType i(1); i < matrixR.rows(); ++i)
        result = result.times(block(i,0,1,matrixR.cols()));

    return result;
}

// column-wise product b = prod(a)
GmpEigenMatrix& GmpEigenMatrix::colProd_new() const
{
    GmpEigenMatrix& result(*(new GmpEigenMatrix(block(0,0,1,matrixR.cols()))));

    for (IndexType i(1); i < matrixR.rows(); ++i)
        result = result.times(block(i,0,1,matrixR.cols()));

    return result;
}

// row-wise product b = prod(a,2)
GmpEigenMatrix GmpEigenMatrix::rowProd() const
{
    GmpEigenMatrix result(block(0,0,matrixR.rows(),1));

    for (IndexType j(1); j < matrixR.cols(); ++j)
        result = result.times(block(0,j,matrixR.rows(),1));

    return result;
}

// row-wise product b = prod(a,2)
GmpEigenMatrix& GmpEigenMatrix::rowProd_new() const
{
    GmpEigenMatrix& result(*(new GmpEigenMatrix(block(0,0,matrixR.rows(),1))));

    for (IndexType j(1); j < matrixR.cols(); ++j)
        result = result.times(block(0,j,matrixR.rows(),1));

    return result;
}



// sorting elements along columns or rows
GmpEigenMatrix GmpEigenMatrix::sort(const int& dim, const int& type, vector < vector < IndexType > >& index) const
{
    GmpEigenMatrix result(*this);

    for (IndexType i(0); i < matrixR.rows(); ++i)
        index.push_back(vector < IndexType > (matrixR.cols(), 0));

    if (dim == 0) {
        for (IndexType j(0); j < matrixR.cols(); ++j) {
            vector < IndexType > indexMap(matrixR.rows(), 0);
            for (IndexType i(0); i < indexMap.size(); ++i)
                indexMap[i] = i;

            if (isComplex)
                std::sort(indexMap.begin(), indexMap.end(), compareThroughIndicesComplex< Matrix< mpreal, Dynamic, 1 >, Matrix< mpreal, Dynamic, 1 > >(matrixR.col(j), matrixI.col(j)));
            else
                std::sort(indexMap.begin(), indexMap.end(), compareThroughIndices< Matrix< mpreal, Dynamic, 1 > >(matrixR.col(j)));

            if (type == 1)
                std::reverse(indexMap.begin(), indexMap.end());

            for (IndexType i(0); i < matrixR.rows(); ++i) {
                result.matrixR(i,j) = matrixR(indexMap[i],j);
                if (isComplex)
                    result.matrixI(i,j) = matrixI(indexMap[i],j);
                index[i][j] = indexMap[i];
            }
        }
    } else {
        for (IndexType i(0); i < matrixR.rows(); ++i) {
            vector < IndexType > indexMap(matrixR.cols(), 0);
            for (IndexType j(0); j < indexMap.size(); ++j)
                indexMap[j] = j;

            if (isComplex)
                std::sort(indexMap.begin(), indexMap.end(), compareThroughIndicesComplex< Matrix< mpreal, Dynamic, 1 >, Matrix< mpreal, Dynamic, 1 > >(matrixR.row(i), matrixI.row(i)));
            else
                std::sort(indexMap.begin(), indexMap.end(), compareThroughIndices< Matrix< mpreal, Dynamic, 1 > >(matrixR.row(i)));

            if (type == 1)
                std::reverse(indexMap.begin(), indexMap.end());

            for (IndexType j(0); j < matrixR.cols(); ++j) {
                result.matrixR(i,j) = matrixR(i,indexMap[j]);
                if (isComplex)
                    result.matrixI(i,j) = matrixI(i,indexMap[j]);
                index[i][j] = indexMap[j];
            }
        }
    }

    return result;
}

// sorting elements along columns or rows
GmpEigenMatrix& GmpEigenMatrix::sort_new(const int& dim, const int& type, vector < vector < IndexType > >& index) const
{
    GmpEigenMatrix& result(*(new GmpEigenMatrix(*this)));

    for (IndexType i(0); i < matrixR.rows(); ++i)
        index.push_back(vector < IndexType > (matrixR.cols(), 0));

    if (dim == 0) {
        for (IndexType j(0); j < matrixR.cols(); ++j) {
            vector < IndexType > indexMap(matrixR.rows(), 0);
            for (IndexType i(0); i < indexMap.size(); ++i)
                indexMap[i] = i;

            if (isComplex)
                std::sort(indexMap.begin(), indexMap.end(), compareThroughIndicesComplex< Matrix< mpreal, Dynamic, 1 >, Matrix< mpreal, Dynamic, 1 > >(matrixR.col(j), matrixI.col(j)));
            else
                std::sort(indexMap.begin(), indexMap.end(), compareThroughIndices< Matrix< mpreal, Dynamic, 1 > >(matrixR.col(j)));

            if (type == 1)
                std::reverse(indexMap.begin(), indexMap.end());

            for (IndexType i(0); i < matrixR.rows(); ++i) {
                result.matrixR(i,j) = matrixR(indexMap[i],j);
                if (isComplex)
                    result.matrixI(i,j) = matrixI(indexMap[i],j);
                index[i][j] = indexMap[i];
            }
        }
    } else {
        for (IndexType i(0); i < matrixR.rows(); ++i) {
            vector < IndexType > indexMap(matrixR.cols(), 0);
            for (IndexType j(0); j < indexMap.size(); ++j)
                indexMap[j] = j;

            if (isComplex)
                std::sort(indexMap.begin(), indexMap.end(), compareThroughIndicesComplex< Matrix< mpreal, Dynamic, 1 >, Matrix< mpreal, Dynamic, 1 > >(matrixR.row(i), matrixI.row(i)));
            else
                std::sort(indexMap.begin(), indexMap.end(), compareThroughIndices< Matrix< mpreal, Dynamic, 1 > >(matrixR.row(i)));

            if (type == 1)
                std::reverse(indexMap.begin(), indexMap.end());

            for (IndexType j(0); j < matrixR.cols(); ++j) {
                result.matrixR(i,j) = matrixR(i,indexMap[j]);
                if (isComplex)
                    result.matrixI(i,j) = matrixI(i,indexMap[j]);
                index[i][j] = indexMap[j];
            }
        }
    }

    return result;
}


/* This function returns the recipe to sort the rows of a matrix */
vector < IndexType > GmpEigenMatrix::sortrowsc(const std::vector < int > ascending) const
{
    vector < IndexType > indexMap(matrixR.rows(), 0);

    for (IndexType i(0); i < indexMap.size(); ++i)
        indexMap[i] = i;

    std::sort(indexMap.begin(), indexMap.end(), compareThroughIndicesRows< Matrix< mpreal, Dynamic, Dynamic > >(matrixR, ascending));

    return indexMap;
}






/* This function transforms complex matrices into twice as big real matrices */
GmpEigenMatrix GmpEigenMatrix::complexIsometry() const
{
    GmpEigenMatrix result(*this);

    result.matrixR.conservativeResize(2*matrixR.rows(), 2*matrixR.cols());
    result.matrixR.block(matrixR.rows(), matrixR.cols(), matrixR.rows(), matrixR.cols()) = matrixR;
    if (isComplex) {
        result.matrixI.resize(0,0);
        result.isComplex = false;
        result.matrixR.block(0, matrixR.cols(), matrixR.rows(), matrixR.cols()) = matrixI;
        result.matrixR.block(matrixR.rows(), 0, matrixR.rows(), matrixR.cols()) = -matrixI;
    }

    return result;
}

/* This function restores the complex matrix corresponding to a big real
   matrix created with the complexIsometry function. */
GmpEigenMatrix GmpEigenMatrix::complexIsometryInverse() const
{
    GmpEigenMatrix result;

    if (isComplex) {
        mexErrMsgTxt("Error in complexIsometryInverse : the provided matrix is not real.");
    }

    result.matrixR = matrixR.block(0, 0, matrixR.rows()/2, matrixR.cols()/2);
    result.matrixI = matrixR.block(0, matrixR.cols()/2, matrixR.rows()/2, matrixR.cols()/2);
    result.checkComplexity();

    return result;
}





// Random matrix
GmpEigenMatrix gemRand(const IndexType& m, const IndexType& n)
{
    GmpEigenMatrix result;

    result.matrixR.resize(m,n);
    for (IndexType i(0); i < m; ++i)
        for (IndexType j(0); j < n; ++j)
            result.matrixR(i,j) = mpfr::random();

    return result;
}

// Random matrix
GmpEigenMatrix& gemRand_new(const IndexType& m, const IndexType& n)
{
    GmpEigenMatrix& result(*(new GmpEigenMatrix));

    result.matrixR.resize(m,n);
    for (IndexType i(0); i < m; ++i)
        for (IndexType j(0); j < n; ++j)
            result.matrixR(i,j) = mpfr::random();

    return result;
}
