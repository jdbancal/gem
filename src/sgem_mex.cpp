#include "mex.h"
#include "class_handle.hpp"
#include "gem.hpp"
#include "sgem.hpp"

/*
  This files receives instructions from matlab about which operation should be
  performed with the c++ library, and it calls the class to execute these
  instructions.

  Note that classes object going to matlab need to be allocated dynamically
  (otherwise they won't survive between two calls to this interface). Some
  special care is thus taken when dealing with these objects (such as calling
  functions ending with "_new" that allocate memory for such objects). Doing so
  allows one to rely on different parts of the code for memory management. Apart
  from the object creation, the only parts of the code which require manipulating
  pointers are in the class_handle.hpp file. We thus choose to deal here with
  references to these objects, rather than pointers. This allows for a cleaner
  code and more reliable memory management.
*/

/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


using namespace std;
using namespace Eigen;


// The following function is supposed to deal with all the memory allocation
// by itself.
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* Get the command string */
    char cmd[64];
    if (nrhs < 1 || mxGetString(prhs[0], cmd, sizeof(cmd)))
        mexErrMsgTxt("First input should be a command string less than 64 characters long.");

    /* New */
    if (!strcmp("new", cmd)) {
        // Check parameters
        if (nlhs != 1)
            mexErrMsgTxt("New: One output expected.");
        if (nrhs == 1)
        {
            // Called with no option, return the identifier to a new C++ instance:
            // for this, we send a reference to a new instance to the createMatlabIdFromObj
            // function.
            plhs[0] = createMatlabIdFromObj<SparseGmpEigenMatrix>(*(new SparseGmpEigenMatrix));
        }
        else if (nrhs == 2)
        {
            // Constructor called with one additional argument; we will call
            // the copy constructor: return a copy of the provided C++ instance

            // First, we restore the object to be copied from the matlab integer
            SparseGmpEigenMatrix& SparseGmpEigenMatrix_instance = recoverObjFromMatlabId<SparseGmpEigenMatrix>(prhs[1]);

            // Now we create a copy of this object
            SparseGmpEigenMatrix& copiedInstance(*(new SparseGmpEigenMatrix(SparseGmpEigenMatrix_instance)));

            //And return it (in the form of a matlab reference to it)
            plhs[0] = createMatlabIdFromObj<SparseGmpEigenMatrix>(copiedInstance);
        }
        else
            mexErrMsgTxt("New: Too many arguments.");
        return;
    }


    /* New sparse object from a dense gem object */
    if (!strcmp("sparse", cmd)) {
        // We check that all the necessary arguments are provided
        if ((nlhs != 1) || (nrhs != 3))
            mexErrMsgTxt("sparse: Unexpected arguments.");

        // We recover the dense object from the second argument
        GmpEigenMatrix& GmpEigenMatrix_instance = recoverObjFromMatlabId<GmpEigenMatrix>(prhs[1]);

        // We extract the required precision
        GmpEigenMatrix& threshold = recoverObjFromMatlabId<GmpEigenMatrix>(prhs[2]);

        // Now we create the new sparse object
        plhs[0] = createMatlabIdFromObj<SparseGmpEigenMatrix>(*(new SparseGmpEigenMatrix(GmpEigenMatrix_instance, threshold.getFirstRealValue())));

        return;
    }


    /* New from Matlab array */
    if (!strcmp("newFromMatlab", cmd)) {
        // We check that all the necessary arguments are provided
        if ((nlhs != 1) || (nrhs != 7))
            mexErrMsgTxt("newFromMatlab: Unexpected arguments.");

        if (!(mxIsDouble(prhs[1])) || !(mxIsDouble(prhs[2])) || !(mxIsDouble(prhs[4])) || !(mxIsDouble(prhs[5])) || !(mxIsDouble(prhs[6])))
            mexErrMsgTxt("newFromMatlab: Unexpected arguments.");

        // We extract the size and required precision
        IndexType m = mxGetScalar(prhs[4]);
        IndexType n = mxGetScalar(prhs[5]);
        int precision = mxGetScalar(prhs[6]);
        if (mxIsDouble(prhs[3])) {
            // The values have been provided in double format

            // We now forward the matlab pointer to the data to the constructor
            // of GmpEigenMatrix to let it copy the data into a new object. We then
            // send back the identifier of the newly created object to matlab.
            plhs[0] = createMatlabIdFromObj<SparseGmpEigenMatrix>(*(new SparseGmpEigenMatrix(prhs[1], prhs[2], prhs[3], m, n, precision)));
        } else {
            // The values have been provided in a gem variable

            // We extract this gem object
            GmpEigenMatrix& GmpEigenMatrix_values = recoverObjFromMatlabId<GmpEigenMatrix>(prhs[3]);

            // We now forward the matlab pointer to the data to the constructor
            // of GmpEigenMatrix to let it copy the data into a new object. We then
            // send back the identifier of the newly created object to matlab.
            plhs[0] = createMatlabIdFromObj<SparseGmpEigenMatrix>(*(new SparseGmpEigenMatrix(prhs[1], prhs[2], GmpEigenMatrix_values, m, n, precision)));
        }

        return;
    }

    /* Delete */
    if (!strcmp("delete", cmd)) {
        // Check parameters
        if (nrhs != 2)
            mexErrMsgTxt("Delete: Wrong number of arguments.");

        // Destroy the C++ object
        destroyObject<SparseGmpEigenMatrix>(prhs[1]);

        return;
    }


    /* Save to a Matlab structure */
    if (!strcmp("saveobj", cmd)) {
        // We check that the parameters are correct
        if ((nlhs != 1) || (nrhs != 2))
            mexErrMsgTxt("saveobj: Unexpected arguments.");

        // Now return a pointer to the created structure
        plhs[0] = recoverObjFromMatlabId<SparseGmpEigenMatrix>(prhs[1]).saveobj();

        return;
    }


    /* Load from a Matlab structure */
    if (!strcmp("loadobj", cmd)) {
        // We check that a second object is given
        if ((nlhs != 1) || (nrhs != 2))
            mexErrMsgTxt("loadobj: Unexpected arguments.");

        // We check that the second object is a matlab struct
        if ((!(mxIsStruct(prhs[1]))) || (mxGetNumberOfElements(prhs[1]) != 1))
            mexErrMsgTxt("loadobj: Argument should be a matlab struct array with 1 element");

        // Now we forward the pointer to the structure to the gem constructor which is
        // tailored for structure arrays.
        plhs[0] = createMatlabIdFromObj<SparseGmpEigenMatrix>(*(new SparseGmpEigenMatrix(prhs[1])));

        return;
    }


    /* Sets the working precision */
    if (!strcmp("setWorkingPrecision", cmd)) {
        // We check that the parameters are correct
        if ((nlhs != 0) || (nrhs != 2))
            mexErrMsgTxt("setWorkingPrecision: Unexpected arguments.");

        // We extract the required precision
        int precision = mxGetScalar(prhs[1]);

        // And set it as default for all mpfr manipulations
        setDefaultPrecForAllThreads(mpfr::digits2bits(precision));

        return;
    }


    // If we reached here, then there must be a second input parameter, and the first
    // input parameter refers to a gem object.
    if (nrhs < 2)
	   mexErrMsgTxt("Second input not found.");

    // We reactivate the class instance referenced by the second input
    // All procedures below don't need it, but it costs virtually nothing.
    SparseGmpEigenMatrix& SparseGmpEigenMatrix_instance = recoverObjFromMatlabId<SparseGmpEigenMatrix>(prhs[1]);


    /* isValid */
    if (!strcmp("isValid", cmd)) {
        // Check parameters
        if ((nlhs != 1) || (nrhs != 2))
            mexErrMsgTxt("isValid: Unexpected arguments.");

        // We allocate space for the result
        plhs[0] = mxCreateNumericMatrix(1, 1, mxUINT8_CLASS, mxREAL);

        // We check where the output data should be places
        int* outputMatrix = (int*)mxGetData(plhs[0]);

        // Call the method
        bool validity(checkValidity<SparseGmpEigenMatrix>(prhs[1]));

        // And stock it at the right place
        if (validity)
            outputMatrix[0] = 1;
        else
            outputMatrix[0] = 0;
        return;
    }


    /* display */
    if (!strcmp("display", cmd)) {
        // Check parameters
        if ((nlhs != 0) || (nrhs != 3))
            mexErrMsgTxt("display: Unexpected arguments.");

        // We extract the required precision
        int precision = mxGetScalar(prhs[2]);

        // This function prints each component individually.
        SparseGmpEigenMatrix_instance.displayIndividual(precision);

        return;
    }


    /* Extract a matlab table */
    if (!strcmp("double", cmd)) {
        // Check parameters
        if ((nlhs != 1) || (nrhs != 2))
            mexErrMsgTxt("double: Unexpected arguments.");

        plhs[0] = SparseGmpEigenMatrix_instance.toDouble();
        return;
    }


    /* Extract a cell array with strings describing each number */
    if (!strcmp("toStrings", cmd)) {
        // Check parameters
        if ((nlhs != 1) || (nrhs != 3))
            mexErrMsgTxt("toStrings: Unexpected arguments.");

        // We extract the required precision
        int precision = mxGetScalar(prhs[2]);

        // ... and call the function
        plhs[0] = SparseGmpEigenMatrix_instance.toStrings(precision);
        return;
    }


    /* Extracts the number of significant digits for each matrix element */
    if (!strcmp("precision", cmd)) {
        // Check parameters
        if ((nlhs != 1) || (nrhs != 2))
            mexErrMsgTxt("precision: Unexpected arguments.");

        // We return to matlab the precision table
        plhs[0] = SparseGmpEigenMatrix_instance.precision();

        return;
    }


    /* Chops components with a magnitude smaller than a given threshold */
    if (!strcmp("clean", cmd)) {
        // Check parameters
        if ((nlhs != 1) || (nrhs != 3))
            mexErrMsgTxt("clean: Unexpected arguments.");

        // We get the tolrance
        GmpEigenMatrix& tolerance = recoverObjFromMatlabId<GmpEigenMatrix>(prhs[2]);

        // We compute the cleaned matrix
        SparseGmpEigenMatrix& result(SparseGmpEigenMatrix_instance.clean_new(tolerance));

        // We prepare the reference to this object for matlab
        plhs[0] = createMatlabIdFromObj<SparseGmpEigenMatrix>(result);
        return;
    }


    /* Extract the nonzero elements from the matrix */
    if (!strcmp("find", cmd)) {
        // Check parameters
        if ((nlhs != 3) || (nrhs != 2))
            mexErrMsgTxt("find: Unexpected arguments.");

        vector < double > rows, cols;
        GmpEigenMatrix& values(SparseGmpEigenMatrix_instance.find_new(rows, cols));

        plhs[0] = vectorToMatlabDoubles(rows);
        plhs[1] = vectorToMatlabDoubles(cols);
        plhs[2] = createMatlabIdFromObj<GmpEigenMatrix>(values);
        return;
    }


    /* Extract a the size of the matrix */
    if (!strcmp("size", cmd)) {
        // Check parameters
        if ((nlhs != 1) || (nrhs != 2))
            mexErrMsgTxt("size: Unexpected arguments.");

        plhs[0] = SparseGmpEigenMatrix_instance.size();
        return;
    }


    /* Resize a matrix */
    if (!strcmp("resize", cmd)) {
        // Check parameters
        if ((nlhs != 0) || (nrhs != 4))
            mexErrMsgTxt("resize: Unexpected arguments.");

        // We extract the new dimension
        IndexType m2(mxGetScalar(prhs[2]));
        IndexType n2(mxGetScalar(prhs[3]));

        // And resize the current object
        SparseGmpEigenMatrix_instance.resize(m2, n2);

        return;
    }


    /* Extract a sub-matrix */
    if (!strcmp("subsref", cmd)) {
        // Check parameters
        if ((nlhs != 1) || (nrhs < 3) || (nrhs > 4))
            mexErrMsgTxt("subsref: Unexpected arguments.");

        // Extract the sub-matrix
        if (nrhs == 3) {
            // We extract the individual indices
            vector < vector < IndexType > > indices(matlabDoublesToMatrixVV<IndexType>(prhs[2]));

            SparseGmpEigenMatrix& result(SparseGmpEigenMatrix_instance.subsref_new(indices));

            // We prepare the reference to this object for matlab
            plhs[0] = createMatlabIdFromObj<SparseGmpEigenMatrix>(result);
        } else if (nrhs == 4) {
            // We extract the individual indices
            vector < IndexType > indicesA(matlabDoublesToVector<IndexType>(prhs[2]));
            vector < IndexType > indicesB(matlabDoublesToVector<IndexType>(prhs[3]));

            SparseGmpEigenMatrix& result(SparseGmpEigenMatrix_instance.subsref_new(indicesA, indicesB));

            // We prepare the reference to this object for matlab
            plhs[0] = createMatlabIdFromObj<SparseGmpEigenMatrix>(result);
        }

        return;
    }


    /* Reshape the matrix */
    if (!strcmp("reshape", cmd)) {
        // Check parameters
        if ((nlhs != 1) || (nrhs != 4))
            mexErrMsgTxt("reshape: Unexpected arguments.");

        // We extract the new dimension
        IndexType m2(mxGetScalar(prhs[2]));
        IndexType n2(mxGetScalar(prhs[3]));

        SparseGmpEigenMatrix& result(SparseGmpEigenMatrix_instance.reshape_new(m2, n2));

        // We prepare the reference to this object for matlab
        plhs[0] = createMatlabIdFromObj<SparseGmpEigenMatrix>(result);

        return;
    }


    /* Create a matrix from its diagonal coefficients */
    if (!strcmp("diagCreate", cmd)) {
        // Check parameters
        if ((nlhs != 1) || (nrhs != 3))
            mexErrMsgTxt("diagCreate: Unexpected arguments.");

        // We extract the diagonal index
        IndexType k(mxGetScalar(prhs[2]));

        SparseGmpEigenMatrix& result(SparseGmpEigenMatrix_instance.diagCreate_new(k));

        // We prepare the reference to this object for matlab
        plhs[0] = createMatlabIdFromObj<SparseGmpEigenMatrix>(result);

        return;
    }


    /* Extracts the diagonal from a matrix */
    if (!strcmp("diagExtract", cmd)) {
        // Check parameters
        if ((nlhs != 1) || (nrhs != 3))
            mexErrMsgTxt("diagExtract: Unexpected arguments.");

        // We extract the diagonal index
        IndexType k(mxGetScalar(prhs[2]));

        SparseGmpEigenMatrix& result(SparseGmpEigenMatrix_instance.diagExtract_new(k));

        // We prepare the reference to this object for matlab
        plhs[0] = createMatlabIdFromObj<SparseGmpEigenMatrix>(result);

        return;
    }


    /* Call the class method "uminus" */
    if (!strcmp("uminus", cmd)) {
        // Check parameters
        if ((nlhs != 1) || (nrhs != 2))
            mexErrMsgTxt("uminus: Unexpected arguments.");

        // Compute the difference
        SparseGmpEigenMatrix& result(SparseGmpEigenMatrix_instance.uminus_new());

        // We return the reference to this object to matlab
        plhs[0] = createMatlabIdFromObj<SparseGmpEigenMatrix>(result);

        return;
    }


    /* Call the class method "transpose" */
    if (!strcmp("transpose", cmd)) {
        // Check parameters
        if ((nlhs != 1) || (nrhs != 2))
            mexErrMsgTxt("transpose: Unexpected arguments.");

        // Compute the transposition
        SparseGmpEigenMatrix& result(SparseGmpEigenMatrix_instance.transpose_new());

        // We return the reference to this object to matlab
        plhs[0] = createMatlabIdFromObj<SparseGmpEigenMatrix>(result);

        return;
    }


    /* Call the class method "conj" */
    if (!strcmp("conj", cmd)) {
        // Check parameters
        if ((nlhs != 1) || (nrhs != 2))
            mexErrMsgTxt("conj: Unexpected arguments.");

        // Compute the complex conjugate
        SparseGmpEigenMatrix& result(SparseGmpEigenMatrix_instance.conj_new());

        // We return the reference to this object to matlab
        plhs[0] = createMatlabIdFromObj<SparseGmpEigenMatrix>(result);

        return;
    }


    /* Call the class method "ctranspose" */
    if (!strcmp("ctranspose", cmd)) {
        // Check parameters
        if ((nlhs != 1) || (nrhs != 2))
            mexErrMsgTxt("ctranspose: Unexpected arguments.");

        // Compute the hermitian conjugate
        SparseGmpEigenMatrix& result(SparseGmpEigenMatrix_instance.ctranspose_new());

        // We return the reference to this object to matlab
        plhs[0] = createMatlabIdFromObj<SparseGmpEigenMatrix>(result);

        return;
    }


    /* Call the class method "real" */
    if (!strcmp("real", cmd)) {
        // Check parameters
        if ((nlhs != 1) || (nrhs != 2))
            mexErrMsgTxt("real: Unexpected arguments.");

        // Extracts the real part
        SparseGmpEigenMatrix& result(SparseGmpEigenMatrix_instance.real_new());

        // We return the reference to this object to matlab
        plhs[0] = createMatlabIdFromObj<SparseGmpEigenMatrix>(result);

        return;
    }


    /* Call the class method "imag" */
    if (!strcmp("imag", cmd)) {
        // Check parameters
        if ((nlhs != 1) || (nrhs != 2))
            mexErrMsgTxt("imag: Unexpected arguments.");

        // Extracts the imaginary part
        SparseGmpEigenMatrix& result(SparseGmpEigenMatrix_instance.imag_new());

        // We return the reference to this object to matlab
        plhs[0] = createMatlabIdFromObj<SparseGmpEigenMatrix>(result);

        return;
    }


    /* Call the class method "isreal" */
    if (!strcmp("isreal", cmd)) {
        // Check parameters
        if ((nlhs != 1) || (nrhs != 2))
            mexErrMsgTxt("isreal: Unexpected arguments.");

        // We return the answer to matlab
        plhs[0] = mxCreateLogicalScalar(SparseGmpEigenMatrix_instance.isreal());

        return;
    }


    /* Call the class method "nnz" */
    if (!strcmp("nnz", cmd)) {
        // Check parameters
        if ((nlhs != 1) || (nrhs != 2))
            mexErrMsgTxt("nnz: Unexpected arguments.");

        // We allocate space for the result
        plhs[0] = mxCreateNumericMatrix(1, 1, mxUINT64_CLASS, mxREAL);

        // We check where the output data should be places
        int* output = (int*)mxGetData(plhs[0]);

        // We return the answer to matlab
        output[0] = SparseGmpEigenMatrix_instance.nnz();

        return;
    }


    /* Call the class method "round" */
    if (!strcmp("round", cmd)) {
        // Check parameters
        if ((nlhs != 1) || (nrhs != 2))
            mexErrMsgTxt("round: Unexpected arguments.");

        // Computes the round
        SparseGmpEigenMatrix& result(SparseGmpEigenMatrix_instance.round_new());

        // We return the reference to this object to matlab
        plhs[0] = createMatlabIdFromObj<SparseGmpEigenMatrix>(result);

        return;
    }


    /* Call the class method "floor" */
    if (!strcmp("floor", cmd)) {
        // Check parameters
        if ((nlhs != 1) || (nrhs != 2))
            mexErrMsgTxt("floor: Unexpected arguments.");

        // Computes the floor
        SparseGmpEigenMatrix& result(SparseGmpEigenMatrix_instance.floor_new());

        // We return the reference to this object to matlab
        plhs[0] = createMatlabIdFromObj<SparseGmpEigenMatrix>(result);

        return;
    }


    /* Call the class method "ceil" */
    if (!strcmp("ceil", cmd)) {
        // Check parameters
        if ((nlhs != 1) || (nrhs != 2))
            mexErrMsgTxt("ceil: Unexpected arguments.");

        // Computes the ceil
        SparseGmpEigenMatrix& result(SparseGmpEigenMatrix_instance.ceil_new());

        // We return the reference to this object to matlab
        plhs[0] = createMatlabIdFromObj<SparseGmpEigenMatrix>(result);

        return;
    }


    /* Call the class method "trunc" */
    if (!strcmp("trunc", cmd)) {
        // Check parameters
        if ((nlhs != 1) || (nrhs != 2))
            mexErrMsgTxt("trunc: Unexpected arguments.");

        // Computes the trunc
        SparseGmpEigenMatrix& result(SparseGmpEigenMatrix_instance.trunc_new());

        // We return the reference to this object to matlab
        plhs[0] = createMatlabIdFromObj<SparseGmpEigenMatrix>(result);

        return;
    }


    /* Call the class method "abs" */
    if (!strcmp("abs", cmd)) {
        // Check parameters
        if ((nlhs != 1) || (nrhs != 2))
            mexErrMsgTxt("abs: Unexpected arguments.");

        // Compute the absolute value (or complex magnitude)
        SparseGmpEigenMatrix& result(SparseGmpEigenMatrix_instance.abs_new());

        // We return the reference to this object to matlab
        plhs[0] = createMatlabIdFromObj<SparseGmpEigenMatrix>(result);

        return;
    }


    /* Call the class method "angle" */
    if (!strcmp("angle", cmd)) {
        // Check parameters
        if ((nlhs != 1) || (nrhs != 2))
            mexErrMsgTxt("angle: Unexpected arguments.");

        // Compute the phase angle
        SparseGmpEigenMatrix& result(SparseGmpEigenMatrix_instance.angle_new());

        // We return the reference to this object to matlab
        plhs[0] = createMatlabIdFromObj<SparseGmpEigenMatrix>(result);

        return;
    }


    // // The two following methods are not needed by the matlab interface :
    //
    // /* Call the class method "exp" */
    // if (!strcmp("exp", cmd)) {
    //     // Check parameters
    //     if ((nlhs != 1) || (nrhs != 2))
    //         mexErrMsgTxt("exp: Unexpected arguments.");
    //
    //     // Compute the exponential
    //     SparseGmpEigenMatrix& result(SparseGmpEigenMatrix_instance.exp_new());
    //
    //     // We return the reference to this object to matlab
    //     plhs[0] = createMatlabIdFromObj<SparseGmpEigenMatrix>(result);
    //
    //     return;
    // }
    //
    //
    // /* Call the class method "log" */
    // if (!strcmp("log", cmd)) {
    //     // Check parameters
    //     if ((nlhs != 1) || (nrhs != 2))
    //         mexErrMsgTxt("log: Unexpected arguments.");
    //
    //     // Compute the natural logarithm
    //     SparseGmpEigenMatrix& result(SparseGmpEigenMatrix_instance.log_new());
    //
    //     // We return the reference to this object to matlab
    //     plhs[0] = createMatlabIdFromObj<SparseGmpEigenMatrix>(result);
    //
    //     return;
    // }


    /* Call the class method "sqrt" */
    if (!strcmp("sqrt", cmd)) {
        // Check parameters
        if ((nlhs != 1) || (nrhs != 2))
            mexErrMsgTxt("sqrt: Unexpected arguments.");

        // Compute the square root
        SparseGmpEigenMatrix& result(SparseGmpEigenMatrix_instance.sqrt_new());

        // We return the reference to this object to matlab
        plhs[0] = createMatlabIdFromObj<SparseGmpEigenMatrix>(result);

        return;
    }


    /* Call the class method "cbrt" */
    if (!strcmp("cbrt", cmd)) {
        // Check parameters
        if ((nlhs != 1) || (nrhs != 2))
            mexErrMsgTxt("cbrt: Unexpected arguments.");

        // Compute the cubic root
        SparseGmpEigenMatrix& result(SparseGmpEigenMatrix_instance.cbrt_new());

        // We return the reference to this object to matlab
        plhs[0] = createMatlabIdFromObj<SparseGmpEigenMatrix>(result);

        return;
    }


    /* Call the class method "sin" */
    if (!strcmp("sin", cmd)) {
        // Check parameters
        if ((nlhs != 1) || (nrhs != 2))
            mexErrMsgTxt("sin: Unexpected arguments.");

        // Compute the cubic root
        SparseGmpEigenMatrix& result(SparseGmpEigenMatrix_instance.sin_new());

        // We return the reference to this object to matlab
        plhs[0] = createMatlabIdFromObj<SparseGmpEigenMatrix>(result);

        return;
    }


    /* Call the class method "tan" */
    if (!strcmp("tan", cmd)) {
        // Check parameters
        if ((nlhs != 1) || (nrhs != 2))
            mexErrMsgTxt("tan: Unexpected arguments.");

        // Compute the cubic root
        SparseGmpEigenMatrix& result(SparseGmpEigenMatrix_instance.tan_new());

        // We return the reference to this object to matlab
        plhs[0] = createMatlabIdFromObj<SparseGmpEigenMatrix>(result);

        return;
    }


    /* Call the class method "asin" */
    if (!strcmp("asin", cmd)) {
        // Check parameters
        if ((nlhs != 1) || (nrhs != 2))
            mexErrMsgTxt("asin: Unexpected arguments.");

        // Compute the cubic root
        SparseGmpEigenMatrix& result(SparseGmpEigenMatrix_instance.asin_new());

        // We return the reference to this object to matlab
        plhs[0] = createMatlabIdFromObj<SparseGmpEigenMatrix>(result);

        return;
    }


    /* Call the class method "atan" */
    if (!strcmp("atan", cmd)) {
        // Check parameters
        if ((nlhs != 1) || (nrhs != 2))
            mexErrMsgTxt("atan: Unexpected arguments.");

        // Compute the cubic root
        SparseGmpEigenMatrix& result(SparseGmpEigenMatrix_instance.atan_new());

        // We return the reference to this object to matlab
        plhs[0] = createMatlabIdFromObj<SparseGmpEigenMatrix>(result);

        return;
    }






    /* Call the class method "rank" */
    if (!strcmp("rank", cmd)) {
        // Check parameters
        if ((nlhs != 1) || (nrhs != 2))
            mexErrMsgTxt("rank: Unexpected arguments.");

        // We allocate space for the result
        plhs[0] = mxCreateNumericMatrix(1, 1, mxUINT64_CLASS, mxREAL);

        // We check where the output data should be places
        int* output = (int*)mxGetData(plhs[0]);

        // We return the answer to matlab
        output[0] = SparseGmpEigenMatrix_instance.rank();

        return;
    }


    /* Call the class method "mldivide_sf" */
    if (!strcmp("mldivide_sf", cmd)) {
        // Check parameters
        if ((nlhs != 1) || (nrhs != 3))
            mexErrMsgTxt("mldivide_sf: Unexpected arguments.");

        // Get the second argument
        GmpEigenMatrix& b = recoverObjFromMatlabId<GmpEigenMatrix>(prhs[2]);

        // Solve the linear system
        GmpEigenMatrix& result(SparseGmpEigenMatrix_instance.mldivide_sf_new(b));

        // We return the reference to these objects to matlab
        plhs[0] = createMatlabIdFromObj<GmpEigenMatrix>(result);

        return;
    }


    /* Call the class method "mldivide" */
    if (!strcmp("mldivide", cmd)) {
        // Check parameters
        if ((nlhs != 1) || (nrhs != 3))
            mexErrMsgTxt("mldivide: Unexpected arguments.");

        // Get the second argument
        SparseGmpEigenMatrix& b = recoverObjFromMatlabId<SparseGmpEigenMatrix>(prhs[2]);

        // Solve the linear system
        SparseGmpEigenMatrix& result(SparseGmpEigenMatrix_instance.mldivide_new(b));

        // We return the reference to these objects to matlab
        plhs[0] = createMatlabIdFromObj<SparseGmpEigenMatrix>(result);

        return;
    }


    /* Call the class method "inv" */
    if (!strcmp("inv", cmd)) {
        // Check parameters
        if ((nlhs != 1) || (nrhs != 2))
            mexErrMsgTxt("inv: Unexpected arguments.");

        // Compute the matric inverse
        SparseGmpEigenMatrix& result(SparseGmpEigenMatrix_instance.inv_new());

        // We return the reference to these objects to matlab
        plhs[0] = createMatlabIdFromObj<SparseGmpEigenMatrix>(result);

        return;
    }


    /* Call the class method "eigs" */
    if (!strcmp("eigs", cmd)) {
        // Check parameters
        if ((nlhs != 2) || (nrhs != 5))
            mexErrMsgTxt("eigs: Unexpected arguments.");

        // We extract the number of eigenvalues requested
        long int nbEigenvalues = mxGetScalar(prhs[2]);
        long int type = mxGetScalar(prhs[3]);
        GmpEigenMatrix& sigma = recoverObjFromMatlabId<GmpEigenMatrix>(prhs[4]);

        // Compute the eigen decomposition
        GmpEigenMatrix& Vmatrix(*(new GmpEigenMatrix));
        GmpEigenMatrix& result(SparseGmpEigenMatrix_instance.eigs_new(nbEigenvalues, Vmatrix, type, sigma));

        // We return the reference to these objects to matlab
        plhs[0] = createMatlabIdFromObj<GmpEigenMatrix>(Vmatrix);
        plhs[1] = createMatlabIdFromObj<GmpEigenMatrix>(result);

        return;
    }







    /* Call the class method "isnan" */
    if (!strcmp("isnan", cmd)) {
        // Check parameters
        if ((nlhs != 1) || (nrhs != 2))
            mexErrMsgTxt("isnan: Unexpected arguments.");

        // Compute the element-wise comparison
        SparseMatrix<bool> result(SparseGmpEigenMatrix_instance.isnan());

        // We return the reference to this object to matlab
        plhs[0] = matrixToMatlabDoubles(result);

        return;
    }


    /* Call the class method "isinf" */
    if (!strcmp("isinf", cmd)) {
        // Check parameters
        if ((nlhs != 1) || (nrhs != 2))
            mexErrMsgTxt("isinf: Unexpected arguments.");

        // Compute the element-wise comparison
        SparseMatrix<bool> result(SparseGmpEigenMatrix_instance.isinf());

        // We return the reference to this object to matlab
        plhs[0] = matrixToMatlabDoubles(result);

        return;
    }


    /* Call the class method "issymmetric" */
    if (!strcmp("issymmetric", cmd)) {
        // Check parameters
        if ((nlhs != 1) || (nrhs != 2))
            mexErrMsgTxt("issymmetric: Unexpected arguments.");

        // We allocate space for the result
        plhs[0] = mxCreateNumericMatrix(1, 1, mxUINT8_CLASS, mxREAL);

        // We check where the output data should be places
        int* outputMatrix = (int*)mxGetData(plhs[0]);

        // Check whether the matrix is symmetric
        bool result(SparseGmpEigenMatrix_instance.issymmetric());

        // And stock the result at the right place
        if (result)
            outputMatrix[0] = 1;
        else
            outputMatrix[0] = 0;

        return;
    }


    /* Call the class method "ishermitian" */
    if (!strcmp("ishermitian", cmd)) {
        // Check parameters
        if ((nlhs != 1) || (nrhs != 2))
            mexErrMsgTxt("ishermitian: Unexpected arguments.");

            // We allocate space for the result
            plhs[0] = mxCreateNumericMatrix(1, 1, mxUINT8_CLASS, mxREAL);

            // We check where the output data should be places
            int* outputMatrix = (int*)mxGetData(plhs[0]);

            // Check whether the matrix is symmetric
            bool result(SparseGmpEigenMatrix_instance.ishermitian());

            // And stock the result at the right place
            if (result)
                outputMatrix[0] = 1;
            else
                outputMatrix[0] = 0;

        return;
    }





    /* Call the class method "colMin" */
    if (!strcmp("colMin", cmd)) {
        // Check parameters
        if ((nlhs != 2) || (nrhs != 2))
            mexErrMsgTxt("colMin: Unexpected arguments.");

        // We prepare a table to gather the optimal indices
        vector<IndexType> indices;

        // Compute the column minimum
        SparseGmpEigenMatrix& result(SparseGmpEigenMatrix_instance.colMin_new(indices));

        // We return the reference to this object to matlab
        plhs[0] = createMatlabIdFromObj<SparseGmpEigenMatrix>(result);
        plhs[1] = vectorToMatlabDoubles(indices);

        return;
    }


    /* Call the class method "rowMin" */
    if (!strcmp("rowMin", cmd)) {
        // Check parameters
        if ((nlhs != 2) || (nrhs != 2))
            mexErrMsgTxt("rowMin: Unexpected arguments.");

        // We prepare a table to gather the optimal indices
        vector<IndexType> indices;

        // Compute the row minimum
        SparseGmpEigenMatrix& result(SparseGmpEigenMatrix_instance.rowMin_new(indices));

        // We return the reference to this object to matlab
        plhs[0] = createMatlabIdFromObj<SparseGmpEigenMatrix>(result);
        plhs[1] = vectorToMatlabDoubles(indices);

        return;
    }


    /* Call the class method "colMax" */
    if (!strcmp("colMax", cmd)) {
        // Check parameters
        if ((nlhs != 2) || (nrhs != 2))
            mexErrMsgTxt("colMax: Unexpected arguments.");

        // We prepare a table to gather the optimal indices
        vector<IndexType> indices;

        // Compute the column maximum
        SparseGmpEigenMatrix& result(SparseGmpEigenMatrix_instance.colMax_new(indices));

        // We return the reference to this object to matlab
        plhs[0] = createMatlabIdFromObj<SparseGmpEigenMatrix>(result);
        plhs[1] = vectorToMatlabDoubles(indices);

        return;
    }


    /* Call the class method "rowMax" */
    if (!strcmp("rowMax", cmd)) {
        // Check parameters
        if ((nlhs != 2) || (nrhs != 2))
            mexErrMsgTxt("rowMax: Unexpected arguments.");

        // We prepare a table to gather the optimal indices
        vector<IndexType> indices;

        // Compute the row maximum
        SparseGmpEigenMatrix& result(SparseGmpEigenMatrix_instance.rowMax_new(indices));

        // We return the reference to this object to matlab
        plhs[0] = createMatlabIdFromObj<SparseGmpEigenMatrix>(result);
        plhs[1] = vectorToMatlabDoubles(indices);

        return;
    }


//     // The two following methods are not needed by the matlab interface :
//
//     /* Call the class method "colSum" */
//     if (!strcmp("colSum", cmd)) {
//         // Check parameters
//         if ((nlhs != 1) || (nrhs != 2))
//             mexErrMsgTxt("colSum: Unexpected arguments.");
//
//         // Compute the column sum
//         GmpEigenMatrix& result(GmpEigenMatrix_instance.colSum_new());
//
//         // We return the reference to this object to matlab
//         plhs[0] = createMatlabIdFromObj<GmpEigenMatrix>(result);
//
//         return;
//     }
//
//
//     /* Call the class method "rowSum" */
//     if (!strcmp("rowSum", cmd)) {
//         // Check parameters
//         if ((nlhs != 1) || (nrhs != 2))
//             mexErrMsgTxt("rowSum: Unexpected arguments.");
//
//         // Compute the row sum
//         GmpEigenMatrix& result(GmpEigenMatrix_instance.rowSum_new());
//
//         // We return the reference to this object to matlab
//         plhs[0] = createMatlabIdFromObj<GmpEigenMatrix>(result);
//
//         return;
//     }


    /* Call the class method "colProd" */
    if (!strcmp("colProd", cmd)) {
        // Check parameters
        if ((nlhs != 1) || (nrhs != 2))
            mexErrMsgTxt("colProd: Unexpected arguments.");

        // Compute the column product
        SparseGmpEigenMatrix& result(SparseGmpEigenMatrix_instance.colProd_new());

        // We return the reference to this object to matlab
        plhs[0] = createMatlabIdFromObj<SparseGmpEigenMatrix>(result);

        return;
    }


    /* Call the class method "rowProd" */
    if (!strcmp("rowProd", cmd)) {
        // Check parameters
        if ((nlhs != 1) || (nrhs != 2))
            mexErrMsgTxt("rowProd: Unexpected arguments.");

        // Compute the row product
        SparseGmpEigenMatrix& result(SparseGmpEigenMatrix_instance.rowProd_new());

        // We return the reference to this object to matlab
        plhs[0] = createMatlabIdFromObj<SparseGmpEigenMatrix>(result);

        return;
    }


    /* Sorts out the elements of a matrix */
    if (!strcmp("sort", cmd)) {
        // Check parameters
        if ((nlhs != 3) || (nrhs != 4))
            mexErrMsgTxt("sort: Unexpected arguments.");

        // We extract the sorting dimension
        int dim(mxGetScalar(prhs[2]));

        // We extract the sorting type
        int type(mxGetScalar(prhs[3]));

        // We prepare a table to gather the indices
        vector < vector < IndexType > > indices;

        // We prepare a list of positions for the zeros
        vector < IndexType > nbNegatives;

        SparseGmpEigenMatrix& result(SparseGmpEigenMatrix_instance.sort_new(dim, type, indices, nbNegatives));

        // We prepare the reference to this object for matlab
        plhs[0] = createMatlabIdFromObj<SparseGmpEigenMatrix>(result);
        plhs[1] = vectorVectorToMatlabDoublesCell(indices);
        plhs[2] = vectorToMatlabDoubles(nbNegatives);

        return;
    }


    // If we reached here, then there must be a third input parameter
    if (nrhs < 3)
       mexErrMsgTxt("Third input not found.");


    // Here first come the methods which also involve full matrices

    /* Call the class method "times_sf" (element-wise multiplication between a sparse and a full matrix) */
    if (!strcmp("times_sf", cmd)) {
        // Check parameters
        if ((nlhs != 1) || (nrhs != 3))
            mexErrMsgTxt("times_sf: Unexpected arguments.");

        // We reactivate the class instance referenced by the second input
        // All procedures below don't need it, but it costs virtually nothing.
        GmpEigenMatrix& GmpEigenMatrix_instance2 = recoverObjFromMatlabId<GmpEigenMatrix>(prhs[2]);

        // Compute the multiplication
        SparseGmpEigenMatrix& result(SparseGmpEigenMatrix_instance.times_sf_new(GmpEigenMatrix_instance2));

        // We return the reference to this object to matlab
        plhs[0] = createMatlabIdFromObj<SparseGmpEigenMatrix>(result);

        return;
    }


    /* Call the class method "mtimes_sf" (matrix multiplication between a sparse and a full matrix) */
    if (!strcmp("mtimes_sf", cmd)) {
        // Check parameters
        if ((nlhs != 1) || (nrhs != 3))
            mexErrMsgTxt("mtimes_sf: Unexpected arguments.");

        // We reactivate the class instance referenced by the second input
        // All procedures below don't need it, but it costs virtually nothing.
        GmpEigenMatrix& GmpEigenMatrix_instance2 = recoverObjFromMatlabId<GmpEigenMatrix>(prhs[2]);

        // Compute the multiplication
        GmpEigenMatrix& result(SparseGmpEigenMatrix_instance.mtimes_sf_new(GmpEigenMatrix_instance2));

        // We return the reference to this object to matlab
        plhs[0] = createMatlabIdFromObj<GmpEigenMatrix>(result);

        return;
    }

    /* Call the class method "kron_sf" (kronecker tensor product between a sparse and a full matrix) */
    if (!strcmp("kron_sf", cmd)) {
        // Check parameters
        if ((nlhs != 1) || (nrhs != 3))
            mexErrMsgTxt("kron_sf: Unexpected arguments.");

        // We reactivate the class instance referenced by the second input
        // All procedures below don't need it, but it costs virtually nothing.
        GmpEigenMatrix& GmpEigenMatrix_instance2 = recoverObjFromMatlabId<GmpEigenMatrix>(prhs[2]);

        // Compute the element-wise multiplication
        SparseGmpEigenMatrix& result(SparseGmpEigenMatrix_instance.kron_sf_new(GmpEigenMatrix_instance2));

        // We return the reference to this object to matlab
        plhs[0] = createMatlabIdFromObj<SparseGmpEigenMatrix>(result);

        return;
    }


    /* Call the class method "power" (element-wise power) */
    if (!strcmp("power", cmd)) {
        // Check parameters
        if ((nlhs != 1) || (nrhs != 3))
            mexErrMsgTxt("power: Unexpected arguments.");

        // We reactivate the class instance referenced by the second input
        // All procedures below don't need it, but it costs virtually nothing.
        GmpEigenMatrix& GmpEigenMatrix_instance2 = recoverObjFromMatlabId<GmpEigenMatrix>(prhs[2]);

        // Compute the element-wise power
        SparseGmpEigenMatrix& result(SparseGmpEigenMatrix_instance.power_new(GmpEigenMatrix_instance2));

        // We return the reference to this object to matlab
        plhs[0] = createMatlabIdFromObj<SparseGmpEigenMatrix>(result);

        return;
    }

    // Now come the functions involving a second instance as well.
    // We reactivate this object from identifier present in the third input.
    SparseGmpEigenMatrix& SparseGmpEigenMatrix_instance2 = recoverObjFromMatlabId<SparseGmpEigenMatrix>(prhs[2]);


    /* Call the class method "plus" */
    if (!strcmp("plus", cmd)) {
        // Check parameters
        if ((nlhs != 1) || (nrhs != 3))
            mexErrMsgTxt("plus: Unexpected arguments.");

        /* The object to be returned to matlab needs to be created explicitely
        by hand (otherwise its memory would be dis-allocated upon exiting of
        this function). Therefore, we need to do things in a slightly unusual
        way here. Below are three examples of how to do that. */
        /* Here we use the the special addition function defined in our
        class, that also creates a new object. We keep track of the result
        produced with the help of a reference object.*/
        SparseGmpEigenMatrix& result(SparseGmpEigenMatrix_instance.plus_new(SparseGmpEigenMatrix_instance2));

        // Now we can return the reference to this object to matlab
        plhs[0] = createMatlabIdFromObj<SparseGmpEigenMatrix>(result);

        return;
    }

    /* Call the class method "minus" */
    if (!strcmp("minus", cmd)) {
        // Check parameters
        if ((nlhs != 1) || (nrhs != 3))
            mexErrMsgTxt("minus: Unexpected arguments.");

        // Compute the difference
        SparseGmpEigenMatrix& result(SparseGmpEigenMatrix_instance.minus_new(SparseGmpEigenMatrix_instance2));

        // We return the reference to this object to matlab
        plhs[0] = createMatlabIdFromObj<SparseGmpEigenMatrix>(result);

        return;
    }

    /* Call the class method "times" (element-wise multiplication) */
    if (!strcmp("times", cmd)) {
        // Check parameters
        if ((nlhs != 1) || (nrhs != 3))
            mexErrMsgTxt("times: Unexpected arguments.");

        // Compute the element-wise multiplication
        SparseGmpEigenMatrix& result(SparseGmpEigenMatrix_instance.times_new(SparseGmpEigenMatrix_instance2));

        // We return the reference to this object to matlab
        plhs[0] = createMatlabIdFromObj<SparseGmpEigenMatrix>(result);

        return;
    }

    /* Call the class method "rdivide" (element-wise division on the right) */
    if (!strcmp("rdivide", cmd)) {
        // Check parameters
        if ((nlhs != 1) || (nrhs != 3))
            mexErrMsgTxt("times: Unexpected arguments.");

        // Compute the element-wise ratio
        SparseGmpEigenMatrix& result(SparseGmpEigenMatrix_instance.rdivide_new(SparseGmpEigenMatrix_instance2));

        // We return the reference to this object to matlab
        plhs[0] = createMatlabIdFromObj<SparseGmpEigenMatrix>(result);

        return;
    }

    /* Call the class method "mtimes" (matrix multiplication) */
    if (!strcmp("mtimes", cmd)) {
        // Check parameters
        if ((nlhs != 1) || (nrhs != 3))
            mexErrMsgTxt("mtimes: Unexpected arguments.");

        // Compute the element-wise multiplication
        SparseGmpEigenMatrix& result(SparseGmpEigenMatrix_instance.mtimes_new(SparseGmpEigenMatrix_instance2));

        // We return the reference to this object to matlab
        plhs[0] = createMatlabIdFromObj<SparseGmpEigenMatrix>(result);

        return;
    }

    /* Call the class method "kron" (kronecker tensor product) */
    if (!strcmp("kron", cmd)) {
        // Check parameters
        if ((nlhs != 1) || (nrhs != 3))
            mexErrMsgTxt("kron: Unexpected arguments.");

        // Compute the element-wise multiplication
        SparseGmpEigenMatrix& result(SparseGmpEigenMatrix_instance.kron_new(SparseGmpEigenMatrix_instance2));

        // We return the reference to this object to matlab
        plhs[0] = createMatlabIdFromObj<SparseGmpEigenMatrix>(result);

        return;
    }

//     /* Call the class method "mpower" (matrix power) */
// /*    if (!strcmp("mpower", cmd)) {
//         // Check parameters
//         if ((nlhs != 1) || (nrhs != 3))
//             mexErrMsgTxt("mpower: Unexpected arguments.");
//
//         // Compute the element-wise power
//         GmpEigenMatrix& result(GmpEigenMatrix_instance.mpower_new(GmpEigenMatrix_instance2));
//
//         // We return the reference to this object to matlab
//         plhs[0] = createMatlabIdFromObj<GmpEigenMatrix>(result);
//
//         return;
//     }*/







    /* Call the class method "subsasgn" (assignment of elements) */
    if (!strcmp("subsasgn", cmd)) {
        // Check parameters
        if ((nlhs != 0) || (nrhs < 4) || (nrhs > 5))
            mexErrMsgTxt("subsasgn: Unexpected arguments.");

        if (nrhs == 4) {
            // We extract the individual indices
            vector < IndexType > indices(matlabDoublesToVector<IndexType>(prhs[3]));

            // We assign the new values
            SparseGmpEigenMatrix_instance.subsasgn(indices, SparseGmpEigenMatrix_instance2);
        } else {
            // We extract the individual indices
            vector < IndexType > indicesA(matlabDoublesToVector<IndexType>(prhs[3]));
            vector < IndexType > indicesB(matlabDoublesToVector<IndexType>(prhs[4]));

            // We assign the new values
            SparseGmpEigenMatrix_instance.subsasgn(indicesA, indicesB, SparseGmpEigenMatrix_instance2);
        }

        return;
    }


    /* Call the class method "vertcat" */
    if (!strcmp("vertcat", cmd)) {
        // Check parameters
        if ((nlhs != 1) || (nrhs != 3))
            mexErrMsgTxt("vertcat: Unexpected arguments.");

        // Concatenates the two matrices together
        SparseGmpEigenMatrix& result(SparseGmpEigenMatrix_instance.vertcat_new(SparseGmpEigenMatrix_instance2));

        // We return the reference to this object to matlab
        plhs[0] = createMatlabIdFromObj<SparseGmpEigenMatrix>(result);

        return;
    }


    /* Call the class method "horzcat" */
    if (!strcmp("horzcat", cmd)) {
        // Check parameters
        if ((nlhs != 1) || (nrhs != 3))
            mexErrMsgTxt("horzcat: Unexpected arguments.");

        // Concatenates the two matrices together
        SparseGmpEigenMatrix& result(SparseGmpEigenMatrix_instance.horzcat_new(SparseGmpEigenMatrix_instance2));

        // We return the reference to this object to matlab
        plhs[0] = createMatlabIdFromObj<SparseGmpEigenMatrix>(result);

        return;
    }







    /* Call the class method "lt" (element-wise less than) */
    if (!strcmp("lt", cmd)) {
        // Check parameters
        if ((nlhs != 1) || (nrhs != 3))
            mexErrMsgTxt("lt: Unexpected arguments.");

        // Compute the element-wise comparison
        SparseMatrix<bool> result(SparseGmpEigenMatrix_instance < SparseGmpEigenMatrix_instance2);

        // We return the reference to this object to matlab
        plhs[0] = matrixToMatlabDoubles(result);

        return;
    }


    /* Call the class method "le" (element-wise less than or equal) */
    if (!strcmp("le", cmd)) {
        // Check parameters
        if ((nlhs != 1) || (nrhs != 3))
            mexErrMsgTxt("le: Unexpected arguments.");

        // Compute the element-wise comparison
        SparseMatrix<bool> result(SparseGmpEigenMatrix_instance <= SparseGmpEigenMatrix_instance2);

        // We return the reference to this object to matlab
        plhs[0] = matrixToMatlabDoubles(result);

        return;
    }


    /* Call the class method "gt" (element-wise greater than) */
    if (!strcmp("gt", cmd)) {
        // Check parameters
        if ((nlhs != 1) || (nrhs != 3))
            mexErrMsgTxt("gt: Unexpected arguments.");

        // Compute the element-wise comparison
        SparseMatrix<bool> result(SparseGmpEigenMatrix_instance > SparseGmpEigenMatrix_instance2);

        // We return the reference to this object to matlab
        plhs[0] = matrixToMatlabDoubles(result);

        return;
    }


    /* Call the class method "ge" (element-wise greater than or equal) */
    if (!strcmp("ge", cmd)) {
        // Check parameters
        if ((nlhs != 1) || (nrhs != 3))
            mexErrMsgTxt("ge: Unexpected arguments.");

        // Compute the element-wise comparison
        SparseMatrix<bool> result(SparseGmpEigenMatrix_instance >= SparseGmpEigenMatrix_instance2);

        // We return the reference to this object to matlab
        plhs[0] = matrixToMatlabDoubles(result);

        return;
    }


    /* Call the class method "eq" (element-wise equality test) */
    if (!strcmp("eq", cmd)) {
        // Check parameters
        if ((nlhs != 1) || (nrhs != 3))
            mexErrMsgTxt("eq: Unexpected arguments.");

        // Compute the element-wise comparison
        SparseMatrix<bool> result(SparseGmpEigenMatrix_instance.eq(SparseGmpEigenMatrix_instance2));

        // We return the reference to this object to matlab
        plhs[0] = matrixToMatlabDoubles(result);

        return;
    }


    /* Call the class method "ne" (element-wise not equal) */
    if (!strcmp("ne", cmd)) {
        // Check parameters
        if ((nlhs != 1) || (nrhs != 3))
            mexErrMsgTxt("ne: Unexpected arguments.");

        // Compute the element-wise comparison
        SparseMatrix<bool> result(SparseGmpEigenMatrix_instance.ne(SparseGmpEigenMatrix_instance2));

        // We return the reference to this object to matlab
        plhs[0] = matrixToMatlabDoubles(result);

        return;
    }


    /* Call the class method "identicalValues" (checks whether all numerical values match, assuming all other characteristics already match) */
    if (!strcmp("identicalValues", cmd)) {
        // Check parameters
        if ((nlhs != 1) || (nrhs != 3))
            mexErrMsgTxt("identicalValues: Unexpected arguments.");

        // We allocate space for the result
        plhs[0] = mxCreateNumericMatrix(1, 1, mxUINT8_CLASS, mxREAL);

        // We check where the output data should be places
        int* outputMatrix = (int*)mxGetData(plhs[0]);

        // Check whether values match
        bool result(SparseGmpEigenMatrix_instance.identicalValues(SparseGmpEigenMatrix_instance2));

        // And stock the result at the right place
        if (result)
            outputMatrix[0] = 1;
        else
            outputMatrix[0] = 0;

        return;
    }


    /* Call the class method "identicalValuesNaNok" (checks whether all numerical values match, assuming all other characteristics already match) */
    if (!strcmp("identicalValuesNaNok", cmd)) {
        // Check parameters
        if ((nlhs != 1) || (nrhs != 3))
            mexErrMsgTxt("identicalValuesNaNok: Unexpected arguments.");

        // We allocate space for the result
        plhs[0] = mxCreateNumericMatrix(1, 1, mxUINT8_CLASS, mxREAL);

        // We check where the output data should be places
        int* outputMatrix = (int*)mxGetData(plhs[0]);

        // Check whether values match
        bool result(SparseGmpEigenMatrix_instance.identicalValuesNaNok(SparseGmpEigenMatrix_instance2));

        // And stock the result at the right place
        if (result)
            outputMatrix[0] = 1;
        else
            outputMatrix[0] = 0;

        return;
    }


    /* Call the class method "min" (element-wise minimum) */
    if (!strcmp("ewMin", cmd)) {
        // Check parameters
        if ((nlhs != 1) || (nrhs != 3))
            mexErrMsgTxt("ewMin: Unexpected arguments.");

        // Compute the element-wise multiplication
        SparseGmpEigenMatrix& result(SparseGmpEigenMatrix_instance.ewMin_new(SparseGmpEigenMatrix_instance2));

        // We return the reference to this object to matlab
        plhs[0] = createMatlabIdFromObj<SparseGmpEigenMatrix>(result);

        return;
    }


    /* Call the class method "max" (element-wise maximum) */
    if (!strcmp("ewMax", cmd)) {
        // Check parameters
        if ((nlhs != 1) || (nrhs != 3))
            mexErrMsgTxt("ewMax: Unexpected arguments.");

        // Compute the element-wise multiplication
        SparseGmpEigenMatrix& result(SparseGmpEigenMatrix_instance.ewMax_new(SparseGmpEigenMatrix_instance2));

        // We return the reference to this object to matlab
        plhs[0] = createMatlabIdFromObj<SparseGmpEigenMatrix>(result);

        return;
    }


    // Got here, so command not recognized
    mexErrMsgTxt("Command not recognized.");
}
