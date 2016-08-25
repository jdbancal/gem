#ifndef __CLASS_HANDLE_HPP__
#define __CLASS_HANDLE_HPP__
#include "mex.h"
#include <stdint.h>
#include <string>
#include <cstring>
#include <typeinfo>
#include <iostream>

/*
  This file provides three functions :
    - createMatlabIdFromObj
    - recoverObjFromMatlabId
    - destroyObject
  which allow one to easily keep track of c++ class objects allocated
  dynamically between two calls from matlab.

  A typical usage starts by calling the function "createMatlabIdFromObj<T>" with some
  dynamically allocated object of type T (i.e. created with the 'new' function).
  This function creates a class_handle instance to keep track of the object of
  type T and returns a pointer to an integer (meant to be used by matlab). This
  integer is a recasted pointer to the class_handle instance that manages the
  object of type T. It can thus be used as an identifier for this object (no two
  objects can have the same identifier).

  Recovering the object of type T corresponding to a given identifier at a later
  time is easily accomplished by calling the function recoverObjFromMatlabId with a
  pointer to the object's identifier. After checking the validity of the
  identifier, this function returns a reference to the corresponding object of
  type T. This object can then be easily manipulated.

  Finally, calling the function destroyObject, again with a pointer to the
  corresponding identifier, allows one to safely delete the object, together
  with its corresponding class_handle instance. Note that once an object has
  been assigned an identifier, it should only be deleted with this function.

  Example :
  ---------
    // We create a new object dynamically
    int& object = *(new int);
    object = 2;

    // Now we return the
    plhs[0] = createMatlabIdFromObj<T>(object);

    // Here we can exit the library and go back to matlab: our object remains in
    // memory and matlab has its identifier.

    ...

    // Now, matlab asks us to do something with our object. It sends us a
    // pointer to the object's identifier in prhs[1]. We start by reactivating
    // the object :
    int& object2 = recoverObjFromMatlabId<T>(prhs[1]);

    // Now we can apply operations on this object:
    object2 += 1;

    // Matlab still has a copy of the identifier, so we don't need to call
    // convertObjb2Mat anymore before handing back the hand to matlab. In fact,
    // we should not call converObj2Mat two times with the same object: this
    // would create several instances of the class_handle for the same object,
    // but only one can later clear the memory.

    ...

    // Now, matlab asks us to delete the our object. It sends us a pointer to
    // the object's identifier in prhs[1]. We call the dedicated delete
    // procedure to do the job. This deletes both our object, and the
    // class_handle instance which kept track of it.
    deleteObject(prhs[1]);


  More details on this file :
  ---------------------------
  The class used here to keep track of c++ object contains three fields:
    - a pointer
    - a signature
    - a name

  When constructed from an object of some type, a pointer to this object is
  copied internally and the two other variables are set as follows :
    - the signature is assigned a value proper to the current implementation
    - the type of the object the pointer points to is remembered in 'name_m'

  These informations allow one to safely reactivate a zone of memory containing
  such an object by just using a pointer to a class instance defined here.
  Indeed, if this pointer points to a zone of memory with a valid signature, and
  a name describing the expected variable type, then one can safely assume that
  the pointer within this class is also a valid one. This last point, however
  cannot be truly checked: for instance, it could be that the memory has been
  freed in the mean time. To avoid this problem, one should only pass here
  objects allocated dynamically, and on should use the dedicated function
  deleteObject to delete the object.

               ________________           _____________
               |              |           | an object |
               | - pointer ---|---------->|  of some  |
      pointer  |              |           |   type    |
     --------->| - signature  |           |___________|
               |              |          A zone of memory
               | - name       |       that can safely remain
               |______________|        unactive for a while
              class_handle object

*/

/* Note : In order to allow different mex files to create and destroy objects
   of the same type, we need to bypass the mexLock feature, and define our own
   counting procedure that counts individually all objects by type... So for
   now we just desactivate the mexLock and mexUnlock features. */


// Here we define the class used by the following functions
#define CLASS_HANDLE_SIGNATURE 0xB20F1D58
template<class T> class class_handle
{
public:
    class_handle(T * const ptr) : ptr_m(ptr), name_m(typeid(T).name()) { signature_m = CLASS_HANDLE_SIGNATURE; }
    ~class_handle() { signature_m = 0; delete ptr_m; }
    bool isValid() {//std::cout << "remembered type = " << name_m << ", calling type = " << typeid(T).name() << std::endl;
        return ((signature_m == CLASS_HANDLE_SIGNATURE) && !strcmp(name_m.c_str(), typeid(T).name())); }
    T *ptr() { return ptr_m; }

private:
    uint64_t signature_m;
    std::string name_m;
    T *ptr_m;
};


// This function creates a class_handle object pointing to the provided
// pointer (which points to an object of type T), and returns a matlab
// integer whose translation to a pointer points to this class_handle,
// which itself contains the provided pointer, i.e. points to the
// original object of type T
template<class T> inline mxArray *createMatlabIdFromPtr(T * const ptr)
{
    //mexLock();
    //mexPrintf("One call to mexLock() for an object of type %s\n", typeid(T).name());
    mxArray *out = mxCreateNumericMatrix(1, 1, mxUINT64_CLASS, mxREAL);
    *((uint64_t *)mxGetData(out)) = reinterpret_cast<uint64_t>(new class_handle<T>(ptr));
    return out;
}

// This function prepares an identifier corresponding to the object obj for
// matlab. The object given here as argument should continue to exist after
// exiting the c++ code and returning to matlab (i.e. is should have been
// allocated manually with 'new')
template<class T> inline mxArray *createMatlabIdFromObj(const T& obj)
{
    return createMatlabIdFromPtr(&obj);
}


// This function finds again the class_handle object created when the
// provided integer was created, and it checks that it is a valid object
// (if the values of signature_m corresponds to the global constant, and
// name_m corresponds to the type 'T' with which this function is
// called).
template<class T> inline class_handle<T> *convertMat2HandlePtr(const mxArray *in)
{
    //std::cout << "convertMat2HandlePtr stage" << std::endl;
    if (mxGetNumberOfElements(in) != 1 || mxGetClassID(in) != mxUINT64_CLASS || mxIsComplex(in))
        mexErrMsgTxt("Input must be a real uint64 scalar.");
    class_handle<T> *ptr = reinterpret_cast<class_handle<T> *>(*((uint64_t *)mxGetData(in)));
    //std::cout << "Validity : " << ptr->isValid() << std::endl;
    if (!ptr->isValid())
        mexErrMsgTxt("Handle not valid.");
    return ptr;
}

// This function returns the pointer contained in the class_handle object
// pointed to by the integer provided.
// In doing so, it re-activates the class_handle object created by the
// createMatlabIdFromPtr which produced the integer, and checks that this
// object is valid (In this way, an integer corresponding to nothing can
// be detected).
template<class T> inline T *recoverPtrFromMatlabId(const mxArray *in)
{
    return convertMat2HandlePtr<T>(in)->ptr();
}


// This functions returns a reference to the object of type T referred to by a
// matlab integer.
template<class T> inline T& recoverObjFromMatlabId(const mxArray *in)
{
    T* pointer;

    // We get the pointer to our object
    pointer = recoverPtrFromMatlabId<T>(in);

    // We place an object of type T at the place in memory pointed by the pointer
    T& result = *pointer;

    // ... and return that object
    return result;
}


// This function deletes and frees the memory of the class_handle object
// associated to a matlab integer. Deleting this class_handle object has the
//  consequence of also deleting the object of type 'T' that the class_handle
// points to.
// Before deletion, the validity of the provided integer as an identifier is
// verified.
template<class T> inline void destroyObject(const mxArray *in)
{
    //std::cout << "Deleting one object of type " << typeid(T).name() << " (destroyObject stage)" << std::endl;
    delete convertMat2HandlePtr<T>(in);
    //mexUnlock();
    //mexPrintf("One call to mexUnlock() for an object of type %s\n", typeid(T).name());
}


// This function just checks whether a matlab integer points to a valid
// object. It does so by constructing the associated class_handle object
// and checking its validity.
template<class T> inline bool checkValidity(const mxArray *in)
{
    if (mxGetNumberOfElements(in) != 1 || mxGetClassID(in) != mxUINT64_CLASS || mxIsComplex(in))
        mexErrMsgTxt("Input must be a real uint64 scalar.");
    class_handle<T> *ptr = reinterpret_cast<class_handle<T> *>(*((uint64_t *)mxGetData(in)));
    return ptr->isValid();
}


#endif // __CLASS_HANDLE_HPP__
