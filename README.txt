                My favorite double
               --------------------

... is a c++ double. And this library gives you access to them.

Well, ok, there is not much difference between matlab and c++ doubles, but this project shows how a c++ class can be seemlessly interfaced with matlab.

This project heavily relies on the library "C++ class interface" by Oliver Woodford. It does go a few steps further however, in terms of readiness to be adapted to a different class (it includes data exchange from and to matlab), and hopefully also in terms of documentation (unless you prefer to read the forum messages here: http://www.mathworks.com/matlabcentral/newsreader/view_thread/278243 , or unless you think that nothing explains c++ code better than itseld of course ,-)

So here is a simple explanation of the working of the interface, and of the role of each component. You can find below an explanation for the installation of the interface, its usage, and the steps that you can try in order to transpose this interface to bring your own c++ class to matlab.


Notes:
------
  - We implement complex numbers as two real tables (a bit like matlab), and not as a table of complex numbers. Both approaches have advantages and inconvenients...


Plans for the future:
---------------------

  - Extend the format to support sparse matrices




The global picture:
-------------------

The interface between matlab and c++ relies of four main components:

          |1        |2        |3
          |         |         |
   pure   |/________|________\|  pure
  matlab  |\   4    |        /|  C++
          |         |         |
          |         |         |

1. myDouble.m
   This file constitutes the matlab frontend of the library.
   It implements a matlab class which contains all the information in relation with the c++ class.
   Since the class defined here is a fully matlab component, it can be manipulated seemlessly in matlab. In particular, it leaves to matlab the job of managing which object should be kept in memory, and which ones can be cleared: when matlab judges that a class instances is not needed anymore, matlab calls its delete function, which is transmitted to the c++ interface to clear the corresponding memory.
   This file takes care of transmiting all calls from matlab onto the c++ file at level 2.

2. myDouble_mex.cpp
   This c++ file communicates with the level 1 through the matlab c++ interface, and with the c++ class at level 3. It translates the information received from matlab in a c-like format into calls to c++ classes, and back.

3. myClass.hpp
   This file describes the c++ class that we want to access from matlab. It is pure c++, with almost no relation to matlab. The only dependency appears in the data encoding and decoding, which needs to be done in c-style (we do it here rather than on level 2 in order to guaratee optimal performance).

4. class_handle.hpp
   This file is crucial to link all three interface levels 1, 2 and 3 above. Indeed, any object created in 2 or 3 (i.e. on the c++ side of the interface) is automatically deleted from memory when the execution returns to level 1 (matlab). Therefore, the approach used here to keep c++ data in memory between two calls to the interface is to allocate memory for every persistent object dynamically in 2 (and sometimes in 3 for optimal performance). When doing so, the objects stay in memory until they are manually deleted, or until matlab stops, so they can in principle be accessed again when needed.
   However, in order to reactivate an object that was created some time ago by the c++ library, matlab needs to provide some kind of identifier to the library, to let it know which object he is interested in reactivating. This should be a number (matlab only accepts numerical data) which was given to matlab by the library in the first place when the object was created. Establishing a correspondance between c++ class objects and integers identifiers is what the class_handle program takes care of. Note that this translation tool is implemented in quite a clever way. In particular, when called to restore a c++ object from its corresponding matlab identifier, this class_handle library performs some checks in order to make sure that the provided integer is a valid identifier.



Installation : Compile by running the command 'mex myDouble_mex.cpp' from within matlab.
--------------


Running : The following commands create two c++ numbers a and b, and adds them up
---------
  a = myDouble(2);
  b = myDouble(3+1i);
  c = a+b

Note : Since the main purpose of this c++ double library is to demonstrate a working interface between matlab and c++, only addition of two numbers is currently implemented. More operations can be easily added by following some of the steps explained below.



Adapting the interfact to your own c++ class :
----------------------------------------------

It should be relatively straightforward to transpose the current interface to your own favorite c++ class by following the following steps:
 - Replace myClass.hpp with the c++ class that you want to access in matlab
 - Make sure the methods that you want to be accessible from matlab are implemented in the following files :
    - myDouble.m
    - myDoule_mex.cpp

To see how to do that more precisly, check out the example of the plus operation in those files. It demonstrates several ways in which a class operation can be implemented. Note that it should not be necessary to change the code in class_handle.hpp.



Credits :
---------

This project is an adaptation by Jean-Daniel Bancal of the work of Oliver Woodford called "C++ class interface". It was finished on 12 April 2016.



