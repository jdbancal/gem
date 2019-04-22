#include "utils.hpp"

/*
  This file contains the implementation of the functions defined in utils.hpp.
  These are several utilities which are used by the GmpEigenMatrix library.
*/

/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


using namespace std;
using namespace mpfr;
using namespace Eigen;


// Here is a function which attempts to extract the n first digits of an mpreal
// number (with some rounding). If the precision of the number is smaller, it
// pads the rest with spaces. The function also takes a flag (modified by
// reference) that signals whether the rounding procedure created one more
// leading digit (e.g. if 99.999 is rounded to 100.0 at n=4). When this is the
// case, we know that if we had asked one more digit, it would have been a 0
// (i.e. the next digit in 100.0 as an approximation of 99.999 is 0, thus giving
// 100.00. We only get a non-zero digit if we go to one more digit.)
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
string firstDigitsToString(const mpreal& x, const int& n, bool& oneMoreDigitFlag)
{
    ostringstream format;

    if (n == 0)
    {
        return "";
    }

    // We determine how many digits we will actually request (not more than the
    // number's precision)
    int digits(n);
    if ((n < 0) || (n > bits2digits(x.get_prec())))
    {
        digits = 1+bits2digits(x.get_prec());
    }

    format << "%." << digits-1 << "RNe";

    // Now we can call the mpfr printing function
    // We need to pass through a c-style string to get our string...
    char *s = NULL;
    string out;

    // We ask the library for the requested digits
    if(!(mpfr_asprintf(&s, format.str().c_str(), x.mpfr_srcptr()) < 0))
    {
        out = string(s);

        mpfr_free_str(s);
    }

    // We determine whether the rounding process created one more leading digit.
    // Here is the exponent we would have without rounding:
    string initialExponent(iszero(x) ? "0" : ToString<double>(mpfr::floor(mpfr::log10(mpfr::abs(x))).toDouble()));
    // Here is the exponent that we found
    size_t exponentSign = out.find_last_of("+");
    if (exponentSign == string::npos)
        exponentSign = out.find_last_of("-")-1;
    string finalExponent(out.substr(exponentSign+1, out.length()-exponentSign+1));
    // If needed, we pad the initial exponent with zeros for faithful comparison
    if (initialExponent.length() < finalExponent.length()) {
        if (initialExponent.substr(0,1) == "-")
            initialExponent = "-" + nZeros(finalExponent.length()-initialExponent.length()) + initialExponent.substr(1,initialExponent.length()-1);
        else
            initialExponent = nZeros(finalExponent.length()-initialExponent.length()) + initialExponent;
    }
    // Now we can compare the two exponents
    if (finalExponent.compare(initialExponent) != 0)
        oneMoreDigitFlag = true;
    else
        oneMoreDigitFlag = false;

//    cout << endl << x << " : " << initialExponent << " versus " << out.substr((sgn(x)<0)+1+digits,out.length()-((sgn(x)<0)+1+digits)) << endl;
//    cout << endl << x << " : " << initialExponent << " versus " << finalExponent << endl;

//    cout << "The " << n << " extracted digits : " << out << endl;
    // Now extract the digits only from the string (i.e. remove the dot, but
    // keep sign
    out = out.substr((sgn(x)<0),1) + out.substr((sgn(x)<0)+2,digits-1);
//    cout << "and after removing the dot : " << out << endl;

    return out;
}


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
string mprealToString(const mpreal& x, int width, const int& exponentShift, int dotPosition, const bool& includeSign)
{

    // This function assumes that the dot is within the printing range (or just
    // out on the right hand side) (it also assumes other things like width > 0,
    // but we don't check for all possible misuse)
    if ((dotPosition < 0) || (dotPosition > width))
        cout << "Error : The dot is not in the printing range (!)" << endl;

    // This variable tells us whether we will be printing a dot (it is not the
    // case when the function is called with dotPosition = width, which we do
    // when printing integers)
    int dotInRange(((dotPosition >= 0) && (dotPosition < width)) ? 1 : 0);

//    cout << endl << x.toDouble() << endl;

    string number("");

    // We first set the sign
    if (includeSign)
    {
        if (sgn(x) < 0)
        {
            number = "-";
        } else {
            number = " ";
        }
        width = width-1;
        dotPosition = dotPosition-1;
    }

    // Let's find out how many digits do we need to compute for this number...
    // The display width is covering the range of powers of 10 between
    int displayPower1(exponentShift + (dotPosition-1));
    // and
    int displayPower2(displayPower1 - (width - dotInRange) + 1);

    // We are interested in the intersection between this range and the ones on
    // which the number has significant digits. Here are the exponents of the
    // first and last digits of our number
    int numberPower1(iszero(x) ? 0 : mpfr::floor(mpfr::log10(mpfr::abs(x))).toLong());
    int numberPower2(numberPower1 - bits2digits(x.get_prec()) + 1);

    // WARNING : If the first digit is rounded up, it may create one more digit higher up (!)
    // We should take care of that : 0.9 is

//    cout << "Printing on         10^[" << displayPower1 << ", " << displayPower2 << "]" << endl;
//    cout << "A number defined on 10^[" << numberPower1 << ", " << numberPower2 << "]" << endl;

    // First, we get rid of simple cases
    if (displayPower2 > numberPower1+1) { // +1 is because rounding off can potentially increase numberPower by one...
//        cout << "Case 1" << endl;
        number = number + nSpaces(max(0,dotPosition-1)) + nZeros(width-dotInRange-max(0,dotPosition-1));
    } else if (numberPower2 > displayPower1) {
        // Note : under expected usage, this case should never happen...
//        cout << "Case 2" << endl;
        number = number + nSpaces(width-dotInRange);
    } else {
//        cout << "Case 3" << endl;
        // Now we obtain the digits we need
        bool oneMoreDigitFlag(false);
//        cout << "We choose the minimum of " << numberPower1-numberPower2+1 << " and " << numberPower1-displayPower2+1 << endl;
        string someDigits(firstDigitsToString(x, min(numberPower1-numberPower2+1, numberPower1-displayPower2+1), oneMoreDigitFlag));
        if (oneMoreDigitFlag) {
            numberPower1 += 1;
            someDigits += "0"; // We add the following known digit
        }

        if (displayPower2 > numberPower1) {
            // In case the rounding process didn't bring a single digit back into
            // the printing field
            number = number + nSpaces(max(0,dotPosition-1)) + nZeros(width-dotInRange-max(0,dotPosition-1));
        } else {
//          cout << "First extracted digits are : '" << someDigits << "'" << endl;

            // We cut this string on the left hand side if needed
            if (numberPower1 > displayPower1){
//                cout << "Case 3.1" << endl;
                someDigits = someDigits.substr(numberPower1-displayPower1, someDigits.length()-(numberPower1-displayPower1));
            } else if (displayPower1 > numberPower1)
            {
//                cout << "Case 3.2" << endl;
                nSpaces(0);
//                cout << "dotPosition, displayPower1-numberPower1 = " << dotPosition << "," << displayPower1-numberPower1 << endl;
//                cout << "displayPower1-numberPower1 = " << displayPower1-numberPower1 << endl;
//                cout << "max(0,dotPosition-1) = " << max(0,dotPosition-1) << endl;
//                cout << "displayPower1-numberPower1-max(0,dotPosition-1) = " << displayPower1-numberPower1-max(0,dotPosition-1) << endl;
                // ... or maybe we complete this string on the left
                if (dotPosition > displayPower1-numberPower1) // Then we just need to add spaces
                    someDigits = nSpaces(displayPower1-numberPower1) + someDigits;
                    //someDigits = nSpaces(width-1-someDigits.length()) + someDigits; // should give the same
                else // Then we add zeros and maybe spaces
                    someDigits = nSpaces(max(0,dotPosition-1)) + nZeros(displayPower1-numberPower1-max(0,dotPosition-1)) + someDigits;
            }

//            cout << "After dealing with the left side : '" << someDigits << "'" << endl;

            // We complete the string on the right if necessary
            if (numberPower2 > displayPower2)
            {
                if (exponentShift == 0) {
                    if (numberPower2 < 0) // We have enough zeros already, we just add spaces
                        someDigits = someDigits + nSpaces(numberPower2-displayPower2);
                    else {
                        // We add zeros before the dot, and spaces afterwards
                        someDigits = someDigits + nZeros(numberPower2);
                        someDigits = someDigits + nSpaces(-displayPower2);
                    }
                } else
                    // We just add spaces
                    someDigits = someDigits + nSpaces(numberPower2-displayPower2);
            }

//            cout << "After dealing with the right side : '" << someDigits << "'" << endl;

            // Now we have our string of digits that we can add the the sign
            number = number + someDigits;
        }
    }

//    cout << "Before adding the dot we have the string '" << number << "'" << endl;
    // Finally, we add the dot where it should be
    if (dotInRange == 1)
        number = number.substr(0,includeSign+dotPosition) + "." + number.substr(includeSign+dotPosition,width-1);

    // Any space between the last digit and the dot should be turned into a zero

//cout << number.length() << " " << width << " " << includeSign << endl;
//    if (number.length() != width + includeSign)
//        cout << "Error : the output has a wrong number of digits(!)" << endl;

//    cout << " -> " << number << " which contains " << number.length() << " digits"<< endl;

    // Actually, we want the sign to be right before the first digit
    if ((includeSign) && (sgn(x) < 0)) {
        // We look for the last space in the beginning
        size_t firstDigit = number.find_first_not_of(" ", 1);
        number = nSpaces(firstDigit-1) + "-" + number.substr(firstDigit, number.length()-firstDigit);
    }

    return number;
}


// Here is another function which produces a formatted output. It produces a
// string of exactly width characters that best represent the number. The first
// character is reserved for a sign, all the rest is free. The procedure
// automatically switches to scientific notation if needed, so as to best
// represent the number.
// Warning : the behavior is not guaranteed when width is too small (e.g. < 3)
// WARNING : The behavior was not tested in presence of rounding ! (the exponent
//           could have changed... -> we need to check that)
string mprealToString2(const mpreal& x, const int& width, const bool& padOnTheRight)
{
    string number;

//    cout << endl << x.toDouble() << endl;

    // First, we check what is the exponent of the number
    int exponent(iszero(x) ? 0 : mpfr::floor(mpfr::log10(mpfr::abs(x))).toLong());
    int lengthExponent(ToString(exponent).length());

//    cout << exponent << " " << lengthExponent << endl;
//    cout << isint(x) << " " << iszero(x) << endl;

    // Based on this information, we deduce how many digits of the number we are
    // able to print
    int n(0);

    int lowestExponentForNonScientificNotation(-4);

    // Remember that two characters may be needed for the sign and the dot
    if ((exponent >= lowestExponentForNonScientificNotation) && (exponent < width-1))
    {
        // Then we don't need to use the scientific notation
        n = width - 2 + (exponent == width-2);

        // We get those digits
        bool oneMoreDigitFlag(false);
        string digits(firstDigitsToString(x, n, oneMoreDigitFlag));
        if (oneMoreDigitFlag) {
            // We update the exponent
            exponent += 1;
            lengthExponent = ToString(exponent).length();
            digits += "0";
        }
//        cout << "The first n digits are : " << digits << endl;
        // Add the sign
        if (x < 0)
        {
            digits = "-" + digits;
        } else {
            digits = " " + digits;
        }

        // We still need to add a dot if there are more digits than exponent-1
        if (digits.length()-1 > exponent+1)
            digits = digits.substr(0,exponent+2) + "." + digits.substr(exponent+2,digits.length()-(exponent+2));

        // Pad them if needed
        for (int i(0); digits.length() < width; ++i)
        {
            if (padOnTheRight)
            {
                digits = digits + " ";
            } else {
                digits = " " + digits;
            }
        }

        // The result is ready
        number = digits;
    } else {
        // Then we need the scientific notation
        n = width - 1 - 1 - 1 - lengthExponent;

        // We get those digits
        bool oneMoreDigitFlag(false);
        string digits(firstDigitsToString(x, n, oneMoreDigitFlag));
        if (oneMoreDigitFlag) {
            // We update the exponent
            exponent += 1;
            lengthExponent = ToString(exponent).length();
            digits += "0";
        }

        // Add the sign
        if (x < 0)
        {
            digits = "-" + digits;
        } else {
            digits = " " + digits;
        }

        // Add the dot
        digits = digits.substr(0,2) + "." + digits.substr(2,digits.length()-2);

        // Add the exponent
        digits = digits + "e" + ToString<int>(exponent);

        // Finally, pad the string if needed
        for (int i(0); digits.length() < width; ++i)
        {
            if (padOnTheRight)
            {
                digits = digits + " ";
            } else {
                digits = " " + digits;
            }
        }

        // The result is ready
        number = digits;
    }

    return number;
}
