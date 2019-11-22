/*                                                     mconf.h
 *
 *     Common include file for math routines
 *
 *
 *
 * SYNOPSIS:
 *
 * #include "mconf.h"
 *
 *
 *
 * DESCRIPTION:
 *
 * This file contains definitions for error codes that are
 * passed to the common error handling routine mtherr()
 * (which see).
 *
 * The file also includes a conditional assembly definition
 * for the type of computer arithmetic (IEEE, DEC, Motorola
 * IEEE, or UNKnown).
 * 
 * For Digital Equipment PDP-11 and VAX computers, certain
 * IBM systems, and others that use numbers with a 56-bit
 * significand, the symbol DEC should be defined.  In this
 * mode, most floating point constants are given as arrays
 * of octal integers to eliminate decimal to binary conversion
 * errors that might be introduced by the compiler.
 *
 * For little-endian computers, such as IBM PC, that follow the
 * IEEE Standard for Binary Floating Point Arithmetic (ANSI/IEEE
 * Std 754-1985), the symbol IBMPC should be defined.  These
 * numbers have 53-bit significands.  In this mode, constants
 * are provided as arrays of hexadecimal 16 bit integers.
 *
 * Big-endian IEEE format is denoted MIEEE.  On some RISC
 * systems such as Sun SPARC, double precision constants
 * must be stored on 8-byte address boundaries.  Since integer
 * arrays may be aligned differently, the MIEEE configuration
 * may fail on such machines.
 *
 * To accommodate other types of computer arithmetic, all
 * constants are also provided in a normal decimal radix
 * which one can hope are correctly converted to a suitable
 * format by the available C language compiler.  To invoke
 * this mode, define the symbol UNK.
 *
 * An important difference among these modes is a predefined
 * set of machine arithmetic constants for each.  The numbers
 * MACHEP (the machine roundoff error), MAXNUM (largest number
 * represented), and several other parameters are preset by
 * the configuration symbol.  Check the file const.c to
 * ensure that these values are correct for your computer.
 *
 * Configurations NANS, INFINITIES, MINUSZERO, and DENORMAL
 * may fail on many systems.  Verify that they are supposed
 * to work on your computer.
 */

/*
 * Cephes Math Library Release 2.3:  June, 1995
 * Copyright 1984, 1987, 1989, 1995 by Stephen L. Moshier
 */


#ifndef CEPHES_MCONF_H
#define CEPHES_MCONF_H


#if defined(_MSC_VER)
        #define NPY_INLINE __inline
#elif defined(__GNUC__)
	#if defined(__STRICT_ANSI__)
		#define NPY_INLINE __inline__
	#else
		#define NPY_INLINE inline
	#endif
#else
        #define NPY_INLINE
#endif

        
#include "protos.h"
#include "mex.h"

/* Constant definitions for math error conditions
 */

#define DOMAIN		1	/* argument domain error */
#define SING		2	/* argument singularity */
#define OVERFLOW	3	/* overflow range error */
#define UNDERFLOW	4	/* underflow range error */
#define TLOSS		5	/* total loss of precision */
#define PLOSS		6	/* partial loss of precision */
#define TOOMANY         7	/* too many iterations */
#define MAXITER        500

#define EDOM		33
#define ERANGE		34

NPY_INLINE static float __npy_inff(void)
{
    const union { unsigned long __i; float __f;} __bint = {0x7f800000UL};
    return __bint.__f;
}

NPY_INLINE static float __npy_nanf(void)
{
    const union { unsigned long __i; float __f;} __bint = {0x7fc00000UL};
    return __bint.__f;
}

/*#define npy_isfinite(x) _finite((x))*/
#define NPY_NAN ((double)__npy_nanf())
#define NPY_INFINITY ((double)__npy_inff())
#define NPY_PI 3.141592653589793238462643383279502884
#define NPY_EULER 0.577215664901532860606512090082402431


#ifdef _MSC_VER
    #define npy_isfinite(x) _finite((x))
#else
    #include <math.h>
    #define npy_isfinite(x) isfinite((x))
#endif 

/* Long double complex numeral.  */
/*
 * typedef struct
 * {
 * long double r;
 * long double i;
 * } cmplxl;
 */

/* Type of computer arithmetic */

/* UNKnown arithmetic, invokes coefficients given in
 * normal decimal format.  Beware of range boundary
 * problems (MACHEP, MAXLOG, etc. in const.c) and
 * roundoff problems in pow.c:
 * (Sun SPARCstation)
 */

/* SciPy note: by defining UNK, we prevent the compiler from
 * casting integers to floating point numbers.  If the Endianness
 * is detected incorrectly, this causes problems on some platforms.
 */
#define UNK 1

/* Define this `volatile' if your compiler thinks
 * that floating point arithmetic obeys the associative
 * and distributive laws.  It will defeat some optimizations
 * (but probably not enough of them).
 *
 * #define VOLATILE volatile
 */
#define VOLATILE

/* For 12-byte long doubles on an i386, pad a 16-bit short 0
 * to the end of real constants initialized by integer arrays.
 *
 * #define XPD 0,
 *
 * Otherwise, the type is 10 bytes long and XPD should be
 * defined blank (e.g., Microsoft C).
 *
 * #define XPD
 */
#define XPD 0,

/* Define to support tiny denormal numbers, else undefine. */
#define DENORMAL 1

/* Define to distinguish between -0.0 and +0.0.  */
#define MINUSZERO 1

/* Define 1 for ANSI C atan2() function
 * See atan.c and clog.c. */
#define ANSIC 1

/* Variable for error reporting.  See mtherr.c.  */
/* extern int merror;*/

#define gamma Gamma

#endif				/* CEPHES_MCONF_H */
