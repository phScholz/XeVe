/* This file defines the current dgf format. If Xia change their format 
 * again, this is the file to modify... It should NOT contain anything
 * which depends on parameter settings like the trace length, and should
 * only need to be modified if the dsp code is changed */

#ifndef _DGF_FORMAT
#define _DGF_FORMAT

/* Number of channels per module */

#define DGF_NCHAN 4            /* Number of channels in module */
#define DGF_BITS_IN_PATTERN 16 /* Number of bits in hit pattern */

/* Length of a DGF format timestamp (3 short integers = 6 bytes) */

#define DGF_TIMELEN 3

/* Maximum size of entire dgf buffer */

#define DGF_MAXBUF 8192

/* Minimumn and maximum valid DGF serial number */

#define DGF_MIN_SN 1100
#define DGF_MAX_SN 1499

/*.......................................................................*/
 
/* Positions of items within buffer header */

#define DGF_H_NWORDS 0
#define DGF_H_MODULE 1
#define DGF_H_FORMAT 2
#define DGF_H_TIMESTAMP 3

/* Length of buffer header */

#define DGF_BUFHEADLEN 6

/* Mask for format */

#define DGF_FORMAT_MASK 0x01FF

/*.......................................................................*/

/* Positions of items within event header */

#define DGF_E_PATTERN 0
#define DGF_E_TIMESTAMP 1

/* Length of event header */

#define DGF_EVENTHEADLEN 3

/*.......................................................................*/

/* Positions within channel header in long format */

#define DGF_C_LNWORDS       0
#define DGF_C_LTRIGGER_TIME 1
#define DGF_C_LENERGY       2
#define DGF_C_LXIAPSA       3
#define DGF_C_LUSERPSA      4
#define DGF_C_LGSLTHI       5
#define DGF_C_LGSLTMI       6
#define DGF_C_LGSLTLO       7
#define DGF_C_LREALA        8

/* Length of long channel header */

#define DGF_LCHANHEADLEN 9

/* Formats (after masking for long format) */

#define DGF_FORMAT_LONG1 0x0100
#define DGF_FORMAT_LONG2 0x0101

/* Positions within channel header in psa format */

#define DGF_C_PTRIGGER_TIME 0
#define DGF_C_PENERGY       1
#define DGF_C_PXIAPSA       2
#define DGF_C_PUSERPSA      3

/* Length of PSA channel header */

#define DGF_PCHANHEADLEN 4

/* Format (after masking for PSA format) */

#define DGF_FORMAT_PSA 0x0102

/* Positions within channel header in short format */

#define DGF_C_STRIGGER_TIME 0
#define DGF_C_SENERGY       1

/* Length of short channel header */

#define DGF_SCHANHEADLEN 2

/* Format (after masking for SHORT format) */

#define DGF_FORMAT_SHORT 0x0103

/* Positions within channel header in unix timestamp format */

#define DGF_C_UTIME_LO    0
#define DGF_C_UTIME_HI    1

/* Length of short channel header */

#define DGF_UCHANHEADLEN 2

/* Format (after masking for unix timestamp format) */

#define DGF_FORMAT_UNIX 0x01FF

#endif
