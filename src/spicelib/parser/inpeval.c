/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1985 Thomas L. Quarles
Modified: 2026- replaced hand-rolled mantissa accumulator with strtod() so
the parsed value is IEEE-754 correctly-rounded. The previous implementation
accumulated digits into a `double` integer mantissa, which lost the 16th
decimal digit at the edge of double precision and then multiplied by an
inexact `pow(10, expo)`, producing 1-ULP errors for values like
"0.0009999990000000001". Engineering suffix handling (T/G/K/U/N/P/F/M/Meg/Mil)
and FORTRAN-style D/d exponent markers are preserved.
**********/

#include "ngspice/ngspice.h"
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include "ngspice/inpdefs.h"
#include "inpxx.h"


double
INPevaluate(char **line, int *error, int gobble)
/* gobble: non-zero to gobble rest of token, zero to leave it alone */
{
    char *token;
    char *here;
    char *tmpline;
    char *endptr;
    char buf[256];
    double value;
    int expo_suffix;
    int len;

    /* setup */
    tmpline = *line;

    if (gobble) {
        /* MW. INPgetUTok should be called with gobble=0 or it make
         * errors in v(1,2) exp */
        *error = INPgetUTok(line, &token, 0);
        if (*error)
            return (0.0);
    } else {
        token = *line;
        *error = 0;
    }

    here = token;

    /* Validity check: after optional sign, must see a digit or '.' */
    {
        char *probe = here;
        if (*probe == '+' || *probe == '-')
            probe++;
        if ((*probe == '\0') ||
            ((!(isdigit((unsigned char)*probe))) && (*probe != '.'))) {
            /* number looks like just a sign! */
            *error = 1;
            if (gobble) {
                FREE(token);
                /* back out the 'gettok' operation */
                *line = tmpline;
            }
            return (0);
        }
    }

    /* Copy the numeric prefix into `buf`, translating FORTRAN-style
     * D/d exponent markers to e/E so strtod() can consume them.
     * Stops at the first non-numeric character (engineering suffix,
     * ':', whitespace, end of token, etc.).
     *
     * Original semantics preserved:
     *   - leading sign accepted once
     *   - at most one '.' (decimal point)
     *   - at most one exponent marker (e/E/d/D)
     *   - sign accepted only at start, or immediately after exponent marker
     *   - ':' immediately after the integer part terminates parsing,
     *     leaving ':' for the caller (matches ternary / subcircuit node use).
     */
    len = 0;
    {
        char *src = here;
        int saw_e = 0;
        int saw_dot = 0;
        int after_e = 0;

        if ((*src == '+' || *src == '-') && len < (int)sizeof(buf) - 1) {
            buf[len++] = *src;
            src++;
        }
        while (*src && len < (int)sizeof(buf) - 1) {
            char c = *src;
            if (isdigit((unsigned char)c)) {
                buf[len++] = c;
                after_e = 0;
            } else if (c == '.' && !saw_dot && !saw_e) {
                buf[len++] = c;
                saw_dot = 1;
            } else if ((c == 'e' || c == 'E') && !saw_e) {
                buf[len++] = c;
                saw_e = 1;
                after_e = 1;
            } else if ((c == 'd' || c == 'D') && !saw_e) {
                buf[len++] = (c == 'd') ? 'e' : 'E';
                saw_e = 1;
                after_e = 1;
            } else if ((c == '+' || c == '-') && after_e) {
                buf[len++] = c;
                after_e = 0;
            } else {
                break;
            }
            src++;
        }
        buf[len] = '\0';
    }

    value = strtod(buf, &endptr);
    /* Advance the original token pointer past the consumed numeric portion. */
    here += (endptr - buf);

    /* Engineering scale suffix (matches original switch). */
    expo_suffix = 0;
    switch (*here) {
    case 't': case 'T':
        expo_suffix = 12;
        here++;
        break;
    case 'g': case 'G':
        expo_suffix = 9;
        here++;
        break;
    case 'k': case 'K':
        expo_suffix = 3;
        here++;
        break;
    case 'u': case 'U':
        expo_suffix = -6;
        here++;
        break;
    case 'n': case 'N':
        expo_suffix = -9;
        here++;
        break;
    case 'p': case 'P':
        expo_suffix = -12;
        here++;
        break;
    case 'f': case 'F':
        expo_suffix = -15;
        here++;
        break;
    case 'a': case 'A':
        expo_suffix = -18;
        here++;
        break;
    case 'm': case 'M':
        if (((here[1] == 'E') || (here[1] == 'e')) &&
            ((here[2] == 'G') || (here[2] == 'g'))) {
            expo_suffix = 6;        /* Meg */
            here += 3;
        } else if (((here[1] == 'I') || (here[1] == 'i')) &&
                   ((here[2] == 'L') || (here[2] == 'l'))) {
            expo_suffix = -6;
            value *= 25.4;          /* Mil */
            here += 3;
        } else {
            expo_suffix = -3;       /* m, milli */
            here++;
        }
        break;
    default:
        break;
    }

    if (expo_suffix != 0) {
        value *= pow(10.0, (double) expo_suffix);
    }

    if (gobble) {
        FREE(token);
    } else {
        *line = here;
    }

    return value;
}


/* In addition to fcn INPevaluate() above, allow values like 4k7,
   similar to the RKM code (used by inp2r) */
double
INPevaluateRKM_R(char** line, int* error, int gobble)
/* gobble: non-zero to gobble rest of token, zero to leave it alone */
{
    char* token;
    char* here;
    double mantis;
    double deci;
    int expo1;
    int expo2;
    int expo3;
    int sign;
    int expsgn;
    char* tmpline;
    bool hasmulti = FALSE;

    /* setup */
    tmpline = *line;

    if (gobble) {
        /* MW. INPgetUTok should be called with gobble=0 or it leads to
         * errors in v(1,2) expression */
        *error = INPgetUTok(line, &token, 0);
        if (*error)
            return (0.0);
    }
    else {
        token = *line;
        *error = 0;
    }

    mantis = 0;
    deci = 0;
    expo1 = 0;
    expo2 = 0;
    expo3 = 0;
    sign = 1;
    expsgn = 1;

    /* loop through all of the input token */
    here = token;

    if (*here == '+')
        here++;                 /* plus, so do nothing except skip it */
    else if (*here == '-') {    /* minus, so skip it, and change sign */
        here++;
        sign = -1;
    }

    if ((*here == '\0') || ((!(isdigit_c(*here))) && (*here != '.') && (*here != 'r'))) {
        /* number looks like just a sign! */
        *error = 1;
        if (gobble) {
            FREE(token);
            /* back out the 'gettok' operation */
            *line = tmpline;
        }
        return (0);
    }

    while (isdigit_c(*here)) {
        /* digit, so accumulate it. */
        mantis = 10 * mantis + *here - '0';
        here++;
    }

    if (*here == '\0') {
        /* reached the end of token - done. */
        if (gobble) {
            FREE(token);
        }
        else {
            *line = here;
        }
        return ((double)mantis * sign);
    }

    if (*here == ':') {
        /* ':' is no longer used for subcircuit node numbering
           but is part of ternary function a?b:c
           FIXME : subcircuit models still use ':' for model numbering
           Will this hurt somewhere? */
        if (gobble) {
            FREE(token);
        }
        else {
            *line = here;
        }
        return ((double)mantis * sign);
    }

    /* after decimal point! */
    if (*here == '.') {
        /* found a decimal point! */
        here++;                 /* skip to next character */

        if (*here == '\0') {
            /* number ends in the decimal point */
            if (gobble) {
                FREE(token);
            }
            else {
                *line = here;
            }
            return ((double)mantis * sign);
        }

        while (isdigit_c(*here)) {
            /* digit, so accumulate it. */
            mantis = 10 * mantis + *here - '0';
            expo1 = expo1 - 1;
            here++;
        }
    }

    /* now look for "E","e",etc to indicate an exponent */
    if ((*here == 'E') || (*here == 'e') || (*here == 'D') || (*here == 'd')) {

        /* have an exponent, so skip the e */
        here++;

        /* now look for exponent sign */
        if (*here == '+')
            here++;             /* just skip + */
        else if (*here == '-') {
            here++;             /* skip over minus sign */
            expsgn = -1;        /* and make a negative exponent */
            /* now look for the digits of the exponent */
        }

        while (isdigit_c(*here)) {
            expo2 = 10 * expo2 + *here - '0';
            here++;
        }
    }

    /* now we have all of the numeric part of the number, time to
     * look for the scale factor (alphabetic)
     */
    switch (*here) {
    case 't':
    case 'T':
        expo1 = expo1 + 12;
        hasmulti = TRUE;
        break;
    case 'g':
    case 'G':
        expo1 = expo1 + 9;
        hasmulti = TRUE;
        break;
    case 'k':
    case 'K':
        expo1 = expo1 + 3;
        hasmulti = TRUE;
        break;
    case 'u':
    case 'U':
        expo1 = expo1 - 6;
        hasmulti = TRUE;
        break;
    case 'r':
    case 'R':
        /* This should be R150, i.e. R followed by an integer number */
        {
            int num;
            char ch;
            if (sscanf(here + 1, "%i%c", &num, &ch) == 1) {
                expo1 = expo1;
                hasmulti = TRUE;
            }
            else {
                *error = 1;
                if (gobble) {
                    FREE(token);
                    /* back out the 'gettok' operation */
                    *line = tmpline;
                }
                return (0);
            }
        }
        break;
    case 'n':
    case 'N':
        expo1 = expo1 - 9;
        hasmulti = TRUE;
        break;
    case 'p':
    case 'P':
        expo1 = expo1 - 12;
        hasmulti = TRUE;
        break;
    case 'm':
    case 'M':
        if (((here[1] == 'E') || (here[1] == 'e')) &&
            ((here[2] == 'G') || (here[2] == 'g')))
        {
            expo1 = expo1 + 6;  /* Meg */
            here += 2;
            hasmulti = TRUE;
        }
        else if (((here[1] == 'I') || (here[1] == 'i')) &&
            ((here[2] == 'L') || (here[2] == 'l')))
        {
            expo1 = expo1 - 6;
            mantis *= 25.4;     /* Mil */
        }
        else {
            expo1 = expo1 - 3;  /* m, M for milli */
            hasmulti = TRUE;
        }
        break;
    case 'l':
    case 'L':
        expo1 = expo1 - 3;  /* m, milli */
        hasmulti = TRUE;
        break;
    default:
        break;
    }

    /* read a digit after multiplier */
    if (hasmulti) {
        here++;
        while (isdigit_c(*here)) {
            deci = 10 * deci + *here - '0';
            expo3 = expo3 - 1;
            here++;
        }
        mantis = mantis + deci * pow(10.0, (double)expo3);
    }

    if (gobble) {
        FREE(token);
    }
    else {
        *line = here;
    }

    return (sign * mantis *
        pow(10.0, (double)(expo1 + expsgn * expo2)));
}

/* In addition to fcn INPevaluate() above, allow values like 4k7,
   similar to the RKM code (used by inp2r) */
double
INPevaluateRKM_C(char** line, int* error, int gobble)
/* gobble: non-zero to gobble rest of token, zero to leave it alone */
{
    char* token;
    char* here;
    double mantis;
    double deci;
    int expo1;
    int expo2;
    int expo3;
    int sign;
    int expsgn;
    char* tmpline;
    bool hasmulti = FALSE;

    /* setup */
    tmpline = *line;

    if (gobble) {
        /* MW. INPgetUTok should be called with gobble=0 or it make
         * errors in v(1,2) exp */
        *error = INPgetUTok(line, &token, 0);
        if (*error)
            return (0.0);
    }
    else {
        token = *line;
        *error = 0;
    }

    mantis = 0;
    deci = 0;
    expo1 = 0;
    expo2 = 0;
    expo3 = 0;
    sign = 1;
    expsgn = 1;

    /* loop through all of the input token */
    here = token;

    if (*here == '+')
        here++;                 /* plus, so do nothing except skip it */
    else if (*here == '-') {    /* minus, so skip it, and change sign */
        here++;
        sign = -1;
    }

    if ((*here == '\0') || ((!(isdigit_c(*here))) && (*here != '.') && (*here != 'r'))) {
        /* number looks like just a sign! */
        *error = 1;
        if (gobble) {
            FREE(token);
            /* back out the 'gettok' operation */
            *line = tmpline;
        }
        return (0);
    }

    while (isdigit_c(*here)) {
        /* digit, so accumulate it. */
        mantis = 10 * mantis + *here - '0';
        here++;
    }

    if (*here == '\0') {
        /* reached the end of token - done. */
        if (gobble) {
            FREE(token);
        }
        else {
            *line = here;
        }
        return ((double)mantis * sign);
    }

    if (*here == ':') {
        /* ':' is no longer used for subcircuit node numbering
           but is part of ternary function a?b:c
           FIXME : subcircuit models still use ':' for model numbering
           Will this hurt somewhere? */
        if (gobble) {
            FREE(token);
        }
        else {
            *line = here;
        }
        return ((double)mantis * sign);
    }

    /* after decimal point! */
    if (*here == '.') {
        /* found a decimal point! */
        here++;                 /* skip to next character */

        if (*here == '\0') {
            /* number ends in the decimal point */
            if (gobble) {
                FREE(token);
            }
            else {
                *line = here;
            }
            return ((double)mantis * sign);
        }

        while (isdigit_c(*here)) {
            /* digit, so accumulate it. */
            mantis = 10 * mantis + *here - '0';
            expo1 = expo1 - 1;
            here++;
        }
    }

    /* now look for "E","e",etc to indicate an exponent */
    if ((*here == 'E') || (*here == 'e') || (*here == 'D') || (*here == 'd')) {

        /* have an exponent, so skip the e */
        here++;

        /* now look for exponent sign */
        if (*here == '+')
            here++;             /* just skip + */
        else if (*here == '-') {
            here++;             /* skip over minus sign */
            expsgn = -1;        /* and make a negative exponent */
            /* now look for the digits of the exponent */
        }

        while (isdigit_c(*here)) {
            expo2 = 10 * expo2 + *here - '0';
            here++;
        }
    }

    /* now we have all of the numeric part of the number, time to
     * look for the scale factor (alphabetic)
     */
    switch (*here) {
    case 't':
    case 'T':
        expo1 = expo1 + 12;
        hasmulti = TRUE;
        break;
    case 'g':
    case 'G':
        expo1 = expo1 + 9;
        hasmulti = TRUE;
        break;
    case 'k':
    case 'K':
        expo1 = expo1 + 3;
        hasmulti = TRUE;
        break;
    case 'u':
    case 'U':
        expo1 = expo1 - 6;
        hasmulti = TRUE;
        break;
    case 'r':
    case 'R':

        expo1 = expo1;
        hasmulti = TRUE;
        break;
    case 'n':
    case 'N':
        expo1 = expo1 - 9;
        hasmulti = TRUE;
        break;
    case 'p':
    case 'P':
        expo1 = expo1 - 12;
        hasmulti = TRUE;
        break;
    case 'f':
    case 'F':
        expo1 = expo1 - 15;
        hasmulti = TRUE;
        break;
    case 'a':
    case 'A':
        expo1 = expo1 - 18;
        break;
    case 'm':
    case 'M':
        if (((here[1] == 'E') || (here[1] == 'e')) &&
            ((here[2] == 'G') || (here[2] == 'g')))
        {
            expo1 = expo1 + 6;  /* Meg */
            here += 2;
            hasmulti = TRUE;
        }
        else if (((here[1] == 'I') || (here[1] == 'i')) &&
            ((here[2] == 'L') || (here[2] == 'l')))
        {
            expo1 = expo1 - 6;
            mantis *= 25.4;     /* Mil */
        }
        else {
            expo1 = expo1 - 3;  /* Meg as well */
            hasmulti = TRUE;
        }
        break;
    case 'l':
    case 'L':
        expo1 = expo1 - 3;  /* m, milli */
        hasmulti = TRUE;
        break;
    default:
        break;
    }

    /* read a digit after multiplier */
    if (hasmulti) {
        here++;
        while (isdigit_c(*here)) {
            deci = 10 * deci + *here - '0';
            expo3 = expo3 - 1;
            here++;
        }
        mantis = mantis + deci * pow(10.0, (double)expo3);
    }

    if (gobble) {
        FREE(token);
    }
    else {
        *line = here;
    }

    return (sign * mantis *
        pow(10.0, (double)(expo1 + expsgn * expo2)));
}

/* In addition to fcn INPevaluate() above, allow values like 4k7,
   similar to the RKM code (used by inp2r) */
double
INPevaluateRKM_L(char** line, int* error, int gobble)
/* gobble: non-zero to gobble rest of token, zero to leave it alone */
{
    char* token;
    char* here;
    double mantis;
    double deci;
    int expo1;
    int expo2;
    int expo3;
    int sign;
    int expsgn;
    char* tmpline;
    bool hasmulti = FALSE;

    /* setup */
    tmpline = *line;

    if (gobble) {
        /* MW. INPgetUTok should be called with gobble=0 or it make
         * errors in v(1,2) exp */
        *error = INPgetUTok(line, &token, 0);
        if (*error)
            return (0.0);
    }
    else {
        token = *line;
        *error = 0;
    }

    mantis = 0;
    deci = 0;
    expo1 = 0;
    expo2 = 0;
    expo3 = 0;
    sign = 1;
    expsgn = 1;

    /* loop through all of the input token */
    here = token;

    if (*here == '+')
        here++;                 /* plus, so do nothing except skip it */
    else if (*here == '-') {    /* minus, so skip it, and change sign */
        here++;
        sign = -1;
    }

    if ((*here == '\0') || ((!(isdigit_c(*here))) && (*here != '.') && (*here != 'r'))) {
        /* number looks like just a sign! */
        *error = 1;
        if (gobble) {
            FREE(token);
            /* back out the 'gettok' operation */
            *line = tmpline;
        }
        return (0);
    }

    while (isdigit_c(*here)) {
        /* digit, so accumulate it. */
        mantis = 10 * mantis + *here - '0';
        here++;
    }

    if (*here == '\0') {
        /* reached the end of token - done. */
        if (gobble) {
            FREE(token);
        }
        else {
            *line = here;
        }
        return ((double)mantis * sign);
    }

    if (*here == ':') {
        /* ':' is no longer used for subcircuit node numbering
           but is part of ternary function a?b:c
           FIXME : subcircuit models still use ':' for model numbering
           Will this hurt somewhere? */
        if (gobble) {
            FREE(token);
        }
        else {
            *line = here;
        }
        return ((double)mantis * sign);
    }

    /* after decimal point! */
    if (*here == '.') {
        /* found a decimal point! */
        here++;                 /* skip to next character */

        if (*here == '\0') {
            /* number ends in the decimal point */
            if (gobble) {
                FREE(token);
            }
            else {
                *line = here;
            }
            return ((double)mantis * sign);
        }

        while (isdigit_c(*here)) {
            /* digit, so accumulate it. */
            mantis = 10 * mantis + *here - '0';
            expo1 = expo1 - 1;
            here++;
        }
    }

    /* now look for "E","e",etc to indicate an exponent */
    if ((*here == 'E') || (*here == 'e') || (*here == 'D') || (*here == 'd')) {

        /* have an exponent, so skip the e */
        here++;

        /* now look for exponent sign */
        if (*here == '+')
            here++;             /* just skip + */
        else if (*here == '-') {
            here++;             /* skip over minus sign */
            expsgn = -1;        /* and make a negative exponent */
            /* now look for the digits of the exponent */
        }

        while (isdigit_c(*here)) {
            expo2 = 10 * expo2 + *here - '0';
            here++;
        }
    }

    /* now we have all of the numeric part of the number, time to
     * look for the scale factor (alphabetic)
     */
    switch (*here) {
    case 't':
    case 'T':
        expo1 = expo1 + 12;
        hasmulti = TRUE;
        break;
    case 'g':
    case 'G':
        expo1 = expo1 + 9;
        hasmulti = TRUE;
        break;
    case 'k':
    case 'K':
        expo1 = expo1 + 3;
        hasmulti = TRUE;
        break;
    case 'u':
    case 'U':
        expo1 = expo1 - 6;
        hasmulti = TRUE;
        break;
    case 'r':
    case 'R':

        expo1 = expo1;
        hasmulti = TRUE;
        break;
    case 'n':
    case 'N':
        expo1 = expo1 - 9;
        hasmulti = TRUE;
        break;
    case 'p':
    case 'P':
        expo1 = expo1 - 12;
        hasmulti = TRUE;
        break;
    case 'f':
    case 'F':
        expo1 = expo1 - 15;
        hasmulti = TRUE;
        break;
    case 'a':
    case 'A':
        expo1 = expo1 - 18;
        break;
    case 'm':
    case 'M':
        if (((here[1] == 'E') || (here[1] == 'e')) &&
            ((here[2] == 'G') || (here[2] == 'g')))
        {
            expo1 = expo1 + 6;  /* Meg */
            here += 2;
            hasmulti = TRUE;
        }
        else if (((here[1] == 'I') || (here[1] == 'i')) &&
            ((here[2] == 'L') || (here[2] == 'l')))
        {
            expo1 = expo1 - 6;
            mantis *= 25.4;     /* Mil */
        }
        else {
            expo1 = expo1 - 3;  /* Meg as well */
            hasmulti = TRUE;
        }
        break;
    case 'l':
    case 'L':
        expo1 = expo1 - 3;  /* m, milli */
        hasmulti = TRUE;
        break;
    default:
        break;
    }

    /* read a digit after multiplier */
    if (hasmulti) {
        here++;
        while (isdigit_c(*here)) {
            deci = 10 * deci + *here - '0';
            expo3 = expo3 - 1;
            here++;
        }
        mantis = mantis + deci * pow(10.0, (double)expo3);
    }

    if (gobble) {
        FREE(token);
    }
    else {
        *line = here;
    }

    return (sign * mantis *
        pow(10.0, (double)(expo1 + expsgn * expo2)));
}

