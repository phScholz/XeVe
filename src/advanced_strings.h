#ifndef ADVANCEDSTRING_H
#define ADVANCEDSTRING_H
extern char *as_trim(char *str);

// trimmed line starts with a '#' or '$'
extern int as_is_comment(char *str);

// checks if line only constists of spaces
extern int as_is_empty(char *str);

// replace $d
extern void as_replace_detector(char *dest, const char *src, int det);

// replace $d and #..#
extern void as_replace_wildcards(char *dest, const char *src, int run, int det);
#endif
