typedef int IDL_LONG;

typedef int IDL_STRING_SLEN_T;
#define IDL_STRING_MAX_SLEN 2147483647


typedef struct {		/* Define string descriptor */
  IDL_STRING_SLEN_T slen;	/* Length of string, 0 for null */
  short stype;			/* type of string, static or dynamic */
  char *s;			/* Addr of string */
} IDL_STRING;

