/*
 * GETPARAM.H: command-line processing definitions
 *
 * Mar. 1987  Joshua Edward Barnes, Princeton NJ.
 * Sep. 1990  -- added finiparam() to list              PJT
 * Nov. 1991  -- added hasvalue()                       PJT
 */

typedef char *string;
typedef char bool;
typedef unsigned char byte;

#define TRUE  (1)
#define FALSE (0)

#define streq(x,y) (strcmp((x), (y)) == 0)

void    initparam(/* argv, defv */);
void    finiparam(/* */);
void    stop(/* lev */);
void	putparam(/* key, value */);
void	promptparam(/* par, prompt */);
string  getparam(/* name */);
int     getiparam(/* name */);
long    getlparam(/* name */);
bool    getbparam(/* name */);
double  getdparam(/* name */);

/*
 * Macro used to obtain name of program.
 */

#define getargv0()      (getparam("argv0"))
