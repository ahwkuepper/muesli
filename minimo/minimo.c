#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "getparam.h"

string *bindvec = NULL;  /* arr. of NULL term'd "name=value" strings */
int help_level = 0;
int debug_level = 0;
int error_level=0;
string yapp_string = NULL;
int yapp_dev = 0;
int review_flag = 0;
#define HELP_PROMPT 0x02
#define DIR_SEP  '/'

char namebuf[64];

void scan_environment(),error();

string parname(),parvalue(),parhelp(),scopy(),bindpar();
int xstrlen(),scanbind();
byte *allocate();
char *getmem(),*copxstr();
bool matchname();

void initparam(argv, defv)
string argv[];                   /* arg vector passed from main */
string defv[];                   /* default parameter values */
{
    int xstrlen(), ndef, i, j;
    byte *allocate();
    char *cp;
    bool posflag, useflag;
    string scopy(), tail();
    string help, name;
    string version_i="*";

    scan_environment();              /* scan env.var and set what is needed */

    ndef = xstrlen(defv, sizeof(string));       /* count params plus NULL   */
    bindvec = (string *) allocate((ndef + 1) * sizeof(string)); /* PPAP */
    for (i = 0; i <= ndef; i++)                 /* init bindings to NULL    */
        bindvec[i] = NULL;			/* note it's NULL terminated*/
    bindvec[0] = bindpar("argv0", tail(argv[0])); /* argv0 = progr name */

#if defined(HELPVEC)
    helpvec = (string *) allocate((ndef + 1) * sizeof(string)); /* PPAP */ 
    for (i = 0; i <= ndef; i++)                 /* init bindings to NULL */   
        helpvec[i] = NULL; 		        /* note it's NULL terminated */
    for (i=0; defv[i] != NULL; i++)
        helpvec[i] = parhelp(defv[i]);
#endif    
    for (i = 0; i < ndef-1; i++)                  /* loop over defaults */
        if (streq(parname(defv[i]),"VERSION")) {  /* until VERSION found */
            version_i = scopy(parvalue(defv[i]));   /* and store */
            break;
        }

    /*
     * the first thing to do is parse over the command line arguments
     * and install them into bindvec[], saves some work for later
     */

    help = NULL;                                /* turn help processing off */
    posflag = TRUE;                             /* start scan by position   */
    for (i = 1; argv[i] != NULL; i++) {         /* loop over stuff in argv  */
        name = parname(argv[i]);                /*   get param name, if any */
        posflag = posflag && (name == NULL);    /*   see how to match args  */
        if (posflag) {                          /*   match by position?     */
            if (i >= ndef)
                error("Too many args");
            bindvec[i] = bindpar(parname(defv[i-1]), argv[i]);
                                                /*     save name=val pair   */
        } else {                                /*   match by name?         */
            if (name == NULL)
                error("Arg \"%s\" must be named", argv[i]);
            j = 1 + scanbind(defv, name);       /*     find index in defv   */
            if (j > 0) {                        /*     was name found?      */
                if (bindvec[j] != NULL) {       /*       check slot is free */
                    error("Parameter \"%s\" duplicated", name);
		    exit(-1);
		}    
                bindvec[j] = argv[i];           /*       store name=value   */
            } else {                            /*     not listed in defv?  */
                if (streq(name, "help")) {      /*       got help keyword?  */
                    help = parvalue(argv[i]);
                    if ((cp = strpbrk(help,"0123456789"))!=NULL)  /* isnum? */
                        help_level = atoi(cp);  /* if so, change help_level */
                }
                else if (streq(name, "debug"))  /*       got debug keyword? */
                    debug_level = atoi(parvalue(argv[i]));
                else if (streq(name, "yapp")) { /*       got yapp keyword?  */
                    yapp_string = scopy(parvalue(argv[i]));
                    yapp_dev = atoi(yapp_string);
                }
                else if (streq(name, "review")) /*     go into review ?   */
                    review_flag = atoi(parvalue(argv[i]));
                else if (streq(name, "error")) /*     allowed error count? */
                    error_level = atoi(parvalue(argv[i]));
#if defined(REMOTE)
                else if (streq(name, "host"))   /*       got host keyword?  */
                    rsh(parvalue(argv[i]), argv);
#endif /* REMOTE */
                else {
                    printf("\nParameter \"%s\" unknown\n\n", name);
		    exit(-1);
		}    
            } /* j>0 (i.e. program/system keyword) */
        } /* if(posflag) */
    } /* for (i) */

#if defined(INTERRUPT)
    if (review_flag) {
        signal(SIGTERM, review);                /* catch interrupts */
        signal(SIGQUIT, review);                /* for review section */
        signal(SIGINT, review);                 /* and ^C also */
    }
#endif

    useflag =  FALSE;                         /* set if usage message needed */

#if defined(INTERACT)
    sprintf(&key_filename[strlen(key_filename)],"%s.def", tail(argv[0]));
    dprintf(2,"Keyword file will be %s\n",key_filename);

    if (help_level & HELP_DEFIO) {           /* read any existing keys-file */
	(void) readkeys("iniparam(1)",useflag,FALSE,defv);
    } /* if (help_level & HELP_DEFIO) */

    if (help_level & HELP_GLOBAL) {          /* read any existing keys-file */
        warning("HELP_GLOBAL level not implemented");
    }

#endif /* INTERACT */

    if (!(help_level & HELP_PROMPT)) {          /* if no help was needed */
        for (i=1; i<ndef; i++)                   /* loop installing defaults */
            if (bindvec[i]==NULL) {                /* empty slot to fill ? */
                useflag = useflag || streq(parvalue(defv[i-1]), "???");
                bindvec[i] = defv[i-1];            /* copy default entry */
            }
    } else if (help_level & HELP_PROMPT) {
#if defined(INTERACT)
        char  *cp, *scopy();
        string key, val;
        int    go=0;            /* set immediate go to false */
        int    v;
        
        if (!isatty(fileno(stdin)) || !isatty(fileno(stdout))) {
            printf("\n(initparam) Can only INTERACT with terminal"\n);
	    exit(-1);
	}    

        for (i=1; i<ndef; i++) {                 /* loop installing defaults */
            v = 0;
            do {                                 /*       loop for input     */
                if (bindvec[i]==NULL) {        /* had not yet been installed */
                    strcpy(keybuf,defv[i-1]);
                    if (streq(parname(keybuf), "VERSION")) {
                        v = 1;
                        version_i = scopy(parvalue("VERSION"));
                        key = parname(keybuf);
                        val = parvalue(keybuf);
                        break;
                    }
                } else {                            /* was already installed */
                    strcpy(keybuf,bindvec[i]);     
                    if (streq(parname(keybuf), "VERSION")) {
                        v = 1;
                        version_e = scopy(parvalue("VERSION"));
                        if (!streq(version_i,version_e)) 
                            warning("version_i=%s version_e=%s",version_i,version_e);
                    }
                }

                if (v) {             /* if keyword was "VERSION", break out here */
                    printf("%s: %s\n", argv[0], keybuf);
                    break;
                }
                if (streq(parvalue(keybuf), "???")) {  /* unanswered */
                    defflag = TRUE;
                    cp=strchr(keybuf,'=');
                    if (cp==NULL) {
                        printf("\n(initparam) Improper keyword (%s)\n",keybuf);
			exit(-1);
	            } else
                        cp++;
                    *cp = '\0';
                    go = readparam(keybuf, "No default allowed");  /* get it */
                } else {                  /* keyword with default value */
                    defflag = FALSE;
                    go = readparam(keybuf, "");     /* or use NULL for "" */
                }
                if (go) break;          /* exit loop when premature done */
                key = parname(keybuf);
                val = parvalue(keybuf);
            } while ((key == NULL) ||
                     ((val == NULL || streq(val, "")) && defflag) );
            bindvec[i] = bindpar(key, val);
            if (go) break;
        } /* for (i) */
       /************* GO still has error - not all par's are bount yet ********/
#else   
    printf("\n(initparam) INTERACTIVE input not allowed; check help=?\n");
    exit(-1);
#endif /* INTERACT */     
    }

#if defined(INTERACT)
    if (help_level&HELP_DEFIO || help_level&HELP_EDIT) {     /* write keyfile */
	writekeys("iniparam(1)");
    }
    if (help_level&HELP_EDIT) {    /* into editor mode, edit and read keyfile */
        char *getenv(), *edtcmd;
                                
        edtcmd = getenv("EDITOR");
        if (edtcmd != NULL)
            sprintf(keybuf,"%s %s",edtcmd,key_filename);       /* your editor */
        else
            sprintf(keybuf,"vi %s",key_filename);           /* default editor */

        if (system(keybuf))                         /*  edit the keyword file */
            error("(initparam) Error running editor (%s)",keybuf);
                        
        useflag = FALSE;                   /* will be set if usage msg needed */
	useflag=readkeys("iniparam(2)",useflag,TRUE,defv);

    }

    if (help_level&HELP_GLOBAL) {                /* write global keyword file */
        warning("HELP_GLOBAL: Not implemented");
    }
#endif /* INTERACT */

#if defined(USEDMM)
    if (help_level&HELP_DMM) {                          /* dmm user interface */
        int startfield = -1;
        int retcode, ncols, nrows, icol, irow;
        extern DMMRETTYPE       dmmRetType;
        DMMDRAWFLAG             drawflag = dmmNew;
        DMMVALFLAG              valflag = dmmInitVal;

        static DMMMESSAGE       msg1 =  {"Writing keyword file",
                                        5, 10, 0};
        static DMMMESGBLK       wmsgblk = {&msg1, 0};
        static DMMMENUITEM  runitem = {"Run", "Run the program"},
			    writeitem = {"Write", "Write the keyword file"},
			    quititem = {"Quit", "Quit without saving"};
        static DMMMENU 	    menu = { &runitem, &writeitem, &quititem, 0};
        DMMFIELD **screen, *fptr;
        /* some NEMO scope needs */
        extern string usage;
        char *scopy(), *sconc();
        byte *allocate();
        string *burststring();

        screen = (DMMFIELD **) allocate((ndef+1) * sizeof(DMMFIELD *));
        for (i=0; i<ndef; i++) {
            fptr =  (DMMFIELD *) allocate(sizeof(DMMFIELD));

            if (i==0) {          /* header with program name, version & usage */
                fptr->labelValue = sconc(sconc(parvalue(bindvec[0]),"  V"),
                            sconc(parvalue(bindvec[ndef-1]),sconc(" :",usage)));
                fptr->labelCol          = 5;
                fptr->inputCol          = 0;
                fptr->labelAttrib       = A_REVERSE|A_BLINK;
                fptr->inputEditValue    = 0;
                fptr->inputLength       = 0;
                fptr->inputInstructVec  = 0;
            } else {           /* Editable fields with keyword and def values */
                if(bindvec[i]==NULL) {
                    fptr->labelValue        = scopy(parname(defv[i-1]));
                    fptr->inputInitValue    = parvalue(defv[i-1]);
                } else {
                    fptr->labelValue        = scopy(parname(bindvec[i]));
                    fptr->inputInitValue    = parvalue(bindvec[i]);
                }
                fptr->labelCol          = 0;
                fptr->inputCol          = 12;
                fptr->labelAttrib       = 0;
                fptr->inputEditValue    = allocate(80);
                fptr->inputLength       = 60;      /* should be 80-12 ??*/
                fptr->inputInstructVec  = burststring(helpvec[i-1],"\n");
            }
            fptr->labelLine         = i;
            fptr->inputLine         = i;
            fptr->inputMaxCols      = 0;
            fptr->inputFullLines    = 0;
            fptr->inputShortCols    = 0;
            fptr->inputAttrib       = A_REVERSE;
            fptr->inputEditAttrib   = A_UNDERLINE;
            fptr->inputNumInstructs = 0;
            fptr->inputAudit        = NULL;
            screen[i] = fptr;
        }
        screen[ndef] = NULL;          /* terminate DMMFIELD array with a NULL */

        dmmInit();                                           /* set up screen */
        for(;;) {                                         /* loop until happy */
            retcode = dmmRun(startfield, screen, menu, drawflag, valflag);
            if(retcode==0) {                                   /* run program */
                       /* should actually check that keywords are still '???' */
                break;
            } else if(retcode==1) {                     /* write keyword file */
                dmmClear();
                dmmPutMsgsHit(wmsgblk);
            } else if(retcode==2)                       /* unconditional quit */
                break;            
            else
                break;
        } /* for(;;) */
        dmmClose();                              /* return to a normal screen */

        if(retcode == 2) {                                      /* catch Quit */
            exit(0);
            /*NOTREACHED*/
	}
        else if(retcode == -1)                /* and other bizarre situations */
            error("bad return code from dmmRun using help=8");

        for(i=1; i<ndef; i++) {            /* save new values back in bindvec */
            bindvec[i]=bindpar(screen[i]->labelValue,screen[i]->inputEditValue);
	    dprintf(1,"Bound par %d as %s\n",i,bindvec[i]);
        }
    }
#endif /* USEDMM */
} /* initparam */


string getparam(name)
string name;
{
    int i;
    byte *allocate();

    if (bindvec == NULL)
        error("(getparam) called before initparam");
    i = scanbind(bindvec, name);
    if (i < 0)
        error("(getparam) \"%s\" unknown", name);
#if defined(MACROREAD)
    cp = parvalue(bindvec[i]);
    if (*cp == '@') {
        cp = get_macro(cp);         /* get the macro in memory  */
        cp1 = (char *) allocate(strlen(cp)+strlen(name)+1);
        strcpy(cp1,name);
        strcat(cp1,"=");
        strcat(cp1,cp);
        bindvec[i]=cp1;
        free(cp);
        cp = strchr(cp1,'=');       /* look again for '=' - is always there */
        cp++;                       /* advance one to point to value part */
    } 
    return cp;                     /* return value part */
#else
    return parvalue(bindvec[i]);
#endif
}


/*
 * GETIPARAM, ..., GETDPARAM: get int, long, bool, or double parameters.
 *          and allow some optional parsing if -DNEMOINP
 */

int getiparam(par)
string par;                     /* name of parameter to look for */
{
    string val;
    int    atoi();

    val = getparam(par);                        /* obtain value of param */
#if !defined(NEMOINP)
    return (atoi(val));                         /* convert to an integer */
#else
    nret = nemoinpi(val,&ipar, 1);
    if (nret < 0)
        warning("getiparam(%s=%s) parsing error %d, assumed %d\n",
                    par,val,nret,ipar);
    if (nret==0)
        return(0);
    else
        return(ipar);
#endif /* !NEMOINP */
}

long getlparam(par)
string par;                     /* name of parameter to look for */
{
    string val;
    long atol();

    val = getparam(par);                        /* obtain value of param */
#if !defined(NEMOINP)
    return (atol(val));                         /* convert to an integer */
#else
    nret = nemoinpl(val,&lpar,1);
    if (nret < 0)
        warning("getlparam(%s=%s) parsing error %d assumed %l\n",
                    par,val,nret,lpar);
    if (nret==0)
        return(0);
    else
        return(lpar);
#endif /* !NEMOINP */
}

bool getbparam(par)
string par;                     /* name of parameter to look for */
{
    char *val;

    val = getparam(par);                        /* obtain value of param */
#if !defined(NEMOINP)
    if (strchr("tTyY1", *val) != NULL)          /* is value true? */
        return (TRUE);
    if (strchr("fFnN0", *val) != NULL)          /* is value false? */
        return (FALSE);
    error("getbparam: %s=%s not bool", par, val);
    return(0);   /*turboc*/
#else
    nret = nemoinpb(val,&bpar,1);
    if (nret < 0)
        warning("getbparam(%s=%s) parsing error %d, assumed %d\n",
                        par,val,nret,bpar);
    if (nret==0)
        return(FALSE);
    else
        return(bpar);
#endif /* !NEMOINP */
}


double getdparam(par)
string par;                     /* name of parameter to look for */
{
    string val;
    double atof();

    val = getparam(par);                        /* obtain value of param */
    if (val == NULL)                            /* insure value is given */
        error("getdparam: %s unknown", par);
#if !defined(NEMOINP)
    return (atof(val));                         /* convert to a double */
#else
    nret = nemoinpd(val,&dpar,1);
    if (nret < 0)
        warning("getdparam(%s=%s) parsing error %d, assumed %g\n",
                            par,val,nret,dpar);
    if (nret==0)
        return(0.0);
    else
        return(dpar);
#endif /* !NEMOINP */
}

void scan_environment()
{
    int i;
    int nemobin;
    char *getenv(), *date_id(), *scopy();
    extern char **environ;      /* (extern) global database for environment */

    for (i = 0; environ[i] != NULL; i++) {
        if (streq("DEBUG", parname(environ[i])))
            debug_level = atoi(parvalue(environ[i]));
        else if (streq("YAPP", parname(environ[i]))) {
            yapp_string = scopy(parvalue(environ[i]));
            yapp_dev = atoi(yapp_string);
        } else if (streq("HELP", parname(environ[i])))
            help_level = atoi(parvalue(environ[i]));
        else if (streq("REVIEW", parname(environ[i])))
            review_flag = atoi(parvalue(environ[i]));
        else if (streq("ERROR", parname(environ[i])))
            error_level = atoi(parvalue(environ[i]));
	else if (streq("NEMOBIN", parname(environ[i])))
            nemobin = 1;
    }
}


string parname(arg)
string arg;                     /* string of the form "name=value" */
{
    char *ap, *np;

    ap = (char *) arg;
    while (*ap == ' ')          /* skip initial whitespace */
        ap++;
    np = &namebuf[0];               /* point to spot where to copy to */
    while ((*np = *ap) != NULL) {       /* copy until '=' or ' ' */
#if 0
        if (*np == '=' || *np == ' ') {
#else
        if (*np == '=') {
#endif
            *np = NULL;
            return ((string) namebuf);
        }
        np++;
        ap++;
    }
    return (NULL);
}

string parvalue(arg)
string arg;                     /* string of the form "name=value" */
{
    char *ap;

    ap = (char *) arg;
    while (*ap != NULL)
        if (*ap++ == '=') {
            while (*ap == ' ')      /* skip whitespace after '=' */
                ap++;
            return ((string) ap);
        }
    return (NULL);
}

string parhelp(arg)
string arg;
{
    char *cp;

    cp = arg;    
    while (*cp != NULL && *cp != '\n')      /* scan until EOS or BOH */
        cp++;
    if (*cp == '\n')
        *cp = '\0';   /* patch BOH with EOS - DANGER with constant strings */
                      /* -fwriteable-strings is required in GNU compiler */
    cp++;
    while (*cp == ' ' || *cp == '\t')       /* skip white spaces */
        cp++;
    return ( (string) cp);      /* return NULL (no help) or valid help */
}

string scopy(s)
string s;
{
    string result;

    result = (string) getmem(strlen(s)+1);
    strcpy(result, s);
    return (result);
}

int xstrlen(xspt, nbyt)         /* returns count of values (w. NULL) */
char *xspt;                     /* ptr to ext. string to copy */
int nbyt;                       /* bytes per value */
{
    int nval, i;
    bool lpflg;

    nval = 0;                                   /* init count of values */
    do {                                        /* loop over values */
        nval++;                                 /*   count one more */
        lpflg = FALSE;                          /*   init loop flag */
        for (i = 0; i < nbyt; i++)              /*   loop over bytes */
            if (*xspt++ != NULL)                /*     a byte of data? */
                lpflg = TRUE;                   /*       set loop flag */
    } while (lpflg);                            /* until a NULL value */
    return (nval);                              /* return total count */
}


int scanbind(bvec, name)
string bvec[];  /* NULL terminated array of "key=val" strings */
string name;    /* 'key' to look for */
{
    int i;

    for (i = 0; bvec[i] != NULL; i++)
        if (matchname(bvec[i], name)) return i;
    return -1;
}


string bindpar(name, value)
string name;
string value;
{
    char bindbuf[128], *copxstr(), *scopy();

    if (strlen(name)+strlen(value)+1 > 128) 
        error("bindpar: not enough space for key=val");

    sprintf(bindbuf, "%s=%s", name, value);
    return (copxstr(bindbuf, sizeof(char)));
}

byte *allocate(nb)
int nb;
{
    byte *mem;

    mem = (byte *) malloc((unsigned)nb);
    if (mem == NULL)
	error("allocate: not enuf memory (%d bytes)", nb);
    return (mem);
}


string tail(filename)
string filename;
{
    char *slashpos;

    slashpos = strrchr(filename, DIR_SEP);
    if (slashpos == NULL)
	return (scopy(filename));
    else
        return (scopy(slashpos + 1));
}

char *callocate(nb)
int nb;
{
    char *mem;

    mem = (char *) malloc((unsigned)nb);
    if (mem == NULL)
        error("callocate: not enuf memory (%d bytes)", nb);
    return (mem);
}

char *getmem(nbytes)
int nbytes;
{
    char *callocate();

    return ( callocate(nbytes) );
}

bool matchname(bind, name)
string bind, name;
{
    char *bp, *np;
    
    if (name==NULL) return FALSE;

    bp = bind;
    np = name;
    while (*bp == ' ')          /* skip initial blanks */
        bp++;
    while (*np == ' ')          /* skip initial blanks */
        np++;
    while (*bp == *np) {
        bp++;
        np++;
    }
    return (*bp == '=' || *bp == ' ') && *np == NULL;
}

char *copxstr(xspt, nbyt)       /* returns ptr to new copy */
char *xspt;                     /* ptr to ext. string to copy */
int nbyt;                       /* bytes per value */
{
    int xstrlen(), n;
    char *dest, *dp, *callocate();

    n = nbyt * xstrlen(xspt, nbyt);             /* get length in bytes */
    dp = dest = callocate(n);                   /* allocate new storage */
    while (--n >= 0)                            /* loop over bytes */
        *dp++ = *xspt++;                        /*   copy each in turn */
    return (dest);                              /* return copy string */
}

void error(string s) {
  printf("\n%s\n\n",s);
  exit(-1);
}
