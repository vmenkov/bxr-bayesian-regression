// This header includes code to handle compiler dependencies.

#ifndef _COMPILER_H
#define _COMPILER_H


#ifndef _MSC_VER
#define stricmp(a, b)   strcasecmp(a,b)
#else
#define finite(x) (_finite(x))
#define _stricmp stricmp
#define sscanf_s sscanf
#define _getch getch
#endif


// TODO: Move to util.h
bool readFromStdin(const char* fname); 

typedef unsigned YType;


#endif // _COMPILER_H

