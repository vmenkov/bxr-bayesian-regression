/*
SL comment the code for Linux version;
 */

//#include <windows.h>


/*

typedef struct _SYSTEMTIME {  WORD wYear;  WORD wMonth;  WORD wDayOfWeek;  WORD wDay;  WORD wHour;  WORD wMinute;  WORD wSecond;  WORD wMilliseconds;
} SYSTEMTIME,  *PSYSTEMTIME;
*/


double TimeMilliseconds() {
/*	SYSTEMTIME st;
	GetSystemTime(&st);
	double t = (double)st.wMilliseconds + 1000.0*(double)st.wSecond + 60000.0*(double)st.wMinute + 3600000.0*(double)st.wHour;
	return t;
*/
    return 0;
}

