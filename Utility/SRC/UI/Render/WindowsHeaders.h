#ifndef _WINDOWSHEADERS_H_
#define _WINDOWSHEADERS_H_

#ifdef WIN32
#define _interlockedbittestandset fk_m$_set_WINDOWSHEADERS
#define _interlockedbittestandreset fk_m$_reset_WINDOWSHEADERS
#include <windows.h>
#undef _interlockedbittestandset
#undef _interlockedbittestandreset
#endif

#endif
