#ifndef _QGLVIEWERHEADERS_H_
#define _QGLVIEWERHEADERS_H_

#ifdef WIN32
#define _interlockedbittestandset fk_m$_set_QGLVIEWERHEADERS
#define _interlockedbittestandreset fk_m$_reset_QGLVIEWERHEADERS
#include <QGLViewer/qglviewer.h>
#undef _interlockedbittestandset
#undef _interlockedbittestandreset
#else
#include <QGLViewer/qglviewer.h>
#endif

#endif
