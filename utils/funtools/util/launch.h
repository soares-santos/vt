/*
 *	Copyright (c) 1999-2003 Smithsonian Astrophysical Observatory
 */

/*
 *
 * launch.h -- declarations for launching a program
 *
 */

#ifndef	__launch_h
#define	__launch_h

#if HAVE_CONFIG_H
#include <conf.h>
#endif

#include <stdio.h>
#include <errno.h>
#include <signal.h>
#include <fcntl.h>                                                       
#include <sys/time.h>
#include <sys/stat.h>
#if HAVE_STRING_H
#include <string.h>
#endif
#if HAVE_MALLOC_H
#include <malloc.h>
#endif
#if HAVE_UNISTD_H
#include <unistd.h>
#endif
#include <xport.h>
#include <word.h>
#include <xalloc.h>
#include <prsetup.h>

#define LAUNCH_SPACE '\001'

_PRbeg

int launch _PRx((char *cmdstring, int wait, char **stdfiles));
pid_t launchpid _PRx((void));

_PRend

#endif
