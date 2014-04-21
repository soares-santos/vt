/*
 *	Copyright (c) 1999-2003 Smithsonian Astrophysical Observatory
 */

#include <launch.h>

#define LAUNCHARGS 1024

/* we one of these must be defined ... */
#if !defined(USE_PIPE) && !defined(USE_WAITPID)
#define USE_PIPE 1
#endif
/* ... but not both */
#if defined(USE_PIPE) && defined(USE_WAITPID)
#error "USE_PIPE and USE_WAITPID are mutually exclusive"
#endif

#ifdef USE_WAITPID
#define WAIT_TRIES  100
#define WAIT_MSEC  5000
#endif

#ifndef WAIT_MSEC
#define WAIT_MSEC  5000
#endif

static pid_t _launchpid=0;

/* spawnvp seems to be broken on cygwin as of 1/06, so just use fork/exec */
#if HAVE_CYGWIN
#define HAVE_CYGWIN_USE_SPAWNVP 0
#endif
/*
 *----------------------------------------------------------------------------
 *
 *
 * 			Private Routines and Data
 *
 *
 *----------------------------------------------------------------------------
 */

#if HAVE_CYGWIN||HAVE_MINGW32
#if HAVE_CYGWIN_USE_SPAWNVP

#ifdef ANSI_FUNC
static int launch_win32(char *cmdstring, int attach, char **stdfiles)
#else
static int launch_win32(cmdstring, attach, stdfiles)
     char *cmdstring;
     int attach;
     char **stdfiles;
#endif
{
  int i, j;
  int len;
  int got;
  int status;
  char *argv[LAUNCHARGS+1];
  char *path=NULL;
  char *s=NULL, *t=NULL;
  struct timeval tv;

  /* for now, we can't support stdfiles */
  if( stdfiles )
    return(-1);

  /* package up the arguments for new process */
  t = (char *)xstrdup(cmdstring);
  for(i=0, got=0, s=(char *)strtok(t, " \t"); s;
      i++, s=(char *)strtok(NULL," \t")){
    if( i < LAUNCHARGS ){ 
      /* save argument */
      argv[i] = xstrdup(s);
      /* change back special char to spaces, if necessary */
      len = strlen(argv[i]);
      for(j=0; j<len; j++){
	if( argv[i][j] == LAUNCH_SPACE){
	  argv[i][j] = ' ';
	}
      }
      argv[i+1] = NULL;
      /* save program name */
      if( i == 0 ) path = (char *)argv[i];
      got++;
    }
  }
  if( t ) xfree(t);
  if( attach )
    i = _P_WAIT;
  else
    i = _P_NOWAIT;
  if((status = spawnvp(i, path, (void *)argv)) != -1){
    status = 0;
    /* wait for child to start */
    tv.tv_sec = 0;
    tv.tv_usec = WAIT_MSEC;
    xselect(1, NULL, NULL, NULL, &tv);
  }
  for(i=0; i<got; i++)
    if( argv[i] ) xfree((char *)argv[i]);
  return(status);
}

#endif
#endif

/*
 *----------------------------------------------------------------------------
 *
 *
 * 			Public Routines and Data
 *
 *
 *----------------------------------------------------------------------------
 */

/*
 *
 * launchpid() -- return pid of last  launched process 
 *
 */
#ifdef ANSI_FUNC
pid_t launchpid(void)
#else
pid_t launchpid()
#endif
{
  return _launchpid;
}

#if HAVE_MINGW32==0

/* 
 *  adapted from the system() code in:
 *  W. Richard Stevens
 *  "Advanced Programming in the Unix Environment" 
 *  Addison-Wesley Publishing Co, 1992
 *  p. 314
 */
#ifdef ANSI_FUNC
int launch(char *cmdstring, int attach, char **stdfiles)
#else
int launch(cmdstring, attach, stdfiles)
     char *cmdstring;
     int attach;
     char **stdfiles;
#endif
{
  int			status;
  pid_t			pid;
  struct sigaction	ignore, saveintr, savequit;
  sigset_t		chldmask, savemask;
#ifdef USE_PIPE
  int			fd[2];
#endif

  /* return false if no command is specified */
  if( !cmdstring || !*cmdstring )
    return(-1);

  ignore.sa_handler = SIG_IGN;	/* ignore SIGINT and SIGQUIT */
  sigemptyset(&ignore.sa_mask);
  ignore.sa_flags = 0;
  if (sigaction(SIGINT, &ignore, &saveintr) < 0)
    return(-1);
  if (sigaction(SIGQUIT, &ignore, &savequit) < 0)
    return(-1);
  
  sigemptyset(&chldmask);	/* now block SIGCHLD */
  sigaddset(&chldmask, SIGCHLD);
  if (sigprocmask(SIG_BLOCK, &chldmask, &savemask) < 0)
    return(-1);
  
#if HAVE_CYGWIN_USE_SPAWNVP
  /* if we are on the Cygwin platform, use fork/exec only if we are
     redirecting stdfiles. Otherwise use spawnvp(), which works better. */
  if( stdfiles ){
#endif

#ifdef USE_PIPE
  /* open a pipe so parent can hear if the child fails to exec */
  if( !attach ){
    if( pipe(fd) < 0 )
      return(-1);
    xfcntl(fd[0], F_SETFD, FD_CLOEXEC);
    xfcntl(fd[1], F_SETFD, FD_CLOEXEC);
  }
#endif

  /* start new process */
  if( (pid = fork()) < 0 ){
#ifdef USE_PIPE
    if( !attach ){
      close(fd[0]);
      close(fd[1]);
    }
#endif
    status = -1;		/* ERROR: probably out of processes */

  } else if( pid == 0 ){	/* child */
    int i, j, len;
    char *argv[LAUNCHARGS+1];
    char *path=NULL;
    char *s=NULL, *t=NULL;

    /* close and reopen stdio files, if necessary */
    if( stdfiles ){
      for(i=0; i<3; i++){
	if( stdfiles[i] ){
	  close(i);
	  switch(i){
	  case 0:
	    if( open(stdfiles[i], O_RDONLY) < 0){
	      _exit(-1);
	    }
	    break;
	  case 1:
	    if( open(stdfiles[i], O_CREAT|O_WRONLY|O_TRUNC, 0600) < 0){
	      _exit(-1);
	    }
	    break;
	  case 2:
	    /* if stderr is the same as stdout, just dup */
	    if( stdfiles[1] && !strcmp(stdfiles[1], stdfiles[i]) ){
	      dup(1);
	    }
	    else{
	      if( open(stdfiles[i], O_CREAT|O_WRONLY|O_TRUNC, 0600) < 0){
		_exit(-1);
	      }
	    }
	    break;
	  }
	}
      }
    }

    /* restore previous signal actions & reset signal mask, but only if
       parent is waiting for completion (i.e., we are "attached") */
    if( attach ){
      sigaction(SIGINT, &saveintr, NULL);
      sigaction(SIGQUIT, &savequit, NULL);
      sigprocmask(SIG_SETMASK, &savemask, NULL);
    }
#ifdef USE_PIPE
    /* child closes reader -- only writes status */
    else{
      close(fd[0]);
    }
#endif

    /* package up the arguments for new process */
    t = (char *)xstrdup(cmdstring);
    for(i=0, s=(char *)strtok(t, " \t"); s;
	i++, s=(char *)strtok(NULL," \t")){
      if( i < LAUNCHARGS ){ 
	/* save argument */
	argv[i] = xstrdup(s);
	/* change back special char to spaces, if necessary */
	len = strlen(argv[i]);
	for(j=0; j<len; j++){
	  if( argv[i][j] == LAUNCH_SPACE){
	    argv[i][j] = ' ';
	  }
	}
	argv[i+1] = NULL;
	/* save program name */
	if( i == 0 ) path = argv[i];
      }
    }
    if( t ) xfree(t);
#ifndef HAVE_CYGWIN
    /* this call is broken in cygwin */
    /* for unattached processes, start a new session (with new process id),
       so that we do not inherit signals from parent (particularly SIGTERM) */
    if( !attach )
      setsid();
#endif
    /* start up the new program */
    if( execvp(path, argv) ){
      status = 127;
#ifdef USE_PIPE
      if( !attach ){
	write(fd[1], &status, 4);
	close(fd[1]);
      }
#endif
      _exit(status);		/* exec error */
    }
  } else {			/* parent */
    _launchpid = pid;
    /* wait for program termination from attached process */
    if( attach ){
      while( waitpid(pid, &status, 0) < 0 ){
	if( xerrno != EINTR ){
	  status = -1; /* error other than EINTR from waitpid() */
	  break;
	}
      }
    }
    else{
#ifdef USE_WAITPID
      int i, got;
      struct timeval tv;
      /* we wait up to WAIT_TRIES millisecs to make sure the child started;
         but if we get an error, we can exit immediately */
      for(i=0; i<WAIT_TRIES; i++){
	errno = 0;
	got=waitpid(pid, &status, WNOHANG);
	/* look for error termination */
	if( (got < 0) || ((got == 0) && xerrno) ){
	  got = -1;
	  /* make sure status shows error */
	  if( status == 0 )
	    status = -1;
	  break;
	}
	/* look for normal termination */
	else if( got > 0 ){
	  break;
	}
	/* no termination, sleep and wait some more */
	else{
	  tv.tv_sec = 0;
	  tv.tv_usec = WAIT_MSEC;
	  xselect(1, NULL, NULL, NULL, &tv);
	}
      }
      /* no termination means the child is still running */
      if( got == 0 )
	status = 0;
#endif
#ifdef USE_PIPE
      close(fd[1]);
      if( read(fd[0], &status, 4) == 0 ){
	status = 0;
      }
      close(fd[0]);
#endif
    }
  }
  
#if HAVE_CYGWIN_USE_SPAWNVP
  }
  /* for Cygwin, call their spawnvp() routine instead of fork()/exec() */
  else{
    status = launch_win32(cmdstring, attach, stdfiles);
  }
#endif

  /* restore previous signal actions & reset signal mask */
  if( sigaction(SIGINT, &saveintr, NULL) < 0 )
    return(-1);
  if( sigaction(SIGQUIT, &savequit, NULL) < 0 )
    return(-1);
  if( sigprocmask(SIG_SETMASK, &savemask, NULL) < 0 )
    return(-1);
  
  return(status);
}

#else

#ifdef ANSI_FUNC
int launch(char *cmdstring, int attach, char **stdfiles)
#else
int launch(cmdstring, attach, stdfiles)
     char *cmdstring;
     int attach;
     char **stdfiles;
#endif
{
  return launch_win32(cmdstring, attach, stdfiles);
}

#endif
