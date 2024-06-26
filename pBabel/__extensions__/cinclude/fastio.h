/***************************************************************************
 *cr
 *cr            (C) Copyright 1995-2009 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ***************************************************************************/
/***************************************************************************
 * RCS INFORMATION:
 *
 *      $RCSfile: fastio.h,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.25 $       $Date: 2011/10/14 03:52:46 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *   This is a simple abstraction layer for system-dependent I/O calls
 * that allow plugins to do binary I/O using the fastest possible method.
 *
 * This code is intended for use by binary trajectory reader plugins that
 * work with multi-gigabyte data sets, reading only binary data.
 *
 ***************************************************************************/

/*--------------------------------------------------------------------------------------------------------------------------------*/
/* . MJF: included from "largefiles.h" as need for O_DIRECT. */
#if defined(_AIX)

/* Define to enable large file extensions on AIX */
#define _LARGE_FILE
#define _LARGE_FILES

#else

/* Defines which enable LFS I/O interfaces for large (>2GB) file support
 * on 32-bit machines.  These must be defined before inclusion of any
 * system headers.
 */

#ifndef _LARGEFILE_SOURCE
#define _LARGEFILE_SOURCE
#endif

#define _FILE_OFFSET_BITS 64
#endif /* . MJF _AIX. */
/*--------------------------------------------------------------------------------------------------------------------------------*/

#ifndef FASTIO_H
#define FASTIO_H
#define FIO_READ    0x01
#define FIO_WRITE   0x02
#define FIO_DIRECT  0x04 /* emulate Unix O_DIRECT flag */
 
/*==================================================================================================================================
! . Windows.
!=================================================================================================================================*/
#if defined(_MSC_VER)

/*--------------------------------------------------------------------------------------------------------------------------------*/
#if 1 /* . Native Windows I/O calls. */
/*--------------------------------------------------------------------------------------------------------------------------------*/
#define FASTIO_NATIVEWIN32 1

#include <stdio.h>
#include <string.h>
#include <windows.h>

typedef HANDLE fio_fd;
typedef LONGLONG fio_size_t;
typedef void * fio_caddr_t;

typedef struct {
  fio_caddr_t iov_base;
  int iov_len;
} fio_iovec;


#define FIO_SEEK_CUR  FILE_CURRENT
#define FIO_SEEK_SET  FILE_BEGIN
#define FIO_SEEK_END  FILE_END

static int fio_win32convertfilename(const char *filename, char *newfilename, int maxlen) {
  int i;
  int len=strlen(filename);
 
  if ((len + 1) >= maxlen)
    return -1;
   
  for (i=0; i<len; i++) {
    if (filename[i] == '/')
      newfilename[i] = '\\';
    else
      newfilename[i] = filename[i];
  }
  newfilename[len] = '\0'; /* NUL terminate the string */

  return 0;
}

static int fio_open(const char *filename, int mode, fio_fd *fd) {
  HANDLE fp;
  char winfilename[8192];
  DWORD access;
  DWORD sharing;
  LPSECURITY_ATTRIBUTES security;
  DWORD createmode;
  DWORD flags;

  if (fio_win32convertfilename(filename, winfilename, sizeof(winfilename)))
    return -1;  

  access = 0;
  if (mode & FIO_READ)
    access |= GENERIC_READ;
  if (mode & FIO_WRITE)
    access |= GENERIC_WRITE;
#if 0
  access = FILE_ALL_ACCESS; /* XXX hack if above modes fail */
#endif
#if 1
  if (mode & FIO_DIRECT)
    flags = FILE_FLAG_NO_BUFFERING;
  else
    flags = FILE_ATTRIBUTE_NORMAL;
#else
  if (mode & FIO_DIRECT)
    return -1; /* not supported yet */
#endif

  sharing = 0;       /* disallow sharing with other processes  */
  security = NULL;   /* child processes don't inherit anything */

  /* since we never append, blow away anything that's already there */
  if (mode & FIO_WRITE)
    createmode = CREATE_ALWAYS;
  else 
    createmode = OPEN_EXISTING;

  fp = CreateFile(winfilename, access, sharing, security, 
                  createmode, flags, NULL);

  if (fp == NULL) {
    return -1;
  } else {
    *fd = fp;
    return 0;
  }
}


static int fio_fclose(fio_fd fd) {
  BOOL rc;
  rc = CloseHandle(fd);
  if (rc) 
    return 0;
  else 
    return -1;
}

static fio_size_t fio_fread(void *ptr, fio_size_t size, 
                            fio_size_t nitems, fio_fd fd) {
  BOOL rc;
  DWORD len;
  DWORD readlen;

  len = size * nitems;

  rc = ReadFile(fd, ptr, len, &readlen, NULL);
  if (rc) {
    if (readlen == len)
      return nitems;
    else 
      return 0;
  } else {
    return 0;
  }
}

static fio_size_t fio_readv(fio_fd fd, const fio_iovec * iov, int iovcnt) {
  int i;
  fio_size_t len = 0; 

  for (i=0; i<iovcnt; i++) {
    fio_size_t rc = fio_fread(iov[i].iov_base, iov[i].iov_len, 1, fd);
    if (rc != 1)
      break;
    len += iov[i].iov_len;
  }

  return len;
}

static fio_size_t fio_fwrite(void *ptr, fio_size_t size, 
                             fio_size_t nitems, fio_fd fd) {
  BOOL rc;
  DWORD len;
  DWORD writelen;

  len = size * nitems; 
 
  rc = WriteFile(fd, ptr, len, &writelen, NULL);
  if (rc) {
    if (writelen == len)
      return nitems;
    else
      return 0;
  } else {
    return 0;
  }
}

static fio_size_t fio_fseek(fio_fd fd, fio_size_t offset, int whence) {
#if 1
  /* code that works with older MSVC6 compilers */
  LONGLONG finaloffset;
  LARGE_INTEGER bigint;
  LARGE_INTEGER finalint;

  bigint.QuadPart = offset;
  finalint = bigint;      /* set the high part, which will be overwritten */
  finalint.LowPart = SetFilePointer(fd, bigint.LowPart, &finalint.HighPart, whence);
  if (finalint.LowPart == -1) {
    /* if (finalint.LowPart == INVALID_SET_FILE_POINTER) { */
    /* INVALID_SET_FILE_POINTER is a possible "ok" low order result when */
    /* working with 64-bit offsets, so we have to also check the system  */
    /* error value for this thread to be sure */
    if (GetLastError() != ERROR_SUCCESS) {
      return -1;
    }
  } 

  finaloffset = finalint.QuadPart;
  return 0;
#else
  BOOL rc;
  LONGLONG finaloffset;

  /* SetFilePointerEx() only exists with new .NET compilers */
  rc = SetFilePointerEx(fd, offset, &finaloffset, whence);

  if (rc) 
    return 0;
  else
    return -1;
#endif
}

static fio_size_t fio_ftell(fio_fd fd) {
  /* code that works with older MSVC6 compilers */
  LONGLONG finaloffset;
  LARGE_INTEGER bigint;
  LARGE_INTEGER finalint;

  bigint.QuadPart = 0;
  finalint = bigint;      /* set the high part, which will be overwritten */

  finalint.LowPart = SetFilePointer(fd, bigint.LowPart, &finalint.HighPart, FILE_CURRENT);
  if (finalint.LowPart == -1) {
    /* if (finalint.LowPart == INVALID_SET_FILE_POINTER) { */
    /* INVALID_SET_FILE_POINTER is a possible "ok" low order result when */
    /* working with 64-bit offsets, so we have to also check the system  */
    /* error value for this thread to be sure */
    if (GetLastError() != ERROR_SUCCESS) {
      return -1;
    }
  }

  finaloffset = finalint.QuadPart;

  return finaloffset;
}

/*--------------------------------------------------------------------------------------------------------------------------------*/
#else /* . Windows version for machines with plain old ANSI C. */
/*--------------------------------------------------------------------------------------------------------------------------------*/

#include <stdio.h>
#include <string.h>

typedef FILE * fio_fd;
typedef size_t fio_size_t;  /* MSVC doesn't uinversally support ssize_t */
typedef void * fio_caddr_t; /* MSVC doesn't universally support caddr_t */

typedef struct {
  fio_caddr_t iov_base;
  int iov_len;
} fio_iovec;

#define FIO_SEEK_CUR SEEK_CUR
#define FIO_SEEK_SET SEEK_SET
#define FIO_SEEK_END SEEK_END

static int fio_open(const char *filename, int mode, fio_fd *fd) {
  char * modestr;
  FILE *fp;

  if (mode & FIO_READ) 
    modestr = "rb";

  if (mode & FIO_WRITE) 
    modestr = "wb";

  if (mode & FIO_DIRECT)
    return -1; /* not supported yet */
 
  fp = fopen(filename, modestr);
  if (fp == NULL) {
    return -1;
  } else {
    *fd = fp;
    return 0;
  }
}

static int fio_fclose(fio_fd fd) {
  return fclose(fd);
}

static fio_size_t fio_fread(void *ptr, fio_size_t size, 
                            fio_size_t nitems, fio_fd fd) {
  return fread(ptr, size, nitems, fd);
}

static fio_size_t fio_readv(fio_fd fd, const fio_iovec * iov, int iovcnt) {
  int i;
  fio_size_t len = 0; 

  for (i=0; i<iovcnt; i++) {
    fio_size_t rc = fread(iov[i].iov_base, iov[i].iov_len, 1, fd);
    if (rc != 1)
      break;
    len += iov[i].iov_len;
  }

  return len;
}

static fio_size_t fio_fwrite(void *ptr, fio_size_t size, 
                             fio_size_t nitems, fio_fd fd) {
  return fwrite(ptr, size, nitems, fd);
}

static fio_size_t fio_fseek(fio_fd fd, fio_size_t offset, int whence) {
  return fseek(fd, offset, whence);
}

static fio_size_t fio_ftell(fio_fd fd) {
  return ftell(fd);
}

#endif /* plain ANSI C */

/*==================================================================================================================================
! . Version for UNIX machines.
!=================================================================================================================================*/
#else

#if defined(__linux)
#ifndef _GNU_SOURCE
#define _GNU_SOURCE            /* required for O_DIRECT */
#endif
#endif
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h> /* . MJF: includes O_DIRECT. */
#include <string.h>

typedef int fio_fd;
typedef off_t fio_size_t;      /* off_t is 64-bits with LFS builds */

/* enable use of kernel readv() if available */
#if defined(__sun) || defined(__APPLE_CC__) || defined(__linux)
#define USE_KERNEL_READV 1
#endif

typedef void * fio_caddr_t;

#if defined(USE_KERNEL_READV)
#include <sys/uio.h>
typedef struct iovec fio_iovec;
#else

typedef struct {
  fio_caddr_t iov_base;
  int iov_len;
} fio_iovec;
#endif

#define FIO_SEEK_CUR SEEK_CUR
#define FIO_SEEK_SET SEEK_SET
#define FIO_SEEK_END SEEK_END

static int fio_open(const char *filename, int mode, fio_fd *fd) {
  int nfd;
  int oflag = 0;
  
  if (mode & FIO_READ) 
    oflag = O_RDONLY;

  if (mode & FIO_WRITE) 
    oflag = O_WRONLY | O_CREAT | O_TRUNC;

#if defined(__linux)
  /* enable direct I/O, requires block-aligned buffers and I/O sizes */
  if (mode & FIO_DIRECT)
    oflag |= O_DIRECT;
#else
  if (mode & FIO_DIRECT)
    return -1; /* not supported yet */
#endif

  nfd = open(filename, oflag, 0666);
  if (nfd < 0) {
    return -1;
  } else {
    *fd = nfd;
    return 0;
  }
}

static int fio_fclose(fio_fd fd) {
  return close(fd);
}

static fio_size_t fio_fread(void *ptr, fio_size_t size, 
                            fio_size_t nitems, fio_fd fd) {
  int i;
  fio_size_t len = 0; 
  int cnt = 0;

  for (i=0; i<nitems; i++) {
    fio_size_t rc = read(fd, (void*) (((char *) ptr) + (cnt * size)), size);
    if (rc != size)
      break;
    len += rc;
    cnt++;
  }

  return cnt;
}

static fio_size_t fio_readv(fio_fd fd, const fio_iovec * iov, int iovcnt) {
#if defined(USE_KERNEL_READV)
  return readv(fd, iov, iovcnt);
#else
  int i;
  fio_size_t len = 0; 

  for (i=0; i<iovcnt; i++) {
    fio_size_t rc = read(fd, iov[i].iov_base, iov[i].iov_len);
    if (rc != iov[i].iov_len)
      break;
    len += iov[i].iov_len;
  }

  return len;
#endif
}

static fio_size_t fio_fwrite(void *ptr, fio_size_t size, 
                             fio_size_t nitems, fio_fd fd) {
  int i;
  fio_size_t len = 0; 
  int cnt = 0;

  for (i=0; i<nitems; i++) {
    fio_size_t rc = write(fd, ptr, size);
    if (rc != size)
      break;
    len += rc;
    cnt++;
  }

  return cnt;
}

static fio_size_t fio_fseek(fio_fd fd, fio_size_t offset, int whence) {
 if (lseek(fd, offset, whence) >= 0)
   return 0;  /* success (emulate behavior of fseek) */
 else 
   return -1; /* failure (emulate behavior of fseek) */
}

static fio_size_t fio_ftell(fio_fd fd) {
  return lseek(fd, 0, SEEK_CUR);
}

#endif /* . _MSC_VER. */
/*--------------------------------------------------------------------------------------------------------------------------------*/

/* higher level routines that are OS independent */

/* . MJF: unnecessary.
static int fio_read_int32(fio_fd fd, int *i) {
  return (fio_fread(i, 4, 1, fd) != 1);
}
*/

/* . MJF: unnecessary.
static int fio_write_str(fio_fd fd, const char *str) {
  int len = strlen(str);
  return (fio_fwrite((void *) str, len, 1, fd) != 1);
}
*/

#endif /* . FASTIO_H. */
/*--------------------------------------------------------------------------------------------------------------------------------*/
