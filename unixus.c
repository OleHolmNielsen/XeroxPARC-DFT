
/*
 * Returns # of major+minor page faults in npages,
 * and user CPU time in time.
 * Callable from Fortran as CALL UNIXUS(NPAGES,TIME)
 * 
 * Read manual pages for getrusage(2)
 * and the file /usr/include/sys/time.h for "timeval" definition
 */

#include <sys/time.h>
#include <sys/resource.h>
#include <stdio.h>

#ifdef hppa
#include <sys/times.h>
#include <unistd.h>
#endif

/* Hewlett-Packard Precision Architecture */
#ifdef hppa
unixus(npages,time)
#endif

/* Sun Microsystems SunOS */
#ifdef sun
unixus_(npages,time)
#endif

/* IBM RS/6000 AIX 3.1 */
#ifdef aix
unixus(npages,time)
#endif

int *npages;
float *time;
{
#ifdef hppa
	struct tms buffer;
	float ticks = (float) sysconf(_SC_CLK_TCK);

	if (times(&buffer) < 0) {
		perror ("usage - bad times return");
		return;
	}
	*time = ((float)buffer.tms_utime + (float)buffer.tms_stime) / ticks;
	
#else
	int who = RUSAGE_SELF;
	struct rusage *rusage;

	rusage = (struct rusage *) malloc(sizeof(struct rusage));
	getrusage(who, rusage);
	*npages = rusage->ru_majflt + rusage->ru_minflt;
	*time   =  (float)rusage->ru_utime.tv_sec  +
	  1.0e-6 * (float)rusage->ru_utime.tv_usec ;
	/*
	*time   =  (float)rusage.ru_utime->tv_sec  +
	  1.0e-6 * (float)rusage.ru_utime->tv_usec +
	           (float)rusage.ru_stime->tv_sec  +
	  1.0e-6 * (float)rusage.ru_stime->tv_usec; */
#endif
}
