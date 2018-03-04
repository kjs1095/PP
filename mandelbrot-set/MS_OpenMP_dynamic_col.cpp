 /* 
   Mandelbort sort OpenMP
   dynamic version, partition by rows
 */

#include <X11/Xlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <unistd.h>
#include <string.h>
#include <omp.h>

#define SAME 0
const double eps = 1e-7;
const int ChunkSize = 10;

typedef struct complexType
{	double real, imag;
} Compl;

int main(int argc, char* argv[])
{	int thread_n;
    double x_left, x_right, y_upper, y_lower;
    int x_point_n, y_point_n;					// x window size
    char cmd[10];								// command to enable x window
    
	int **canvas;								// for draw mandelbort set

	if (argc == 9) {
        sscanf(argv[1], "%d", &thread_n);
        sscanf(argv[2], "%lf", &x_left);
        sscanf(argv[3], "%lf", &x_right);
        sscanf(argv[4], "%lf", &y_lower);
        sscanf(argv[5], "%lf", &y_upper);
        sscanf(argv[6], "%d", &x_point_n);
		sscanf(argv[7], "%d", &y_point_n);
		sscanf(argv[8], "%s", cmd);
	} else {
		thread_n = 1;
		x_point_n = y_point_n = 400;
		x_left = y_lower = -2.0;
		x_right = y_upper = 2.0;
		strcpy(cmd, "enable");
	}

	int draw;
	if (strcmp(cmd, "enable") == SAME)	draw = 1;
	else								draw = 0;

	canvas = new int*[x_point_n];
	for (int r = 0; r < x_point_n; ++r)
		canvas[r] = new int[y_point_n];

	Display *display;
	Window window;      //initialization for a window
	int screen;         //which screen 

	/* create graph */
	GC gc;
	XGCValues values;
	long valuemask = 0;

	if (draw) {
		/* open connection with the server */ 
		display = XOpenDisplay(NULL);
		if(display == NULL) {
			fprintf(stderr, "cannot open display\n");
			return 0;//exit(1);
		}

		screen = DefaultScreen(display);

		/* set window position */
		int window_x = 0;
		int window_y = 0;

		/* border width in pixels */
		int border_width = 0;

		/* create window */
		window = XCreateSimpleWindow(display, RootWindow(display, screen), window_x, window_y, x_point_n, y_point_n, border_width,
					BlackPixel(display, screen), WhitePixel(display, screen));
	
		gc = XCreateGC(display, window, valuemask, &values);
		//XSetBackground (display, gc, WhitePixel (display, screen));
		XSetForeground (display, gc, BlackPixel (display, screen));
		XSetBackground(display, gc, 0X0000FF00);
		XSetLineAttributes (display, gc, 1, LineSolid, CapRound, JoinRound);
	
		/* map(show) the window */
		XMapWindow(display, window);
		XSync(display, 0);
	}

	/* Theorem : If |c| <= 1/4, then c belongs to M */

	/* draw points */
	Compl z, c;
	int repeats;
	double temp, lengthsq;
	int i, j;
	int job_n = 0;

	struct timeval tv, tv2;
	unsigned long long start_utime, end_utime;
	gettimeofday(&tv, NULL);
	start_utime = tv.tv_sec * 1000000 + tv.tv_usec;

	int tid;
	#pragma omp parallel private(tid, i, j, repeats, z, c, lengthsq, temp, job_n) num_threads(thread_n) 
	{	tid = omp_get_thread_num();
		job_n = 0;
		
		struct timeval t_tv, t_tv2;
		unsigned long long t_start_utime, t_end_utime;
		gettimeofday(&t_tv, NULL);
		t_start_utime = t_tv.tv_sec * 1000000 + t_tv.tv_usec;

		#pragma omp for schedule(dynamic) nowait
		for (i = 0; i < y_point_n; ++i) {
			++job_n;
			//#pragma omp for schedule(dynamic, ChunkSize) nowait 
			for (j = 0; j < x_point_n; ++j) {
				//printf("tid = %d, do %d.%d\n", tid, i, j);
				z.real = 0.0;
				z.imag = 0.0;
				c.real = x_left + (double)j * ((x_right - x_left)/(double)x_point_n);
				c.imag = y_lower + (double)i * ((y_upper - y_lower)/(double)y_point_n);
				repeats = 0;
				lengthsq = 0.0;

				while(repeats < 10000 && lengthsq < 4.0) { /* Theorem : If c belongs to M, then |Zn| <= 2. So Zn^2 <= 4 */
					temp = z.real*z.real - z.imag*z.imag + c.real;
					z.imag = 2*z.real*z.imag + c.imag;
					z.real = temp;
					lengthsq = z.real*z.real + z.imag*z.imag; 
					repeats++;
				}

				canvas[j][i] = 1024*1024*(repeats%256);
			}
		}

		gettimeofday(&t_tv2, NULL);
		t_end_utime = t_tv2.tv_sec * 1000000 + t_tv2.tv_usec;
		//printf("tid %2d's runtime = %llu.%03llu ms, process %d rows\n", tid, (t_end_utime - t_start_utime)/1000, (t_end_utime - t_start_utime)%1000, job_n);
		printf("%d\t%llu.%03llu\t%d\n", tid, (t_end_utime - t_start_utime)/1000, (t_end_utime - t_start_utime)%1000, job_n);
	}

	 gettimeofday(&tv2, NULL);
	 end_utime = tv2.tv_sec * 1000000 + tv2.tv_usec;

	printf("-1\t%llu.%03llu\t%d\n", (end_utime - start_utime)/1000, (end_utime - start_utime)%1000, thread_n);

	if (draw) {
		for (i = 0;i < x_point_n; ++i)
			for (j = 0; j < y_point_n; ++j) {
				XSetForeground (display, gc, canvas[i][j]);
				XDrawPoint (display, window, gc, i, y_point_n - j -1);
			}
		XFlush(display);
		sleep(5);
	}

	for (int r = 0; r < x_point_n; ++r)
		delete[] canvas[r];
	delete[] canvas;

	return 0;
}
