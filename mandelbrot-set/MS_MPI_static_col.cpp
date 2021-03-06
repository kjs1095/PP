#include <stdio.h>
#include "mpi.h"
#include <X11/Xlib.h>
#include <sys/time.h>
#include <unistd.h>
#include <string.h>
#define SAME 		0
#define MASTER 		0
#define MSG_TAG		0
#define DATA_TAG	1
#define STANDBY		-1
#define NO_RESULT	-2

typedef struct complexType
{	double real, imag;
} Compl;

int main(int argc, char *argv[])
{	int thread_n;
	double x_left, x_right, y_upper, y_lower;
	int x_point_n, y_point_n;					// x window size
	char cmd[10];								// command to enable x window
	int **canvas;								// for draw mandelbort set
	int draw;

	MPI_Status status;
	int myRank, nProc;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	MPI_Comm_size(MPI_COMM_WORLD, &nProc);

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
		thread_n = 2;
		x_point_n = y_point_n = 400;
		x_left = y_lower = -2.0;
		x_right = y_upper = 2.0;
		strcpy(cmd, "enable");
	}

	if (strcmp(cmd, "enable") == SAME)	draw = 1;
	else								draw = 0;

	int v_y_point_n;
	if (y_point_n % nProc != 0)	v_y_point_n = (y_point_n/nProc +1) * nProc;
	else						v_y_point_n = y_point_n;

	canvas = new int*[v_y_point_n];
	for (int r = 0; r < v_y_point_n; ++r)
		canvas[r] = new int[x_point_n];

	Display *display;
	Window window;      //initialization for a window
	int screen;         //which screen 

	// create graph 
	GC gc;
	XGCValues values;
	long valuemask = 0;

	if (draw && myRank == MASTER) {
		// open connection with the server  
		display = XOpenDisplay(NULL);
		if(display == NULL) {
			fprintf(stderr, "cannot open display\n");
			return 0;//exit(1);
		}

		screen = DefaultScreen(display);

		// set window position 
		int window_x = 0;
		int window_y = 0;

		// border width in pixels 
		int border_width = 0;

		// create window 
		window = XCreateSimpleWindow(display, RootWindow(display, screen), window_x, window_y, x_point_n, y_point_n, border_width,
					BlackPixel(display, screen), WhitePixel(display, screen));
	
		gc = XCreateGC(display, window, valuemask, &values);
		//XSetBackground (display, gc, WhitePixel (display, screen));
		XSetForeground (display, gc, BlackPixel (display, screen));
		XSetBackground(display, gc, 0X0000FF00);
		XSetLineAttributes (display, gc, 1, LineSolid, CapRound, JoinRound);
	
		// map(show) the window 
		XMapWindow(display, window);
		XSync(display, 0);
	}

	// draw points 
	Compl z, c;
	int repeats;
	double temp, lengthsq;
	int msg, slave_rank;

	struct timeval tv, tv2;
	unsigned long long start_utime, end_utime;
	gettimeofday(&tv, NULL);
	start_utime = tv.tv_sec * 1000000 + tv.tv_usec;

	int gap = 0;
	gap = v_y_point_n / nProc;
	for (int task_col = gap*myRank; task_col < gap*(myRank +1); ++task_col) {
		for (int row = 0; row < x_point_n; ++row) {
			z.real = 0.0;
			z.imag = 0.0;
			c.real = x_left + (double)row * ((x_right - x_left)/(double)x_point_n);
			c.imag = y_lower + (double)task_col * ((y_upper - y_lower)/(double)y_point_n);
			repeats = 0;
			lengthsq = 0.0;

			while (repeats < 10000 && lengthsq < 4.0) { // Theorem : If c belongs to M, then |Zn| <= 2. So Zn^2 <= 4 
				temp = z.real*z.real - z.imag*z.imag + c.real;
				z.imag = 2.0 * z.real * z.imag + c.imag;
				z.real = temp;
				lengthsq = z.real*z.real + z.imag*z.imag;
				repeats++;
			}

			canvas[task_col][row] = 1024*1024*(repeats%256);
		}
	}

	if (myRank == MASTER) {
		int working_n = nProc -1;
		while (working_n > 0) {
			MPI_Recv(&msg, 1, MPI_INT, MPI_ANY_SOURCE, MSG_TAG, MPI_COMM_WORLD, &status);
			slave_rank = status.MPI_SOURCE;
			
			if (msg == NO_RESULT)
				--working_n;
			else
				MPI_Recv(canvas[msg], x_point_n, MPI_INT, slave_rank, DATA_TAG, MPI_COMM_WORLD, &status);
		} 
	} else {
		for (int c = gap * myRank; c < gap*(myRank+1); ++c) {
			MPI_Send(&c, 1, MPI_INT, MASTER, MSG_TAG, MPI_COMM_WORLD);
			MPI_Send(canvas[c], x_point_n, MPI_INT, MASTER, DATA_TAG, MPI_COMM_WORLD);
		}
		msg = NO_RESULT;
		MPI_Send(&msg, 1, MPI_INT, MASTER, MSG_TAG, MPI_COMM_WORLD);
	}

	gettimeofday(&tv2, NULL);
	end_utime = tv2.tv_sec * 1000000 + tv2.tv_usec;
	//printf("rank %d's runtime = %llu.%03llu ms, process %d rows\n", myRank, (end_utime - start_utime)/1000, (end_utime - start_utime)%1000, gap);
	printf("%d\t%llu.%03llu\t%d\n", myRank, (end_utime - start_utime)/1000, (end_utime - start_utime)%1000, gap);
	if (draw && myRank == MASTER) {
		for (int i = 0; i < x_point_n; ++i)
			for (int j = 0; j < y_point_n; ++j) {
				XSetForeground (display, gc, canvas[j][i]);
				XDrawPoint (display, window, gc, i, y_point_n -j -1);
			}
		
		XFlush(display);
		sleep(10);
	}
	
	for (int r = 0; r < y_point_n; ++r)
		delete[] canvas[r];
	delete[] canvas;
	
	MPI_Finalize();

	return 0;
}
