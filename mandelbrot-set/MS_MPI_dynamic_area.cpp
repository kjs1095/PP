#include <stdio.h>
#include "mpi.h"
#include <X11/Xlib.h>
#include <sys/time.h>
#include <unistd.h>
#include <string.h>
#include <algorithm>
using std::min;
#define SAME 		0
#define MASTER 		0
#define MSG_TAG		0
#define DATA_TAG	1
#define STANDBY		-1
#define NO_JOB		-2
//const int W = 1;
//const int H = 1;

typedef struct complexType
{	double real, imag;
} Compl;

typedef struct 
{	int top;
	int left;
	int bottom;
	int right;
}Area;

int main(int argc, char *argv[])
{	int thread_n;
	double x_left, x_right, y_upper, y_lower;
	int x_point_n, y_point_n;					// x window size
	char cmd[10];								// command to enable x window
	int **canvas;								// for draw mandelbort set
	int draw;
	Area areaInfo;

	MPI_Status status;
	int myRank, nProc;

	MPI_Datatype MPI_AREA;
	MPI_Datatype type[4] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT};
	int blocklen[4] = {1, 1, 1, 1};
	MPI_Aint disp[4];

	MPI_Init(&argc, &argv);

	/*disp[0] = &areaInfo.top - &areaInfo;
	disp[1] = &areaInfo.left - &areaInfo.top;
	disp[2] = &areaInfo.bottom - &areaInfo.left;
	disp[3] = &areaInfo.right - &areaInfo.bottom;
	*/
	//disp[0] = disp[1] = disp[2] = disp[3] = sizeof(int);
	disp[0] = offsetof(Area, top);
	disp[1] = offsetof(Area, left);
	disp[2] = offsetof(Area, bottom);
	disp[3] = offsetof(Area, right);
	MPI_Type_create_struct(4, blocklen, disp, type, &MPI_AREA);
	MPI_Type_commit(&MPI_AREA);

	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	MPI_Comm_size(MPI_COMM_WORLD, &nProc);

	int W, H;
	if (argc == 11) {
		sscanf(argv[1], "%d", &thread_n);
		sscanf(argv[2], "%lf", &x_left);
		sscanf(argv[3], "%lf", &x_right);
		sscanf(argv[4], "%lf", &y_lower);
		sscanf(argv[5], "%lf", &y_upper);
		sscanf(argv[6], "%d", &x_point_n);
		sscanf(argv[7], "%d", &y_point_n);
		sscanf(argv[8], "%s", cmd);
		sscanf(argv[9], "%d", &H);
		sscanf(argv[10], "%d", &W);
	} else {
		thread_n = 1;
		x_point_n = y_point_n = 400;
		x_left = y_lower = -2.0;
		x_right = y_upper = 2.0;
		strcpy(cmd, "enable");
	}

	if (strcmp(cmd, "enable") == SAME)	draw = 1;
	else								draw = 0;

	canvas = new int*[x_point_n];
	for (int r = 0; r < x_point_n; ++r)
		canvas[r] = new int[y_point_n];

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
	int slave_rank;
	int task_row;
	int job_n = 0;

	struct timeval tv, tv2;
	unsigned long long start_utime, end_utime;
	gettimeofday(&tv, NULL);
	start_utime = tv.tv_sec * 1000000 + tv.tv_usec;

	struct timeval ctv, ctv2;
	unsigned long long c_start_utime, c_end_utime, compu = 0;

	if (myRank == MASTER) {
		int working_n = 0;
		int task_x = 0, task_y = -W;
		while (task_x != NO_JOB || working_n > 0) {
			MPI_Recv(&areaInfo, 1, MPI_AREA, MPI_ANY_SOURCE, MSG_TAG, MPI_COMM_WORLD, &status);
			slave_rank = status.MPI_SOURCE;
			//printf("m get %d %d %d %d\n", areaInfo.top, areaInfo.bottom, areaInfo.left, areaInfo.right);	
			if (areaInfo.top != STANDBY) {
				for (int i = areaInfo.top; i < areaInfo.bottom; ++i)
					//for (int j = areaInfo.left; j < areaInfo.right; ++j)
					//	MPI_Recv(&canvas[i][j], 1, MPI_INT, slave_rank, DATA_TAG, MPI_COMM_WORLD, &status);
					MPI_Recv(&canvas[i][areaInfo.left], areaInfo.right - areaInfo.left, MPI_INT, slave_rank, DATA_TAG, MPI_COMM_WORLD, &status);
				--working_n;
		//		printf("update from %d, remain %d\n", slave_rank, working_n);
			}

			if (task_x != NO_JOB) {
				if (task_y + W >= y_point_n) {
					task_y = 0;
					task_x += H;
				} else {
					task_y += W;
				}

				if (task_x >= x_point_n) 	task_x = NO_JOB;
				else 						++working_n;
			}

			areaInfo.top = task_x;
			areaInfo.left = task_y;
			areaInfo.bottom = min(x_point_n, task_x + H);
			areaInfo.right = min(y_point_n, task_y + W);

			//printf("m send %d %d %d %d to %d, remain %d\n", areaInfo.top, areaInfo.bottom, areaInfo.left, areaInfo.right, slave_rank, working_n);
			MPI_Send(&areaInfo, 1, MPI_AREA, slave_rank, MSG_TAG, MPI_COMM_WORLD);
			//printf("%d %d\n", task_x, working_n);
		}
	} else {
		areaInfo.top = STANDBY;
		while (1) {
			MPI_Send(&areaInfo, 1, MPI_AREA, MASTER, MSG_TAG, MPI_COMM_WORLD);
		//	printf("s%d send %d %d %d %d\n", myRank, areaInfo.top, areaInfo.bottom, areaInfo.left, areaInfo.right);	
			if (areaInfo.top != STANDBY) {
				for (int i = areaInfo.top; i < areaInfo.bottom; ++i)
					//for (int j = areaInfo.left; j < areaInfo.right; ++j)	
					//	MPI_Send(&canvas[i][j], 1, MPI_INT, slave_rank, DATA_TAG, MPI_COMM_WORLD);
					MPI_Send(&canvas[i][areaInfo.left], areaInfo.right - areaInfo.left, MPI_INT, MASTER, DATA_TAG, MPI_COMM_WORLD);
			}

			MPI_Recv(&areaInfo, 1, MPI_AREA, MASTER, MSG_TAG, MPI_COMM_WORLD, &status);
		//	printf("s%d get %d %d %d %d\n", myRank, areaInfo.top, areaInfo.bottom, areaInfo.left, areaInfo.right);
			if (areaInfo.top == NO_JOB)	break;
		
			gettimeofday(&ctv, NULL);
			c_start_utime = ctv.tv_sec * 1000000 + ctv.tv_usec;
			for (int row = areaInfo.top; row < areaInfo.bottom; ++row) {
				for (int col = areaInfo.left; col < areaInfo.right; ++col) {
					z.real = 0.0;
					z.imag = 0.0;
					c.real = x_left + (double)row * ((x_right - x_left)/(double)x_point_n);
					c.imag = y_lower + (double)col * ((y_upper - y_lower)/(double)y_point_n);
					repeats = 0;
					lengthsq = 0.0;

					while (repeats < 10000 && lengthsq < 4.0) { // Theorem : If c belongs to M, then |Zn| <= 2. So Zn^2 <= 4 
						temp = z.real*z.real - z.imag*z.imag + c.real;
						z.imag = 2.0 * z.real * z.imag + c.imag;
						z.real = temp;
						lengthsq = z.real*z.real + z.imag*z.imag;
						repeats++;
						job_n++;
					}

					canvas[row][col] = 1024*1024*(repeats%256);
				}
			}
			gettimeofday(&ctv2, NULL);
			c_end_utime = ctv2.tv_sec * 1000000 + ctv2.tv_usec;
			compu += (c_end_utime - c_start_utime);
		}
	}
	
	gettimeofday(&tv2, NULL);
	end_utime = tv2.tv_sec * 1000000 + tv2.tv_usec;
	//printf("rank %d's runtime = %llu.%03llu ms, process %d rows\n", myRank, (end_utime - start_utime)/1000, (end_utime - start_utime)%1000, job_n);
	if (myRank != MASTER)
		printf("%d\t%llu.%03llu\t%llu.%03llu\t%d\n", myRank, (end_utime - start_utime)/1000, (end_utime - start_utime)%1000, compu/1000, compu%1000, job_n);

	if (draw && myRank == MASTER) {
		for (int i = 0; i < x_point_n; ++i)
			for (int j = 0; j < y_point_n; ++j) {
				XSetForeground (display, gc, canvas[i][j]);
				XDrawPoint (display, window, gc, i, y_point_n -j -1);
			}
		
		XFlush(display);
		sleep(5);
	}
	
	for (int r = 0; r < x_point_n; ++r)
		delete[] canvas[r];
	delete[] canvas;
	
	MPI_Finalize();

	return 0;
}
