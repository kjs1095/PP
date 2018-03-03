#include <stdio.h>
#include <sys/time.h>
#include <algorithm>
#include "mpi.h"
using std::sort;
#define INF 2147483647
void MergeSplit(int local_n, int *arr, int *recv_arr, int smallPart)
{	int *tmp = new int[local_n];
	for (int i = 0; i < local_n; ++i)
		tmp[i] = arr[i];
	int i, j, k;
	if (smallPart) {
		i = j = k = 0;
		while (k < local_n) 
			if (tmp[i] < recv_arr[j])	arr[k++] = tmp[i++];
			else						arr[k++] = recv_arr[j++];
	} else {
		i = j = k = local_n -1;
		while (k >= 0) {
			if (tmp[i] > recv_arr[j])	arr[k--] = tmp[i--];
			else						arr[k--] = recv_arr[j--];
		}
	}

	delete [] tmp;
}
int main(int argc, char *argv[])
{	int rank, nProc;
	int n;
	int local_n;
	int *local_arr;
	int *recv_arr;
	int bufSize;
	int oddPhaseRank;
	int evenPhaseRank;
	MPI_Status status;
	MPI_File fh;

	struct timeval tv, tv2;
	unsigned long long start_utime, end_utime, io_utime = 0;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nProc);

	sscanf(argv[1], "%d", &n);
	local_n = (n + nProc -1) / nProc;
	local_arr = new int[local_n];
	recv_arr = new int[local_n];
	
	for (int i = 0; i < local_n; ++i)
		local_arr[i] = INF;

	if (rank * local_n < n) {
		bufSize = local_n;
		if (n - local_n * rank < local_n)
			bufSize = n - local_n * rank;
	} else {
		bufSize = 0;
	}
	
	MPI_File_open(MPI_COMM_WORLD, argv[2], MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
	MPI_File_read_at_all(fh, rank * local_n * sizeof(int), local_arr, bufSize, MPI_INT, &status); 
	MPI_File_close(&fh);

	sort(local_arr, local_arr + local_n);

	if (rank %2 == 0) {
		evenPhaseRank = rank +1;
		oddPhaseRank = rank -1;
	} else {
		evenPhaseRank = rank -1;
		oddPhaseRank = rank +1;
	}

	if (oddPhaseRank < 0 || oddPhaseRank == nProc)
		oddPhaseRank = MPI_PROC_NULL;
	if (evenPhaseRank < 0 || evenPhaseRank == nProc)
		evenPhaseRank = MPI_PROC_NULL;

	for (int phase = 0; phase < nProc; ++phase) {
		if (phase %2 == 0) {	// even phase
			if (rank %2 == 0) {
				MPI_Send(local_arr, local_n, MPI_INT, evenPhaseRank, 0, MPI_COMM_WORLD);
				MPI_Recv(recv_arr, local_n, MPI_INT, evenPhaseRank, 0, MPI_COMM_WORLD, &status);
			} else {
				MPI_Recv(recv_arr, local_n, MPI_INT, evenPhaseRank, 0, MPI_COMM_WORLD, &status);
				MPI_Send(local_arr, local_n, MPI_INT, evenPhaseRank, 0, MPI_COMM_WORLD);
			}
		} else {				// odd phase
			if (rank %2 == 0) {
				MPI_Recv(recv_arr, local_n, MPI_INT, oddPhaseRank, 0, MPI_COMM_WORLD, &status);
				MPI_Send(local_arr, local_n, MPI_INT, oddPhaseRank, 0, MPI_COMM_WORLD);
			} else {
				MPI_Send(local_arr, local_n, MPI_INT, oddPhaseRank, 0, MPI_COMM_WORLD);
				MPI_Recv(recv_arr, local_n, MPI_INT, oddPhaseRank, 0, MPI_COMM_WORLD, &status);
			}
		}

		if (status.MPI_SOURCE != MPI_PROC_NULL)
			MergeSplit(local_n, local_arr, recv_arr, rank < status.MPI_SOURCE); 
	}

	MPI_File_open(MPI_COMM_WORLD, argv[3], MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
	MPI_File_write_at_all(fh, rank * local_n * sizeof(int), local_arr, bufSize, MPI_INT, &status);
	MPI_File_close(&fh);

	MPI_Finalize();
	
	delete[] local_arr;
	delete[] recv_arr;
	return 0;
}
