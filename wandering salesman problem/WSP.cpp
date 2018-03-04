#include <stdio.h>
#include <sys/time.h>
#include <pthread.h>
#include <algorithm>
#include <queue>
using std::min;
using std::queue;
#define NONE -1
#define SOURCE 0

const int N = 20 +1;
const int INF = 0x7fffffff;
const int THREAD_N = 100;
const int MODE_BOUND = 20;

typedef struct
{	int walked[N];
	int cost;	
}Path;

/* global protect!!! */
int path[N];
int ans;
queue<Path> src;
/* global protect!!! */

int Map[N][N];
int n;

pthread_mutex_t global_ans_mtx;
pthread_mutex_t global_queue_mtx;

/******************* BB part ********************/
void tsp_bf(int n_last, int tmp_length, int tmp_path[], int undo, int status)
{	if (undo == 0) {
		/* critical section */
		pthread_mutex_lock( &global_ans_mtx );
		if (tmp_length < ans || (tmp_length == ans && path[0] < tmp_path[0])) {
			ans = tmp_length;
			for (int i = 0; i < n; ++i)
				path[i] = tmp_path[i];
		}	
		pthread_mutex_unlock( &global_ans_mtx );
		/* critical section */
	}

	if (ans <= tmp_length)
		return ;

	for (int nxt_last = 1; nxt_last <= n; ++nxt_last) {
		if ((status & (1<<(nxt_last-1))) && Map[n_last][nxt_last] != NONE) {
			status -= (1<<(nxt_last-1));
			tmp_path[n - undo] = nxt_last;

			tsp_bf(nxt_last, tmp_length + Map[n_last][nxt_last], tmp_path, undo -1, status);

			status += (1<<(nxt_last-1));
		}
	}	
}
void *run(void *param)
{	unsigned long long t_start_utime, t_end_utime;
	struct timeval t_tv1, t_tv2;
	gettimeofday(&t_tv1, NULL);
	t_start_utime = t_tv1.tv_sec * 1000000 + t_tv1.tv_usec;

	int *tid = (int*)param;
	Path start;
	while (1) {
		pthread_mutex_lock( &global_queue_mtx );
		if (src.empty()) {
			pthread_mutex_unlock( &global_queue_mtx );
			break;
		}
		start = src.front();
		src.pop();
		pthread_mutex_unlock( &global_queue_mtx );

		int status = (1<<n)-1;
		int tmp_path[N];

		tmp_path[0] = start.walked[0];
		status -= (1<<((tmp_path[0]) -1));
		tsp_bf(tmp_path[0], start.cost, tmp_path, n-1, status);
	}

	gettimeofday(&t_tv2, NULL);
	t_end_utime = t_tv2.tv_sec * 1000000 + t_tv2.tv_usec;
	//printf("thread %d'sruntime = %llu\n", *tid, t_end_utime - t_start_utime);
	pthread_exit(NULL);
}
/****************** BB part *********************/
int main(int argc, char *argv[])
{	struct timeval tv1, tv2;
	unsigned long long start_utime, end_utime;
	gettimeofday(&tv1, NULL);
	start_utime = tv1.tv_sec * 1000000 + tv1.tv_usec;

	freopen(argv[2], "r", stdin);

	int thread_n;
	/* Input file */
	scanf("%d", &n);
	for (int i = 1; i <= n; ++i)
		for (int j = 1; j <= n; ++j) 
			scanf("%d", &Map[i][j]);
	/* End of Input file */

	if (argc == 2)	thread_n = n;
	else			sscanf(argv[1], "%d", &thread_n);

	//printf("%d\n", thread_n);
	pthread_t thread_id[THREAD_N];

	ans = INF;
	pthread_mutex_init( &global_ans_mtx, NULL);
	pthread_mutex_init( &global_queue_mtx, NULL);

	for (int i = 1; i <= n; ++i) {
		Path tmp;
		tmp.walked[0] = i;
		tmp.cost = 0;
		src.push(tmp);
	}

	int tida[THREAD_N];
	for (int i = 0; i < thread_n; ++i) {
		tida[i] = i;
		pthread_create( &thread_id[i], NULL, run, (void*)&tida[i]);
	}

	for (int i = 0; i < thread_n; ++i)
		pthread_join( thread_id[i], NULL);

	pthread_mutex_destroy( &global_ans_mtx );
	pthread_mutex_destroy( &global_queue_mtx );

	printf("%d:", ans);
	for (int i = 0; i < n; ++i) {
		if (i)
			printf(",");
		printf("%d", path[i]);
	}

	gettimeofday(&tv2, NULL);
	end_utime = tv2.tv_sec * 1000000 + tv2.tv_usec;

	//printf(" #thread = %d runtime = %llu\n", thread_n, end_utime - start_utime);
	pthread_exit(NULL);
	return 0;
}

