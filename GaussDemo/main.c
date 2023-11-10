#include "gausslib.h"

#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <string.h>

#define N1 4
#define N2 5

#define EPSILON 0.0001
#define RUN_COUNT 5

#define USAGE "Usage: [v]alidate | [g]enerate n:int r:int outname:str | [p]rint inputname:str | [s]ingle-thread inputname:str is_verbose:int | [m]ulti-thread inputname:str is_verbose:int thread_count:int"

static void validate1(void(*echelon_form)(int, double **));
static void validate2(void(*echelon_form)(int, double **));
static void validate(void(*echelon_form)(int, double **),
					 int n, double ** matrix, double * x, double * expected_x);

static void performance(void(*echelon_form)(int, double **),
						int n, double ** original_matrix, double ** copy_matrix, double * x, int is_verbose);

static void matrix_generate(int n, int r, FILE * f);
static void matrix_read(int * n, double *** matrix, FILE * f);

static double ** matrix_malloc(int n);
static void matrix_print(int n, double ** matrix);
static void matrix_free(int n, double ** matrix);

int main(int argc, char ** argv)
{
	if (argc < 2)
	{
		puts(USAGE);
		return EXIT_SUCCESS;
	}

	int n, is_verbose;
	FILE * f = NULL;
	double ** matrix; double * x; double ** copy_matrix;
	switch (argv[1][0])
	{
		case 'v':

			validate1(echelon_form_1t);
			putchar('\n');
			validate2(echelon_form_1t);
			putchar('\n');
			validate1(echelon_form_mt);
			putchar('\n');
			validate2(echelon_form_mt);

			break;

		case 'g':

			if (argc != 5)
			{
				puts(USAGE);
				return EXIT_SUCCESS;
			}

			n = atoi(argv[2]);
			if (n <= 0)
			{
				printf("Invalid matrix dimension: %s\n", argv[2]);
				return EXIT_SUCCESS;
			}

			int r = atoi(argv[3]);
			if (r <= 0)
			{
				printf("Invalid matrix dimension: %s\n", argv[3]);
				return EXIT_SUCCESS;
			}

			fopen_s(&f, argv[4], "w");
			if (f == NULL)
			{
				printf("Can't open file: %s\n", argv[4]);
				return EXIT_SUCCESS;
			}

			matrix_generate(n, r, f);

			fflush(f);
			fclose(f);

			break;

		case 'p':

			if (argc != 3)
			{
				puts(USAGE);
				return EXIT_SUCCESS;
			}

			fopen_s(&f, argv[2], "r");
			if (f == NULL)
			{
				printf("Can't open file: %s\n", argv[2]);
				return EXIT_SUCCESS;
			}

			matrix_read(&n, &matrix, f);
			matrix_print(n, matrix);
			matrix_free(n, matrix);

			break;

		case 's':

			if (argc != 4)
			{
				puts(USAGE);
				return EXIT_SUCCESS;
			}

			fopen_s(&f, argv[2], "r");
			if (f == NULL)
			{
				printf("Can't open file: %s\n", argv[2]);
				return EXIT_SUCCESS;
			}

			is_verbose = atoi(argv[3]);
			if (is_verbose < 0)
			{
				printf("Invalid verbose flag: %s\n", argv[3]);
				return EXIT_SUCCESS;
			}

			matrix_read(&n, &matrix, f);
			x = (double *)malloc(n * sizeof(double));
			copy_matrix = matrix_malloc(n);

			performance(echelon_form_1t, n, matrix, copy_matrix, x, is_verbose);
			
			free(x);
			matrix_free(n, matrix);
			matrix_free(n, copy_matrix);

			break;

		case 'm':

			if (argc != 5)
			{
				puts(USAGE);
				return EXIT_SUCCESS;
			}

			fopen_s(&f, argv[2], "r");
			if (f == NULL)
			{
				printf("Can't open file: %s\n", argv[2]);
				return EXIT_SUCCESS;
			}

			is_verbose = atoi(argv[3]);
			if (is_verbose < 0)
			{
				printf("Invalid verbose flag: %s\n", argv[3]);
				return EXIT_SUCCESS;
			}

			int thread_count = atoi(argv[4]);
			if (thread_count <= 0)
			{
				printf("Invalid thread count: %s\n", argv[4]);
				return EXIT_SUCCESS;
			}

			matrix_read(&n, &matrix, f);
			x = (double *)malloc(n * sizeof(double));
			copy_matrix = matrix_malloc(n);

			omp_set_num_threads(thread_count);
			performance(echelon_form_mt, n, matrix, copy_matrix, x, is_verbose);
			
			free(x);
			matrix_free(n, matrix);
			matrix_free(n, copy_matrix);

			break;

		default:
			puts(USAGE);
			break;
	}

	return EXIT_SUCCESS;
}

static void validate1(void(*echelon_form)(int, double **))
{
	double ** matrix = matrix_malloc(N1);
	matrix[0][0] = 1; matrix[0][1] = 2; matrix[0][2] = 3; matrix[0][3] = 4; matrix[0][4] = 0;
	matrix[1][0] = 7; matrix[1][1] = 14; matrix[1][2] = 20; matrix[1][3] = 27; matrix[1][4] = 0;
	matrix[2][0] = 5; matrix[2][1] = 10; matrix[2][2] = 16; matrix[2][3] = 19; matrix[2][4] = -2;
	matrix[3][0] = 3; matrix[3][1] = 5; matrix[3][2] = 6; matrix[3][3] = 13; matrix[3][4] = 5;
	double x[N1];
	double expected_x[N1] = { 1, -1, -1, 1 };

	validate(echelon_form, N1, matrix, x, expected_x);

	matrix_free(N1, matrix);
}

static void validate2(void(*echelon_form)(int, double **))
{
	double ** matrix = matrix_malloc(N2);
	matrix[0][0] = 2; matrix[0][1] = -1; matrix[0][2] = -1; matrix[0][3] = -4; matrix[0][4] = -1; matrix[0][5] = 2;
	matrix[1][0] = -1; matrix[1][1] = 2; matrix[1][2] = -1; matrix[1][3] = -1; matrix[1][4] = -1; matrix[1][5] = 0;
	matrix[2][0] = 4; matrix[2][1] = 1; matrix[2][2] = -5; matrix[2][3] = -8; matrix[2][4] = -5; matrix[2][5] = 1;
	matrix[3][0] = 1; matrix[3][1] = 1; matrix[3][2] = 2; matrix[3][3] = 1; matrix[3][4] = 1; matrix[3][5] = 0;
	matrix[4][0] = 1; matrix[4][1] = 1; matrix[4][2] = 1; matrix[4][3] = 2; matrix[4][4] = 1; matrix[4][5] = 0.5;
	double x[N2];
	double expected_x[N2] = { 5.0/18.0, 4.0/9.0, (-4.0)/3.0, (-5.0)/6.0, 25.0/9.0 };

	validate(echelon_form, N2, matrix, x, expected_x);

	matrix_free(N2, matrix);
}

static void validate(void(*echelon_form)(int, double **),
					 int n, double ** matrix, double * x, double * expected_x)
{
	echelon_form(n, matrix);
	back_substitution(n, matrix, x);
	for (int i = 0; i < n; i++)
	{
		printf("| %lf ", x[i]);
		if (fabs(x[i] - expected_x[i]) > EPSILON)
		{
			printf("Assertion failed at %d", i);
			break;
		}
	}
	putchar('|');
}

static void performance(void(*echelon_form)(int, double **),
						int n, double ** original_matrix, double ** copy_matrix, double * x, int is_verbose)
{
	double s_time, e_time, time, min_time = DBL_MAX;
	for (int i = 0; i < RUN_COUNT; i++)
	{
		for (int j = 0; j < n; j++)
			memcpy_s(copy_matrix[j], (n + 1) * sizeof(double), original_matrix[j], (n + 1) * sizeof(double));

		s_time = omp_get_wtime();
		echelon_form(n, copy_matrix);
		back_substitution(n, copy_matrix, x);
		e_time = omp_get_wtime();
		time = e_time - s_time;

		printf("run %d: time = %lf;\n", i, time);

		if (is_verbose)
		{
			for (int j = 0; j < n; j++)
				printf("| %lf ", x[j]);
			putchar('\n');
		}

		if (min_time > time)
			min_time = time;
	}
	printf("Best time: %lf\n", min_time);
}

static void matrix_generate(int n, int r, FILE * f)
{
	fwrite(&n, sizeof(int), 1, f);

	double d;
	for (int i = 0; i < n * (n + 1); i++)
	{
		d = rand() % r;
		fwrite(&d, sizeof(double), 1, f);
	}
}

static void matrix_read(int * n, double *** matrix, FILE * f)
{
	int _n;
	fread(&_n, sizeof(int), 1, f);

	double ** _matrix = matrix_malloc(_n);

	for (int i = 0; i < _n; i++)
	{
		fread(_matrix[i], sizeof(double), _n + 1, f);
	}

	*n = _n;
	*matrix = _matrix;
}

static double ** matrix_malloc(int n)
{
	double ** m = (double **)malloc(n * sizeof(double *));
	for (int i = 0; i < n; i++)
	{
		m[i] = (double *)malloc((n + 1) * sizeof(double));
	}
	return m;
}

static void matrix_print(int n, double ** matrix)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n + 1; j++)
		{
			printf("| %lf ", matrix[i][j]);
		}
		puts("|\n");
	}
}

static void matrix_free(int n, double ** matrix)
{
	for (int i = 0; i < n; i++)
	{
		free(matrix[i]);
	}
	free(matrix);
}