#include "pch.h"

#include "gausslib.h"

static int swap_rows_ith_nonzero(int n, int i, double ** matrix);

void echelon_form_1t(int n, double ** matrix)
{
	double c; 
	for (int i = 0; i < n; i++) {
		if (matrix[i][i] == 0)
			if (swap_rows_ith_nonzero(n, i, matrix) == n) continue;

		c = matrix[i][i];
		
		for (int j = n; j >= i; j--) {
			matrix[i][j] /= c;
		}

		for (int j = i + 1; j < n; j++) {
			c = matrix[j][i];
			for (int k = n; k >= i; k--)
			{
				matrix[j][k] -= c * matrix[i][k];
			}
		}
	}
}

void back_substitution(int n, double ** matrix, double * x)
{
	x[n - 1] = matrix[n - 1][n];
	for (int i = n - 2; i >= 0; i--)
	{
		x[i] = matrix[i][n];
		for (int j = i + 1; j < n; j++)
		{
			x[i] -= matrix[i][j] * x[j];
		}
	}
}

static int swap_rows_ith_nonzero(int n, int i, double ** matrix)
{
	int j; double * tmp;
	for (j = i + 1; (j < n) && (matrix[j][i] == 0); j++);
	if (j == n)
	{
		return j;
	}
	else
	{
		tmp = matrix[i];
		matrix[i] = matrix[j];
		matrix[j] = tmp;
		return j;
	}
}
