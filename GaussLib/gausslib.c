#include "pch.h"

#include "gausslib.h"

#include <omp.h>

static int swap_rows_ith_nonzero(int n, int i, double ** matrix);

/// <summary>
/// Single thread row echelon form implementation.
/// Mutates augmented matrix (NxN+1) into row echelon form.
/// </summary>
/// <param name="n"> Augmented matrix row count. </param>
/// <param name="matrix"> Augmented matrix. </param>
void echelon_form_1t(int n, double ** matrix)
{
	double c; 
	for (int i = 0; i < n; i++) // row cycle
	{
		// Check for zero division
		if (matrix[i][i] == 0)
			if (swap_rows_ith_nonzero(n, i, matrix) == 0) continue;

		// Coefficient to divide ith row with
		c = matrix[i][i];
		
		// Divide current row
		for (int j = n; j >= i; j--) // element cycle
		{
			matrix[i][j] /= c;
		}

		// Subtract ith row from all below
		for (int j = i + 1; j < n; j++) // row cycle
		{
			// Coefficient to multiply kth row with
			c = matrix[j][i];

			for (int k = n; k >= i; k--) // element cycle
			{
				matrix[j][k] -= c * matrix[i][k];
			}
		}
	}
}

/// <summary>
/// Multi thread row echelon form implementation.
/// Mutates augmented matrix (NxN+1) into row echelon form.
/// </summary>
/// <param name="n"> Augmented matrix row count. </param>
/// <param name="matrix"> Augmented matrix. </param>
void echelon_form_mt(int n, double ** matrix)
{
	int i, j;
	for (i = 0; i < n; i++) // row cycle
	{
		// Check for zero division
		if (matrix[i][i] == 0)
			if (swap_rows_ith_nonzero(n, i, matrix) == 0) continue;

		// Coefficient to divide ith row with
		double c = matrix[i][i];

		// Divide current row
		for (j = n; j >= i; j--) // element cycle
		{
			matrix[i][j] /= c;
		}

		// Subtract ith row from all below
		#pragma omp parallel for
		for (j = i + 1; j < n; j++)
		{
			// Coefficient to multiply kth row with
			double d = matrix[j][i];

			for (int k = n; k >= i; k--) // element cycle
			{
				matrix[j][k] -= d * matrix[i][k];
			}
		}
	}
}

/// <summary>
/// Runs back substitution algorithm on matrix.
/// </summary>
/// <param name="n"> Augmented matrix row count. </param>
/// <param name="matrix"> Augmented matrix. </param>
/// <param name="x"> Array to save solutions in. </param>
void back_substitution(int n, double ** matrix, double * x)
{
	x[n - 1] = matrix[n - 1][n];
	for (int i = n - 2; i >= 0; i--) // row cycle
	{
		x[i] = matrix[i][n];

		int j;
		for (j = i + 1; j < n; j++)
		{
			x[i] -= matrix[i][j] * x[j];
		}
	}
}

/// <summary>
/// Attempts to find fst row below ith, where e[j][i] != 0.
/// If e[j][i] != 0 row exists, swaps ith row with jth row, returns 1;
/// Else returns 0;
/// </summary>
/// <param name="n"> Matrix row count. </param>
/// <param name="i"> Index of row to swap with. </param>
/// <param name="matrix"> Matrix of equation. </param>
/// <returns> Flag to report if swap was performed. </returns>
static int swap_rows_ith_nonzero(int n, int i, double ** matrix)
{
	int nonzero_flag = 0;
	for (int j = i + 1; j < n; j++) // row cycle
	{
		// Check if element is non zero
		if (matrix[j][i] != 0)
		{
			// Non zero element is found, perform swap, set flag
			nonzero_flag = 1;
			double * tmp = matrix[i];
			matrix[i] = matrix[j];
			matrix[j] = tmp;
			break;
		}
	}
	return nonzero_flag;
}
