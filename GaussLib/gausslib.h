#pragma once

/// <summary>
/// Single thread row echelon form implementation.
/// Mutates augmented matrix (NxN+1) into row echelon form.
/// </summary>
/// <param name="n"> Augmented matrix row count. </param>
/// <param name="matrix"> Augmented matrix. </param>
void echelon_form_1t(int n, double ** matrix);

/// <summary>
/// Multi thread row echelon form implementation.
/// Mutates augmented matrix (NxN+1) into row echelon form.
/// </summary>
/// <param name="n"> Augmented matrix row count. </param>
/// <param name="matrix"> Augmented matrix. </param>
void echelon_form_mt(int n, double ** matrix);

/// <summary>
/// Runs back substitution algorithm on matrix.
/// </summary>
/// <param name="n"> Augmented matrix row count. </param>
/// <param name="matrix"> Augmented matrix. </param>
/// <param name="x"> Array to save solutions in. </param>
void back_substitution(int n, double ** matrix, double * x);
