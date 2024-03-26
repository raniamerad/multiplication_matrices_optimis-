#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

void mat_alloc(float *** A, int rows, int cols)
{
        int i;

        *A = malloc(sizeof(float *) * rows);
        for (i = 0 ; i < cols ; ++i)
        {
                (*((*A)+i)) = malloc(sizeof(float) * cols);
                //(*A)[i] = malloc(sizeof(float) * cols);
        }
}

void mat_dealloc(float ** A)
{
        free(A);
}

void mat_populate_val(float ** A, int rows, int cols, float val)
{
        int i, j;

        for(i = 0 ; i < rows ; ++i)
        {
                for(j = 0 ; j < cols ; ++j)
                {
                        *(*(A+i)+j) = val;
                }
        }
}

void mat_populate_rand(float ** A, int rows, int cols)
{
        int i, j;

        srand((unsigned int) time(NULL));

        for(i = 0 ; i < rows ; ++i)
        {
                for(j = 0 ; j < cols ; ++j)
                {
                        *(*(A+i)+j) = (float)rand()/RAND_MAX;
                }
        }
}

void mat_display(float ** A, int rows, int cols)
{
        int i, j;

        for(i = 0 ; i < rows ; ++i)
        {
                for(j = 0 ; j < cols ; ++j)
                {
                        printf("[%d][%d] = ", i, j);
                        printf("%.6f \n", A[i][j]);
                }
                printf("\n");
        }
}

float ** mat_add(float ** A, float ** B, int rows, int cols)
{
        int i, j, ii;
        float tmp;

        float ** R;
        mat_alloc(&R, rows, rows);

        for(i = 0 ; i < rows ; ++i)
        {
                for(ii = 0 ; ii < rows ; ++ii)
                {
                        R[i][ii] = 0.0;
                        for(j = 0 ; j < cols ; ++j)
                        {
                        	R[i][ii] = A[ii][j] + B[j][ii];
                       	}
             	}
                        
        }
       	return R;
}

float ** mat_sub(float ** A, float ** B, int rows, int cols)
{
        int i, j, ii;
        float tmp;

        float ** R;
        mat_alloc(&R, rows, rows);

        for(i = 0 ; i < rows ; ++i)
        {
                for(ii = 0 ; ii < rows ; ++ii)
                {
                        R[i][ii] = 0.0;
                        for(j = 0 ; j < cols ; ++j)
                        {
                        	R[i][ii] = A[ii][j] - B[j][ii];
                       	}
             	}
                        
        }
       	return R;
}
                                              
float ** mat_mul(float ** A, float ** B, int rows, int cols)
{
        int i, j, ii;
        float tmp;

        float ** R;
        mat_alloc(&R, rows, rows);

        for(i = 0 ; i < rows ; ++i)
        {
                for(ii = 0 ; ii < rows ; ++ii)
                {
                        R[i][ii] = 0.0;
                        for(j = 0 ; j < cols ; ++j)
                        {
                                R[i][ii] += A[ii][j] * B[j][ii];
                        }
                }
        }

        return R;
}

// mat_split : splits A into four submatrices A0;A1;A2;A3
void mat_split(float ** A,
	float ** A0, float ** A1, float ** A2, float ** A3,
	int sub_rows)
{
	int i, j;

	for (i = 0 ; i < sub_rows ; i++)
	{
		for (j = 0 ; j < sub_rows ; j++)
		{
			A0[i][j] = A[i][j];
			A1[i][j] = A[i+sub_rows][j];
			A2[i][j] = A[i][j+sub_rows];
			A3[i][j] = A[i+sub_rows][j+sub_rows];
		}	
	}
}

// mat_copy : copies a subset of A in B
void mat_copy(float ** A, float ** B, int x, int y, int len)
{
	int i, j;
	for (i = 0 ; i < len ; i++)
	{
		for (j = 0 ; j < len ; j++)
		{
			B[x+i][y+j] = A[i][j];
		}
		
	}
}

float ** mat_multiply_strassen_elements(float ** A, float ** B, int rows)
{
        int i, j, k;
        float m1, m2, m3, m4, m5, m6, m7;

        float ** C;
        mat_alloc(&C, rows, rows);

	// 1 - sub-matrix 2x2 - divide-and-conquer
	for (i = 0 ; i < rows ; i+=2)
	{
		for (j = 0 ; j < rows ; j+=2)
		{
			// Strassen kernels
			m1 = (A[i][j]+A[i+1][j+1])*(B[i][j]+B[i+1][j+1]);
			m2 = (A[i+1][j]+A[i+1][j+1])*B[i][j];
			m3 = A[i][j]*(B[i][j+1]-B[i+1][j+1]);
			m4 = A[i+1][j+1]*(B[i+1][j]-B[i][j]);
			m5 = (A[i][j]+A[i][j+1])*B[i+1][j+1];
			m6 = (A[i+1][j]-A[i][j])*(B[i][j]+B[i][j+1]);
			m7 = (A[i][j+1]-A[i+1][j+1])*(B[i+1][j]+B[i+1][j+1]);  
			C[i][j] = m1 + m4 - m5 + m7;
			C[i+1][j] = m2 + m4;
			C[i][j+1] = m3 + m5;
			C[i+1][j+1] = m1 - m2 + m3 + m6;
		}	
	}
        return C;
}

float ** mat_multiply_strassen_matrices(float ** A, float ** B, int rows)
{
        float ** C;
	    mat_alloc(&C, rows, rows);	

        if (rows > 4) 
        {
                int i, j, k, sub_row;
                sub_row = rows/2;
                float ** m1;
                float ** m2;
                float ** m3;
                float ** m4;
                float ** m5;
                float ** m6;
                float ** m7;

                float ** mA0;
                float ** mA1;
                float ** mA2;
                float ** mA3;

                float ** mB0;
                float ** mB1;
                float ** mB2;
                float ** mB3;

                float ** mC0;
                float ** mC1;
                float ** mC2;
                float ** mC3;

                mat_alloc(&m1, sub_row, sub_row);
                mat_alloc(&m2, sub_row, sub_row);
                mat_alloc(&m3, sub_row, sub_row);
                mat_alloc(&m4, sub_row, sub_row);
                mat_alloc(&m5, sub_row, sub_row);
                mat_alloc(&m6, sub_row, sub_row);
                mat_alloc(&m7, sub_row, sub_row);

                mat_alloc(&mC0, sub_row, sub_row);
                mat_alloc(&mC1, sub_row, sub_row);
                mat_alloc(&mC2, sub_row, sub_row);
                mat_alloc(&mC3, sub_row, sub_row);

                mat_alloc(&mA0, sub_row, sub_row);
                mat_alloc(&mA1, sub_row, sub_row);
                mat_alloc(&mA2, sub_row, sub_row);
                mat_alloc(&mA3, sub_row, sub_row);

                mat_alloc(&mB0, sub_row, sub_row);
                mat_alloc(&mB1, sub_row, sub_row);
                mat_alloc(&mB2, sub_row, sub_row);
                mat_alloc(&mB3, sub_row, sub_row);

                mat_split(A, mA0, mA1, mA2, mA3, sub_row);
                mat_split(B, mB0, mB1, mB2, mB3, sub_row);

                m1 = mat_multiply_strassen_matrices(mat_add(mA0,mA3,sub_row,sub_row), mat_add(mB0,mB2,sub_row,sub_row), sub_row);
	            m2 = mat_multiply_strassen_matrices(mat_add(mA1,mA3,sub_row,sub_row), mB0, sub_row);
	            m3 = mat_multiply_strassen_matrices(mA0, mat_sub(mB1, mB3,sub_row,sub_row), sub_row);
	            m4 = mat_multiply_strassen_matrices(mA3, mat_sub(mB1, mB0,sub_row,sub_row), sub_row);
	            m5 = mat_multiply_strassen_matrices(mat_add(mA0,mA2,sub_row,sub_row), mB3, sub_row);
	            m6 = mat_multiply_strassen_matrices(mat_sub(mA1, mA0,sub_row,sub_row),mat_add(mB0, mB2,sub_row,sub_row), sub_row);
	            m7 = mat_multiply_strassen_matrices(mat_sub(mA2, mA3,sub_row,sub_row),mat_add(mB2, mB3,sub_row,sub_row), sub_row);	

	            mC0 = mat_add(mat_sub(mat_add(m1, m4, sub_row, sub_row), m5, sub_row, sub_row), m7, sub_row, sub_row);
	            mC1 = mat_add(m2, m4, sub_row, sub_row);
	            mC2 = mat_add(m3, m5, sub_row, sub_row);
	            mC3 = mat_add(mat_add(mat_sub(m1, m2, sub_row, sub_row), m3, sub_row, sub_row), m6, sub_row, sub_row);
	
	            mat_copy(mC0, C, 0, 0, sub_row);
	            mat_copy(mC1, C, sub_row, 0, sub_row);
	            mat_copy(mC2, C, 0, sub_row, sub_row);
	            mat_copy(mC3, C, sub_row, sub_row, sub_row);	

                return C;
        }
        else
        {
                C = mat_mul(A, B, rows, rows);

                return C;
        }
}


float ** algo_strassen(float ** A, float ** B, int rows)
{
	//1 - mult matrices A and B @t highest resolution (256x256 matrices of size 2x2)
	//2 - downgrade resolution by 4 
	//	(e.g. 128x128 matrices of size 2x2 -> 128x128 matrices of size 4x4)
	//3 - repeat until resolution is a 1x1 matrices of size 512x512 

	float ** C;
    mat_alloc(&C, rows, rows);

	C = mat_multiply_strassen_matrices(A, B, rows);

	return C;
}

int main(int argc, char ** argv)
{	
        float ** matA;
        float ** matB;
        int rA, cA, rB, cB;
        rA = 32;
        cA = rA;
        rB = cA;
        cB = rA;

        mat_alloc(&matA, rA, cA);
        mat_populate_val(matA, rA, cA, 1.0);
        //mat_populate_rand(matA, rA, cA);
        //printf("matrix A :\n");
        //mat_display(matA, rA, cA);

        mat_alloc(&matB, rB, cB);
        //mat_populate_rand(matB, rB, cB);
        mat_populate_val(matB, rB, cB, 1.0);
        //printf("matrix B :\n");
        //mat_display(matB, rB, cB);

        float ** matR;
	    mat_alloc(&matR, rA, cA);
	    matR = algo_strassen(matA, matB, rA);
        //matR = mat_mul(matA, matB, rA, cA);

        printf("matrix R :\n");
       	mat_display(matR, rA, cA);

        mat_dealloc(matA);
        mat_dealloc(matB);
        mat_dealloc(matR);

        return 0;
}
