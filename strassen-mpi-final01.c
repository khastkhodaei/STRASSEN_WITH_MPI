#include <stdio.h>
#include <mpi.h>
#include <time.h>
#include <stdlib.h>


const int n = 2000;
const int M = n / 2;
int rank, size, start_row, end_row;

int A[2000][2000], B[2000][2000], C[2000][2000];
int m1[1000][1000], m2[1000][1000], m3[1000][1000], m4[1000][1000], m5[1000][1000], m6[1000][1000], m7[1000][1000];

int A11[1000][1000], A12[1000][1000], A21[1000][1000], A22[1000][1000], B11[1000][1000], B12[1000][1000], B21[1000][1000], B22[1000][1000];
int A11A22[1000][1000], B11B22[1000][1000], A21A22[1000][1000], A11A12[1000][1000], B12B22[1000][1000], B11B12[1000][1000], B21B22[1000][1000], B21B11[1000][1000], A12A22[1000][1000], A11A21[1000][1000];
int C11[1000][1000], C12[1000][1000], C21[1000][1000], C22[1000][1000];

void multiplym1(int mySize, int AA[M][M], int BB[M][M], int resfinal[M][M], int mystart_row, int myend_row)
{
    // printf("multiplym1 mystart_row: %d--------- myend_row: %d\n", mystart_row, myend_row);

    for (int i = mystart_row; i < myend_row; i++)
    {
        for (int k = 0; k < mySize; k++)
        {
            for (int j = 0; j < mySize; j += 8)
            // for (int j = 0; j < mySize; j++)
            {
                // resfinal[i][j] += (AA[i][k] * BB[k][j]);

                resfinal[i][j] += (AA[i][k] * BB[k][j]);
                resfinal[i][j + 1] += (AA[i][k] * BB[k][j + 1]);
                resfinal[i][j + 2] += (AA[i][k] * BB[k][j + 2]);
                resfinal[i][j + 3] += (AA[i][k] * BB[k][j + 3]);
                resfinal[i][j + 4] += (AA[i][k] * BB[k][j + 4]);
                resfinal[i][j + 5] += (AA[i][k] * BB[k][j + 5]);
                resfinal[i][j + 6] += (AA[i][k] * BB[k][j + 6]);
                resfinal[i][j + 7] += (AA[i][k] * BB[k][j + 7]);
            }
        }
    }
}
/*

int multiplym2(int mySize, int AA[M][M], int BB[M][M], int resfinal[M][M], int mystart_row, int myend_row)
{

    for (int i = mystart_row; i < myend_row; i++)
    {
        for (int j = 0; j < mySize; j++)
        {
            for (int k = 0; k < mySize; k++)
            {
                resfinal[i][j] += (AA[i][k] * BB[k][j]);
            }
        }
    }
}



int multiplym3(int mySize, int AA[M][M], int BB[M][M], int resfinal[M][M], int mystart_row, int myend_row)
{

    for (int i = mystart_row; i < myend_row; i++)
    {
        for (int j = 0; j < mySize; j++)
        {
            for (int k = 0; k < mySize; k++)
            {
                resfinal[i][j] += (AA[i][k] * BB[k][j]);
            }
        }
    }
}

int multiplym4(int mySize, int AA[M][M], int BB[M][M], int resfinal[M][M], int mystart_row, int myend_row)
{
    for (int i = mystart_row; i < myend_row; i++)
    {
        for (int j = 0; j < mySize; j++)
        {
            for (int k = 0; k < mySize; k++)
            {
                resfinal[i][j] += (AA[i][k] * BB[k][j]);
            }
        }
    }
}

*/

int main(int argc, char **argv)
{

    srand(time(NULL));
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            A[i][j] = 2; // rand() % 5;
            B[i][j] = 2; // rand() % 5;
            C[i][j] = 0;
        }
    }

    for (int i = 0; i < n / 2; i++)
    {
        for (int j = 0; j < n / 2; j++)
        {

            A11[i][j] = A[i][j];
            A12[i][j] = A[i][j + M];
            A21[i][j] = A[i + M][j];
            A22[i][j] = A[i + M][j + M];

            B11[i][j] = B[i][j];
            B12[i][j] = B[i][j + M];
            B21[i][j] = B[i + M][j];
            B22[i][j] = B[i + M][j + M];
        }
    }

    for (int i = 0; i < n / 2; i++)
    {
        for (int j = 0; j < n / 2; j++)
        {
            A11A22[i][j] = A11[i][j] + A22[i][j];
            B11B22[i][j] = B11[i][j] + B22[i][j];
            A21A22[i][j] = A21[i][j] + A22[i][j];
            A11A12[i][j] = A11[i][j] + A12[i][j];
            B12B22[i][j] = B12[i][j] - B22[i][j];
            B21B11[i][j] = B21[i][j] - B11[i][j];
            A11A21[i][j] = A11[i][j] - A21[i][j];
            A12A22[i][j] = A12[i][j] - A22[i][j];
            B11B12[i][j] = B11[i][j] + B12[i][j];
            B21B22[i][j] = B21[i][j] + B22[i][j];
        }
    }

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    MPI_Barrier(MPI_COMM_WORLD);

    start_row = rank * (M / size);
    if (rank + 1 == size)
    {
        end_row = M;
    }
    else
    {
        end_row = (rank + 1) * (M / size);
    }

    multiplym1(n / 2, A11A22, B11B22, m1, start_row, end_row);

    int *counts = malloc(size * sizeof(int));
    int *displs = malloc(size * sizeof(int));
    for (int i = 0; i < size; i++)
    {
        counts[i] = (M / size) * M;
        displs[i] = i * (M / size) * M;
    }
    counts[size - 1] = ((M / size) + (M % size)) * M;

    MPI_Gatherv(&m1[start_row][0], counts[rank], MPI_INT, m1, counts, displs, MPI_INT, 0, MPI_COMM_WORLD);

    // MPI_Barrier(MPI_COMM_WORLD);

    multiplym1(n / 2, A21A22, B11, m2, start_row, end_row);
    // multiplym2(n / 2, A21A22, B11, m2, start_row, end_row);

    MPI_Gatherv(&m2[start_row][0], counts[rank], MPI_INT, m2, counts, displs, MPI_INT, 0, MPI_COMM_WORLD);
    // MPI_Barrier(MPI_COMM_WORLD);

    multiplym1(n / 2, A11A12, B22, m5, start_row, end_row);

    // multiplym2(n / 2, A11A12, B22, m5, start_row, end_row);

    MPI_Gatherv(&m5[start_row][0], counts[rank], MPI_INT, m5, counts, displs, MPI_INT, 0, MPI_COMM_WORLD);
    // MPI_Barrier(MPI_COMM_WORLD);

    multiplym1(n / 2, A11, B12B22, m3, start_row, end_row);
    // multiplym3(n / 2, A11, B12B22, m3, start_row, end_row);

    MPI_Gatherv(&m3[start_row][0], counts[rank], MPI_INT, m3, counts, displs, MPI_INT, 0, MPI_COMM_WORLD);
    // MPI_Barrier(MPI_COMM_WORLD);

    multiplym1(n / 2, A22, B21B11, m4, start_row, end_row);
    // multiplym3(n / 2, A22, B21B11, m4, start_row, end_row);

    MPI_Gatherv(&m4[start_row][0], counts[rank], MPI_INT, m4, counts, displs, MPI_INT, 0, MPI_COMM_WORLD);
    // MPI_Barrier(MPI_COMM_WORLD);

    // multiplym4(n / 2, A21, A11, B11, B12, m6, start_row, end_row);

    multiplym1(n / 2, A11A21, B11B12, m6, start_row, end_row);
    // multiplym4(n / 2, A11A21, B11B12, m6, start_row, end_row);

    MPI_Gatherv(&m6[start_row][0], counts[rank], MPI_INT, m6, counts, displs, MPI_INT, 0, MPI_COMM_WORLD);
    // MPI_Barrier(MPI_COMM_WORLD);

    multiplym1(n / 2, A12A22, B21B22, m7, start_row, end_row);
    // multiplym4(n / 2, A12A22, B21B22, m7, start_row, end_row);

    MPI_Gatherv(&m7[start_row][0], counts[rank], MPI_INT, m7, counts, displs, MPI_INT, 0, MPI_COMM_WORLD);
    // MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 0)
    {
        for (int i = 0; i < M; i++)
        {
            for (int j = 0; j < M; j++)
            {
                /*
                C11[i][j] = m1[i][j] + m4[i][j] - m5[i][j] + m7[i][j];
                C12[i][j] = m3[i][j] + m5[i][j];
                C21[i][j] = m2[i][j] + m4[i][j];
                //C22[i][j] = m1[i][j] - m2[i][j] + m3[i][j] + m6[i][j];
                C22[i][j] = m1[i][j] - m2[i][j] + m3[i][j] - m6[i][j];
                //C22[i][j] = m1[i][j] + m3[i][j] - m2[i][j] + m6[i][j];
                //C22[i][j] = m1[i][j] + m2[i][j] - m3[i][j] + m6[i][j];


                -------------------------------------------------*/

                C11[i][j] = m1[i][j] + m4[i][j] - m5[i][j] + m7[i][j];
                C12[i][j] = m3[i][j] + m5[i][j];
                C21[i][j] = m2[i][j] + m4[i][j];
                C22[i][j] = m1[i][j] - m2[i][j] + m3[i][j] + m6[i][j];
            }
        }

        // Combine the four quadrants of the resulting matrix C
        for (int i = 0; i < M; i++)
        {
            for (int j = 0; j < M; j++)
            {
                C[i][j] = C11[i][j];
                C[i][j + M] = C12[i][j];
                C[i + M][j] = C21[i][j];
                C[i + M][j + M] = C22[i][j];
            }
        }

        // printf("Matrix m1,m2, m5:\n");

       

        //Print Resault of C
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                // printf("%d: %d\t", rank, m1[i][j]);
                // printf("%d: %d\t", rank, m2[i][j]);
                // printf("%d: m1[%d][%d]--%d\n", rank, i, j, m1[i][j]);
                // printf("%d: m2[%d][%d]--%d\n", rank, i, j, m2[i][j]);
                // printf("%d: m3[%d][%d]--%d\n", rank, i, j, m3[i][j]);
                // printf("%d: m4[%d][%d]--%d\n", rank, i, j, m4[i][j]);
                // printf("%d: m5[%d][%d]--%d\n", rank, i, j, m5[i][j]);
                // printf("%d: m6[%d][%d]--%d\n", rank, i, j, m6[i][j]);
                // printf("%d: m7[%d][%d]--%d\n", rank, i, j, m7[i][j]);

                printf("C[%d][%d]--%d\n", i, j, C[i][j]);
            }
            // printf("\n");
        }

        /*
printf("Matrix A:\n");
for (int i = 0; i < n; i++)
{
    for (int j = 0; j < n; j++)
    {
        printf("%d: %d\t", rank, A[i][j]);
    }
    printf("\n");
}

printf("\nMatrix B:\n");
for (int i = 0; i < n; i++)
{
    for (int j = 0; j < n; j++)
    {
        printf("%d: %d\t", rank, B[i][j]);
    }
    printf("\n");
}

printf("\nMatrix C (A * B):\n");
for (int i = 0; i < n; i++)
{
    for (int j = 0; j < n; j++)
    {
        printf("%d: %d\t", rank, C[i][j]);
    }
    printf("\n");
}
*/
    }

    MPI_Finalize();
    return 0;
}
