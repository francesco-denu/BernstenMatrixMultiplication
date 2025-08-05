/* cannon.c -- uses Cannon's algorithm to multiply two (square or rectangular) matrices
 *
 * Notes:  
 *     1.  Assumes the number of processes is p^2q
 *     2.  The array member of the matrices is statically allocated
 *     3.  The size of the matrix must be divisible by the square_root(p)
 */

#include <stdio.h>
#include "mpi.h"
#include <math.h>
#include <stdlib.h>
//#include <unistd.h>

#define MAX 8192*8192 // Possiamo moltiplicare massimo matrici di cui la matrice locale ha questa dimensione (nota che dipende dal numero di processori la dimensione massima delle matrici A,B,C)

typedef unsigned int dtype; // Definisco un tipo in questo modo cosi che sia piu facile cambiarlo dopo
#define MY_MPI_TYPE MPI_UNSIGNED // Definisco anche il corrispondente MPI_TYPE

typedef struct {
    int       p;         // Total number of processes    
    MPI_Comm  comm;      // Communicator for entire grid   
    int       q;         // Order of grid                
    int       my_row;    // My row number                
    int       my_col;    // My column number             
    int       my_rank;   // My rank in the grid comm     
} GRID_INFO_T;

typedef struct {
    int     n_rows;         // Number of rows of the matrix
    int     n_cols;         // Number of cols of the matrix
    dtype   entries[MAX];   // Salviamo in maniera contigua gli elementi della matrice
} LOCAL_MATRIX_T;

#define Entry(A,i,j) (*(((A)->entries) + ((A)->n_cols)*(i) + (j))) // Data una LOCAL_MATRIX_T estrae gli elementi a partire da indici i,j

/* Function Declarations */
void             Read_matrix(char* prompt, LOCAL_MATRIX_T* local_A, GRID_INFO_T* grid, char* matrix_filename, int *n, int*m);
void             Print_matrix(char* filename, LOCAL_MATRIX_T* local_A, GRID_INFO_T* grid, int total_rows);
void             Print_local_matrix(char* title, LOCAL_MATRIX_T* local_A, GRID_INFO_T* grid);
void             print_matrix_to_file(const char *filename, LOCAL_MATRIX_T *matrix);

void             Set_to_zero(LOCAL_MATRIX_T* local_A);

void             Build_matrix_type(LOCAL_MATRIX_T* local_A);

void             Setup_grid(GRID_INFO_T*  grid);

void             Local_matrix_multiply(LOCAL_MATRIX_T* local_A, LOCAL_MATRIX_T* local_B, LOCAL_MATRIX_T* local_C);
void             Cannon_Setup(GRID_INFO_T* grid, LOCAL_MATRIX_T* local_A, LOCAL_MATRIX_T* local_B);
void             Cannon_Multiplication(int n, GRID_INFO_T* grid, LOCAL_MATRIX_T* local_A, LOCAL_MATRIX_T* local_B, LOCAL_MATRIX_T* local_C);

MPI_Datatype     local_matrix_mpi_t;

LOCAL_MATRIX_T*  temp_mat;

/*********************************************************/

int main(int argc, char* argv[]) {
    int              my_rank;
    GRID_INFO_T      grid;
    LOCAL_MATRIX_T*  local_A;
    LOCAL_MATRIX_T*  local_B;
    LOCAL_MATRIX_T*  local_C;
    int              n;
    int              n_A, n_B, m_A, m_B;
    double seconds; // Memorizza il tempo per eseguire l'algoritmo su ogni processo
    double total_seconds; // Memorizza il tempo del processore che ha lavorato per più tempo


    if (argc != 4) {
        fprintf(stderr, "Usage: %s <matrix1file> <matrix2file> <outputfile>\n", argv[0]);
        MPI_Finalize();
        return 1;
    }

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    Setup_grid(&grid);
    
    local_A = (LOCAL_MATRIX_T*) malloc(sizeof(LOCAL_MATRIX_T));
    local_B = (LOCAL_MATRIX_T*) malloc(sizeof(LOCAL_MATRIX_T));
    local_C = (LOCAL_MATRIX_T*) malloc(sizeof(LOCAL_MATRIX_T));
    
    Read_matrix("Reading A...", local_A, &grid, argv[1], &n_A, &m_A);
    //Print_matrix("matrixA", local_A, &grid, n_A);

    Read_matrix("Reading B...", local_B, &grid, argv[2], &n_B, &m_B);
    //Print_matrix("matrixB", local_B, &grid, n_B);

    local_C->n_rows = local_A->n_rows;
    local_C->n_cols = local_B->n_cols;
    Set_to_zero(local_C);

    Build_matrix_type(local_A);

    MPI_Barrier (MPI_COMM_WORLD); // Da qui cominciamo a prendere i tempi
    seconds = - MPI_Wtime();

    Cannon_Setup(&grid, local_A, local_B);

    Cannon_Multiplication(n, &grid, local_A, local_B, local_C);

    MPI_Barrier (MPI_COMM_WORLD); // Fin qui prendiamo i tempi
    seconds += MPI_Wtime();

    MPI_Allreduce (&seconds, &total_seconds, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    if (!grid.my_rank) {
        printf ("CANNON Algorithm N = %d M = %d, Processes = %d, Time = %12.6f sec\n", n_A, m_A, grid.p, total_seconds);
        //printf ("%f \n", total_seconds);

    }

    Print_matrix(argv[3], local_C, &grid, n_A);
    
    free(local_A);
    free(local_B);
    free(local_C);

    free(temp_mat);

    MPI_Finalize();

    return 0;
}  /* main */


/*********************************************************/
void Setup_grid(GRID_INFO_T*  grid  /* out */) {
    int dimensions[2];
    int wrap_around[2];
    int coordinates[2];
    //int free_coords[2];

    MPI_Comm_size(MPI_COMM_WORLD, &(grid->p)); // Estraggo numero di processori nel comunicatore generale

    grid->q = (int) sqrt((double) grid->p); // Supponiamo p sia quadrato perfetto
    if (grid->q * grid->q != grid->p){
            perror("p is not perfect square");
            printf("p is not perfect square");
            MPI_Abort(grid->comm, 1);
            MPI_Finalize();
            return;
    }

    dimensions[0] = dimensions[1] = grid->q;
    wrap_around[0] = wrap_around[1] = 1; // Vogliamo circular shift 
    
    MPI_Cart_create(MPI_COMM_WORLD, 2, dimensions, wrap_around, 1, &(grid->comm));

    MPI_Comm_rank(grid->comm, &(grid->my_rank));
    MPI_Cart_coords(grid->comm, grid->my_rank, 2, coordinates);
    grid->my_row = coordinates[0];
    grid->my_col = coordinates[1];

} /* Setup_grid */


/*********************************************************/
/* Read and distribute matrix:  
 *     foreach global row of the matrix,
 *         foreach cpu in the grid 
 *             read a block of row of the corresponding cpu
 *             and send them to the appropriate process.
 */
void Read_matrix(char* prompt /* in  */, LOCAL_MATRIX_T* local_A /* out */, GRID_INFO_T* grid /* in  */, char* matrix_filename, int* n, int* m) {
    int        mat_row, mat_col;
    int        grid_row, grid_col;
    int        dest;
    int        coords[2];
    dtype*     temp;
    MPI_Status status;
    FILE *finptr;

    if (grid->my_rank == 0) {

        finptr = fopen(matrix_filename, "r");

        if (finptr == NULL) {
            perror("Error opening file");
            printf("Matrix file not found");
            MPI_Abort(grid->comm, 1);
            MPI_Finalize();

        }        
        
        fscanf(finptr, "%d %d", n, m);
    }

        MPI_Bcast(n, 1, MPI_INT, 0, MPI_COMM_WORLD); // Invio in broadcast l'ordine delle matrici A,B,C
        MPI_Bcast(m, 1, MPI_INT, 0, MPI_COMM_WORLD); // Invio in broadcast l'ordine delle matrici A,B,C

        // Supposing the number of processes for each submatrix is perfect square 2^q with q even
        local_A->n_rows = *n / grid->q;
        local_A->n_cols = *m / grid->q;


    if (grid->my_rank == 0) {
        //printf("%s\n", prompt);
        temp = (dtype*) malloc((local_A->n_cols)*sizeof(dtype));
        fflush(stdout);
        
        for (mat_row = 0;  mat_row < *n; mat_row++) {
            grid_row = mat_row/local_A->n_rows;
            coords[0] = grid_row;
            for (grid_col = 0; grid_col < grid->q; grid_col++) {
                coords[1] = grid_col;
                MPI_Cart_rank(grid->comm, coords, &dest);
                for (mat_col = 0; mat_col < local_A->n_cols; mat_col++){
                    if (dest == 0){
                        fscanf(finptr ,"%d", &Entry(local_A, mat_row, mat_col)); // Siamo nel caso del processo root, quindi mat_row coincide con la riga di local_A
                    }else{
                        fscanf(finptr ,"%d", temp + mat_col);
                    }
                }
                if (dest != 0){
                MPI_Send(temp, local_A->n_cols, MY_MPI_TYPE, dest, 0, grid->comm);
                }
            }
        }
        free(temp);
        fclose(finptr);

    } else {
        for (mat_row = 0; mat_row < local_A->n_rows; mat_row++){
           MPI_Recv(&Entry(local_A, mat_row, 0), local_A->n_cols, MY_MPI_TYPE, 0, 0, grid->comm, &status);
       }
    }

}  /* Read_matrix */

void Print_matrix(char* filename /* in  */, LOCAL_MATRIX_T* local_A /* out */, GRID_INFO_T* grid /* in  */, int total_rows) {
    int        mat_row, mat_col;
    int        grid_row, grid_col;
    int        source;
    int        coords[2];
    dtype*     temp;
    MPI_Status status;

    if (grid->my_rank == 0) {
        
        FILE *file = fopen(filename, "w+");
        if (file == NULL) {
            perror("Error opening file");
            return;
        }

        fprintf(file, "%d %d\n", local_A->n_rows * grid->q, local_A->n_cols* grid->q);

        temp = (dtype*) malloc(local_A->n_rows*local_A->n_cols);
        for (mat_row = 0;  mat_row < total_rows; mat_row++) {
            grid_row = mat_row/local_A->n_rows;
            coords[0] = grid_row;
            for (grid_col = 0; grid_col < grid->q; grid_col++) {
                coords[1] = grid_col;
                MPI_Cart_rank(grid->comm, coords, &source);
                if (source == 0) {
                    for(mat_col = 0; mat_col < local_A->n_cols; mat_col++)
                        fprintf(file, "%d ", Entry(local_A, mat_row, mat_col));
                } else {
                    MPI_Recv(temp, local_A->n_cols, MY_MPI_TYPE, source, 0, grid->comm, &status);
                    for(mat_col = 0; mat_col < local_A->n_cols; mat_col++)
                        fprintf(file, "%d ", temp[mat_col]);
                }
            }
            fprintf(file,"\n");
        }
        free(temp);
        fclose(file);
    } else {
        for (mat_row = 0; mat_row < local_A->n_rows; mat_row++) 
            MPI_Send(&Entry(local_A, mat_row, 0), local_A->n_cols, MY_MPI_TYPE, 0, 0, grid->comm);
    }
                     
}  /* Print_matrix */


/*********************************************************/
void Set_to_zero(LOCAL_MATRIX_T*  local_A  /* out */) {
    int i, j;
    for (i = 0; i < local_A->n_rows; i++)
        for (j = 0; j < local_A->n_cols; j++)
            Entry(local_A,i,j) = 0.0;
}  /* Set_to_zero */


/*********************************************************/

void Build_matrix_type(LOCAL_MATRIX_T*  local_A) {
    MPI_Datatype  temp_mpi_t;
    int           block_lengths[3];
    MPI_Aint      displacements[3];
    MPI_Datatype  typelist[3];
    MPI_Aint      start_address;
    MPI_Aint      address;

    // Crea un MPI_DATATYPE (temp_mpi_t) che corrisponde a un blocco contiguo di FLOAT 
    MPI_Type_contiguous(local_A->n_rows*local_A->n_cols, MY_MPI_TYPE, &temp_mpi_t);

    MPI_Get_address(local_A, &start_address); // Indirizzo di inizio local_A (struttura di tipo LOCAL_MATRIX_T)

    MPI_Get_address(&(local_A->n_rows), &address); // Indirizzo di n_bar all'interno di LOCAL_MATRIX_T
    displacements[0] = address - start_address; // Dopo quanti byte rispetto all'inizio della struttrura c'è n_bar

    MPI_Get_address(&(local_A->n_cols), &address); // Indirizzo di n_bar all'interno di LOCAL_MATRIX_T
    displacements[1] = address - start_address; // Dopo quanti byte rispetto all'inizio della struttrura c'è n_bar

    MPI_Get_address(local_A->entries, &address); // Indirizzo di entries all'interno di LOCAL_MATRIX_T
    displacements[2] = address - start_address; // Dopo quanti byte rispetto all'inizio della struttrura c'è entries

    typelist[0] = MPI_INT;
    typelist[1] = MPI_INT;
    typelist[2] = temp_mpi_t;
    block_lengths[0] = block_lengths[1] = block_lengths[2] = 1; // Specifica che ci sarà 1 intero (n_bar) e 1 array (entries)

    // Creiamo adesso un MPI_Datatype per la matrice, che ha 2 elementi, ciascuno con 1 elemento (block_lengths), con i displacement in displacements, e di tipo in typelist (quindi un MPI_INT e un temp_mpi_t costruito poco fa, ovvero il blocco contiguo di float) 
    MPI_Type_create_struct(3, block_lengths, displacements, typelist, &local_matrix_mpi_t);
    MPI_Type_commit(&local_matrix_mpi_t); // Rendiamo il tipo local_matrix_mpi_t utilizzabile nelle comunicazioni
} 


/*********************************************************/

void Local_matrix_multiply(LOCAL_MATRIX_T*  local_A , LOCAL_MATRIX_T*  local_B, LOCAL_MATRIX_T*  local_C) {
    int i, j, k;

    for (i = 0; i < local_C->n_rows; i++)
        for (j = 0; j < local_C->n_cols; j++)
            for (k = 0; k < local_A->n_cols; k++)
                Entry(local_C,i,j) = Entry(local_C,i,j) + Entry(local_A,i,k)*Entry(local_B,k,j);

}  


/*********************************************************/
void Print_local_matrix(char* title /* in */, LOCAL_MATRIX_T* local_A  /* in */, GRID_INFO_T* grid /* in */) {
    
    printf("Process %d > grid_row = %d, grid_col = %d\n", grid->my_rank, grid->my_row, grid->my_col);
    for (int i = 0; i <local_A->n_rows ; i++) {
        for (int j = 0; j < local_A->n_cols; j++){
            printf("%d ", Entry(local_A,i,j));
        }
        printf("\n");
    }
    
}  /* Print_local_matrix */


void Cannon_Setup(GRID_INFO_T* grid, LOCAL_MATRIX_T* local_A, LOCAL_MATRIX_T* local_B){
    MPI_Status status;
    int source_rank, destination_rank;
    int source_position[2];
    int destination_position[2]; 

    // A
    MPI_Cart_shift(grid->comm, 1, -grid->my_row, &source_rank, &destination_rank); // sulla base della mia posizione trovo il rank che è nella posizione in cui devo andare
    MPI_Cart_coords(grid->comm, source_rank, 2, source_position);
    MPI_Cart_coords(grid->comm, destination_rank, 2, destination_position);
    //printf("Io sono in A in %d %d e mi sposto di %d in %d %d mentre al posto mio viene %d %d\n", grid->my_row, grid->my_col, -grid->my_row, destination_position[0], destination_position[1], source_position[0], source_position[1]);
    MPI_Sendrecv_replace(local_A, 1, local_matrix_mpi_t, destination_rank, 0, source_rank, 0, grid->comm, &status);

    // B
    MPI_Cart_shift(grid->comm, 0, -grid->my_col, &source_rank, &destination_rank); // sulla base della mia posizione trovo il rank che è nella posizione in cui devo andare
    MPI_Cart_coords(grid->comm, source_rank, 2, source_position);
    MPI_Cart_coords(grid->comm, destination_rank, 2, destination_position);
    //printf("Io sono in B in %d %d e mi sposto di %d in %d %d mentre al posto mio viene %d %d\n", grid->my_row, grid->my_col, -grid->my_row, destination_position[0], destination_position[1], source_position[0], source_position[1]);
    MPI_Sendrecv_replace(local_B, 1, local_matrix_mpi_t, destination_rank, 0, source_rank, 0, grid->comm, &status);
}

void Cannon_Multiplication(int n, GRID_INFO_T* grid, LOCAL_MATRIX_T* local_A, LOCAL_MATRIX_T* local_B, LOCAL_MATRIX_T* local_C){
    MPI_Status status;
    int source_rank, destination_rank;

    for (int i=0; i < grid->q - 1; i++){
        Local_matrix_multiply(local_A, local_B, local_C);

        MPI_Cart_shift(grid->comm, 1, -1, &source_rank, &destination_rank); // sulla base della mia posizione trovo il rank che è nella posizione in cui devo andare
        MPI_Sendrecv_replace(local_A, 1, local_matrix_mpi_t, destination_rank, 0, source_rank, 0, grid->comm, &status);

        MPI_Cart_shift(grid->comm, 0, -1, &source_rank, &destination_rank); // sulla base della mia posizione trovo il rank che è nella posizione in cui devo andare
        MPI_Sendrecv_replace(local_B, 1, local_matrix_mpi_t, destination_rank, 0, source_rank, 0, grid->comm, &status);
    }

    Local_matrix_multiply(local_A, local_B, local_C);


}
