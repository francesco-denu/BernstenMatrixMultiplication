/* bernsten.c -- uses Bernsten's algorithm to multiply two matrices
 *
 * Notes:  
 *     1.  Assumes the number of processes is p^3*i
 *     2.  The array member of the matrices is statically allocated
 *     3.  The size of the matrices must be divisible by p^1/3, and the size of the submatrices A(i)/B(i) must be divisible by p^1/3 in order to apply cannon
 */

#include <stdio.h>
#include "mpi.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define MAX 8192*8192 // Possiamo moltiplicare massimo matrici di cui la matrice locale ha questa dimensione (nota che dipende dal numero di processori la dimensione massima delle matrici A,B,C)

// Definisco un tipo, e il corrispondente MPI_TYPE, in questo modo in modo che sia immediato cambiarlo dopo se necessario
typedef unsigned int dtype; 
#define MY_MPI_TYPE MPI_UNSIGNED 

// La seguente struttura dati memorizza tutte le informazioni relative alla topologia dei processi
// full_comm è un comunicatore cartesiano di dimensione q*q*h, che viene diviso in h comunicatori cartesiani bidimensionali q*q memorizzati in grid_comm
typedef struct {
    int       p;            // Total number of processes    
    MPI_Comm  full_comm;    // 3D communicator for all the processes
    MPI_Comm  grid_comm;    // Communicator for grid where each proces belongs
    MPI_Comm  vertical_comm;// Verical communicator with processes of fixed i,j and free z
    int       q;            // Order of bidimensional grid (q*q)
    int       h;            // height of the 3d communicator (2^k)         
    int       my_row;       // Row of the current process       
    int       my_col;       // Column of the current process
    int       my_z;         // Height of the current process
    int       full_rank;    // Rank of the current process in the full_comm communicator   
} GRID_INFO_T;

// La seguente struttura definisce una matrice e ne memorizza i dati in maniera contigua
typedef struct {
    int     n_rows;         // Number of rows of the matrix
    int     n_cols;         // Number of cols of the matrix
    dtype   entries[MAX];  // The actual contiguous matrix
} LOCAL_MATRIX_T;

#define Entry(A,i,j) (*(((A)->entries) + ((A)->n_cols)*(i) + (j))) // Data una LOCAL_MATRIX_T estrae gli elementi a partire da indici i,j
MPI_Datatype     local_matrix_mpi_t; // MPI type per LOCAL_MATRIX_T

/* Function Declarations */
void            Read_matrix_by_columns(LOCAL_MATRIX_T* full_A , GRID_INFO_T* grid , char* matrix_filename, int* n, int* m);
void            Print_matrix_by_column(char* filename, LOCAL_MATRIX_T* local_A, GRID_INFO_T* full_grid , int total_rows) ;
void            Read_matrix_by_rows(LOCAL_MATRIX_T* full_B , GRID_INFO_T* grid , char* matrix_filename, int* n, int* m);
void            Print_matrix_by_rows(char* filename, LOCAL_MATRIX_T* local_B, GRID_INFO_T* full_grid , int total_rows);
void            Distribute_local_matrices(LOCAL_MATRIX_T* full_A, LOCAL_MATRIX_T* local_A , GRID_INFO_T* grid);
void            Print_level_matrix(char* filename, LOCAL_MATRIX_T* local_A, GRID_INFO_T* grid);
void            Print_local_matrix(char* title, LOCAL_MATRIX_T* local_A, GRID_INFO_T* grid);

void            Setup_communicator(GRID_INFO_T*  comm);
void            Build_matrix_type(LOCAL_MATRIX_T* local_A);

void            Local_matrix_multiply(LOCAL_MATRIX_T* local_A, LOCAL_MATRIX_T* local_B, LOCAL_MATRIX_T* local_C);
void            Cannon_Setup(GRID_INFO_T* grid, LOCAL_MATRIX_T* local_A, LOCAL_MATRIX_T* local_B);
void            Cannon_Multiplication(GRID_INFO_T* grid, LOCAL_MATRIX_T* local_A, LOCAL_MATRIX_T* local_B, LOCAL_MATRIX_T* local_C);

void            Set_to_zero(LOCAL_MATRIX_T* local_A);


/*********************************************************/

int main(int argc, char* argv[]) {

    // Informazioni sui comunicatori
    GRID_INFO_T full_comm;

    // Le matrici A/B vengono divise in 2^k colonne/righe, e un processore per livello memorizzerà A(i) e B(i), per poi distribuire le varie sottomatrici (di A(i), B(i)) a tutti gli altri processi dello stesso livello 
    // Seguono quindi le matrici dove si memorizzeranno A(i), B(i)
    LOCAL_MATRIX_T* full_A; 
    LOCAL_MATRIX_T* full_B; 
    
    // Ogni processore moltiplicherà local_A * local_B per ottenere local_C.
    LOCAL_MATRIX_T* local_A; 
    LOCAL_MATRIX_T* local_B;
    LOCAL_MATRIX_T* local_C;

    // I local_C corrispondenti faranno una reduction per ottenere Cln definitivo 
    LOCAL_MATRIX_T* Cln; 

    int n_A, n_B, m_A, m_B; // Dimensioni di A,B

    double seconds; // Memorizza il tempo per eseguire l'algoritmo su ogni processo
    double total_seconds; // Memorizza il tempo del processore che ha lavorato per più tempo

    if (argc != 4) {
        fprintf(stderr, "Usage: %s <matrix_A_file> <matrix_B_file> <output_file>\n", argv[0]);
        MPI_Finalize();
        return -1;
    }

    MPI_Init(&argc, &argv);

    Setup_communicator(&full_comm);

    full_A = (LOCAL_MATRIX_T*) malloc(sizeof(LOCAL_MATRIX_T));
    full_B = (LOCAL_MATRIX_T*) malloc(sizeof(LOCAL_MATRIX_T));
    local_A = (LOCAL_MATRIX_T*) malloc(sizeof(LOCAL_MATRIX_T));
    local_B = (LOCAL_MATRIX_T*) malloc(sizeof(LOCAL_MATRIX_T));
    local_C = (LOCAL_MATRIX_T*) malloc(sizeof(LOCAL_MATRIX_T));
    Cln = (LOCAL_MATRIX_T*) malloc(sizeof(LOCAL_MATRIX_T));

    // Memorizziamo A(i), B(i)
    Read_matrix_by_columns(full_A, &full_comm, argv[1], &n_A, &m_A);
    //Print_matrix_by_column("testA", full_A, &full_comm, n_A);
    Read_matrix_by_rows(full_B, &full_comm, argv[2], &n_B, &m_B);
    //Print_matrix_by_rows("testB", full_B, &full_comm, n_B);

    // Ogni A(i), B(i) viene successivamente diviso in slice per effettuare cannon
    Distribute_local_matrices(full_A, local_A, &full_comm);
    //Print_level_matrix("Ax", local_A, &full_comm);
    Distribute_local_matrices(full_B, local_B, &full_comm);
    //Print_level_matrix("Bx", local_B, &full_comm);

    free(full_A);
    free(full_B);

    local_C->n_rows = local_A->n_rows;
    local_C->n_cols = local_B->n_cols;
    Cln->n_rows = local_C->n_rows;
    Cln->n_cols = local_C->n_cols;

    Set_to_zero(local_C);
    Set_to_zero(Cln);

    Build_matrix_type(local_A); // Costruisco un MPI type per la matrice local_matrix

    MPI_Barrier (MPI_COMM_WORLD); // Da qui cominciamo a prendere i tempi
    seconds = - MPI_Wtime();

    Cannon_Setup(&full_comm, local_A, local_B);
    //Print_level_matrix("Ax_after_setup", local_A, &full_comm);
    //Print_level_matrix("Ax_after_setup", local_B, &full_comm);

    Cannon_Multiplication(&full_comm, local_A, local_B, local_C);
    //Print_level_matrix("Cx", local_C, &full_comm);

    // Reduction per sommare i vari local_C corrispondenti di ogni livello in modo da ottenere Cln finale
    MPI_Reduce(local_C->entries, Cln->entries, local_C->n_rows*local_C->n_cols, MY_MPI_TYPE, MPI_SUM, 0, full_comm.vertical_comm);
    //MPI_Allreduce(local_C->entries, Cln->entries, local_C->n_rows*local_C->n_cols, MY_MPI_TYPE, MPI_SUM, full_comm.vertical_comm); // Se vogliamo la all reduce invece della reduce

    MPI_Barrier (MPI_COMM_WORLD); // Fin qui prendiamo i tempi
    seconds += MPI_Wtime();

    MPI_Allreduce (&seconds, &total_seconds, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    if (!full_comm.full_rank) {
      printf ("BERNSTEN Algorithm N = %d, Processes = %d, Time = %12.6f sec\n", n_A, full_comm.p, total_seconds);
      //printf ("%f\n", total_seconds);
    }

    if (full_comm.my_z == 0) // Se sei nel livello 0, avrai Cln finale, quindi stampiamo matrice finale
        Print_level_matrix(argv[3], Cln, &full_comm);

    // NOTA: potrei anche fare la free di local_A e local_B poco prima
    free(local_A); 
    free(local_B);
    free(local_C);
    free(Cln);

    MPI_Type_free(&local_matrix_mpi_t);

    MPI_Finalize();

    return 0;
}  /* main */


/*********************************************************/
/*
La seguente funzione crea:
    - comunicatore cartesiano 3d
    - un comunicatore 2d per ogni livello 
    - un comunicatore unidimensionale verticale che verrà usato per la reduction dei vari Cln dei vari livelli 
*/
void Setup_communicator(GRID_INFO_T*  comm  /* out */) {
    int dimensions[3];
    int wrap_around[3];
    int coordinates[3];

    MPI_Comm_size(MPI_COMM_WORLD, &(comm->p)); // Estraggo numero di processori totali

    int n_parts = cbrt(comm->p); // Il numero di parti in cui si dividono le matrici A,B (ovvero 2^k) è la radice cubica di p
    if (n_parts * n_parts * n_parts != comm->p){
        printf("Il numero di processori non è adatto, non posso calcolarne radice cubica. \n");
        MPI_Abort(MPI_COMM_WORLD, 1);
        return;
    }

    comm->h = n_parts; // L'altezza del comunicatore tridimensionale sarà quante sono le parti Ai/Bi, ogni sottogriglia bidimensionale infatti si occuperà di A(i)*B(i)
    comm->q = (int) sqrt(comm->p / n_parts); //  Il numero di processori per ogni "parte" è p^2/3, ovvero p / n_parti. Ogni griglia di processori sarà allora qxq, dove q = radice(p/n_parti)

    // Il comunicatore tridimensionale sarà q*q*n_parts: ogni "parte" calcolerà A(i)*B(i) tramite cannnon
    dimensions[0] = dimensions[1] = comm->q;
    dimensions[2] = n_parts;
    wrap_around[0] = wrap_around[1] = 1; // Vogliamo circular shift 
    wrap_around[2] = 0; 

    MPI_Cart_create(MPI_COMM_WORLD, 3, dimensions, wrap_around, 1, &(comm->full_comm));

    MPI_Comm_rank(comm->full_comm, &(comm->full_rank)); // Estraggo il rank del processore corrente nel comunicatore generale
    MPI_Cart_coords(comm->full_comm, comm->full_rank, 3, coordinates); // Estraggo le coordinate del processo all'interno del comunicatore generale
    comm->my_row = coordinates[0];
    comm->my_col = coordinates[1];
    comm->my_z = coordinates[2];

    // Creo un comunicatore bidimensionale a griglia per ogni "parte" (che eseguirà cannon su ogni livello)
    int grid_dimensions[3] = {1,1,0};
    MPI_Cart_sub(comm->full_comm, grid_dimensions, &(comm->grid_comm));

    // Creo un comunicatore monodimensionale sull'altezza che useremo alla fine per effettuare la reduction
    grid_dimensions[0] = grid_dimensions[1] = 0;
    grid_dimensions[2] = 1;
    MPI_Cart_sub(comm->full_comm, grid_dimensions, &(comm->vertical_comm));
} /* Setup_grid */

/*********************************************************/
/*
La seguente funzione crea il tipo di dato MPI per la matrice di tipo LOCAL_MATRIX_T, cosi che possa essere inviata tra i processi 
*/
void Build_matrix_type(LOCAL_MATRIX_T*  local_A) {
    MPI_Datatype  temp_mpi_t;
    int           block_lengths[3];
    MPI_Aint      displacements[3];
    MPI_Datatype  typelist[3];
    MPI_Aint      start_address;
    MPI_Aint      address;

    // Crea un MPI_DATATYPE (temp_mpi_t) che corrisponde a un blocco contiguo di MY_MPI_TYPE
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
/*
La seguente funzione fa leggere la matrice A da 0,0,0 che la distribuisce le colonne (i vari A(i)) ai processi 0,0,i (i=0,..,h-1), che poi le divideranno ulteriormente sul proprio piano per procedere con algoritmo di cannon per A(i)*B(i)
Per ogni riga
    per per ogni altezza
        calcola destinatario di quella riga 
        se sono io
            memorizza nella mia matrice
        se non sono io
            invia al destinatario
*/
void Read_matrix_by_columns(LOCAL_MATRIX_T* full_A , GRID_INFO_T* comm , char* matrix_filename, int* n, int* m) {
    int        mat_row, mat_col;
    int        grid_row, grid_col;
    int        coords[3] = {0,0,-1}; // i processori 0,0,h di ogni griglia ricevono la matrice Ai, il valore -1 è un vaore sentinella che indicherà l'altezza
    int        dest;
    MPI_Status status;
    FILE *finptr;

    if (comm->full_rank == 0) {
        finptr = fopen(matrix_filename, "r");
        if (finptr == NULL) {
            perror("Error opening file");
            printf("Matrix file not found");
            MPI_Abort(MPI_COMM_WORLD, 1);
            return;
        }        
        fscanf(finptr, "%d %d", n, m); // Leggo n m
    }

    // Invio in broadcast le dimensioni della matrice
    MPI_Bcast(n, 1, MPI_INT, 0, MPI_COMM_WORLD); 
    MPI_Bcast(m, 1, MPI_INT, 0, MPI_COMM_WORLD);

    full_A->n_rows = *n;
    full_A->n_cols = *m /comm->h;

    if (comm->full_rank == 0) { // Il processo 0,0,0 si occupa della distribuzione
        dtype* temp = (dtype*) malloc((full_A->n_cols)*sizeof(dtype)); // array temporaneo che si occupa di memorizzare le righe per poi inviarle
        for (mat_row = 0;  mat_row < full_A->n_rows; mat_row++) { // per ogni riga della matrice
            for (int grid_height = 0; grid_height < comm->h; grid_height++) { // per ognuna delle matrice A(i)
                coords[2] = grid_height; // la destinazione di questa riga sarà 0,0, grid_height (ovvero la root di ogni livello di altezza)
                MPI_Cart_rank(comm->full_comm, coords, &dest); // calcolo rank a partire da destinazione
                for (mat_col = 0; mat_col < full_A->n_cols; mat_col++){ // leggo ogni elemento della riga
                    if (dest == 0){
                        fscanf(finptr ,"%d", &Entry(full_A, mat_row, mat_col)); // Siamo nel caso del processo root, quindi mat_row coincide con la riga di full_A
                    }else{
                        fscanf(finptr ,"%d", temp + mat_col); // memorizzo la riga un elemento per volta, si potrebbe fare meglio
                    }
                }
                if (dest != 0){
                    MPI_Send(temp, full_A->n_cols, MY_MPI_TYPE, dest, 0, comm->full_comm);
                }
            }
        }
        free(temp);
        fclose(finptr);

    } else if (comm->my_row == 0 && comm->my_col == 0){ // Se non sono root, ma sono 0,0,i mi preparo a ricevere A(i)
            for (mat_row = 0; mat_row < full_A->n_rows; mat_row++){
                MPI_Recv(&Entry(full_A, mat_row, 0), full_A->n_cols, MY_MPI_TYPE, 0, 0, comm->full_comm, &status);
            }
    }

}  // Read_matrix_by_columns

// La seguente funzione è utile per il debug per controllare se la matrice A è stata correttamente divisa per colonne
// Se è stata divisa correttamente, A verrà anche riunificata e stampata
void Print_matrix_by_column(char* filename, LOCAL_MATRIX_T* local_A, GRID_INFO_T* full_grid , int total_rows) {
    int        mat_row, mat_col;
    int        grid_row, grid_col;
    int        source;
    int        coords[3] = {0,0,-1};
    dtype*     temp;
    MPI_Status status;

    if (full_grid->full_rank == 0) {

        FILE *file = fopen(filename, "w+");
        if (file == NULL) {
            perror("Error opening file");
            return;
        }

        fprintf(file, "%d %d\n", local_A->n_rows, local_A->n_cols* full_grid->h);

        temp = (dtype*) malloc(sizeof(dtype) * local_A->n_cols);
        for (mat_row = 0;  mat_row < total_rows; mat_row++) {
            for (grid_col = 0; grid_col < full_grid->h; grid_col++) {
                coords[2] = grid_col;
                MPI_Cart_rank(full_grid->full_comm, coords, &source);
                if (source == 0) {
                    for(mat_col = 0; mat_col < local_A->n_cols; mat_col++)
                        fprintf(file, "%d ", Entry(local_A, mat_row, mat_col));
                } else {
                    MPI_Recv(temp, local_A->n_cols, MY_MPI_TYPE, source, 0, full_grid->full_comm, &status);
                    for(mat_col = 0; mat_col < local_A->n_cols; mat_col++)
                        fprintf(file, "%d ", temp[mat_col]);
                }
            }
            fprintf(file,"\n");
        }
        free(temp);
        fclose(file);
    } else {
        if (full_grid->my_row == 0 && full_grid->my_col == 0){
            for (mat_row = 0; mat_row < local_A->n_rows; mat_row++) 
                MPI_Send(&Entry(local_A, mat_row, 0), local_A->n_cols, MY_MPI_TYPE, 0, 0, full_grid->full_comm);
            }
    }
                     
}  // Print_matrix_by_column 

/*********************************************************/
/*
La seguente funzione fa leggere la matrice B da 0,0,0 che la distribuisce per righe (i vari B(i)) ai processi 0,0,i (i=0,..,h-1), che poi le divideranno ulteriormente sul proprio piano per procedere con algoritmo di cannon per A(i)*B(i)
Per ogni riga
    per per ogni altezza
        calcola destinatario di quella riga 
        se sono io
            memorizza nella mia matrice
        se non sono io
            invia al destinatario
*/
void Read_matrix_by_rows(LOCAL_MATRIX_T* full_B , GRID_INFO_T* comm , char* matrix_filename, int* n, int* m) {
    int        mat_row, mat_col;
    int        grid_row, grid_col;
    int        coords[3] = {0,0,-1}; // i processori 0,0 di ogni griglia ricevono la matrice Ai, il valore -1 indicherà l'altezza ed ora è un valore sentinella
    int        dest;
    MPI_Status status;
    FILE *finptr;

    if (comm->full_rank == 0) {
        finptr = fopen(matrix_filename, "r");

        if (finptr == NULL) {
            perror("Error opening file");
            printf("Matrix file not found");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }        
        
        fscanf(finptr, "%d %d", n, m);
    }

    // Invio in broadcast le dimensioni della matrice
    MPI_Bcast(n, 1, MPI_INT, 0, MPI_COMM_WORLD); 
    MPI_Bcast(m, 1, MPI_INT, 0, MPI_COMM_WORLD); 

    // Supposing the number of processes for each submatrix is perfect square 2^q with q even
    full_B->n_rows = *n / comm->h;
    full_B->n_cols = *m;
    
    if (comm->full_rank == 0) { // se sono root
        dtype* temp = (dtype*) malloc((full_B->n_cols)*sizeof(dtype)); // array temporaneo che memorizzerà le righe di B(i) prima di mandarle   
        for (mat_row = 0;  mat_row < *n; mat_row++) { // per ogni riga
            coords[2] = mat_row/full_B->n_rows; // calcolo destinazione della riga 
            MPI_Cart_rank(comm->full_comm, coords, &dest); // calcolo rank della destinazione della riga
            for (mat_col = 0; mat_col < full_B->n_cols; mat_col++){ 
                if (dest == 0){
                    fscanf(finptr ,"%d", &Entry(full_B, mat_row, mat_col)); // Siamo nel caso del processo root, quindi mat_row coincide con la riga di full_A
                }else{
                    fscanf(finptr ,"%d", temp + mat_col);
                }
            }
            if (dest != 0){
                MPI_Send(temp, full_B->n_cols, MY_MPI_TYPE, dest, 0, comm->full_comm);
            }
        }
        free(temp);
        fclose(finptr);

    } else {
        if (comm->my_row == 0 && comm->my_col == 0){
            for (mat_row = 0; mat_row < full_B->n_rows; mat_row++){
                MPI_Recv(&Entry(full_B, mat_row, 0), full_B->n_cols, MY_MPI_TYPE, 0, 0, comm->full_comm, &status);
            }
       }
    }

}  // Read_matrix_by_rows

// La seguente funzione è utile per il debug per controllare se la matrice B è stata correttamente divisa per righe
// Se è stata divisa correttamente, B verrà anche riunificata e stampata
void Print_matrix_by_rows(char* filename, LOCAL_MATRIX_T* local_B, GRID_INFO_T* full_grid , int total_rows) {
    int        mat_row, mat_col;
    int        grid_row, grid_col;
    int        source;
    int        coords[3] = {0,0,-1};
    dtype*     temp;
    MPI_Status status;

    if (full_grid->full_rank == 0) {
        FILE *file = fopen(filename, "w+");
        if (file == NULL) {
            perror("Error opening file");
            return;
        }
        
        fprintf(file, "%d %d\n", local_B->n_rows *full_grid->h , local_B->n_cols);
        
        temp = (dtype*) malloc(sizeof(dtype) * local_B->n_cols);
        for (mat_row = 0;  mat_row < local_B -> n_rows * full_grid->h; mat_row++) {
            for (mat_col = 0; mat_col < local_B->n_cols; mat_col++) {
                coords[2] = mat_row / local_B -> n_rows;
                MPI_Cart_rank(full_grid->full_comm, coords, &source);
                if (source == 0) {
                    for(mat_col = 0; mat_col < local_B->n_cols; mat_col++)
                        fprintf(file, "%d ", Entry(local_B, mat_row, mat_col));
                } else {
                    MPI_Recv(temp, local_B->n_cols, MY_MPI_TYPE, source, 0, full_grid->full_comm, &status);
                    for(mat_col = 0; mat_col < local_B->n_cols; mat_col++)
                        fprintf(file, "%d ", temp[mat_col]);
                }
            }
            fprintf(file,"\n");
        }
        free(temp);
        fclose(file);
    } else {
        if (full_grid->my_row == 0 && full_grid->my_col == 0){
            for (mat_row = 0; mat_row < local_B->n_rows; mat_row++) 
                MPI_Send(&Entry(local_B, mat_row, 0), local_B->n_cols, MY_MPI_TYPE, 0, 0, full_grid->full_comm);
            }
    }
                     
}  // Print_matrix_by_rows 

/*********************************************************/
/*
La seguente funzione viene lanciata da ogni processo 0,0,i (i=0,..,h-1) sulla propria griglia e distribuisce le sottomatrici di A(i) o B(i) in modo da poter calcolare C(i) tramite Cannon
*/
void Distribute_local_matrices(LOCAL_MATRIX_T* full_A, LOCAL_MATRIX_T* local_A , GRID_INFO_T* comm) {
    int        mat_row, mat_col;
    int        grid_row, grid_col;
    int        coords[3];
    int        dest;
    MPI_Status status;

    local_A->n_rows = full_A->n_rows / comm->q;
    local_A->n_cols = full_A->n_cols / comm->q;

    if (comm->my_col == 0 && comm-> my_row == 0) { // se sono il "root" del mio livello comm->my_z, questo codice lo eseguirà ogni 0,0,i
        for (mat_row = 0;  mat_row < full_A->n_rows; mat_row++) { // per ogni riga
            coords[0] = mat_row/local_A->n_rows;
            for (grid_col = 0; grid_col < comm->q; grid_col++) { // per ogni processore nella griglia
                coords[1] = grid_col;
                coords[2] = comm->my_z; // ogni "local root" invierà dati solo coloro che sono sul suo stesso piano
                if (comm->my_row == coords[0] && comm->my_col == coords[1] && comm->my_z == coords[2]){ // "se è per me stesso",  vale "per ogni root"
                    for (mat_col = 0; mat_col < local_A->n_cols; mat_col++){
                        Entry(local_A, mat_row, mat_col) = Entry(full_A, mat_row, mat_col ); // Siamo nel caso del processo root, quindi mat_row coincide con la riga di local_A e anche matcol
                    }
                } else {
                    MPI_Cart_rank(comm->full_comm, coords, &dest); // calcolo rank di destinazione
                    MPI_Send(&Entry(full_A, mat_row , grid_col*local_A->n_cols), local_A->n_cols, MY_MPI_TYPE, dest, 0, comm->full_comm); // invio riga al processo interessato
                }
            }
        }
    } else { // se non sono "una root"
        int source;
        int source_coords[3] = {0,0,-1}; // rappresenterà le coordinate del root del proprio piano che gli invierà i dati
        source_coords[2] = comm->my_z;
        MPI_Cart_rank(comm->full_comm, source_coords, &source);
        for (int local_mat_row = 0; local_mat_row < local_A->n_rows; local_mat_row++){
           MPI_Recv(&Entry(local_A, local_mat_row, 0), local_A->n_cols, MY_MPI_TYPE, source, 0, comm->full_comm, &status);
       }
    }

}  // Distribute_local_matrices

// La seguente funzione è usata per stampare la matrice composta da tutte le sottomatrici di un determinato livello
// La seguente funzione stamperà quindi per esempio, la matrice C(i) (che è composta da tante sottomatrici quadrate quanti sono i processori del livello)
// Questa funzione si può usare per stampare per esempio i vari C(i), che si compongono di tutti i vari local_C di ogni processore dello stesso piano... 
// Oppure in particolare è usata per stampare il risultato finale (perchè anche esso si compone dei vari Cln che sono spalmati sui procesori di uno stesso piano)
void Print_level_matrix(char* filename, LOCAL_MATRIX_T* local_A, GRID_INFO_T* grid) {
    int        mat_row, mat_col;
    int        grid_row, grid_col;
    int        source;
    int        coords[2];
    dtype*     temp;
    MPI_Status status;

    /*
    // Le seguenti due righe di codice sono utili se si vogliono stampare i risultati parziali C(i), in modo che i file che li contengono abbiano nomi diversi in base al "livello" (z) a cui corrispondono
    if (strlen(filename) < 2 || strlen(filename) > 64) {
        fprintf(stderr, "Filename of Ai/Bi non adeguato\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
        return;
    }
    char new_filename[12]; // Adjust size as necessary
    strcpy(new_filename, filename);
    char new_char = '0' + (grid->my_z % 10); // Ensure it's a single digit
    new_filename[1] = new_char;
    */

    int my_grid_rank;
    int my_coords[2] = {grid->my_row, grid->my_col};
    MPI_Cart_rank(grid->grid_comm, my_coords, &my_grid_rank);

    if (my_grid_rank == 0) {
        FILE *file = fopen(filename, "w+");
        if (file == NULL) {
            perror("Error opening file");
            return;
        }

        fprintf(file, "%d %d\n", local_A->n_rows * grid->q, local_A->n_cols* grid->q);

        temp = (dtype*) malloc(local_A->n_rows*local_A->n_cols);
        for (mat_row = 0;  mat_row < local_A->n_rows * grid->q; mat_row++) {
            grid_row = mat_row/local_A->n_rows;
            coords[0] = grid_row;
            for (grid_col = 0; grid_col < grid->q; grid_col++) {
                coords[1] = grid_col;
                MPI_Cart_rank(grid->grid_comm, coords, &source);
                if (source == 0) {
                    for(mat_col = 0; mat_col < local_A->n_cols; mat_col++)
                        fprintf(file, "%d ", Entry(local_A, mat_row, mat_col));
                } else {
                    MPI_Recv(temp, local_A->n_cols, MY_MPI_TYPE, source, 0, grid->grid_comm, &status);
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
            MPI_Send(&Entry(local_A, mat_row, 0), local_A->n_cols, MY_MPI_TYPE, 0, 0, grid->grid_comm);
    }
}  // Print_level_matrix 

/*********************************************************/
// La seguente funzione setta tutti i valori della matrice local_A a 0
void Set_to_zero(LOCAL_MATRIX_T*  local_A  /* out */) {
    int i, j;
    for (i = 0; i < local_A->n_rows; i++)
        for (j = 0; j < local_A->n_cols; j++)
            Entry(local_A,i,j) = 0;
}  /* Set_to_zero */

/*********************************************************/
/*
Funzione che effettua moltiplicazione sequenziale standard di due matrici local_A e local_B
*/
void Local_matrix_multiply(LOCAL_MATRIX_T*  local_A , LOCAL_MATRIX_T*  local_B, LOCAL_MATRIX_T*  local_C) {
    int i, j, k;
    for (i = 0; i < local_C->n_rows; i++)
        for (j = 0; j < local_C->n_cols; j++)
            for (k = 0; k < local_A->n_cols; k++)
                Entry(local_C,i,j) = Entry(local_C,i,j) + Entry(local_A,i,k)*Entry(local_B,k,j);

}  


/*********************************************************/
/*
La seguente funzione ha solo scopo di debug per stampare una matrice "locale" di un processore
// NOTA: se si chiama questa funzione su tutti i processori l'output sarà un disastro
*/
void Print_local_matrix(char* title /* in */, LOCAL_MATRIX_T* local_A  /* in */, GRID_INFO_T* grid /* in */) {
    printf("Process %d > comm_row = %d, comm_col = %d, comm_height = %d\n", grid->full_rank, grid->my_row, grid->my_col, grid->my_z);
    for (int i = 0; i <local_A->n_rows ; i++) {
        for (int j = 0; j < local_A->n_cols; j++){
            printf("%d ", Entry(local_A,i,j));
        }
        printf("\n");
    }
}  /* Print_local_matrix */

/******************************************************************************************************************/

/*
La seguente funzione effettua il setup per l'algoritmo di Cannon
    - shift verso sinistra di i colonne gli elementi di A, dove i è la riga dell'elemento
    - shift verso l'alto di j righe gli elementi di B, dove j è la colonna dell'elemento
*/
void Cannon_Setup(GRID_INFO_T* grid, LOCAL_MATRIX_T* local_A, LOCAL_MATRIX_T* local_B){
    MPI_Status status;
    int source_rank, destination_rank;
    int source_position[2];
    int destination_position[2]; 

    // A
    MPI_Cart_shift(grid->grid_comm, 1, -grid->my_row, &source_rank, &destination_rank); // sulla base della mia posizione trovo il rank che è nella posizione in cui devo andare
    MPI_Cart_coords(grid->grid_comm, source_rank, 2, source_position);
    MPI_Cart_coords(grid->grid_comm, destination_rank, 2, destination_position);
    //printf("Io sono in A in %d %d e mi sposto di %d in %d %d mentre al posto mio viene %d %d\n", grid->my_row, grid->my_col, -grid->my_row, destination_position[0], destination_position[1], source_position[0], source_position[1]);
    MPI_Sendrecv_replace(local_A, 1, local_matrix_mpi_t, destination_rank, 0, source_rank, 0, grid->grid_comm, &status);

    // B
    MPI_Cart_shift(grid->grid_comm, 0, -grid->my_col, &source_rank, &destination_rank); // sulla base della mia posizione trovo il rank che è nella posizione in cui devo andare
    MPI_Cart_coords(grid->grid_comm, source_rank, 2, source_position);
    MPI_Cart_coords(grid->grid_comm, destination_rank, 2, destination_position);
    //printf("Io sono in B in %d %d e mi sposto di %d in %d %d mentre al posto mio viene %d %d\n", grid->my_row, grid->my_col, -grid->my_row, destination_position[0], destination_position[1], source_position[0], source_position[1]);
    MPI_Sendrecv_replace(local_B, 1, local_matrix_mpi_t, destination_rank, 0, source_rank, 0, grid->grid_comm, &status);
}

/*
La seguente funzione implementa l'algoritmo di cannon, in cui dopo il setup, ogni processore moltiplica le matrici lolcal_A e local_B che possiede, e invia local_A a sinistra e local_B in alto... per poi ricevere da destra e da giù le nuove matrici local_A e local_B
*/
void Cannon_Multiplication(GRID_INFO_T* comm, LOCAL_MATRIX_T* local_A, LOCAL_MATRIX_T* local_B, LOCAL_MATRIX_T* local_C){
    MPI_Status status;
    int source_rank, destination_rank;

    for (int i=0; i < comm->q - 1; i++){
        Local_matrix_multiply(local_A, local_B, local_C);

        MPI_Cart_shift(comm->grid_comm, 1, -1, &source_rank, &destination_rank); // sulla base della mia posizione trovo il rank che è nella posizione in cui devo andare
        MPI_Sendrecv_replace(local_A, 1, local_matrix_mpi_t, destination_rank, 0, source_rank, 0, comm->grid_comm, &status);

        MPI_Cart_shift(comm->grid_comm, 0, -1, &source_rank, &destination_rank); // sulla base della mia posizione trovo il rank che è nella posizione in cui devo andare
        MPI_Sendrecv_replace(local_B, 1, local_matrix_mpi_t, destination_rank, 0, source_rank, 0, comm->grid_comm, &status);
    }

    Local_matrix_multiply(local_A, local_B, local_C);
}
