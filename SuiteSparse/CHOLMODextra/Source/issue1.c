/*
 * Issue #1
 *
 * Problem with getting the inverse
 *
 * https://github.com/jluttine/cholmod-extra/issues/1
 */

#include "cholmod_extra.h"
#include <cholmod.h>
#include <cholmod_internal.h>
#include <math.h>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main(void)
{
    int N = 2 ;
    int i, j, n ;
    int nz = 0;
    double *Ax ;
    double x, error ;
    cholmod_dense *A, *invK, *spinvK, *I ;
    cholmod_sparse *K, *V ;
    cholmod_factor *L ;
    cholmod_common Common,*cm ;
    clock_t start, end;
    double cpu_time_used;
    cm=&Common;
    // Start using CHOLMOD
    cholmod_start(cm) ;
    cm->print=5;
    /* SPARSE COVARIANCE MATRIX CONSTRUCTION */

    // Generate random symmetric positive (semi)definite matrix
    A = cholmod_zeros(N, N, CHOLMOD_REAL, &Common) ;
    Ax =(double*) A->x ;
    nz = N ;

    // Make positive-definite by adding something positive to the
    // diagonal
    for (n = 0; n < N; n++)
    {
        Ax[n+n*N] += 5;
    }

    // Make the matrix sparse
    K = cholmod_dense_to_sparse(A, TRUE, &Common) ;
    K->stype = 1 ; // NEED TO MAKE THE MATRIX SYMMETRIC

    // Identity matrix
    I = cholmod_eye(N,N,CHOLMOD_REAL,&Common) ;

    /* SIMPLICIAL */

    // Factorize
    Common.supernodal = CHOLMOD_SIMPLICIAL ;
    L = cholmod_analyze(K, &Common) ;
    cholmod_factorize(K, L, &Common) ;
    invK = cholmod_solve(CHOLMOD_A, L, I, &Common) ;

    // Compute the sparse inverse and the full inverse
    start = clock();
    V = cholmod_spinv(L, &Common) ;
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

    // Show results
    cholmod_print_sparse(K,"Original",&Common) ;
    cholmod_print_factor(L,"Factor",&Common) ;
    cholmod_print_sparse(V,"Sparse inverse",&Common) ;
    cholmod_print_dense(invK,"Dense inverse",&Common);

    // Free memory
    cholmod_free_factor(&L, &Common) ;
    cholmod_free_sparse(&K, &Common) ;
    cholmod_free_dense(&I, &Common) ;
    cholmod_free_dense(&A, &Common) ;
    cholmod_finish(&Common) ;
    return 0 ;
}
