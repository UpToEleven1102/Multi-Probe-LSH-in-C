/****** File: LSH_L2.c ******/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "LSH_L2.h"

#define min(x,y)  	((x) < (y) ? (x) : (y))
#define max(x,y)  	((x) > (y) ? (x) : (y))


/**
 * guassian_rand
 */


double gaussian_rand(char phase) /*** phase = 0 or 1 ***/
{
    double U, V, Z , r;

    U = 2.0*rand()/RAND_MAX - 1.0 ;    V = 2.0*rand()/RAND_MAX - 1.0 ;
    r = U*U + V*V ;
    r = sqrt(-2.0 * log(r) / r) ;

    if(phase == 0) {
//        U = 2.0*rand()/RAND_MAX - 1.0 ;    V = 2.0*rand()/RAND_MAX - 1.0 ;
//        r = U*U + V*V ;
//        r = sqrt(-2.0 * log(r) / r) ;
        Z = U*r;
    }
    else {  Z = V*r ;  }

    return Z;
}



/**************************************************************************/
int init_LSHparameters(int dim, int i0, int im, double *data,              // input: dataset
                       double *centroid, double *dimMinMax, double *dimVariance,         // intermediate data
                       double *datum, double *cluster_center[2], int cluster_size[2],   // buffers
                       struct LSH_Parameters *param_ptr)                                // output
/*** Array sizes: datum[dim], cluster_center[2][dim], cluster_size[2] ***/
{
    srand(1);

    int     i, j, k, dim_maxvar ;
    char    membership, m, m_max, L, L_max ; // L is num_tables
    double  tmp, maxvar ;

/******** Scan dataset to find centroid, variance, and min & max of each dimension ********
 ********                      and find initial values for W, m, L                 ********/
    for(j=0; j<dim; j++) {
        tmp = data[i0*dim+j] ;
        centroid[j] = tmp ;     dimVariance[j] = tmp*tmp ;
        dimMinMax[j] = tmp ;    dimMinMax[j+dim]= tmp ;  // Min & Max pts. in dimension j
    }
    for(i=i0+1; i<im; i++) {// scan dataset to find centroid and min & max of each dimension
        memcpy(datum, data+i*dim, dim*sizeof(double)) ;
        for(j=0; j<dim; j++) {
            centroid[j] += datum[j] ;
            if(datum[j] > dimMinMax[j+dim]) dimMinMax[j+dim] = datum[j] ;
            else { if(datum[j] < dimMinMax[j]) dimMinMax[j] = datum[j] ; }
            dimVariance[j] += datum[j]*datum[j] ;
        }
    }
    for(j=0; j<dim; j++) {
        centroid[j] /= (im-i0) ;
        dimVariance[j] -= (im-i0)*centroid[j]*centroid[j] ;
    }

    maxvar=dimVariance[0];  dim_maxvar=0; ///// Find dimension with the max variance
    for(j=1; j<dim; j++)
        if(dimVariance[j] > maxvar) { maxvar=dimVariance[j];   dim_maxvar=j; }
    for(i=i0; i<im; i++) { /*** Partition whole dataset at centroid along dim_maxvar into 2 subsets ***/
        memcpy( datum, data+dim*i, dim*sizeof(double) ) ;
        membership = ( (datum[dim_maxvar] < centroid[dim_maxvar])? 0 : 1 ) ;
        cluster_size[membership]++ ;
        for(j=0; j<dim; j++)
            cluster_center[membership][j] +=  datum[j] ;
    }
    for(k=0; k<2; k++) {
        if( cluster_size[k] > 0 ) {
            for(j=0; j<dim; j++) cluster_center[k][j] /= (double)cluster_size[k] ;
        }
        else {printf("ERROR in init_LSHparameters(): cluster[%d] is empty\n", k);   return 0; }
    }
/*** End of Partition dataset into two subsets ***/


/****** Enter initial values for some of the LSH parameters ******/
    param_ptr->m_max = min( (int) (0.25*dim+0.5) , 12) ;
    param_ptr->L_max = 20 ;
    param_ptr->W_init = 0.0 ;
    for(j=0; j<dim; j++) {
        tmp = cluster_center[0][j]-cluster_center[1][j] ;
        param_ptr->W_init += tmp*tmp ;
    }
    param_ptr->W_init = sqrt(param_ptr->W_init) ;
    param_ptr->b = dimMinMax ; // Choose min value of each dim. May also choose centroid.

/****** Allocate memory for hash functions & generate hash functions ******/
    param_ptr->hfunction = (double **) calloc( (m_max*L_max), sizeof(double *) ) ;
    for(i=0; i<(m_max*L_max); i++)
        param_ptr->hfunction[i] = (double *) calloc(dim, sizeof(double)) ;

    char phase = 0 ;
    for(l=0; l<=L_max; l++)
        for(i=0; i<=m_max; i++) { // Generate L_max tables, each with m_max hash functions
            length = 0.0 ;
            for(j=0; j<dim; j++) {
                tmp = gaussian_rand(phase) ;    phase = 1 - phase ;
                param_ptr->hfunction[l*m_max + i][j] = tmp ;
                length += tmp*tmp ;
            }
            length = 1.0/sqrt(length) ;
            for(j=0; j<dim; j++) param_ptr->hfunction[l*m_max + i][j] *= length ;
        }


    return 1 ;
}



/******************************************************************************/
int choose_LSHparameters(int dim, int i0, int im, double *data,   // input: small dataset
                         int nqueries, double *queries,           // input
                         struct LSH_Parameters *param_ptr)        // input & output
{
    int     i, j, k,  *datum_hashval;  // datum_hashVal[L_max*m_max], uses only [L*m]
    char    membership, m, m_max, L, L_max, mm, ll ;
    double  tmp, length, W_init, W ;

    m_max = param_ptr->m_max ;
    L_max = param_ptr->L_max ;
    W_init= param_ptr->W_init ;
    datum_hashval = (double *)calloc( (L_max*m_max), sizeof(double*) ) ; // buffer

/****** Loop thru values of L, M, W to find the best performance parameters ******/
    int W_count, num_Ws, W_min, W_max ;
    struct LSH_Performance  ***performances ; // Array performances[L_max][m_max][num_Ws]
    struct LSH_Buckets *buckets_ptr = (struct LSH_Buckets *) malloc(sizeof(struct LSH_Buckets)) ;

    W_min = 0.7*W_init ;
    W_max = 1.5*W_init ;
    num_Ws = (int) (0.5 + (0.1 + W_max - W_min)/0.1) ;
    performances = (struct LSH_Performance ***) malloc(L_max*sizeof(struct LSH_Performance **));
    for(i=0; i<L_max; i++) {
        performances[i] = (struct LSH_Performance **) malloc(m_max*sizeof(struct LSH_Performance *));
        for(j=0; j<m_max; j++)
            performances[i][j] = (struct LSH_Performance *) malloc(num_Ws*sizeof(struct LSH_Performance));
    }
    for(L=3; L<=L_max; L+=3)
        for(m=2; m<=m_max; m+=2) {
            W_count = 0 ;
            for(W=W_min; W<=W_max; W+=0.1*W_init) {
                ////// Generate buckets using parameters m, L, W and produce search performance results
                param_ptr->m = m ;
                param_ptr->L = L ;
                param_ptr->W = W ;

                //pass performance to applyLSH
                applyLSH(dim, i0, im, data, param_ptr,  // input
                         datum, datum_hashval,          // buffers
                         buckets_ptr) ;                 // output
                searchLSH(dim, i0, im, data, nqueries, queries, param_ptr, buckets_ptr, // input
                          datum, datum_hashval,                                         // buffers
                          &performances[L][m][W_count]) ;                               // output
                W_count ++ ;
            }
        }

    // Choose the m, L, W that produce the best performance
    // NEED TO DEFINE HOW TO EVALUATE PERFORMANCE
//    param_ptr-> m = ;
//    param_ptr-> L = ;
//    param_ptr-> W = ;

/****** Deallocate memories ******/
    for(i=0; i<L_max; i++) {
        for(j=0; j<m_max; j++) free(performances[i][j]) ;
        free(performances[i]) ;
    }
    free(performances) ;   free(datum_hashval);

    return 1 ;
}



/*****************************************************************************/
int applyLSH(int dim, int i0, int im, double *data,    // input: dataset
             struct LSH_Parameters *param_ptr,         // input
             double *datum, int *datum_hashval,        // buffers
             struct LSH_Buckets *buckets_ptr)          // output
/*** datum_hashVal[L_max][m_max], uses only [L][m] ***/
{
    int     i, j, k, max_nclusters, max_clustersize ;
    char    membership, m, m_max, L, L_max, ll, mm ;
    double  tmp, W ;

    m = param_ptr->m ;    m_max = param_ptr->m_max ;
    L = param_ptr->L ;    L_max = param_ptr->L_max ;
    W = param_ptr->W ;

#if 0   ////// Already defined. Only for reference
    struct LSH_Buckets {
   int   nclusters, max_nclusters, max_clustersize,
         *clustersize, *clustersize_limit; // clustersize_limit[max_nclusters]
                                           // cluster_size[max_nclusters]
   int **cluster_hashval,  // cluster_hashval[max_nclusters][L_max*m_max]
       **data_indices;     // data_indices[max_nclusters][max_clustersize]
} ;
#endif

    max_nclusters   = 8 * (int)sqrt(im-i0) ;
    max_clustersize = 1024 ;

    buckets_ptr->clustersize_limit= (int *) calloc(max_nclusters, sizeof(int)) ;
    buckets_ptr->clustersize     = (int *) calloc(max_nclusters, sizeof(int)) ;

    buckets_ptr->cluster_hashval = (int **) calloc(max_nclusters, sizeof(int *)) ;
    buckets_ptr->data_indices    = (int **) calloc(max_nclusters, sizeof(int *)) ;
    for(i=0; i<max_nclusters; i++) {
        buckets_ptr->cluster_hashval[i]= (int *) calloc((L_max*m_max), sizeof(int)) ;
        buckets_ptr->data_indices[i]   = (int *) calloc(max_clustersize, sizeof(int)) ;
    }

/**************** Calc hash val for each datum and put to a bucket **************/
    int isEqual, cluster_size, nclusters ;
    buckets_ptr->nclusters = 0 ;
    for(k=0; k<(buckets_ptr->nclusters); k++) buckets_ptr->cluster_size[k] = 0 ;
    fo(k=0; k<max_nclusters; k++) buckets_ptr->clusterzie_limit[k] = max_clustersize ;

    for(i=0; i<(im-i0); i++) {/*** Calc hash val of each datum and put to a bucket ***/
        memcpy(datum, data+i*dim, dim*sizeof(double)) ;

        for(ll=0; ll<L; ll++)
            for(mm=0; mm<m; mm++) {
                tmp = 0.0 ;
                for(j=0; j<dim; j++)
                    tmp += (datum[j] - param_ptr->b[j])*(param_ptr->hfunctions[ll*m_max+mm][j]) ;
                datum_hashval[ll*m_max + mm] = (int) (tmp/W) ;
            }
        for(k=0; k<(buckets_ptr->nclusters); k++) {// Compare datum_hashval with cluster hashvals
            isEqual = compareHashVals(param_ptr, datum_hashval, buckets_ptr->cluster_hashval[k]) ;
            if(isEqual) {
                buckets_ptr->cluster_size[k] ++ ;
                cluster_size = buckets_ptr->cluster_size[k] ;
                if(cluster_size >= buckets_ptr->clusterzie_limit[k]) {// Grow ptr->data_indices
                    buckets_ptr->clusterzie_limit[k] += 1024 ; /////////////// Increase by 1024
                    buckets_ptr->data_indices[k] =
                            (int *) realloc(buckets_ptr->clusterzie_limit[k], sizeof(int)) ;
                    if(max_clustersize < buckets_ptr->clusterzie_limit[k])
                        max_clustersize = buckets_ptr->clusterzie_limit[k] ;
                }
                buckets_ptr->data_indices[k][cluster_size] = i ;
                break ;
            }
        }
        if(k > (buckets_ptr->nclusters)) {// datum not in any existing bucket. Create new bucket.
            buckets_ptr->nclusters ++ ;
            nclusters = buckets_ptr->nclusters ;
            if(nclusters > max_nclusters) { // max_nclusters too small
                printf("ERROR in applyLSH (): max_nclusters too small.\n\n");
                return 0;
            }
            for(ll=0; ll<L; ll++)
                for(mm=0; mm<m; mm++) { // buckets_ptr->cluster_hashval[k] = datum_hashval
                    buckets_ptr->cluster_hashval[k][ll*m_max+mm] = datum_hashval[ll*m_max + mm] ;
                }
        }
    }
    buckets_ptr->max_nclusters   = max_nclusters ;
    buckets_ptr->max_clustersize = max_clustersize ;


    return 1 ;
} /****** End of function applyLSH() ******/


/*********************************************************/
int searchLSH(int dim, int i0, int im, double *data, int nqueries, double *queries, // input
              struct LSH_Paramters *param_ptr, struct LSH_Buckets *buckets_ptr,     // input
              double *datum, int *datum_hashval,                                     // buffers
              struct LSH_Performance *performance_ptr)                              // output
{



    return 1 ;
}


