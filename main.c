/****** File: LSH_L2.c ******/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "lib/utils.h"

#define min(x, y)    ((x) < (y) ? (x) : (y))
#define max(x, y)    ((x) > (y) ? (x) : (y))


//line 432  can't be moved outside

// apply LSH buckets - hash vals = 0; print result at the end

// ai + bi * q/ w?? -bi ?? ai*q + bi

struct LSH_Parameters {
    char m_max, m,
            L_max, L;
    double W_init, W;
    double **hfunction, *b; // Array sizes: hfunction[m_max*L_max][dim], b[dim]
    // datum_hashval = <datum - b, hfunction> / W
};


struct LSH_Buckets {
    int nclusters, max_nclusters, max_clustersize,
            *clustersize, *clustersize_limit;
    char **cluster_hashval;  // cluster_hashval[max_nclusters][L_max*m_max]
    int **data_indices;     // data_indices[max_nclusters][max_clustersize]
};


struct LSH_Performance {
    double ClusteringTime,
            avg_SearchingTime, wrst_SearchingTime,
            avg_PtsChecked, wrst_PtsChecked,
            avg_Dist, wrst_Dist, avg_RelativeDist, wrst_RelativeDist, // RelativeDist = approx_distance/exact_distance
            avg_rho,   // The recall rho = 1 if closest pt in multi-probed buckets, 0 otherwise
            tau,       // The selectivity tau = (num_PtsChecked)/ndata
            tau_other; // = (m * l * num_buckets)/(dim*ndata)
    int num_buckets;
};

int applyLSH(int dim, int i0, int im, double *data, struct LSH_Parameters *param_ptr,  // input
             double *datum, char *datum_hashval,                                        // buffers
             struct LSH_Buckets *buckets_ptr, struct LSH_Performance *performance_ptr);

int searchLSH(int dim, int i0, int im, double *data, int nqueries, double *queries, // input
              struct LSH_Parameters *param_ptr, struct LSH_Buckets *buckets_ptr,     // input
              double *datum, char *datum_hashval,                                     // buffers
              struct LSH_Performance *performance_ptr);

//double gaussian_rand(char phase) /*** phase = 0 or 1 ***/
//{
//    double U, V, Z, r;
//
//    U = 2.0 * rand() / RAND_MAX - 1.0;
//    V = 2.0 * rand() / RAND_MAX - 1.0;
//    r = U * U + V * V;
//    r = sqrt(-2.0 * log(r) / r);
//
//    if (phase == 0) {
//        Z = U * r;
//    } else { Z = V * r; }
//
//    return Z;
//}
//print dataset
double z1;

double gaussian_rand(int phase) {
    const double epsilon = 2.22507e-308;
    const double two_pi = 2 * M_PI;

    if (phase == 1)
        return z1;

    double u1, u2;
    do {
        u1 = (rand() + 1.) / (RAND_MAX + 2.);
        u2 = rand() / (RAND_MAX + 1.);
    } while (u1 <= epsilon);

    double z0;
    z0 = sqrt(-2.0 * log(u1)) * cos(two_pi * u2);
    z1 = sqrt(-2.0 * log(u1)) * sin(two_pi * u2);

    return z0;
}

/**************************************************************************/
/***  Array sizes: centroid[dim],  dimMinMax[2*dim],  dimVariance[dim]  ***/
/***  Array sizes: datum[dim], cluster_center[2][dim], cluster_size[2]  ***/
/**************************************************************************/
int init_LSHparameters(int dim, int i0, int im, double *data,             // input: dataset
                       double *centroid, double *dimMinMax, double *dimVariance,         // intermediate data
                       double *datum, double *cluster_center[2], int cluster_size[2],   // buffers
                       struct LSH_Parameters *param_ptr,
                       struct LSH_Buckets *buckets_ptr)                                // output
{
    srand(1);
    int i, j, k, dim_maxvar;
    char membership, m, m_max, L, L_max; // L is num_tables
    double tmp, maxvar, length;

    cluster_size[0] = 0;
    cluster_size[1] = 0;

/******** Scan dataset to find centroid, variance, and min & max of each dimension ********
 ********                      and find initial values for W, m, L                 ********/
    for (j = 0; j < dim; j++) {
        tmp = data[i0 * dim + j];
        centroid[j] = tmp;
        dimVariance[j] = tmp * tmp;
        dimMinMax[j] = tmp;
        dimMinMax[j + dim] = tmp;  // Min & Max pts. in dimension j
    }

    for (i = i0 + 1; i < im; i++) {// scan dataset to find centroid and min & max of each dimension
        memcpy(datum, data + i * dim, dim * sizeof(double));
        for (j = 0; j < dim; j++) {
            centroid[j] += datum[j];
            if (datum[j] > dimMinMax[j + dim]) dimMinMax[j + dim] = datum[j];
            else if (datum[j] < dimMinMax[j]) dimMinMax[j] = datum[j];
            dimVariance[j] += datum[j] * datum[j];
        }
    }


    for (j = 0; j < dim; j++) {
        centroid[j] /= (im - i0);
        dimVariance[j] -= (im - i0) * centroid[j] * centroid[j];
    }

    maxvar = dimVariance[0];
    dim_maxvar = 0; ///// Find dimension with the max variance
    for (j = 1; j < dim; j++)
        if (dimVariance[j] > maxvar) {
            maxvar = dimVariance[j];
            dim_maxvar = j;
        }


    for (i = i0; i < im; i++) {/*** Partition dataset at centroid in dim_maxvar into 2 subsets ***/
        memcpy(datum, data + dim * i, dim * sizeof(double));
        membership = ((datum[dim_maxvar] < centroid[dim_maxvar]) ? 0 : 1);
        cluster_size[membership]++;
        for (j = 0; j < dim; j++) cluster_center[membership][j] += datum[j];
    }

    for (k = 0; k < 2; k++) {
        if (cluster_size[k] > 0) {
            for (j = 0; j < dim; j++) cluster_center[k][j] /= (double) cluster_size[k];
        } else {
            printf("ERROR in init_LSHparameters(): cluster[%d] is empty\n", k);
            return 0;
        }
    }
/*** End of Partition dataset into two subsets ***/

/****** Enter initial values for some of the LSH parameters ******/
    param_ptr->m_max = min((int) (0.25 * dim + 0.5), 12);
    param_ptr->L_max = 20;
    param_ptr->W_init = 0.0;
    for (j = 0; j < dim; j++) {
        tmp = cluster_center[0][j] - cluster_center[1][j];
        param_ptr->W_init += tmp * tmp;
    }
    param_ptr->W_init = sqrt(param_ptr->W_init);
    param_ptr->b = dimMinMax; // Choose min value of each dim. May also choose centroid.

/****** Allocate memory for hash functions & generate hash functions ******/

    L_max = param_ptr->L_max;
    m_max = param_ptr->m_max;

    param_ptr->hfunction = (double **) calloc((L_max * m_max), sizeof(double *));
    for (i = 0; i < (L_max * m_max); i++)
        param_ptr->hfunction[i] = (double *) calloc(dim, sizeof(double));

    char ll, mm, phase = 0;
    for (ll = 0; ll < L_max; ll++)
        for (mm = 0; mm < m_max; mm++) { // Generate L_max tables, each with m_max hash functions
            length = 0.0;
            for (j = 0; j < dim; j++) {
                tmp = gaussian_rand(phase);
                phase = 1 - phase;
                param_ptr->hfunction[ll * m_max + mm][j] = tmp;
                length += tmp * tmp;
            }
            length = 1.0 / sqrt(length);
            for (j = 0; j < dim; j++) {
                param_ptr->hfunction[ll * m_max + mm][j] *= length;
            }
        }


    int max_nclusters = 128 * (int) sqrt(im - i0);
    int max_clustersize = 1024;

    buckets_ptr->max_nclusters = max_nclusters;
    buckets_ptr->max_clustersize = max_clustersize;


    return 1;
}

/*************** helper functions *********************/

double calculateDistance(int dim, int idx1, int idx2, const double *data) {
    double distance = 0;
    for (int i = 0; i < dim; ++i) {
        distance += (data[dim * idx1 + i] - data[dim * idx2 + i]) * (data[dim * idx1 + i] - data[dim * idx2 + i]);
    }

    return sqrt(distance);
}

int compareHashVals(struct LSH_Parameters *param_ptr, const char *datum_hashval, const char *cluster_hashval) {
    for (int i = 0; i < param_ptr->L; ++i) {
        for (int j = 0; j < param_ptr->m; ++j) {
            if (datum_hashval[i * param_ptr->m_max + j] != cluster_hashval[i * param_ptr->m_max + j]) {
                return 0;
            }
        }
    }

    return 1;
}

double exactClosestDistance(int dim, int i0, int im, int idx, double *data) {
    double minDistance = RAND_MAX;
    int min_idx = -1;
    for (int i = i0; i < im; ++i) {
        double distance = calculateDistance(dim, idx, i, data);
        if (distance < minDistance) {
            minDistance = distance;
            min_idx = i;
        }
//        printf("i: %d - distance: %f \n", i, distance);
    }

//    printf("closest idx: %d \n", min_idx);

    return minDistance;
}

double
distanceToBucket(int dim, struct LSH_Parameters *param_ptr, const double *datum, const char *datum_hashval,
                 const char *bucket_hashval) {
    double distance = 0, w, leftDistance, tmp; // temp = ai * q + bi
    int l, m, n_step, ll, mm, m_max;
    l = param_ptr->L;
    m = param_ptr->m;
    w = param_ptr->W;
    m_max = param_ptr->m_max;

    for (ll = 0; ll < l; ++ll) {
        for (mm = 0; mm < m; ++mm) {
            tmp = 0;
            n_step = bucket_hashval[ll * m_max + mm] - datum_hashval[ll * m_max + mm];

//            printf("%d, %d \n", bucket_hashval[ll * m_max + mm] ,datum_hashval[ll * m_max + mm]);

            if (n_step != 0) {
                distance += (abs(n_step) - 1) * w;

                for (int j = 0; j < dim; ++j) {
                    tmp += (datum[j] - param_ptr->b[j]) * param_ptr->hfunction[ll * m_max + mm][j];

                }

                leftDistance = tmp - (datum_hashval[ll * m_max + mm] * w);

                distance += n_step < 0 ? leftDistance : w - leftDistance;
            }

        }
    }

    return distance;
}

/******************************************************************************/
int choose_LSHparameters(int dim, int i0, int im, double *data,   // input: small dataset
                         double *datum,                           // buffer
                         int nqueries, double *queries,           // input
                         struct LSH_Parameters *param_ptr,
                         struct LSH_Buckets *buckets_ptr)        // input & output
{
    int i, j, W_count, num_Ws; //datum_hashVal[L_max*m_max], uses only [L*m]
    char m, m_max, L, L_max, *datum_hashval;
    double tmp, W_init, W_min, W_max, W;


    m_max = param_ptr->m_max; //7
    L_max = param_ptr->L_max; //20
    W_init = param_ptr->W_init;

    datum_hashval = (char *) calloc((L_max * m_max), sizeof(char)); // buffer

    struct LSH_Performance ***performances; // Array performances[L_max][m_max][num_Ws]

    W_min = 0.8 * W_init;
    W_max = 1.5 * W_init;
    num_Ws = (int) (0.5 + (0.1 * W_init + W_max - W_min) / (0.1 * W_init)); //47 // *W_init = 10

    performances = (struct LSH_Performance ***) malloc((L_max + 1) * sizeof(struct LSH_Performance **));
    for (i = 0; i <= L_max; i++) {
        performances[i] = (struct LSH_Performance **) malloc((m_max + 1) * sizeof(struct LSH_Performance *));
        for (j = 0; j <= m_max; j++)
            performances[i][j] = (struct LSH_Performance *) malloc((num_Ws + 1) * sizeof(struct LSH_Performance));
    }

/****** Loop thru values of L, M, W to find the best performance parameters ******/
    for (L = 3; L <= L_max; L += 3)
        for (m = 2; m <= m_max; m += 1) {
            W_count = 0;
            for (W = W_min; W <= W_max; W += 0.1 * W_init) {
                ////// Generate buckets using parameters m, L, W and produce search performance results
                param_ptr->m = m;
                param_ptr->L = L;
                param_ptr->W = W;
                applyLSH(dim, i0, im, data, param_ptr,                        // input
                         datum, datum_hashval,                                // buffers
                         buckets_ptr, &performances[L][m][W_count]);         // output

                searchLSH(dim, i0, im, data, nqueries, queries, param_ptr, buckets_ptr, // input
                          datum, datum_hashval,                                         // buffers
                          &performances[L][m][W_count]);                               // output

                W_count++;
            }
        }

    double wrst_relative_dist_threshold = 1.3;
    double avg_rho_threshold = .8;

    int _L = -1, _m = -1;
    double _W = -1, min_clustering_time = RAND_MAX, min_search_time = RAND_MAX;

    for (L = 3; L <= L_max; L += 3)
        for (m = 2; m <= m_max; m += 1) {
            W_count = 0;
            for (W = W_min; W <= W_max; W += 0.1 * W_init) {
                if (performances[L][m][W_count].wrst_RelativeDist < wrst_relative_dist_threshold
                    && performances[L][m][W_count].avg_rho > avg_rho_threshold
                    && performances[L][m][W_count].ClusteringTime < min_clustering_time
                    && performances[L][m][W_count].avg_SearchingTime < min_search_time) {
                    min_clustering_time = performances[L][m][W_count].ClusteringTime;
                    min_search_time = performances[L][m][W_count].avg_SearchingTime;
                    _L = L;
                    _m = m;
                    _W = W;
                }

                printf(" L : %d, m: %d, W: %f, performances: ", L, m, W);

                printf("wrst_PtsChecked %f, avg_PtsChecked %f, wrst_Dist %f, avg_Dist %f, wrst_RelativeDist %f, avg_RelativeDist %f, avg_rho %f, tau %f, tau_other %f, num_buckets %d, ClusteringTime %f, avg_SearchingTime %f, wrst_SearchingTime %f \n",
                       performances[L][m][W_count].wrst_PtsChecked,
                       performances[L][m][W_count].avg_PtsChecked,
                       performances[L][m][W_count].wrst_Dist,
                       performances[L][m][W_count].avg_Dist,
                       performances[L][m][W_count].wrst_RelativeDist,
                       performances[L][m][W_count].avg_RelativeDist,
                       performances[L][m][W_count].avg_rho,
                       performances[L][m][W_count].tau,
                       performances[L][m][W_count].tau_other,
                       performances[L][m][W_count].num_buckets,
                       performances[L][m][W_count].ClusteringTime,
                       performances[L][m][W_count].avg_SearchingTime,
                       performances[L][m][W_count].wrst_SearchingTime);
                W_count++;
            }
        }

    printf("final L %d m %d w %f \n", _L, _m, _W);

    getchar();
    // HERE: Choose the m, L, W that produce the best performance
    param_ptr->m = _m;
    param_ptr->L = _L;
    param_ptr->W = _W;

/****** Deallocate memories ******/
    for (i = 0; i < L_max; i++) {
        for (j = 0; j < m_max; j++) free(performances[i][j]);
        free(performances[i]);
    }
    free(performances);
    free(datum_hashval);

    // first deallocate memories inside buckets
    free(buckets_ptr);

    return 1;
}


/**************************************************************************************/
int applyLSH(int dim, int i0, int im, double *data, struct LSH_Parameters *param_ptr,  // input
             double *datum, char *datum_hashval,                                        // buffers
             struct LSH_Buckets *buckets_ptr, struct LSH_Performance *performance_ptr) // output
/*** datum_hashVal[L_max][m_max], uses only [L][m] ***/
{
    int i, j, k, max_nclusters, max_clustersize;
    char membership, m, m_max, L, L_max, ll, mm;
    double tmp, W, time_start, time_end;

    time_start = clock();
    m = param_ptr->m;
    m_max = param_ptr->m_max;
    L = param_ptr->L;
    L_max = param_ptr->L_max;
    W = param_ptr->W;

#if 0   /*** Already defined. Only for reference ***/
    struct LSH_Buckets {
   int   nclusters, max_nclusters, max_clustersize,
         *clustersize, *clustersize_limit; // clustersize_limit[max_nclusters]
                                           // cluster_size[max_nclusters]
   int **cluster_hashval,  // cluster_hashval[max_nclusters][L_max*m_max]
       **data_indices;     // data_indices[max_nclusters][max_clustersize]
} ;
#endif

    max_nclusters = buckets_ptr->max_nclusters;
    max_clustersize = buckets_ptr->max_clustersize;

    buckets_ptr->clustersize_limit = (int *) calloc(max_nclusters, sizeof(int));
    buckets_ptr->clustersize = (int *) calloc(max_nclusters, sizeof(int));

    buckets_ptr->cluster_hashval = (char **) calloc(max_nclusters, sizeof(char *));
    buckets_ptr->data_indices = (int **) calloc(max_nclusters, sizeof(int *));
    for (i = 0; i < max_nclusters; i++) {
        buckets_ptr->cluster_hashval[i] = (char *) calloc((L_max * m_max), sizeof(char));
        buckets_ptr->data_indices[i] = (int *) calloc(max_clustersize, sizeof(int));
    }

/**************** Calculate hash val for each datum and put to a bucket **************/
    int isEqual, new_limit, cluster_size, nclusters;
    buckets_ptr->nclusters = 0;
    for (k = 0; k < (buckets_ptr->nclusters); k++) buckets_ptr->clustersize[k] = 0;
    for (k = 0; k < max_nclusters; k++) buckets_ptr->clustersize_limit[k] = max_clustersize;

    for (i = i0; i < im; i++) {/*** Calc hash val of each datum and put to a bucket ***/
        memcpy(datum, data + i * dim, dim * sizeof(double));

//        printf("datum hashval: \n");

        for (ll = 0; ll < L; ll++)
            for (mm = 0; mm < m; mm++) {
                tmp = 0.0;
                for (j = 0; j < dim; j++)
                    tmp += (datum[j] - param_ptr->b[j]) * (param_ptr->hfunction[ll * m_max + mm][j]);
                datum_hashval[ll * m_max + mm] = (int) floor(tmp / W);
//                printf("%d \n",datum_hashval[ll * m_max + mm]);
            }

//        getchar();

        for (k = 0; k < (buckets_ptr->nclusters); k++) {// Compare datum_hashval with cluster hashvals
            isEqual = compareHashVals(param_ptr, datum_hashval, buckets_ptr->cluster_hashval[k]);
            if (isEqual) {
                buckets_ptr->clustersize[k]++;
                cluster_size = buckets_ptr->clustersize[k];
                if (cluster_size >= buckets_ptr->clustersize_limit[k]) {// Grow ptr->data_indices
                    buckets_ptr->clustersize_limit[k] += 1024; /////////////// Increase by 1024
                    new_limit = buckets_ptr->clustersize_limit[k];
                    buckets_ptr->data_indices[k] =
                            (int *) realloc(buckets_ptr->data_indices[k], new_limit * sizeof(int));
                    if (max_clustersize < buckets_ptr->clustersize_limit[k])
                        max_clustersize = buckets_ptr->clustersize_limit[k];
                }
                buckets_ptr->data_indices[k][cluster_size - 1] = i;

                break;
            }
        }
        if (k == (buckets_ptr->nclusters)) {// datum not in any existing bucket. Create new bucket.
            if (buckets_ptr->nclusters == max_nclusters) { // max_nclusters too small
                printf("ERROR in applyLSH (): max_nclusters too small.n data: %d n max clusters: %d\n\n", im - i0,
                       max_nclusters);
                break;
//                return 0;
            }
            buckets_ptr->nclusters++;

            cluster_size = buckets_ptr->clustersize[k];
            buckets_ptr->clustersize[k]++;

            buckets_ptr->data_indices[k][cluster_size] = i;

            for (ll = 0; ll < L; ll++)
                for (mm = 0; mm < m; mm++) { // buckets_ptr->cluster_hashval[k] = datum_hashval
                    buckets_ptr->cluster_hashval[k][ll * m_max + mm] = datum_hashval[ll * m_max + mm];
                }
        }
    }

    buckets_ptr->max_nclusters = max_nclusters;
    buckets_ptr->max_clustersize = max_clustersize;

    time_end = clock();

    if (performance_ptr) performance_ptr->ClusteringTime = 0.000001 * (time_end - time_start);

    return 1;
} /****** End of function applyLSH() ******/


/*********************************************************/
int searchLSH(int dim, int i0, int im, double *data, int nqueries, double *queries, // input
              struct LSH_Parameters *param_ptr, struct LSH_Buckets *buckets_ptr,     // input
              double *datum, char *datum_hashval,                                     // buffers
              struct LSH_Performance *performance_ptr)                              // output
{
    int counter = 0, n_probing_buckets = (int) (buckets_ptr->nclusters * 5 / 100.0) + 1;

    n_probing_buckets = max(n_probing_buckets, 2);

//    printf("n cluster %d : cluster sizes: \n", buckets_ptr->nclusters);
//
//    for (int l = 0; l < buckets_ptr->nclusters; ++l) {
//        printf("%d \n", buckets_ptr->clustersize[l]);
//        counter += buckets_ptr->clustersize[l];
//    }
//
//    printf("n data : %d", counter);
//    getchar();

    int i, j, k, ll, mm, m, m_max, L, result_idx;
    char isEqual;
    double tmp, exact_distance, W, time_start, time_end, min_distance, search_time;

    m = param_ptr->m;
    m_max = param_ptr->m_max;
    L = param_ptr->L;
    W = param_ptr->W;

    int *bucket_indices = (int *) calloc(buckets_ptr->nclusters, sizeof(int));
    double *bucket_distances = (double *) calloc(buckets_ptr->nclusters, sizeof(double));

    for (i = 0; i < i0; ++i) {
        time_start = clock();

        bool notFound = true;
        min_distance = RAND_MAX;
        result_idx = -1;
        counter = 0;

        //calculate hashvals
        memcpy(datum, data + i * dim, dim * sizeof(double));


        for (ll = 0; ll < L; ll++)
            for (mm = 0; mm < m; mm++) {
                tmp = 0.0;
                for (j = 0; j < dim; j++)
                    tmp += (datum[j] - param_ptr->b[j]) * (param_ptr->hfunction[ll * m_max + mm][j]);
                datum_hashval[ll * m_max + mm] = (int) floor(tmp / W);
            }


//        printf("%d n bucket\n", buckets_ptr->nclusters);
        for (k = 0; k < (buckets_ptr->nclusters); k++) {
            bucket_indices[k] = k;
            bucket_distances[k] = distanceToBucket(dim, param_ptr, datum, datum_hashval,
                                                   buckets_ptr->cluster_hashval[k]);

//            printf("idx : %d, bucket distance : %f \n", bucket_indices[k], bucket_distances[k]);
            if (notFound) {
                isEqual = compareHashVals(param_ptr, datum_hashval, buckets_ptr->cluster_hashval[k]);

                if (isEqual) {
                    double distance;
                    notFound = false;
                    for (j = 0; j < buckets_ptr->clustersize[k]; j++) {
                        counter++;
                        distance = calculateDistance(dim, i, buckets_ptr->data_indices[k][j], data);
                        if (distance < min_distance) {
                            min_distance = distance;
                            result_idx = buckets_ptr->data_indices[k][j];
                        }
                    }
                }
            }
        }

        exact_distance = exactClosestDistance(dim, i0, im, i, data);

//        printf("checked: %d, result idx: %d, distance: %f, exact closest: %f \n", counter, result_idx, min_distance,
//               exact_distance);

        int tmp_idx;
        double tmp_dis;
        for (k = 0; k < buckets_ptr->nclusters; ++k) {
            for (int l = k + 1; l < buckets_ptr->nclusters; ++l) {
                if (bucket_distances[k] > bucket_distances[l]) {
                    tmp_dis = bucket_distances[k];
                    bucket_distances[k] = bucket_distances[l];
                    bucket_distances[l] = tmp_dis;

                    tmp_idx = bucket_indices[k];
                    bucket_indices[k] = bucket_indices[l];
                    bucket_indices[l] = tmp_idx;
                }
            }
        }

        for (int n = 1; n < n_probing_buckets; ++n) {
            k = bucket_indices[n];
            double distance;
            for (j = 0; j < buckets_ptr->clustersize[k]; j++) {
                counter++;
                distance = calculateDistance(dim, i, buckets_ptr->data_indices[k][j], data);
                if (distance < min_distance) {
                    min_distance = distance;
                    result_idx = buckets_ptr->data_indices[k][j];
                }
            }

//            printf("idx : %d, distance: %f \n", bucket_indices[n], bucket_distances[n]);
        }
//        printf("after probing: checked: %d, idx: %d, distance: %f \n", counter, result_idx, min_distance);

        if (performance_ptr) {
            time_end = clock();

            search_time = 0.000001 * (time_end - time_start);

            performance_ptr->avg_SearchingTime += search_time;
            if (performance_ptr->wrst_SearchingTime < search_time) performance_ptr->wrst_SearchingTime = search_time;

            performance_ptr->avg_PtsChecked += counter;
            if (performance_ptr->wrst_PtsChecked < counter) performance_ptr->wrst_PtsChecked = counter;

            performance_ptr->avg_Dist += min_distance - exact_distance;
            if (performance_ptr->wrst_Dist < min_distance - exact_distance)
                performance_ptr->wrst_Dist = min_distance - exact_distance;

            performance_ptr->avg_RelativeDist += min_distance / exact_distance;
            if (performance_ptr->wrst_RelativeDist < min_distance / exact_distance)
                performance_ptr->wrst_RelativeDist = min_distance / exact_distance;

            performance_ptr->avg_rho += min_distance == exact_distance ? 1 : 0;

            performance_ptr->tau += counter / (double) (im - i0);
        }
    }

    if (performance_ptr) {
        performance_ptr->num_buckets = buckets_ptr->nclusters;
        performance_ptr->avg_SearchingTime /= i0;
        performance_ptr->avg_PtsChecked /= i0;
        performance_ptr->avg_Dist /= i0;
        performance_ptr->avg_RelativeDist /= i0;
        performance_ptr->avg_rho /= i0;

//        performance_ptr->tau /= i0;

        performance_ptr->tau_other = (m * L * buckets_ptr->nclusters) / (double) (dim * (im - i0));
    }

    return 1;
}


#define DATASET        1

int main() {
    int dim, ndata, i0, im, nqueries, cluster_size[2], n_batches, data_size;

    double **data, *centroid, *dimMinMax, *dimVariance, *queries, *datum;

    double *cluster_center[2];
    char *datum_hashval;

    struct LSH_Parameters *param_ptr;
    struct LSH_Buckets *buckets_ptr;

    param_ptr = (struct LSH_Parameters *) malloc(sizeof(struct LSH_Parameters));
    buckets_ptr = (struct LSH_Buckets *) malloc(sizeof(struct LSH_Buckets));

#if (DATASET == 1)/*** Read data from HIGGS binary file of double floating-pt data ***/
    FILE *fp = fopen("../data_sets/tr_HIGGS.dat", "rb");


    dim = 29;

//    data_size = 9000000;
//    ndata = 1000000;

    data_size = 90000;

    ndata = 10000;
    n_batches = data_size / ndata;

    nqueries = ndata / 10;

    data = (double **) malloc(n_batches * sizeof(double *));
    for (int i = 0; i < n_batches; ++i) {
        data[i] = (double *) calloc(dim * ndata, sizeof(double));
    }

    queries = data[0]; // 1st nqueries data in data[] are queries, i.e. i0=nqueries, im=ndata
    if (fp != NULL) {
        for (int i = 0; i < n_batches; ++i) {
            fseek(fp, sizeof(double) * dim * ndata * i, SEEK_SET);
            fread(data[i], sizeof(double), dim * ndata, fp);
        }
        fclose(fp);
    } else {
        printf("failed to open file");
    }
#endif

#if (DATASET == 2)/*** bio_train dataset from http://osmot.cs.cornell.edu/kddcup/  ***/
    FILE *fp = fopen("../data_sets/bio_train.dat","rb");
     dim = 77 ;   ndata = 1000 ;
//     dim = 74 ;   ndata = 145751 ;
     data = (double *) calloc(dim * ndata, sizeof(double)) ;
     if(fp!= NULL) {
         fread(data, sizeof(double), dim * ndata, fp);
         fclose(fp);
     } else {
         printf("failed to open file");
     }
#endif

#if (DATASET == 3)
    FILE *fp = fopen("../data_sets/tlc_nyc2016_norm_41M_dim16.dat", "rb");

    dim = 16;
//    ndata = 41000000;
    ndata = 100000;
    data = (double *) calloc(dim * ndata, sizeof(double));
    if (fp != NULL) {
        fread(data, sizeof(double), dim * ndata, fp);
        fclose(fp);
    } else {
        printf("failed to open file");
    }
#endif

#if (DATASET == 4)
    FILE *fp = fopen("../data_sets/heterogeneity_activity_norm.dat", "rb");

    dim = 24;
//    ndata = 16000000;
    ndata = 10000;
    data = (double *) calloc(dim * ndata, sizeof(double));
    if (fp != NULL) {
        fread(data, sizeof(double), dim * ndata, fp);
        fclose(fp);
    } else {
        printf("failed to open file");
    }
#endif
    // add more data sets here

    centroid = (double *) calloc(dim, sizeof(double));
    dimMinMax = (double *) calloc((2 * dim), sizeof(double));
    dimVariance = (double *) calloc(dim, sizeof(double));
    datum = (double *) calloc(dim, sizeof(double));
    cluster_center[0] = (double *) calloc(dim, sizeof(double));
    cluster_center[1] = (double *) calloc(dim, sizeof(double));


    /*** Determine LSH parameters from a small dataset ***/
//    int n_smalldata = ndata / 100, n_smallqueries = n_smalldata / 10;
//    i0 = n_smallqueries;
//    im = n_smalldata;

    im = ndata;
    i0 = ndata/10;


    init_LSHparameters(dim, i0, im, data[0], centroid, dimMinMax, dimVariance,
                       datum, cluster_center, cluster_size, param_ptr, buckets_ptr);


#if 1
    choose_LSHparameters(dim, i0, im, data[0], datum, nqueries, queries, param_ptr, buckets_ptr);
    /*** End of determining LSH parameters ***/


    /*** Apply LSH with determined parameters to whole dataset ***/
    char L_max, m_max;

    nqueries = ndata / 10;
    i0 = nqueries;
    im = ndata;
    L_max = param_ptr->L_max;
    m_max = param_ptr->m_max;
    datum_hashval = (char *) calloc((L_max * m_max), sizeof(char *));

    applyLSH(dim, i0, im, data[0], param_ptr, datum, datum_hashval,
             buckets_ptr, NULL);
    searchLSH(dim, i0, im, data[0], nqueries, queries, param_ptr, buckets_ptr,
              datum, datum_hashval, NULL);

#endif


/****** Deallocate memory space ******/
//    free(datum_hashval);
    free(datum);
    free(data);

    // first deallocate hfucntion[] inside param_ptr, then deallocate param_ptr
    // first deallocate memories inside buckets_ptr, then deallocate buckets_ptr


} /************ End of main() ************/


