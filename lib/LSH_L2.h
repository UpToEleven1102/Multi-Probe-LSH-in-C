
#ifndef LSH_PROBING_LSH_L2_H
#define LSH_PROBING_LSH_L2_H

struct LSH_Parameters {
    char    m_max, m,
            L_max, L;
    double  W_init, W ;
    double  **hfunction, *b ; // Array sizes: hfunction[m_max*L_max][dim], b[dim]
    // datum_hashval = <datum - b, hfunction> / W ////////*****
} ;


struct LSH_Buckets {
    int   nclusters, max_nclusters, max_clustersize,
            *clustersize, *clustersize_limit ;
//   struct LSH_Parameters *param_ptr ;
    int **cluster_hashval ; // cluster_hashval[max_nclusters][L_max*m_max]
    int **data_indices;     // data_indices[max_nclusters][max_clustersize]
} ;


struct LSH_Performance {
    double avg_ClusteringTime, wrst_ClusteringTime,
            avg_SearchingTime, wrst_SearchingTime,
            avg_RelativeDist, wrst_RelativeDist, // RelativeDist = approx_distance/exact_distance
            avg_rho , // The recall rho = 1 if closest pt in multi-probed buckets, 0 otherwise
            tau,// The selectivity tau = (num pts checked) / ndata
            ion_metric;  // num buckets /n_data

} ;

int init_LSHparameters(int dim, int i0, int im, double *data,              // input: dataset
                       double *centroid, double *dimMinMax, double *dimVariance,         // intermediate data
                       double *datum, double *cluster_center[2], int cluster_size[2],   // buffers
                       struct LSH_Parameters *param_ptr);                                // output

#endif //LSH_PROBING_LSH_L2_H
