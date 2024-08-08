#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <time.h>
#include <omp.h>


#define ERR_NO_NUM -1
#define ERR_NO_MEM -2
#define FREED_RAND -3


int *clusters_sizes;
/*As funções de print não serão paralelizadas, não há sentido em fazer isso*/
void print_vector(double *vector, int vector_size) {
	printf("(");

	
	for (int i = 0; i < vector_size; ++i) {
		if (i > 0)
		
		
			printf(", ");
		
		printf("%.2f", vector[i]);
	}
	
	printf(")");
}

void print_observations(double **observations, int observations_size, int vector_size) {
	printf("[");

	
	for (int i = 0; i < observations_size; ++i) {
		if (i > 0)
		
			printf(", ");

		print_vector(observations[i], vector_size);
	}
	
	printf("]");
}

void print_clusters(double ***clusters, int k, int observations_size, int vector_size) {
	printf("{");

	
	for (int i = 0; i < k; ++i) {
		if (i > 0)
		
			printf(", ");
		print_observations(clusters[i], clusters_sizes[i], vector_size);
	}
	
	free(clusters_sizes);
	printf("}");
}


int compare_clusters(const int *clusters_map1, const int *clusters_map2, int clusters_size) {
    int i;
    int result = 1;

	#pragma omp parallel for reduction(&& : result)
    for (i = 0; i < clusters_size; ++i) {
        if (clusters_map1[i] != clusters_map2[i]) {
 
            result = 0;
        }
    }
    
    return result;
}

/*Se alguém conseguir paralelizar eficientemente essa criação lovecraftiana merece um prêmio nobel*/
double ***km(double **observations, int k, int observations_size, int vector_size) {
	clusters_sizes = (int *) calloc(k, sizeof(int));
	int *clusters_map = (int *) calloc(observations_size, sizeof(int));
	double **cs = initialize(observations, k, observations_size, vector_size);
	
	if (observations_size < k) {
		printf("Could not compute clusters.");

		#pragma omp parallel for
		for (int i = 0; i < k; ++i)
			free(cs[i]);
		free(cs);
		free(clusters_map);
		free(clusters_sizes);
		
		exit(1);
	}
	
	while (1) {
		int *new_clusters_map = partition(observations, cs, k, observations_size, vector_size);
		
		if (compare_clusters(clusters_map, new_clusters_map, observations_size)) {
			double ***clusters = map_clusters(clusters_map, observations, k, observations_size, vector_size);
	
      for (int i = 0; i < k; ++i)
				free(cs[i]);
			free(cs);
			free(clusters_map);
			free(new_clusters_map);
			
			return clusters;
		}
		for (int i = 0; i < k; ++i)
			free(cs[i]);
		free(cs);
		free(clusters_map);
		clusters_map = new_clusters_map;
		cs = re_centroids(clusters_map, observations, k, observations_size, vector_size);
	}
}

/*Na função Re_Centroids explica melhor as infelicidades que tive tentando resolver o segfault que essa função causa
tirar o free do looop não resolve, usar single não resolve, single trava o programa e eu não sei por que */

double *centroid(double **observations, int observations_size, int vector_size) {
	double *vector = (double *) calloc(vector_size, sizeof(double));
	
	for (int i = 0; i < observations_size; ++i) {
		double *temp = vsum(vector, observations[i], vector_size);
		free(vector);
		vector = temp;
	}
	 
	#pragma omp parallel for
	for (int i = 0; i < vector_size; ++i)
		
		#pragma omp critical
		vector[i] /= observations_size;
	
	
	return vector;
}

/*Apesar de vsum ser chamado dentro de centroid, vale mais a pena paralelizar aqui do que em centroids
centroids é uma função do diabo que causa segfault */

double *vsum(const double *vector1, const double *vector2, int vector_size) {
	double *vector = (double *) malloc(sizeof(double) * vector_size);
	
	#pragma omp parallel for 
	for (int i = 0; i < vector_size; ++i)
		vector[i] = vector1[i] + vector2[i];
	
	return vector;
}

/*Não se deve paralelizar vsub e innerprod, já que ambas as funções são chamadas
em partition, que está paralelizada.*/

double *vsub(const double *vector1, const double *vector2, int vector_size) {
	double *vector = (double *) malloc(sizeof(double) * vector_size);
		 
	for (int i = 0; i < vector_size; ++i)
		vector[i] = vector1[i] - vector2[i];
	
	return vector;
}

double innerprod(const double *vector1, const double *vector2, int vector_size) {
    double prod = 0;
    
    for (int i = 0; i < vector_size; ++i)
        prod += vector1[i] * vector2[i];
    
    return prod;
}

double norm(const double *vector, int vector_size) {
	return sqrt(innerprod(vector, vector, vector_size));
}

/* Loved this shuffling random algorithm
 * Source: http://stackoverflow.com/a/5064432
 */
int rand_num(int size) {
	static int *numArr = NULL;
	static int numNums = 0;
	int i, n;
	
	if (size == -22) {
		free(numArr);
		return FREED_RAND;
	}
	
	if (size >= 0) {
		if (numArr != NULL)
			free(numArr);
		
		if ((numArr = (int *) malloc(sizeof(int) * size)) == NULL)
			return ERR_NO_MEM;
		
		#pragma omp parallel for
		for (i = 0; i < size; ++i)
			numArr[i] = i;
		
		numNums = size;
	}
	
	if (numNums == 0)
		return ERR_NO_NUM;
	
	n = rand() % numNums;
	i = numArr[n];
	numArr[n] = numArr[numNums - 1];
	numNums--;
	
	if (numNums == 0) {
		free(numArr);
		numArr = 0;
	}
	
	return i;
}

double **initialize(double **observations, int k, int observations_size, int vector_size) {
	double **centroids = (double **) malloc(sizeof(double *) * k);
	
	srand(time(NULL));
	int r = rand_num(observations_size);
	
	#pragma omp parallel for 
	for (int i = 0; i < k; ++i) {
		centroids[i] = (double *) malloc(sizeof(double) * vector_size);
		for (int j = 0; j < vector_size; ++j) {
			centroids[i][j] = observations[r][j];
			r = rand_num(-1);
		}
	}
	
	rand_num(-22);
	
	return centroids;
}

int *partition(double **observations, double **cs, int k, int observations_size, int vector_size) {
    int *clusters_map = (int *) malloc(sizeof(int) * observations_size);
	double *temp;

	/*vsub e norm(innerprod) são chamadas aqui dentro, não precisamos paralelizar nenhuma das duas funções acima
	há um problema, que vsub está realizando alocação de memória a cada loop, e isso é ineficiente, 
	alguém me ajuda pelo amor de deus*/

	#pragma omp parallel for private(temp)
    for (int i = 0; i < observations_size; ++i) {
        double min_distance = DBL_MAX;
        int centroid = -1;

        for (int c = 0; c < k; ++c) {
            temp = vsub(observations[i], cs[c], vector_size);
            double curr_distance = norm(temp, vector_size);

            if (curr_distance < min_distance) {
                min_distance = curr_distance;
                centroid = c;
            }

            free(temp);
        }

        clusters_map[i] = centroid;
    }

    return clusters_map;
}

/*Eu não consigo paralelizar isso daqui de jeito nenhum, tentei single, critical, tudo,*/

double **re_centroids(int *clusters_map, double **observations, int k, int observations_size, int vector_size) {
	double **centroids = (double **) malloc(sizeof(double *) * k);
	double **temp_arr = (double **) malloc(sizeof(double *) * observations_size);
	int count = 0;

	for (int c = 0; c < k; ++c) {
		for (int i = 0; i < observations_size; ++i) {
			int curr = clusters_map[i];
			
			if (curr == c) {
				temp_arr[count] = observations[i];
				++count;
			}
		}
		centroids[c] = centroid(temp_arr, count, vector_size);
		count = 0;
	}
	
	free(temp_arr);
	
	return centroids;
}

double ***map_clusters(int *clusters_map, double **observations, int k, int observations_size, int vector_size) {
	double ***clusters = (double ***) malloc(sizeof(double **) * k);
	int i;
	#pragma omp parallel for private(i)
	for (i= 0; i < k; ++i)
		clusters[i] = map_cluster(clusters_map, observations, i, observations_size, vector_size);
	
	return clusters;
}

double **map_cluster(const int *clusters_map, double **observations, int c, int observations_size, int vector_size) {
	int count = 0;
	int i;
	int *temp_arr = (int *) malloc(sizeof(int) * observations_size);
	#pragma omp parallel for private(count, i)
	for (i = 0; i < observations_size; ++i) {
		if (clusters_map[i] == c) {
			temp_arr[count] = i;
			++count;
		}
	}
	
	double **cluster = (double **) malloc(sizeof(double *) * count);
	#pragma omp parallel for private(i)
	for (i = 0; i < count; ++i)
		cluster[i] = observations[temp_arr[i]];
	
	free(temp_arr);
	clusters_sizes[c] = count;
	
	return cluster;
}