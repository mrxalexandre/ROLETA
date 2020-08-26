/*
 * SampleDecoder.cpp
 *
 *  Created on: Jan 14, 2011
 *      Author: rtoso
 */

#include "SampleDecoder.h"
#include <iostream>
#include <cmath>

using namespace std;

extern std::vector<std::vector<double> > points;

SampleDecoder::SampleDecoder() { }

SampleDecoder::~SampleDecoder() { }


inline double d(unsigned i, unsigned j) {

	unsigned dim = points[0].size();

	double sum = 0.0;
	for(unsigned d = 0; d<dim; ++d) {
		sum += pow(points[i][d] - points[j][d], 2);
	}

	return sqrt(sum);
}

inline double b(unsigned i, unsigned w, std::vector<std::vector<unsigned> >& cluster) {

	unsigned tam = cluster.size();


	double sum_b = 10000000.0; //FIXME
	for(unsigned t=0; t<tam; ++t) {
		unsigned tam_w = cluster[t].size();

		if(tam_w == 0 or t == w) continue; // Se Cluster w for vazio e t !=w

		double sum = 0.0;
		for(unsigned &j : cluster[t]) {
			sum += d(i,j);
		}
		sum = sum/tam_w;

		if(sum_b > sum and sum != 0.0)
			sum_b = sum;

	}

	return sum_b;
}


inline double a(unsigned i, std::vector<unsigned>& Cw) {

	unsigned tam = Cw.size();
	if(tam == 1) return 0.0;

	double sum = 0.0;
	for(unsigned& j : Cw) {

		if(i == j) continue;

		sum += d(i,j);

	}

	return sum/(tam-1);


}

int closerst_cluster(std::vector<std::vector<unsigned>>& clusters,unsigned point, std::vector< std::pair< double, unsigned > > &ranking, unsigned k) {

	unsigned bestcluster = 0;
	double bestdist = 99999999999.0;

	for(unsigned j = 0; j < k; j++){
		double sum = 0.0;

		for(unsigned c: clusters[j]){
			sum += d(point, c);

		}

		if(sum/clusters[j].size() < bestdist){
			bestdist = sum/clusters[j].size();
			bestcluster = j;

		}

	}

	return bestcluster;


}


double SampleDecoder::decode(const std::vector< double >& chromosome) const {

	unsigned n_points = points.size();

	unsigned k = number_of_cluster;

	std::vector<unsigned> cluster_of_i(n_points + 1,0);
	std::vector<std::vector<unsigned>> clusters(k);
	std::vector< std::pair< double, unsigned > > ranking(chromosome.size());
	for(unsigned i = 0; i < chromosome.size(); ++i) {
		ranking[i] = std::pair< double, unsigned >(chromosome[i], i);

	}
	// Here we sort 'permutation', which will then produce a permutation of [n] in pair::second:
	std::sort(ranking.begin(), ranking.end());

	for(unsigned i=0; i < k; ++i) {
		unsigned point = ranking[i].second; //ranking second é o indice do pnt
		clusters[i].push_back(point);
		cluster_of_i[point] = i;
	}

	for(unsigned i = k; i < n_points; ++i) { //FIXME


		unsigned point = ranking[i].second;
		unsigned c_i = closerst_cluster(clusters, point, ranking,k);
		clusters[c_i].push_back(point);
		cluster_of_i[point] = c_i;


	}

	double silhueta = 0.0;
	for(unsigned i=0; i<n_points; ++i) {

		unsigned c_i = cluster_of_i[i];
		double b_i = b(i, c_i, clusters);
		double a_i = a(i, clusters[c_i]);

		double max = a_i > b_i ? a_i : b_i;

		silhueta += (b_i - a_i)/max;

	}


//		std::cerr << "V:" << (-1)*silhueta/n_points << std::endl;
//
	return (-1)*silhueta/n_points;


/*
	std::vector< std::pair< double, unsigned > > ranking(chromosome.size());
	for(unsigned i = 0; i < chromosome.size(); ++i) {
		ranking[i] = std::pair< double, unsigned >(chromosome[i], i);
	}

	// Here we sort 'permutation', which will then produce a permutation of [n] in pair::second:
	std::sort(ranking.begin(), ranking.end());

	// permutation[i].second is in {0, ..., n - 1}; a permutation can be obtained as follows
	std::list< unsigned > permutation;
	for(std::vector< std::pair< double, unsigned > >::const_iterator i = ranking.begin();
			i != ranking.end(); ++i) {
		permutation.push_back(i->second);
	}

	// sample fitness is the first allele
	return chromosome.front();
*/
}
