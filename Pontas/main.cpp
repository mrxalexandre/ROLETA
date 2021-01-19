#include <iostream>
#include <sstream>
#include <vector>
#include <queue>
#include "../Decoders/SampleDecoder.hpp"
#include "../BRKGA_API/MTRand.hpp"
#include "../BRKGA_API/BRKGA.hpp"
#include "../Timer/BossaTimer.hpp"
#include "ArgPack.hpp"
#include <fstream>
#include <cmath>

std::vector<std::vector<double> > points;

int main(int argc, char* argv[]) {


	ArgPack single_ap(argc, argv);

//ArgPack::ap().time
	unsigned n = 0;			// size of chromosomes
	const unsigned p = ArgPack::ap().population;		// size of population
	const double pe = ArgPack::ap().populationElite;		// fraction of population to be the elite-set
	const double pm = ArgPack::ap().populationMutants;		// fraction of population to be replaced by mutants
	const double rhoe = ArgPack::ap().rhoe;	// probability that offspring inherit an allele from elite parent
	const unsigned K = ArgPack::ap().K;		// number of independent populations
	const unsigned MAXT = ArgPack::ap().threads;	// number of threads for parallel decoding

	const double cutoff_time = ArgPack::ap().time;

	// Reading instance
	std::string s;
	std::string file = ArgPack::ap().inputFile;

	ifstream f(file);

	int n_points, dim;
	f >> n_points >> dim;

	points = std::vector<std::vector<double> > (n_points);
	for(int i=0; i<n_points; ++i)
		points[i] = std::vector<double> (dim);

	for(int i=0;i<n_points;++i) {
		for(int d=0; d<dim; ++d) {
			f >> points[i][d];
		}
	}

	SampleDecoder decoder;			// initialize the decoder


	const long unsigned rngSeed = ArgPack::ap().rngSeed;	// seed to the random number generator
	MTRand rng(rngSeed);				// initialize the random number generator

	n = n_points;

	// initialize the BRKGA-based heuristic


	BossaTimer timer;
	timer.start();

	double bestValue = -1;
	double timerToBest;
	bool verbose = ArgPack::ap().verbose;

	int k_max = sqrt(n_points); // 2 <= k <= sqrt(n)
	int number_pop = k_max - 1;
	std::vector<BRKGA<SampleDecoder, MTRand>*> populations(number_pop);
	vector<SampleDecoder*> vec_decoders(number_pop);
	vector<double> best_values(number_pop);

	int best_population = 0;
	double sum_best= 0.0;

	if(verbose)
		cout << "Initializing populations\n";
	for(int i=0; i< number_pop; ++i) {
		vec_decoders[i] = new SampleDecoder();
		vec_decoders[i]->set_k(i+2);
		populations[i] = new BRKGA<SampleDecoder, MTRand> (n, p, pe, pm, rhoe, *vec_decoders[i], rng, K, MAXT);
		best_values[i] = (-1)*populations[i]->getBestFitness();
		sum_best += best_values[i];

		if(best_values[i] > bestValue) { // bestValue é o melhor global
			bestValue = best_values[i];
			best_population = i;
		}

		if(verbose) {
			std::cout << (-1)*populations[i]->getBestFitness() << " ";
			if( i == number_pop- 1)
				cout << endl;
		}
	}
	cout << "Best solution " << bestValue << " from pop " << best_population << endl;

	unsigned generation = 0;		// current generation
	// const unsigned X_INTVL =  ArgPack::ap().exchangeBest;	// exchange best individuals at every 100 generations
	// const unsigned X_NUMBER = ArgPack::ap().exchangeTop;	// exchange top 2 best
	const unsigned MAX_GENS = ArgPack::ap().generations;	// run for 1000 gens

	int k_idx = (rng.randInt(1)) ? 1 : number_pop - 1;
	int not_improvement = 0;
	// Initial evolution
	do {
		++generation;

		bool improvement = false;
		for (int i = 0; i < 3; ++i) {
			float value = populations[ k_idx - 1 + i ]->getBestFitness();
			populations[ k_idx - 1 + i ]->evolve();
			if ( fabs(populations[ k_idx - 1 + i ]->getBestFitness() - value) > 1e-4 ){
				improvement = true;
				if ( bestValue < -populations[ k_idx - 1 + i ]->getBestFitness() ){
					best_population = k_idx - 1 + i;
					bestValue = -populations[ k_idx - 1 + i ]->getBestFitness();
					timerToBest = timer.getTime();
				}
			}
		}

		if ( !improvement ){
			not_improvement++;
		}

		if ( not_improvement > 50 ){ // TODO: parameter
			not_improvement = 0;
			k_idx += rng.randExc() > 0.5 ? 1 : -1;
			k_idx = max(k_idx, 1);
			k_idx = min(k_idx, number_pop - 1);
		}

	} while (generation < MAX_GENS and timer.getTime() < cutoff_time); // TODO: Verificar isso (quantidade de evoluções)
	
	double ans = 0;
	int clusters = 0;
	for (int i = 0; i < number_pop ; ++i) {
		if( ans < -populations[i]->getBestFitness() ){
			ans = -populations[i]->getBestFitness();
			clusters = i + 2;
		}
	}

	cout << "Best solution " << ans << " with k = " << clusters << endl;
	cout << "Total time = " << timer.getTime() << std::endl;
	cout << "Time to Best ttb = " << timerToBest << std::endl;

	for(int i=0;i<number_pop; ++i) {
		delete populations[i];
		delete vec_decoders[i];
	}


	return 0;
}
