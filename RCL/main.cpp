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
			timerToBest = timer.getTime();
		}

		if(verbose) {
			std::cout << (-1)*populations[i]->getBestFitness() << " ";
			if( i == number_pop- 1)
				cout << endl;
		}
	}
	cout << "Best solution " << bestValue << " from pop " << best_population << endl;

	unsigned generation = 0;		// current generation
	const unsigned X_INTVL =  ArgPack::ap().exchangeBest;	// exchange best individuals at every 100 generations
	const unsigned X_NUMBER = ArgPack::ap().exchangeTop;	// exchange top 2 best
	const unsigned MAX_GENS = ArgPack::ap().generations;	// run for 1000 gens
	// Initial evolution
	do {
		++generation;
		if(verbose)	cout << "Envolving generation " << generation << " ";
		for (int i = 0; i < number_pop; ++i) {
			populations[i]->evolve();	// evolve the population for one generation
			// Pegar o melhor valor até então e setar o tempo
			if( (-1)*populations[i]->getBestFitness() > bestValue ){
				bestValue = (-1)*populations[i]->getBestFitness();
				best_population = i;
				timerToBest = timer.getTime();
			}
		}
	} while (generation < MAX_GENS/10 and timer.getTime() < cutoff_time/10); // TODO: Verificar isso (quantidade de evoluções)
	
	int Kp = populations.size()/2;

	priority_queue<pair<double,int>> ranking_population;

	// Ordering the population
	for (int i = 0; i < number_pop ; ++i) {
		ranking_population.push( { populations[i]->getBestFitness(), i } );
	}

	// Get K best population
	vector<int> KBestPopulations(Kp);
	for (int i = 0; i < Kp; ++i) {
		pair<double,int> value = ranking_population.top();
		ranking_population.pop();
		KBestPopulations[i] = value.second;
	}
	// ranking_population.clear();

	// Evolving the K best population
	do {
		++generation;
		if(verbose)	cout << "Envolving generation " << generation << " ";
		for (int i = 0; i < Kp ; ++i) {
			populations[ KBestPopulations[i] ]->evolve();	// evolve the population for one generation
			// Pegar o melhor valor até então e setar o tempo
			if( (-1)*populations[i]->getBestFitness() > bestValue ){
				bestValue = (-1)*populations[i]->getBestFitness();
				best_population = i;
				timerToBest = timer.getTime();
			}
		}
	} while (generation < MAX_GENS and timer.getTime() < cutoff_time); // TODO: Verificar isso (quantidade de evoluções)

	timer.pause();

	double ans = 0;
	int clusters = 0;
	for (int i = 0; i < number_pop ; ++i) {
		if( ans < -populations[i]->getBestFitness() ){
			ans = -populations[i]->getBestFitness();
			clusters = i + 2;
		}
	}

	cout << "Best solution " << ans << " k = " << clusters << endl;
	std::cout << "Total time = " << timer.getTime() << std::endl;
	std::cout << "Time to Best ttb = " << timerToBest << std::endl;

	// Show answer:
	vector<double> ch = populations[clusters-2]->getBestChromosome();
	// vector<vector<int>> allocate = vec_decoders[clusters-2]->answer(ch);

	// for (unsigned i = 0; i < allocate.size() ; ++i) {
	// 	cout << "Cluter " << i+1 << " is : " << allocate[i][0] << endl;
	// 	cout << "Points:";
	// 	for (unsigned j = 1; j < allocate[i].size(); ++j) {
	// 		cout << " " << allocate[i][j];
	// 	}
	// 	cout << endl;
	// }

	for(int i=0;i<number_pop; ++i) {
		delete populations[i];
		delete vec_decoders[i];
	}


	return 0;
}
