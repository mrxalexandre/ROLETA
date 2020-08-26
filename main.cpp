#include <iostream>
#include <sstream>
#include <vector>
#include "SampleDecoder.h"
#include "MTRand.h"
#include "BRKGA.h"
#include "bossa_timer.h"
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
	//dstd::getline(std::cin, s);
	//std::stringstream st(s);
	f >> n_points >> dim;

	points = std::vector<std::vector<double> > (n_points);
	for(int i=0; i<n_points; ++i)
		points[i] = std::vector<double> (dim);


	//for(int d=0;d<dim;++d) {
		//CLUSTER_L[d] = 1000;
		//CLUSTER_U[d] = -1000;
	//}

	for(int i=0;i<n_points;++i) {

		//std::getline(std::cin, s);
		//std::stringstream st(s);
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
	const unsigned X_INTVL =  ArgPack::ap().exchangeBest;	// exchange best individuals at every 100 generations


	const unsigned X_NUMBER = ArgPack::ap().exchangeTop;	// exchange top 2 best
	const unsigned MAX_GENS = ArgPack::ap().generations;	// run for 1000 gens
	do {

		++generation;


		if(verbose)
			cout << "Envolving generation " << generation << " ";
	/*   Metodo da roleta, escolho um random e vejo em que parte ele
		 se encontra de acordo com os best de cada população

							| random number
							|
							X
			+--------------+-----+----+--+--++---+-----+
			|              |     |    |  |  ||   |     |
			+--------------+-----+----+--+--++---+-----+
			0             b0     b1                    sum_best         */
		double r = rng.randDblExc(sum_best);
		int pop_pick = -1;

		double acum = 0.0;

		// 1, 3, 2, 4, 8, 5, 9           0 1 4 6 10 18 23 32

		for (int i=0;i<number_pop; ++i) {
			acum += best_values[i];
			if(r < acum) {
				pop_pick = i;
				break;
			}
		}

		if(pop_pick == -1) // In case of bug
			pop_pick = number_pop -1;

		if(verbose)
			cout << "We pick pop number " << pop_pick  << " ";


		populations[pop_pick]->evolve();	// evolve the population for one generation

		double best_temp = (-1)* populations[pop_pick]->getBestFitness();
		if(verbose)
			cout << best_temp << endl;


		if( best_temp > best_values[pop_pick]) // If improve the pop best
		{
			sum_best = sum_best + best_temp - best_values[pop_pick]; // update sum of best (sum new and subtract old)
			best_values[pop_pick] = best_temp;
		}
		if(best_temp > bestValue) {
			timerToBest = timer.getTime();
			bestValue = best_temp;
			best_population = pop_pick;

			cout << "New best " << bestValue << " from pop " << best_population << endl;
		}

		//if((++generation) % X_INTVL == 0) { // VEJA se usa isso
		//	algorithm.exchangeElite(X_NUMBER);	// exchange top individuals
		//}



	} while (generation < MAX_GENS and timer.getTime() < cutoff_time);
	timer.pause();

		cout << "Best solution " << bestValue << " from pop " << best_population << " k = " << best_population +2 << endl;
		std::cout << "Total time = " << timer.getTime() << std::endl;
		std::cout << "Time to Best ttb = " << timerToBest << std::endl;


	for(int i=0;i<number_pop; ++i) {
		delete populations[i];
		delete vec_decoders[i];
	}


	return 0;
}
