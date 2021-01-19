#include <iostream>
#include <sstream>
#include <vector>
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

	for(int i = 0; i < n_points; ++i) {
		for(int d = 0; d < dim; ++d) {
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
	double sum_best = 0.0;

	if(verbose)	cout << "Initializing populations\n";
	
	for(int i = 0; i < number_pop; ++i) {
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
			cout << -populations[i]->getBestFitness() << " ";
			if( i == number_pop - 1) cout << endl;
		}
	}
	
	const unsigned MAX_GENS = ArgPack::ap().generations;	// run for 1000 gens

	unsigned generation = 0;		// current generation
	do {
		++generation;
		if(verbose)	cout << "Envolving generation " << generation << " ";
		
		double q_sum = 0;
		vector<double> q(number_pop, 0);
		for (unsigned i = 0; i < (unsigned)number_pop ; ++i) {
			for (unsigned j = 0; j < populations[i]->getPe() ; ++j) {
				q[i] -= populations[i]->getPopulation().getFitness(j);
			}
			q_sum += q[i]; 
		}
		
		double r = rng.randExc() * q_sum;
		int pop_pick = -1;
		double acum = 0;
		
		for (int i = 0; i < number_pop ; ++i) {
			acum += q[i];
			if(r < acum) {
				pop_pick = i;
				break;
			}
		}
		if(pop_pick == -1) pop_pick = number_pop -1;

		if(verbose)	cout << "We pick pop number " << pop_pick  << " ";

		populations[pop_pick]->evolve();	// evolve the population for one generation
		
		if( -populations[pop_pick]->getBestFitness() > bestValue) {
			timerToBest = timer.getTime();
			bestValue = -populations[pop_pick]->getBestFitness();
			best_population = pop_pick;
			cout << "New best " << bestValue << " from pop " << best_population << endl;
		}

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
