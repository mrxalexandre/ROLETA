#ifndef ARGPACK_H_
#define ARGPACK_H_

#include <string>
#include <assert.h>

using namespace std;

class ArgPack {
public:
	string inputFile;
	string outputFile;

	unsigned population; // Size of population
	double populationElite; // Fraction of population to be the elite-set
	double populationMutants; // Fraction of population to be replaced by mutants
	double rhoe; // Probability offspring inherits elite parent allele
	unsigned K; // Number of independent populations
	unsigned threads; // Number of threads for parallel decoding
	double time; // Time limit

	unsigned rngSeed; // Seed random number generator
	unsigned exchangeBest; // Exchange best individuals
	unsigned exchangeTop; // Exchange top 2 best
	unsigned generations; // Run for 1000 generations

	bool verbose;

	//------------
	// singleton functions
	//------------

	static const ArgPack &ap() { assert(def_ap_); return *def_ap_; }

	void print() const;

	bool Ok() const;

	ArgPack(int argc, char * const argv []);

	~ArgPack() { assert(def_ap_); def_ap_ = 0; }

private:

	//------------
	// singleton instance
	//------------

	static ArgPack *def_ap_;

	string value;
	string programName;
	bool ok;
	// static void show_usage(string name);
};



#endif /* ARGPACK_H_ */
