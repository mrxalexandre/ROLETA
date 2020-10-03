#include "ArgPack.hpp"

#include <string>
#include <cstring>
#include <unistd.h> // for getopt
#include <iostream>

using namespace std;

extern int optind;
ArgPack *ArgPack::def_ap_ = 0;

static void show_usage(string name)
{
    cerr << "Usage: " << name << " <option(s)> SOURCES"
			<< "Options:\n"
			<< "\t-h,--help\t\tShow this help message\n"
			<< "\t-i,--input Input file\tSpecify the destination path\n"
			<< "\t-v,--verbose \n"
			<< "\t-s,--seed Pseudo number generator seed\n"
			<< "\t-p,--population Population\tSize of population\n"
			<< "\t-e,--population-elite Population elite-set\tFraction of population to be the elite-set\n"
			<< "\t-m,--population-mutants Population mutants\tFraction of population to be replaced by mutants\n"
			<< "\t-r,--rhoe Inherits allele\tProbability offspring inherits elite parent allele\n"
			<< "\t-K,--K Independent populations\tNumber of independent populations\n"
			<< "\t-t,--threads Threads\tNumber of threads for parallel decoding\n"
			<< "\t-b,--exchangeBest\texchange the best individuals \n"
			<< "\t-x,--exchangeTop \t exchange the top individuals\n"
			<< "\t-T,--time Time\tCPU-time (in secondes)\n"
			<< "\t-g,--generations Generations\tNumber of threads for parallel decoding\n"
			;
}

ArgPack::ArgPack(int argc, char * const argv []) :
	inputFile(""),
	outputFile(""),

	population(100),
	populationElite(0.2),
	populationMutants(0.1),
	rhoe(0.7),
	K(1),
	threads(1),
	time(100000.0),

	rngSeed(0),
	exchangeBest(100),
	exchangeTop(2),
	generations(1000),
	verbose(false){


	assert(!def_ap_);
	def_ap_ = this;
	programName = argv[0];
	ok = false;
	bool error = false;

	for (int i = 1; i < argc; ++i) {
		string arg = argv[i];
		if( arg == "-h" or arg == "--help" ){
			show_usage(programName);
			exit(0);
		} else if(arg == "-i" or arg == "--input") {
			try {
				inputFile = argv[++i];
			}
			catch (...)	{
				cerr << "Unknown exception caught!" << endl;
				error = true;
			}
		} else if (arg == "-v" or arg == "--verbose"){
			try {
				verbose = true;
			}
			catch (...)	{
				cerr << "Unknown exception caught!" << endl;
				error = true;
			}
		} else if (arg == "-o" or arg == "--output"){
			try {
				outputFile = argv[++i];
			}
			catch (...)	{
				cerr << "Unknown exception caught!" << endl;
				error = true;
			}
		} else if(arg == "-p" or arg == "--population"){
			try {
				value = argv[++i];
				population = stoi(value, nullptr, 10);
			}
			catch (...)	{
				error = true;
				cerr << "Unknown exception caught!" << endl;
			}
		} else if(arg == "-e" or arg == "--population-elite"){
			try {
				value = argv[++i];
				populationElite = stof(value);
			}
			catch (...)	{
				error = true;
				cerr << "Unknown exception caught!" << endl;
			}
		} else if(arg == "-m" or arg == "--population-mutants"){
			try {
				value = argv[++i];
				populationMutants = stof(value);
			}
			catch (...)	{
				error = true;
				cerr << "Unknown exception caught!" << endl;
			}
		} else if(arg == "-r" or arg == "--rhoe"){
			try {
				value = argv[++i];
				rhoe = stof(value);
			}
			catch (...)	{
				cerr << "Unknown exception caught!" << endl;
				error = true;
			}
		} else if(arg == "-s" or arg == "--seed"){
			try {
				value = argv[++i];
				rngSeed = strtoul(value.c_str(),NULL,10);
			}
			catch (...)	{
				cerr << "Unknown exception caught!" << endl;
				error = true;
			}
		} else if(arg == "-t" or arg == "--threads"){
			try {
				value = argv[++i];
				threads = stof(value);
			}
			catch (...)	{
				error = true;
				cerr << "Unknown exception caught!" << endl;
			}
		} else if(arg == "-T" or arg == "--time"){
			try {
				value = argv[++i];
				time = stof(value);
			}
			catch (...)	{
				error = true;
				cerr << "Unknown exception caught!" << endl;
			}
		} else if(arg == "-K" or arg == "--K"){
			try {
				value = argv[++i];
				K = stoi(value, nullptr, 10);
			}
			catch (...)	{
				error = true;
				cerr << "Unknown exception caught!" << endl;
			}
		} else if(arg == "-T" or arg == "--threads"){
			try {
				value = argv[++i];
				threads = stoi(value, nullptr, 10);
			}
			catch (...)	{
				error = true;
				cerr << "Unknown exception caught!" << endl;
			}
		} else if(arg == "-g" or arg == "--generations"){
			try {
				value = argv[++i];
				generations = stoi(value, nullptr, 10);
			}
			catch (...)	{
				error = true;
				cerr << "Unknown exception caught!" << endl;
			}
		} else if(arg == "-b" or arg == "--exchange-best"){
			try {
				value = argv[++i];
				exchangeBest = stoi(value, nullptr, 10);
			}
			catch (...)	{
				error = true;
				cerr << "Unknown exception caught!" << endl;
			}
		} else if(arg == "-x" or arg == "--exchange-top"){
			try {
				value = argv[++i];
				exchangeTop = stoi(value, nullptr, 10);
			}
			catch (...)	{
				error = true;
				cerr << "Unknown exception caught!" << endl;
			}
		} else {
			error = true;
			cerr << "Unexpected parameter!" << endl;
		}
	}

	if( !error ){
		if( inputFile == "" ){
			error = true;
			//cerr << "Missing parameters!" << endl;
		}
		if( !error ){
			if( outputFile == "" ){
				int pos = inputFile.find_last_of('.');
				if( pos < 0 ){
					outputFile = inputFile + "_H.out";
				} else {
					outputFile = inputFile.substr(0,pos) + "_H.out";
				}
			}
			ok = true;
		}
	}
}

bool ArgPack::Ok() const {
	return ok;
}

void ArgPack::print() const{
	// cout << "Parameters:\n";
	// cout << "Instance: " << input_name << endl;
	// if(primal_bound == 0)
	// 	cout << "\tNo primal bound\n";
	// else
	// 	cout << "Primal bound = " << primal_bound << endl;
	// if(seed == 0)
	// 	cout << "\tRandom seed\n";
	// else
	// 	cout << "\tSeed = " << seed << endl;
	// if(time_limit == 0)
	// 	cout << "\tNo time limit\n";
	// else
	// 	cout << "\tTime limit = " << time_limit << endl;
	// cout << "\tI_R = " << I_R << endl;
	// cout << "\tI_ILS = " << I_ILS << endl;
	// cout << "\tnMovesPert   = " << nMovesPert << endl;
	// cout << "\tnPert   = " << nPert << endl;
	// cout << "\ttypePert = " << typePert << endl;
	// cout << "\tneighSelector = " << neighSelector << endl;
	// cout << "\tgreedyAlg = " << greedyAlg << endl;
	// cout << "\tILS_DP = " << ILS_DP << endl;
	// cout << "\tnaiveEval = " << naiveEval << endl;
}
