/*
 * SNPs.h
 *
 *  Created on: Jan 14, 2013
 *      Author: pickrell
 */

#ifndef SNPS_H_
#define SNPS_H_

#include "SNP.h"
#include "fgwas_params.h"
using namespace std;

typedef double LLKFunction(const gsl_vector *, void *);

class SNPs{
public:
	SNPs();
	SNPs(Fgwas_params *);
	Fgwas_params *params;
	vector<SNP> d;

	//snp annotations
	vector<double> snppri;
	vector<double> snppost;
	vector<double> lambdas;
	vector<string> annotnames;
	vector<string> chrnames;
	int nannot;
	vector<vector<pair<int, int> > > dmodels; // hold the distance models
	double condlambda; //for conditional analysis
	
	//segment annotations
	double segpi;
	int nsegannot;
	vector<string> segannotnames;
	vector<pair<int, int> > segments;
	vector<pair<int, int> > chrsegments;
	vector<double> seglambdas;
	vector<vector<bool> > segannot;
	vector<double> segpriors;
	void init_segpriors();
	void set_segpriors();

	//funtions for working with distance annotations
	void append_dannotnames(string, vector<pair<int, int> >);
	vector<pair<int, int> > read_dmodel(string);

	//10-fold cross-validation
	double cross10(bool);
	vector<set<int> > make_cross10();

	double phi, resphi;
	void load_snps(string, double, vector<string>);
	void load_snps_z(string, vector<double>, vector<string>, vector<string>, vector<string>);
	void print();
	void print(string, string);

	void make_segments(int);
	void make_segments(string);
	map<string, vector<pair<int, int> > > read_bedfile(string);
	void make_chrsegments();
	void make_segments_finemap();
	void print_segments();
	void print_chrsegments();
	void set_priors(int);
	void set_priors();
	void set_priors_cond(int);
	void set_priors_cond();

	void set_post(int);
	void set_post();
	void GSL_optim();
	void GSL_optim_ridge();
	void GSL_xv_optim(set<int> toskip, bool penalize);
	void GSL_optim(LLKFunction* pLLKFunc, set<int> toskip, bool penalize);
	
	double llk();
	double llk(set<int> skip, bool penalize);
	double llk(int which);
	
	double data_llk;
	double sumlog(double, double);
	void print_segprobs(string);
	void optimize_segpi();
	void optimize_condlambda();
	void optimize_l0();
	vector<pair<pair<int, int>, pair<double, double> > > get_cis();
	pair<pair<int, int>, pair<double, double> > get_cis_segpi();
	pair<pair<int, int>, pair<double, double> > get_cis_condlambda();
	pair<pair<int, int>, pair<double, double> > get_cis_lambda(int);
	pair<pair<int, int>, pair<double, double> > get_cis_seglambda(int);

	int golden_section_segpi(double, double, double, double);
	int golden_section_condlambda(double, double, double, double);
	int golden_section_segpi_ci(double, double, double, double, double, int *);
	int golden_section_lambda_ci(double, double, double, double, int, double);
	int golden_section_condlambda_ci(double, double, double, double, double);
	int golden_section_seglambda_ci(double, double, double, double, int, double);
	int golden_section_l0(double, double, double, double);

	void check_input();
	void check_string2digit(string);
	
	//**************************** OPTIMIZATION ****************************//
	// Vector that stores a cache of the sum of lambda values returned by SNP::get_x.
	// Because combinations of annotations re-occur across SNPs, rather than recomputing 
	// this sum for each SNP, we cache. Each SNP has an index into this cache that depends
	// on the annotations of the particular SNP.
	// This vector has length equal to the number of distinct annotation combinations
	// across all SNPs.
	vector<double> snpXcache;
	// This vector has length equal to the number of distinct annotation combinations.
	// It stores the unique combination of annotation 1's/0's.
	vector<vector<bool> > cacheAnnots;
	
	// A vector of length (number of SNPs), with the index into the cache for the SNP.
	// Indices in this correspond with indices in the vector "d" above.
	vector<int> snpCacheIndex;

	// A vector of the Bayes factor for each SNP. These are d[i].BF.
	vector<double> BFs;
	
	// A vector of snppri + BF - computed in setpriors as an optimization.
	vector<double> snppriBFsum;
	
	void init_priors_optimization();
	void init_BF_optimization();
	void precompute_prior_sums();
	//**********************************************************************//
};

struct GSL_params{
        SNPs *d;
        int which;
        set<int> toskip;
        bool penalize;
};
extern double GSL_llk(const gsl_vector *, void *GSL_params);
extern double GSL_llk_ridge(const gsl_vector *, void *GSL_params);
extern double GSL_llk_fine(const gsl_vector *, void *GSL_params);
extern double GSL_llk_ridge_fine(const gsl_vector *, void *GSL_params);
extern double GSL_llk_xv(const gsl_vector *, void *GSL_params);

#endif /* SNPS_H_ */
