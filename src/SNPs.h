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
	double llk_xv(set<int>, bool);

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
	void GSL_xv_optim(set<int>, bool);
	void GSL_optim_fine();
	void GSL_optim_ridge();
	void GSL_optim_ridge_fine();
	double llk(int);
	double llk();
	double llk_ridge();

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
