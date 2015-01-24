/*
 * SNP.h
 *
 *  Created on: Jan 14, 2013
 *      Author: pickrell
 */

#ifndef SNP_H_
#define SNP_H_

#include <string>
#include <cmath>
#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <sys/stat.h>
#include "gzstream.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>
#include <boost/algorithm/string.hpp>
#include "CmdLine.h"
using namespace std;

class SNP{
public:
	SNP();
	SNP(string, string, int, int, double, double, vector<double>, vector<bool>, vector<int>, vector<vector<pair<int, int> > >);
	SNP(string, string, int, int, int, double, double, vector<double>, vector<bool>, vector<int>, vector<vector<pair<int, int> > >);
	//SNP(string, string, int, double, double, double, double, vector<bool>);
	string id;
	string chr;
	int pos;
	double Z; // Z-score
	double BF; // Bayes factor
	double f; // Minor allele frequency
	double N; // Sample size
	double Ncase, Ncontrol; //Sample sizes for case/control study
	double V; // Variance in estimate of effect size
	vector<double> W; // Prior variances on effect size (averaged);
	int chunknumber;
	float dens;
	vector<bool> annot;
	vector<float> annot_weight;
	vector<int> dists;
	bool condannot;
	void append_distannots(vector<vector<pair<int, int> > >); // convert distances to annotations according to distance models
	int nannot;
	double calc_logBF();
	double calc_logBF_ind(double);
	double approx_v(); // approximate V using f and N
	double approx_v_cc(); //approximate V in case control setting
	double get_x(vector<double>);
	double get_x_cond(vector<double>, double);
	double sumlog(double, double);
};


#endif /* SNP_H_ */
