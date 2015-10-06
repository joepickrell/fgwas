/*
 * test.cpp
 *
 *  Created on: Jan 14, 2013
 *      Author: pickrell
 */

#include "SNPs.h"
#include "fgwas_params.h"
#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

int main(){

	Fgwas_params p;
	//p.wannot.push_back("wgEncodeDukeDnase8988T");
	//p.wannot.push_back("wgEncodeDukeDnaseAoSMC");
	p.cc = true;
	vector<double> tmp;
	tmp.push_back(0.05);
	p.infile = "/Users/jkpickrell/Documents/workspace/fgwas/testing/testin.gz";
	p.V = tmp;
	SNPs s(&p);
	//s.GSL_optim();
	for (int i = 1; i < 100; i++){
		double tmps = (double) i / (double) 100;
		s.segpi = tmps;
		s.set_segpriors();
		cout << i << " "<< s.llk() << "\n";
	}


	return 0;
}

