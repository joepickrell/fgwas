/*
 * test.cpp
 *
 *  Created on: Jan 14, 2013
 *      Author: pickrell
 */

//#include "SNPs.h"
//#include "fgwas_params.h"
#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

int main(){
	vector<double> tmp;
	tmp.push_back(1.23); tmp.push_back(5.45);
	cout << tmp.size() << " "<< log(tmp.size()) << "\n";
	/*
	Fgwas_params p;
	//p.wannot.push_back("wgEncodeDukeDnase8988T");
	//p.wannot.push_back("wgEncodeDukeDnaseAoSMC");
	//p.cc = false;
	p.K = 5;
	//p.finemap = false;
	p.segment_bedfile = "/Users/jkpickrell/projects/GWAS/test.bed";
	p.bedseg = true;
	p.infile = "/Users/jkpickrell/projects/GWAS/t.gz";
	string bedfile = "/Users/jkpickrell/projects/GWAS/test.bed";
	SNPs s(&p);
	for (vector<pair<int, int> >::iterator it = s.segments.begin(); it != s.segments.end(); it++){
		int st = it->first;
		int sp = it->second;
		cout << st << " "<< sp << "\n";
		cout << s.d.at(st).chr << " "<< s.d.at(st).pos << " "<< s.d.at(sp-1).chr << " "<< s.d.at(sp-1).pos << "\n";
	}
	*/
	//map<string, vector<pair<int, int> > > t = s.read_bedfile(bedfile);
	//for (map<string, vector<pair<int, int> > >::iterator it = t.begin(); it != t.end(); it++){
	//	for (vector<pair<int, int> >::iterator it2 = it->second.begin(); it2 != it->second.end(); it2++){
	//		cout << it->first << " "<< it2->first << " "<< it2->second << "\n";
	//	}
	//}
	//p.dannot.push_back("tssdist");
	//p.segannot.push_back("tssdist");
	//p.loquant = 0.33;
	//p.hiquant = 0.67;
	//p.distmodels.push_back("dist_model");
	//string test = "1340.5";
	//cout << atoi(test.c_str()); cout.flush();
	//p.infile = "/Users/pickrell/Documents/workspace/GWAS/test/hb.wpos.wf.COMBINED.wdist.gz";
	//p.infile = "/Users/pickrell/Documents/workspace/GWAS/test/hb_testin.gz";
	//SNPs s(&p);
	//s.GSL_optim();
	//s.cross10(false);
	//s.GSL_optim();
	/*
	for (int i = 0; i < s.segments.size(); i++){
		pair<int, int> seg = s.segments[i];
		int st = seg.first;
		int sp = seg.second;
		int sum  = 0;
		int total = 0;
		for (int j = st ;j < sp; j++){
			sum+= s.d[j].dens;
			total++;
		}
		double mean = (double) sum / (double) total;
		cout << i << " "<< mean << " "<< s.d[st].pos << " "<< s.d[sp].pos << " "<< s.d[st].chr << "\n";

	}
	*/
	//s.print_segments();
	//s.GSL_optim();
	//s.print("testout.bfs.gz");

	return 0;
}

