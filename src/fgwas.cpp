/*
 * fgwas.cpp
 *
 *  Created on: Apr 2, 2013
 *      Author: pickrell
 */




#include "SNPs.h"
#include "fgwas_params.h"
using namespace std;


void printopts(){
        cout << "\nfgwas v. 0.3.6\n";
        cout << "by Joe Pickrell (jkpickrell@nygenome.org)\n\n";
        cout << "-i [file name] input file w/ Z-scores\n";
        cout << "-o [string] stem for names of output files\n";
        cout << "-w [string] which annotation(s) to use. Separate multiple annotations with plus signs\n";
        cout << "-dists [string:string] the name of the distance annotation(s) and the file(s) containing the distance model(s)\n";
        cout << "-k [integer] block size in number of SNPs (5000)\n";
        cout << "-bed [string] read block positions from a .bed file\n";
        cout << "-v [float] variance of prior on normalized effect size. To average priors, separate with commas (0.01,0.1,0.5)\n";
        cout << "-p [float] penalty on sum of squared lambdas, only relevant for -print (0.2)\n";
        //cout << "-mse input is in mean/standard error format (default is Z-score, sample size)\n";
        cout << "-print print the per-SNP output\n";
        //cout << "-drop [string] chromosome to drop (none)\n";
        cout << "-xv do 10-fold cross-validation\n";
        cout << "-dens [string float float] model segment probability according to quantiles of an annotation\n";
        cout << "-cc this is a case-control study, which implies a different input file format\n";
        cout << "-fine this is a fine mapping study, which implies a different input file format\n";
        cout << "-onlyp only do optimization under penalized likelihood\n";
        cout << "-cond [string] estimate the effect size of an annotation conditional on the others in the model\n";
        cout << "\n";
}


int main(int argc, char *argv[]){
	Fgwas_params p;

    CCmdLine cmdline;
    if (cmdline.SplitLine(argc, argv) < 1){
        printopts();
        exit(1);
    }
    if (cmdline.HasSwitch("-i")) p.infile = cmdline.GetArgument("-i", 0).c_str();
    else{
        printopts();
        exit(1);
    }
    if (cmdline.HasSwitch("-o")) p.outstem = cmdline.GetArgument("-o", 0);
    if (cmdline.HasSwitch("-k")) p.K = atoi(cmdline.GetArgument("-k", 0).c_str());
    if (cmdline.HasSwitch("-bed")) {
    	p.bedseg = true;
    	p.segment_bedfile = cmdline.GetArgument("-bed", 0);
    }
    if (cmdline.HasSwitch("-cc")) {
    	p.cc = true;
    }
    if (cmdline.HasSwitch("-v")) {
    	p.V.clear();
    	vector<string> strs;
    	string s = cmdline.GetArgument("-v", 0);
    	boost::split(strs, s ,boost::is_any_of(","));
    	for (int i  = 0; i < strs.size(); i++) {
    		p.V.push_back( atof(strs[i].c_str()) );
    	}
    }
    if (cmdline.HasSwitch("-p")) p.ridge_penalty = atof(cmdline.GetArgument("-p", 0).c_str());
    if (cmdline.HasSwitch("-xv")) p.xv = true;
    if (cmdline.HasSwitch("-mse")) p.zformat = false;
    if (cmdline.HasSwitch("-print")) p.print = true;
    if (cmdline.HasSwitch("-onlyp")) p.onlyp = true;
    if (cmdline.HasSwitch("-cond")){
    	p.cond = true;
    	p.testcond_annot = cmdline.GetArgument("-cond", 0);
    }
    if (cmdline.HasSwitch("-w")){
    	vector<string> strs;
    	string s = cmdline.GetArgument("-w", 0);
    	boost::split(strs, s ,boost::is_any_of("+"));
    	for (int i  = 0; i < strs.size(); i++) {
    		p.wannot.push_back( strs[i] );
    	}
    }
    if (cmdline.HasSwitch("-dists")){
     	vector<string> strs;
     	string s = cmdline.GetArgument("-dists", 0);
     	boost::split(strs, s ,boost::is_any_of("+"));
     	for (int i  = 0; i < strs.size(); i++) {
     		vector<string> strs2;
     		boost::split(strs2, strs[i], boost::is_any_of(":"));
     		p.dannot.push_back( strs2[0] );
     		p.distmodels.push_back(strs2[1]);
     	}
     }
    if (cmdline.HasSwitch("-drop")){
    	p.dropchr = true;
    	string s = cmdline.GetArgument("-drop", 0);
    	p.chrtodrop = s;
    }

    if (cmdline.HasSwitch("-fine")) p.finemap = true;
    if (cmdline.HasSwitch("-dens")) {
    	p.segannot.push_back(cmdline.GetArgument("-dens", 0));
    	p.loquant = atof(cmdline.GetArgument("-dens", 1).c_str());
    	p.hiquant = atof(cmdline.GetArgument("-dens", 2).c_str());
    }

	SNPs s(&p);

	//if doing unpenalized optimization
	if (!p.onlyp){
		s.GSL_optim();

		// conditional analysis
		if (p.cond){
			string llkoutfile = p.outstem+".llk";
			ofstream lkout(llkoutfile.c_str());
			lkout << "ln(lk) [w/o cond]: "<<  s.llk() << "\n";
			s.optimize_condlambda();
			double cl = s.condlambda;
			lkout << "ln(lk) [w cond]: "<< s.llk() << "\n";

			string outparam = p.outstem+".params";
			ofstream out(outparam.c_str());
			pair<pair<int, int>, pair<double, double> > cond_cis = s.get_cis_condlambda();
			out << "pi " << s.segpi << "\n";
			for (int i = 0; i < s.seglambdas.size(); i++) out << s.segannotnames[i] << " " << s.seglambdas[i]<< "\n";
			for (int i = 0; i < s.lambdas.size(); i++) out << s.annotnames[i] << " " << s.lambdas[i]<< "\n";
			out << p.testcond_annot << " ";
			if (cond_cis.first.first == 0)  out << cond_cis.second.first << " " << cl<< " ";
			else out << "fail "<< cl << " ";
			if (cond_cis.first.second == 0) out << cond_cis.second.second << "\n";
			else out << "fail\n";
		}

		//standard analysis
		else{
			string llkoutfile = p.outstem+".llk";
			ofstream lkout(llkoutfile.c_str());
			int np = s.lambdas.size() + s.seglambdas.size();

			lkout << "ln(lk): "<<  s.llk() << "\n";
			lkout << "nparam: "<< np << "\n";
			lkout << "AIC: "<< 2.0* (double) np - 2* s.llk() << "\n";
			//cout << "getting cis\n"; cout.flush();
			vector<pair<pair<int, int>, pair<double, double> > > cis = s.get_cis();
			//cout << "ok\n"; cout.flush();
			string outparam = p.outstem+".params";
			ofstream out(outparam.c_str());
			out << "parameter CI_lo estimate CI_hi\n";

			//not fine-mapping (there are segment-level annotations)
			if (!p.finemap){
				pair<pair<int, int>,  pair<double, double> > sci = cis[0];
				out << "pi_region ";
				if (sci.first.first == 0) out<< sci.second.first << " " << s.segpi << " ";
				else if (sci.first.first == 2) out << "<"<< sci.second.first << " "<< s.segpi << " ";
				else out << "fail "<< s.segpi << " ";
				if (sci.first.second == 0) out << sci.second.second << "\n";
				else if (sci.first.second == 2) out << ">"<< sci.second.second << "\n";
				else out << "fail\n";


				for (int i = 0; i < s.seglambdas.size(); i++){
					out << s.segannotnames[i] << "_ln ";
					int index = i+1;
					if (cis[index].first.first == 0)  out<< cis[index].second.first << " " << s.seglambdas[i]<< " ";
					else if (cis[index].first.first == 2) out<< "<"<< cis[index].second.first << " " << s.seglambdas[i]<< " ";
					else out << "fail "<< s.seglambdas[i] << " ";
					if (cis[index].first.second == 0) out << cis[index].second.second << "\n";
					else if (cis[index].first.second == 2) out << ">"<<cis[index].second.second << "\n";
					else out << "fail\n";
				}
			}

			//print annotations
			for (int i = 0; i < s.lambdas.size(); i++){
				out << s.annotnames[i] << "_ln ";
				int index = i+1+s.seglambdas.size();
				if (p.finemap) index = i;
				if (cis[index].first.first == 0)  out<< cis[index].second.first << " " << s.lambdas[i]<< " ";
				else if (cis[index].first.first == 2) out<< "<"<< cis[index].second.first << " " << s.lambdas[i]<< " ";
				else out << "fail "<< s.lambdas[i] << " ";
				if (cis[index].first.second == 0) out << cis[index].second.second << "\n";
				else if (cis[index].first.second == 2) out << ">"<<cis[index].second.second << "\n";
				else out << "fail\n";
			}
		}
	}

	// penalized likelihood
	if ( (p.print || p.xv) && !p.cond) {
		s.GSL_optim_ridge();
		string outridge = p.outstem+".ridgeparams";
		ofstream outr(outridge.c_str());
		outr << "ridgeparam: "<< p.ridge_penalty << "\n";
		for (int i = 0; i < s.seglambdas.size(); i++) outr << s.segannotnames[i] << " " << s.seglambdas[i]<< "\n";
		for (int i = 0; i < s.lambdas.size(); i++) outr << s.annotnames[i] << " " << s.lambdas[i]<< "\n";
		if (p.xv){
			double xvlk = s.cross10(true);
			outr << "X-validation penalize ln(lk): "<< xvlk << "\n";
		}
		// print the PPAs
		if (p.print) s.print(p.outstem+".bfs.gz", p.outstem+".segbfs.gz");
	}


	return 0;
}

