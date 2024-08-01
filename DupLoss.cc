/*
Copyright (C) 2024 Mukul S. Bansal (mukul.bansal@uconn.edu).
Based on open-source code originally written by Andre Wehe and 
Mukul S. Bansal for the DupTree software package.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
 
You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


bool LIMIT_LOSSES = false;
float WEIGHTED_RECON_COST = 0;

#include <signal.h>
#include "rmq.h"
#include "rmq.c"
#include "common.h"
#include "argument.h"
#include "input.h"
#include "node.h"
#include "tree.h"

namespace buildtree {
//#include "buildtree-treeset-random.h"
#include "buildtree-treeset-leafadd.h"
#include "heuristic-leafadd.h"
#include "buildtree-heuristicsecondfast.h"
//#include "buildtree-heuristicsecond.h"
}

namespace gtpspr {
#include "gtp-spr-treeset.h"
#include "gtp-spr-heuristic.h"
#include "gtp-spr-simpleheuristicrandom.h"
}


namespace genedupreport {
#include "report-treeset.h"
#include "report.h"
}


int main(int argc, char *argsv[]) {
	initialization();
	signal(SIGINT, interruptFunction);

	Argument::add(argc, argsv);

	// random seed given
	unsigned int randomseed = unsigned(time(NULL)) % 0xFFFF;
	{
		const Argument *arg = Argument::find("--seed");
		if (arg != NULL) arg->convert(randomseed);
		srand(randomseed);
	}



	// score
	bool score_flag = true;

	// report
	bool report_flag = true;

	ReRoot reroot;

	// output format
	Format oformat;
	{
		const Argument *arg = Argument::find("--oformat");
		if (arg == NULL) oformat = NEWICK;
		else {
			string str;
			arg->convert(str);
			if (str == "newick") oformat = NEWICK;
			else
			if (str == "nexus") oformat = NEXUS;
			else EXCEPTION("--oformat has a wrong argument");
		}
	}

	// output input gene trees
	bool genetreesFlag = false;
	{
		const Argument *arg = Argument::find("--genetrees");
		if (arg != NULL) {
			genetreesFlag = true;
			if (arg->hasValue()) EXCEPTION("--genetrees doesn't need a value");
		}
	}


	// Check if user wants to limit species tree to leaf set of species tree when computing losses. Corresponds to version 1 of DupLoss 
      
        {
                const Argument *arg = Argument::find("--limit");
                if (arg != NULL) {
                        LIMIT_LOSSES = true;
                        if (arg->hasValue()) EXCEPTION("--limit doesn't need a value");
                }
        }


	// choose tree generator
	int generator;
	{
		const Argument *arg1 = Argument::find("-g");
		const Argument *arg2 = Argument::find("--generator");
		const Argument *arg = arg1 != NULL ? arg1 : arg2;
		if (arg == NULL) generator = 1;
		else arg->convert(generator);
	}

	// puts the leaf add into fast mode
	bool fastFlag = true;

	// constraints input
	istream *constraints_in;
	ifstream constraints_infile;
	{
		const Argument *arg = Argument::find("--constraints");
		if (arg == NULL) constraints_in = NULL;
		else {
			string filename;
			arg->convert(filename);
			msgout << "Constraints file: " << filename << endl;
			constraints_infile.open(filename.c_str());
			if (!constraints_infile) EXCEPTION("constraints file " << filename << " not found");
			constraints_in = &constraints_infile;
			if ((generator != 1) && (generator != 2)) WARNING("--constraints does not apply to generator " << generator);
		}
	}

	#include "commonarguments.cc"

	// output help message
	if (helpFlag) {
		cout << "Usage: " << argsv[0] << " [ARGUMENT]" << endl;
		cout << endl;
		cout << "  -i, --input                   Input file" << endl;
		cout << "  -o, --output                  Output file" << endl;
 		cout << "      --genetrees               Output input gene trees and their reconciliation costs after the speciestree" << endl;
		cout << "  -g, --generator 0|1           Options for species tree to use to initialize the local search heuristic" << endl;
		cout << "                                0 - use user-provided species tree as initial species tree" << endl;
		cout << "                                1 - use leaf adding heuristic to build initial species tree [default]" << endl;
		cout << "      --constraints <file>      A file containing groupings of species for generator 1." << endl;
		cout << "      --limit                   Limit species tree to leaf set of gene tree when computing losses. Not recommended." << endl;
		cout << "  -q, --quiet                   No processing output." << endl;
		cout << "      --seed <integer number>   Set a user defined random number generator seed." << endl;
		cout << "  -v, --version                 Output the version number." << endl;
		cout << "  -h, --help                    Output brief help message and example commands." << endl;
		cout << endl;
		cout << "Examples:" << endl;
		cout << "  " << argsv[0] << " -i trees.newick -o speciestree.newick" << endl;
		cout << "  cat speciestree.newick genetrees.newick | " << argsv[0] << " -g 0" << endl;
		return 0;
	}


	// random seed
	msgout << "Random seed: " << randomseed << endl;

	// final gene duplications
// 	if (score_flag) msgout << "Printing final gene duplications with resulting species tree(s)" << endl;




	// -------------------------------------------------------------------------------------------
	// run tree generator
	ostringstream ostreamBuildtree;
	switch (generator) {
		case 0: {
				msgout << "Using user-provided initial species tree" << endl;
		} break;
		case 1: {
			// copy input into a string
			ostringstream ostreamTemp;
			input.skipComments = false;
			char ch;
			while (input.nextChar(ch)) ostreamTemp << ch;
			input.skipComments = true;
			istringstream istreamBuildtree;
			istreamBuildtree.str(ostreamTemp.str());
			Input inputBuildtree(&istreamBuildtree);

			buildtree::HeuristicLeafAdd *heuristic = NULL;
//			if(fastFlag)
//			{
				msgout << "Using fast randomized leaf adding heuristic to build initial species tree" << endl;
				heuristic = new buildtree::HeuristicSecondFast();

				// read input trees
				heuristic->readTrees(inputBuildtree);
//			}


			// read constraints
			if (constraints_in!=NULL) {
				Input constraints_input(constraints_in);
				heuristic->readConstraints(constraints_input);
			}

			ReRoot reroot=OPT;
			heuristic->run(ostreamBuildtree, reroot);

			// output genetrees
			ostreamBuildtree << ostreamTemp.str();
			delete heuristic;
		} break;

		default: {
			EXCEPTION("Unknown initial species tree generator " << generator);
		} break;
	}
	// create input for next step - search step
	istringstream istreamSearch;
	istreamSearch.str(ostreamBuildtree.str());
	Input inputSearch2(&istreamSearch);
	Input &inputSearch = (generator==0) ? input : inputSearch2;

	// -------------------------------------------------------------------------------------------
	// run the search heuristic
	ostringstream ostreamSearchSpeciesTree;
	ostringstream ostreamSearchGeneTrees;
	gtpspr::Heuristic *heuristic = NULL;


// Commenting out different heuristic cases below since only going to use simple SPR search heuristic. 
//	switch (heuristic_type) {
//		case 1: {
			msgout << "Using SPR local search heuristic to infer final species tree" << endl;
			heuristic = new gtpspr::SimpleHeuristicRandom(oformat, score_flag);
//		} break;
//		case 2: {
//			msgout << "Using partial queue based heuristic" << endl;
//			heuristic = new gtpspr::AdvancedHeuristic(oformat, score_flag);
//		} break;
//		case 3: {
//			msgout << "Using full queue based heuristic" << endl;
//			heuristic = new gtpspr::BioHeuristic(queue, max_trees, oformat, score_flag);
//		} break;
//		default: {
//			EXCEPTION("Unknown heuristic " << heuristic);
//		} break;
//	}

	// read input trees
	heuristic->readTrees(inputSearch);

	// run timed
	time_t startTime = time(NULL);
	heuristic->run(ostreamSearchSpeciesTree, reroot);
	time_t endTime = time(NULL);

	
 	// output genetrees
	for (int i = 0; i < heuristic->genetree_rooted.size(); i++) {
		heuristic->genetree_rooted.at(i)->tree2newick(ostreamSearchGeneTrees); ostreamSearchGeneTrees << endl;
	}
	for (int i = 0; i < heuristic->genetree_unrooted.size(); i++) {
		ostreamSearchGeneTrees << "[&U]"; heuristic->genetree_unrooted.at(i)->tree2newick(ostreamSearchGeneTrees); ostreamSearchGeneTrees << endl;
	}

	delete heuristic;

	// timing
	{
		long unsigned int duration = endTime - startTime;
		const int days = duration / (60 * 60 * 24); duration -= days * (60 * 60 * 24);
		const int hours = duration / (60 * 60); duration -= hours * (60 * 60);
		const int minutes = duration / 60; duration -= minutes * 60;
		const int seconds = duration;
		ostringstream os;
		bool on = false;
		if (days!=0) {on = true; os << days << "d ";}
		if (on || (hours!=0)) {on = true; os << hours << "h ";}
		if (on || (minutes!=0)) {on = true; os << minutes << "m ";}
		os << seconds << 's';
		msgout <<"Time taken for search heuristic: " << os.str() << endl;
	}






	// -------------------------------------------------------------------------------------------
	// create a report
	{
		// header
		if (oformat == NEXUS) {
			*out << "#nexus" << endl;
		}
		// each set of species tree and gene trees
		istringstream MyStream(ostreamSearchSpeciesTree.str());
		int no=1;
		for (string tree;getline(MyStream, tree);no++) {
			if (no > 1) *out << endl;

			ostringstream ostreamReport;
			ostreamReport << tree << ostreamSearchGeneTrees.str();
			istringstream istreamReport(ostreamReport.str());
			Input inputReport(&istreamReport);

			genedupreport::Report report;

			// read input trees
			report.readTrees(inputReport);

			// create record
			report.createReportInComment(*out, true, no, genetreesFlag, score_flag, oformat);
		}
		// footer
		if (oformat == NEXUS) {
		}
	}


	#ifdef DEBUG
	for (map<string,int>::iterator itr=NodeID::namelist.begin();itr!=NodeID::namelist.end();itr++)
		cout << "Species namelist entry:" << itr->first << ' ' << itr->second << endl;
	msgout << "\nAll done!" << endl;
	#endif

	return 0;
}
