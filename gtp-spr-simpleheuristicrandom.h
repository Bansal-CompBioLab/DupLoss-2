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

#ifndef GTP_SPR_SIMPLEHEURISTICRANDOM_H
#define GTP_SPR_SIMPLEHEURISTICRANDOM_H

class SimpleHeuristicRandom :  public Heuristic {
protected:
        double Best_score;
	bool update;
	bool score_flag;
	struct data
	{
		SpeciesNode *BestSubtreeRoot;
		SpeciesNode *BestNewLocation;
	} temp, old;
	vector<data> queue;
	int j, index;
	Format format;
	unsigned int countTotal;

	
	
public:
        inline void scoreComputed(SpeciesNode &node)
	{
                const double &genedup = node.score;
				const double &geneloss = node.lossScore;
                if(genedup + geneloss < Best_score)
                {
			update = true;
			queue.clear();
			temp.BestSubtreeRoot = speciestree->nodes[j];
                        temp.BestNewLocation = &node;
			queue.push_back(temp);
                        Best_score = genedup + geneloss;

                        msgout<< "\rCurrent best score: " << genedup  << " + " << geneloss << " = " <<Best_score<<"         ";
			flush(cout);
                }
		else if ((update== true) && ((genedup + geneloss) == Best_score))
		{
			temp.BestSubtreeRoot = speciestree->nodes[j];
                        temp.BestNewLocation = &node;
                        queue.push_back(temp);
		}

        }

        SimpleHeuristicRandom(Format format, bool score_flag) : format(format), score_flag(score_flag)
        {
		index=0;
		//srand(time(0));
		Best_score = UINT_MAX;
		update = true;
		temp.BestSubtreeRoot = NULL;
		temp.BestNewLocation = NULL;
		countTotal = 0;
		queue.push_back(temp);
	}

	void run(ostream &os, const ReRoot reroot)
	{
		bool rerooting = (reroot == ALL);
                createLeafMapping();

                int num_nodes, Left_Right;
                num_nodes = speciestree->nodes.size();

//		list<data>:: iterator iter;

		msgout << "Computing...\n";
		do {

                        for ( j=0; j< num_nodes;j++)
                        {

                                SpeciesNode *prnt, *sblng;

                                if (speciestree->nodes[j] == speciestree->root)
                                        continue;

                                prnt = speciestree->nodes[j]->parent();
                                if (prnt->child(0) == speciestree->nodes[j])
                                        Left_Right = 0;
                                else Left_Right =1;

                                sblng = prnt->child(1-Left_Right);

// checkConstraintsStructure(speciestree->root,0);
                                speciestree->moveSubtree(speciestree->nodes[j], speciestree->root);

                                computeGeneDuplications(speciestree->nodes[j], rerooting);
                                speciestree->moveSubtree(speciestree->root->child(Left_Right), sblng);
// checkConstraintsStructure(speciestree->root,0);

                        }
						


			if (update == false) {

				computeBestRooting();
				if (rerooting) break;
				rerooting = true;
	
				queue.clear();
			} else {
//				computeBestRooting();
				rerooting = (reroot == ALL);
			}

			update = false;
		
			
			
			if (!queue.empty()) {
				countTotal++;
				index = rand()%queue.size();
				old.BestSubtreeRoot = queue[index].BestSubtreeRoot;
				old.BestNewLocation = queue[index].BestNewLocation;
				speciestree->moveSubtree(old.BestSubtreeRoot, old.BestNewLocation);
				computeBestRooting();
				
			}
		} while(!interruptFlag);

		msgout << endl;
		msgout << "Number of rSPR tree edit operations: " << countTotal << endl;

		// output score
		msgout << "Final weighted reconciliation cost: " << getCurrentScore() << endl;
		WEIGHTED_RECON_COST = getCurrentScore();
		speciestree->tree2newick(os); os << endl;


	}


	double getCurrentScore() {
		return Best_score;
	}
};

#endif
