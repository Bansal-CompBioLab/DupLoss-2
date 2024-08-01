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


#ifndef REPORT_H
#define REPORT_H

// ------------------------------------------------------------------------------------------------------------------
class Report : public TreeSet {
public:
	Report() {
		#ifdef DEBUG
		cout << "Report created" << endl;
		#endif
	}

	virtual ~Report() {
		#ifdef DEBUG
		cout << "Report destroyed" << endl;
		#endif
	}

	// ------------------------------------------------------------------------------------------------------
	// compute the gene duplications
	inline void computeGeneDuplications(bool reroot = true) {
		speciestree->establishOrder();
		speciestree->preprocessLCA();
//		resetLossStuff(speciestree->root);
		// process all rooted trees
		for(int i=0; i<genetree_rooted.size(); i++) {
			GeneTreeRooted &tree = *genetree_rooted[i];
			createPrimaryMapping(tree);
			unsigned int &score = tree.score;
			score = getScore(tree);

			resetRelevantTree(speciestree->root);
			buildRelevantTree(tree);
			unsigned int &lossScore = tree.lossScore;
			lossScore  = computeGeneLossForRoot(tree);
			
		}
		// process all unrooted trees
		for(int i=0; i<genetree_unrooted.size(); i++) {
			GeneTreeUnrooted &tree = *genetree_unrooted[i];
			if (reroot) { // find the best geneduplication score of all rootings (rerooting of the genetrees)
				createPrimaryMappingUnrooted(tree);
				unsigned int &score = tree.score;
				score = getScore(tree);
				resetRelevantTree(speciestree->root);
				buildRelevantTree(tree);
				unsigned int &lossScore = tree.lossScore;
				lossScore  = computeGeneLossForRoot(tree);
		
				GeneNodeUnrooted *u = tree.root->child(0);
				GeneNodeUnrooted *v = tree.root->child(1);
				best_node[0] = u;
				best_node[1] = v;
				#ifdef DEBUG
				position_counter_debug = 2;
				if ((u == NULL) || (v == NULL)) EXCEPTION("child = NULL in computeGeneDuplications" << endl);
				#endif
				moveRoot(tree, u, score, lossScore);
				moveRoot(tree, v, score, lossScore);
				#ifdef DEBUG
				if (position_counter_debug != tree.nodes.size()-1)
					WARNING("tree traversal failed in computeGeneDuplications" << position_counter_debug << " != " << tree.nodes.size()-1 << endl);
				#endif
				tree.reroot(best_node[0], best_node[1]);
			} else { // find the best genedupication of the current rooting
				createPrimaryMapping(tree);
				unsigned int &score = tree.score;
				score = getScore(tree);
				
				resetRelevantTree(speciestree->root);
				buildRelevantTree(tree);
				unsigned int &lossScore = tree.lossScore;
				lossScore  = computeGeneLossForRoot(tree);
			}
		}
		speciestree->postprocessLCA();
	}

	// calculate the gene duplication for all rootings
	#ifdef DEBUG
	int position_counter_debug;
	#endif
	GeneNodeUnrooted *best_node[2];
	inline void moveRoot(GeneTreeUnrooted &tree, GeneNodeUnrooted *p, unsigned int &best_score1, unsigned int &best_score2) {
		GeneNodeUnrooted *c[2];
		for (int i=0; i<2; i++) c[i] = p->child(i);
		for (int i=0; i<2; i++) {
			if (c[i] != NULL) {
				moveRoot(tree, c[i], best_score1, best_score2);
				tree.reroot(c[i], p);
				#ifdef DEBUG
				tree.checkStructure();
				position_counter_debug++;
				#endif
				const int score = getScore(tree);
				resetRelevantTree(speciestree->root);
				buildRelevantTree(tree);
				const int lossScore  = computeGeneLossForRoot(tree);
				
				if ((score + lossScore) < (best_score1 + best_score2)) {
					best_score1 = score; 
					best_score2 = lossScore;
					best_node[1] = c[i];
					best_node[0] = p;
				}
			}
		}
	}
	
	// Reset the subtreeSize counter at each node to 0, reset the isRelevant flag to false and and set nodeDepth to 0
	void resetRelevantTree(SpeciesNode *node){
		if (node->child(0) != NULL) {
			resetRelevantTree(node->child(0));
		}
		if (node->child(1) != NULL) {
			resetRelevantTree(node->child(1));
		}

		node->subtreeSize = 0;
//		node->lossScoreTemp = 0;
		node->isRelevant = false;
		node->nodeDepth = 0;
	}



	// Mark the useless nodes that are not present in the species tree build on the gene tree's taxa set. 
	void buildRelevantTree(GeneTreeRooted &tree) {
		vector<NamedGeneNodeRooted*> &g = tree.leafnodes;
		// process leaf mappings
		for(vector<NamedGeneNodeRooted*>::iterator itr=g.begin(); itr!=g.end(); itr++) {
			NamedGeneNodeRooted *genenode = *itr;
			SpeciesNode *&mapping = genenode->getMapping();
			mapping->subtreeSize = mapping->subtreeSize + 1;
		}

		// perform post order traversal, to compute the value of subtreeSize at each species node.
		doPostOrder(speciestree->root);

		// Set the isRelevant flag for each node 
		doSecondPostOrder(speciestree->root);
		
		// set the nodeDepth flag to the correct value
		doPreOrder(speciestree->root, 1);



	}
   
	
	void buildRelevantTree(GeneTreeUnrooted &tree) {
		vector<NamedGeneNodeUnrooted*> &g = tree.leafnodes;
		// process leaf mappings
		for(vector<NamedGeneNodeUnrooted*>::iterator itr=g.begin(); itr!=g.end(); itr++) {
			NamedGeneNodeUnrooted *genenode = *itr;
			SpeciesNode *&mapping = genenode->getMapping();
			mapping->subtreeSize = mapping->subtreeSize + 1;
		}

		// perform post order traversal, to compute the value of subtreeSize at each species node.
		doPostOrder(speciestree->root);

		// Set the isRelevant flag for each node 
		doSecondPostOrder(speciestree->root);
		
		// set the nodeDepth flag to the correct value
		doPreOrder(speciestree->root, 1);


	}	
	
	
	// builds the subtreeSize counter by doing a post order traversal of the species tree 
	void doPostOrder(SpeciesNode *node) {
		if (node->child(0) != NULL) {
			doPostOrder(node->child(0));
		}
		if (node->child(1) != NULL) {
			doPostOrder(node->child(1));
		}

		if((node->child(0) != NULL) && (node->child(1) != NULL))
			node->subtreeSize = node->child(0)->subtreeSize + node->child(1)->subtreeSize;
	}

	// sets the isRelevant flag at each node by doing a post order traversal of the species tree
	void doSecondPostOrder(SpeciesNode *node) {
		if (node->child(0) != NULL) {
			doSecondPostOrder(node->child(0));
		}
		if (node->child(1) != NULL) {
			doSecondPostOrder(node->child(1));
		}

		if((node->child(0) != NULL) && (node->child(1) != NULL))
		{
			if((node->child(0)->subtreeSize != 0) && (node->child(1)->subtreeSize != 0))
				node->isRelevant = true;
			else 
			{
                                if (LIMIT_LOSSES == false)
                                        node->isRelevant = true;
                                else
                                        node->isRelevant = false;
                        }


		}
		else {
			if(node->subtreeSize == 0)
			{
                                if (LIMIT_LOSSES == false)
                                        node->isRelevant = true;
                                else
                                        node->isRelevant = false;
                        }

			else node->isRelevant = true;
		}
	}
 
	// do a preOrder traversal of the species tree to label nodes with their depth (used for computing loss values)
	void doPreOrder(SpeciesNode *node, int depthCount) {
		if (node->isRelevant == true)
		{
			node->nodeDepth = depthCount;
			depthCount++;
		}
		
		if (node->child(0) != NULL) {
			doPreOrder(node->child(0), depthCount);
		}
		if (node->child(1) != NULL) {
			doPreOrder(node->child(1), depthCount);
		}
	}


	
	// compute the loss value for the current species tree (which corresponds to "regrafting at the root")
	unsigned int computeGeneLossForRoot(GeneTreeRooted &tree) {
		
		unsigned int c=0;
		vector<GeneNodeRooted*> &g = tree.nodes;
		for(vector<GeneNodeRooted*>::iterator itr=g.begin(); itr!=g.end(); itr++) {
			GeneNodeRooted *genenode = *itr;
			if (!genenode->isLeaf())
			{
				SpeciesNode *&mapping = genenode->getMapping();
				SpeciesNode *&leftchild_mapping = genenode->child(0)->getMapping();
				SpeciesNode *&rightchild_mapping = genenode->child(1)->getMapping();
			
				// start a case by case evaluation of the loss value for the "mapping" node
				if((mapping == leftchild_mapping) && (mapping == rightchild_mapping))
				{
				// do nothing
				}
			
				else if (mapping == leftchild_mapping)
				{
					c = c + rightchild_mapping->nodeDepth -  mapping->nodeDepth;
				
				}
			
				else if (mapping == rightchild_mapping)
				{
			
					c = c + leftchild_mapping->nodeDepth -  mapping->nodeDepth;
				}
			
				else 
				{
					c = c + leftchild_mapping->nodeDepth + rightchild_mapping->nodeDepth -  (2 * mapping->nodeDepth) -2;
				
				}
			}
		}
		return c;
//		node->lossScore = c; 
	}
	
	
	unsigned int computeGeneLossForRoot(GeneTreeUnrooted &tree) {
		
		unsigned int c=0;
		vector<GeneNodeUnrooted*> &g = tree.nodes;
		for(vector<GeneNodeUnrooted*>::iterator itr=g.begin(); itr!=g.end(); itr++) {
			GeneNodeUnrooted *genenode = *itr;
			if (!genenode->isLeaf())
			{
				SpeciesNode *&mapping = genenode->getMapping();
				SpeciesNode *&leftchild_mapping = genenode->child(0)->getMapping();
				SpeciesNode *&rightchild_mapping = genenode->child(1)->getMapping();
			
				// start a case by case evaluation of the loss value for the "mapping" node
				if((mapping == leftchild_mapping) && (mapping == rightchild_mapping))
				{
				// do nothing
				}
			
				else if (mapping == leftchild_mapping)
				{
					c = c + rightchild_mapping->nodeDepth -  mapping->nodeDepth;
				
				}
			
				else if (mapping == rightchild_mapping)
				{
			
					c = c + leftchild_mapping->nodeDepth -  mapping->nodeDepth;
				}
			
				else 
				{
					c = c + leftchild_mapping->nodeDepth + rightchild_mapping->nodeDepth -  (2 * mapping->nodeDepth) -2;
				
				}
			}
		}
		return c;
//		node->lossScore = c; 
	}
	
	// ------------------------------------------------------------------------------------------------------
	// output a report of the trees
	void createReport(ostream &os, const bool rerooting) {
		createLeafMapping();
		computeGeneDuplications(rerooting);
		msgout << endl;
		// species tree
		{
			os << "Species tree" << endl;
			os << "------------" << endl;
			os << "Tree: "; speciestree->tree2newick(os); os << endl;
			os << "Size: " << speciestree->nodes.size() << " nodes " << speciestree->leafnodes.size() << " leaves " << endl;
			os << endl;
		}

		// output gene trees
		{
			os << "Gene tree" << endl;
			os << "---------" << endl;
			for (int i = 0; i < genetree_rooted.size(); i++) {
				os << "Tree: "; genetree_rooted.at(i)->tree2newick(os); os << endl;
				os << "Size: " << genetree_rooted.at(i)->nodes.size() << " nodes " << genetree_rooted.at(i)->leafnodes.size() << " leaves" << endl;
				os << "Gene duplications + Losses: " << genetree_rooted.at(i)->score << " + " << genetree_rooted.at(i)->lossScore << " = " << genetree_rooted.at(i)->score + genetree_rooted.at(i)->lossScore << endl;
				os << endl;
			}
			for (int i = 0; i < genetree_unrooted.size(); i++) {
				os << "Tree: "; os << "[&U]"; genetree_unrooted.at(i)->tree2newick(os); os << endl;
				os << "Size: " << genetree_unrooted.at(i)->nodes.size() << " nodes " << genetree_unrooted.at(i)->leafnodes.size() << " leaves" << endl;
				os << "Gene duplications + Losses: " << genetree_unrooted.at(i)->score << " + " << genetree_unrooted.at(i)->lossScore << " = " << genetree_unrooted.at(i)->score + genetree_unrooted.at(i)->lossScore << endl;
				os << endl;
			}
		}

		// total score
		{
			unsigned int total_score = 0;
			for (int i = 0; i < genetree_rooted.size(); i++) {
				total_score += genetree_rooted.at(i)->score + genetree_rooted.at(i)->lossScore;
			}
			for (int i = 0; i < genetree_unrooted.size(); i++) {
				total_score += genetree_unrooted.at(i)->score + genetree_unrooted.at(i)->lossScore;
			}
			os << "Total Reconciliation Cost: " << total_score << endl;
			os << endl;
		}
	}

	// ------------------------------------------------------------------------------------------------------
	// output a report of the trees
	void createReportInComment(ostream &os, const bool rerooting, const int no, const bool genetreesFlag, const bool score_flag, const Format format) {
		createLeafMapping();
		computeGeneDuplications(rerooting);

		// species tree
		{
			if (format == NEXUS) os << "begin trees;" << endl;
//			os << "[Species tree " << no << "]" << endl;
                        os << "[Species tree]" << endl;
		
			os << "[Size: " << speciestree->nodes.size() << " nodes " << speciestree->leafnodes.size() << " leaves]" << endl;
			// total score
			if (score_flag) {
				unsigned int total_scoreGD = 0;
				unsigned int total_scoreLosses = 0;
				for (int i = 0; i < genetree_rooted.size(); i++) {
					total_scoreGD += genetree_rooted.at(i)->score;
					total_scoreLosses += genetree_rooted.at(i)->lossScore;

				}
				for (int i = 0; i < genetree_unrooted.size(); i++) {
					total_scoreGD += genetree_unrooted.at(i)->score;
					total_scoreLosses += genetree_unrooted.at(i)->lossScore;

				}
				os << "[Weighted reconciliation cost: " << WEIGHTED_RECON_COST << "]" << endl;
				os << "[Total number of gene duplications + losses: " << total_scoreGD << " + " << total_scoreLosses << " = " << total_scoreGD + total_scoreLosses << "]" << endl;
			}
			// species tree
			if (format == NEWICK) {
				speciestree->tree2newick(os); os << endl;
			} else
			if (format == NEXUS) {
				os << "tree speciestree" << no << " = "; speciestree->tree2newick(os); os << endl;
			}
			if (format == NEXUS) os << "end;" << endl;
		}

		// gene trees
		if (genetreesFlag) {
			os << endl;
			os << "[Gene trees]" << endl;
			if (format == NEXUS) os << "begin trees;" << endl;
			int count = 1;
			// rooted gene trees
			for (int i = 0; i < genetree_rooted.size(); i++, count++) {
				os << "[Rooted gene tree " << count << "]" << endl;
				os << "[Size: " << genetree_rooted.at(i)->nodes.size() << " nodes " << genetree_rooted.at(i)->leafnodes.size() << " leaves]" << endl;
				if (score_flag) {
					os << "[Gene duplications + Losses: " << genetree_rooted.at(i)->score << " + " << genetree_rooted.at(i)->lossScore << " = " << genetree_rooted.at(i)->score + genetree_rooted.at(i)->lossScore << "]" << endl;
				}
				if (format == NEWICK) {
					genetree_rooted.at(i)->tree2newick(os); os << endl;
				} else
				if (format == NEXUS) {
					os << "tree genetree" << count << " = "; genetree_rooted.at(i)->tree2newick(os); os << endl;
				}
				if (i != genetree_rooted.size()) os << endl;
			}
			// unrooted gene trees
			for (int i = 0; i < genetree_unrooted.size(); i++, count++) {
				os << "[Rooted gene tree " << count;
				if (rerooting) os << " (unrooted in input)";
				os << "]" << endl;
				os << "[Size: " << genetree_unrooted.at(i)->nodes.size() << " nodes " << genetree_unrooted.at(i)->leafnodes.size() << " leaves]" << endl;
				if (score_flag) {
					os << "[Gene duplications + Losses: " << genetree_unrooted.at(i)->score << " + " << genetree_unrooted.at(i)->lossScore << " = " << genetree_unrooted.at(i)->score + genetree_unrooted.at(i)->lossScore << "]" << endl;
				}
				if (format == NEWICK) {
					genetree_unrooted.at(i)->tree2newick(os); os << endl;
				} else
				if (format == NEXUS) {
					os << "tree genetree" << count << " = "; genetree_unrooted.at(i)->tree2newick(os); os << endl;
				}
				if (i != genetree_unrooted.size()-1) os << endl;
			}
			if (format == NEXUS) os << "end;" << endl;
		}

	}

	// ------------------------------------------------------------------------------------------------------
	// outputs the trees into a string stream
	friend ostream & operator << (ostream &os, Report &t);
};

// outputs the trees into a string stream
ostream & operator << (ostream &os, Report &t) {
	os << (TreeSet&) t;
	return os;
}

#endif
