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

#ifndef GTP_SPR_HEURISTIC_H
#define GTP_SPR_HEURISTIC_H

// ------------------------------------------------------------------------------------------------------------------
class Heuristic : public TreeSet {
public:
	using TreeSet::getScore;

	Heuristic() {
		#ifdef DEBUG
		cout << "Heuristic created" << endl;
		#endif
	}

	virtual ~Heuristic() {
		#ifdef DEBUG
		cout << "Heuristic destroyed" << endl;
		#endif
	}

	// ------------------------------------------------------------------------------------------------------
	// compute the gene duplications
	// subtree is the subtree pruned (P)

	inline void computeGeneDuplications(SpeciesNode *subtree, bool reroot = true) {
		SpeciesNode *sibling = subtree->getSibling();
		SpeciesNode *LossSibling;
		speciestree->establishOrder();
		speciestree->preprocessLCA();
		resetGeneDuplications(sibling);
		
		//resetting stuff related to losses;
		resetLossStuff(speciestree->root);

		// process all rooted trees

//cout << "next iteration....................." << endl;
		for(int i=0; i<genetree_rooted.size(); i++) {
			GeneTreeRooted &tree = *genetree_rooted[i];
			createPrimaryMapping(tree);
			const double score = getScore(tree) * tree.weight;
			createSecondaryMapping(tree, subtree);
			computeGeneDuplicationsTriple();
//			removeSecondaryMapping(tree);
			computeGeneDuplicationsAdd(sibling, score, tree.weight);

				// Add stuff to handle losses
			resetRelevantTree(speciestree->root);
			buildRelevantTree(tree);
			double initLoss  = tree.weight * computeGeneLossForRoot(tree);

			copyInitialLossScoreToAllNodes(sibling, initLoss);

			if(speciestree->root->isRelevant == false)
			{
//				cout << "Am I here by any chance?" << endl;
				// If the gene tree is fully contained in either the pruned subtree or in the other subtree then there is nothing to be done
			}
			else
			{
				computeGeneLossCounters(tree, subtree);
				
				if(speciestree->root->child(0) == subtree)
					LossSibling = speciestree->root->LossChild2;
				else LossSibling = speciestree->root->LossChild1;
				computeLossScores(LossSibling, tree.weight);				
			}
			
			
			removeSecondaryMapping(tree);
			
		}



		// process all unrooted trees
		for(int i=0; i<genetree_unrooted.size(); i++) {
			GeneTreeUnrooted &tree = *genetree_unrooted[i];
			if (reroot) { // find the best geneduplication score of all rootings (rerooting of the genetrees)
				createPrimaryMappingUnrooted(tree);
				double &score = best_score;
				score = getScore(tree) * tree.weight;
				createSecondaryMapping(tree, subtree);
				computeGeneDuplicationsTriple();
//				removeSecondaryMapping(tree);
				computeGeneDuplicationsTempReplace(sibling, score, tree.weight);
				
				
				
				
				resetRelevantTree(speciestree->root);
				buildRelevantTree(tree);
				double initLoss  = tree.weight * computeGeneLossForRoot(tree);

				copyInitialLossScoreToAllNodes_Temp(sibling, initLoss);

				if(speciestree->root->isRelevant == false)
				{
					// If the gene tree is fully contained in either the pruned subtree or in the other subtree then there is nothing to be done
				}
				else
				{
					computeGeneLossCounters(tree, subtree);
			
					if(speciestree->root->child(0) == subtree)
						LossSibling = speciestree->root->LossChild2;
					else LossSibling = speciestree->root->LossChild1;
					computeLossScores_Temp(LossSibling, tree.weight);				
				}
				removeSecondaryMapping(tree);
				
				
				
				
				GeneNodeUnrooted *u = tree.root->child(0);
				GeneNodeUnrooted *v = tree.root->child(1);
				best_node[0] = u;
				best_node[1] = v;
				#ifdef DEBUG
				position_counter_debug = 2;
				if ((u == NULL) || (v == NULL)) EXCEPTION("child = NULL in computeGeneDuplications" << endl);
				#endif
				moveRoot(tree, u, subtree, sibling);
				moveRoot(tree, v, subtree, sibling);
				#ifdef DEBUG
				if (position_counter_debug != tree.nodes.size()-1)
					WARNING("tree traversal failed in computeGeneDuplications" << position_counter_debug << " != " << tree.nodes.size()-1 << endl);
				#endif
				tree.reroot(u, v);
				addTempGeneDuplications(sibling);
			} 
			else 
			{ // find the best genedupication of the current rooting
				createPrimaryMapping(tree);
				const double score = getScore(tree) * tree.weight;
				createSecondaryMapping(tree, subtree);
				computeGeneDuplicationsTriple();
//				removeSecondaryMapping(tree);
				computeGeneDuplicationsAdd(sibling, score, tree.weight);
				
					// Add stuff to handle losses
				resetRelevantTree(speciestree->root);
				buildRelevantTree(tree);
				double initLoss  = tree.weight * computeGeneLossForRoot(tree);

				copyInitialLossScoreToAllNodes(sibling, initLoss);

				if(speciestree->root->isRelevant == false)
				{
					// If the gene tree is fully contained in either the pruned subtree or in the other subtree then there is nothing to be done
				}
				else
				{
					computeGeneLossCounters(tree, subtree);
			
					if(speciestree->root->child(0) == subtree)
						LossSibling = speciestree->root->LossChild2;
					else LossSibling = speciestree->root->LossChild1;
					computeLossScores(LossSibling, tree.weight);				
				}
				removeSecondaryMapping(tree);
				
			}
		}
		speciestree->postprocessLCA();
		
		// 1 line of code to assign loss scores for non-relevant nodes.
//		AssignLossScoresToNonRelevantNodes(sibling);		

/*		// extra 1 line of code for losses
//		AddLossCounterToDupScoreCounter(sibling);
*/		
//		outputSpeciesTreeNodes(speciestree->root);	
//		cout << endl<< endl;
		forEachCallScoreComputed(subtree, sibling);
		
	}


//-------------------Losses stuff begin-------------------------------------------------------------------

/* Function used for testing purposes */	
	void outputSpeciesTreeNodes(SpeciesNode *node)
	{
		if (node->child(0) != NULL) {
			outputSpeciesTreeNodes(node->child(0));
		}
		if (node->child(1) != NULL) {
			outputSpeciesTreeNodes(node->child(1));
		}

//		if (node->child(0) == NULL)
//			cout << (NodeID)node << " ";
		cout << *node;
	}
	

	void resetLossStuff(SpeciesNode *&node) {
		if (node == NULL) return;
		node->lossScore = 0;
		node->lossScoreTemp = 0;
		for (int i=0; i<2; i++)
			resetLossStuff(node->child(i));
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
		node->LossScoreDiff = 0;
//		node->lossScoreTemp = 0;
		node->isRelevant = false;
		node->nodeDepth = 0;
		node->lossCounter1 = 0;
		node->lossCounter2 = 0;
		node->lossCounter3 = 0;
		node->lossCounter4 = 0;
		node->lossCounter5 = 0;
		node->lossCounter6 = 0;
		node->LossParent = NULL;
		node->LossChild1 = NULL;
		node->LossChild2 = NULL;

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

		// Set the LossParent, LossChild[0], and LossChild[1] pointers
		setPointersPostOrder(speciestree->root);

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

		// Set the LossParent, LossChild[0], and LossChild[1] pointers
		setPointersPostOrder(speciestree->root);

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
				if(LIMIT_LOSSES == false)
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

	void setPointersPostOrder(SpeciesNode *node) {
		if (node->child(0) != NULL) {
			setPointersPostOrder(node->child(0));
		}
		if (node->child(1) != NULL) {
			setPointersPostOrder(node->child(1));
		}

		if(node->isRelevant == true)
		{
			SpeciesNode* Temp_Node = node;
			while(Temp_Node!= speciestree->root)
			{
				Temp_Node = Temp_Node->parent();
				if (Temp_Node->isRelevant == true)
				{
					if(Temp_Node->LossChild1 == NULL)
						Temp_Node->LossChild1 = node;
					else Temp_Node->LossChild2 = node;
					node->LossParent = Temp_Node;
					break;
				}
			}
		}
	}
	
	// compute the loss value for the current species tree (which corresponds to "regrafting at the root")
	int computeGeneLossForRoot(GeneTreeRooted &tree) {
		
		int c=0;
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
	
	
	int computeGeneLossForRoot(GeneTreeUnrooted &tree) {
		
		int c=0;
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
	

	// Add the loss score computed at the root throughout the tree, to add and subtract to these values later.
	void copyInitialLossScoreToAllNodes(SpeciesNode *&node, double initLoss)
	{
		if (node == NULL) return;
		node->lossScore = node->lossScore + initLoss;
		for (int i=0; i<2; i++)
			copyInitialLossScoreToAllNodes(node->child(i), initLoss);
	}


	void computeGeneLossCounters(GeneTreeRooted &tree,  SpeciesNode *subtree)
	{
		vector<GeneNodeRooted*> &g = tree.nodes;
		SpeciesNode *temp_mapping, *sibli;

		

		for(vector<GeneNodeRooted*>::iterator itr=g.begin(); itr!=g.end(); itr++) {
			GeneNodeRooted *genenode = *itr;
			if (!genenode->isLeaf())
			{
				SpeciesNode *&mapping = genenode->getMapping();
				SpeciesNode *&leftchild_mapping = genenode->child(0)->getMapping();
				SpeciesNode *&rightchild_mapping = genenode->child(1)->getMapping();
			
				// start a case by case analysis for the different possible cases
				
				// Case 1: g does not map to the root node
				if(mapping != speciestree->root)
				{
					
					temp_mapping = genenode->child(0)->getMapping();
					while(temp_mapping != mapping)
					{
						temp_mapping->lossCounter1 = temp_mapping->lossCounter1 + 1;
						temp_mapping = temp_mapping->LossParent;
					}
					
					temp_mapping = genenode->child(1)->getMapping();
					while(temp_mapping != mapping)
					{
						temp_mapping->lossCounter1 = temp_mapping->lossCounter1 + 1;
						temp_mapping = temp_mapping->LossParent;
					}
				}

				// Case 2: g, g' and g'' all map to the root
				else if((mapping == speciestree->root) && (rightchild_mapping == speciestree->root) && (leftchild_mapping == speciestree->root))
				{
					SpeciesNode *&secmapping = genenode->secondarymapping;
					SpeciesNode *&leftchild_secmapping = genenode->child(0)->secondarymapping;
					SpeciesNode *&rightchild_secmapping = genenode->child(1)->secondarymapping;
					
					if((secmapping!= leftchild_secmapping) && (secmapping!= rightchild_secmapping))
					{
						temp_mapping = genenode->child(0)->secondarymapping;
						while(temp_mapping != secmapping)
						{
							temp_mapping->lossCounter2 = temp_mapping->lossCounter2 + 1;
//							if(temp_mapping->parent() == secmapping)
//								temp_mapping->lossCounter2 = temp_mapping->lossCounter2 + 1;
							temp_mapping = temp_mapping->LossParent;
						}
			
						temp_mapping = genenode->child(1)->secondarymapping;
						while(temp_mapping != secmapping)
						{
							temp_mapping->lossCounter2 = temp_mapping->lossCounter2 + 1;
//							if(temp_mapping->parent() == secmapping)
//								temp_mapping->lossCounter2 = temp_mapping->lossCounter2 + 1;
							temp_mapping = temp_mapping->LossParent;
						}
					}
					
					else if((secmapping == leftchild_secmapping) && (secmapping!= rightchild_secmapping))
					{
						
						temp_mapping = genenode->child(1)->secondarymapping;
						while(temp_mapping != secmapping)
						{
							temp_mapping->lossCounter2 = temp_mapping->lossCounter2 + 1;
//							if(temp_mapping->parent() == secmapping)
//								temp_mapping->lossCounter2 = temp_mapping->lossCounter2 + 1;
							temp_mapping = temp_mapping->LossParent;
						}
					}
					
					
					else if((secmapping != leftchild_secmapping) && (secmapping == rightchild_secmapping))
					{
						
						temp_mapping = genenode->child(0)->secondarymapping;
						while(temp_mapping != secmapping)
						{
							temp_mapping->lossCounter2 = temp_mapping->lossCounter2 + 1;
//							if(temp_mapping->parent() == secmapping)
//								temp_mapping->lossCounter2 = temp_mapping->lossCounter2 + 1;
							temp_mapping = temp_mapping->LossParent;
						}
					}
					
				}
				
				// Case 3: g = root, g' = root, g'' = supporting
				
				else if((mapping == speciestree->root) && (leftchild_mapping == speciestree->root)  && ((rightchild_mapping->no < subtree->begin) || (rightchild_mapping->no > subtree->end)))
				{
					SpeciesNode *&secmapping = genenode->secondarymapping;
					SpeciesNode *&leftchild_secmapping = genenode->child(0)->secondarymapping;
					SpeciesNode *&rightchild_secmapping = genenode->child(1)->getMapping();

					if((secmapping!= leftchild_secmapping) && (secmapping!= rightchild_secmapping))
					{
						temp_mapping = genenode->secondarymapping;
						while(temp_mapping->LossParent != speciestree->root)
						{
							temp_mapping->lossCounter6 = temp_mapping->lossCounter6 + 1;
							//SpeciesNode *sib = temp_mapping->getSibling();
							temp_mapping->getLossSibling()->lossCounter6 = temp_mapping->getLossSibling()->lossCounter6 + 1;
							temp_mapping = temp_mapping->LossParent;
						}
//cout << secmapping->LossChild1 << " " << secmapping->child(0) << endl;
//cout << secmapping->LossChild2 << " " << secmapping->child(1) << endl;
//cout << secmapping->LossParent << " " << secmapping->parent() << endl;
    					secmapping->LossChild1->lossCounter6 = secmapping->LossChild1->lossCounter6 + 1;
						secmapping->LossChild2->lossCounter6 = secmapping->LossChild2->lossCounter6 + 1;
						temp_mapping = genenode->child(0)->secondarymapping;
						while(temp_mapping->LossParent != secmapping)
						{
							temp_mapping->lossCounter2 = temp_mapping->lossCounter2 + 1;
//							temp_mapping->lossCounter1 = temp_mapping->lossCounter1 + 1;
							temp_mapping = temp_mapping->LossParent;
						}
						temp_mapping->lossCounter6 = temp_mapping->lossCounter6 + 1;
						temp_mapping = genenode->child(1)->getMapping();
						while(temp_mapping != secmapping)
						{
							temp_mapping->lossCounter1 = temp_mapping->lossCounter1 + 1;
							temp_mapping = temp_mapping->LossParent;
						}
					}
					
					else if((secmapping == leftchild_secmapping) && (secmapping != rightchild_secmapping))
					{
						
					
						temp_mapping = genenode->secondarymapping;
						while(temp_mapping->LossParent != speciestree->root)
						{
							temp_mapping->lossCounter6 = temp_mapping->lossCounter6 + 1;
							temp_mapping->getLossSibling()->lossCounter6 = temp_mapping->getLossSibling()->lossCounter6 + 1;

							temp_mapping = temp_mapping->LossParent;
							
						}
						
						secmapping->LossChild1->lossCounter6 = secmapping->LossChild1->lossCounter6 + 1;
						secmapping->LossChild2->lossCounter6 = secmapping->LossChild2->lossCounter6 + 1;
						temp_mapping = genenode->child(1)->getMapping();
						
						while(temp_mapping != secmapping)
						{
							temp_mapping->lossCounter1 = temp_mapping->lossCounter1 + 1;
							
							temp_mapping = temp_mapping->LossParent;
						}
					}					
					
					else if((secmapping != leftchild_secmapping) && (secmapping == rightchild_secmapping))
					{
						
						temp_mapping = genenode->secondarymapping;
						while(temp_mapping->LossParent != speciestree->root)
						{
							temp_mapping->lossCounter6 = temp_mapping->lossCounter6 + 1;
							temp_mapping->getLossSibling()->lossCounter6 = temp_mapping->getLossSibling()->lossCounter6 + 1;
							
							temp_mapping = temp_mapping->LossParent;
						}
						
						temp_mapping = genenode->child(0)->secondarymapping;
						sibli = temp_mapping->getLossSibling();
						while(temp_mapping->LossParent != secmapping)
						{
							temp_mapping->lossCounter2 = temp_mapping->lossCounter2 + 1;
							temp_mapping = temp_mapping->LossParent;
							sibli = temp_mapping->getLossSibling();
						}
						
						sibli->lossCounter6 = sibli->lossCounter6 + 1;
						
					}
					
					else if((secmapping == leftchild_secmapping) && (secmapping == rightchild_secmapping))
					{
						
						temp_mapping = genenode->secondarymapping;
						
						while(temp_mapping->LossParent != speciestree->root)
						{
							
							temp_mapping->lossCounter6 = temp_mapping->lossCounter6 + 1;
							temp_mapping->getLossSibling()->lossCounter6 = temp_mapping->getLossSibling()->lossCounter6 + 1;
							
							temp_mapping = temp_mapping->LossParent;
							
						}
					
						
						if(genenode->secondarymapping->LossChild1 != NULL)
						{
							genenode->secondarymapping->LossChild1->lossCounter6 = genenode->secondarymapping->LossChild1->lossCounter6 + 1;
							genenode->secondarymapping->LossChild2->lossCounter6 = genenode->secondarymapping->LossChild2->lossCounter6 + 1;
						}
						
					}
				}

				// Case 3: symmetric case when g' and g'' are switched
				else if((mapping == speciestree->root) && (rightchild_mapping == speciestree->root)  && ((leftchild_mapping->no < subtree->begin) || (leftchild_mapping->no > subtree->end)))
				{
					SpeciesNode *&secmapping = genenode->secondarymapping;
					SpeciesNode *&leftchild_secmapping = genenode->child(0)->getMapping();
					SpeciesNode *&rightchild_secmapping = genenode->child(1)->secondarymapping;

					if((secmapping!= leftchild_secmapping) && (secmapping!= rightchild_secmapping))
					{
						temp_mapping = genenode->secondarymapping;
						while(temp_mapping->LossParent != speciestree->root)
						{
							temp_mapping->lossCounter6 = temp_mapping->lossCounter6 + 1;
							temp_mapping->getLossSibling()->lossCounter6 = temp_mapping->getLossSibling()->lossCounter6 + 1;
							
							temp_mapping = temp_mapping->LossParent;
						}
						secmapping->LossChild1->lossCounter6 = secmapping->LossChild1->lossCounter6 + 1;
						secmapping->LossChild2->lossCounter6 = secmapping->LossChild2->lossCounter6 + 1;
						
						temp_mapping = genenode->child(1)->secondarymapping;
						while(temp_mapping->LossParent != secmapping)
						{
							temp_mapping->lossCounter2 = temp_mapping->lossCounter2 + 1;
//							temp_mapping->lossCounter1 = temp_mapping->lossCounter1 + 1;
							temp_mapping = temp_mapping->LossParent;
						}
						temp_mapping->lossCounter6 = temp_mapping->lossCounter6 + 1;

						temp_mapping = genenode->child(0)->getMapping();
						while(temp_mapping != secmapping)
						{
							temp_mapping->lossCounter1 = temp_mapping->lossCounter1 + 1;
							temp_mapping = temp_mapping->LossParent;
						}
						
					}
					
					else if((secmapping != leftchild_secmapping) && (secmapping == rightchild_secmapping))
					{
						temp_mapping = genenode->secondarymapping;
						while(temp_mapping->LossParent != speciestree->root)
						{
							temp_mapping->lossCounter6 = temp_mapping->lossCounter6 + 1;
							temp_mapping->getLossSibling()->lossCounter6 = temp_mapping->getLossSibling()->lossCounter6 + 1;
							
							temp_mapping = temp_mapping->LossParent;
						}
						
						secmapping->LossChild1->lossCounter6 = secmapping->LossChild1->lossCounter6 + 1;
						secmapping->LossChild2->lossCounter6 = secmapping->LossChild2->lossCounter6 + 1;
						temp_mapping = genenode->child(0)->getMapping();
						while(temp_mapping != secmapping)
						{
							temp_mapping->lossCounter1 = temp_mapping->lossCounter1 + 1;
							temp_mapping = temp_mapping->LossParent;
						}
					}					

					else if((secmapping == leftchild_secmapping) && (secmapping != rightchild_secmapping))
					{
						temp_mapping = genenode->secondarymapping;
						while(temp_mapping->LossParent != speciestree->root)
						{
							temp_mapping->lossCounter6 = temp_mapping->lossCounter6 + 1;
							temp_mapping->getLossSibling()->lossCounter6 = temp_mapping->getLossSibling()->lossCounter6 + 1;
							
							temp_mapping = temp_mapping->LossParent;
						}
						
						temp_mapping = genenode->child(1)->secondarymapping;
						sibli = temp_mapping->getLossSibling();
						while(temp_mapping->LossParent != secmapping)
						{
							temp_mapping->lossCounter2 = temp_mapping->lossCounter2 + 1;
							
							temp_mapping = temp_mapping->LossParent;
							sibli = temp_mapping->getLossSibling();
							
						}
						
						sibli->lossCounter6 = sibli->lossCounter6 + 1;
						
					}
					
					else if((secmapping== leftchild_secmapping) && (secmapping == rightchild_secmapping))
					{
						temp_mapping = genenode->secondarymapping;
						while(temp_mapping->LossParent != speciestree->root)
						{
							temp_mapping->lossCounter6 = temp_mapping->lossCounter6 + 1;
							temp_mapping->getLossSibling()->lossCounter6 = temp_mapping->getLossSibling()->lossCounter6 + 1;
							
							temp_mapping = temp_mapping->LossParent;
						}
						
						if(genenode->secondarymapping->LossChild1 != NULL)
						{
							genenode->secondarymapping->LossChild1->lossCounter6 = genenode->secondarymapping->LossChild1->lossCounter6 + 1;
							genenode->secondarymapping->LossChild2->lossCounter6 = genenode->secondarymapping->LossChild2->lossCounter6 + 1;
						}
					}
				}

				// Case 4: g = root, g' = p, g'' = root
				
				else if((mapping == speciestree->root) && (rightchild_mapping == speciestree->root)  && ((leftchild_mapping->no >= subtree->begin) && (leftchild_mapping->no <= subtree->end)))
				{
					
					if(genenode->secondarymapping != genenode->child(1)->secondarymapping)
					{
						cout << "Problem in case 4!!" << endl;
						exit(1);
					}
					
					if(genenode->secondarymapping->LossParent != speciestree->root)
						genenode->secondarymapping->lossCounter4 = genenode->secondarymapping->lossCounter4 + 1;
					

					if(genenode->secondarymapping->LossChild1 != NULL)
					{
						genenode->secondarymapping->LossChild1->lossCounter3 = genenode->secondarymapping->LossChild1->lossCounter3 + 1;
						genenode->secondarymapping->LossChild2->lossCounter3 = genenode->secondarymapping->LossChild2->lossCounter3 + 1;
					}
				}
				
				// case 4 symmetric when g' and g'' are switched 
				else if((mapping == speciestree->root) && (leftchild_mapping == speciestree->root)  && ((rightchild_mapping->no >= subtree->begin) && (rightchild_mapping->no <= subtree->end)))
				{
					
					if(genenode->secondarymapping != genenode->child(0)->secondarymapping)
					{
						cout << "Problem in case 4!!" << endl;
						exit(1);
					}
					
					if(genenode->secondarymapping->LossParent != speciestree->root)
						genenode->secondarymapping->lossCounter4 = genenode->secondarymapping->lossCounter4 + 1;
					
//					genenode->secondarymapping->lossCounter3 = genenode->secondarymapping->lossCounter3 + 1;

					if(genenode->secondarymapping->LossChild1 != NULL)
					{
						genenode->secondarymapping->LossChild1->lossCounter3 = genenode->secondarymapping->LossChild1->lossCounter3 + 1;
						genenode->secondarymapping->LossChild2->lossCounter3 = genenode->secondarymapping->LossChild2->lossCounter3 + 1;
					}
				}
					
				
				// Case 5: g = root, g' = P, g'' = supporting
				
				else if((mapping == speciestree->root) && ((leftchild_mapping->no >= subtree->begin) && (leftchild_mapping->no <= subtree->end)) && ((rightchild_mapping->no < subtree->begin) || (rightchild_mapping->no > subtree->end)))
				{
					
					if(genenode->secondarymapping != rightchild_mapping)
					{
						cout << "Problem in case 5!!" << endl;
						exit(1);
					}
					
					if(genenode->secondarymapping->LossParent != speciestree->root)
						genenode->secondarymapping->lossCounter5 = genenode->secondarymapping->lossCounter5 + 1;
					
					
					// handle the change in he number of losses because the node g becomes a gene duplication
					if(genenode->secondarymapping->LossChild1 != NULL)
					{
						genenode->secondarymapping->LossChild1->lossCounter3 = genenode->secondarymapping->LossChild1->lossCounter3 + 1;
						genenode->secondarymapping->LossChild2->lossCounter3 = genenode->secondarymapping->LossChild2->lossCounter3 + 1;
						genenode->secondarymapping->LossChild1->lossCounter2 = genenode->secondarymapping->LossChild1->lossCounter2 + 1;
						genenode->secondarymapping->LossChild2->lossCounter2 = genenode->secondarymapping->LossChild2->lossCounter2 + 1;
					}
				}
				
				// case 5 symmetric when g' and g'' are switched 
				else if((mapping == speciestree->root) && ((rightchild_mapping->no >= subtree->begin) && (rightchild_mapping->no <= subtree->end)) && ((leftchild_mapping->no < subtree->begin) || (leftchild_mapping->no > subtree->end)))
				{
					
					
					
					if(genenode->secondarymapping != leftchild_mapping)
					{
						cout << "Problem in case 5!!" << endl;
						exit(1);
					}
					
					if(genenode->secondarymapping->LossParent != speciestree->root)
						genenode->secondarymapping->lossCounter5 = genenode->secondarymapping->lossCounter5 + 1;
					
//					genenode->secondarymapping->lossCounter3 = genenode->secondarymapping->lossCounter3 + 1;
					
					// handle the change in he number of losses because the node g becomes a gene duplication
					if(genenode->secondarymapping->LossChild1 != NULL)
					{
						genenode->secondarymapping->LossChild1->lossCounter3 = genenode->secondarymapping->LossChild1->lossCounter3 + 1;
						genenode->secondarymapping->LossChild2->lossCounter3 = genenode->secondarymapping->LossChild2->lossCounter3 + 1;
						genenode->secondarymapping->LossChild1->lossCounter2 = genenode->secondarymapping->LossChild1->lossCounter2 + 1;
						genenode->secondarymapping->LossChild2->lossCounter2 = genenode->secondarymapping->LossChild2->lossCounter2 + 1;
					}
				}

				else
				{
					cout << "Shouldn't be here at all!!" << endl;
					exit(1);
				}

			} 
		}
				
	}
	
	
	
	
	
	
	
	
	
	
	
	void computeGeneLossCounters(GeneTreeUnrooted &tree,  SpeciesNode *subtree)
	{
		vector<GeneNodeUnrooted*> &g = tree.nodes;
		SpeciesNode *temp_mapping, *sibli;

		

		for(vector<GeneNodeUnrooted*>::iterator itr=g.begin(); itr!=g.end(); itr++) {
			GeneNodeUnrooted *genenode = *itr;
			if (!genenode->isLeaf())
			{
				SpeciesNode *&mapping = genenode->getMapping();
				SpeciesNode *&leftchild_mapping = genenode->child(0)->getMapping();
				SpeciesNode *&rightchild_mapping = genenode->child(1)->getMapping();
			
				// start a case by case analysis for the different possible cases
				
				// Case 1: g does not map to the root node
				if(mapping != speciestree->root)
				{
					
					temp_mapping = genenode->child(0)->getMapping();
					while(temp_mapping != mapping)
					{
						temp_mapping->lossCounter1 = temp_mapping->lossCounter1 + 1;
						temp_mapping = temp_mapping->LossParent;
					}
					
					temp_mapping = genenode->child(1)->getMapping();
					while(temp_mapping != mapping)
					{
						temp_mapping->lossCounter1 = temp_mapping->lossCounter1 + 1;
						temp_mapping = temp_mapping->LossParent;
					}
				}

				// Case 2: g, g' and g'' all map to the root
				else if((mapping == speciestree->root) && (rightchild_mapping == speciestree->root) && (leftchild_mapping == speciestree->root))
				{
					SpeciesNode *&secmapping = genenode->secondarymapping;
					SpeciesNode *&leftchild_secmapping = genenode->child(0)->secondarymapping;
					SpeciesNode *&rightchild_secmapping = genenode->child(1)->secondarymapping;
					
					if((secmapping!= leftchild_secmapping) && (secmapping!= rightchild_secmapping))
					{
						temp_mapping = genenode->child(0)->secondarymapping;
						while(temp_mapping != secmapping)
						{
							temp_mapping->lossCounter2 = temp_mapping->lossCounter2 + 1;
//							if(temp_mapping->parent() == secmapping)
//								temp_mapping->lossCounter2 = temp_mapping->lossCounter2 + 1;
							temp_mapping = temp_mapping->LossParent;
						}
			
						temp_mapping = genenode->child(1)->secondarymapping;
						while(temp_mapping != secmapping)
						{
							temp_mapping->lossCounter2 = temp_mapping->lossCounter2 + 1;
//							if(temp_mapping->parent() == secmapping)
//								temp_mapping->lossCounter2 = temp_mapping->lossCounter2 + 1;
							temp_mapping = temp_mapping->LossParent;
						}
					}
					
					else if((secmapping == leftchild_secmapping) && (secmapping!= rightchild_secmapping))
					{
						
						temp_mapping = genenode->child(1)->secondarymapping;
						while(temp_mapping != secmapping)
						{
							temp_mapping->lossCounter2 = temp_mapping->lossCounter2 + 1;
//							if(temp_mapping->parent() == secmapping)
//								temp_mapping->lossCounter2 = temp_mapping->lossCounter2 + 1;
							temp_mapping = temp_mapping->LossParent;
						}
					}
					
					
					else if((secmapping != leftchild_secmapping) && (secmapping == rightchild_secmapping))
					{
						
						temp_mapping = genenode->child(0)->secondarymapping;
						while(temp_mapping != secmapping)
						{
							temp_mapping->lossCounter2 = temp_mapping->lossCounter2 + 1;
//							if(temp_mapping->parent() == secmapping)
//								temp_mapping->lossCounter2 = temp_mapping->lossCounter2 + 1;
							temp_mapping = temp_mapping->LossParent;
						}
					}
					
				}
				
				// Case 3: g = root, g' = root, g'' = supporting
				
				else if((mapping == speciestree->root) && (leftchild_mapping == speciestree->root)  && ((rightchild_mapping->no < subtree->begin) || (rightchild_mapping->no > subtree->end)))
				{
					SpeciesNode *&secmapping = genenode->secondarymapping;
					SpeciesNode *&leftchild_secmapping = genenode->child(0)->secondarymapping;
					SpeciesNode *&rightchild_secmapping = genenode->child(1)->getMapping();

					if((secmapping!= leftchild_secmapping) && (secmapping!= rightchild_secmapping))
					{
						temp_mapping = genenode->secondarymapping;
						while(temp_mapping->LossParent != speciestree->root)
						{
							temp_mapping->lossCounter6 = temp_mapping->lossCounter6 + 1;
							//SpeciesNode *sib = temp_mapping->getSibling();
							temp_mapping->getLossSibling()->lossCounter6 = temp_mapping->getLossSibling()->lossCounter6 + 1;
							temp_mapping = temp_mapping->LossParent;
						}
    					secmapping->LossChild1->lossCounter6 = secmapping->LossChild1->lossCounter6 + 1;
						secmapping->LossChild2->lossCounter6 = secmapping->LossChild2->lossCounter6 + 1;
						temp_mapping = genenode->child(0)->secondarymapping;
						while(temp_mapping->LossParent != secmapping)
						{
							temp_mapping->lossCounter2 = temp_mapping->lossCounter2 + 1;
//							temp_mapping->lossCounter1 = temp_mapping->lossCounter1 + 1;
							temp_mapping = temp_mapping->LossParent;
						}
						temp_mapping->lossCounter6 = temp_mapping->lossCounter6 + 1;
						temp_mapping = genenode->child(1)->getMapping();
						while(temp_mapping != secmapping)
						{
							temp_mapping->lossCounter1 = temp_mapping->lossCounter1 + 1;
							temp_mapping = temp_mapping->LossParent;
						}
					}
					
					else if((secmapping == leftchild_secmapping) && (secmapping != rightchild_secmapping))
					{
						
					
						temp_mapping = genenode->secondarymapping;
						while(temp_mapping->LossParent != speciestree->root)
						{
							temp_mapping->lossCounter6 = temp_mapping->lossCounter6 + 1;
							temp_mapping->getLossSibling()->lossCounter6 = temp_mapping->getLossSibling()->lossCounter6 + 1;

							temp_mapping = temp_mapping->LossParent;
							
						}
						
						secmapping->LossChild1->lossCounter6 = secmapping->LossChild1->lossCounter6 + 1;
						secmapping->LossChild2->lossCounter6 = secmapping->LossChild2->lossCounter6 + 1;
						temp_mapping = genenode->child(1)->getMapping();
						
						while(temp_mapping != secmapping)
						{
							temp_mapping->lossCounter1 = temp_mapping->lossCounter1 + 1;
							
							temp_mapping = temp_mapping->LossParent;
						}
					}					
					
					else if((secmapping != leftchild_secmapping) && (secmapping == rightchild_secmapping))
					{
						
						temp_mapping = genenode->secondarymapping;
						while(temp_mapping->LossParent != speciestree->root)
						{
							temp_mapping->lossCounter6 = temp_mapping->lossCounter6 + 1;
							temp_mapping->getLossSibling()->lossCounter6 = temp_mapping->getLossSibling()->lossCounter6 + 1;
							
							temp_mapping = temp_mapping->LossParent;
						}
						
						temp_mapping = genenode->child(0)->secondarymapping;
						sibli = temp_mapping->getLossSibling();
						while(temp_mapping->LossParent != secmapping)
						{
							temp_mapping->lossCounter2 = temp_mapping->lossCounter2 + 1;
							temp_mapping = temp_mapping->LossParent;
							sibli = temp_mapping->getLossSibling();
						}
						
						sibli->lossCounter6 = sibli->lossCounter6 + 1;
						
					}
					
					else if((secmapping == leftchild_secmapping) && (secmapping == rightchild_secmapping))
					{
						
						temp_mapping = genenode->secondarymapping;
						
						while(temp_mapping->LossParent != speciestree->root)
						{
							
							temp_mapping->lossCounter6 = temp_mapping->lossCounter6 + 1;
							temp_mapping->getLossSibling()->lossCounter6 = temp_mapping->getLossSibling()->lossCounter6 + 1;
							
							temp_mapping = temp_mapping->LossParent;
							
						}
					
						
						if(genenode->secondarymapping->LossChild1 != NULL)
						{
							genenode->secondarymapping->LossChild1->lossCounter6 = genenode->secondarymapping->LossChild1->lossCounter6 + 1;
							genenode->secondarymapping->LossChild2->lossCounter6 = genenode->secondarymapping->LossChild2->lossCounter6 + 1;
						}
						
					}
				}

				// Case 3: symmetric case when g' and g'' are switched
				else if((mapping == speciestree->root) && (rightchild_mapping == speciestree->root)  && ((leftchild_mapping->no < subtree->begin) || (leftchild_mapping->no > subtree->end)))
				{
					SpeciesNode *&secmapping = genenode->secondarymapping;
					SpeciesNode *&leftchild_secmapping = genenode->child(0)->getMapping();
					SpeciesNode *&rightchild_secmapping = genenode->child(1)->secondarymapping;

					if((secmapping!= leftchild_secmapping) && (secmapping!= rightchild_secmapping))
					{
						temp_mapping = genenode->secondarymapping;
						while(temp_mapping->LossParent != speciestree->root)
						{
							temp_mapping->lossCounter6 = temp_mapping->lossCounter6 + 1;
							temp_mapping->getLossSibling()->lossCounter6 = temp_mapping->getLossSibling()->lossCounter6 + 1;
							
							temp_mapping = temp_mapping->LossParent;
						}
						secmapping->LossChild1->lossCounter6 = secmapping->LossChild1->lossCounter6 + 1;
						secmapping->LossChild2->lossCounter6 = secmapping->LossChild2->lossCounter6 + 1;
						
						temp_mapping = genenode->child(1)->secondarymapping;
						while(temp_mapping->LossParent != secmapping)
						{
							temp_mapping->lossCounter2 = temp_mapping->lossCounter2 + 1;
//							temp_mapping->lossCounter1 = temp_mapping->lossCounter1 + 1;
							temp_mapping = temp_mapping->LossParent;
						}
						temp_mapping->lossCounter6 = temp_mapping->lossCounter6 + 1;

						temp_mapping = genenode->child(0)->getMapping();
						while(temp_mapping != secmapping)
						{
							temp_mapping->lossCounter1 = temp_mapping->lossCounter1 + 1;
							temp_mapping = temp_mapping->LossParent;
						}
						
					}
					
					else if((secmapping != leftchild_secmapping) && (secmapping == rightchild_secmapping))
					{
						temp_mapping = genenode->secondarymapping;
						while(temp_mapping->LossParent != speciestree->root)
						{
							temp_mapping->lossCounter6 = temp_mapping->lossCounter6 + 1;
							temp_mapping->getLossSibling()->lossCounter6 = temp_mapping->getLossSibling()->lossCounter6 + 1;
							
							temp_mapping = temp_mapping->LossParent;
						}
						
						secmapping->LossChild1->lossCounter6 = secmapping->LossChild1->lossCounter6 + 1;
						secmapping->LossChild2->lossCounter6 = secmapping->LossChild2->lossCounter6 + 1;
						temp_mapping = genenode->child(0)->getMapping();
						while(temp_mapping != secmapping)
						{
							temp_mapping->lossCounter1 = temp_mapping->lossCounter1 + 1;
							temp_mapping = temp_mapping->LossParent;
						}
					}					

					else if((secmapping == leftchild_secmapping) && (secmapping != rightchild_secmapping))
					{
						temp_mapping = genenode->secondarymapping;
						while(temp_mapping->LossParent != speciestree->root)
						{
							temp_mapping->lossCounter6 = temp_mapping->lossCounter6 + 1;
							temp_mapping->getLossSibling()->lossCounter6 = temp_mapping->getLossSibling()->lossCounter6 + 1;
							
							temp_mapping = temp_mapping->LossParent;
						}
						
						temp_mapping = genenode->child(1)->secondarymapping;
						sibli = temp_mapping->getLossSibling();
						while(temp_mapping->LossParent != secmapping)
						{
							temp_mapping->lossCounter2 = temp_mapping->lossCounter2 + 1;
							
							temp_mapping = temp_mapping->LossParent;
							sibli = temp_mapping->getLossSibling();
							
						}
						
						sibli->lossCounter6 = sibli->lossCounter6 + 1;
						
					}
					
					else if((secmapping== leftchild_secmapping) && (secmapping == rightchild_secmapping))
					{
						temp_mapping = genenode->secondarymapping;
						while(temp_mapping->LossParent != speciestree->root)
						{
							temp_mapping->lossCounter6 = temp_mapping->lossCounter6 + 1;
							temp_mapping->getLossSibling()->lossCounter6 = temp_mapping->getLossSibling()->lossCounter6 + 1;
							
							temp_mapping = temp_mapping->LossParent;
						}
						
						if(genenode->secondarymapping->LossChild1 != NULL)
						{
							genenode->secondarymapping->LossChild1->lossCounter6 = genenode->secondarymapping->LossChild1->lossCounter6 + 1;
							genenode->secondarymapping->LossChild2->lossCounter6 = genenode->secondarymapping->LossChild2->lossCounter6 + 1;
						}
					}
				}

				// Case 4: g = root, g' = p, g'' = root
				
				else if((mapping == speciestree->root) && (rightchild_mapping == speciestree->root)  && ((leftchild_mapping->no >= subtree->begin) && (leftchild_mapping->no <= subtree->end)))
				{
					
					if(genenode->secondarymapping != genenode->child(1)->secondarymapping)
					{
						cout << "Problem in case 4!!" << endl;
						exit(1);
					}
					
					if(genenode->secondarymapping->LossParent != speciestree->root)
						genenode->secondarymapping->lossCounter4 = genenode->secondarymapping->lossCounter4 + 1;
					
//					genenode->secondarymapping->lossCounter3 = genenode->secondarymapping->lossCounter3 + 1;

					if(genenode->secondarymapping->LossChild1 != NULL)
					{
						genenode->secondarymapping->LossChild1->lossCounter3 = genenode->secondarymapping->LossChild1->lossCounter3 + 1;
						genenode->secondarymapping->LossChild2->lossCounter3 = genenode->secondarymapping->LossChild2->lossCounter3 + 1;
					}
				}
				
				// case 4 symmetric when g' and g'' are switched 
				else if((mapping == speciestree->root) && (leftchild_mapping == speciestree->root)  && ((rightchild_mapping->no >= subtree->begin) && (rightchild_mapping->no <= subtree->end)))
				{
					
					if(genenode->secondarymapping != genenode->child(0)->secondarymapping)
					{
						cout << "Problem in case 4!!" << endl;
						exit(1);
					}
					
					if(genenode->secondarymapping->LossParent != speciestree->root)
						genenode->secondarymapping->lossCounter4 = genenode->secondarymapping->lossCounter4 + 1;
					
//					genenode->secondarymapping->lossCounter3 = genenode->secondarymapping->lossCounter3 + 1;

					if(genenode->secondarymapping->LossChild1 != NULL)
					{
						genenode->secondarymapping->LossChild1->lossCounter3 = genenode->secondarymapping->LossChild1->lossCounter3 + 1;
						genenode->secondarymapping->LossChild2->lossCounter3 = genenode->secondarymapping->LossChild2->lossCounter3 + 1;
					}
				}
					
				
				// Case 5: g = root, g' = P, g'' = supporting
				
				else if((mapping == speciestree->root) && ((leftchild_mapping->no >= subtree->begin) && (leftchild_mapping->no <= subtree->end)) && ((rightchild_mapping->no < subtree->begin) || (rightchild_mapping->no > subtree->end)))
				{
					
					if(genenode->secondarymapping != rightchild_mapping)
					{
						cout << "Problem in case 5!!" << endl;
						exit(1);
					}
					
					if(genenode->secondarymapping->LossParent != speciestree->root)
						genenode->secondarymapping->lossCounter5 = genenode->secondarymapping->lossCounter5 + 1;
					
//					genenode->secondarymapping->lossCounter3 = genenode->secondarymapping->lossCounter3 + 1;
					
					// handle the change in he number of losses because the node g becomes a gene duplication
					if(genenode->secondarymapping->LossChild1 != NULL)
					{
						genenode->secondarymapping->LossChild1->lossCounter3 = genenode->secondarymapping->LossChild1->lossCounter3 + 1;
						genenode->secondarymapping->LossChild2->lossCounter3 = genenode->secondarymapping->LossChild2->lossCounter3 + 1;
						genenode->secondarymapping->LossChild1->lossCounter2 = genenode->secondarymapping->LossChild1->lossCounter2 + 1;
						genenode->secondarymapping->LossChild2->lossCounter2 = genenode->secondarymapping->LossChild2->lossCounter2 + 1;
					}
				}
				
				// case 5 symmetric when g' and g'' are switched 
				else if((mapping == speciestree->root) && ((rightchild_mapping->no >= subtree->begin) && (rightchild_mapping->no <= subtree->end)) && ((leftchild_mapping->no < subtree->begin) || (leftchild_mapping->no > subtree->end)))
				{
					
					
					
					if(genenode->secondarymapping != leftchild_mapping)
					{
						cout << "Problem in case 5!!" << endl;
						exit(1);
					}
					
					if(genenode->secondarymapping->LossParent != speciestree->root)
						genenode->secondarymapping->lossCounter5 = genenode->secondarymapping->lossCounter5 + 1;
					
//					genenode->secondarymapping->lossCounter3 = genenode->secondarymapping->lossCounter3 + 1;
					
					// handle the change in he number of losses because the node g becomes a gene duplication
					if(genenode->secondarymapping->LossChild1 != NULL)
					{
						genenode->secondarymapping->LossChild1->lossCounter3 = genenode->secondarymapping->LossChild1->lossCounter3 + 1;
						genenode->secondarymapping->LossChild2->lossCounter3 = genenode->secondarymapping->LossChild2->lossCounter3 + 1;
						genenode->secondarymapping->LossChild1->lossCounter2 = genenode->secondarymapping->LossChild1->lossCounter2 + 1;
						genenode->secondarymapping->LossChild2->lossCounter2 = genenode->secondarymapping->LossChild2->lossCounter2 + 1;
					}
				}

				else
				{
					cout << "Shouldn't be here at all!!" << endl;
					exit(1);
				}

			} 
		}
				
	}
	
	
	
	
	
	
	
	
	
	
	void computeLossScores(SpeciesNode *&node, const double &weight)
	{
//cout << *node << endl;	
		// convert counter 5 into counter 3 and counter 6
		convertCounterFive(node);
//cout << *node << endl;		
		// handle counter types 2 and 6, and also 1
		PreOrderCounterFirst(node, 0);
//cout << *node << endl;		
		// convert counter 4 into counter 3
		convertCounterFour(node); 
//cout << *node << endl;		
		// handle counter 3 
		PreOrderCounterThreeStep1(node, 0);
//cout << *node << endl;		
		PreOrderCounterThreeStep2(node, 0);
//cout << *node << endl;
		// copies the loss values computed at the relevant nodes to all the other nodes of the species tree.
		TransferLossScoretoFullTree(node, weight);
		
	}
	
	void PreOrderCounterFirst(SpeciesNode *&node, int c)
	{

		if((node->LossParent == speciestree->root) && ((node->lossCounter2 !=0) || (node->lossCounter1 != 0) || (node->lossCounter6 != 0)))
		{
			cout << "root-child counters not 0!!" << endl;
			exit(0);
		}

		if (node->isRelevant == true)
		{
			c =  c + node->lossCounter2 - node->lossCounter6;
			node->LossScoreDiff = node->LossScoreDiff + c + node->lossCounter1;
		}
		else {
			cout << "isRelevant is false!" << endl;
			exit(0);
		}
		
		if (node->LossChild1 != NULL) {
			PreOrderCounterFirst(node->LossChild1, c);
		}
		if (node->LossChild2 != NULL) {
			PreOrderCounterFirst(node->LossChild2, c);
		}
	}
	

	// NOTE: this function does not distinguish between relevant and non-relevant nodes... TO BE FIXED
	void convertCounterFour(SpeciesNode *&node)
	{
	
		if (node->LossChild1 != NULL) {
			convertCounterFour(node->LossChild1);
		}
		if (node->LossChild2 != NULL) {
			convertCounterFour(node->LossChild2);
		}
		
		if(node->LossParent !=speciestree->root)
		{
			
			if(node->LossChild1 == NULL)
			{
				SpeciesNode *sib = node->getLossSibling();
				sib->lossCounter3 = sib->lossCounter3 + node->lossCounter4;
			}
			
			else
			{
				node->lossCounter4 = node->lossCounter4 + node->LossChild1->lossCounter4 + node->LossChild2->lossCounter4;
				SpeciesNode *sib = node->getLossSibling();
				sib->lossCounter3 = sib->lossCounter3 + node->lossCounter4;
			}
		}
				
		
		
		
/*		if ((node->parent() == speciestree->root) && (node->child(0) != NULL) && (node->child(1) != NULL))
		{
//			node->lossCounter4 = node->child(0)->lossCounter4 + node->child(1)->lossCounter4;
//			SpeciesNode *sib = node->getSibling();
//			sib->lossCounter3 = sib->lossCounter3 + node->lossCounter4;
		}
		
		else if((node->child(0) != NULL) && (node->child(1) != NULL))
		{
			
			node->lossCounter4 = node->lossCounter4 + node->child(0)->lossCounter4 + node->child(1)->lossCounter4;
			SpeciesNode *sib = node->getSibling();
			sib->lossCounter3 = sib->lossCounter3 + node->lossCounter4;
		}
		else if(node->child(0) == NULL)
		{
			SpeciesNode *sib = node->getSibling();
			sib->lossCounter3 = sib->lossCounter3 + node->lossCounter4;
		}
*/
	}
	
		// NOTE: this function does not distinguish between relevant and non-relevant nodes... TO BE FIXED
	void convertCounterFive(SpeciesNode *&node)
	{
		if (node->LossChild1 != NULL) {
			convertCounterFive(node->LossChild1);
		}
		if (node->LossChild2 != NULL) {
			convertCounterFive(node->LossChild2);
		}
		
		
		
		if(node->LossParent !=speciestree->root)
		{
			
			if(node->LossChild1 == NULL)
			{
				SpeciesNode *sib = node->getLossSibling();
				sib->lossCounter3 = sib->lossCounter3 + node->lossCounter5;
				sib->lossCounter6 = sib->lossCounter6 + node->lossCounter5;
				node->lossCounter6 = node->lossCounter6 + node->lossCounter5;
			}
			
			else
			{
				node->lossCounter5 = node->lossCounter5 + node->LossChild1->lossCounter5 + node->LossChild2->lossCounter5;
				SpeciesNode *sib = node->getLossSibling();
				sib->lossCounter3 = sib->lossCounter3 + node->lossCounter5;
				sib->lossCounter6 = sib->lossCounter6 + node->lossCounter5;
				node->lossCounter6 = node->lossCounter6 + node->lossCounter5;
			}
		}
		
		
/*		if ((node->parent() == speciestree->root) && (node->child(0) != NULL) && (node->child(1) != NULL))
		{
//			node->lossCounter5 = node->child(0)->lossCounter5 + node->child(1)->lossCounter5;
//			SpeciesNode *sib = node->getSibling();
//			sib->lossCounter3 = sib->lossCounter3 + node->lossCounter5;
//			node->lossCounter6 = node->lossCounter6 + node->lossCounter5;
		}
		
		else if((node->child(0) != NULL) && (node->child(1) != NULL))
		{
			
			node->lossCounter5 = node->lossCounter5 + node->child(0)->lossCounter5 + node->child(1)->lossCounter5;
			SpeciesNode *sib = node->getSibling();
			sib->lossCounter3 = sib->lossCounter3 + node->lossCounter5;
			node->lossCounter6 = node->lossCounter6 + node->lossCounter5;
		}
		else if(node->child(0) == NULL)
		{
			SpeciesNode *sib = node->getSibling();
			sib->lossCounter3 = sib->lossCounter3 + node->lossCounter5;
			node->lossCounter6 = node->lossCounter6 + node->lossCounter5;
		}
*/		
		
	}

		// NOTE: this function does not distinguish between relevant and non-relevant nodes... TO BE FIXED
	void PreOrderCounterThreeStep1(SpeciesNode *&node, int c)
	{
		
		if((node->LossParent == speciestree->root) && (node->lossCounter3 !=0))
		{
			cout << "root-child counter 3 not 0!!" << endl;
			exit(0);
		}
			
		c =  c + node->lossCounter3;
		node->lossCounter3 = c;
		
		if (node->LossChild1 != NULL) {
			PreOrderCounterThreeStep1(node->LossChild1, c);
		}
		if (node->LossChild2 != NULL) {
			PreOrderCounterThreeStep1(node->LossChild2, c);
		}
	}

	// NOTE: this function does not distinguish between relevant and non-relevant nodes... TO BE FIXED
	void PreOrderCounterThreeStep2(SpeciesNode *&node, int c)
	{
		
		if((node->LossParent == speciestree->root) && (node->lossCounter3 !=0))
		{
			cout << "root-child counter 3 not 0!!" << endl;
			exit(0);
		}
		c =  c + node->lossCounter3;
		node->LossScoreDiff = node->LossScoreDiff + c;
		
		if (node->LossChild1 != NULL) {
			PreOrderCounterThreeStep2(node->LossChild1, c);
		}
		if (node->LossChild2 != NULL) {
			PreOrderCounterThreeStep2(node->LossChild2, c);
		}
	}


	void TransferLossScoretoFullTree(SpeciesNode *&node, const double &weight)
	{
		if (node->child(0) != NULL) {
			TransferLossScoretoFullTree(node->child(0), weight);
		}
		if (node->child(1) != NULL) {
			TransferLossScoretoFullTree(node->child(1), weight);
		}
		
		if(node->isRelevant == true)
		{
			node->lossScore = node->lossScore + (weight * node->LossScoreDiff);
//			if(node->LossParent != speciestree->root)
			if(node->parent() != speciestree->root)
			{		

				SpeciesNode* Temp_Node = node->parent();
				while(Temp_Node->isRelevant == false)
				{
					Temp_Node->lossScore = Temp_Node->lossScore + (weight * node->LossScoreDiff);
					if(Temp_Node->child(0)->subtreeSize == 0)
						TransferLossScoreStep2(Temp_Node->child(0), (weight * node->LossScoreDiff));
					else 
						TransferLossScoreStep2(Temp_Node->child(1), (weight * node->LossScoreDiff));
							
					Temp_Node = Temp_Node->parent();
				}
			}
		}
	}
		
	void TransferLossScoreStep2(SpeciesNode *&tmpnode, double diff)
	{
		tmpnode->lossScore = tmpnode->lossScore + diff;
		
		if (tmpnode->child(0) != NULL) {
			TransferLossScoreStep2(tmpnode->child(0), diff);
		}
		if (tmpnode->child(1) != NULL) {
			TransferLossScoreStep2(tmpnode->child(1), diff);
		}
	}
		
	

	
// -----------------Losses stuff end-----------------------------------------------------------------------


	// find the best rooting of the unrooted gene trees
	void computeBestRooting() {
		speciestree->establishOrder();
		speciestree->preprocessLCA();
		// process all unrooted trees
		for(int i=0; i<genetree_unrooted.size(); i++) {
			GeneTreeUnrooted &tree = *genetree_unrooted[i];
			createPrimaryMappingUnrooted(tree);
			
			resetRelevantTree(speciestree->root);
			buildRelevantTree(tree);
			double &score = best_score;
			score = getScore(tree) * tree.weight + computeGeneLossForRoot(tree) * tree.weight;
			GeneNodeUnrooted *u = tree.root->child(0);
			GeneNodeUnrooted *v = tree.root->child(1);
			best_node[0] = u;
			best_node[1] = v;
			#ifdef DEBUG
			position_counter_debug = 2;
			if ((u == NULL) || (v == NULL)) EXCEPTION("child = NULL in computeBestRooting" << endl);
			#endif
			moveRoot2(tree, u, score);
			moveRoot2(tree, v, score);
			#ifdef DEBUG
			if (position_counter_debug != tree.nodes.size()-1)
				WARNING("tree traversal failed in computeBestRooting" << position_counter_debug << " != " << tree.nodes.size()-1 << endl);
			#endif
			tree.reroot(best_node[0], best_node[1]);
		}
		speciestree->postprocessLCA();
	}

	// calculate the gene duplication for all rootings and rSPR operation
	#ifdef DEBUG
	int position_counter_debug;
	#endif
	double best_score;
	GeneNodeUnrooted *best_node[2];
	inline void moveRoot(GeneTreeUnrooted &tree, GeneNodeUnrooted *p, SpeciesNode *subtree, SpeciesNode *sibling) {
		GeneNodeUnrooted *c[2];
		SpeciesNode *LossSibling;
		for (int i=0; i<2; i++) c[i] = p->child(i);
		for (int i=0; i<2; i++) {
			if (c[i] != NULL) {
				moveRoot(tree, c[i], subtree, sibling);
				tree.reroot(c[i], p);
				#ifdef DEBUG
				tree.checkStructure();
				position_counter_debug++;
				#endif
				
				createSecondaryMapping(tree, subtree);
				computeGeneDuplicationsTriple();
				resetRelevantTree(speciestree->root);
				buildRelevantTree(tree);
				const double score = getScore(tree) * tree.weight;
				double initLoss  = tree.weight * computeGeneLossForRoot(tree);
				if (score + initLoss < best_score) {
					best_score = score + initLoss;
					best_node[1] = c[i];
					best_node[0] = p;
				}

				
				
			// Loss stuff begin	
				if(speciestree->root->isRelevant == false)
				{
					// If the gene tree is fully contained in either the pruned subtree or in the other subtree then there is nothing to be done
				}
				else
				{
					computeGeneLossCounters(tree, subtree);
			
					if(speciestree->root->child(0) == subtree)
						LossSibling = speciestree->root->LossChild2;
					else LossSibling = speciestree->root->LossChild1;
					computeLossScores_Temp2(LossSibling);				
				}
				
			// loss stuff end
				
				removeSecondaryMapping(tree);
				
				
				
				
				
				computeGeneDuplicationsTempMin2(sibling, score, initLoss, tree.weight);
			}
		}
	}

	// calculate the gene duplication for all rootings
	inline void moveRoot2(GeneTreeUnrooted &tree, GeneNodeUnrooted *p, double &best_score) {
		GeneNodeUnrooted *c[2];
		for (int i=0; i<2; i++) c[i] = p->child(i);
		for (int i=0; i<2; i++) {
			if (c[i] != NULL) {
				moveRoot2(tree, c[i], best_score);
				tree.reroot(c[i], p);
				#ifdef DEBUG
				tree.checkStructure();
				position_counter_debug++;
				#endif
				
				
				
				const double score = getScore(tree) * tree.weight + computeGeneLossForRoot(tree) * tree.weight;
				

				
				
				if (score < best_score) {
					best_score = score;
					best_node[1] = c[i];
					best_node[0] = p;
				}
			}
		}
	}

	// reset the gene duplication score to 0
	inline void resetGeneDuplications(SpeciesNode *&node) {
		if (node == NULL) return;
		node->score = 0;
		for (int i=0; i<2; i++)
			resetGeneDuplications(node->child(i));
	}

	// compute the gene duplication score and add it to the current score
	inline void computeGeneDuplicationsAdd(SpeciesNode *&node, const double &score, const double &weight) {
		if (node == NULL) return;
		node->score += score;
		for (int i=0; i<2; i++) {
			computeGeneDuplicationsAdd(node->child(i), score + (double(node->gain) - double(node->lost[i])) * weight, weight);
		}
	}

	// compute the gene duplication score into a temporary variable
	inline void computeGeneDuplicationsTempReplace(SpeciesNode *&node, const double &score, const double &weight) {
		if (node == NULL) return;
		node->tempscore = score;
		for (int i=0; i<2; i++)
			computeGeneDuplicationsTempReplace(node->child(i), score + (double(node->gain) - double(node->lost[i])) * weight, weight);
	}
	// compute the gene duplication score into a temporary variable (replace is smaller than current score)
	inline void computeGeneDuplicationsTempMin(SpeciesNode *&node, const double &score, const double &weight) {
		if (node == NULL) return;
		if (score < node->tempscore) node->tempscore = score;
		for (int i=0; i<2; i++)
			computeGeneDuplicationsTempMin(node->child(i), score + (double(node->gain) - double(node->lost[i])) * weight, weight);
	}
	
	
	
	inline void computeGeneDuplicationsTempMin2(SpeciesNode *&node, const double &score, double &initLoss, const double &weight) {
		if (node == NULL) return;
		if ((score + initLoss + (weight * node->LossScoreDiff)) < (node->tempscore + node->lossScoreTemp)) 
		{
			node->tempscore = score;
			node->lossScoreTemp = initLoss + (weight * node->LossScoreDiff);
		}
		for (int i=0; i<2; i++)
			computeGeneDuplicationsTempMin2(node->child(i), score + (double(node->gain) - double(node->lost[i])) * weight, initLoss, weight);
	}
	
	
	
	// add the temporary gene duplication score to the final
	inline void addTempGeneDuplications(SpeciesNode *&node) {
		if (node == NULL) return;
		node->score += node->tempscore;
		node->lossScore = node->lossScore + node->lossScoreTemp;
		for (int i=0; i<2; i++)
			addTempGeneDuplications(node->child(i));
	}

	
//	----------------Some more Loss Stuff -------------------
	
	void copyInitialLossScoreToAllNodes_Temp(SpeciesNode *&node, double initLoss)
	{
		if (node == NULL) return;
		node->lossScoreTemp = initLoss;
		for (int i=0; i<2; i++)
			copyInitialLossScoreToAllNodes_Temp(node->child(i), initLoss);
	}
	
	
	
	
	void computeLossScores_Temp(SpeciesNode *&node, const double &weight)
	{
//cout << *node << endl;	
		// convert counter 5 into counter 3 and counter 6
		convertCounterFive(node);
//cout << *node << endl;		
		// handle counter types 2 and 6, and also 1
		PreOrderCounterFirst(node, 0);
//cout << *node << endl;		
		// convert counter 4 into counter 3
		convertCounterFour(node); 
//cout << *node << endl;		
		// handle counter 3 
		PreOrderCounterThreeStep1(node, 0);
//cout << *node << endl;		
		PreOrderCounterThreeStep2(node, 0);
//cout << *node << endl;
		// copies the loss values computed at the relevant nodes to all the other nodes of the species tree.
		TransferLossScoretoFullTree_Temp(node, weight);
		
	}
	
	
	
		void TransferLossScoretoFullTree_Temp(SpeciesNode *&node, const double &weight)
	{
		if (node->child(0) != NULL) {
			TransferLossScoretoFullTree_Temp(node->child(0), weight);
		}
		if (node->child(1) != NULL) {
			TransferLossScoretoFullTree_Temp(node->child(1), weight);
		}
		
		if(node->isRelevant == true)
		{
			node->lossScoreTemp = node->lossScoreTemp + (weight * node->LossScoreDiff);
			if(node->parent() != speciestree->root)
			{		

				SpeciesNode* Temp_Node = node->parent();
				while(Temp_Node->isRelevant == false)
				{
					Temp_Node->lossScoreTemp = Temp_Node->lossScoreTemp + (weight * node->LossScoreDiff);
					if(Temp_Node->child(0)->subtreeSize == 0)
						TransferLossScoreStep2_Temp(Temp_Node->child(0), (weight * node->LossScoreDiff));
					else 
						TransferLossScoreStep2_Temp(Temp_Node->child(1), (weight * node->LossScoreDiff));
							
					Temp_Node = Temp_Node->parent();
				}
			}
		}
	}
		
	void TransferLossScoreStep2_Temp(SpeciesNode *&tmpnode, double diff)
	{
		tmpnode->lossScoreTemp = tmpnode->lossScoreTemp + diff;
		
		if (tmpnode->child(0) != NULL) {
			TransferLossScoreStep2_Temp(tmpnode->child(0), diff);
		}
		if (tmpnode->child(1) != NULL) {
			TransferLossScoreStep2_Temp(tmpnode->child(1), diff);
		}
	}
	
	
	
	
	void computeLossScores_Temp2(SpeciesNode *&node)
	{
//cout << *node << endl;	
		// convert counter 5 into counter 3 and counter 6
		convertCounterFive(node);
//cout << *node << endl;		
		// handle counter types 2 and 6, and also 1
		PreOrderCounterFirst(node, 0);
//cout << *node << endl;		
		// convert counter 4 into counter 3
		convertCounterFour(node); 
//cout << *node << endl;		
		// handle counter 3 
		PreOrderCounterThreeStep1(node, 0);
//cout << *node << endl;		
		PreOrderCounterThreeStep2(node, 0);
//cout << *node << endl;
		// copies the loss values computed at the relevant nodes to all the other nodes of the species tree.
		TransferLossScoretoFullTree_Temp2(node);
		
	}
	
	
	void TransferLossScoretoFullTree_Temp2(SpeciesNode *&node)
	{
		if (node->child(0) != NULL) {
			TransferLossScoretoFullTree_Temp2(node->child(0));
		}
		if (node->child(1) != NULL) {
			TransferLossScoretoFullTree_Temp2(node->child(1));
		}
		
		if(node->isRelevant == true)
		{
//			node->lossScoreTemp = node->lossScoreTemp + node->LossScoreDiff;
			if(node->parent() != speciestree->root)
			{		

				SpeciesNode* Temp_Node = node->parent();
				while(Temp_Node->isRelevant == false)
				{
					Temp_Node->LossScoreDiff = node->LossScoreDiff;
					if(Temp_Node->child(0)->subtreeSize == 0)
						TransferLossScoreStep2_Temp2(Temp_Node->child(0), node->LossScoreDiff);
					else 
						TransferLossScoreStep2_Temp2(Temp_Node->child(1), node->LossScoreDiff);
							
					Temp_Node = Temp_Node->parent();
				}
			}
		}
	}
		
	void TransferLossScoreStep2_Temp2(SpeciesNode *&tmpnode, int diff)
	{
		tmpnode->LossScoreDiff = diff;
		
		if (tmpnode->child(0) != NULL) {
			TransferLossScoreStep2_Temp2(tmpnode->child(0), diff);
		}
		if (tmpnode->child(1) != NULL) {
			TransferLossScoreStep2_Temp2(tmpnode->child(1), diff);
		}
	}
	
	
// ---------------------Loss stuff end -------------------	
	
	
	
	// travers the tree and call the virtual function scoreComputed for valid rSPR operations
	inline void forEachCallScoreComputed(SpeciesNode *&subtree, SpeciesNode *&node) {
		if (node == NULL) return;
		// call scoreComputed when constraints are not violated
		const unsigned int color = subtree->parent()->constraint;
		if (color == node->constraint) {
			scoreComputed(*node);
		} else
		if ((subtree->parent() != node->parent()) && (color == node->parent()->constraint)) {
			scoreComputed(*node);
		}
		for (int i=0; i<2; i++)
			forEachCallScoreComputed(subtree, node->child(i));
	}

	virtual void scoreComputed(SpeciesNode &node) = 0;
	virtual void run(ostream &os, const ReRoot reroot) = 0;
	virtual double getCurrentScore() = 0;

	// ------------------------------------------------------------------------------------------------------
	// outputs the trees into a string stream
	friend ostream & operator << (ostream &os, Heuristic &t);
};

// outputs the trees into a string stream
ostream & operator << (ostream &os, Heuristic &t) {
	os << (TreeSet&) t;
	return os;
}

#endif
