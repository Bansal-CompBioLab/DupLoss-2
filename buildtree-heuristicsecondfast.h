/*Copyright (C) 2024 Mukul S. Bansal (mukul.bansal@uconn.edu).
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

#ifndef HEURISTICSECONDFAST_H

#define TreeSet LEAFADD::TreeSet
#define SpeciesNode LEAFADD::SpeciesNode
#define NamedSpeciesNode LEAFADD::NamedSpeciesNode
#define GeneNodeUnrooted LEAFADD::GeneNodeUnrooted
#define NamedGeneNodeUnrooted LEAFADD::NamedGeneNodeUnrooted
#define GeneTreeUnrooted LEAFADD::GeneTreeUnrooted
#define GeneNodeRooted LEAFADD::GeneNodeRooted
#define NamedGeneNodeRooted LEAFADD::NamedGeneNodeRooted
#define GeneTreeRooted LEAFADD::GeneTreeRooted

class HeuristicSecondFast :  public HeuristicLeafAdd {
protected:
        unsigned int Best_score;
	bool update, tripFlag;
	struct data
	{
		int BestSubtreeRoot;
		SpeciesNode *BestNewLocation;
	} temp, old;

	struct trip_data
	{
		int i_val;
		int l_val;
		int k_val;
	} trip_temp;
	vector<data> queue;
	vector<trip_data> trip_queue;
	int j, i, k, l, currentconstraint, index, index_extra, index_extra2, index_extra3;



public:
        inline void scoreComputed(SpeciesNode &node)
	{

                const int &genedup = node.genedup;
				const int &geneloss = node.lossScore;
		if((node.constraint == currentconstraint) || (node.parent()->constraint == currentconstraint))
		{
	                if(genedup + geneloss < Best_score)
        	        {

				if(tripFlag)
				{

					update = true;
					trip_queue.clear();
					trip_temp.i_val = i;
        	                	trip_temp.l_val = l;
					trip_temp.k_val = k;
					trip_queue.push_back(trip_temp);
                	        	Best_score = genedup + geneloss;
				}

				else
				{
					update = true;
					queue.clear();
					temp.BestSubtreeRoot = j;
               			        temp.BestNewLocation = &node;
					queue.push_back(temp);
	                        	Best_score = genedup + geneloss;

//                        msgout<< "\rCurrent best score: " <<Best_score<<"      ";
					flush(cout);
				}
                	}
			else if ((update== true) && ((genedup + geneloss) == Best_score))
			{

				if(tripFlag)
				{
					trip_temp.i_val = i;
	                        	trip_temp.l_val = l;
					trip_temp.k_val = k;

                	        	trip_queue.push_back(trip_temp);
				}
				else
				{

					temp.BestSubtreeRoot = j;
        	                	temp.BestNewLocation = &node;
                		     	queue.push_back(temp);
				}
			}

        	}
	}

        HeuristicSecondFast()
        {
		index=0;
		index_extra = 0;
		index_extra2 = 0;
		index_extra3 = 0;
		Best_score = ~0;
		update = true;
		temp.BestSubtreeRoot = -1;
		temp.BestNewLocation = NULL;
		queue.push_back(temp);
		tripFlag=true;

	}

	void run(ostream &os, const ReRoot reroot)
	{
		bool rerooting = (reroot == ALL);
//                createLeafMapping();


		// build list of unique leaf nodes
		vector<string> leafs;
		{
			set<string> templeafs;
			for (vector<GeneTreeRooted*>::iterator itr = genetree_rooted.begin(); itr != genetree_rooted.end(); itr++) {
				GeneTreeRooted &genetree_rooted = **itr;
				for (vector<NamedGeneNodeRooted*>::iterator itr = genetree_rooted.leafnodes.begin(); itr != genetree_rooted.leafnodes.end(); itr++) {
					NamedGeneNodeRooted &node = **itr;
					templeafs.insert(node.getName());
				}
			}
			for (set<string>::iterator itr = templeafs.begin(); itr != templeafs.end(); itr++)
				leafs.push_back(*itr);
		}




		// create and initialize the leafConstraints array
		int leafConstraint[leafs.size()];
		int leafSetSize = leafs.size();
		for(i= 0; i< leafSetSize; i++)
			leafConstraint[i] = -1;






		// sort constraints by constraint size
		multimap<int, vector<string>* > sorted_constraints;
		int last;
		for ( i=0, last=constraints.size(); i<last; i++) {
			vector<string> &c = constraints[i];
			sorted_constraints.insert(pair<int, vector<string>* >(c.size(), &c));
		}


		vector<vector<string>* > myconstraints;

		// traverse the sorted constraints, starting from the smallest constraint
		for (multimap<int, vector<string>* >::iterator itr=sorted_constraints.begin(), last=sorted_constraints.end(); itr!=last; itr++) {
			vector<string> &constraint = *itr->second;


			myconstraints.push_back(&constraint);


//			cout << "size: " << constraint.size() << " - ";
//			for (int i=0, last=constraint.size(); i<last; i++) {
//				cout << constraint[i] << " ";
//			}
//			cout << endl;
		}



	// Label each leaf node with its constraint (if any)



		int numlists =0 , listsize =0 ;
		numlists = myconstraints.size();

		// iterate through each string in constraints vector and search for it in the leaf set
		for(i= 0; i < numlists; i++)
		{
			listsize = myconstraints[i]->size();
			for (l = 0; l < listsize; l++)
			{
				for(k = 0; k < leafSetSize; k++)
				{

					// label the leafs in the leaf set with the relevant index from the constraints vector
					if(leafs[k] == (*myconstraints[i])[l] )
					{
						if( leafConstraint[k] == -1)
							leafConstraint[k] = i;
						// make sure the constraints are not nested
						else
						{
							EXCEPTION("The constraints are either contradictory or nested! This program can only handle non-nested constraints!" << endl << "Terminating...");
						}
					}
				}
			}
		}


		// create a boolean array of size numlists, to mark which lists are represented in the current tree

		bool leafpresent[numlists];

		for(i=0; i<numlists; i++)
			leafpresent[i] = false;


		// Create all possible triplets to get starting species tree
		vector<SpeciesNode*> nodes;
		i=0;
	       	l=1;
		k=2;


		msgout << "Leaf set size: " << leafSetSize <<endl;
		msgout <<  "Building Starting triplet ... ";

		index_extra = rand() % (leafSetSize-2);

		index_extra2 = index_extra+ 1 + (rand() % (leafSetSize - index_extra -2) );
		index_extra3 = index_extra2 + 1 + (rand() % (leafSetSize - index_extra2 -1));

		for (i= index_extra; i< index_extra + 1; i++)
		{
			for(l = index_extra2; l< index_extra2 + 1; l++)
			{
				for(k= index_extra3; k< index_extra3+ 1; k++)
				{
					NamedSpeciesNode *node1 = new NamedSpeciesNode(leafs[i]);
					node1->constraint = leafConstraint[i];
					if(leafConstraint[i] != -1)
						leafpresent[leafConstraint[i]] = true;

					NamedSpeciesNode *node2 = new NamedSpeciesNode(leafs[l]);
					node2->constraint = leafConstraint[l];
					if(leafConstraint[l] != -1)
						leafpresent[leafConstraint[l]] = true;

					NamedSpeciesNode *node3 = new NamedSpeciesNode(leafs[k]);
					node3->constraint = leafConstraint[k];

					currentconstraint = -1;
					if(leafConstraint[k] != -1)
					{
						if(leafpresent[leafConstraint[k]] == false)
							currentconstraint = -1;
						else
							currentconstraint = leafConstraint[k];
					}

					if(leafConstraint[k] != -1)
						leafpresent[leafConstraint[k]] = true;

					nodes.push_back(node1);
					nodes.push_back(node2);
					nodes.push_back(node3);
					speciestree->nodes.push_back(node1);
					speciestree->leafnodes.push_back(node1);
					speciestree->nodes.push_back(node2);
					speciestree->leafnodes.push_back(node2);
					speciestree->nodes.push_back(node3);
					speciestree->leafnodes.push_back(node3);

					SpeciesNode *node = new SpeciesNode();
					node->child(0) = nodes[0];
				       	node->child(1) = nodes[1];
					node->child(0)->parent() = node;
					node->child(1)->parent() = node;
					nodes.push_back(node);
					if(node->child(0)->constraint == node->child(1)->constraint)
						node->constraint = node->child(0)->constraint;
					else node->constraint = -1;
					speciestree->nodes.push_back(node);

					SpeciesNode *Secondnode = new SpeciesNode();
					Secondnode->child(0) = nodes[2];
				       	Secondnode->child(1) = nodes[3];
					Secondnode->child(0)->parent() = Secondnode;
					Secondnode->child(1)->parent() = Secondnode;
					nodes.push_back(Secondnode);
					if(Secondnode->child(0)->constraint == Secondnode->child(1)->constraint)
						Secondnode->constraint = Secondnode->child(0)->constraint;
					else Secondnode->constraint = -1;

					speciestree->nodes.push_back(Secondnode);

					speciestree->root = nodes[4];
					createLeafMapping();

	                                computeGeneDuplications(speciestree->nodes[2], rerooting);


					speciestree->nodes.clear();
					speciestree->leafnodes.clear();
					nodes.clear();
					if(leafConstraint[i] != -1)
						leafpresent[leafConstraint[i]] = false;
					if(leafConstraint[l] != -1)
						leafpresent[leafConstraint[l]] = false;
					if(leafConstraint[k] != -1)
						leafpresent[leafConstraint[k]] = false;
					delete node1;
					delete node2;
					delete node3;
					delete node;
					delete Secondnode;
				}
			}
		}


		index = rand()%trip_queue.size();

		i= trip_queue[index].i_val;
		l= trip_queue[index].l_val;
		k= trip_queue[index].k_val;
		trip_queue.clear();
		tripFlag = false;
		Best_score =  ~0;

		// Recreate the starting species tree

		NamedSpeciesNode *node1 = new NamedSpeciesNode(leafs[i]);
		NamedSpeciesNode *node2 = new NamedSpeciesNode(leafs[l]);
		NamedSpeciesNode *node3 = new NamedSpeciesNode(leafs[k]);
		nodes.push_back(node1);
		nodes.push_back(node2);
		nodes.push_back(node3);

		node1->constraint = leafConstraint[i];
		if(leafConstraint[i] != -1)
			leafpresent[leafConstraint[i]] = true;

		node2->constraint = leafConstraint[l];
		if(leafConstraint[l] != -1)
			leafpresent[leafConstraint[l]] = true;

		node3->constraint = leafConstraint[k];

		currentconstraint = -1;
		if(leafConstraint[k] != -1)
		{
			if(leafpresent[leafConstraint[k]] == false)
				currentconstraint = -1;
			else
				currentconstraint = leafConstraint[k];
		}
		if(leafConstraint[k] != -1)
			leafpresent[leafConstraint[k]] = true;

		speciestree->nodes.push_back(node1);
		speciestree->leafnodes.push_back(node1);
		speciestree->nodes.push_back(node2);
		speciestree->leafnodes.push_back(node2);
		speciestree->nodes.push_back(node3);
		speciestree->leafnodes.push_back(node3);

			// create internal nodes and set their constraint values

		SpeciesNode *node = new SpeciesNode();
		node->child(0) = nodes[0];
	       	node->child(1) = nodes[1];
		node->child(0)->parent() = node;
		node->child(1)->parent() = node;
		if(node->child(0)->constraint == node->child(1)->constraint)
			node->constraint = node->child(0)->constraint;
		else node->constraint = -1;


		nodes.push_back(node);
		speciestree->nodes.push_back(node);
		SpeciesNode *Secondnode = new SpeciesNode();
		Secondnode->child(0) = nodes[2];
	       	Secondnode->child(1) = nodes[3];
		Secondnode->child(0)->parent() = Secondnode;
		Secondnode->child(1)->parent() = Secondnode;
		if(Secondnode->child(0)->constraint == Secondnode->child(1)->constraint)
			Secondnode->constraint = Secondnode->child(0)->constraint;
		else Secondnode->constraint = -1;

		nodes.push_back(Secondnode);
		speciestree->nodes.push_back(Secondnode);
		speciestree->root = nodes[4];
		j=k;
       		createLeafMapping();
		computeGeneDuplications(speciestree->nodes[2], rerooting);




//		speciestree.nodes.clear();
//		speciestree.leafnodes.clear();
//		nodes.clear();
//		delete node1;
//		delete node2;
//		delete node3;
//		delete node;
//		delete Secondnode;

//speciestree->checkStructure();

		index = rand()%queue.size();
		old.BestNewLocation = queue[index].BestNewLocation;
		speciestree->moveSubtree(speciestree->nodes[2], old.BestNewLocation);


		if(old.BestNewLocation->constraint == leafConstraint[k])
			old.BestNewLocation->parent()->constraint = leafConstraint[k];
		else old.BestNewLocation->parent()->constraint = -1;


//cout << *speciestree << endl;
//speciestree->tree2newick(os); os << endl;
//cout << node1->constraint << node2->constraint << node3->constraint << node->constraint << Secondnode->constraint << endl;


		// remove the leafs that have already been added



//for (int count =0; count < leafs.size(); count++)
//cout << leafs[count] << " " << leafConstraint[count] << endl;
//cout << endl << i << " " << l<< " "<< k;

		leafs[k].swap(leafs[leafs.size() -1]);
		leafConstraint[k] = leafConstraint[leafs.size()-1];
		leafs.pop_back();

		leafs[l].swap(leafs[leafs.size() -1]);
		leafConstraint[l] = leafConstraint[leafs.size()-1];
		leafs.pop_back();

		leafs[i].swap(leafs[leafs.size() -1]);
		leafConstraint[i] = leafConstraint[leafs.size()-1];
		leafs.pop_back();
//		leafs.erase(i);
//		leafs.erase(l);
//		leafs.erase(k);



	// Start the additive procedure
		leafSetSize = leafs.size();
		queue.clear();
		Best_score =  ~0;

		msgout << "\rCurrent species tree leaf set size: " << speciestree->leafnodes.size() << "        ";
	while(leafSetSize > 0)
	{



		int nodeSetSize = nodes.size();
		index_extra = rand()%leafSetSize;
		for (i=index_extra; i< index_extra + 1; i++)
		{

			//build temporary tree by adding one leaf and check its scores
			NamedSpeciesNode *node1 = new NamedSpeciesNode(leafs[i]);

			node1->constraint = leafConstraint[i];
			currentconstraint = -1;
			if(leafConstraint[i] != -1)
			{
				if(leafpresent[leafConstraint[i]] == false)
					currentconstraint = -1;
				else
					currentconstraint = leafConstraint[i];
			}





			if(leafConstraint[i] != -1)
				leafpresent[leafConstraint[i]] = true;

			nodes.push_back(node1);
			speciestree->nodes.push_back(node1);
			speciestree->leafnodes.push_back(node1);

			SpeciesNode *node = new SpeciesNode();
			node->child(0) = nodes[nodeSetSize];
		       	node->child(1) = speciestree->root;
			node->child(0)->parent() = node;
			node->child(1)->parent() = node;

			if(node->child(0)->constraint == node->child(1)->constraint)
				node->constraint = node->child(0)->constraint;
			else node->constraint = -1;


			nodes.push_back(node);
			speciestree->nodes.push_back(node);
			speciestree->root = nodes[nodeSetSize + 1];
			j= i;

			createLeafMapping();
//cout << leafs[i] << " " << leafConstraint[i] << " current constraint is: " << currentconstraint <<endl;
	        	computeGeneDuplications(speciestree->nodes[nodeSetSize], rerooting);

			//undo the leaf addition step to try out other leaf additions
			nodes.pop_back();
			nodes.pop_back();
			if(leafConstraint[i] != -1)
				if(currentconstraint == -1)
					leafpresent[leafConstraint[i]] = false;
			speciestree->leafnodes.pop_back();
			speciestree->nodes.pop_back();
			speciestree->nodes.pop_back();
			speciestree->root = node->child(1);
			delete node1;
			delete node;


		}




			// reconstruct the optimal species tree from queue data


			index = rand()%queue.size();
			NamedSpeciesNode *node1 = new NamedSpeciesNode(leafs[queue[index].BestSubtreeRoot]);

			node1->constraint = leafConstraint[queue[index].BestSubtreeRoot];

			currentconstraint = -1;
			if(leafConstraint[queue[index].BestSubtreeRoot] != -1)
			{
				if(leafpresent[leafConstraint[queue[index].BestSubtreeRoot]] == false)
					currentconstraint = -1;
				else
					currentconstraint = leafConstraint[queue[index].BestSubtreeRoot];
			}


			if(leafConstraint[queue[index].BestSubtreeRoot] != -1)
				leafpresent[leafConstraint[queue[index].BestSubtreeRoot]] = true;

			nodes.push_back(node1);
			speciestree->nodes.push_back(node1);
			speciestree->leafnodes.push_back(node1);
			SpeciesNode *node = new SpeciesNode();
			node->child(0) = nodes[nodeSetSize];
		       	node->child(1) = speciestree->root;
			node->child(0)->parent() = node;
			node->child(1)->parent() = node;

			if(node->child(0)->constraint == node->child(1)->constraint)
				node->constraint = node->child(0)->constraint;
			else node->constraint = -1;


			nodes.push_back(node);
			speciestree->nodes.push_back(node);
			speciestree->root = nodes[nodeSetSize + 1];

			createLeafMapping();
			speciestree->moveSubtree(speciestree->nodes[nodeSetSize], queue[index].BestNewLocation);

			if(queue[index].BestNewLocation->constraint == leafConstraint[queue[index].BestSubtreeRoot])
				queue[index].BestNewLocation->parent()->constraint = leafConstraint[queue[index].BestSubtreeRoot];
			else queue[index].BestNewLocation->parent()->constraint = -1;








			leafs[queue[index].BestSubtreeRoot].swap(leafs[leafs.size() -1]);
			leafConstraint[queue[index].BestSubtreeRoot] = leafConstraint[leafs.size()-1];
			leafs.pop_back();

			leafSetSize = leafs.size();

			queue.clear();
	 		msgout << "\rCurrent species tree leaf set size: " << speciestree->leafnodes.size() << "        ";
			if ( leafSetSize == 0)
				msgout << endl << "Initial species tree with " << speciestree->leafnodes.size() << " leaves computed. Starting local search..." << endl << endl;		
//				cout << endl << "Total number of gene duplications and losses for initial species tree: "<< Best_score << endl << endl;
			Best_score =  ~0;
	}



		speciestree->tree2newick(os); os << endl;
	}

	int currentScore() {
		return Best_score;
	}

};


#undef TreeSet
#undef SpeciesNode
#undef NamedSpeciesNode
#undef GeneNodeUnrooted
#undef NamedGeneNodeUnrooted
#undef GeneTreeUnrooted
#undef GeneNodeRooted
#undef NamedGeneNodeRooted
#undef GeneTreeRooted

#endif
