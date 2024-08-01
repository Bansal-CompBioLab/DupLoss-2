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

#ifndef REPORT_TREESET_H
#define REPORT_TREESET_H

#include "rmq.h"
#include "rmq.c"

// ------------------------------------------------------------------------------------------------------------------
// primary mapping (unrooted)
template<class T>
class PrimaryMappingUnrooted {
protected:
	// mapping
	T* mapping[3];

public:
	PrimaryMappingUnrooted() {
		resetMapping();
	}

	// get mapping in direction i
	inline T*& getMapping(const int &i) {
		return mapping[i];
	}

	// set mapping in direction i
	template<class T2>
	inline void setMapping(int &i, T2* &node) {
		mapping[i] = node;
	}

	// reset mappings for all directions
	inline void resetMapping() {
		for (int i=0; i<3; i++) mapping[i] = NULL;
	}

	template<class T2>
	friend ostream & operator << (ostream & os, PrimaryMappingUnrooted<T2> & m);
};

// output the node into a string stream
template<class T>
ostream & operator << (ostream & os, PrimaryMappingUnrooted<T> & m) {
	os << "1st-map:{" << m.mapping[0] << ',' << m.mapping[1] << ','  << m.mapping[2] << '}';
	return os;
}

// ------------------------------------------------------------------------------------------------------------------
// primary mapping (rooted)
template<class T>
class PrimaryMappingRooted {
protected:
	// mapping
	T* mapping;

public:
	PrimaryMappingRooted() {
		resetMapping();
	}

	// get mapping in direction i
	inline T*& getMapping() {
		return mapping;
	}

	// set mapping in direction i
	template<class T2>
	inline void setMapping(T2* &node) {
		mapping = node;
	}

	// reset mappings for all directions
	inline void resetMapping() {
		mapping = NULL;
	}

	template<class T2>
	friend ostream & operator << (ostream & os, PrimaryMappingRooted<T2> & m);
};

// output the node into a string stream
template<class T>
ostream & operator << (ostream & os, PrimaryMappingRooted<T> & m) {
	os << "1st-map:" << m.mapping;
	return os;
}

// ------------------------------------------------------------------------------------------------------------------
// secondary mapping
template<class T>
class SecondaryMapping {
public:
	// mapping
	T* secondarymapping;

	SecondaryMapping() {
		secondarymapping = NULL;
	}

	// check if contained in Gamma-Tree
	inline bool belongs2GammaTree() {
		return secondarymapping != NULL;
	}

	template<class T2>
	friend ostream & operator << (ostream & os, SecondaryMapping<T2> & m);
};

// output the node into a string stream
template<class T>
ostream & operator << (ostream & os, SecondaryMapping<T> & m) {
	os << "2nd-map:" << m.secondarymapping;
	return os;
}

// ------------------------------------------------------------------------------------------------------------------
// order number
class Order {
public:
	unsigned int no, begin, end;

	Order() {
		no = 0;
		begin = 0;
		end = 0;
	}

	friend ostream & operator << (ostream & os, Order & m);
};

// output the node into a string stream
ostream & operator << (ostream & os, Order & m) {
	os << "order:(" << m.begin << '<' << m.no << '<' << m.end << ')';
	return os;
}

// ------------------------------------------------------------------------------------------------------------------
// order number
class GeneDuplication {
public:
	unsigned int gain, lost[2];
	unsigned int genedup, tempgenedup;

	GeneDuplication() {
		gain = 0;
		for (int i=0; i<2; i++) lost[i] = 0;
		genedup = 0;
	}

	friend ostream & operator << (ostream & os, GeneDuplication & m);
};

// output the node into a string stream
ostream & operator << (ostream & os, GeneDuplication & m) {
	os << "genedup:" << m.genedup << '(' << m.lost[0] << '[' << m.gain << ']' << m.lost[1] << ')';
	return os;
}



class GeneLoss {
public:
	int subtreeSize, nodeDepth;
	
	//subtreeSize counts the number leaves in the subtree that have mappings from the gene tree
	// nodeDepth is 1 for root. 
	bool isRelevant;

	GeneLoss() {
		subtreeSize = 0;
//		lossScore = 0;
		nodeDepth = 0;
		isRelevant = false;
//		LossParent = NULL;
//		LossChild[0] = NULL;
//		LossChild[1] = NULL;
	}

	friend ostream & operator << (ostream & os, GeneLoss & m);
};

// output the node into a string stream
ostream & operator << (ostream & os, GeneLoss & m) {
	//os << "Losses:" << m.subtreeSize << " " << m.isRelevant << " " << m.lossCounter1 << " " << m.nodeDepth << " " << m.lossScore;
	os << "Losses:" << m.subtreeSize << " " << m.isRelevant;
	return os;
}


// ==================================================================================================================
// ==================================================================================================================
// ==================================================================================================================

// ------------------------------------------------------------------------------------------------------------------
class SepciesNode;
class NamedSepciesNode;
class GeneNodeUnrooted;
class NamedGeneNodeUnrooted;
class SpeciesTree;
class GeneTree;

// ------------------------------------------------------------------------------------------------------------------
// species node
class SpeciesNode : public TreeNodeRooted<SpeciesNode>, public Order, public GeneDuplication, public GeneLoss {
public:
	int idx;

	SpeciesNode(SpeciesNode  *parent = NULL) : TreeNodeRooted<SpeciesNode>(parent) {}

	friend ostream & operator << (ostream & os, SpeciesNode & m);
};

// output a node into a string stream
ostream & operator << (ostream & os, SpeciesNode & m) {
	os << (TreeNodeRooted<SpeciesNode> &)m << ", "
	   << "idx:" << m.idx << ", "
	   << (Order &)m << ", "
	   << (GeneDuplication &)m;
	return os;
}

// ------------------------------------------------------------------------------------------------------------------
// species node with name
class NamedSpeciesNode : public SpeciesNode, public NodeID {
public:
	NamedSpeciesNode(const string &id, SpeciesNode *parent = NULL) : SpeciesNode(parent), NodeID(id) {}

	friend ostream & operator << (ostream & os, NamedSpeciesNode & m);
};

// output a node into a string stream
ostream & operator << (ostream & os, NamedSpeciesNode & m) {
	os << (NodeID &)m << ", "
	   << (SpeciesNode &)m;
	return os;
}

// ------------------------------------------------------------------------------------------------------------------
// gene node (rooted)
class GeneNodeRooted : public TreeNodeRooted<GeneNodeRooted>, public PrimaryMappingRooted<SpeciesNode>, public SecondaryMapping<SpeciesNode> {
public:
	using PrimaryMappingRooted<SpeciesNode>::getMapping;

	GeneNodeRooted(GeneNodeRooted *parent = NULL) : TreeNodeRooted<GeneNodeRooted>(parent) {}

	// get LCA mapping in rooted direction
	inline SpeciesNode*& getMapping() {
		return mapping;
	}

	// set LCA mapping in rooted direction
	template<class T2>
	inline void setMapping(T2* &node) {
		mapping = node;
	}

	friend ostream & operator << (ostream & os, GeneNodeRooted & m);
};

// output the node into a string stream
ostream & operator << (ostream & os, GeneNodeRooted & m) {
	os << (TreeNodeRooted<GeneNodeRooted> &)m << ", "
	   << (PrimaryMappingRooted<SpeciesNode> &)m << ", "
	   << (SecondaryMapping<SpeciesNode> &)m;
	return os;
}

// ------------------------------------------------------------------------------------------------------------------
// gene node with name
class NamedGeneNodeRooted : public GeneNodeRooted, public NodeID {
public:
	NamedGeneNodeRooted(const string &id, GeneNodeRooted *parent = NULL) : GeneNodeRooted(parent), NodeID(id) {}

	friend ostream & operator << (ostream & os, NamedGeneNodeRooted & m);
};

// output a node into a string stream
ostream & operator << (ostream & os, NamedGeneNodeRooted & m) {
	os << (NodeID &)m << ", " << (GeneNodeRooted &)m;
	return os;
}

// ------------------------------------------------------------------------------------------------------------------
// gene node (unrooted)
class GeneNodeUnrooted : public TreeNodeUnrooted<GeneNodeUnrooted>, public PrimaryMappingUnrooted<SpeciesNode>, public SecondaryMapping<SpeciesNode> {
public:
	using PrimaryMappingUnrooted<SpeciesNode>::getMapping;

	GeneNodeUnrooted(GeneNodeUnrooted *parent = NULL) : TreeNodeUnrooted<GeneNodeUnrooted>(parent) {}

	// get LCA mapping in rooted direction
	inline SpeciesNode*& getMapping() {
		return mapping[parentno];
	}

	// set LCA mapping in rooted direction
	template<class T2>
	inline void setMapping(T2* &node) {
		mapping[parentno] = node;
	}

	// return the LCA mappings (unrooted)
	inline SpeciesNode* &getMappingDirected(GeneNodeUnrooted* &parent) {
		for (int i=0; i<3; i++) {
			if (relative[i] == parent) {
				return mapping[i];
			}
		}
		EXCEPTION("direction not found in getMappingDirected");
		return mapping[0];
	}

	// set the LCA mapping (unrooted)
	template<class T2>
	inline void setMappingDirected(GeneNodeUnrooted* &parent, T2* &node) {
		for (int i=0; i<3; i++) {
			if (relative[i] == parent) {
				mapping[i] = node;
				return;
			}
		}
		EXCEPTION("direction not found in setMappingDirected");
	}

	friend ostream & operator << (ostream & os, GeneNodeUnrooted & m);
};

// output the node into a string stream
ostream & operator << (ostream & os, GeneNodeUnrooted & m) {
	os << (TreeNodeUnrooted<GeneNodeUnrooted> &)m << ", "
	   << (PrimaryMappingUnrooted<SpeciesNode> &)m << ", "
	   << (SecondaryMapping<SpeciesNode> &)m;
	return os;
}

// ------------------------------------------------------------------------------------------------------------------
// gene node with name
class NamedGeneNodeUnrooted : public GeneNodeUnrooted, public NodeID {
public:
	NamedGeneNodeUnrooted(const string &id, GeneNodeUnrooted *parent = NULL) : GeneNodeUnrooted(parent), NodeID(id) {}

	friend ostream & operator << (ostream & os, NamedGeneNodeUnrooted & m);
};

// output a node into a string stream
ostream & operator << (ostream & os, NamedGeneNodeUnrooted & m) {
	os << (NodeID &)m << ", " << (GeneNodeUnrooted &)m;
	return os;
}

// ==================================================================================================================
// ==================================================================================================================
// ==================================================================================================================

// ------------------------------------------------------------------------------------------------------------------
class SpeciesTree : public TreeIO<Tree<SpeciesNode, NamedSpeciesNode> > {
public:
	SpeciesTree() {
		R = NULL; E = NULL; L = NULL;
	}

	virtual ~SpeciesTree() {
		delete [] L, E, R;
	}

	// set the gene duplication triple to 0
	void resetGeneDubTriple() {
		for (vector<SpeciesNode*>::iterator itr = nodes.begin(); itr != nodes.end(); itr++) {
			SpeciesNode &t = **itr;
			t.gain = 0;
			t.lost[0] = 0;
			t.lost[1] = 0;
		}
	}

	// assigns each node an index
	void assignIndex() {
		for (int i = 0; i < nodes.size(); i++) nodes[i]->idx = i;
	}

	// node numbering + range of subtree
	inline void establishOrder() {
		int pos = 0;
		establishOrderDFS(root, pos);
	}
	inline void establishOrderDFS(SpeciesNode *node, int &pos) {
		node->begin = pos;
		if (node->child(0) != NULL) {
			establishOrderDFS(node->child(0), pos);
			pos++;
		}
		node->no = pos;
		if (node->child(1) != NULL) {
			pos++;
			establishOrderDFS(node->child(1), pos);
		}
		node->end = pos;
	}

	// preprocess the LCA computation (build RMQ)
	INT *R; // first occurences in sequence
	VAL *E, *L; // sequence - E:nodes, L:levels
	struct rmqinfo *ri; // lookup table
	void preprocessLCA() {
		const int n = nodes.size();
		if ((E == NULL) || (L == NULL) || (R == NULL)) {
			assignIndex();
			delete [] R; R = new INT[n];
			delete [] E; E = new VAL[n];
			delete [] L; L = new VAL[n];
		}
		VAL *e = &E[0], *l = &L[0];
		createInorderSequence(root, e, l);
		for (int i=n-1; i>=0; i--) R[E[i]] = i;
		ri = rm_query_preprocess(L, n);
	}

	// clean up of LCA computation
	void postprocessLCA() {
		rm_free(ri);
	}

	// return the LCA of 2 species nodes
	SpeciesNode*& getLCA(SpeciesNode* &u, SpeciesNode* &v) {
		const INT uidx = R[u->idx];
		const INT vidx = R[v->idx];
		int y;
		if (uidx < vidx) y = rm_query(ri, uidx, vidx);
		else y = rm_query(ri, vidx, uidx);
		return nodes[E[y]];
	}

	// create sequences of nodes and levels of the species tree in an Euler Tour
	void createInorderSequence(SpeciesNode *&node, VAL *&e, VAL *&l, int lvl = 0) {
		if (node == NULL) return;
		bool const isLeaf = node->isLeaf();
		createInorderSequence(node->child(0), e, l, lvl+1);
		*l = lvl;l++;
		*e = node->idx; e++;
		createInorderSequence(node->child(1), e, l, lvl+1);
	}
};

// ------------------------------------------------------------------------------------------------------------------
class GeneTree {
public:
	unsigned int score, lossScore;
	GeneTree() {
		score = 0;
		lossScore = 0;
	}
};

// ------------------------------------------------------------------------------------------------------------------
// a rooted binary gene tree
class GeneTreeRooted : public GeneTree, public TreeIO<Tree<GeneNodeRooted, NamedGeneNodeRooted> > {
public:
	inline bool isRooted() {
		return true;
	}
};

// a quasi unrooted binary gene tree
// (rooted gene tree, but it allows rerooting)
class GeneTreeUnrooted : public GeneTree, public TreeIO<Tree<GeneNodeUnrooted, NamedGeneNodeUnrooted> > {
public:
	inline bool isRooted() {
		return false;
	}

	// move the root into a different edge
	inline void reroot(GeneNodeUnrooted *u, GeneNodeUnrooted *v) {
		disconnectRoot();
		connectRoot(u, v);
		rerootDFS(root->child(0), root);
		rerootDFS(root->child(1), root);
	}
	inline void rerootDFS(GeneNodeUnrooted *&node, GeneNodeUnrooted *&parent) {
		if (node == NULL) return;
		if (node->parent() == parent) return;
		node->assignParent(parent);
		rerootDFS(node->child(0), node);
		rerootDFS(node->child(1), node);
	}
};

// ==================================================================================================================
// ==================================================================================================================
// ==================================================================================================================

// ------------------------------------------------------------------------------------------------------------------
// the container for a whole set of trees (1 species tree and several gene trees)
class TreeSet {
public:
	// data container for the species tree
	SpeciesTree *speciestree;

	// data container for the gene trees
	vector<GeneTreeRooted*> genetree_rooted;
	vector<GeneTreeUnrooted*> genetree_unrooted;

	TreeSet() {
		#ifdef DEBUG
		cout << "TreeSet created" << endl;
		#endif
		speciestree = NULL;
	}

	virtual ~TreeSet() {
		for(vector<GeneTreeRooted*>::iterator itr=genetree_rooted.begin(); itr!=genetree_rooted.end(); itr++) delete *itr;
		for(vector<GeneTreeUnrooted*>::iterator itr=genetree_unrooted.begin(); itr!=genetree_unrooted.end(); itr++) delete *itr;
		if (speciestree != NULL) delete speciestree;
		#ifdef DEBUG
		cout << "TreeSet destroyed" << endl;
		#endif
	}

	// ------------------------------------------------------------------------------------------------------
	// read all trees from the input
	void readTrees(Input &input) {
		// read the species tree
		speciestree = new SpeciesTree;
		if (speciestree->stream2tree(input)) {
// 			msgout << "input species tree of " << speciestree->leafnodes.size() << " taxa" << endl;
		} else EXCEPTION("missing input for species tree " << input.getLastPos());

		// read all gene trees
		for (;;) {
			string str;
			bool rooted = true;
			while (!(str=input.getComment()).empty()) {
				if (str == "[&U]") rooted = false;
				if (str == "[&R]") rooted = true;
			}

			if (rooted) {
				GeneTreeRooted *tree = new GeneTreeRooted;
				if (tree->stream2tree(input)) {
					genetree_rooted.push_back(tree);
// 					msgout << "rooted input gene tree of " << tree->leafnodes.size() << " taxa" << endl;
				} else {
					delete tree;
					break;
				}
			} else {
				GeneTreeUnrooted *tree = new GeneTreeUnrooted;
				if (tree->stream2tree(input)) {
					genetree_unrooted.push_back(tree);
// 					msgout << "unrooted input gene tree of " << tree->leafnodes.size() << " taxa" << endl;
				} else {
					delete tree;
					break;
				}
			}
		}
		if (genetree_rooted.size() + genetree_unrooted.size() < 1)
			EXCEPTION("missing input for gene trees " << input.getLastPos()); // ERROR no gene trees

// 		msgout << genetree_rooted.size() + genetree_unrooted.size() << " input gene trees total" << endl;
	}

	// ------------------------------------------------------------------------------------------------------
	// establish the initial primary mapping
	// mappings between the leaf nodes of gene trees and the species tree; according to their names
	void createLeafMapping() {
		// get the leaf nodes of the species tree
		// the leaf nodes of the species tree have to be a full set of all species and be unique
		vector<NamedSpeciesNode*> &s = speciestree->leafnodes;
		// sort these node by their name
		map<const string*, NamedSpeciesNode*> sm;
		for(vector<NamedSpeciesNode*>::iterator itr=s.begin(); itr!=s.end(); itr++) {
			NamedSpeciesNode &species = **itr;
			if (sm.find(&species.name->first) != sm.end()) // leafe node is not unique
				EXCEPTION('"' << species.name->first << "\" is not a unique species");
			sm[&species.name->first] = &species;
		}

		// create the leaf node mapping for all rooted gene trees
		for(int i=0; i<genetree_rooted.size(); i++) {
			// get the leaf nodes of the gene tree
			vector<NamedGeneNodeRooted*> &g = genetree_rooted[i]->leafnodes;
			// create links for each leaf node to the species tree
			for(vector<NamedGeneNodeRooted*>::iterator itr=g.begin(); itr!=g.end(); itr++) {
				NamedGeneNodeRooted *genenode = *itr;
				map<const string*, NamedSpeciesNode*>::iterator smt = sm.find(&genenode->name->first);
				if (smt == sm.end()) // leaf node of genetree is not in the species tree!
					EXCEPTION("species \"" << genenode->name->first << "\" is not covered in the species tree");
				// establish mapping
				genenode->setMapping(smt->second);
			}
		}

		// create the leaf node mapping for all unrooted gene trees
		for(int i=0; i<genetree_unrooted.size(); i++) {
			// get the leaf nodes of the gene tree
			vector<NamedGeneNodeUnrooted*> &g = genetree_unrooted[i]->leafnodes;
			// create links for each leaf node to the species tree
			for(vector<NamedGeneNodeUnrooted*>::iterator itr=g.begin(); itr!=g.end(); itr++) {
				NamedGeneNodeUnrooted *genenode = *itr;
				map<const string*, NamedSpeciesNode*>::iterator smt = sm.find(&genenode->name->first);
				if (smt == sm.end()) // leaf node of genetree is not in the species tree!
					EXCEPTION("species \"" << genenode->name->first << "\" is not covered in the species tree");
				// establish mapping
				genenode->setMapping(smt->second);
			}
		}

		// check for complete leaf checking (all species nodes have to map to at least 1 node in any gene tree)
		for(int i=0; i<genetree_rooted.size(); i++) {
			vector<NamedGeneNodeRooted*> &g = genetree_rooted[i]->leafnodes;
			for(vector<NamedGeneNodeRooted*>::iterator itr=g.begin(); itr!=g.end(); itr++) {
				NamedGeneNodeRooted *genenode = *itr;
				sm.erase(&genenode->name->first);
			}
		}
		for(int i=0; i<genetree_unrooted.size(); i++) {
			vector<NamedGeneNodeUnrooted*> &g = genetree_unrooted[i]->leafnodes;
			for(vector<NamedGeneNodeUnrooted*>::iterator itr=g.begin(); itr!=g.end(); itr++) {
				NamedGeneNodeUnrooted *genenode = *itr;
				sm.erase(&genenode->name->first);
			}
		}
		if (!sm.empty()) {
			for(map<const string*, NamedSpeciesNode*>::iterator itr=sm.begin(); itr!=sm.end(); itr++) {
				const string *name = itr->first;
				NamedSpeciesNode *&node = itr->second;
				WARNING("species \"" << *name << "\" is not covered by any gene tree - species ignored");
				speciestree->deleteLeafNode(node);
			}
		}
	}

	// ------------------------------------------------------------------------------------------------------
	// establish the LCA mapping between one gene tree and the species tree
	void createPrimaryMapping(GeneTreeRooted &tree) {
		// reset mappings for all internal nodes
		for (vector<GeneNodeRooted*>::iterator itr=tree.nodes.begin(); itr!=tree.nodes.end(); itr++) {
			GeneNodeRooted &t = **itr;
			if (!t.isLeaf()) t.resetMapping();
		}
		// establish LCA mapping
		getLCA(tree.root, tree.root->parent(), speciestree->E, speciestree->R, speciestree->ri);
	}
	void createPrimaryMapping(GeneTreeUnrooted &tree) {
		// reset mappings for all internal nodes
		for (vector<GeneNodeUnrooted*>::iterator itr=tree.nodes.begin(); itr!=tree.nodes.end(); itr++) {
			GeneNodeUnrooted &t = **itr;
			if (!t.isLeaf()) t.resetMapping();
		}
		// establish LCA mapping
		getLCA(tree.root, tree.root->parent(), speciestree->E, speciestree->R, speciestree->ri);
	}

	// establish the LCA mapping between one unrooted gene tree and the species tree
	void createPrimaryMappingUnrooted(GeneTreeRooted &tree) {
		// reset mappings for all internal nodes
		for (vector<GeneNodeRooted*>::iterator itr=tree.nodes.begin(); itr!=tree.nodes.end(); itr++) {
			GeneNodeRooted &t = **itr;
			if (!t.isLeaf()) t.resetMapping();
		}
		// establish LCA mapping (unrooted)
		getLCA(tree.root, tree.root->parent(), speciestree->E, speciestree->R, speciestree->ri);
		for (vector<NamedGeneNodeRooted*>::iterator itr=tree.leafnodes.begin(); itr!=tree.leafnodes.end(); itr++) {
			GeneNodeRooted *node = (GeneNodeRooted*)*itr;
			getLCA(node->parent(), node, speciestree->E, speciestree->R, speciestree->ri);
		}
	}
	void createPrimaryMappingUnrooted(GeneTreeUnrooted &tree) {
		// reset mappings for all internal nodes
		for (vector<GeneNodeUnrooted*>::iterator itr=tree.nodes.begin(); itr!=tree.nodes.end(); itr++) {
			GeneNodeUnrooted &t = **itr;
			if (!t.isLeaf()) t.resetMapping();
		}
		// establish LCA mapping (unrooted)
		getLCA(tree.root, tree.root->parent(), speciestree->E, speciestree->R, speciestree->ri);
		for (vector<NamedGeneNodeUnrooted*>::iterator itr=tree.leafnodes.begin(); itr!=tree.leafnodes.end(); itr++) {
			GeneNodeUnrooted *node = (GeneNodeUnrooted*)*itr;
			getLCA(node->parent(), node, speciestree->E, speciestree->R, speciestree->ri);
		}
	}

	// returns the LCA mapping (if necessary establish LCA mapping)
	SpeciesNode* &getLCA(GeneNodeRooted* &node, GeneNodeRooted* &node2, VAL E[], INT R[], struct rmqinfo *&ri) {
		#ifdef DEBUG
		if (node == NULL) EXCEPTION("node is NULL in getLCA" << endl);
		#endif
		SpeciesNode *&mapping = node->getMapping();
		if (mapping != NULL) return mapping;
		GeneNodeRooted *gn[2];
		node->getChildren(gn);
		if (gn[0] == NULL) {
			node->setMapping(getLCA(gn[1], node, E, R, ri));
		} else
		if (gn[1] == NULL) {
			node->setMapping(getLCA(gn[0], node, E, R, ri));
		} else {
			SpeciesNode *&u = getLCA(gn[0], node, E, R, ri);
			SpeciesNode *&v = getLCA(gn[1], node, E, R, ri);
			SpeciesNode *&p = speciestree->getLCA(u, v);
			node->setMapping(p);
		}
		return node->getMapping();
	}
	SpeciesNode* &getLCA(GeneNodeUnrooted* &node, GeneNodeUnrooted* &node2, VAL E[], INT R[], struct rmqinfo *&ri) {
		#ifdef DEBUG
		if (node == NULL) EXCEPTION("node is NULL in unrooted getLCA" << endl);
		#endif
		SpeciesNode *&mapping = node->getMappingDirected(node2);
		if (mapping != NULL) return mapping;
		GeneNodeUnrooted *gn[2];
		node->getChildrenDirected(node2, gn);
		if (gn[0] == NULL) {
			node->setMappingDirected(node2, getLCA(gn[1], node, E, R, ri));
		} else
		if (gn[1] == NULL) {
			node->setMappingDirected(node2, getLCA(gn[0], node, E, R, ri));
		} else {
			SpeciesNode *&u = getLCA(gn[0], node, E, R, ri);
			SpeciesNode *&v = getLCA(gn[1], node, E, R, ri);
			SpeciesNode *&p = speciestree->getLCA(u, v);
			node->setMappingDirected(node2, p);
		}
		return node->getMappingDirected(node2);
	}

	// ------------------------------------------------------------------------------------------------------
	// establish the secondary mapping
	template<class GeneTree>
	void createSecondaryMapping(GeneTree &tree, SpeciesNode *&subtreeroot) {
		// establish LCA mapping
		if (tree.root->getMapping() == speciestree->root) {
			get2ndLCA(tree.root, *subtreeroot, speciestree->E, speciestree->R, speciestree->ri);
		}
	}

	// returns the 2nd LCA mapping (if necessary establish mapping)
	template<class GeneNode>
	SpeciesNode* &get2ndLCA(GeneNode* &node, SpeciesNode &subtreeroot, VAL E[], INT R[], struct rmqinfo *&ri) {
		#ifdef DEBUG
		if (node == NULL) EXCEPTION("node is NULL in get2ndLCA" << endl);
		#endif
		SpeciesNode *&mapping2 = node->secondarymapping;
		if (mapping2 != NULL) return mapping2;

		SpeciesNode *&mapping1 = node->getMapping();
		if (mapping1 != speciestree->root) {
			const int &begin = subtreeroot.begin;
			const int &end = subtreeroot.end;
			const int &no = mapping1->no;
			if ((begin>no) || (no>end)) {
				mapping2 = mapping1;
				addSupportNode(node);
			}
		} else {
			GeneNode *gn[2];
			node->getChildren(gn);
			#ifdef DEBUG
			for (int i = 0; i < 2; i++) {
				if (node->child(i) == NULL) {
					EXCEPTION("child node->child(" << i << ") is NULL in get2ndLCA" << endl);
				}
			}
			#endif
			SpeciesNode *&u = get2ndLCA(node->child(0), subtreeroot, E, R, ri);
			SpeciesNode *&v = get2ndLCA(node->child(1), subtreeroot, E, R, ri);
			if (u == NULL) mapping2 = v;
			else
			if (v == NULL) mapping2 = u;
			else
			mapping2 = speciestree->getLCA(u, v);
		}
		return mapping2;
	}

	// remove mapping to ghost nodes in the species tree
	vector<GeneNodeRooted*> supportnode_rooted;
	vector<GeneNodeUnrooted*> supportnode_unrooted;
	inline void addSupportNode(GeneNodeRooted* &node) {
		supportnode_rooted.push_back(node);
	}
	inline void addSupportNode(GeneNodeUnrooted* &node) {
		supportnode_unrooted.push_back(node);
	}
	inline void removeSecondaryMapping(GeneTreeRooted &tree) {
		resetSecondaryMappingDFS(tree.root);
		supportnode_rooted.clear();
	}
	inline void removeSecondaryMapping(GeneTreeUnrooted &tree) {
		resetSecondaryMappingDFS(tree.root);
		supportnode_unrooted.clear();
	}
	template<class GeneNode>
	inline void resetSecondaryMappingDFS(GeneNode *&node) {
		if ((node == NULL) || (node->secondarymapping == NULL)) return;
		node->secondarymapping = NULL;
		resetSecondaryMappingDFS(node->child(0));
		resetSecondaryMappingDFS(node->child(1));
	}

	// ------------------------------------------------------------------------------------------------------
	// compute the gene duplications
	inline void computeGeneDuplicationsTriple() {
		speciestree->resetGeneDubTriple();
		computeGeneDuplicationsTripleAdd();
	}

	inline void computeGeneDuplicationsTripleAdd() {
		for (vector<GeneNodeRooted*>::iterator itr=supportnode_rooted.begin(); itr!=supportnode_rooted.end(); itr++) {
			GeneNodeRooted &node = **itr;
			GeneNodeRooted &parent = *node.parent();
			GeneNodeRooted &sibling = *node.getSibling();
			if (!sibling.belongs2GammaTree()) {
				// both nodes are in the ghost-node-tree
				node.secondarymapping->gain++;
			} else {
				// only one node is in the ghost-node-tree
				const int umapid = parent.secondarymapping->no;
				const int vmapid = sibling.secondarymapping->no;
				const int omapid = node.secondarymapping->no;
				if ((vmapid < umapid) && (umapid < omapid)) {
					parent.secondarymapping->lost[0]++;
				} else
				if ((omapid < umapid) && (umapid < vmapid)) {
					parent.secondarymapping->lost[1]++;
				}
			}
		}
		for (vector<GeneNodeUnrooted*>::iterator itr=supportnode_unrooted.begin(); itr!=supportnode_unrooted.end(); itr++) {
			GeneNodeUnrooted &node = **itr;
			GeneNodeUnrooted &parent = *node.parent();
			GeneNodeUnrooted &sibling = *node.getSibling();
			if (!sibling.belongs2GammaTree()) {
				// both nodes are in the ghost-node-tree
				node.secondarymapping->gain++;
			} else {
				// only one node is in the ghost-node-tree
				const int umapid = parent.secondarymapping->no;
				const int vmapid = sibling.secondarymapping->no;
				const int omapid = node.secondarymapping->no;
				if ((vmapid < umapid) && (umapid < omapid)) {
					parent.secondarymapping->lost[0]++;
				} else
				if ((omapid < umapid) && (umapid < vmapid)) {
					parent.secondarymapping->lost[1]++;
				}
			}
		}
	}

	// ------------------------------------------------------------------------------------------------------
	// calculate gene duplication score
	template<class GeneTree>
	inline unsigned int getScore(GeneTree &genetree) {
		unsigned int score = 0;
		getScoreDFS(score, genetree.root);
		return score;
	}
	template<class GeneNode>
	inline void getScoreDFS(unsigned int &score, GeneNode *&genenode) {
		if (!genenode->isLeaf()) {
			// child left
			GeneNode *&child1 = genenode->child(0);
			// child right
			GeneNode *&child2 = genenode->child(1);

			// if one child mapps to the same node as the parent increase the gene duplication score by 1
			SpeciesNode *&nodemap = genenode->getMapping();
			if ((nodemap == child1->getMapping()) || (nodemap == child2->getMapping())) score++;

			// DFS to the left
			getScoreDFS(score, child1);
			// DFS to the right
			getScoreDFS(score, child2);
		}

	}

	// ------------------------------------------------------------------------------------------------------
	// outputs the trees into a string stream
	friend ostream & operator << (ostream &os, TreeSet &t);
};

// outputs the trees into a string stream
ostream & operator << (ostream &os, TreeSet &t) {
	for(int i = 0; i < t.genetree_rooted.size(); i++ ) os << "gene tree" << i << ':' << endl << *t.genetree_rooted[i] << endl;
	for(int i = 0; i < t.genetree_unrooted.size(); i++ ) os << "gene tree" << i << ':' << endl << *t.genetree_unrooted[i] << endl;
	os << "species tree:" << endl << *t.speciestree << endl;
	return os;
}

#endif
