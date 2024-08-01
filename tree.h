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

#ifndef TREE_H
#define TREE_H

template<class NODE, class LEAFNODE>
class Tree {
public:
	typedef NODE Node;
	typedef LEAFNODE LeafNode;

	NODE *root;
	vector<NODE*> nodes;
	vector<LEAFNODE*> leafnodes;

	Tree() {
		root = NULL;
	}

	Tree(const Tree &r) {
		// copying gives you an empty tree rigth now!!!
		if (r.root != NULL) cerr << "WARNING: Tree is not going to be copied!" << endl;
	}

	virtual ~Tree() {
		#ifdef DEBUG
		checkStructure();
		#endif
		for (typename vector<NODE*>::iterator itr = nodes.begin(); itr != nodes.end(); itr++) delete *itr;
		nodes.clear();
	}

	// performs a NNI operation on the edge above 'node' and it's parent's sibling.
	inline void oprNNI(NODE *node) {
		NODE *a = node;
		NODE *ap = a->parent();
		int ai = ap->child(0) == a ? 0 : 1;
		NODE *b = a->parent()->getSibling();
		NODE *bp = b->parent();
		int bi = bp->child(0) == b ? 0 : 1;
		a->parent() = bp;
		bp->child(bi) = a;
		b->parent() = ap;
		ap->child(ai) = b;
	}

	// moves a subtree to a new location in the tree
	// subroot = root node of the tree to be pruned
	// targetnode and subroot are going to be sibling in the new tree
	inline void moveSubtree(NODE *subroot, NODE *targetnode) {
		// prune out subtree
		#ifdef DEBUG
		if (subroot == NULL) EXCEPTION("moveSubtree: subtree root node is NULL");
		if (subroot->parent() == NULL) EXCEPTION("moveSubtree: no parent node");
		for (NODE *n = targetnode; n != NULL; n = n->parent()) {
			if (n == subroot) EXCEPTION("moveSubtree: cannot move into subtree");
		}
		#endif
		NODE *p = subroot->parent();
		if (targetnode == p) return;
		NODE *pp = p->parent();
		const int pi = p->child(0) != subroot ? 0 : 1;
		NODE *pc = p->child(pi);
		pc->parent() = pp;
		if (pp != NULL) {
			const int ppi = pp->child(0) == p ? 0 : 1;
			pp->child(ppi) = pc;
		} else {
			root = p->child(pi);
			p->child(pi) = NULL;
		}

		// insert subtree
		#ifdef DEBUG
		if (targetnode == NULL) EXCEPTION("moveSubtree: target node is NULL");
		if (targetnode == subroot) EXCEPTION("moveSubtree: invalid target node");
		#endif
		NODE *q = targetnode->parent();
		if (q != NULL) {
			const int qi = q->child(0) == targetnode ? 0 : 1;
			q->child(qi) = p;
		} else {
			root = p;
			p->parent() = NULL;
		}
		p->parent() = q;
		targetnode->parent() = p;
		p->child(pi) = targetnode;
	}

	// delete a leaf node
	void deleteLeafNode(LEAFNODE *node) {
		// prune out subtree
		#ifdef DEBUG
		if (node == NULL) EXCEPTION("removeLeafNode: subtree root node is NULL");
		if (node->parent() == NULL) EXCEPTION("removeLeafNode: no parent node");
		#endif
		NODE *p = node->parent();
		NODE *pp = p->parent();
		const int pi = p->child(0) != node ? 0 : 1;
		NODE *pc = p->child(pi);
		pc->parent() = pp;
		if (pp != NULL) {
			const int ppi = pp->child(0) == p ? 0 : 1;
			pp->child(ppi) = pc;
		} else {
			root = p->child(pi);
			p->child(pi) = NULL;
		}

		// delete leaf node and it's parent node
		// from leaf set
		for (typename vector<LEAFNODE*>::iterator itr = leafnodes.begin(); itr != leafnodes.end(); itr++) {
			if (*itr == node) {
				leafnodes.erase(itr);
				break;
			}
		}
		// from node set
		for (typename vector<NODE*>::iterator itr = nodes.begin(); itr != nodes.end(); itr++) {
			if (*itr == node) {
				nodes.erase(itr);
				break;
			}
		}
		for (typename vector<NODE*>::iterator itr = nodes.begin(); itr != nodes.end(); itr++) {
			if (*itr == p) {
				nodes.erase(itr);
				break;
			}
		}
		// free the memory
		delete p;
		delete node;
		#ifdef DEBUG
		checkStructure();
		#endif
	}

	// take out root node
	inline void disconnectRoot() {
		NODE *u = root->child(0);
		NODE *v = root->child(1);
		u->parent() = v;
		v->parent() = u;
	}

	// insert root node
	inline void connectRoot(NODE *u, NODE *v) {
		root->child(0) = u;
		root->child(1) = v;
		u->getDirection(v) = root;
		v->getDirection(u) = root;
	}

	#ifdef DEBUG
	void checkStructure() {
		// DFS starting from root node
		int size = 0;
		checkStructureDFS(root, NULL, size);
		// check amount of registered nodes
		if (size != nodes.size())
			EXCEPTION("tree structure broken - " << size <<" nodes in the tree, but " << nodes.size() << " registered");
		// are all nodes of known datatype
		for (typename vector<NODE*>::iterator itr = nodes.begin(); itr != nodes.end(); itr++)
			EXCEPTIONDATATYPE2("checkStructure", **itr, NODE, LEAFNODE);
		// are all leaf nodes of datatype LEAFNODE
		for (typename vector<LEAFNODE*>::iterator itr = leafnodes.begin(); itr != leafnodes.end(); itr++)
			EXCEPTIONDATATYPE("checkStructure", **itr, LEAFNODE);
	}
	void checkStructureDFS(NODE *&node, NODE *parent, int &size) {
		if (node == NULL) return;
		size++;
		if (node->parent() != parent) EXCEPTION("tree structure broken - node " << node << " has wrong parent");
		NODE *children[2];
		node->getChildren(children);
		for (int i = 0; i < 2; i++) {
			checkStructureDFS(children[i], node, size);
		}
	}
	#endif

	// outputs a subtree into a string stream
	static void outputSubtree(ostream &os, NODE *&node, int lvl) {
		// indent according to the depth of the tree
		for (int i=0; i<lvl; i++) os << "  ";
		// output the node
		EXCEPTIONDATATYPE2("outputSubtree", *node, LEAFNODE, NODE);
		if (typeid(*node)==typeid(LEAFNODE)) os << *((LEAFNODE *)node);
		else if (typeid(*node)==typeid(NODE)) os << *((NODE *)node);
		else EXCEPTION("treenode of unknown type");

		if ((node->child(0)!=NULL) || (node->child(1)!=NULL)) os << " (";
		os << endl;
		// proceed with the children nodes
		if (node->child(0)!=NULL) outputSubtree(os, node->child(0), lvl+1);
		if (node->child(1)!=NULL) outputSubtree(os, node->child(1), lvl+1);
		// indent according to the depth of the tree
		if ((node->child(0)!=NULL) || (node->child(1)!=NULL)) {
			for (int i=0; i<lvl; i++) os << "  ";
			os << ')' << endl;
		}
	}

	// outputs the tree into a string stream
	template<class N1, class N2> friend ostream & operator << (ostream &os, Tree<N1,N2> &t);
};

// outputs the tree into a string stream
template<class NODE, class LEAFNODE>
ostream & operator << (ostream &os, Tree<NODE, LEAFNODE> &t) {
	t.outputSubtree(os, t.root, 0);
	return os;
}

// ------------------------------------------------------------------------------------------------------------------
// extends a binary tree by input/output functions
template<class TREE>
class TreeIO : public TREE {
public:
	using TREE::root;
	using TREE::nodes;
	using TREE::leafnodes;

	typedef typename TREE::Node NODE;
	typedef typename TREE::LeafNode NAMEDNODE;

	// read the tree from a string stream (newick formatted e.g. ((name1,name2),name3);)
	public: bool stream2tree(Input &is) {
		char c;
		input = &is;
		if (!is.nextAnyChar(c)) return false;
		if (c != '(') EXCEPTION("problem in the tree expression " << input->getPos() << " " << c << "this");
		stream2tree(c, is);
		return true;
	}

	// read the tree from a string stream; given the first character of the stream
	public: void stream2tree(char &c, Input &is) {
		int lvl = 1;
		input = &is;

		// start reading the tree ... begin with root node
		commentcallback.comments.clear();
		input->setCommentCallback(&commentcallback);
		stream2nextNode(c, -1);
		input->setCommentCallback(NULL);

		// ';' terminator found?
		if (c != ';') EXCEPTION("missing ';' in tree expression " << input->getPos());
		buildTree();
		nodeorder.clear();
	}

	// called when a comment was explored
	class CommentCallBack : public InputCommentCallbackBase {
		void commentFound(string &str) {
			comments.push_back(str);
		}
		public: vector<string> comments;
	} commentcallback;

	// outputs the tree in newick format
	public: void tree2newick(ostream &os) {
		subtree2newickDFS(os, root);
		os << ';';
	}
	protected: void subtree2newickDFS(ostream &os, NODE *&node) {
		if (node == NULL) return;

		// if leaf node then output the node
		if (!node->isLeaf()) {
			EXCEPTIONDATATYPE2("tree2newickDFS", *node, NAMEDNODE, NODE);
			if ((node->child(0)!=NULL) || (node->child(1)!=NULL)) os << '(';
			// proceed with the children nodes
			subtree2newickDFS(os, node->child(0));
			if (node->child(1)!=NULL) os << ',';
			subtree2newickDFS(os, node->child(1));
			// close the subtree
			if ((node->child(0)!=NULL) || (node->child(1)!=NULL)) os << ')';
		}

		// output node name and properties
		if (typeid(*node)==typeid(NAMEDNODE)) {
			namednode2newick(os, *((NAMEDNODE *)node));
		} else {
			node2newick(os, *node);
		}
	}

	// output a named tree node in newick format
	public: virtual void namednode2newick(ostream &os, NAMEDNODE &node) {
		os << name2newick(node.getName());
	}

	// output a unnamed tree node in newick format
	public: virtual void node2newick(ostream &os, NODE &node) {
	}

	protected: Input *input;

	// temporary all-purpose node
	protected: class TempNode {
	public:
		TempNode() {
			parent = -1;
			name = "";
			weight = 0;
			node = NULL;
		}
		TempNode(int parent) : parent(parent) {
			name = "";
			weight = 0;
			node = NULL;
		}
		TempNode(int parent, string name) : parent(parent), name(name) {
			weight = 0;
			node = NULL;
		}
		TempNode(const TempNode& p) {
			parent = p.parent;
			name = p.name;
			node = p.node;
			weight = p.weight;
			comments = p.comments;
		}

		int parent;
		string name;
		double weight;
		vector<string> comments;
		NODE* node;
	};
	vector<TempNode> nodeorder;

	// build the tree from all-purpose nodes
	protected: void buildTree() {
		// resolve multifurcations by making them arbitrarily binary
		vector<vector<int> > children(nodeorder.size());
		for (int i=0; i<nodeorder.size(); i++) {
			TempNode *n = &nodeorder[i];
			if (n->parent != -1) children[n->parent].push_back(i);
		}
		for (int i=0; i<children.size(); i++) {
			while (children[i].size() > 2) {
				// create a common parent
				int newid = nodeorder.size();
				nodeorder.push_back(TempNode(i));

				// pick 1st node randomly
				int l = (int) (children[i].size() * (rand() / (RAND_MAX + 1.0)));
				nodeorder[children[i][l]].parent = newid;
				children[i].erase(children[i].begin()+l);

				// pick 2nd node randomly
				int r = (int) (children[i].size() * (rand() / (RAND_MAX + 1.0)));
				nodeorder[children[i][r]].parent = newid;
				children[i].erase(children[i].begin()+r);

				// add new node as child
				children[i].push_back(newid);
			}
		}

		// create all nodes
		for (int i=0, last=nodeorder.size(); i<last; i++) {
			TempNode &n = nodeorder[i];
			createNode(n);
			nodes.push_back(n.node);
		}

		// connect all nodes
		for (int i=0; i<nodeorder.size(); i++) {
			TempNode &n = nodeorder.at(i);
			if (n.parent == -1) {
				root = n.node;
				n.node->parent() = NULL;
			} else {
				TempNode &p = nodeorder.at(n.parent);
				n.node->parent() = p.node;
				if (p.node->child(0) == NULL) {
					p.node->child(0) = n.node;
				} else
				if (p.node->child(1) == NULL) {
					p.node->child(1) = n.node;
				} else {
					cout<<"  node:"<<i<<" parent:"<<n.parent<<endl;
					EXCEPTION("buildTree: cannot connect more than 2 children to one parent");
				}
			}
		}

		// build leaf set
		for (int i=0; i<nodeorder.size(); i++) {
			TempNode &n = nodeorder.at(i);
			if (n.node->isLeaf()) leafnodes.push_back((NAMEDNODE*)n.node);
		}
	}

	// creates a tree node from a temporary node
	public: virtual void createNode(TempNode &n) {
		if (n.name.empty()) {
			n.node = new NODE(NULL);
		} else {
			n.node = new NAMEDNODE(n.name, NULL);
		}
	}

	// process an internal or leaf node
	protected: void stream2nextNode(char &c, const int &parentid) {
		// create a temporary node
		const int nodeid = nodeorder.size();
		nodeorder.push_back(TempNode(parentid));
		#define tempnode nodeorder[nodeid]

		// add comments prior to the node expression
		tempnode.comments = commentcallback.comments;
		commentcallback.comments.clear();

		if (c == '(') { // internal node
			int childCount = 0;
			do {
				// check for multifurcations
				childCount++;
				if (childCount > 2) WARNING("multifurcation detected at " << input->getLastPos())

				// node/subtree
				if (!input->nextAnyChar(c)) EXCEPTION("unexpected end in the tree expression " << input->getPos());
				stream2nextNode(c, nodeid);
			} while (c == ',');
			if (c != ')') EXCEPTION("missing ')' in input tree expression " << input->getLastPos());

			// read name
			if (!input->nextAnyChar(c)) EXCEPTION("unexpected end in the tree expression " << input->getPos());
			if ((c != ')') && (c != ',') && (c != ':')) {
				input->pushBack(c);
				tempnode.name = input->getName();

				if (!input->nextAnyChar(c)) EXCEPTION("unexpected end in the tree expression " << input->getPos());
			}
		} else { // leaf node
			input->pushBack(c);
			tempnode.name = input->getName();

			if (!input->nextAnyChar(c)) EXCEPTION("unexpected end in the tree expression " << input->getPos());
		}
		// read weight
		if (c == ':') {
			double weight = input->readNumber(); // read weight + get next character
			tempnode.weight = weight;

			if (!input->nextAnyChar(c)) EXCEPTION("unexpected end in the tree expression " << input->getPos());
		}
		// add comments post to the node expression
		tempnode.comments.insert(tempnode.comments.end(), commentcallback.comments.begin(), commentcallback.comments.end());
		commentcallback.comments.clear();
		#undef tempnode
	}

	void callbackComment(string &str) {
		cout << endl << str << endl;
	}
};

#endif
