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


#ifndef NODE_H
#define NODE_H

// ------------------------------------------------------------------------------------------------------------------
// this is a general node of a binary tree
// it contains links to the parent and childen nodes
template<class T>
class TreeNode {
public:
	// access parent
	virtual T* &parent() = 0;

	// access child
	virtual T* &child(const int i) = 0;

	// access child
	inline void swapChildren() {
		T *n = child(0);
		child(0) = child(1);
		child(1) = n;
	}

	// return true if this node is a leaf (has no children)
	inline bool isLeaf() {
		return ((child(0) == NULL) && (child(1) == NULL));
	}

	// return true if this node is the root node (has no parent)
	inline bool isRoot() {
		return (parent() == NULL);
	}

	// returns the child nodes
	inline void getChildren(T *child[]) {
		child[0] = this->child(0);
		child[1] = this->child(1);
	}

	// returns the silbing node
	// if the node has no sibling it return NULL;
	inline T *&getSibling() {
		#ifdef DEBUG
		if (parent() == NULL) EXCEPTION("getSibling: node has no sibling");
		#endif
		return parent()->child(0) == this ? parent()->child(1) : parent()->child(0);
	}

	// connect to a node as a child node
// 	inline void connectAsChild(T *&parent) {
// 		this->parent = parent;
// 		if (parent == NULL) return;
// 		for (int i=0; i<2; i++) {
// 			if (parent->child[i] == NULL) {
// 				parent->child[i] = this;
// 				return;
// 			}
// 		}
// 		EXCEPTION("connectAsChild: node " << this << " already has 2 child nodes connected");
// 	}
//
// 	// connect to a node as a certain child node
// 	inline void connectAsChild(T *parent, const int child) {
// 		this->parent() = parent;
// 		parent->child(child) = this;
// 	}

	template<class T2>
	friend ostream & operator << (ostream & os, TreeNode<T2> & m);
};

// output a node into a string stream
template<class T>
ostream & operator << (ostream & os, TreeNode<T> & m) {
	os << "addr:" << &m
	   << ", parent:" << m.parent()
	   << ", children:{" << m.child(0) << ", " << m.child(1) << '}';
	return os;
}

// ------------------------------------------------------------------------------------------------------------------
// node of a rooted binary tree
template<class T>
class TreeNodeRooted : public TreeNode<T> {
protected:
	T* parentlink;
	T* childlink[2];

public:
	// create a node (with a given parent)
	TreeNodeRooted(T *parent = NULL) : parentlink(parent) {
		for (int i = 0; i < 2; i++) childlink[i] = NULL;
	}

	// destroy a node
	virtual ~TreeNodeRooted() {}

	// access parent
	inline T* &parent() {
		return parentlink;
	}

	// access child
	inline T* &child(const int i) {
		return childlink[i];
	}

	template<class T2>
	friend ostream & operator << (ostream & os, TreeNodeRooted<T2> & m);
};

// output a node into a string stream
template<class T>
ostream & operator << (ostream & os, TreeNodeRooted<T> & m) {
	os << "addr:" << &m
	   << ", parent:" << m.parent()
	   << ", children:{" << m.child(0) << ", " << m.child(1) << '}';
	return os;
}

// ------------------------------------------------------------------------------------------------------------------
// node of a rooted binary tree
template<class T>
class TreeNodeUnrooted : public TreeNode<T> {
protected:
	T *relative[3];

public:
	unsigned short parentno;
	#define DEFAULT_ROOT_ID 0

	// create a node (with a given parent)
	TreeNodeUnrooted(T *parent = NULL) {
		for (int i = 0; i < 3; i++) relative[i] = NULL;
		parentno = DEFAULT_ROOT_ID;
		relative[parentno] = parent;
	}

	// destroy a node
	virtual ~TreeNodeUnrooted() {}

	// access parent
	inline T* &parent() {
		return relative[parentno];
	}

	// assign parent node
	inline void assignParent(T *&parent) {
		for (int i=0; i < (sizeof(relative)/sizeof(T*)); i++) {
			if (relative[i] == parent) {
				parentno = i;
				return;
			}
		}
		EXCEPTION("parent not found in assignParent" << endl);
	}

	// access child
	inline T* &child(const int i) {
		return i < parentno ? relative[i] : relative[1+i];
	}

	// access an edge
	inline T* &direction(const int i) {
		return relative[i];
	}

	// return the id of an edge (unrooted)
	inline int getDirectionID(const T *node) {
		for (int i=0; i < (sizeof(relative)/sizeof(T*)); i++) {
			if (relative[i] == node) return i;
		}
		EXCEPTION("node not found in getDirection" << endl);
	}

	// return the direction of an edge (unrooted)
	inline T *&getDirection(const T *node) {
		for (int i=0; i < (sizeof(relative)/sizeof(T*)); i++) {
			if (relative[i] == node) return relative[i];
		}
		EXCEPTION("node not found in getDirection" << endl);
	}

	// returns the child nodes (all nodes that are not parent)
	inline void getChildrenDirected(T *&parent, T *child[]) {
		T **r = &relative[0];
		for (int i=0; i < (sizeof(relative)/sizeof(T*)-1); i++) {
			if (relative[i] == parent) r++;
			child[i] = *r++;
		}
	}

	template<class T2>
	friend ostream & operator << (ostream & os, TreeNodeUnrooted<T2> & m);
};

// output a node into a string stream
template<class T>
ostream & operator << (ostream & os, TreeNodeUnrooted<T> & m) {
	os << "addr:" << &m
	   << ", parent:" << m.parent()
	   << ", children:{" << m.child(0) << ", " << m.child(1) << '}';
	return os;
}

// ------------------------------------------------------------------------------------------------------------------
// contains the ID (name) of a node
// nodes with the same name are automatically grouped together
class NodeID {
public:
	static map<string,int> namelist;
	map<string,int>::iterator name;

	// create a node (with a given parent)
	NodeID(const string &id) {
		name = namelist.find(id);
		if (name == namelist.end()) { // name doesn't exist
			namelist[id] = 0;
			name = namelist.find(id);
		}
		name->second++;
	}

	// destroy a node
	virtual ~NodeID() {
		name->second--;
		if (name->second == 0) namelist.erase(name);
	}

	// return the name
	inline const string &getName() {
		return name->first;
	}

	friend ostream & operator << (ostream & os, NodeID & m);
};
map<string,int> NodeID::namelist;

// output NodeID into a string stream
ostream & operator << (ostream & os, NodeID & m) {
	os << "name:" << m.name->first;
	return os;
}

#endif
