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

#ifndef ARGUMENT_H

class Argument {
public:
	static map<string, Argument> args;
	string key;
	string value;
	bool visited;

	Argument() {
		key = "";
		value = "";
		visited = false;
	}
	Argument(string key) {
		this->key = key;
		value = "";
		visited = false;
	}
	Argument(string key, string value) {
		this->key = key;
		this->value = value;
		visited = false;
	}

	inline bool operator < (const Argument& ref) const {
		return key < ref.key;
	}
	inline bool operator > (const Argument& ref) const {
		return key > ref.key;
	}
	inline bool operator == (const Argument& ref) const {
		return key == ref.key;
	}

	// adds arguments from the command line
	static void add(int argc, char *argsv[]) {
		#define INSERT(key,value) { \
			if (args.find(key) != args.end()) EXCEPTION("more than one argument " << key); \
			args[key] = Argument(key, value); \
			key = ""; \
		}
		if (argc <= 1) return;
		string key = "";
		for (int i=1; i<argc; i++) {
			string str(argsv[i]);
			if (str[0] == '-') {
				if (!key.empty()) INSERT(key,"")
				const int pos = str.find("=");
				if (pos != string::npos) {
					key = str.substr(0,pos);
					string value = str.substr(pos+1);
					INSERT(key,value);
				} else {
					key = str;
				}
			} else {
				string value = str;
				INSERT(key,value);
			}
		}
		if (!key.empty()) INSERT(key,"")
		#undef INSERT
	}

	// searches for a particular argument
	inline static Argument *find(const string &key) {
		map<string, Argument>::iterator itr = args.find(key);
		if (itr == args.end()) return NULL;
		Argument &arg = itr->second;
		arg.visited = true;
		return &arg;
	}

	// searches for any argument
	inline static Argument *findAny(vector<string> &list) {
		for (vector<string>::iterator itr = list.begin(); itr != list.end(); itr++) {
			Argument *arg = find(*itr);
			if (arg != NULL) return arg;
		}
		return NULL;
	}

	// convert a character string into another datatype
	template<class T> void convert(T &var) const;

	// true if argument has a value
	bool hasValue() const {
		return (!value.empty());
	}

	// return unused arguments
	static vector<Argument*> unusedArgs() {
		vector<Argument*> unused;
		for (map<string, Argument>::iterator itr = args.begin(); itr != args.end(); itr++) {
			if (!itr->second.visited) unused.push_back(&itr->second);
		}
		return unused;
	}
};
map<string, Argument> Argument::args;

// convert a character string into another datatype
template<class T>
void Argument::convert(T &var) const {
	istringstream ist(value);
	if (!(ist >> var)) EXCEPTION("invalid value of argument " << key << ' ' << value);
}
// specialization for string
template<>
void Argument::convert(string &var) const {
	var = value;
}

#endif
