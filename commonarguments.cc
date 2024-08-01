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

bool helpFlag = false;
{
	const Argument *arg1 = Argument::find("-h");
	const Argument *arg2 = Argument::find("--help");
	const Argument *arg = arg1 != NULL ? arg1 : arg2;
	if (arg != NULL) {
		helpFlag = true;
		if (arg->hasValue()) EXCEPTION("--help doesn't need a value");
	}
}

quiet = false;
{
	const Argument *arg1 = Argument::find("-q");
	const Argument *arg2 = Argument::find("--quiet");
	const Argument *arg = arg1 != NULL ? arg1 : arg2;
	if (arg != NULL) {
		quiet = true;
		if (arg->hasValue()) EXCEPTION("--quiet doesn't need a value");
	}
}

bool versionFlag = false;
{
	const Argument *arg1 = Argument::find("-v");
	const Argument *arg2 = Argument::find("--version");
	const Argument *arg = arg1 != NULL ? arg1 : arg2;
	if (arg != NULL) {
		versionFlag = true;
		if (arg->hasValue()) EXCEPTION("--version doesn't need a value");
	}
}

// bool helpFlag = (Argument::find("-h") != NULL) || (Argument::find("--help") != NULL);
// quiet = (Argument::find("-q") != NULL) || (Argument::find("--quiet") != NULL);
// bool versionFlag = (Argument::find("-v") != NULL) || (Argument::find("--version") != NULL);

// choose input
istream *in;
ifstream infile;
{
	const Argument *arg1 = Argument::find("-i");
	const Argument *arg2 = Argument::find("--input");
	const Argument *arg = arg1 != NULL ? arg1 : arg2;
	if (arg == NULL) in = &cin;
	else {
		string filename;
		arg->convert(filename);
		msgout << "Input file: " << filename << endl;
		infile.open(filename.c_str());
		if (!infile) EXCEPTION("input file " << filename << " not found");
		in = &infile;
	}
}
Input input(in);

// choose output
ostream *out;
ofstream outfile;
{
	const Argument *arg1 = Argument::find("-o");
	const Argument *arg2 = Argument::find("--output");
	const Argument *arg = arg1 != NULL ? arg1 : arg2;
	if (arg == NULL) out = &cout;
	else {
		string filename;
		arg->convert(filename);
		msgout << "Output file: " << filename << endl;
		outfile.open(filename.c_str());
		if (!outfile) EXCEPTION("Output file " << filename << " not created");
		out = &outfile;
	}
}

// output version number
if (versionFlag) {
	cout << version << endl;
	return 0;
}

// complain about unknown arguments
vector<Argument*> unused = Argument::unusedArgs();
for (vector<Argument*>::iterator itr = unused.begin(); itr != unused.end(); itr++) {
	EXCEPTION("Unknown argument " << (*itr)->key << endl << "Try " << argsv[0] << " --help for a list of arguments");
}
