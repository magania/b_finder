/*
 * VrtSaver.cpp
 *
 *  Created on: Mar 17, 2009
 *      Author: magania
 */
#include <VrtSaver.h>

VrtSaver::VrtSaver(const char* name, TTree &tree){
	tree.Branch(strccat(name, "_x"), &x, strccat(name, "_x/D"));
	tree.Branch(strccat(name, "_y"), &y, strccat(name, "_y/D"));
	tree.Branch(strccat(name, "_z"), &z, strccat(name, "_z/D"));

	tree.Branch(strccat(name, "_v11"), &v11, strccat(name, "_v11/D"));
	tree.Branch(strccat(name, "_v12"), &v12, strccat(name, "_v12/D"));
	tree.Branch(strccat(name, "_v13"), &v13, strccat(name, "_v13/D"));
	tree.Branch(strccat(name, "_v22"), &v22, strccat(name, "_v22/D"));
	tree.Branch(strccat(name, "_v23"), &v23, strccat(name, "_v23/D"));
	tree.Branch(strccat(name, "_v33"), &v33, strccat(name, "_v33/D"));

	tree.Branch(strccat(name, "_chi2"), &chi2, strccat(name, "_chi2/D"));
}

VrtSaver::~VrtSaver(){
	return;
}

void VrtSaver::fill(const AA::Vrt &v) {
	x=v.x(1);
	y=v.x(2);
	z=v.x(3);

	v11=v.v(1,1);
	v12=v.v(1,2);
	v13=v.v(1,3);
	v22=v.v(2,2);
	v23=v.v(2,3);
	v33=v.v(3,3);

	chi2 = v.chi2();
}

const char* VrtSaver::strccat(const char *a, const char *b){
	char *temp = new char[20];
	names.push_back(temp);
	strcpy(temp, a);
	const char *result = strcat(temp,b);
	return result;
}
