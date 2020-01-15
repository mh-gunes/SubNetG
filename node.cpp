/*
 * node.cpp
 *
 *  Created on: Nov 26, 2010
 *      Author: mba
 */

#include "node.h"

node::node(int idNo) {
	id=idNo;
	desiredAlias=0;
	currentAlias=0;
	currentDegree=0;
}

node::~node() {
	// TODO Auto-generated destructor stub
}
