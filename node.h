/*
 * node.h
 *
 *  Created on: Nov 26, 2010
 *      Author: mba
 */

#ifndef NODE_H_
#define NODE_H_
#include "subnet.h"
#include <vector>
class subnet;
class node {
public:
	int id;					//unique id
	std::vector<subnet * > subnetList;  //a pointer to the Linked list which contains the subnets this node is attached
	long int desiredAlias;
	long int currentAlias;
	long int currentDegree;
	node(int idNo=-1);
	virtual ~node();
};

#endif /* NODE_H_ */
