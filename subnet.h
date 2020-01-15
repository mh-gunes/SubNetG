/*
 * subnet.h
 *
 *  Created on: Nov 26, 2010
 *      Author: mba
 */

#ifndef SUBNET_H_
#define SUBNET_H_
#include "node.h"
#include<vector>
class node; //forward declaration, to avoid circular dependency;
class subnet {
public:
	int id;
	std::vector<node *> nodeList;
	long int desiredSize;
	long int currentSize;
	subnet();
	virtual ~subnet();
};

#endif /* SUBNET_H_ */
