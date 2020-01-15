/*
 * network.h
 *
 *  Created on: Jan 2, 2012
 *      Author: mba
 */

#ifndef NETWORK_H_
#define NETWORK_H_

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <time.h>
#include  <vector>
#include "math.h"
#include "node.h" //standard network node
#include "subnet.h"//subnet
using namespace std;

#define maxDistributionIndex 100000

class network
{
public:
	network(float nadSLopeUser, float ssdSlopeUser, int NN_input);
	virtual ~network();
	float nadSlope;
	float ssdSlope;
	long long int NN; //desired total number of nodes//global
	long long int IpCount;//total number of interfaces//global
	long long int subnetCount;//total number of subnets//global

	vector<int> NAD;
	vector<int> SSD;
	vector<long long int> NDD;
	vector<int> currentNDD;

	vector<subnet> subnetRoot; //vector of all subnets in the network.
	vector<node> nodeRoot;  //vector of all nodes in the network.
	void generate(void);
	void calculateNAD(void);
	void calculateNADHT(void);
	void calculateSSD(void);
	void calculateSSDHT(void);
	bool isConnectable(void);
	void generateNodes(void);
	void generateSubnets(void);
	void reverseSubnetOrder(void);
	void reverseNodeOrder(void);
	void interfacePreferentialAttachment(void);
	void degreePreferentialAttachment(void);
	void mixSubnets(void);
	void mixNodes(void);
	void connectMST(void);
	void connectRestofMST(void);
	void inOrderConnect(void);
	void randomConnect(void);
	void randomAttach(void);
	bool isAttachable(subnet &s, node &n);
	void attach(subnet & s, node & n);
	void randomizeNetwork(void);
	void calculateDegree(void);
	bool isSwappable(node &n1,int s1, node &n2, int s2);
	void calculateCurrentNDD(void);
	void printCurrentNDD(void);
	void printNDD(void);
	void powerifyNDD(void);
	bool isConverge(node &n1, int s1, node &n2, int s2);
	void printMap(void);
	void printSSD(void);
	void printNAD(void);
	void printEdgeList(void);
	void setClusterIds(void);
	void printClusterIds(void);
};

#endif /* NETWORK_H_ */
