/*
 * network.cpp
 *
 *  Created on: Jan 2, 2012
 *      Author: mba
 */

#include "network.h"
using namespace std;

ofstream logfile;
ofstream alpha;
network::network(float nadSlopeUser, float ssdSlopeUser, int NN_input)
{
// TODO Auto-generated constructor stub
	NN=NN_input;
	subnetCount=0;
	nadSlope=nadSlopeUser;
	ssdSlope=ssdSlopeUser;
	NAD.resize(maxDistributionIndex,0);
	SSD.resize(maxDistributionIndex,0);
	logfile.open("log.txt");
	alpha.open("alpha.txt",fstream::app);
}

network::~network()
{
	logfile<<"exiting...\n";
	logfile.close();
	alpha.close();
}

void network::generate(void)
{
		calculateNADHT();
		calculateSSDHT();
//		if( isConnectable())
//			cout<< "connectable"<<endl;
//		else
//			cout<< "NOT connectable"<<endl;

		generateNodes();
		generateSubnets();

//		mixSubnets();
//		mixNodes();

		reverseSubnetOrder();
//		reverseNodeOrder();
		inOrderConnect();
//		interfacePreferentialAttachment();

		calculateDegree();
		calculateCurrentNDD();
		printCurrentNDD();
		printNAD();
		printSSD();
		printEdgeList();
		printMap();
		//
//		mixSubnets();
//
//		mixNodes();
//
//		connectMST();
//
//		connectRestofMST();
//
//		randomConnect();
//
//		randomAttach();
//
//		printMap();
////
//		calculateDegree();
//
//		calculateCurrentNDD();
//
//		printCurrentNDD();
////
//		printSSD();
//
//		printNAD();
////
//		printEdgeList();
//
//		cout<<"done!"<<endl;
//		cout<<"ssd:"<<ssdSlope<<"     nad:"<<nadSlope<<"      ndd:"<<nddSlope<<endl;
}

void network::calculateNAD(void)
{
	//Need to scale the NAD curve such that its integral(sum) exactly matches the NN.
	//to do this we first scale it with normal ratio. Due to integerization we come up with a smaller integral than desired NN.
	//Than we slowly scale the curve(floating point) with 1.000001 until integral exactly matches NN.
//	cout<<"calculating NAD"<<endl;
	//float nadRef[maxDistributionIndex];
	vector<float> nadRef(maxDistributionIndex,0);

	long long int nadRefCount=0;
	for (unsigned int i=1; i<nadRef.size(); i++)
	{
		nadRef.at(i)=10000000*pow(i,nadSlope);
		nadRefCount+=nadRef.at(i);
	}
	cout<<nadRefCount<<endl;
	long long int integral=0;
	for(unsigned int i=1; i<nadRef.size(); i++)
	{
		nadRef.at(i)=(  nadRef.at(i) * ( float(NN)/float(nadRefCount) )    );
		integral+= int(nadRef.at(i));
	}
//	cout<<"first attempt. NAD integral:"<<integral<<endl;
	while(integral<NN)
	{
		integral=0;
		for(unsigned int i=1; i<nadRef.size(); i++)
		{
			nadRef.at(i)*=1.000001;
			integral+= int (nadRef.at(i));
		}
	}
	//cout<<"integral of NAD curve="<<integral<<endl;
	//in some cases scaling further overflows the integral. Than we make little reduction
	if (integral>NN and (integral-NN)<maxDistributionIndex)
	{
		for(int i=1; i<=(integral-NN); i++)
		{
			nadRef[i]--;
		}
	}
	integral=0;
	IpCount=0;
	for (unsigned int i=0; i<NAD.size(); i++)
	{
		NAD.at(i)=int (nadRef.at(i));
		integral+=NAD.at(i);
		IpCount+=NAD.at(i)*i;
	}
//	cout<<"final NAD integral:"<<integral<<endl;
	if(integral!=NN)
		cerr<<"NAD integral cannot match NN";
	cout<<"Total IPCount=    "<<IpCount<<endl;

//	logfile<<"NAD calculated"<<endl;
}

void network::calculateNADHT(void){
	long long int Count=0;
	vector<int> referenceCurve(maxDistributionIndex,0);
	for (unsigned int i=0; i<referenceCurve.size(); i++)
		referenceCurve.at(i)=0;
	for (unsigned int i=1; i<referenceCurve.size(); i++)
	{
		referenceCurve.at(i)=1000000000*pow(i,nadSlope);
		Count+=referenceCurve.at(i);
	}

	for(int i=0; i<NN; i++)
	{
		long long int random=int (float(rand()) *  (float(Count)/float(RAND_MAX)));
		int j=1;
		while (random > 0)
		{
			random -= referenceCurve.at(j);
			j++;
		}
		NAD.at(j-1)++;
	}
	IpCount=0;
	for (unsigned int i=0; i<NAD.size(); i++)
	{
		IpCount+=NAD.at(i)*i;
	}
	cout<<IpCount<<endl;
//	int acc=0;
//	for (int i=0;i<maxDistributionIndex; i++)
//		if (NAD[i]>0)
//			acc+=NAD[i];
//			cout<<i<<"\t"<<NAD[i]<<endl;
//	cout<<acc;
}

void network::calculateSSD(void)
{
//	cout<<"calculating SSD"<<endl;
	//integral of (SSD[i]*i) curve must give the total number of interfaces
	//which is also the integral of (NAD[i]*i) curve.
	vector<float> ssdRef(maxDistributionIndex,0);
	//float ssdRef[maxDistributionIndex];
	//minimum subnet size is 2.
	long long int ssdRefCount=0;
	for (unsigned int i=2; i<ssdRef.size(); i++)
	{
		ssdRef.at(i)=10000000*pow(i,ssdSlope);
		ssdRefCount+=ssdRef.at(i)*i;
	}
//	cout<<ssdRefCount<<endl;
	long long int integral=0;
	for(unsigned int i=2; i<ssdRef.size(); i++)
	{
		ssdRef.at(i)=(  ssdRef.at(i) * ( float(IpCount)/float(ssdRefCount) )    );
		integral+= int(ssdRef.at(i))*i;
	}
//	cout<<"first attempt, SSD integral:"<<integral<<endl;
	while(integral<IpCount)
	{
		integral=0;
		for(unsigned int i=1; i<ssdRef.size(); i++)
		{
			ssdRef.at(i)*=1.000001;
			integral+= int (ssdRef.at(i))*i;
		}
	}
//	cout<<"integral of SSD curve="<<integral<<endl;
	//in some cases scaling further overflows the integral. Than we make little reduction
	if (integral>IpCount and (integral-IpCount)<maxDistributionIndex)
	{
		for(int i=2; i<=(integral-IpCount)+1; i++)
		{
			ssdRef.at(i)++;//both increment
			ssdRef.at(i+1)--;//and decrement should be done to reduce total integral by 1.
		}
	}
	integral=0;
	subnetCount=0;
	for (unsigned int i=0; i<SSD.size(); i++)
	{
		SSD.at(i)=int (ssdRef.at(i));
		integral+=SSD.at(i)*i;
		subnetCount+=SSD.at(i);
	}
//	cout<<"final SSD integral:"<<integral<<endl;
	if(integral!=IpCount)
		cerr<<"SSD[i]*i integral cannot match IPCount"<<endl;
	cout<<"Total subnet Count="<<subnetCount<<endl;

//	logfile<<"SSD calculated"<<endl;
}

void network::calculateSSDHT(void){

	long long int Count=0;
	vector<int> referenceCurve(maxDistributionIndex,0);
	for (unsigned int i=0; i<referenceCurve.size(); i++)
		referenceCurve.at(i)=0;
	for (unsigned int i=2; i<referenceCurve.size(); i++)
	{
		referenceCurve.at(i)=1000000000*pow(i,ssdSlope);
		Count+=referenceCurve.at(i);
	}
	long long int currentIPCount=0;
	subnetCount=0;
	while( currentIPCount <= IpCount+1  )
	{
		//normalize rand to our range,  % remainder would not be accurate enough
		long long int random=int (float(rand()) *  (float(Count)/float(RAND_MAX)));
		int j=1;
		while (random > 0)
		{
			random -= referenceCurve.at(j);
			j++;
		}
		if ((currentIPCount + j -1) > IpCount )
			break;
		else
			SSD.at(j-1)++;
			currentIPCount += j-1;
			subnetCount++;
	}
		if (SSD.at(2)<=0)
			cerr<<"Check network size"<<endl;
		SSD.at(2)--;
		SSD.at( 2+ IpCount - currentIPCount )++;


	int acc=0;
	for (unsigned int i=2;i<SSD.size(); i++)
		if (SSD.at(i)>0)
			acc+=SSD.at(i)*i;
//			cout<<i<<"\t"<<SSD[i]<<endl;
	cout<<acc<<endl;

}

bool network::isConnectable(void)
{//True if it is possible to generate a connected
 //graph with given alpha's and calculated SSD and NAD; False otherwise
	//if sum(i, NAD[i]*(i-1)) + 1 > subnetCount ; for i>1;  return True
	int accumulator=1;
	for (int i=2; i<maxDistributionIndex;i++)
	{
		accumulator+= NAD[i] * (i-1);
	}
	if (accumulator > subnetCount)
	{
		logfile<<"graph connectable!"<<endl;
		return true;
	}
	else
	{
		logfile<<"NOT connectable!"<<endl;
		return false;
	}
}

void network::generateNodes(void)
{//generate an array of nodes. total NN nodes
	nodeRoot.resize(NN,0);
	for(unsigned int i=0; i< nodeRoot.size(); i++)//assign their id
	{
		nodeRoot.at(i).id=i;
	}

	//make a temporary copy of NAD[i]
	vector <int> nadTemp(NAD.size(),0);
	for(unsigned int i=0; i<nadTemp.size(); i++)
	{
		nadTemp.at(i)=NAD.at(i);
	}
//	for(int i=0; i< maxDistributionIndex; i++)
//	{
//		cout<<i<<"   "<<NAD[i]<<endl;
//	}
	int i=0;//assign desired aliasCount for each node.
	for(unsigned int j=0; j< nadTemp.size(); j++)
	{
		while(nadTemp.at(j)>0)
		{
			nodeRoot.at(i).desiredAlias=j;
			i++;
			nadTemp.at(j)--;
		}
	}
//	for(int i=0;i<NN;i++)
//	{
//		cout<<i<<"   "<<nodeRoot[i].desiredAlias<<endl;
//	}
	//generate an array of subnet pointers.
	//Size of the array equals desired alias count for that node
	for (int i=0; i<NN;i++)
	{
		nodeRoot.at(i).subnetList.resize(nodeRoot.at(i).desiredAlias,NULL);
	}
	logfile<<"nodes generated"<<endl;
}

void network::generateSubnets(void)
{//generate an array of subnets. total subnetCount ones.
	subnetRoot.resize(subnetCount);
	//subnetRoot= new subnet[subnetCount];
	for(unsigned int i=0; i< subnetRoot.size(); i++)
	{
		subnetRoot.at(i).id=i;
	}

	//make a temporary copy of SSD[i]
	vector<int> ssdTemp(SSD.size(),0);
	for(unsigned int i=0; i<ssdTemp.size(); i++)
	{
		ssdTemp.at(i)=SSD.at(i);
	}
//	for(int i=0; i< maxDistributionIndex; i++)
//	{
//		cout<<i<<"   "<<SSD[i]<<endl;
//	}
	int i=0;//assign desiredSize for each subnet.
	for(unsigned int j=0; j< ssdTemp.size(); j++)
	{
		while(ssdTemp.at(j)>0)
		{
			subnetRoot.at(i).desiredSize=j;
			i++;
			ssdTemp.at(j)--;
		}
	}
//	for(int i=0;i<subnetCount;i++)
//	{
//		cout<<i<<"   "<<subnetRoot[i].desiredSize<<endl;
//	}

	//generate an array of node pointers.
	//Size of the array equals desireSize for that subnet
	for (int i=0; i<subnetCount;i++)
	{
		subnetRoot[i].nodeList.resize(subnetRoot.at(i).desiredSize,NULL);
	}

	logfile<<"subnets generated"<<endl;
}

void network::reverseSubnetOrder(void)
{ //assumes already sorted in increasing order wrt size
  //using shallow copy only, nodeList members should not have been allocated !!
	subnet temp;
	for(int i=0; i<subnetCount-i-1; i++)
	{
		temp=subnetRoot.at(i);
		subnetRoot.at(i)=subnetRoot.at(subnetCount-i-1);
		subnetRoot.at(subnetCount-i-1)=temp;
	}
	logfile<<"subnets reversed"<<endl;
}

void network::degreePreferentialAttachment(void)
{
//	vector<int> weightedNodeList;
//	long long int totalDegree=0;
//	for(unsigned int i=0; i<subnetCount; i++)
//		totalDegree+=subnetRoot.at(i).desiredSize * (subnetRoot.at(i).desiredSize-1);
//	weightedNodeList.reserve(totalDegree+NN);
//	for(unsigned int=0;i<NN;i++)
//		weightedNodeList.push_back(i);
//	for (unsigned int i=0 ; i< subnetRoot.size(); i++)
//	{
//		for (int j=0; j<subnetRoot.at(i).desiredSize; j++ )
//		{
//			int randomNode=weightedNodeList.at(rand() % weightedNodeList.size());
//			if ( nodeRoot.at(randomNode).desiredAlias > nodeRoot.at(randomNode).currentAlias)
//
//		}
//	}
}

void network::interfacePreferentialAttachment(void)
{
	//This function assumes node with id=i is the i-th entry in nodeRoot vector.
	//This assumption makes it much faster
	//Do not call this function after functions that modify nodeRoot (like reverseNodeOrder)

	//a list of node ids where a node appears number of Interface times
	vector<int> weightedNodeList;
	weightedNodeList.reserve(IpCount);
	for (unsigned int i=0;i< nodeRoot.size(); i++)
		for(int j=0; j<nodeRoot.at(i).desiredAlias;j++)
			weightedNodeList.push_back(nodeRoot.at(i).id);

	//randomize the order
	for(unsigned int i=0 ; i<weightedNodeList.size() ; i++)
	{
		int j= rand() % weightedNodeList.size();
		int temp = weightedNodeList.at(i);
		weightedNodeList.at(i) = weightedNodeList.at(j);
		weightedNodeList.at(j) = temp;
	}

	int k=0;
	for (unsigned int i=0 ; i< subnetRoot.size(); i++)
	{
		for (int j=0; j<subnetRoot.at(i).desiredSize; j++ )
		{
			if (   isAttachable(  subnetRoot.at(i),nodeRoot.at(weightedNodeList.at(k) )  )  )
			{
				attach(subnetRoot.at(i), nodeRoot.at(weightedNodeList.at(k)));
				k++;
			}
			else
				cout<<"subnet:"<<i<<"node:"<<weightedNodeList.at(k)<<"faield to attach"<<endl;
		}
	}
}

void network::reverseNodeOrder(void)
{ //assumes already sorted in increasing order wrt size
//using shallow copy only, subnetList should not have been allocated !!
	node temp;
	for (int i=0; i<NN-i-1; i++)
	{
		temp=nodeRoot.at(i);
		nodeRoot.at(i) = nodeRoot.at(NN-i-1);
		nodeRoot.at(NN-i-1)=temp;
	}
}

void network::mixSubnets(void)
{//shallow copy is safe since no dynamic memory is allocated yet
	subnet temp;
	for(int i=0; i<subnetCount; i++)
	{
		long long int random = rand()%subnetCount;
		temp=subnetRoot[random];
		subnetRoot[random]=subnetRoot[i];
		subnetRoot[i]=temp;
	}

	for(int i=0; i<subnetCount; i++)
	{
		subnetRoot[i].id=i;
	}

	logfile<<"subnets mixed"<<endl;
}

void network::mixNodes(void)
{//shallow copy is safe since no dynamic memory is allocated yet
	node temp;
	for(int i=0; i<NN; i++)
	{
		long long int random = rand()%NN;
		temp=nodeRoot[random];
		nodeRoot[random]=nodeRoot[i];
		nodeRoot[i]=temp;
	}
	for(int i=0; i<NN; i++)
	{
		nodeRoot[i].id=i;
	}
	logfile<<"nodes mixed"<<endl;
}

void network::connectMST(void)
{//Minimum Spanning Tree which generates at least one path between any two subnets
	int subnetIndex=1;
	for (int i=0;i<NN;i++)//for each node
	{
		if (nodeRoot[i].desiredAlias >= 2)//nodes with degree<2 does not involve in MST
		{
			subnetIndex--;//go back to the previous subnet to continue attaching giant component
			while ( nodeRoot[i].desiredAlias > nodeRoot[i].currentAlias )
			{
				if (subnetIndex >= subnetCount)
					break;

				attach( subnetRoot[subnetIndex] , nodeRoot[i] );
				subnetIndex++;
			}
		}
		if (subnetIndex >= subnetCount) //If done with subnets, exit.
			break;
	}
}

void network::connectRestofMST(void)
{
	//push the IDs of each subnet to a vector the number of available interfaces times
	vector<int> availableSubnets;
	for (int i=0; i<subnetCount; i++)
	{
		for (int j=0; j < subnetRoot[i].desiredSize-subnetRoot[i].currentSize; j++)
		{
			availableSubnets.push_back(i);
		}
	}

	//randomize the vector
	vector <int>::iterator it;
	for(it=availableSubnets.begin() ; it<availableSubnets.end() ; it++)
	{
		int randomIndex= rand()%availableSubnets.size();
		int temp=availableSubnets.at(randomIndex);
		availableSubnets.at(randomIndex)=*it;
		*it=temp;
	}

	//attach nodes to the subnets in the randomized available list
	for(int i=0; i<NN; i++)
	{
		while(nodeRoot[i].desiredAlias-nodeRoot[i].currentAlias > 0)
		{
			attach(subnetRoot[availableSubnets.back()], nodeRoot[i]);
			availableSubnets.pop_back();
		}
	}
	cout<<"remaining count"<<availableSubnets.size()<<endl;
}

void network::inOrderConnect(void)
{//iterates over the nodes. For each node iterates over the subnets in order.
	int firstAvailableSubnet=0;
	cout<<"subnetCount:"<<subnetCount<<endl;
	for(int i=0;i<NN;i++)
	{
		//cout<<i<<"\t"<<nodeRoot.at(i).desiredAlias<<endl;

		int k=firstAvailableSubnet;
		while(nodeRoot[i].currentAlias < nodeRoot.at(i).desiredAlias)
		{
			//cout<<"k:"<<k<<endl;
			if (   isAttachable(  subnetRoot.at(k),nodeRoot.at(i)  )  )
			{
				attach(subnetRoot.at(k), nodeRoot.at(i));
				//cout<<"-"<<endl;
			}
			else
				firstAvailableSubnet=k;
			k=(k+1)%subnetCount;
		}
	}
}

void network::randomConnect(void)
{
	long int disconnectCount=0;

	attach( subnetRoot[0] , nodeRoot[0] );//initialize giant component
	int i=1;
	int j=1;
	long int freeNodeInterfaces, freeSubnetInterfaces;
	while(   i<NN || j< subnetCount)
	{
		freeNodeInterfaces=0;
		freeSubnetInterfaces=0;

		for(int dummy=0;dummy<i;dummy++)
			freeNodeInterfaces+=nodeRoot[dummy].desiredAlias-nodeRoot[dummy].currentAlias;
		for(int dummy=0;dummy<j;dummy++)
				freeSubnetInterfaces+=subnetRoot[dummy].desiredSize-subnetRoot[dummy].currentSize;

//		cout<<i<<"\t"<<j<<endl;
		if(   i<NN && freeSubnetInterfaces>0    )
		{
//			cout<<"x";
			int temp= rand()%freeSubnetInterfaces;
			int k=0;
			for( ; k<subnetCount ; k++)
			{
				temp-=subnetRoot[k].desiredSize-subnetRoot[k].currentSize;
				if (temp<0)
					break;
			}
			attach(   subnetRoot[ k  ] , nodeRoot[i] );
			i++;
			freeNodeInterfaces+=nodeRoot[i].desiredAlias-1;
			freeSubnetInterfaces--;
		}
		if(j<subnetCount && freeNodeInterfaces>0   )
		{
//			cout<<"y";
			int temp = rand()%freeNodeInterfaces;
			int k=0;
			for(; k<NN ; k++)
			{
				temp-=nodeRoot[k].desiredAlias-nodeRoot[k].currentAlias;
				if(temp<0)
					break;
			}
			attach(subnetRoot[j],nodeRoot[k]);
			j++;
			freeSubnetInterfaces+=subnetRoot[j].desiredSize-1;
			freeNodeInterfaces--;
		}
		if (   freeNodeInterfaces==0  && freeSubnetInterfaces==0   )
		{
			disconnectCount++;
	//		cout<<"z";
			attach(subnetRoot[j] , nodeRoot[i]);
			i++;
			j++;
		}
	}
	cout<<freeNodeInterfaces<<"\t"<<freeSubnetInterfaces<<endl;
	long int sum=0;
	long int sum2=0;
	for(int i=0;i<NN;i++)
	{
		sum+=nodeRoot[i].desiredAlias;
		sum2+=nodeRoot[i].currentAlias;
	}
	cout<<"sum, sum2     "<<sum<<"\t"<<sum2<<endl;
	sum=0;
	sum2=0;
	for(int i=0;i<subnetCount;i++)
	{
		sum+=subnetRoot[i].desiredSize;
		sum2+=subnetRoot[i].currentSize;
	}
	cout<<"sum, sum2:    "<<sum<<"\t"<<sum2<<endl;
	cout<<disconnectCount<<endl;
}

void network::randomAttach(void)
{
	long int freeInterfaceCount=0;
	for (int i=0;i<NN;i++)//calculate number of free Interfaces in the network
	{
		freeInterfaceCount+=nodeRoot[i].desiredAlias-nodeRoot[i].currentAlias;
	}
	int * nodeInterfaces=new int [freeInterfaceCount];
	int * subnetInterfaces=new int [freeInterfaceCount];

	int index=0;
	for(int i=0; i<NN; i++)//create a list of size= freeInterfaceCount
	{                      //Write the ID of each node depending as many free interfaces as each have.
		for(int j=0; j< (nodeRoot[i].desiredAlias-nodeRoot[i].currentAlias) ;j++)
		{
			nodeInterfaces[index]=nodeRoot[i].id;
			index++;
		}
	}
	index=0;
	for(int i=0; i<subnetCount;i++)
	{//Write the ID of subnet node depending as many free interfaces as each have.
		for (int j=0; j<subnetRoot[i].desiredSize- subnetRoot[i].currentSize ; j++)
		{
			subnetInterfaces[index]=i;
			index++;
		}
	}
	for(int i=0; i<freeInterfaceCount;i++)//mix-randomize both lists
	{
		int temp=0;

		int random=rand()%freeInterfaceCount;
		temp=nodeInterfaces[i];
		nodeInterfaces[i]=nodeInterfaces[random];
		nodeInterfaces[random]=temp;

		random=rand()%freeInterfaceCount;
		temp=subnetInterfaces[i];
		subnetInterfaces[i]=subnetInterfaces[random];
		subnetInterfaces[random]=temp;
	}
	for (int i=0;i<freeInterfaceCount;i++)
	{
		attach(subnetRoot[   subnetInterfaces[i]  ], nodeRoot[   nodeInterfaces[i]   ] );
	}
}

bool network::isAttachable(subnet &s, node &n)
{
	//cout<<"-"<<endl;
	if ( n.currentAlias > n.desiredAlias )	//check if node is in problem
	{
		cerr<<"node alias overflow"<<endl;
		return false;
	}
	if (n.currentAlias == n.desiredAlias )	//check if node is full
		return false;
	if (s.currentSize > s.desiredSize)		//check if subnet is in problem
	{
		//cout<<s.id<<"\t"<<n.id<<"\tsubnet size overflow"<<endl;
		//cout<<"ids:\t"<<s.currentSize<<"\t"<<s.desiredSize<<endl;
		return false;
	}
	if (s.currentSize == s.desiredSize)		//check if subnet is full
		return false;
	return true;
}

void network::attach(subnet &s, node &n)
{
	//cout<<"+"<<endl;
	s.nodeList[s.currentSize]=&n;
	s.currentSize++;
	n.subnetList[n.currentAlias]=&s;
	n.currentAlias++;
	if (n.currentAlias>n.desiredAlias)
		cout<<"ERROR! Alias Count Overflow.NodeID:"<<n.id<<endl;
	if (s.currentSize>s.desiredSize)
		cout<<"ERROR! Subnet Size Overflow. SubnetID"<<s.id<<endl;
}

void network::randomizeNetwork()
{
	long long int randomNode=0, randomSubnet=0;
	for( int i=0;i<NN; i++)
	{
		for(int j=0; j<nodeRoot[i].currentAlias;j++)
		{
			randomNode=rand() % NN;
			randomSubnet= rand() % nodeRoot[randomNode].currentAlias;
			if (   isSwappable(nodeRoot[i] , j , nodeRoot[randomNode], randomSubnet)   )
			;
			else
			{
				j--;
			}
		}
	}
}

void network::calculateDegree(void)
{
	for (int i=0; i<NN; i++)
	{
		nodeRoot[i].currentDegree=0;
		for(  int j=0  ;    j<nodeRoot[i].currentAlias   ; j++ )
		{
			nodeRoot[i].currentDegree+= nodeRoot[i].subnetList[j]->currentSize-1;
		}
//	cout<<i<<"   "<<nodeRoot[i].currentDegree<<endl;
	}
}

bool network::isSwappable(node &n1,int s1, node &n2, int s2)//This function checks if the argument subnets can be swapped
{// if true, it swaps directly.
	for(int i=0; i<n1.currentAlias;i++)//check if swapping causes another duplicate in n1
	{
		if(n1.subnetList[i]->id  ==  n2.subnetList[s2]->id)
			return 0;
	}
	for(int i=0; i<n2.currentAlias;i++)//check if swapping causes another duplicate in n2
	{
		if(n2.subnetList[i]->id  ==  n1.subnetList[s1]->id)
			return 0;
	}

	//otherwise swap them.

	for (int i=0; i<n1.subnetList[s1]->currentSize;i++)//update the pointers in subnet structures.
	{
		if (n1.subnetList[s1]->nodeList[i]->id  ==  n1.id )
		{
			n1.subnetList[s1]->nodeList[i]= &n2;
			break;
		}
	}
	for (int i=0; i<n2.subnetList[s2]->currentSize; i++ )
	{
		if (n2.subnetList[s2]->nodeList[i]->id  ==  n2.id )
		{
			n2.subnetList[s2]->nodeList[i]= &n1;
			break;
		}
	}
	//update currentDegree of the nodes
	n1.currentDegree-=n1.subnetList[s1]->currentSize-1;
	n1.currentDegree+=n2.subnetList[s2]->currentSize-1;
	n2.currentDegree-=n2.subnetList[s2]->currentSize-1;
	n2.currentDegree+=n1.subnetList[s1]->currentSize-1;

	//finally swap the subnet pointers in the node structures
	subnet * tempSubnet=0;
	tempSubnet=n1.subnetList[s1];
	n1.subnetList[s1]=n2.subnetList[s2];
	n2.subnetList[s2]=tempSubnet;

	return 1;
}

void network::calculateCurrentNDD(void)
{
	currentNDD.clear();
	currentNDD.resize(maxDistributionIndex,0);
	for (int i=0;i<NN;i++)
		if (nodeRoot.at(i).currentDegree > maxDistributionIndex)
			cerr<<"node exceed maxDistributionIndex"<<i<<":  "<<nodeRoot.at(i).currentDegree<<"\tignoring"<<endl;
	for (int i=0;i<maxDistributionIndex;i++)
		currentNDD.at(i)=0;
	for(int i=0;i<NN;i++)
	{
		if (unsigned (nodeRoot[i].currentDegree) < currentNDD.size())
			currentNDD.at(   nodeRoot[i].currentDegree    )++;
	}
}

void network::printCurrentNDD(void)
{
	ofstream NDDfile;
	NDDfile.open("NDD.txt");
	cout<<"Printing Current NDD"<<endl;
	for (int i=0;i<maxDistributionIndex;i++)
		if (currentNDD[i]>0)
			NDDfile<<i<<"\t"<<currentNDD[i]<<endl;

	NDDfile.close();
}

void network::printNDD(void)
{
	cout<<"Printing Desired NDD:"<<endl;
	for (int i=0;i<maxDistributionIndex;i++)
	{
		cout<<i<<"   "<<NDD[i]<<endl;
	}
}

void network::powerifyNDD(void)
{
	long long int node1=0,node2=0,subnet1=0,subnet2=0;
	for(int i=0;i<10000;i++)
	{
		node1=rand()%NN;
		node2=rand()%NN;
		subnet1=rand() % (nodeRoot[node1].currentAlias);
		subnet2=rand() % (nodeRoot[node2].currentAlias);

		if (   isConverge(nodeRoot[node1],subnet1,nodeRoot[node2],subnet2)  )
		{
			if (   isSwappable(nodeRoot[node1],subnet1,nodeRoot[node2],subnet2)   )
			{//update currentNDD

				currentNDD[ nodeRoot[node1].currentDegree ]++;
				currentNDD[ nodeRoot[node2].currentDegree ]++;
				currentNDD[ nodeRoot[node1].currentDegree - nodeRoot[node1].subnetList[subnet1]->currentSize + nodeRoot[node2].subnetList[subnet2]->currentSize]--;
				currentNDD[ nodeRoot[node2].currentDegree - nodeRoot[node2].subnetList[subnet2]->currentSize + nodeRoot[node1].subnetList[subnet1]->currentSize]--;
			}
		}
	}
}

bool network::isConverge(node &n1, int s1, node &n2, int s2)
{
	long long signed int difference = (n1.subnetList[s1]->currentSize-1)   -   (n2.subnetList[s2]->currentSize-1);
	if (difference ==  0 )
		return 0;

	long long int degree1=0, degree2=0;

	degree1= n1.currentDegree;
	degree2= n2.currentDegree;

	float effect1= ( currentNDD[degree1] > NDD[degree1] ) ? 1:-1;
	float effect2= ( currentNDD[degree2] > NDD[degree2] ) ? 1:-1;
	float effect3= (NDD[degree1-difference]   >  currentNDD[degree1-difference])  ? 1:-1;
	float effect4= (NDD[degree2+difference]   >  currentNDD[degree2+difference])  ? 1:-1;

//	float effect1=float(currentNDD[degree1]-NDD[degree1]) /  float(NDD[degree1]);
//	float effect2=float(currentNDD[degree2]-NDD[degree2])/float(NDD[degree2]);
//	float effect3=0;  //we have to consider cases where denominator is zero.
//	if 	(NDD[degree1-difference]  == 0 )
//		effect3= -1;
//	else
//		effect3= float(NDD[degree1-difference]-currentNDD[degree1-difference])/float(NDD[degree1-difference]);
//	float effect4=0;
//	if  (NDD[degree2+difference]  ==  0)
//		effect4=0;
//	else
//		effect4=float(NDD[degree2+difference]-currentNDD[degree2+difference])/float(NDD[degree2+difference]);

	float netEffect=  effect1+effect2+effect3+effect4;

//    cout<<"currentNDD[degree1]:"<<currentNDD[degree1]<<"  NDD[degree1] "<<NDD[degree1]<<endl;
//    cout<<"degree1:   "<<degree1<<"   degree2:  "<<degree2<<"   difference:   "<<difference<<"  effects:"<<effect1<<"  "<<effect2<<"  "<<effect3<<"   "<<effect4<<"  neteffect:  "<<netEffect<<endl;

    if( netEffect > 0)
    	return 1;
    else return 0;

}

void network::printMap(void)
{
	ofstream subnetMap;
	subnetMap.open("subnetMap.txt");

	for (int i=0; i<subnetCount;i++)
	{
		for(int j=0;j<subnetRoot[i].currentSize; j++)
		{
			subnetMap<<subnetRoot[i].nodeList[j]->id;
			subnetMap<<"\t";
		}
		subnetMap<<endl;
	}
	subnetMap.close();
}

void network::printSSD(void)
{
	ofstream SSDfile;
	SSDfile.open("SSD.txt");
	cout<<"printing subnet size distribution"<<endl;
	for(int i=0;i<maxDistributionIndex;i++)
	{
		if (SSD.at(i)>0)
			SSDfile<<i<<"  "<<SSD[i]<<endl;
	}
	SSDfile.close();
}

void network::printNAD(void)
{
	cout<<"printing node alias distribution"<<endl;
	ofstream NADfile;
	NADfile.open("NAD.txt");
	for (int i=0;i<maxDistributionIndex;i++)
	{
		if (NAD.at(i)>0)
			NADfile<<i<<"  "<<NAD[i]<<endl;
	}
	NADfile.close();
}

void network::printEdgeList(void)
{
	ofstream EdgeListFile;
	EdgeListFile.open("edgeList.txt");

	for(int i=0;i<subnetCount; i++)
	{
		for (int j=0; j< subnetRoot[i].currentSize - 1 ;  j++)
		{
			for (int k=j+1; k< subnetRoot[i].currentSize ; k++)
			{
				EdgeListFile<<subnetRoot[i].nodeList[j]->id<<"\t"<<subnetRoot[i].nodeList[k]->id<<endl;
			}
		}
	}
	EdgeListFile.close();
}
