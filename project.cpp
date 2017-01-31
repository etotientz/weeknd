#include <bits/stdc++.h>
#define ll long long
using namespace std;

double Eelc = 0.0000000001;
double Eemp = 0.0000000000001;
double Eagg = 0.0000000000000001;
double lengthOfPacket = 1000, n;

struct node {
	double xcor, ycor, transmissionRange;
	int id, trackNumber, clusterNumber, sectorNumber, parentId;
	double residualEnergy = 0.5;
	bool isClusterHead = false;
	bool istrackHead = false;
	vector <node> neighbourhoodNodes;
};

double dist(struct node a, struct node b) {
	return sqrt((a.xcor - b.xcor) * (a.xcor - b.xcor) + (a.ycor - b.ycor) * (a.ycor - b.ycor));
}

bool mySort(node a, node b) {
	return a.clusterNumber < b.clusterNumber;	
}

bool mySort2(node a, node b) {
	if (a.sectorNumber < b.sectorNumber) {
		return true;
	} else if (a.sectorNumber == b.sectorNumber) {
		return (a.trackNumber > b.trackNumber);
	} else {
		return false;
	}
}

vector<node>  getNodes(double xc, double yc, double r, int N) {
	vector<node> points;
	set<double> xAxis;
	set<double> yAxis;

	ll i = 0;
	while (i < N) {
		double lower_bound = 0;
	    double upper_bound = 2 * r + 1;
	    std::random_device rd;
    	std::mt19937 gen(rd());
	    std::uniform_real_distribution<> dis(lower_bound,upper_bound);
	    double x = dis(gen) - r;
	    double y = dis(gen) - r;
		if (xAxis.find(x) != xAxis.end() && yAxis.find(y) != yAxis.end()) {
			continue;
		}
		double temp1 = abs(x - xc) * abs(x - xc);
		double temp2 = abs(y - yc) * abs(y - yc);

		if((temp1 + temp2) <= (r * r)){
			node temp;
			temp.xcor = x;
			temp.ycor = y;
			xAxis.insert(x);
			yAxis.insert(y);
			temp.id = i + 1;
			points.push_back(temp);
			i++;
		}
	}
	return points;
}

void divideInToTracks(vector<node> &nodeList, int noOfTracks, int r) {
	double temp = (double)(r / noOfTracks);

	for (int i = 0; i < nodeList.size(); i++) {
		double num = (abs(nodeList[i].xcor) * abs(nodeList[i].xcor)) + (abs(nodeList[i].ycor) * abs(nodeList[i].ycor));
		double distanceFromBS = sqrt(num);
		//cout<<nodeList[i].xcor<<" "<<nodeList[i].ycor<<" "<<num<<" "<<distanceFromBS<<"\n";
		double zero = 0.0;
		if (fmod(distanceFromBS, temp) == zero) {
			nodeList[i].trackNumber = (int)(distanceFromBS / temp);
		} else {
			nodeList[i].trackNumber = (int)(distanceFromBS / temp) + 1;
		}
	}
}

int main() {
	int noOfTracks, noOfSectors, N;
	double radius, areaOfASegment;

	cout<<"Enter radius of circle: ";
	cin>>radius;
	cout<<"Enter the total number of nodes: ";
	cin>>N;
	cout<<"Enter no of tracks and no of sectors respectively: ";
	cin>>noOfTracks>>noOfSectors;
	
	double sectorAngle = 360 / noOfSectors;
	areaOfASegment = ((M_PI * (sectorAngle / 180)) / 2) * radius * radius;
	cout<<"Area of each and every Sector: "<<areaOfASegment<<"\n\n";

	vector <node> nodeList = getNodes(0, 0, radius, N);
	
	divideInToTracks(nodeList, noOfTracks, radius);


	for (int i = 0; i < nodeList.size(); i++) {

		if (nodeList[i].ycor == 0 && nodeList[i].xcor > 0) {
			nodeList[i].sectorNumber = 1;
			nodeList[i].clusterNumber = 1;
		} else {
			double tempAngle = atan2(nodeList[i].ycor, nodeList[i].xcor) * (180 / M_PI);
			if (tempAngle < 0)
				tempAngle = 360 + tempAngle;
			nodeList[i].sectorNumber = tempAngle / sectorAngle + 1;
			//cout<<nodeList[i].xcor<<" "<<nodeList[i].ycor<<" "<<tempAngle<<"\n";

			if (fmod(tempAngle, sectorAngle) == 0) {
				nodeList[i].clusterNumber = 1;
			} else {
				double temp = sectorAngle / (2 * nodeList[i].trackNumber - 1);
				double temp2 = ((nodeList[i].sectorNumber - 1) * sectorAngle);
				int ans = 0;
				while (temp2 < tempAngle) {
					ans++;
					temp2 += temp;
				}
				nodeList[i].clusterNumber = ans;
			}
		}
	}
/*
	for (int i = 0; i < nodeList.size(); i++) {
		cout<<nodeList[i].xcor<<" "<<nodeList[i].ycor<<" ";
		cout<<nodeList[i].sectorNumber<<" "<<nodeList[i].trackNumber<<" "<<nodeList[i].clusterNumber<<"\n";
	}
*/

	// changing enery inside every cluster
	for (int i = 1; i <= noOfSectors; i++) {
		for (int j = 1; j <= noOfTracks; j++) {
			for (int k = 1; k <= (2 * j - 1); k++) {
				
				double maxEnergy = INT_MIN;
				int maxId = -1;
				for (int l = 0; l < nodeList.size(); l++) {
					if (nodeList[l].sectorNumber == i && nodeList[l].trackNumber == j &&
						nodeList[l].clusterNumber == k) {
						if (nodeList[l].residualEnergy > maxEnergy) {
							maxEnergy = nodeList[l].residualEnergy;
							maxId = nodeList[l].id;
						}
					}
				}
				for (int l = 0; l < nodeList.size(); l++) {
					if (nodeList[l].id == maxId) {
						nodeList[l].isClusterHead = true;
					}
				}

				int count = 0;
				for (int l = 0; l < nodeList.size(); l++) {
					if (nodeList[l].sectorNumber == i && nodeList[l].trackNumber == j &&
						nodeList[l].clusterNumber == k) {
						if (!nodeList[l].isClusterHead) {
							count++;
							nodeList[l].parentId = maxId;
							// subtracting the energy of non-clusterhead node
							nodeList[l].residualEnergy -= (Eelc + Eemp * (radius / noOfTracks)) * lengthOfPacket;
						} else {
							nodeList[l].residualEnergy -= Eelc * lengthOfPacket * count;
						}
					}
				}
			}
		}
	}


	vector <node> trackHeads;

	//choose track head
	for (int i = 1; i <= noOfSectors; i++) {
		for (int j = 1; j <= noOfTracks; j++) {
			bool found = false;
			for (int l = 0; l < nodeList.size(); l++) {
				if (nodeList[l].isClusterHead && nodeList[l].sectorNumber == i 
					&& nodeList[l].trackNumber == j && nodeList[l].clusterNumber == nodeList[l].trackNumber) {
					//cout<<nodeList[l].sectorNumber<<" "<<nodeList[l].trackNumber<<" "<<nodeList[l].xcor<<" "<<nodeList[l].ycor<<endl;
					found = true;
					nodeList[l].istrackHead = true;
					//trackHeads.push_back(nodeList[l]);
					break;
				}
			}
			if (!found) {
				for (int l = 0; l < nodeList.size(); l++) {
					if (nodeList[l].isClusterHead && nodeList[l].sectorNumber == i 
						&& nodeList[l].trackNumber == j) {
						//cout<<nodeList[l].sectorNumber<<" "<<nodeList[l].trackNumber<<" "<<nodeList[l].xcor<<" "<<nodeList[l].ycor<<endl;
						nodeList[l].istrackHead = true;
						//trackHeads.push_back(nodeList[l]);
						break;
					}
				}
			}
		}
	}

	int p = 1;

	// it will store all the tracks
	vector<vector<node> >storage(noOfSectors * noOfTracks + 1);

	for (int i = 1; i <= noOfSectors; i++) {
		for (int j = 1; j <= noOfTracks; j++) {
			//cout<<"sector number = "<<i<<" track number = "<<j<<"\n";
			for (int l = 0; l < nodeList.size(); l++) {
				if (nodeList[l].isClusterHead && nodeList[l].sectorNumber == i 
					&& nodeList[l].trackNumber == j) {
					storage[p].push_back(nodeList[l]);
				}
			}
			p++;
		}
	}

	
	for (int i = 1; i < storage.size(); i++) {
		sort(storage[i].begin(), storage[i].end(), mySort);
	}

	for (int i = 1; i < storage.size(); i++) {
		cout<<"Sector Number = "<<storage[i][0].sectorNumber<<"  Track Number = "<<storage[i][0].trackNumber<<"\n";
		for (int j = 0; j < storage[i].size(); j++) {
			if (storage[i][j].istrackHead)
				cout<<"Track Head is:";
			cout<<storage[i][j].xcor<<" "<<storage[i][j].ycor<<"  ";
		}
		cout<<"\n\n";
	}


	//sending energy from clusterheads to trackheads in every track(track in middle)
	for (int i = 1; i < storage.size(); i++) {
		n = 1;
		int l = 0, h = storage[i].size() - 1;
		while ((l + 1) < storage[i][l].trackNumber) {
			double Dab = sqrt((storage[i][l + 1].xcor - storage[i][l].xcor) * (storage[i][l + 1].xcor - storage[i][l].xcor)
			 + (storage[i][l + 1].ycor - storage[i][l].ycor) * (storage[i][l + 1].ycor - storage[i][l].ycor));
			storage[i][l].residualEnergy -= (Eelc + Eemp * Dab * Dab) * n * lengthOfPacket;
			storage[i][l + 1].residualEnergy -= (Eelc) * lengthOfPacket * n;

			//storage[i][l + 1].residualEnergy -= (storage[i].size() * Eagg);

			storage[i][l].parentId = storage[i][l + 1].id;
			n++;
			l++;
		}
		n = 1;
		while ((h + 1) > storage[i][l].trackNumber) {
			double Dab = sqrt((storage[i][h - 1].xcor - storage[i][h].xcor) * (storage[i][h - 1].xcor - storage[i][h].xcor)
			 + (storage[i][h - 1].ycor - storage[i][h].ycor) * (storage[i][h - 1].ycor - storage[i][h].ycor));
			storage[i][h].residualEnergy -= (Eelc + Eemp * Dab * Dab) * n * lengthOfPacket;
			storage[i][h - 1].residualEnergy -= (Eelc) * lengthOfPacket * n;

			//storage[i][h - 1].residualEnergy -= (storage[i].size() * Eagg);

			storage[i][h].parentId = storage[i][h - 1].id;
			n++;
			h--;
		}

		//Aggregation Energy on track head
		storage[i][storage[i].size() / 2].residualEnergy -= (storage[i].size() * Eagg);

		double distanceOfTrackHeadFromBS = sqrt((abs(storage[i][storage[i].size() / 2].xcor) * abs(storage[i][storage[i].size() / 2].xcor))
		 + (abs(storage[i][storage[i].size() / 2].ycor) * abs(storage[i][storage[i].size() / 2].ycor)));

		//Transmission Range
		storage[i][storage[i].size() / 2].transmissionRange = 2 * (radius / noOfTracks);
		cout<<storage[i][storage[i].size() / 2].transmissionRange<<"\n";
	}


	for (int i = 1; i < storage.size(); i++) {
		trackHeads.push_back(storage[i][storage[i].size() / 2]);
	}

	sort(trackHeads.begin(), trackHeads.end(), mySort2);


	for (int i = 0; i < trackHeads.size(); i ++) {
		cout<<trackHeads[i].trackNumber<<" "<<trackHeads[i].sectorNumber<<"\n";
	}

	for (int i = 0; i < trackHeads.size(); i += noOfTracks) {
		//cout<<trackHeads[i].transmissionRange<<"\n";
		for (int j = i; j < i + noOfTracks; j++) {
			for (int k = j + 1; k < i + noOfTracks; k++) {
				cout<<dist(trackHeads[j], trackHeads[k])<<" k "<<trackHeads[j].transmissionRange<<"\n";
				if (dist(trackHeads[j], trackHeads[k]) < trackHeads[j].transmissionRange) {
					trackHeads[j].neighbourhoodNodes.push_back(trackHeads[k]);
					cout<<"found for "<<trackHeads[j].id<<"\n";
				}
			}
		}

	}

	double tempResidualA, tempResidualB;
	for (int i = 0; i < trackHeads.size(); i++) {
		//cout<<"grgref\n";
		double tempMax = INT_MIN;
		int tempParentId;
		for (int j = 0; j < trackHeads[i].neighbourhoodNodes.size(); j++) {
			tempResidualA = trackHeads[i].residualEnergy;
			tempResidualB = trackHeads[i].neighbourhoodNodes[j].residualEnergy;

			tempResidualA -= (Eelc + Eemp * dist(trackHeads[i], trackHeads[i].neighbourhoodNodes[j]) * 
								dist(trackHeads[i], trackHeads[i].neighbourhoodNodes[j]));
			tempResidualB -= (lengthOfPacket * Eelc);

			if (tempMax < min(tempResidualA, tempResidualB)) {
				tempMax = min(tempResidualA, tempResidualB);
				tempParentId = trackHeads[i].neighbourhoodNodes[j].id;
			}
		}
		trackHeads[i].parentId = tempParentId;
	}

	for (int i = 0; i < trackHeads.size(); i++) {
		cout<<trackHeads[i].id<<" "<<trackHeads[i].neighbourhoodNodes.size()<<"\n";
	}

	for (int i = 0; i < trackHeads.size(); i++) {
		cout<<"parent of track head in sector:"<<trackHeads[i].sectorNumber<<" and track number:"<<
		trackHeads[i].trackNumber<<" with id:"<<trackHeads[i].id<<" is: "<<trackHeads[i].parentId<<"\n";
	}


	return 0;
}