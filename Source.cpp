#include<stdio.h>
#include<fstream>
#include<sstream>
#include<iostream>
#include<string>
#include<vector>
#include<map>
#include<math.h>
#include<algorithm>
//#include<conio.h>

using namespace std;

int inside = 0;

class BWT
{
public:
	vector<char> BWTString;
	vector<int> LCPVal;
	vector<unsigned int> componentIds;
	vector<int> bucketSizes;
	int numOfBuckets;
	int compSize;
	int kmerLength;
	int minLCPLength;
	int numEltsInEachBucket; // Bucket size.

	string origString;
	vector<string> LCPArray;
	void QuickSort(int start, int n);
	int Partition(int start, int n);
	BWT()
	{
		//bucketSizes = NULL;
	}
	void swap(unsigned int &s1, unsigned int &s2);
	void findSuperMaximalRepeats(string outputFName);
	void findBWT();
	void findLCPArray();
	void fillUpComponentIds(int bucketId);
	unsigned int convertCharacter(char);
	unsigned int getKmerMask();
};
void BWT::findBWT()
{
	for (int i = 0; i < componentIds.size(); i++)
	{
		if (componentIds[i]>0)
			BWTString.push_back(origString[componentIds[i] - 1]);
		else
			BWTString.push_back(origString[origString.size() - 1]);
	}
}

void BWT::findLCPArray()
{
	LCPArray.push_back("\0");
	LCPVal.push_back(-1);
	for (int i = 1; i < componentIds.size(); i++)
	{
		string s1 = origString.substr(componentIds[i - 1]);
		string s2 = origString.substr(componentIds[i]);
		int j = 0;
		string LCPStr = "";
		while (s1[j] == s2[j])
		{
			LCPStr += s1[j];
			j++;
		}
		LCPStr[j] += '\0';
		LCPArray.push_back(LCPStr);
		LCPVal.push_back(j);
	}
}

void BWT::QuickSort(int start, int n)
{
	if (n > 1){
		int pivotIndex = Partition(start, n);
		int n1 = pivotIndex - start;
		int n2 = n - n1 - 1;
		QuickSort(start, n1);
		QuickSort(pivotIndex + 1, n2);
		inside--;
	}
}

int BWT::Partition(int start, int n)
{
	inside++;
	string pivot = origString.substr(componentIds[start], compSize);
	int i = start + 1;
	int j = start + n - 1;
	while (i <= j){
		string iComp = origString.substr(componentIds[i], compSize);
		string jComp = origString.substr(componentIds[j], compSize);
		string iPivot = pivot, jPivot = pivot;
		int pivotPoint = componentIds[start];
		int iInc = componentIds[i], jInc = componentIds[j];
		while (iComp == iPivot)
		{
			iInc += compSize;
			pivotPoint += compSize;

			iComp = origString.substr(iInc, compSize);
			iPivot = origString.substr(pivotPoint, compSize);
		}
		//iPivot = pivot;
		pivotPoint = componentIds[start];
		while (jComp == jPivot)
		{
			jInc += compSize;
			pivotPoint += compSize;
			jPivot = origString.substr(pivotPoint, compSize);
			jComp = origString.substr(jInc, compSize);
		}
		if (iComp < iPivot)
			i++;
		else if (jComp > jPivot)
			j--;
		else{
			swap(componentIds[i], componentIds[j]);
			i++;
			j--;
		}//else
	}//while
	swap(componentIds[start], componentIds[j]);
	return j;
}//end Partition

unsigned int BWT::convertCharacter(char c)
{
	switch (c)
	{
	case 'A':
		return 0x00;
		break;
	case 'N':
		return 0x00;
		break;
	case 'C':
		return 0x01;
		break;
	case 'G':
		return 0x02;
		break;
	case 'T':
		return 0x03;
		break;
	default:
		return 0xffffffff;
		break;
	}
}

void BWT::swap(unsigned int &s1, unsigned int &s2)
{
	unsigned int temp = s1;
	s1 = s2;
	s2 = temp;
}

void BWT::fillUpComponentIds(int bucketId)
{
	componentIds.clear();
	LCPVal.clear();
	LCPArray.clear();
	BWTString.clear();
	unsigned int kmer = 0;
	unsigned int maxKmer = 0;
	unsigned int mask = getKmerMask();
	int minIndex = bucketId*numEltsInEachBucket;
	int maxIndex = (bucketId + 1)*numEltsInEachBucket - 1;
	for (int i = 0; i < origString.length() - kmerLength; i++)
	{
		kmer = 0;
		for (int j = 0; j < kmerLength; j++)
		{
			kmer <<= 2;
			kmer &= mask;
			unsigned int convChar = convertCharacter(origString[i+j]);
			if (convChar == 0xffffffff)
			{
				cout << "i=" << i << endl;
				cout << "origString[" << i+j << "]=" << origString[i+j] << endl;
				throw "Invalid character encountered!!";
			}
			kmer |= convChar;
			if (kmer > maxKmer)
			{
				maxKmer = kmer;
			}
		}
		if (bucketId == 0)
		{
			int buckSizeIndex = kmer / numEltsInEachBucket;
			bucketSizes[buckSizeIndex]++;
		}
		if(kmer >= minIndex && kmer<maxIndex)
			componentIds.push_back(i);
		

	}
}

unsigned int BWT::getKmerMask()
{
	unsigned int x = 0xffffffff;
	x <<= (kmerLength*2);
	x = ~x;
	return x;
}

void BWT::findSuperMaximalRepeats(string outputFName)
{
	ofstream fout;
	fout.open(outputFName, ios::app);
	bool currentUp = false, currDown = false;
	int startInt = -1, endInt = -1;

	for (int i = 0; i < LCPVal.size() - 1; i++)
	{
		//if (!currentUp && bwt.LCPVal[i+1] < 3)
		//break;
		if (i + 1 != LCPVal.size() && LCPVal[i]<LCPVal[i + 1])
		{
			currentUp = true;
			startInt = i;
			endInt = i + 1;
		}
		if (currentUp)
		{
			if (LCPVal[i] == LCPVal[i + 1])
			{
				endInt = i + 1;
			}
			else if (LCPVal[i] > LCPVal[i + 1])
			{
				currentUp = false;
				currDown = true;
			}
		}
		if (!currentUp && currDown)
		{
			//put stint and endint in file.
			if (endInt - startInt + 1 <= 4 && LCPVal[endInt]>minLCPLength) // possiblity for supermaximal repeat 4 and 15
			{
				// check for pairwise distinct
				bool pairWiseDistinct = true;
				for (int j = startInt; j < endInt; j++)
				{
					for (int k = j + 1; k <= endInt; k++)
					{
						if (BWTString[j] == BWTString[k])
						{
							pairWiseDistinct = false;
							break;
						}
					}
				}

				if (fout.is_open() && pairWiseDistinct)
				{
					string superMaxRep = origString.substr(componentIds[startInt], LCPVal[endInt]);
					fout << componentIds[startInt] << "\t" << LCPVal[endInt] <<"\t"<<endInt-startInt + 1<< endl;
					cout << componentIds[startInt] << "\t" << LCPVal[endInt] << "\t" << endInt - startInt + 1<<"\t"<<superMaxRep<< endl;
					pairWiseDistinct = false;
					currentUp = false;
					currDown = false;
				}
			}
		}
	}
	fout.close();
}

int main(int argc, char *argv[])
{
	//start: Assign all the command line arguments
	BWT bwt;
	if (argc != 6)
	{
		cout << "Invalid number of arguments!";
			return 1;
	}
	bwt.kmerLength = atoi(argv[1]);
	bwt.minLCPLength = atoi(argv[2]);
	bwt.compSize = atoi(argv[3]);
	string inputFName = argv[4];
	string outputFName = argv[5];
	// end: Assign all the command line arguments
	if (bwt.kmerLength == 0 || bwt.minLCPLength == 0 || bwt.compSize == 0)
	{
		cout << "Error : Invalid kmer length or component size of minimum LCP length!" << endl;
		return 1;
	}
	ifstream myFile(inputFName.c_str());
	stringstream ss;
	string s;
	if (myFile.is_open())
	{
		string line;
		int i = 0;
		int t = 0;
		while (getline(myFile, line))
		{
			if (line[0] != '>')
			{
				line.erase(line.size() - 1);
				for (int k = 0; k<line.size(); k++)
				{
					if (line[k] == 'A' || line[k] == 'C' || line[k] == 'G' || line[k] == 'T'){
						ss << line[k];
					}
				}
			}
		}
		myFile.close();
	}
	else
	{
		cout << "Error : Input file does not exist";
		return 1;
	}
	s = ss.str() + "$";
	//cout<<s;
	/*fstream fout3;
	fout3.open("lcpval.txt", ios::out);
	fout3 << s;
	fout3.close();*/
	bwt.origString = s;
	bwt.numOfBuckets = pow(4,bwt.kmerLength);
	bwt.numEltsInEachBucket = 1024;

	//Start: Initialize bucketSizes Array
	vector<int> tempBuckSizes(bwt.numOfBuckets, 0);
	bwt.bucketSizes = tempBuckSizes;
	//end: Inistialize bucketsizes array
	int cnt1 = 0;
	try
	{
		fstream fout2;
		fout2.open("lcpval.txt", ios::out);
		fstream fout;
		fout.open(outputFName.c_str(), ios::out);
		fout.close();
		
		for (int i = 0; i < bwt.numOfBuckets; i++)
		{
			// fill up the components
			if (i==0 || bwt.bucketSizes[i] != 0)
			{
				bwt.fillUpComponentIds(i);

				// Sort the components
				bwt.QuickSort(0, bwt.componentIds.size());
				// find BWT from sorted component ids
				bwt.findBWT();
				// find LCPs from sorted component ids
				bwt.findLCPArray();
				fout2 << "--------LCPVAL ARRAY   BUCKET   " << i << "   -------" << endl;
				for (int j = 0; j < bwt.LCPVal.size() && j < bwt.componentIds.size() && j < bwt.BWTString.size(); j++)
				{
					fout2 << "pos : " << j << "\tIndex : " << bwt.componentIds[j] << "\t LCP val :" << bwt.LCPVal[j] << "\t BWT : " << bwt.BWTString[j] << endl;
				}
				// Start : compute Super maximal repeats 
				bwt.findSuperMaximalRepeats(outputFName);
				// End : compute Super maximal repeats
				cnt1++;
			}
		}
	}
	catch (const char* msg)
	{
		cout << "Exception : " << msg << endl;
	}
	cout << "over";
	return 0;

}	