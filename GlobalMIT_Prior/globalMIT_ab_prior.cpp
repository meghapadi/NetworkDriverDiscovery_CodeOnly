//////////////////////////////////////////////////////////
//                                                      //
//     globalMIT algorithm for DBN Learning             //
//     Author          : Vinh Nguyen                    //
//                                                      //
//	   Multi time series version						//
//     Project Members :					            //
//                                                      //
//     Last Updated    : 3 Feb, 2010                    //
//     Contact info    : n.x.vinh@unsw.edu.au           //
//                       vinh.nguyen@monash.edu         //
//                       vthesniper@yahoo.com           //
//                                                      //
//////////////////////////////////////////////////////////




#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <time.h>
#include <limits>
#include <vector>
#include <string>
#include "combination.h"

using namespace std;
using namespace stdcomb;
template<class BidIt>
void display(BidIt begin,BidIt end)
{
  for (BidIt it=begin;it!=end;++it)
      cout<<*it<<" ";
  cout<<endl;
}


typedef unsigned long long int UL;
double optimized_double_C(UL n,UL k);//nchoosek function double version
int nchoosek(int n,int k);//nchoosek function
int Factorial(int num);

	//general functions
    void help();
    void optionInfo(int argc, int* argv[],char *inputFile);

    // I/O
    int readInputFile(char *inputFile, char *priorFile);
    int writeResult(char *outputFile);
 
	//DBN related
	float getScore(int **net);		   //get the MIT score of a DBN
    int getNPa(int* &Pa,int **net, int a);	   //get the number of parents of X_i	
	
	float conditional_MI_DBN(int **data,int **datb,int i, int j,int nPa,int *Pa, int n_state);
	
	void Contingency(int Mem1,int Mem2,int n_state);
	float Mutu_Info(int **T, int n_state)		;
	void ClearT(int n_state);			//clear the share contingency table
	
	int compare_Parent_config(int **data,int nPa,int* Pa,int a,int posi, int posj); 
	
	float myEntropy(int x);			//calculate entropy of a 1-unit-shifted in time vector
	int findPstar(float *g_score,float g);
	void updateBestPa(int * &best_Pa,int * Pa,int p);
	int getPosition(int* powo,int p, int * Pa);
	int findLexicalIndex(int n, int p, int * Pa); //Find lexical index of a parent set, variable IDs are from 0->dim-1

	time_t time_=time(NULL);
    void tic(){time_=time(NULL);}
    float toc(){return difftime(time(NULL),time_);}
    void cleanUp();

//////////////////////Global variables//////////////////////////////////
//Data & Data structure
    int            N = 0;				//number of data items
	int			   Ne =0;				//number of effective samples
    int            dim=0;				//number of variables
    int**          data=NULL;			
	int**          datb=NULL;			//1-unit shifted data

	int			   n_state=0;			//number of states, densely mapped [1:n_state], i.e Matlab style
	float*		   chi=NULL;		    //chi values
	int**	       net=NULL;			//the DBN
	int**		prior=NULL;          //MPchange the prior
	int *parset = NULL;
	int**		   T=NULL;				//the shared contingency table
	float*		   g_score=NULL;		//the g-score
// I/O & timing
    ofstream       outLog;              //log report
    ofstream       timeLog;
    time_t		   Start_t=time(NULL);
    time_t         End_t  =time(NULL);
    float          setupTime=0;

//Parameters
    float          alpha=0.999;         //Description: significance level for MI test
    int			   maxFanIn=-1;
	int			   maxTime=-1;
	int			   allowSelfLink=1;     //allow self loop or not?

//debug variables
	int			   readChi=1;			//read in the chi value 
	int			   readNet=0;			//read in the initial net
//MPchange

//-------------------------------------------Main program start----------------------
    int main(int argc, char* argv[]) {

    outLog.open("GlobalMIT_version_0.txt",ios::app);

    outLog<<endl<<"-----------Global RSC release version 1--------------"<<endl;
    cout  <<endl<<"-----------Global RSC release version 1--------------"<<endl;

    if (argc<2){
        help();
        return 0;
    }

    for(int i = 0; i < argc; i++){
        cout   << "argv[" << i << "] = " << argv[i] << endl;
        outLog << "argv[" << i << "] = " << argv[i] << endl;
    }

    
	//char* inputFile="./myDBNinput.txt";
    //char* outputFile="./myDBNoutput.txt";
	//alpha=0.999;

	alpha=atof(argv[1]);
	allowSelfLink=atoi(argv[2]);
    char*  inputFile=argv[3];
    char*  outputFile=argv[4];
    char*  priorFile=argv[5];  //MPchange adding another argument, the prior

    //parsing parameters
    //optionInfo(argc,argv,inputFile);       


    //reading input
    if (readInputFile(inputFile,priorFile)!=0){exit(-1);}
    
	//pre-calculate the g-score
	g_score=new float[16]; //MPchange g_score only goes up to max parents
	g_score[0]=0;	   //g_score[0]->0 parent, g_score[1]->1 parent...
	g_score[1]=chi[1]; //chi[0]-> 0 Parent, chi[1]-> 1 parent...
	for(int i=2;i<=16;i++){   //MPchange dim=6
		g_score[i]=g_score[i-1]+chi[i];		
	}

	//the shared contingency table
	T=new int* [n_state];
	for (int i=0;i<n_state;i++){
		T[i]=new int[n_state];	
		for (int j=0;j<n_state;j++){
			T[i][j]=0;
		}
	}

	//initialize the DBN
	net=new int* [dim];
	for (int i=0;i<dim;i++){
		net[i]=new int[dim];
		for(int j=0;j<dim;j++){
			net[i][j]=0;
		}
	}

//////////////////Main program//////////////////////////////////////////////////
	
	float *best_score_arr=new float [dim];//%best score for each individual node
	float *HX_arr=new float[dim];		//entropy of each (shifted) variable

	cout<<"GlobalMIT Started... "<<endl;
	
	Ne=N; //number of effective samples

	//main loop
	//for(int i=0;i<dim;i++){//loop through the variables

	for(int i=0;i<dim;i++){

//MPchange find all potential parents

		int numpar = 0;
		for(int j=0;j<dim;j++){
			if ((allowSelfLink!=0) || (i!=j)){
			numpar = numpar+prior[j][i];}
		}
		cout<<"Node "<<i <<" Number of parents= "<< numpar<<endl;

		int *parset = new int [numpar];

		if (numpar>0){
		
		int counter = 0;
		for(int j=0;j<dim;j++){
			if ((allowSelfLink!=0) || (i!=j)){
			if (prior[j][i]>0){
				parset[counter]=j;
				counter+=1;}}
		}
		} else { int parset= -1;}


		cout<<"Firstparent "<<parset[0]<<endl;
		HX_arr[i]=myEntropy(i);	
		//investigate all set from 1->P*-1 elements
		//int Pstar=findPstar(g_score,2*Ne*HX_arr[i]);	
		int Pstar = numpar+1;
		cout<<"Node "<<i <<" Pstar= "<< Pstar<<endl;

		int* powo= new int[Pstar-1];
		for(int j=0;j<Pstar-1;j++) {powo[j]= pow(float(dim),j);}//base, for finding the corresponding index in the 1-d array
		
		int *best_Pa=NULL;//empty set
		float best_s_MIT=2*Ne*HX_arr[i]; //score of the empty network
		int	 best_nPa=0;

		//co-ordinate coding: [Pai] -> p*-1 elements
		float * score=new float[dim]; //1-d array to store the scores of all parent combination
		float * new_score=NULL; //1-d array to store the scores of all parent combination
		
		for(int p=1;p<Pstar;p++) { //loop through parent sets of increasing size

		
			int* ca=new int[numpar];//set of all parents
			int* cb=new int[p];  //subset of parent of cardinality p
			int* cbpar = new int[p];  //MPchange parents corresponding to cb indices
			
			for(int j=0;j<numpar;j++)     {ca[j]=j;} //MPchange from dim to numpar 
			for(int j=0;j<p;j++)	   {cb[j]=j;}  //

			if(p>1) {new_score=new float[nchoosek(numpar,p)];}  //MPchange from dim to numpar

			int combi_count=0; //combination count
			do  //generate all parent combination and score
			{
				for(int j=0;j<p;j++){
					cbpar[j]=parset[cb[j]];}
				
				//for (int k=0;k<p;k++){
				//	cout<<cb[k]<<endl;}

				if(allowSelfLink==0){//check self-link
					int selfLink=0;
					for(int j=0;j<p;j++){
						if(parset[cb[j]]==i){
							selfLink=1;
							break;
						}
					}
					if(selfLink) continue;
				}
				
				//int pos=findLexicalIndex(dim,p,cb);   MPchange
				int pos = findLexicalIndex(numpar,p,cb);
				//cout<<"combi_count="<<combi_count<<"Position= "<<pos<<endl;
				combi_count++;
				//score this set of parents
				if (p==1){ //only canculate the score for the 1st level
					 float CMI=conditional_MI_DBN(data,datb,cbpar[0], i, 0,cbpar,n_state);
					 float d_MIT=2*Ne*(HX_arr[i]-CMI);
					 float s_MIT=g_score[p]+d_MIT;
					 if(best_s_MIT>s_MIT){
						 best_s_MIT=s_MIT;
						 updateBestPa(best_Pa,cbpar,p);
						 best_nPa=p;
					 }
					 int pos=cb[0];
					 score[pos]=CMI; //store the score
				}else{
					 float score_i=0;
					 //get from cache
					 //int pos=getPosition(powo,p-1,cb);
					 int pos=findLexicalIndex(numpar,p-1,cb); //MPchange from dim to numpar

					 score_i=score[pos];

					 
					 //calculate the last score and store
					 float CMI=conditional_MI_DBN(data,datb,cbpar[p-1],i, p-1 ,cbpar,n_state);
					 score_i+=CMI;	
					 
					 //pos=getPosition(powo,p,cb);
					 pos=findLexicalIndex(numpar,p,cb);

					 new_score[pos]=score_i; //store the last calculated score
					 //cout<<"hello"<<endl;
					 float d_MIT=2*Ne*(HX_arr[i]-score_i);
					 float s_MIT=g_score[p]+d_MIT;
					 if(best_s_MIT>s_MIT){
						 best_s_MIT=s_MIT;
						 updateBestPa(best_Pa,cbpar,p);
						 best_nPa=p;
					 }

				}
			
			//}while(next_combination(ca,ca+dim,cb,cb+p)); //MPchange take dim to numpar
			}while(next_combination(ca,ca+numpar,cb,cb+p));
		  
			delete[] ca;
			delete[] cb;
			if(p>1) {delete[] score; score=new_score;}
		}// of p loop


		cout<<"Complete scoring all parent sets of node " << i <<"!"<<endl;
		cout<<"Best score: "<<best_s_MIT<< "best Pa=["; for(int k=0;k<best_nPa;k++) {cout<<best_Pa[k]<<" ";}cout<<"]"<<endl;
		best_score_arr[i]=best_s_MIT;
		for(int k=0;k<best_nPa;k++) {net[best_Pa[k]][i]=1;}
		delete[] powo;
		delete[] score;
	}// of i loop


	float best_score=0;
	for(int i=0;i<dim;i++) best_score+=2*Ne*HX_arr[i]-best_score_arr[i];
	cout<<"Best score =" << best_score<<endl;

	float score=getScore(net);
	cout<<"Best score (check)= "<<score<<endl;

    //writeResult(outputFile);
	ofstream  outFile;
	outFile.open(outputFile);
	cout.precision(20);
	outFile<<fixed <<score<<endl;
	for(int i=0;i<dim;i++){
		for(int j=0;j<dim;j++){
			outFile<<net[i][j]<<" ";
		}
		outFile<<endl;
	}
	outFile.close();

    //cleaning up
    cleanUp();
	delete[] best_score_arr;
	delete[] HX_arr;

    return 0;

}
///////////////////////////////////END of MAIN PROGRAM////////////////////////////////////////
int getPosition(int* powo,int p, int * Pa){ //get position in the 1-d array
	int position=0;
	for(int i=0;i<p;i++){
		position+=powo[i]*Pa[i];
	}
	return position;
}

void updateBestPa(int * &best_Pa,int * Pa,int p){ //update the best Parent set
	if (best_Pa!=NULL) delete[] best_Pa;
	best_Pa=new int[p];
	for(int i=0;i<p;i++){
		best_Pa[i]=Pa[i];
	}
}

int findPstar(float *g_score,float g){ //search for the max-fan-in
	int p=0;
	while (g_score[p]<g && p<dim){
		p++;
	}
	return p;
}

//======================Global RSC functions==========================

void optionInfo(int argc, char* argv[],char *inputFile){
    for(int i=3;i<argc;++i){
        if     (!strcmp(argv[i],"-maxFanIn"     )) maxFanIn    = atoi(argv[++i]);
        else if(!strcmp(argv[i],"-maxTime "     )) maxTime     = atoi(argv[++i]);
        else {cout << "Argument not recognized! " << argv[i] << endl; exit(0);}
    }
    cout  <<"---------------"<<endl;outLog<<"---------------"<<endl;
    cout<<"Options:"<<endl;outLog<<"Options:"<<endl;
    if (maxFanIn==1){
        cout  <<"Max-Fan-In= "<<maxFanIn<<endl;
        outLog<<"Max-Fan-In= "<<maxFanIn<<endl;
    }
    cout  <<"Maximum time: "<<maxTime<<endl;
    outLog<<"Maximum time: "<<maxTime<<endl;

    cout  <<"---------------"<<endl;outLog<<"---------------"<<endl;
}

int readInputFile(char *inputFile, char *priorFile){
//MPchange adding prior as a variable
	//read the data file
	cout  <<"Reading numeric input file: "<< inputFile<<endl;
	outLog<<"Reading numeric input file: "<< inputFile<<endl;
   
	tic();   
	FILE * inFile;
	inFile = fopen (inputFile,"r");

	if (inFile==NULL) {
		cout << "Unable to open Data file";outLog << "Unable to open Data file";
		outLog.close();
		exit(-1); // terminate with error
	}
	fscanf (inFile, "%d %d", &N,&dim);
	cout<< "N= "<< N << " Dim= " <<dim<<endl;
   data=new int* [N];
   for(int i=0;i<N;i++){
	   data[i]=new int[dim];
		for (int j=0;j<dim;j++){
			fscanf (inFile, "%d",&data[i][j]);
			if (data[i][j]>n_state) {n_state=data[i][j];}
			data[i][j]--;  //remapping the data to 1-> n_state-1;
			//cout<<data[i][j]<<" ";
		}
		//cout<<endl;
   }

   datb=new int* [N];
   for(int i=0;i<N;i++){
	   datb[i]=new int[dim];
		for (int j=0;j<dim;j++){
			fscanf (inFile, "%d",&datb[i][j]);
			if (datb[i][j]>n_state) {n_state=datb[i][j];}
			datb[i][j]--;  //remapping the data to 1-> n_state-1;
			//cout<<datb[i][j]<<" ";
		}
		//cout<<endl;
   }


   cout<<"N_state= "<<n_state<<endl;

   fclose (inFile);

   if (readChi){

	   inFile = fopen ("myChiValue_nstate_3_alpha_0.999.txt","r");

	   chi=new float[16]; //MPchange maximum chi is not the total number of features!

	   chi[0]=0;  //chi[0]-> 0 Parent, chi[1]-> 1 parent...
	   for (int i=1;i<=16;i++){

			fscanf (inFile, "%f",&chi[i]);

			cout<<"Chi["<<i<<"]="<<chi[i]<<endl;
	   }
	   fclose (inFile);
	}


//MPchange reading in the prior
	inFile = fopen (priorFile,"r");
	if (inFile!=NULL){
	prior=new int* [dim];
	for (int i=0;i<dim;i++){
		prior[i]=new int[dim];
		for (int j=0;j<dim;j++){
		fscanf (inFile, "%d",&prior[i][j]);
		}
	   }
	   fclose (inFile);
	}
	cout<<"Net00"<<prior[0][0]<<endl;

   if (readNet){
	   inFile = fopen ("./myNet.txt","r");
	   net=new int* [dim];
	    for(int i=0;i<dim;i++){
		   net[i]=new int[dim];
				for (int j=0;j<dim;j++){
					fscanf (inFile, "%d",&net[i][j]);	
					//cout<<net[i][j]<<" ";
				}
			//cout<<endl;
		}
	}

   cout  <<"DONE in " << toc() << " seconds."<<endl;
   outLog<<"DONE in " << toc() << " seconds."<<endl;

return (0);
}

void cleanUp(){
    if (data!=NULL){
        for (int i=0;i<N;i++){delete data[i];}
        delete[] data;
    }

	if (datb!=NULL){
        for (int i=0;i<N;i++){delete datb[i];}
        delete[] datb;
    }

	if (T!=NULL){
        for (int i=0;i<n_state;i++){delete T[i];}
        delete[] T;
    }

	if (net!=NULL){
        for (int i=0;i<dim;i++){delete net[i];}
        delete[] net;
    }

	if (chi!=NULL)             delete[] chi;
}

void help(){
    using namespace std;
    cout << "Usage: globalMIT.exe <alpha> <allowSelfLink> <inputFile> <outputFile> [options]"<<endl;
	cout << "<alpha>				: Significance level for the Mutual Information test\n";
	cout << "<allowSelfLink>	    : 1 for allowing self-link, 0 otherwise\n";
    cout << "<inputFile>            : Ascii file containing the data matrix. The first rows contain the number of rows and collumns\n";
    cout << "<outputFile>           : Ascii file containing the graph\n";
    cout << "-----------------------------------------------------------------------------------------------------\n";
    //cout << "options :\n";
    //cout << "-maxFanIn              : set the max-fan-in parameter\n";
    //cout << "-maxTime               : set max run time\n";
    //cout << "-----------------------------------------------------------------------------------------------------\n";
}

float getScore(int **net){
	float score=0;
	for(int i=0;i<dim;i++){
		int* Pa=NULL;
		int nPa=getNPa(Pa,net,i);
		
		float score_i=0;
		//cout<<"Node "<<i << " nPa= "<<nPa<<endl;  //MPchange commented out
		for(int j=0;j<nPa;j++){
			score_i+=2*Ne*conditional_MI_DBN(data, datb, Pa[j],i,j,Pa,n_state)-chi[j+1];
		}
		//cout<<"Node "<< i <<" score " << score_i<<endl;  //MPchange commented out
		score+=score_i;
	}	
	cout<<"Total score: "<<score <<endl;
	return score;
}

int getNPa(int* &Pa,int **net, int a){  //get the number of parents of X_a and put the set in Pa
	int nPa=0;
	for(int j=0;j<dim;j++){
		nPa=nPa+net[j][a];
	}

	if(Pa!=NULL)  {delete [] Pa;}

	Pa=new int[nPa];
	int pos=0;
	for(int j=0;j<dim;j++){
		if (net[j][a]){ //add this parent
			Pa[pos]=j;
			pos++;
		}
	}
	return nPa;
}


//conditional MI between node a-> node b given other parent Pa
float conditional_MI_DBN(int **data,int **datb,int a, int b,int nPa, int* Pa, int n_state){
float MI=0;

if (nPa==0){ //no parent
	 Contingency(a,b,n_state);
	 return Mutu_Info(T, n_state);
}
else {	//with some parents?
	int  * scanned=new int[N];

	for(int i=0;i<N;i++){scanned[i]=0;}

	for(int i=0;i<N-1;i++){ //scan all rows of data
		if(scanned[i]==0){  //a new  combination of Pa found		
			scanned[i]=1;
			float count=1;
			ClearT(n_state);
			T[data[i][a]][datb[i][b]]++;

			for(int j=i+1;j<N;j++){
				if(scanned[j]==0 && compare_Parent_config(data,nPa,Pa,b,i,j)){
					scanned[j]=1;	 				
					T[data[j][a]][datb[j][b]]++;
					count++;
				}
			}
			MI+=(count/Ne)*Mutu_Info(T,n_state);
		}
	}
	delete[] scanned;	
}

return MI;
}


//compare a parent set configuration of node a at two position in the data
int compare_Parent_config(int **data,int nPa,int* Pa,int a,int posi, int posj){
	int	isSame=1;
	for (int i=0;i<nPa;i++){ //scan through the list of parents
		if(data[posi][Pa[i]]!=data[posj][Pa[i]]){//check this parent value at posi & posj
			return 0;
		}
	}
return isSame;
}


//calculate the unconditional contingency table between node a=> node b (shifted, store to the global var T
void Contingency(int a,int b,int n_state){
	//clear up T
	ClearT(n_state);

	//build table
	for(int i =0;i<N;i++){
		T[data[i][a]][datb[i][b]]++;  //note: b is one unit shifted in time
	}
}

void ClearT(int n_state){
	for(int i=0;i<n_state;i++){
		for(int j=0;j<n_state;j++){
			T[i][j]=0;
		}
	}
}

float Mutu_Info(int **T, int n_state){  //get the mutual information from a contingency table
	float MI=0;
	int *a = new int[n_state];
	int *b = new int[n_state];
	int N=0;

	for(int i=0;i<n_state;i++){ //row sum
		a[i]=0;
		for(int j=0;j<n_state;j++)
		{a[i]+=T[i][j];}
	}

	for(int i=0;i<n_state;i++){ //col sum
		b[i]=0;
		for(int j=0;j<n_state;j++)
		{b[i]+=T[j][i];}
	}

	for(int i=0;i<n_state;i++) {N+=a[i];}
	
	for(int i=0;i<n_state;i++){
		for(int j=0;j<n_state;j++){
			if(T[i][j]>0){
				MI+= T[i][j]*log(float(T[i][j])*N/a[i]/b[j]);
			}
		}
	}
	delete []a;
	delete []b;

	return MI/N;
}

float myEntropy(int x){
	float *H =new float[n_state];
	for(int i=0;i<n_state;i++) {H[i]=0;}
	//entropy of a 1-unit-shifted in time vector
	for(int i=0;i<N;i++){// i run from 1
		H[datb[i][x]]++;	
	}
	float e=0;
	for(int i=0;i<n_state;i++) {
		H[i]/=Ne;
		if (H[i]!=0) {e-=H[i]*log(H[i]);}		
	}
	return e;
}

int findLexicalIndex(int n, int p, int * Pa){
if(p==0) return -1;
if(p==1) return Pa[0];

int pos=0;
int last_pos=0;

for(int i=0;i<p-1;i++){
	if(i==0) {last_pos=0;}
	else{
		last_pos=Pa[i-1]+1;
	}
	for(int j=last_pos;j<Pa[i];j++){
		pos=pos+nchoosek(n-(j+1),p-(i+1));
	}

}//for i

pos=pos+Pa[p-1]-Pa[p-2]-1;
return pos;

}


int nchoosek(int n,int k){
    int i,temp = 1;
    if(k > (n/2))
         k = n-k;
    for(i = n; i >= (n-k+1); i--){
       temp = temp * i;
    }
    return (temp/Factorial(k));
}


double optimized_double_C(UL n,UL k){
    double answer = 1.0;
    UL i;
    if(k > (n/2))
        k = n-k;
    for(i = 0; i < k; i++){
         answer = answer * ((double)(n-i)/(double)(i+1));
      }
    return answer;
}

int Factorial(int num) {
   int res = 1;
    while(num > 0){
        res = res * num;
        num = num - 1;
    }
   return res;
} 


