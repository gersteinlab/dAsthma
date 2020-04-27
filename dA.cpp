/*
refer:
http://deeplearning.net/tutorial/
https://github.com/yusugomori/DeepLearning

*/


#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string.h>
#include <string>
#include <fstream>
#include <sstream>
#include <math.h>
#include <ctime>
#include <vector>
#include "dA.h"
using namespace std;


double uniform(double min, double max) {
  return rand() / (RAND_MAX + 1.0) * (max - min) + min;
}

int binomial(int n, double p) {
  if(p < 0 || p > 1) return 0;
  
  int c = 0;
  double r;
  
  for(int i=0; i<n; i++) {
    r = rand() / (RAND_MAX + 1.0);
    if (r < p) c++;
  }

  return c;
}

double sigmoid(double x) {
  return 1.0 / (1.0 + exp(-x));
}


dA::dA(int size, int n_v, int n_h, double **w, double *hb, double *vb) {
  N = size;
  n_visible = n_v;
  n_hidden = n_h;

  if(w == NULL) {
    W = new double*[n_hidden];
    for(int i=0; i<n_hidden; i++) W[i] = new double[n_visible];
    double a = 1.0 / n_visible;

    for(int i=0; i<n_hidden; i++) {
      for(int j=0; j<n_visible; j++) {
        W[i][j] = uniform(-a, a);
      }
    }
  } else {
    W = w;
  }

  if(hb == NULL) {
    hbias = new double[n_hidden];
    for(int i=0; i<n_hidden; i++) hbias[i] = 0;
  } else {
    hbias = hb;
  }

  if(vb == NULL) {
    vbias = new double[n_visible];
    for(int i=0; i<n_visible; i++) vbias[i] = 0;
  } else {
    vbias = vb;
  }
}

dA::~dA() {
    /* */
    for(int i=0; i<n_hidden; i++) delete[] W[i];
    delete[] W;
    delete[] hbias;
    delete[] vbias;
   
}

void dA::get_corrupted_input(int *x, int *tilde_x, double p) {
  for(int i=0; i<n_visible; i++) {
    if(x[i] == 0) {
      tilde_x[i] = 0;
    } else {
      tilde_x[i] = binomial(1, p);
    }
  }
}
void dA::get_corrupted_input(double *x, double *tilde_x, double p) {
    //reset if the value 
    for(int i=0; i<n_visible; i++) {
        double tmp= rand()/(RAND_MAX+1.0);
        if (tmp > p){
            tilde_x[i]=uniform(0.0, 1.0);
        }else{
            tilde_x[i]=x[i];
        }
        
    }
}


// Encode
void dA::get_hidden_values(int *x, double *y) {
  for(int i=0; i<n_hidden; i++) {
    y[i] = 0;
    for(int j=0; j<n_visible; j++) {
      y[i] += W[i][j] * x[j];
    }
    y[i] += hbias[i];
    y[i] = sigmoid(y[i]);
  }
}

// Encode
void dA::get_hidden_values(double *x, double *y) {
  for(int i=0; i<n_hidden; i++) {
    y[i] = 0;
    for(int j=0; j<n_visible; j++) {
      y[i] += W[i][j] * x[j];
    }
    y[i] += hbias[i];
    y[i] = sigmoid(y[i]);
  }
}



// Decode
void dA::get_reconstructed_input(double *y, double *z) {
  for(int i=0; i<n_visible; i++) {
    z[i] = 0;
    for(int j=0; j<n_hidden; j++) {
      z[i] += W[j][i] * y[j];
    }
    z[i] += vbias[i];
    z[i] = sigmoid(z[i]);
  }
}



void dA::train(double *x, double lr, double corruption_level) {
  double *tilde_x = new double[n_visible];

  
  //sgd
  double *y = new double[n_hidden];
  double *z = new double[n_visible];
  // vbias
  // hbias 
  double *L_vbias = new double[n_visible];
  double *L_hbias = new double[n_hidden];

  double p = 1 - corruption_level;

  
  get_corrupted_input(x, tilde_x, p);
  get_hidden_values(tilde_x, y);
  get_reconstructed_input(y, z);

  //get neg-likelihood
  

  
    for(int i=0; i< n_visible; i++){
        L_vbias[i] = x[i] - z[i];
        neglik += x[i]*log(z[i])+(1-x[i])*log(1-z[i]);
        vbias[i] += lr * L_vbias[i] / N;

    }

    for(int i=0; i<n_hidden; i++) {
      L_hbias[i] = 0;
      for(int j=0; j<n_visible; j++) {
          L_hbias[i] += W[i][j] * L_vbias[j];
      }
      L_hbias[i] *= y[i] * (1 - y[i]);

      hbias[i] += lr * L_hbias[i] / N;
    }

  // W
    for(int i=0; i<n_hidden; i++) {
        for(int j=0; j<n_visible; j++) {
            W[i][j] += lr * (L_hbias[i] * tilde_x[j] + L_vbias[j] * y[i]) / N;
        }
    }
  
  delete[] L_hbias;
  delete[] L_vbias;
  delete[] z;
  delete[] y;
  delete[] tilde_x;
}


void dA::reconstruct(int *x, double *z) {
  double *y = new double[n_hidden];

  get_hidden_values(x, y);
  get_reconstructed_input(y, z);

  delete[] y;
}

void dA::reconstruct(double *x, double *z) {
  double *y = new double[n_hidden];

  get_hidden_values(x, y);
  get_reconstructed_input(y, z);

  delete[] y;
}

void usage(char* prog){
    std::cout<< "Usage:" << prog << " [Options]\n"
             <<" -c tranning cycle[30]\n"
             <<" -0 number of hidden[30]\n"
             <<" -i Input file contains only gene exprs convertedinto[0,1]\n"
             <<" -r learnning rate[0.1]\n"
             <<" -e Error rate, corruption_level [0.2]\n"
             <<" -o output file[tmp,tmp.hidden.txt,tmp.log.txt,tmp.weight.txt]\n"
             <<" -w ncol, width of input file\n"
             <<" -h nrow, height of input file\n";
    

    exit(1);
}

int main(int argc, char* argv[]) {
    //test_dA();
    char* input=NULL;
    //string input=string("/home2/sl2373/data4dA.txt");
    unsigned int nrow=112;
    unsigned int ncol=22148;
    unsigned int training_epochs=30;
    double lr=0.1;
    double err=0.1;
    double p=0.3; //the number of sample to be left-out
    string line,tmp;
    char* output="tmp";
    //nsample,

    //double train_X[train_row][ncol]={};
    //double test_X[test_row][ncol]={};

    int n_visible=ncol, n_hidden=20;
    int c;
    /* */
    while((c = getopt(argc, argv,"0:i:r:e:c:o:h:w:u"))!=-1)
    {
        switch(c)
        {
            case 'c': training_epochs=atoi(optarg); break;
            case '0': n_hidden=atoi(optarg);break;
            case 'i': input=optarg;break;
            case 'r': lr=atof(optarg);break;  //learning rate
            case 'e': err=atof(optarg);break;
            case 'o': output=optarg;break;
            case 'h': nrow=atoi(optarg); break; //nrow
            case 'w': ncol=atoi(optarg); break; //ncol
            case 'u': usage(argv[0]); exit(1);
            default:  usage(argv[0]); exit(1);
        }
    }

    if (nrow < 20 | lr >1 | lr < 0 | err < 0 | err > 1 | p > 1 | p < 0 | ncol < 0 | nrow < 0 | (!input)){
        usage(argv[0]);
        exit(1);
        
    }
    
    int train_row= (int)ceil(nrow * (1-p));
    int test_row = nrow -train_row;
    n_visible=ncol;
    //int train_row = nrow;
    double ** train_X=NULL;
    double ** test_X=NULL;

    train_X = (double**) malloc(sizeof(double*) * train_row);
    for(int i=0; i < train_row; i++){
        train_X[i] = (double*) malloc(sizeof(double)* ncol);
        //memset(train_X[i], 0, sizeof(double)*ncol);
    }
    //memset(train_X, 0, sizeof(double*)*train_row);
    //generate testing set    
    test_X = (double**) malloc(sizeof(double*) * test_row);
    for(int i=0; i < test_row; i++){
       test_X[i] = (double*) malloc(sizeof(double)* ncol);
    }
    //memset(test_X, 0, sizeof(double)*test_row*ncol);
        
    
    //read data to train and testing set
    ifstream ifs;
    ifs.open(input);
    int train_idx=0, test_idx=0;
    cout << "Randomize" << endl;
    srand(time(NULL));
    while(getline(ifs,line)){
        double pp= uniform(0.0,1.0);
        bool test_flag=false;
        
        if (pp<=(p*1.1)){
            test_flag=true;
        }
        istringstream ss(line);
        //vector<double> expr;
        int colidx=0;
        cout << pp <<",";
        if (train_idx >= train_row){
            test_flag=true;
        }
        while(getline(ss,tmp,'\t')){
            // expr.push_back(atof(tmp.c_str()));
            if(test_flag && test_idx<test_row){// 
                test_X[test_idx][colidx] = atof(tmp.c_str());
            }
            else{
                train_X[train_idx][colidx] = atof(tmp.c_str());
                //train_idx++;
            }
            colidx++;
        }
        
        if(test_flag && test_idx < test_row){
            test_idx++;
        }else{
            train_idx++;
        }
    }
    ifs.close();
    cout << endl;
    cout << "Reading data done..., start trainning"<< endl;
    //
    dA da(train_row, n_visible,n_hidden, NULL, NULL, NULL);
    dA da2(test_row, n_visible, n_hidden, NULL, NULL, NULL); //for testing
    
    
    char* logfilename=NULL;
    logfilename=(char*)malloc(sizeof(char)*(strlen(output)+20));
    strcpy(logfilename,output);
    strcpy(logfilename+strlen(output), ".log.txt");
    FILE *logfile=fopen(logfilename,"w");
    fprintf(logfile, ">>>%d\t%d\t%d\n",n_visible, train_row,test_row);

    for(int epoch=0; epoch < training_epochs; epoch++){
        da.neglik=0.0;
        for(int i=0; i < train_row; i++){
            da.train(train_X[i], lr, err);
        }

 
        for (int m =0; m< n_hidden;m++){
            da2.hbias[m] = da.hbias[m];
            for (int n=0; n< n_visible; n++){
                da2.W[m][n]= da.W[m][n];
                if(m==0){
                    da2.vbias[n]=da.vbias[n];
                }
            }
        }


        da2.neglik=0.0;
        for (int i=0; i < test_row; i++){
            double *y = new double[n_hidden];
            double *z = new double[n_visible];
            da2.get_hidden_values(test_X[i], y); //x n by 2k; 2kxn_hidden, n by n-hidden; x n-hidden x 2k = n by 2k
            da2.get_reconstructed_input(y, z);
            for (int j=0; j < n_visible; j++){
                da2.neglik += test_X[i][j] * log(z[j])+ (1-test_X[i][j])*log(1-z[j]);
            }

            delete[] y;
            delete[] z;

            
        }

        


        
        
        fprintf(logfile, ">>>%d\t%f\t%f\n",epoch,-1.0*da.neglik, -da2.neglik);
        
        cout <<"Run("<< epoch<<")"<<", cross_entropy="<< -da.neglik << ", test loss=" << -da2.neglik << endl;

        
        //output likelihood#ifdef DEBUG
        
    }
    fclose(logfile);

    //output value
    //output W and y and Z
    char* hidfile=NULL;
    hidfile=(char*)malloc(sizeof(char)*(strlen(output)+20));
    strcpy(hidfile,output);
    strcpy(hidfile+strlen(output), ".hidden.txt");
    
    FILE * hiddenfs=fopen(hidfile,"a");
    for(int i=0; i < train_row; i++){
        double* y= new double[n_hidden];
        da.get_hidden_values(train_X[i], y);
        //output y;
        fprintf(hiddenfs,"%d",i);
        for(int j=0; j<n_hidden; j++){
            fprintf(hiddenfs, "\t%f", y[j]);
        }
        fprintf(hiddenfs,"\n");
        
    }
    fclose(hiddenfs);

    char* weightfile=NULL;
    weightfile = (char*) malloc(sizeof(char)*(strlen(output)+20));
    strcpy(weightfile, output);
    strcpy(weightfile+strlen(output),".weight.txt");
    
    FILE * wfs=fopen(weightfile,"a");
    for(int i=0; i < n_hidden; i++){
        fprintf(wfs,"%d",i);
        for(int j=0; j<n_visible; j++){
            //
            fprintf(wfs,"\t%f",da.W[i][j]);

        }
        fprintf(wfs,"\n");
    }
    fclose(wfs);
    


    return 0;
}
