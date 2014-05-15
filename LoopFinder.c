//---------------------------------------------------------------------------------
// Brute force search for closed motif detection in undirected simple graphs.
// Matheus Viana, 19.04.2014 - vianamp@gmail.com
//---------------------------------------------------------------------------------

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <list>

int   N;// Number of nodes

int lmax;// Length of the paths to be analysed.

bool *A;// Adjacenty Matrix
        // Element Aij should be accessed as A[i+j*N]. A is assummed to be symmetric
        // and "multiple" or "auto" connections are not allowed. Aij = 1 means that
        // there is a connection between nodes i and j. Aij = 0 means no connection.

bool *U;// Vector of size N that flags nodes that are being used in the current path.
        // U[i] = 1 means that node i is currently being used and shouldn't be used
        // any more for the current path.

int  *L;// Vector of size N that counts the number of loops that each node belongs
        // to.

std::list<bool*> UsageList;// List that will store all motifis found.
std::list<bool*>::iterator itbool;// Iterator for the list.

std::list<int*> Motifs;// List that will store all motifis found.
std::list<int*>::iterator itint;// Iterator for the list.

FILE *f;// Pointer for the network gnet file.

// Binary Jaccard distance between two vectors. The input should be U-type vectors
// that represent the vertices used in two different paths. If the paths are the
// same, function should return 1, or 0 otherwise.

bool _are_similar(bool *A, bool *B) {
  int s(0);
  for (int i=N;i--;) {
    s += (A[i]==B[i]) ? 1 : 0;
  }
  return (s==N) ? 1 : 0;
}

void _add_new_motif(int *V) {
  int *newM = new int[lmax];
  bool *newU = new bool[N];
  for ( int i=N; i--; ) newU[i]=U[i];
  for ( int i=lmax; i--; ) newM[i]=V[i];
  UsageList.insert(UsageList.begin(),newU);  
  Motifs.insert(Motifs.begin(),newM);  
}

// Check whether the new motif just found is already present in the list "UsageList"
// or not. If not, add the new motif.

void _check_new_motif(int *V) {
  if (UsageList.size()>0) {
    bool *addedU;
    for ( itbool=UsageList.begin(); itbool!=UsageList.end(); ++itbool) {
      addedU = *itbool;
      if (_are_similar(U,addedU)) return;
    }
  }
  _add_new_motif(V);
}

void _print_motifs() {
  int i, *W;
  printf("@Number of motifs of length %d found: %d.\n",lmax,(int)Motifs.size());
  printf("@List of motifs:\n");
  for (itint=Motifs.begin(); itint!=Motifs.end(); ++itint) {
    for (i=0;i<lmax;i++) {
      W = *itint;
      printf("%d ",W[i]);
      L[W[i]]++;
    }
    printf("\n");
  }
  printf("@Participation of each node:\n");
  for (i=0;i<N;i++) {
    printf("%d\n",L[i]);
  }
}

// Recursive function to take one step further in the current path. This function is
// actually an implementation of a number "length_max" of nested loops of type
// for (int i = N; i--;). The idea is to generate all possible sequences of "length"
// interger numbers without repetition and to check if the given sequence corresponds
// to a closed path in the network.

void _advance(int length, int *V) {
  int i, j = V[length];
  if ( length ) {
    for ( i = N; i--; ) {
      if ( !U[i] && A[i+j*N] ) {
        V[length-1] = i;
        U[i] = true;
        _advance(length-1, V);
        U[i] = false;
      }
    }
  } else {
      if ( V[0] == V[lmax] ) {
        _check_new_motif(V);
      }
      U[V[0]] = false;
  }
}

// Reading the gnet file and initializing the both adjacency matrix and vector U.
// Edges weights are not being used in the current version (third column of gnet
// file), but in future we might be interested in using it in order to calculate
// motif properties, such as total real length.

void _read_network(const char *filename){
  float w;
  int i, j;
  char path[32];
  sprintf(path,"%s.gnet",filename);
  f = fopen(path,"r");
  fscanf(f,"%d",&N);
  L = new int[N];
  U = new bool[N];
  A = new bool[N*N];
  for ( i=N; i--; ) {
    L[i] = 0;
    U[i] = false;
    for ( j=N; j--; ) A[i+j*N] = false;
  }
  while ( fscanf(f,"%d %d %f", &i, &j, &w) != EOF ) {
    if (i!=j) A[i+j*N] = A[j+i*N] = true;
  }  
  fclose(f);
}

// Export network as gml file

void _export_gml(const char *filename) {

  _print_motifs();

  int i, j;
  char path[32];
  sprintf(path,"%s.xgmml",filename);
  f = fopen(path,"w");

  fprintf(f,"<?xml version=\"1.0\"?>\n");
  fprintf(f,"<graph label=\"%s\" xmlns:dc=\"http://purl.org/dc/elements/1.1/\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" xmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#\" xmlns:cy=\"http://www.cytoscape.org\" xmlns=\"http://www.cs.rpi.edu/XGMML\" directed=\"%d\">",filename,0);

  for ( i = 0; i < N; i++ ){
    fprintf(f,"<node id=\"%d\" label=\"%d\" name=\"%d\">\n",i,i,i);
    fprintf(f,"\t<att type=\"real\" name=\"cycles\" value=\"%f\"/>\n",(float)L[i]);
    fprintf(f,"\t<graphics fill=\"GREEN\" outline=\"BLACK\" h=\"%f\" w=\"%f\" x=\"0.0\" y=\"0.0\" type=\"ELLIPSE\"/>\n",1+sqrt(10*(float)L[i]),1+sqrt(10*(float)L[i]));
    fprintf(f,"</node>\n");
  }

  int ecount = 0;
  for ( i = 0; i < N; i++ ) {
    for ( j = i+1; j < N; j++ ) {
      if (A[i+j*N]) {
        fprintf(f,"<edge source=\"%d\" target=\"%d\" label=\"%d\">\n",i,j,ecount);
        fprintf(f,"</edge>\n");
        ecount++;
      }
    }
  }

  fprintf(f,"</graph>\n");
  fclose(f);
}

// Main function: calculate all paths of specified length in the specified network.

int main(int argc, char *argv[]) {

  _read_network(argv[1]);
  lmax = atoi(argv[2]);
  int V[lmax+1];

  for (int i=N; i--; ) {
    V[lmax] = i;
    _advance(lmax,V);
  }
 
  _export_gml(argv[1]);

  return 0;
}