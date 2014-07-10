//---------------------------------------------------------------------------------
// Brute force search for closed loop in undirected simple graphs.
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

std::list<bool*> UsageList;// List that will store all loops found.
std::list<bool*>::iterator itbool;// Iterator for the list.

std::list<int*> Loops;// List that will store all loops found.
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

void _add_new_loop(int *V) {
  int *newM = new int[lmax];
  bool *newU = new bool[N];
  for ( int i=N; i--; ) newU[i]=U[i];
  for ( int i=lmax; i--; ) newM[i]=V[i];
  UsageList.insert(UsageList.begin(),newU);  
  Loops.insert(Loops.begin(),newM);  
}

// Check whether the new loop just found is already present in the list "UsageList"
// or not. If not, add the new loop.

void _check_new_loop(int *V) {
  if (UsageList.size()>0) {
    bool *addedU;
    for ( itbool=UsageList.begin(); itbool!=UsageList.end(); ++itbool) {
      addedU = *itbool;
      if (_are_similar(U,addedU)) return;
    }
  }
  _add_new_loop(V);
}

void _print_loops() {
  int i, *W;
  printf("@Number of loops of length %d found: %d.\n",lmax,(int)Loops.size());
  printf("@List of loops:\n");
  for (itint=Loops.begin(); itint!=Loops.end(); ++itint) {
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
        _check_new_loop(V);
      }
      U[V[0]] = false;
  }
}

// Reading the gnet file and initializing the both adjacency matrix and vector U.
// Edges weights are not being used in the current version (third column of gnet
// file), but in future we might be interested in using it in order to calculate
// loop properties, such as total real length.

void _read_network(const char filename[]){
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

  _print_loops();

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

  int i, n;
  char _impath[128];
  sprintf(_impath,"");

  // Collecting input parameters
  for (i = 0; i < argc; i++) {
      if (!strcmp(argv[i],"-path")) {
          sprintf(_impath,"%s//",argv[i+1]);
      }
      if (!strcmp(argv[i],"-n")) {
          n = atoi(argv[i+1]);
      }
  }

  // Generating list of files to run
  char _cmd[256];
  sprintf(_cmd,"ls %s*.gnet | sed -e 's/.gnet//' > %sLoopFinder.files",_impath,_impath);
  system(_cmd);

  // Iterating over many files in the list
  char _gnetfilename[256];
  char _gnetlistpath[128];
  char _xgmlfilename[256];
  sprintf(_gnetlistpath,"%sLoopFinder.files",_impath);
  FILE *f = fopen(_gnetlistpath,"r");
  while (fgets(_gnetfilename,256, f) != NULL) {
      printf(">>%s",_gnetfilename);
      _gnetfilename[strcspn(_gnetfilename, "\n" )] = '\0';

      // Finding loops
      _read_network(_gnetfilename);
      lmax = n;
      int V[lmax+1];

      for (int i=N; i--; ) {
        V[lmax] = i;
        _advance(lmax,V);
      }

      // Saving GML file
    _export_gml(_gnetfilename);

  }
  fclose(f);

  return 1;
}