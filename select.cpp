#include <gravray.cpp>
using namespace std;

#define VERBOSE 0

int main(int argc,char* argv[])
{
  /*
    Example: ./program.exe 10 1
    Where:
    10: Number of particles
    1: Use diagonal covariance

    Input: 

    Output: 
    * File
      Rows: 
      Cols:
          0:i
	  1:tdb (terminal)
  */

  ////////////////////////////////////////////////////
  //INITIALIZE CSPICE
  ////////////////////////////////////////////////////
  initSpice();
  
  ////////////////////////////////////////////////////
  //PARAMETERS
  ////////////////////////////////////////////////////
  int qdiagonal=0;
  int Npart=1;
  if(argc>1){
    Npart=atoi(argv[1]);
  }
  if(argc>2){
    qdiagonal=atoi(argv[2]);
  }
  fprintf(stdout,"Running with Npart = %d, qdiagonal = %d\n",Npart,qdiagonal);

  return 0;
}
