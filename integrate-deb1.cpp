#include <gravray.cpp>
using namespace std;

#define VERBOSE 1

int main(int argc,char* argv[])
{
  /*
    Example: ./select.exe 

    Input: 
    * cloud.data
    * RVGaia.data
    Output: 
  */
  ////////////////////////////////////////////////////
  //INITIALIZE CSPICE
  ////////////////////////////////////////////////////
  int Npart=1;
  if(argc>1){
    Npart=atoi(argv[1]);
  }
  VPRINT(stdout,"Calculating with %d test particles\n",Npart);

  ////////////////////////////////////////////////////
  //INITIALIZE CSPICE
  ////////////////////////////////////////////////////
  initSpice();
  
  ////////////////////////////////////////////////////
  //GENERAL VARIABLES
  ////////////////////////////////////////////////////
  double tmp;
  char ctmp[100],line[10000];

  ////////////////////////////////////////////////////
  //UNITS
  ////////////////////////////////////////////////////
  /*
  //GRAVITATIONAL SYSTEM
  UL=8.4e3*PARSEC;
  UM=1e11*MSUN;
  UT=sqrt(UL*UL*UL/(GCONST*UM));
  */

  //*
  //ASTRONOMICAL
  UM=MSUN;
  UL=PARSEC;
  UT=YEAR;
  //*/

  /*
  //KM,KM/S
  UL=1e3;
  UT=1;
  //*/

  //UNITS OF VELOCITY
  UV=UL/UT;

  VPRINT(stdout,"Units:\n\tUM = %.5e kg=%.5e Msun\n\tUL = %.17e m = %.17e pc\n\tUT = %.5e s = %.5e yr\n\tUV = %.5e m/s = %.5e km/s\n",
	 UM,UM/MSUN,UL,UL/PARSEC,UT,UT/YEAR,UV,UV/1e3);

  ////////////////////////////////////////////////////
  //TEST CONVERSION
  ////////////////////////////////////////////////////
  double v[6];
  double xLSR_s[]={-3.16028e+02*PARSEC, 1.84156e+01*PARSEC, -3.58527e+02*PARSEC, 
                   1.33746e+01, -1.72967e+01, -1.54271e+01};
  copyVec(v,xLSR_s,6);
  vscl_c(1/PARSEC,xLSR_s,v);
  VPRINT(stdout,"LSR (star) = %s\n",vec2strn(v,6,"%.5e "));
  double xGC_s[6];
  LSR2GC(xLSR_s,xGC_s);
  copyVec(v,xGC_s,6);
  vscl_c(1/PARSEC,xGC_s,v);
  VPRINT(stdout,"GC (star) = %s (pc,km/s)\n",vec2strn(v,6,"%.5e "));

  ////////////////////////////////////////////////////
  //READING TEST PARTICLES
  ////////////////////////////////////////////////////
  int Ntimes=2;
  int nsys=6*Npart;
  
  FILE *fc=fopen("cloud.data","r");
  double **xLSR,**xGC,xGC0[nsys],ts[Ntimes],tbody;
  double x0[nsys];
  
  xLSR=(double**)malloc(Ntimes*sizeof(double*));
  for(int j=0;j<Ntimes;j++) xLSR[j]=(double*)malloc(nsys*sizeof(double));

  xGC=(double**)malloc(Ntimes*sizeof(double*));
  for(int j=0;j<Ntimes;j++) xGC[j]=(double*)malloc(nsys*sizeof(double));

  VPRINT(stdout,"Size of each state vector: %d\n",nsys);

  for(int i=0;i<Npart;i++){
    VPRINT(stdout,"Reading particle %d...\n",i);

    //STRIDE
    int ip=6*i;

    //READ POSITIONS
    for(int j=0;j<=39;j++) fscanf(fc,"%lf",&tmp);
    for(int j=0;j<6;j++){
      fscanf(fc,"%lf",&xLSR[0][ip+j]);
      if(j<3) xLSR[0][ip+j]*=1e3;//TRANSFORM TO METERS
    }
    for(int j=46;j<=58;j++) fscanf(fc,"%lf",&tmp);

    //TRANSFORM FROM LSR->GC
    LSR2GC(xLSR[0]+ip,xGC0+ip);
    copyVec(v,xLSR[0]+ip,6);
    vscl_c(1/PARSEC,xLSR[0]+ip,v);
    VPRINT(stdout,"LSR = %s\n",vec2strn(v,6,"%.5e "));
    copyVec(v,xGC0+ip,6);
    vscl_c(1/PARSEC,xGC0+ip,v);
    VPRINT(stdout,"GC = %s\n",vec2strn(v,6,"%.5e "));

    //CONVERT TO DYNAMIC UNITS
    for(int j=0;j<6;j++){
      if(j<3) xLSR[0][ip+j]/=UL;
      else xLSR[0][ip+j]/=UV/1e3;
    }
    for(int j=0;j<6;j++){
      if(j<3) xGC0[ip+j]/=UL;
      else xGC0[ip+j]/=UV/1e3;
    }
    //exit(0);
    //copyVec(xGC0+ip,xLSR[0]+ip,6); //SIMPLE
  }
  fclose(fc);

  //STATE VECTOR OF PARTICLE 5
  /*
  int ip=6*5;
  fprintf(stdout,"%s\n",vec2strn(xLSR[0]+ip,6,"%.17e "));
  */

  ////////////////////////////////////////////////////
  //INTEGRATING TEST PARTICELS
  ////////////////////////////////////////////////////
  double dydt[nsys];
  double params[]={nsys,Npart};
  double tini=0.0;
  double h=1e4*YEAR/UT;
  double duration=+6.22197e+06*YEAR/UT;

  VPRINT(stdout,"Integrating from t = %.5e during %.5e UT (h = %.5e)...\n",
	 tini,duration,h);
  EoMGalactic(0,xGC0,dydt,params);
  VPRINT(stdout,"Initial Conditions:\n");
  VPRINT(stdout,"\ty = %s\n",vec2strn(xGC0,nsys,"%.7e "));
  VPRINT(stdout,"\tdydt = %s\n",vec2strn(dydt,nsys,"%.17e "));

  //INTEGRATE
  //*
  integrateEoM(0,xGC0,h,Ntimes,duration,
	       nsys,EoMGalactic,params,
	       ts,xGC);
  //*/
  
  VPRINT(stdout,"SOLUTION:\n");
  for(int i=0;i<Ntimes;i++){
    VPRINT(stdout,"t = %.5e\n",ts[i]);
    VPRINT(stdout,"\ty = %s\n",vec2strn(xGC[i],nsys,"%.5e "));
  }

  //TRANSFORM GC->LSR

  //COMPUTE DIFFERENCE
  double dvec[3],dist;
  vsub_c(xGC[1],xGC[0],dvec);
  dist=vnorm_c(dvec);
  fprintf(stdout,"dvec = %s, D = %.5e\n",vec2str(dvec,"%.5e"),dist);
  exit(0);

  ////////////////////////////////////////////////////
  //INTEGRATE TEST PARTICLES
  ////////////////////////////////////////////////////

  ////////////////////////////////////////////////////
  //READING SELECTED OBJECTS
  ////////////////////////////////////////////////////
  fc=fopen("candidates.csv","r");

  int Nsur=100;
  int nfields;

  char **fields=(char**)malloc(MAXCOLS*sizeof(char*));
  for(int i=0;i<MAXCOLS;i++) fields[i]=(char*)malloc(MAXTEXT*sizeof(char));
  while(fscanf(fc,"%s",line)==1){
    //PARSE FIELDS
    parseLine(line,fields,&nfields);

    //SPATIAL COORDINATES AND LSR VELOCITY
    
    //INTEGRATE 

  }
  
  return 0;
}
