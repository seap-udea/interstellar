#include <gravray.cpp>
using namespace std;

#define VERBOSE 1

int main(int argc,char* argv[])
{
  /*
    Example: ./integrate.exe 10

    Where: 
    10: Number of particles

    Input: 
    * cloud.data: test particles compatible with A/2017U1 orbit.

    Output:
    * cloud-int.csv
  */

  ////////////////////////////////////////////////////
  //INPUT PARAMETERS
  ////////////////////////////////////////////////////
  int Npart=1;
  if(argc>1){
    Npart=atoi(argv[1]);
  }
  VPRINT(stdout,"Analysing %d test particles against candidate stars\n",Npart);
  VPRINT(stdout,"\n");

  ////////////////////////////////////////////////////
  //INITIALIZE CSPICE
  ////////////////////////////////////////////////////
  initSpice();
  
  ////////////////////////////////////////////////////
  //GENERAL VARIABLES
  ////////////////////////////////////////////////////
  double tmp;
  char ctmp[100],line[MAXLINE];
  int Ntimes,Ntimesp,nsys,nsysp; 
  int ip,n;
  int nfields;
  double params[10];
  double hstep,duration;
  //MATRICES WITH INTEGRATIONS
  double *xIntp0,**xIntp,**xIntc;
  double *xInt0,**xInt;
  double *tsp,*ts;
  //INITIAL CONDITIONS
  double *dxIntdt,*x,*xg,*x0;
  double G;

  ////////////////////////////////////////////////////
  //UNITS
  ////////////////////////////////////////////////////
  UM=MSUN;
  UL=PARSEC;
  UT=YEAR;
  UV=UL/UT;
  G=GCONST/(UL*UL*UL/(UM*UT*UT));
  
  VPRINT(stdout,"Units:\n\tUM = %.5e kg=%.5e Msun\n\tUL = %.17e m = %.17e pc\n\tUT = %.5e s = %.5e yr\n\tUV = %.5e m/s = %.5e km/s\n\tG = %.5e\n",
	 UM,UM/MSUN,UL,UL/PARSEC,UT,UT/YEAR,UV,UV/1e3,G);
  VPRINT(stdout,"\n");

  ////////////////////////////////////////////////////
  //GLOBAL ALLOCATION
  ////////////////////////////////////////////////////
  nsysp=6*Npart;
  nsys=6;

  Ntimesp=10000;
  Ntimes=2;

  hstep=1e4*YEAR/UT;

  x=(double*)malloc(6*sizeof(double));//LSR STATE VECTOR
  xg=(double*)malloc(6*sizeof(double));//GC STATE VECTOR

  xIntp0=(double*)malloc(nsysp*sizeof(double));
  xInt0=(double*)malloc(nsys*sizeof(double));

  xInt=(double**)malloc(Ntimes*sizeof(double*));
  for(int j=0;j<Ntimes;j++) xInt[j]=(double*)malloc(nsys*sizeof(double));

  xIntp=(double**)malloc(Ntimesp*sizeof(double*));
  for(int j=0;j<Ntimesp;j++) xIntp[j]=(double*)malloc(nsysp*sizeof(double));

  xIntc=(double**)malloc(Ntimesp*sizeof(double*));
  for(int j=0;j<Ntimesp;j++) xIntc[j]=(double*)malloc(nsysp*sizeof(double));
  
  char **fields=(char**)malloc(MAXCOLS*sizeof(char*));
  for(int i=0;i<MAXCOLS;i++) fields[i]=(char*)malloc(MAXTEXT*sizeof(char));

  tsp=(double*)malloc(Ntimesp*sizeof(double));
  ts=(double*)malloc(Ntimesp*sizeof(double));

  ////////////////////////////////////////////////////
  //GLOBAL PROPERTIES FOR INTEGRATION
  ////////////////////////////////////////////////////
  ip=1;
  params[ip++]=G*MDISK*MSUN/UM;
  params[ip++]=ADISK*PARSEC/UL;
  params[ip++]=BDISK*PARSEC/UL;
  params[ip++]=G*MBULGE*MSUN/UM;
  params[ip++]=0.0;
  params[ip++]=BBULGE*PARSEC/UL;
  params[ip++]=G*MHALO*MSUN/UM;
  params[ip++]=0.0;
  params[ip++]=BHALO*PARSEC/UL;

  ////////////////////////////////////////////////////
  //READING TEST PARTICLES
  ////////////////////////////////////////////////////
  VPRINT(stdout,"Reading test particles\n",Npart);
  VPRINT(stdout,"\n");
  FILE *fc=fopen("cloud.data","r");
  for(int i=0;i<Npart;i++){
    ip=6*i;
    for(int j=0;j<=39;j++) fscanf(fc,"%lf",&tmp);
    for(int j=0;j<6;j++) fscanf(fc,"%lf",&x[j]);//READ LSR COORDINATES
    for(int j=46;j<=58;j++) fscanf(fc,"%lf",&tmp);
    //TRANSFORM COORDINATES FROM LSR->GC
    LSR2GC(x,xg);
    vscl_c(1e3/UL,xg,xg);//SET UNITS
    vscl_c(1e3/UV,xg+3,xg+3);
    //CONVERT TO CYLINDRICAL GALACTIC COORDINATES
    cart2polar(xg,xIntp0+ip,1.0);
  }
  fclose(fc);

  ////////////////////////////////////////////////////
  //INTEGRATING TEST PARTICLES FOR A LONG PERIOD 
  ////////////////////////////////////////////////////
  //TIME CONDITIONS
  duration=-1e8*YEAR/UT;

  //INTEGRATE
  params[0]=nsysp;
  integrateEoM(0,xIntp0,hstep,Ntimesp,duration,
	       nsysp,EoMGalactic,params,
	       tsp,xIntp);

  //SAVE POSITIONS IN FILE
  fc=fopen("cloud-int.csv","w");

  //HEADER
  fprintf(fc,"t,");
  for(int j=0;j<Npart;j++){
    fprintf(fc,"part%d-R,part%d-phi,part%d-Z,part%d-vR,part%d-dphi,part%d-vZ,",j,j,j,j,j,j);
    fprintf(fc,"part%d-x,part%d-y,part%d-z,part%d-vx,part%d-xy,part%d-vz,",j,j,j,j,j,j);
  }
  fprintf(fc,"dummy\n");

  //SAVE PARTICLE POSITIONS IN POLAR AND CARTESIAN
  for(int i=0;i<Ntimesp;i++){
    fprintf(fc,"%.5e,",tsp[i]);
    for(int j=0;j<Npart;j++){
      ip=6*j;
      fprintf(fc,"%s",vec2strn(xIntp[i]+ip,6,"%.17e,"));
      polar2cart(xIntp[i]+ip,xIntc[i]+ip,1.0);
      fprintf(fc,"%s",vec2strn(xIntc[i]+ip,6,"%.17e,"));
    }
    fprintf(fc,"\n");
  }
  fclose(fc);

  return 0;
}
