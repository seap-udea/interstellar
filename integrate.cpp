#include <gravray.cpp>
using namespace std;

#define VERBOSE 1

//FIELDS OF THE CANDIDATES.CSV FILE
enum {N,POSTARX,POSTARY,POSTARZ,VELSTARX,VELSTARY,VELSTARZ,POSBODYPERIX,POSBODYPERIY,POSBODYPERIZ,POSTARPERIX,POSTARPERIY,POSTARPERIZ,DMIN,TMIN,VRELX,VRELY,VRELZ,VREL,HIP,TYCHO2_ID,SOLUTION_ID,SOURCE_ID,RANDOM_INDEX,REF_EPOCH,RA,RA_ERROR,DEC,DEC_ERROR,PARALLAX,PARALLAX_ERROR,PMRA,PMRA_ERROR,PMDEC,PMDEC_ERROR,RA_DEC_CORR,RA_PARALLAX_CORR,RA_PMRA_CORR,RA_PMDEC_CORR,DEC_PARALLAX_CORR,DEC_PMRA_CORR,DEC_PMDEC_CORR,PARALLAX_PMRA_CORR,PARALLAX_PMDEC_CORR,PMRA_PMDEC_CORR,ASTROMETRIC_N_OBS_AL,ASTROMETRIC_N_OBS_AC,ASTROMETRIC_N_GOOD_OBS_AL,ASTROMETRIC_N_GOOD_OBS_AC,ASTROMETRIC_N_BAD_OBS_AL,ASTROMETRIC_N_BAD_OBS_AC,ASTROMETRIC_DELTA_Q,ASTROMETRIC_EXCESS_NOISE,ASTROMETRIC_EXCESS_NOISE_SIG,ASTROMETRIC_PRIMARY_FLAG,ASTROMETRIC_RELEGATION_FACTOR,ASTROMETRIC_WEIGHT_AL,ASTROMETRIC_WEIGHT_AC,ASTROMETRIC_PRIORS_USED,MATCHED_OBSERVATIONS,DUPLICATED_SOURCE,SCAN_DIRECTION_STRENGTH_K1,SCAN_DIRECTION_STRENGTH_K2,SCAN_DIRECTION_STRENGTH_K3,SCAN_DIRECTION_STRENGTH_K4,SCAN_DIRECTION_MEAN_K1,SCAN_DIRECTION_MEAN_K2,SCAN_DIRECTION_MEAN_K3,SCAN_DIRECTION_MEAN_K4,PHOT_G_N_OBS,PHOT_G_MEAN_FLUX,PHOT_G_MEAN_FLUX_ERROR,PHOT_G_MEAN_MAG,PHOT_VARIABLE_FLAG,L,B,ECL_LON,ECL_LAT,RAJ2000,DEJ2000,RV,ERV,CAT};

int main(int argc,char* argv[])
{
  /*
    Example: ./integrate.exe 10

    Where: 
    10: Number of particles

    Input: 
    * cloud.data: test particles compatible with A/2017U1 orbit.
    * candidates.csv: candidates to close encounters.

    Output:
    *
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

  Ntimesp=100;
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
  duration=-1e7*YEAR/UT;

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
  exit(0);

  ////////////////////////////////////////////////////
  //READING SELECTED OBJECTS
  ////////////////////////////////////////////////////
  params[0]=nsys;
  fc=fopen("candidates.csv","r");
  fgets(line,MAXLINE,fc);

  n=0;
  while(fgets(line,MAXLINE,fc)!=NULL){

    //PARSE FIELDS
    parseLine(line,fields,&nfields);
    n++;
    if(strcmp(fields[HIP],"40170.0")!=0) continue;
    fprintf(stdout,"Star %d,%s:\n",n,fields[HIP]);

    //LMA ENCOUNTER TIME
    duration=atof(fields[TMIN]);

    //SPATIAL COORDINATES AND LSR VELOCITY
    x[0]=atof(fields[POSTARX])*PARSEC/1e3;
    x[1]=atof(fields[POSTARY])*PARSEC/1e3;
    x[2]=atof(fields[POSTARZ])*PARSEC/1e3;
    x[3]=atof(fields[VELSTARX]);
    x[4]=atof(fields[VELSTARY]);
    x[5]=atof(fields[VELSTARZ]);

    //CONVERT TO GC
    LSR2GC(x,xg);
    vscl_c(1e3/UL,xg,xg);//SET UNITS
    vscl_c(1e3/UV,xg+3,xg+3);

    //CONVERT TO CYLINDRICAL COORDINATES
    cart2polar(xg,xInt0,1.0);

    //INTEGRATE 
    integrateEoM(0,xInt0,hstep,Ntimes,duration,
		 nsys,EoMGalactic,params,
		 ts,xInt);

    //RETURN TO CARTESIAN
    double p2[6];
    polar2cart(xInt[1],p2,1.0);
    fprintf(stdout,"\tFinal position:%s\n",vec2strn(p2,6,"%.5e "));

    /*
    double dp[6];
    vsubg_c(p2,p1,6,dp);
    fprintf(stdout,"\tDifference: [%s] (pos:%.5e,vel:%.5e]\n",
    vec2strn(dp,6,"%.5e "),vnorm_c(dp),vnorm_c(dp+3)*UV);
    */

    break;
  }

  return 0;
}
