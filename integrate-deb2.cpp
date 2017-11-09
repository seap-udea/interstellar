#include <gravray.cpp>
using namespace std;

#define VERBOSE 1

enum {N,POSTARX,POSTARY,POSTARZ,VELSTARX,VELSTARY,VELSTARZ,POSBODYPERIX,POSBODYPERIY,POSBODYPERIZ,POSTARPERIX,POSTARPERIY,POSTARPERIZ,DMIN,TMIN,VRELX,VRELY,VRELZ,VREL,HIP,TYCHO2_ID,SOLUTION_ID,SOURCE_ID,RANDOM_INDEX,REF_EPOCH,RA,RA_ERROR,DEC,DEC_ERROR,PARALLAX,PARALLAX_ERROR,PMRA,PMRA_ERROR,PMDEC,PMDEC_ERROR,RA_DEC_CORR,RA_PARALLAX_CORR,RA_PMRA_CORR,RA_PMDEC_CORR,DEC_PARALLAX_CORR,DEC_PMRA_CORR,DEC_PMDEC_CORR,PARALLAX_PMRA_CORR,PARALLAX_PMDEC_CORR,PMRA_PMDEC_CORR,ASTROMETRIC_N_OBS_AL,ASTROMETRIC_N_OBS_AC,ASTROMETRIC_N_GOOD_OBS_AL,ASTROMETRIC_N_GOOD_OBS_AC,ASTROMETRIC_N_BAD_OBS_AL,ASTROMETRIC_N_BAD_OBS_AC,ASTROMETRIC_DELTA_Q,ASTROMETRIC_EXCESS_NOISE,ASTROMETRIC_EXCESS_NOISE_SIG,ASTROMETRIC_PRIMARY_FLAG,ASTROMETRIC_RELEGATION_FACTOR,ASTROMETRIC_WEIGHT_AL,ASTROMETRIC_WEIGHT_AC,ASTROMETRIC_PRIORS_USED,MATCHED_OBSERVATIONS,DUPLICATED_SOURCE,SCAN_DIRECTION_STRENGTH_K1,SCAN_DIRECTION_STRENGTH_K2,SCAN_DIRECTION_STRENGTH_K3,SCAN_DIRECTION_STRENGTH_K4,SCAN_DIRECTION_MEAN_K1,SCAN_DIRECTION_MEAN_K2,SCAN_DIRECTION_MEAN_K3,SCAN_DIRECTION_MEAN_K4,PHOT_G_N_OBS,PHOT_G_MEAN_FLUX,PHOT_G_MEAN_FLUX_ERROR,PHOT_G_MEAN_MAG,PHOT_VARIABLE_FLAG,L,B,ECL_LON,ECL_LAT,RAJ2000,DEJ2000,RV,ERV,CAT};

int main(int argc,char* argv[])
{
  /*
    Example: ./integrate.exe 

    Input: 

    * cloud.data: test particles compatible with A/2017U1 orbit.

    * candidates.csv: candidates to close encounters.

    Output: 

    
  */
  ////////////////////////////////////////////////////
  //INPUT PARAMETERS
  ////////////////////////////////////////////////////
  int Npart=1;
  if(argc>1){
    Npart=atoi(argv[1]);
  }
  VPRINT(stdout,"Propagating %d test particles against the candidate stars\n",Npart);
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
  
  int Ntimes,nsys; 
  int ip,n;
  int nfields;
  double params[10];
  double tini,h,duration;

  //MATRICES WITH INTEGRATIONS
  double **xLSR,**xGC,**xGCp,**xInt,*ts;

  //INITIAL CONDITIONS
  double *xGC0,*xGCp0,*xInt0,*dxGCdt,*x0;

  ////////////////////////////////////////////////////
  //UNITS
  ////////////////////////////////////////////////////
  double G;
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
  G=GCONST/(UL*UL*UL/(UM*UT*UT));
  //*/

  //UNITS OF VELOCITY
  UV=UL/UT;
  
  //CORRECTION FACTOR FOR POLAR->CARTESIAN TRANSFORMATION
  double units=1; //UV/UL;
  
  VPRINT(stdout,"Units:\n\tUM = %.5e kg=%.5e Msun\n\tUL = %.17e m = %.17e pc\n\tUT = %.5e s = %.5e yr\n\tUV = %.5e m/s = %.5e km/s\n\tG = %.5e\n",
	 UM,UM/MSUN,UL,UL/PARSEC,UT,UT/YEAR,UV,UV/1e3,G);
  VPRINT(stdout,"\n");

  ////////////////////////////////////////////////////
  //TEST TRANSFORMATIONS
  ////////////////////////////////////////////////////
  /*
  double xLSR_s[]={-3.16028e+02*PARSEC/1e3,
		   1.84156e+01*PARSEC/1e3,
		   -3.58527e+02*PARSEC/1e3, 
                   1.33746e+01, -1.72967e+01, -1.54271e+01};
  VPRINT(stdout,"LSR (star) = %s\n",vec2strn(xLSR_s,6,"%.17e "));
  double xGC_s[6];
  
  //CONVERT FROM LSR TO GC
  LSR2GC(xLSR_s,xGC_s);
  vscl_c(1e3/PARSEC,xGC_s,xGC_s);
  VPRINT(stdout,"GC (star) = %s (pc,km/s)\n",vec2strn(xGC_s,6,"%.17e "));

  //CONVERT FROM CARTESIAN TO POLAR
  double xGC_p[6],units=1e3/PARSEC;
  cart2polar(xGC_s,xGC_p,units);
  VPRINT(stdout,"GC (polar) = %s (pc,km/s)\n",vec2strn(xGC_p,6,"%.17e "));

  polar2cart(xGC_p,xGC_s,units);
  VPRINT(stdout,"GC (cart.) = %s (pc,km/s)\n",vec2strn(xGC_s,6,"%.17e "));
  VPRINT(stdout,"GC (polar) = %s (pc,km/s)\n",vec2strn(xGC_p,6,"%.17e "));
  exit(0);
  */

  ////////////////////////////////////////////////////
  //READING TEST PARTICLES
  ////////////////////////////////////////////////////
  Ntimes=2;
  nsys=6*Npart;
  ts=(double*)malloc(Ntimes*sizeof(double));
  
  FILE *fc=fopen("cloud.data","r");
  
  //ALLOCATE INTEGRATION RESULTS
  xLSR=(double**)malloc(Ntimes*sizeof(double*));
  for(int j=0;j<Ntimes;j++) xLSR[j]=(double*)malloc(nsys*sizeof(double));
  xGC0=(double*)malloc(nsys*sizeof(double));
  xGCp0=(double*)malloc(nsys*sizeof(double));

  //READ PARTICLES
  for(int i=0;i<Npart;i++){
    ip=6*i;
    for(int j=0;j<=39;j++) fscanf(fc,"%lf",&tmp);
    for(int j=0;j<6;j++) fscanf(fc,"%lf",&xLSR[0][ip+j]);
    for(int j=46;j<=58;j++) fscanf(fc,"%lf",&tmp);
    
    //TRANSFORM COORDINATES FROM LSR->GC
    LSR2GC(xLSR[0]+ip,xGC0+ip);
    //copyVec(xGC0+ip,xLSR[0]+ip,6);//Use this line to test

    vscl_c(1e3/UL,xGC0+ip,xGC0+ip);//SET UNITS
    vscl_c(1e3/UV,xGC0+ip+3,xGC0+ip+3);

    //CONVERT TO CYLINDRICAL COORDINATES
    cart2polar(xGC0+ip,xGCp0+ip,units);
  }
  fclose(fc);

  ////////////////////////////////////////////////////
  //GALACTIC POTENTIAL
  ////////////////////////////////////////////////////
  params[0]=nsys;
  ip=1;
  params[ip++]=G*MDISK*MSUN/UM;params[ip++]=ADISK*PARSEC/UL;params[ip++]=BDISK*PARSEC/UL;
  params[ip++]=G*MBULGE*MSUN/UM;params[ip++]=0.0;params[ip++]=BBULGE*PARSEC/UL;
  params[ip++]=G*MHALO*MSUN/UM;params[ip++]=0.0;params[ip++]=BHALO*PARSEC/UL;
  
  /*
  double R=ROSUN/UL*PARSEC,z=ZSUN*PARSEC/UL,dphidR,dphidq,dphidz;
  gradGalacticPotential(R,z,&dphidR,&dphidq,&dphidz,params);
  VPRINT(stdout,"Parameters : %s\n",vec2strn(params,10,"%.5e "));
  fprintf(stdout,"\tR = %.17e, z = %.17e, dphidR = %.17e, dphidq = %.17e, dphidz = %.17e\n",
	  R,z,dphidR,dphidq,dphidz);
  //CONVERT GRADIENT OF POTENTIAL TO SI
  fprintf(stdout,"dphi/dz = %.17e UL/UT^2\n",dphidz);
  dphidz=dphidz*UL/(UT*UT);
  fprintf(stdout,"dphi/dz = %.17e m/s^2\n",dphidz);
  exit(0);
  */

  ////////////////////////////////////////////////////
  //INTEGRATING TEST PARTICLES
  ////////////////////////////////////////////////////
  tini=0.0;
  h=1e4*YEAR/UT;
  duration=-5.629030e6*YEAR/UT;

  //ALLOCATE INTEGRATION RESULTS
  xGC=(double**)malloc(Ntimes*sizeof(double*));
  for(int j=0;j<Ntimes;j++) xGC[j]=(double*)malloc(nsys*sizeof(double));

  xGCp=(double**)malloc(Ntimes*sizeof(double*));
  for(int j=0;j<Ntimes;j++) xGCp[j]=(double*)malloc(nsys*sizeof(double));

  dxGCdt=(double*)malloc(nsys*sizeof(double));
  x0=(double*)malloc(nsys*sizeof(double));

  //CHOOSE COORDINATE SYSTEM
  int qpolar=1;
  if(qpolar){xInt0=xGCp0;xInt=xGCp;}
  else{xInt0=xGC0;xInt=xGC;}

  //GRADIENT AT INITIAL CONDITIONS
  EoMGalactic(0,xInt0,dxGCdt,params);
  fprintf(stdout,"\tInitial vector (cart.):%s\n",vec2strn(xGC0,nsys,"%.5e "));
  fprintf(stdout,"\tInitial vector (polar):%s\n",vec2strn(xGCp0,nsys,"%.5e "));
  fprintf(stdout,"\tInitial gradient:%s\n",vec2strn(dxGCdt,nsys,"%.5e "));

  //INTEGRATE
  fprintf(stdout,"Integrating %d test particles...\n",Npart);
  integrateEoM(0,xInt0,h,Ntimes,duration,
	       nsys,EoMGalactic,params,
	       ts,xInt);
  fprintf(stdout,"\tResult:%s\n",vec2strn(xInt[1],nsys,"%.5e "));

  if(!qpolar){
    //TRANSFORM TO POLAR
    double polar[nsys];
    for(int i=0;i<Npart;i++){
      ip=6*i;
      cart2polar(xInt[1]+ip,polar+ip,units);
    }
    fprintf(stdout,"\tResult (polar):%s\n",vec2strn(polar,nsys,"%.5e "));
  }

  //STORE RESULT
  //Convert from polar to cartesian
  double p1[6];
  if(qpolar) polar2cart(xInt[1],p1,units);
  else copyVec(p1,xGC[1],6);
  fprintf(stdout,"\tFinal position:%s\n",vec2strn(p1,6,"%.5e "));
  fprintf(stdout,"\n");

  ////////////////////////////////////////////////////
  //READING SELECTED OBJECTS
  ////////////////////////////////////////////////////
  Ntimes=2;
  nsys=6;

  //COMMONG SETTINGS INTEGRATION
  params[0]=nsys;
  tini=0.0;
  h=1e4*YEAR/UT;
  duration=duration;

  ts=(double*)malloc(Ntimes*sizeof(double));

  fc=fopen("candidates.csv","r");
  fgets(line,MAXLINE,fc);

  char **fields=(char**)malloc(MAXCOLS*sizeof(char*));
  for(int i=0;i<MAXCOLS;i++) fields[i]=(char*)malloc(MAXTEXT*sizeof(char));

  xLSR=(double**)malloc(Ntimes*sizeof(double*));
  for(int j=0;j<Ntimes;j++) xLSR[j]=(double*)malloc(nsys*sizeof(double));
  xGC0=(double*)malloc(nsys*sizeof(double));
  xGCp0=(double*)malloc(nsys*sizeof(double));

  xGC=(double**)malloc(Ntimes*sizeof(double*));
  for(int j=0;j<Ntimes;j++) xGC[j]=(double*)malloc(nsys*sizeof(double));
  
  xGCp=(double**)malloc(Ntimes*sizeof(double*));
  for(int j=0;j<Ntimes;j++) xGCp[j]=(double*)malloc(nsys*sizeof(double));

  n=0;
  while(fgets(line,MAXLINE,fc)!=NULL){
    //PARSE FIELDS
    parseLine(line,fields,&nfields);
    n++;
    if(strcmp(fields[HIP],"40170.0")!=0) continue;
    fprintf(stdout,"Star %d,%s:\n",n,fields[HIP]);

    //SPATIAL COORDINATES AND LSR VELOCITY
    xLSR[0][0]=atof(fields[POSTARX])*PARSEC/1e3;
    xLSR[0][1]=atof(fields[POSTARY])*PARSEC/1e3;
    xLSR[0][2]=atof(fields[POSTARZ])*PARSEC/1e3;
    xLSR[0][3]=atof(fields[VELSTARX]);
    xLSR[0][4]=atof(fields[VELSTARY]);
    xLSR[0][5]=atof(fields[VELSTARZ]);
    VPRINT(stdout,"\tInitial contitions (LSR):%s\n",vec2strn(xLSR[0],6,"%.5e "));

    //CONVERT TO GC
    LSR2GC(xLSR[0],xGC0);
    //copyVec(xGC0,xLSR[0],6);//Use this line to test

    vscl_c(1e3/UL,xGC0,xGC0);//SET UNITS
    vscl_c(1e3/UV,xGC0+3,xGC0+3);

    //CONVERT TO CYLINDRICAL COORDINATES
    cart2polar(xGC0,xGCp0,units);

    VPRINT(stdout,"\tInitial contitions (cart.):%s\n",vec2strn(xGC0,6,"%.5e "));
    VPRINT(stdout,"\tInitial contitions (polar):%s\n",vec2strn(xGCp0,6,"%.5e "));

    if(qpolar){xInt0=xGCp0;xInt=xGCp;}
    else{xInt0=xGC0;xInt=xGC;}

    //INTEGRATE 
    integrateEoM(0,xInt0,h,Ntimes,duration,
		 nsys,EoMGalactic,params,
		 ts,xInt);

    //RESULT
    VPRINT(stdout,"\tResult (GC):%s\n",vec2strn(xInt[1],nsys,"%.5e "));

    double p2[6];
    if(qpolar) polar2cart(xInt[1],p2,units);
    else copyVec(p2,xGC[1],6);
    fprintf(stdout,"\tFinal position:%s\n",vec2strn(p2,6,"%.5e "));

    double dp[6];
    vsubg_c(p2,p1,6,dp);
    fprintf(stdout,"\tDifference: [%s] (pos:%.5e,vel:%.5e]\n",
	    vec2strn(dp,6,"%.5e "),vnorm_c(dp),vnorm_c(dp+3)*UV);

    break;
  }

  return 0;
}
