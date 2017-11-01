#include <gravray.cpp>
using namespace std;

#define VERBOSE 0

int main(int argc,char* argv[])
{
  ////////////////////////////////////////////////////
  //INITIALIZE CSPICE
  ////////////////////////////////////////////////////
  initSpice();
  FDEBUG=fopen("debug.dat","w");

  int Npart=atoi(argv[1]);

  ////////////////////////////////////////////////////
  //INITIAL CONDITIONS
  ////////////////////////////////////////////////////
#include <initials.hpp>

  ////////////////////////////////////////////////////
  //GLOBAL DECLARATIONS
  ////////////////////////////////////////////////////
  //UNITS
  UL=GSL_CONST_MKSA_ASTRONOMICAL_UNIT;
  UM=MSUN;
  GGLOBAL=1.0;
  UT=sqrt(UL*UL*UL/(GCONST*UM));
  UV=UL/UT;
  //GRAVITATIONAL PARAMETER
  double n=ini_n*DEG/DAY;
  double a=ini_a*(-AU/1E3);
  double mu=n*n*a*a*a;
  //EPHEMERIS TIME
  double dt;
  double to=unitim_c(ini_to_jed,"JDTDB","TDB");
  double tp;
  //ELEMENTS
  double q,e,inc,W,w,Mo;
  //POSITION
  SpiceDouble position[6],sun[6],ltmp;
  SpiceDouble M_Eclip_J2000[3][3],M_Eclip_Galactic[3][3];
  SpiceDouble posJ2000[6],RA,DEC,d;
  SpiceDouble posGalactic[6],l,b;
  SpiceDouble posFutureJ2000[6],posFuture[6],posFutureGalactic[6];
  double RAfut,DECfut,lfut,bfut,dfut;
  //INTEGRATION
  int npoints=2;
  double tini=to;
  double duration=-100.0*YEAR;
  double direction=duration/fabs(duration);
  double params[]={6,0};
  duration/=UT;
  double h,t_start,t_step,tend,t_stop,t;
  double h_used,h_next,h_adjust,deltat;
  int i,status;
  //INITIAL CONDITIONS
  double X0[6],X[6],Xu[6],as,tdyn;
  //ELEMENTS
  double elts[8];
  //FUTURE
  double tfut=to-5e5*YEAR;

  ////////////////////////////////////////////////////
  //PREPARE ELEMENTS GENERATION
  ////////////////////////////////////////////////////
  SpiceDouble elements[6];
  gsl_vector* ielements=gsl_vector_alloc(6);
  gsl_vector_set(ielements,0,ini_e);
  gsl_vector_set(ielements,1,ini_q);
  gsl_vector_set(ielements,2,ini_tp);
  gsl_vector_set(ielements,3,ini_W);
  gsl_vector_set(ielements,4,ini_w);
  gsl_vector_set(ielements,5,ini_i);
  gsl_vector* relements=gsl_vector_alloc(6);
  gsl_matrix* Lo=gsl_matrix_alloc(6,6);
  gsl_matrix* L=gsl_matrix_alloc(6,6);
  for(int i=0;i<6;i++) for(int j=0;j<6;j++) gsl_matrix_set(Lo,i,j,ini_cov[i][j]);

  ////////////////////////////////////////////////////
  //LOOP OVER PARTICLES
  ////////////////////////////////////////////////////
  FILE* fc=fopen("cloud.data","w");
  int Nfreq=ceil(Npart/10);
  Nfreq=Nfreq==0?1:Nfreq;
  
  for(int j=0;j<Npart;j++){
    
    if ((j%Nfreq)==0)
      fprintf(stdout,"Particle %d...\n",j+1);

    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //GENERATE INITIAL ELEMENTS (MULTIVAR)
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //GENERATE ELEMENTS BY GAUSSAN MULTIVARIATE
    gsl_matrix_memcpy(L,Lo);
    gsl_linalg_cholesky_decomp1(L);
    gsl_ran_multivariate_gaussian(RAND,ielements,L,relements);
    n=ini_n+gsl_ran_gaussian(RAND,ini_dn);
    tp=gsl_vector_get(relements,2);
    Mo=n*(ini_to_jed-tp);

    //EXTRACT ELEMENTS
    /*q=*/elements[0]=gsl_vector_get(relements,1)*AU/1e3;
    /*e=*/elements[1]=gsl_vector_get(relements,0);
    /*i=*/elements[2]=gsl_vector_get(relements,5)*DEG;
    /*W=*/elements[3]=gsl_vector_get(relements,3)*DEG;
    /*w=*/elements[4]=gsl_vector_get(relements,4)*DEG;
    /*M=*/elements[5]=Mo*DEG;
    /*to=*/elements[6]=to;
    /*mu=*/elements[7]=mu;

    /*
    fprintf(stdout,"q = %.17e\n",elements[0]);
    fprintf(stdout,"e = %.17e\n",elements[1]);
    fprintf(stdout,"i = %.17e\n",elements[2]*RAD);
    fprintf(stdout,"W = %.17e\n",elements[3]*RAD);
    fprintf(stdout,"w = %.17e\n",elements[4]*RAD);
    fprintf(stdout,"M = %.17e\n",elements[5]*RAD);
    fprintf(stdout,"to = %.17e\n",elements[6]);
    fprintf(stdout,"mu = %.17e\n",elements[7]);
    exit(0);
    */

    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //INITIAL POSITION @ SSB
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    conics_c(elements,to,position);
    if(VERBOSE) fprintf(stdout,"\tInitial position @ Sun (EJ2000) : %s\n",
			vec2strn(position,3,"%.17e "));
    if(VERBOSE) fprintf(stdout,"\tInitial velocity @ Sun (EJ2000) : %s\n",
			vec2strn(position+3,3,"%.17e "));
    
    spkezr_c("SUN",to,"ECLIPJ2000","NONE",SSB,sun,&ltmp);
    vaddg_c(position,sun,6,position);
    if(VERBOSE) fprintf(stdout,"\tInitial position @ SSB (EJ2000) : %s\n",
	    vec2strn(position,3,"%.17e "));
    if(VERBOSE) fprintf(stdout,"\tInitial velocity @ SSB (EJ2000) : %s\n",
	    vec2strn(position+3,3,"%.17e "));

    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //INITIAL CONDITIONS
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    vpack_c(position[0]*1E3/UL,position[1]*1E3/UL,position[2]*1E3/UL,X0);
    vpack_c(position[3]*1E3/UV,position[4]*1E3/UV,position[5]*1E3/UV,X0+3);
    as=vnorm_c(X0);
    tdyn=2*M_PI*sqrt(as*as*as/(GGLOBAL*MSUN/UM));
    h=direction*tdyn/1000.0;
    t_start=tini/UT;
    t_step=duration/(npoints-1);
    tend=t_start+duration;
    t_stop=tend;
    t=t_start;

    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //INTEGRATION
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    h_used=h;
    for(i=0;i<npoints;i++) {
      deltat=(t-tini/UT)*UT/YEAR;
      if(direction*((t_start+t_step)-tend)>0) t_step=(tend-t_start);
      t_stop = t_start + t_step;
      h_used = h;
      do {
	while(1){
	  status=Gragg_Bulirsch_Stoer(EoM,X0,X,t,h_used,&h_next,1.0,
				      TOLERANCE,EXTMET,params);
	  if(status) h_used/=4.0;
	  else break;
	}
	t+=h_used;
	copyVec(X0,X,6);
	if(direction*(t+h_next-t_stop)>0) h_used=t+h_next-t_stop;
	else h_used=h_next;
      }while(direction*(t-(t_stop-direction*1.e-10))<0);
      if(direction*(t-t_stop)>0){
	h_adjust=(t_stop-t);
	status=Gragg_Bulirsch_Stoer(EoM,X0,X,t,h_adjust,&h_next,1.0,
				    TOLERANCE,EXTMET,params);
	copyVec(X0,X,6);
	t=t_stop;
      }
      t_start = t;
      if(direction*(t_start-tend)>0) break;
    }

    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //FINAL POSITION
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    vscl_c(UL/1E3,X0,Xu);vscl_c(UV/1E3,X0+3,Xu+3);
    if(VERBOSE) fprintf(stdout,"\tFinal position : %s\n",
	    vec2strn(Xu,3,"%.17e "));

    //ASYMPTOTIC ELEMENTS
    oscelt_c(Xu,t,MUTOT,elts);
    
    //J2000
    pxform_c("ECLIPJ2000","J2000",t,M_Eclip_J2000);
    mxv_c(M_Eclip_J2000,Xu,posJ2000);
    mxv_c(M_Eclip_J2000,Xu+3,posJ2000+3);
    recrad_c(posJ2000,&d,&RA,&DEC);
    if(VERBOSE) fprintf(stdout,"\tRA(+HH:MM:SS) = %s, DEC(DD:MM:SS) = %s, d = %.3lf AU\n",
	    dec2sex(RA*RAD/15.0),dec2sex(DEC*RAD),d*1e3/AU);

    //GALACTIC
    pxform_c("ECLIPJ2000","GALACTIC",t,M_Eclip_Galactic);
    mxv_c(M_Eclip_Galactic,Xu,posGalactic);
    mxv_c(M_Eclip_Galactic,Xu+3,posGalactic+3);
    recrad_c(posGalactic,&d,&l,&b);
    if(VERBOSE) fprintf(stdout,"\tl(+DD:MM:SS) = %s, b(DD:MM:SS) = %s\n",
	    dec2sex(l*RAD),dec2sex(b*RAD));

    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //POSITION IN THE DISTANT FUTURE
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    conics_c(elts,tfut,posFuture);
    if(VERBOSE) fprintf(stdout,"\tDistant future position : %s\n",
			vec2strn(posFuture,3,"%.17e "));

    //J2000
    pxform_c("ECLIPJ2000","J2000",tfut,M_Eclip_J2000);
    mxv_c(M_Eclip_J2000,posFuture,posFutureJ2000);
    recrad_c(posFutureJ2000,&dfut,&RAfut,&DECfut);
    if(VERBOSE) fprintf(stdout,"\tFuture: RA(+HH:MM:SS) = %s, DEC(DD:MM:SS) = %s, d = %.3lf pc\n",
	    dec2sex(RAfut*RAD/15.0),dec2sex(DECfut*RAD),dfut*1e3/PARSEC);

    //GALACTIC
    pxform_c("ECLIPJ2000","GALACTIC",tfut,M_Eclip_Galactic);
    mxv_c(M_Eclip_Galactic,posFuture,posFutureGalactic);
    mxv_c(M_Eclip_Galactic,posFuture+3,posFutureGalactic+3);
    recrad_c(posFutureGalactic,&dfut,&lfut,&bfut);
    if(VERBOSE) fprintf(stdout,"\tFuture: l(+DD:MM:SS) = %s, b(DD:MM:SS) = %s\n",
	    dec2sex(lfut*RAD),dec2sex(bfut*RAD));

    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //SAVE POSITION
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //PARTICLE TEND
    fprintf(fc,"%-10d %-27.17e %-27.17e",j,tend*UT,tfut);
    //POSITION ECLIPJ2000
    fprintf(fc,"%s ",vec2strn(Xu,6,"%.17e "));
    //POSITION J2000
    fprintf(fc,"%s ",vec2strn(posJ2000,6,"%.17e "));
    //POSITION GALACTIC
    fprintf(fc,"%s ",vec2strn(posGalactic,6,"%.17e "));
    //J2000
    fprintf(fc,"%-27.17e %-27.17e ",RA*RAD/15.0,DEC*RAD);
    //GALACTIC
    fprintf(fc,"%-27.17e %-27.17e ",l*RAD,b*RAD);
    //DISTANCE
    fprintf(fc,"%-27.17e ",d*1E3/AU);
    //ASYMPTOTIC ELEMENTS
    fprintf(fc,"%s ",vec2strn(elts,8,"%.17e "));
    //DISTANT FUTURE
    fprintf(fc,"%s ",vec2strn(posFuture,6,"%.17e "));
    //DISTANT FUTURE
    fprintf(fc,"%s ",vec2strn(posFutureGalactic,6,"%.17e "));
    //J2000
    fprintf(fc,"%-27.17e %-27.17e ",RAfut*RAD/15.0,DECfut*RAD);
    //GALACTIC
    fprintf(fc,"%-27.17e %-27.17e ",lfut*RAD,bfut*RAD);
    //DISTANCE
    fprintf(fc,"%-27.17e ",dfut*1E3/PARSEC);
    //ELEMENTOS
    fprintf(fc,"%s ",vec2strn(elements,8,"%.17e "));

    fprintf(fc,"\n");
  }
  fclose(fc);
}
