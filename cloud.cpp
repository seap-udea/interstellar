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
  //ELEMENTS
  double q,e,inc,W,w,Mo;
  //POSITION
  SpiceDouble position[6],sun[6],ltmp;
  SpiceDouble M_Eclip_J2000[3][3],M_Eclip_Galactic[3][3];
  SpiceDouble posJ2000[6],RA,DEC,d;
  SpiceDouble posGalactic[6],l,b;
  SpiceDouble posFuture[6],RAfut,DECfut,lfut,bfut,dfut;
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
  //LOOP OVER PARTICLES
  ////////////////////////////////////////////////////
  FILE* fc=fopen("cloud.data","w");
  int Npart=50;
  int Nfreq=ceil(Npart/10);
  int j;
  for(j=0;j<Npart;j++){

    if ((j%Nfreq)==0)
      fprintf(stdout,"Particle %d...\n",j+1);

    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //GENERATE INITIAL ELEMENTS
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    double fac=1.0;
    //q
    q=ini_q+gsl_ran_gaussian(RAND,ini_dq*fac);
    if(VERBOSE) fprintf(stdout,"\tq = %.17e (mean = %.17e, sigma = %.17e)\n",q,ini_q,ini_dq);
    //e
    e=ini_e+gsl_ran_gaussian(RAND,ini_de*fac);
    if(VERBOSE) fprintf(stdout,"\te = %.17e (mean = %.17e, sigma = %.17e)\n",e,ini_e,ini_de);
    //i
    inc=ini_i+gsl_ran_gaussian(RAND,ini_di*fac);
    if(VERBOSE) fprintf(stdout,"\ti = %.17e (mean = %.17e, sigma = %.17e)\n",inc,ini_i,ini_di);
    //W
    W=ini_W+gsl_ran_gaussian(RAND,ini_dW*fac);
    if(VERBOSE) fprintf(stdout,"\tW = %.17e (mean = %.17e, sigma = %.17e)\n",W,ini_W,ini_dW);
    //w
    w=ini_w+gsl_ran_gaussian(RAND,ini_dw*fac);
    if(VERBOSE) fprintf(stdout,"\tw = %.17e (mean = %.17e, sigma = %.17e)\n",w,ini_w,ini_dw);
    //Mo
    Mo=ini_M+gsl_ran_gaussian(RAND,ini_dM*fac);
    if(VERBOSE) fprintf(stdout,"\tM = %.17e (mean = %.17e, sigma = %.17e)\n",Mo,ini_M,ini_dM); 
    SpiceDouble elements[]={q*AU/1e3,e,inc*DEG,W*DEG,w*DEG,Mo*DEG,to,mu};

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
    recrad_c(posJ2000,&d,&RA,&DEC);
    if(VERBOSE) fprintf(stdout,"\tRA(+HH:MM:SS) = %s, DEC(DD:MM:SS) = %s, d = %.3lf AU\n",
	    dec2sex(RA*RAD/15.0),dec2sex(DEC*RAD),d*1e3/AU);

    //GALACTIC
    pxform_c("ECLIPJ2000","GALACTIC",t,M_Eclip_Galactic);
    mxv_c(M_Eclip_Galactic,Xu,posGalactic);
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
    mxv_c(M_Eclip_J2000,posFuture,posJ2000);
    recrad_c(posJ2000,&dfut,&RAfut,&DECfut);
    if(VERBOSE) fprintf(stdout,"\tFuture: RA(+HH:MM:SS) = %s, DEC(DD:MM:SS) = %s, d = %.3lf pc\n",
	    dec2sex(RAfut*RAD/15.0),dec2sex(DECfut*RAD),dfut*1e3/PARSEC);

    //GALACTIC
    pxform_c("ECLIPJ2000","GALACTIC",tfut,M_Eclip_Galactic);
    mxv_c(M_Eclip_Galactic,posFuture,posGalactic);
    recrad_c(posGalactic,&dfut,&lfut,&bfut);
    if(VERBOSE) fprintf(stdout,"\tFuture: l(+DD:MM:SS) = %s, b(DD:MM:SS) = %s\n",
	    dec2sex(lfut*RAD),dec2sex(bfut*RAD));

    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //SAVE POSITION
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //PARTICLE TEND
    fprintf(fc,"%-10d %-27.17e ",j,tend);
    //POSITION ECLIPJ2000
    fprintf(fc,"%s ",vec2strn(Xu,6,"%.17e "));
    //POSITION J2000
    fprintf(fc,"%s ",vec2strn(posJ2000,3,"%.17e "));
    //POSITION GALACTIC
    fprintf(fc,"%s ",vec2strn(posGalactic,3,"%.17e "));
    //J2000
    fprintf(fc,"%-27.17e %-27.17e ",RA*RAD/15.0,DEC*RAD);
    //GALACTIC
    fprintf(fc,"%-27.17e %-27.17e ",l*RAD,b*RAD);
    //DISTANCE
    fprintf(fc,"%-27.17e ",d*1E3/AU);
    //ASYMPTOTIC ELEMENTS
    fprintf(fc,"%s ",vec2strn(elts,8,"%.17e "));
    //DISTANT FUTURE
    fprintf(fc,"%s ",vec2strn(posFuture,8,"%.17e "));
    //J2000
    fprintf(fc,"%-27.17e %-27.17e ",RAfut*RAD/15.0,DECfut*RAD);
    //GALACTIC
    fprintf(fc,"%-27.17e %-27.17e ",lfut*RAD,bfut*RAD);
    //DISTANCE
    fprintf(fc,"%-27.17e ",dfut*1E3/PARSEC);

    fprintf(fc,"\n");
    //break;
  }
  fclose(fc);
}
