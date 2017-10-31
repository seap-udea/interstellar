#include <gravray.cpp>
using namespace std;

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
  //COMPUTE GALACTIC TIDAL FORCE
  ////////////////////////////////////////////////////
  
  ////////////////////////////////////////////////////
  //INITIAL POSITIONS
  ////////////////////////////////////////////////////
  fprintf(stdout,"Initial conditions:\n");
  //GRAVITATIONAL CONSTANT
  double n=ini_n*DEG/DAY;
  double a=ini_a*(-AU/1E3);
  double mu=n*n*a*a*a;
  fprintf(stdout,"\tmu = %.17e\n",mu);

  //EPHEMERIS TIME
  double dt;
  double to=unitim_c(ini_to_jed,"JDTDB","TDB");
  deltet_c(to,"et",&dt);
  fprintf(stdout,"\tto = %.17e\n",to);
  
  //INITIAL ELEMENTS
  SpiceDouble elements[]={ini_q*AU/1e3,ini_e,ini_i*DEG,ini_W*DEG,ini_w*DEG,
			  ini_M*DEG,to,mu};

  //INITIAL POSITION OF THE OBJECT RELATIVE TO THE SUN
  SpiceDouble position[6];
  conics_c(elements,to,position);

  fprintf(stdout,"\tInitial position @ Sun (EJ2000) : %s\n",
	  vec2strn(position,3,"%.17e "));
  fprintf(stdout,"\tInitial velocity @ Sun (EJ2000) : %s\n",
	  vec2strn(position+3,3,"%.17e "));

  //INITIAL POSITION RELATIVE TO THE SSB
  SpiceDouble sun[6],ltmp;
  spkezr_c("SUN",to,"ECLIPJ2000","NONE",SSB,sun,&ltmp);
  vaddg_c(position,sun,6,position);

  fprintf(stdout,"\tInitial position @ SSB (EJ2000) : %s\n",
	  vec2strn(position,3,"%.17e "));
  fprintf(stdout,"\tInitial velocity @ SSB (EJ2000) : %s\n",
	  vec2strn(position+3,3,"%.17e "));

  //COORDINATES IN J2000
  double d,RA,DEC;
  SpiceDouble M[3][3],posJ2000[3];
  pxform_c("ECLIPJ2000","J2000",to,M);
  mxv_c(M,position,posJ2000);
  recrad_c(posJ2000,&d,&RA,&DEC);
  
  fprintf(stdout,"\tD = %.5f\n",d*1e3/AU);
  fprintf(stdout,"\tRA(+HH:MM:SS) = %s, DEC(DD:MM:SS) = %s\n",dec2sex(RA*RAD/15.0),dec2sex(DEC*RAD));

  ////////////////////////////////////////////////////
  //CONVERSION TO GALACTIC COORDINATES
  ////////////////////////////////////////////////////
  double posGalactic[6];
  double lam,phi;
  pxform_c("ECLIPJ2000","GALACTIC",to,M);
  mxv_c(M,position,posGalactic);
  recrad_c(posGalactic,&d,&lam,&phi);
  fprintf(stdout,"\tD = %.5f AU\n",d*1e3/AU);
  fprintf(stdout,"\tz (galactic) = %.5f AU\n",posGalactic[2]*1e3/AU);
  fprintf(stdout,"\tLong.(+DD:MM:SS) = %s, Lat.(DD:MM:SS) = %s\n",dec2sex(lam*RAD),dec2sex(phi*RAD));

  ////////////////////////////////////////////////////
  //INITIAL CONDITIONS FOR INTEGRATION
  ////////////////////////////////////////////////////
  int npoints=2;
  double tini=to;
  double duration=-100.0*365.25*GSL_CONST_MKSA_DAY;
  double direction=duration/fabs(duration);
  double params[]={6,0};
  int eparams[]={6,1};
  
  //UNITS
  UL=GSL_CONST_MKSA_ASTRONOMICAL_UNIT;
  UM=MSUN;
  GGLOBAL=1.0;
  UT=sqrt(UL*UL*UL/(GCONST*UM));
  UV=UL/UT;

  //INITIAL CONDITIONS
  double X0[6],X[6],Xu[6],E[8],as;
  vpack_c(position[0]*1E3/UL,position[1]*1E3/UL,position[2]*1E3/UL,X0);
  vpack_c(position[3]*1E3/UV,position[4]*1E3/UV,position[5]*1E3/UV,X0+3);
  as=vnorm_c(X0);

  //DYNAMICAL TIMESCALE
  double tdyn=2*M_PI*sqrt(as*as*as/(GGLOBAL*MSUN/UM));

  //TIME LIMITS
  duration/=UT;
  double h=direction*tdyn/1000.0,h_used,h_next,h_adjust,deltat;
  double t_start=tini/UT;
  double t_step=duration/(npoints-1);
  double tend=t_start+duration;
  double t_stop=tend;
  double t=t_start;

  ////////////////////////////////////////////////////
  //COMPUTE SOLAR SYSTEM FORCES
  ////////////////////////////////////////////////////
  double dydt[6],uF[3],F,ur[3],r;
  EoM(t_start,X0,dydt,eparams);
  unorm_c(X0,ur,&r);
  fprintf(stdout,"\tSolar direction : %s, r = %.17e\n",
	  vec2strn(ur,3,"%.17e "),r);
  unorm_c(dydt+3,uF,&F);
  fprintf(stdout,"\tForces : %s\n",
	  vec2strn(dydt+3,3,"%.17e "));
  fprintf(stdout,"\tForce direction : %s, F = %.17e\n",
	  vec2strn(uF,3,"%.17e "),F);

  ////////////////////////////////////////////////////
  //GALACTIC TIDAL FORCE
  ////////////////////////////////////////////////////
  double rhodisk=0.10*MSUN/(PARSEC*PARSEC*PARSEC);
  double zdisk,adisk;
  zdisk=posGalactic[2]*1e3;
  adisk=-4*M_PI*GCONST*rhodisk*zdisk;
  adisk=adisk/(UL/(UT*UT));
  fprintf(stdout,"\trhodisk = %.17e\n",rhodisk);
  fprintf(stdout,"\tadisk = %.17e\n",adisk);

  ////////////////////////////////////////////////////
  //COMPUTE SOLAR SYSTEM BARYCENTER ORBITAL ELEMENTS
  ////////////////////////////////////////////////////
  double elts[8];
  vscl_c(UL/1E3,X0,Xu);vscl_c(UV/1E3,X0+3,Xu+3);
  oscelt_c(Xu,t*UT,MUTOT,elts);
  fprintf(stdout,"\tOrbital elements :\n");
  fprintf(stdout,"\t\tq = %.17e\n",elts[0]*1E3/AU); 
  fprintf(stdout,"\t\te = %.17e\n",elts[1]); 
  fprintf(stdout,"\t\ti = %.17e\n",elts[2]*RAD); 
  fprintf(stdout,"\t\tW = %.17e\n",elts[3]*RAD); 
  fprintf(stdout,"\t\tw = %.17e\n",elts[4]*RAD); 
  fprintf(stdout,"\t\tM = %.17e\n",elts[5]*RAD); 
  fprintf(stdout,"\t\tto = %.17e\n",elts[6]); 
  fprintf(stdout,"\t\tmu = %.17e\n",elts[7]); 
  
  ////////////////////////////////////////////////////
  //BACKWARD INTEGRATION
  ////////////////////////////////////////////////////
  //INTEGRATION
  fprintf(stdout,"Final conditions:\n");
  int i,status;
  h_used=h;
  for(i=0;i<npoints;i++) {
    deltat=(t-tini/UT)*UT/YEAR;
    //STOP INTEGRATION
    if(direction*((t_start+t_step)-tend)>0) t_step=(tend-t_start);
    t_stop = t_start + t_step;
    //STEP INTEGRATION
    h_used = h;
    do {
      while(1){
	status=Gragg_Bulirsch_Stoer(EoM,X0,X,t,h_used,&h_next,1.0,TOLERANCE,EXTMET,params);
	if(status) h_used/=4.0;
	else break;
      }
      t+=h_used;
      copyVec(X0,X,6);
      if(direction*(t+h_next-t_stop)>0) h_used=t+h_next-t_stop;
      else h_used=h_next;
    }while(direction*(t-(t_stop-direction*1.e-10))<0);
    //LAST INTEGRATION
    if(direction*(t-t_stop)>0){
      h_adjust=(t_stop-t);
      status=Gragg_Bulirsch_Stoer(EoM,X0,X,t,h_adjust,&h_next,1.0,TOLERANCE,EXTMET,params);
      copyVec(X0,X,6);
      t=t_stop;
    }
    t_start = t;
    if(direction*(t_start-tend)>0) break;
  }
  //ORBITAL UNITS
  double te=tini+deltat*YEAR;
  SpiceChar utc[100];
  deltet_c(te,"et",&dt);
  et2utc_c(te+dt,"C",3,100,utc);
  vscl_c(UL/1E3,X0,Xu);vscl_c(UV/1E3,X0+3,Xu+3);

  fprintf(stdout,"\tte = %.17e\n",te);
  fprintf(stdout,"\tte (UTC) = %s\n",utc);
  fprintf(stdout,"\tFinal position (EJ2000) : %s\n",
	  vec2strn(Xu,3,"%.17e "));
  fprintf(stdout,"\tFinal velocity (EJ2000) : %s\n",
	  vec2strn(Xu+3,3,"%.17e "));

  ////////////////////////////////////////////////////
  //FORCES
  ////////////////////////////////////////////////////
  EoM(t_start,X0,dydt,eparams);
  unorm_c(X0,ur,&r);
  fprintf(stdout,"\tSolar direction : %s, r = %.17e\n",
	  vec2strn(ur,3,"%.17e "),r);
  unorm_c(dydt+3,uF,&F);
  fprintf(stdout,"\tForces : %s\n",
	  vec2strn(dydt+3,3,"%.17e "));
  fprintf(stdout,"\tForce direction : %s, F = %.17e\n",
	  vec2strn(uF,3,"%.17e "),F);
  
  ////////////////////////////////////////////////////
  //COMPUTE SOLAR SYSTEM BARYCENTER ORBITAL ELEMENTS
  ////////////////////////////////////////////////////
  oscelt_c(Xu,te,MUTOT,elts);
  fprintf(stdout,"\tOrbital elements @ SSB :\n");
  fprintf(stdout,"\t\tq = %.17e\n",elts[0]*1E3/AU); 
  fprintf(stdout,"\t\te = %.17e\n",elts[1]); 
  fprintf(stdout,"\t\ti = %.17e\n",elts[2]*RAD); 
  fprintf(stdout,"\t\tW = %.17e\n",elts[3]*RAD); 
  fprintf(stdout,"\t\tw = %.17e\n",elts[4]*RAD); 
  fprintf(stdout,"\t\tM = %.17e\n",elts[5]*RAD); 
  fprintf(stdout,"\t\tto = %.17e\n",elts[6]); 
  fprintf(stdout,"\t\tmu = %.17e\n",elts[7]); 

  ////////////////////////////////////////////////////
  //PREDICT FUTURE POSITION USING THE ASYMPTOTIC EL.
  ////////////////////////////////////////////////////
  double tfut=to-100*YEAR;
  //double tfut=to+duration;
  conics_c(elts,tfut,position);
  fprintf(stdout,"\t\tto - tfut= %.17e\n",elts[6]-tfut); 
  fprintf(stdout,"\tFuture position @ SSB (EJ2000), t = %.17e : %s\n",
	  tfut,vec2strn(position,3,"%.17e "));
  fprintf(stdout,"\tFuture velocity @ SSB (EJ2000), t = %.17e : %s\n",
	  tfut,vec2strn(position+3,3,"%.17e "));

  ////////////////////////////////////////////////////
  //CONVERSION TO J2000 COORDINATES
  ////////////////////////////////////////////////////
  pxform_c("ECLIPJ2000","J2000",to,M);
  mxv_c(M,Xu,posJ2000);
  recrad_c(posJ2000,&d,&RA,&DEC);
  fprintf(stdout,"\tD = %.5f AU\n",d*1e3/AU);
  fprintf(stdout,"\tz = %.5f AU\n",Xu[2]*1e3/AU);
  fprintf(stdout,"\tRA(+HH:MM:SS) = %s, DEC(DD:MM:SS) = %s\n",dec2sex(RA*RAD/15.0),dec2sex(DEC*RAD));

  ////////////////////////////////////////////////////
  //CONVERSION TO GALACTIC COORDINATES
  ////////////////////////////////////////////////////
  pxform_c("ECLIPJ2000","GALACTIC",to,M);
  mxv_c(M,position,posGalactic);
  recrad_c(posGalactic,&d,&lam,&phi);
  fprintf(stdout,"\tD = %.5f AU\n",d*1e3/AU);
  fprintf(stdout,"\tz (galactic) = %.5f AU\n",posGalactic[2]*1e3/AU);
  fprintf(stdout,"\tLong.(+DD:MM:SS) = %s, Lat.(DD:MM:SS) = %s\n",dec2sex(lam*RAD),dec2sex(phi*RAD));

  ////////////////////////////////////////////////////
  //GALACTIC TIDAL FORCE
  ////////////////////////////////////////////////////
  zdisk=posGalactic[2]*1e3;
  adisk=-4*M_PI*GCONST*rhodisk*zdisk;
  adisk=adisk/(UL/(UT*UT));
  fprintf(stdout,"\trhodisk = %.17e\n",rhodisk);
  fprintf(stdout,"\tadisk = %.17e\n",adisk);

}
