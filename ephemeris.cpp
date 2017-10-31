#include <gravray.cpp>
using namespace std;
#define MU 1.32712440017987106e+11

int main(int argc,char* argv[])
{
  ////////////////////////////////////////////////////
  //INITIALIZE CSPICE
  ////////////////////////////////////////////////////
  initSpice();
  FDEBUG=fopen("debug.dat","w");

  ////////////////////////////////////////////////////
  //DATE
  ////////////////////////////////////////////////////
  SpiceDouble et,dt,tjd,tjdb,tdb;
  str2et_c("10/30/2017 00:00:00.000 UTC",&et);
  deltet_c(et,"et",&dt);
  tjd=t2jd(et);
  fprintf(stdout,"Julian Date at ET = %.6lf\n",tjd);
  tdb=et-dt;
  tjdb=t2jd(tdb);
  fprintf(stdout,"Julian Date at TDB = %.6lf\n",tjdb);

  ////////////////////////////////////////////////////
  //ELEMENTS OF THE COMET
  ////////////////////////////////////////////////////
  //Solution 2017-Oct-29
  double todb=2458052.5;//2017-Oct-26.0
  double e=1.197188708990351,de=0.0016908;
  double a=-1.290422401198905,da=0.0077688;//AU
  double q=0.2544567273446408,dq=0.0006499;//AU
  double i=122.604016492383,di=0.060915;//deg
  double W=24.60295925659798,dW=0.0027382;//deg
  double w=241.5429238264925,dw=0.11715;//deg
  double M=31.28328035977056,dM=0.27399;//deg
  double tp=2458005.972892582579;
  double n=0.6723667577077285,dn=0.0060718;//deg/d
  str2et_c("10/26/2017 00:00:00.000",&et);
  deltet_c(et,"et",&dt);
  tdb=et-dt;
  tjdb=t2jd(tdb);
  fprintf(stdout,"Julian Date at TDB = %.6lf\n",tjdb);
  double to=tdb;

  fprintf(stdout,"mu=%.17e\n",n*DEG/DAY*n*DEG/DAY*a*a*(-a)*(AU/1e3)*(AU/1e3)*(AU/1e3));
  fprintf(stdout,"mu = %.17e\n",GMASSES[0]);
  fprintf(stdout,"to = %.17e\n",to);

  SpiceDouble position[6],sun[6],ltmp;
  SpiceDouble elements[]={q*AU/1e3,e,i*DEG,W*DEG,w*DEG,
			  M*DEG,to,/**/MU/*GMASSES[0]*/};

  //MOON
  /*
  spkezr_c("301",et,"ECLIPJ2000","NONE",SSB,position,&ltmp);
  spkezr_c("301",et+365.25*DAY,"ECLIPJ2000","NONE",SSB,position,&ltmp);
  */

  spkezr_c("8",et,"ECLIPJ2000","NONE",SSB,sun,&ltmp);

  //COMET
  spkezr_c("SUN",et,"ECLIPJ2000","NONE",SSB,sun,&ltmp);
  fprintf(stdout,"Sun = %s\n",vec2strn(sun,6,"%.17e "));
  conics_c(elements,to,position);
  vaddg_c(position,sun,6,position);
  /*
    1.89270050566111952e+08 8.36792936384781599e+07 4.24613912213196978e+06 4.10740430517903334e+01 8.67467036631901145e+00 1.44039404989131761e+01
    1.32712440040944000e+11
  */

  //MOON
  fprintf(stdout,"%.10e %s 1.0 2\n",et,vec2strn(position,6,"%.17e "));

  //VELOCITY
  fprintf(stdout,"v = %.17e\n",vnorm_c(position+3));

  et+=100*365.25*DAY;
  deltet_c(et,"et",&dt);

  SpiceChar utc[100];
  et2utc_c(et,"C",3,100,utc);
  fprintf(stdout,"UTC=%s\n",utc);

  return 0;
}
