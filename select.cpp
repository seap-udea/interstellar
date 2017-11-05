#include <gravray.cpp>
using namespace std;

#define VERBOSE 1

int main(int argc,char* argv[])
{
  /*
    Example: ./select.exe 

    Input: 
    * cloud.data
    * stars_mot.data
      Rows: First are names of columns
      Cols:
          0: ID
	  1: ref_epoch 
	  2: ra (deg, ref_epoch) 
	  3: ra_error (mas) 
	  4: dec (deg, ref_epoch)
	  5: dec_error (mas)
	  6: parallax (mas)
	  7: parallax_error (mas)
	  8: pmra (mas/year)
	  9: pmra_error (mas/year)
	  10: pmdec (mas/year)
	  11: pmdec_error (mas/year)
	  12: phot_g_mean_mag (gband, mag)
	  13: l (Galactic longitude)
	  14: b (Galactic latitude)
	  15: distance (pc)

    Output: 
  */

  ////////////////////////////////////////////////////
  //INITIALIZE CSPICE
  ////////////////////////////////////////////////////
  int Npart=1;
  if(argc>1){
    Npart=atoi(argv[1]);
  }

  ////////////////////////////////////////////////////
  //INITIALIZE CSPICE
  ////////////////////////////////////////////////////
  initSpice();
  
  ////////////////////////////////////////////////////
  //GENERAL VARIABLES
  ////////////////////////////////////////////////////
  double tmp;
  char ctmp[100];

  ////////////////////////////////////////////////////
  //READING DATA
  ////////////////////////////////////////////////////

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //READING NOMINAL SOLUTION
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  FILE *fc=fopen("cloud.data","r");
  double posbody[6],tbody;
  for(int i=0;i<=0;i++) fscanf(fc,"%lf",&tmp);
  fscanf(fc,"%lf",&tbody);
  VPRINT(stdout,"Epoch body: %.17e\n",tbody);
  for(int i=2;i<=8;i++) fscanf(fc,"%lf",&tmp);
  for(int j=0;j<6;j++) fscanf(fc,"%lf",&posbody[j]);
  VPRINT(stdout,"Position body: %s\n",vec2strn(posbody,6,"%.17e "));
  fclose(fc);

  ////////////////////////////////////////////////////
  //COMPUTING MINIMUM DISTANCE
  ////////////////////////////////////////////////////
  fc=fopen("stars_mot.data","r");
  char id[100];
  double tstar,dt;
  double raep,decep;
  double muraep,dmura,mudecep,dmudec;
  double ra,dra,dec,ddec;
  double postar[6];
  double gmag,gMag;
  double lep,bep,l,b;
  double par,dpar,d,dd;
  double M_Epoch_J2000[3][3];
  double M_J2000_Galactic[3][3];

  for(int i=1;i<=15;i++) fscanf(fc,"%s",&ctmp);
  int n=0;
  int Nstars;
  while(fscanf(fc,"%s",&id)==1){

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //PRIMARY
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //EPOCH
    fscanf(fc,"%lf",&tstar);
    sprintf(ctmp,"01/01/%d 00:00:00.000",(int)tstar);
    str2et_c(ctmp,&tstar);
    deltet_c(tstar,"et",&dt);
    tstar-=dt;

    VPRINT(stdout,"id %d: %s, epoch = %s, t = %lf\n",n,id,ctmp,tstar);

    //COORDINATES AT EPOCH
    fscanf(fc,"%lf",&ra);fscanf(fc,"%lf",&dra);
    fscanf(fc,"%lf",&dec);fscanf(fc,"%lf",&ddec);
    VPRINT(stdout,"\tRA(epoch) = %.17lf +/- %.3lf mas\n",ra,dra);
    VPRINT(stdout,"\tDEC(epoch) = %.17lf +/- %.3lf mas\n",dec,ddec);
    VPRINT(stdout,"\tRA(epoch) = %s, DEC(epoch) = %s\n",dec2sex(ra/15.0),dec2sex(dec));

    //PARALLAX
    fscanf(fc,"%lf",&par);
    fscanf(fc,"%lf",&dpar);
    VPRINT(stdout,"\tParallax = %.17lf +/- %.3lf mas\n",par,dpar);

    //PROPER MOTION
    fscanf(fc,"%lf",&muraep);fscanf(fc,"%lf",&dmura);
    fscanf(fc,"%lf",&mudecep);fscanf(fc,"%lf",&dmudec);
    VPRINT(stdout,"\tmu_RA(epoch) = %.17lf +/- %.3lf mas\n",muraep,dmura);
    VPRINT(stdout,"\tmu_DEC(epoch) = %.17lf +/- %.3lf mas\n",mudecep,dmudec);

    //GMAG
    fscanf(fc,"%lf",&gmag);
    VPRINT(stdout,"\tmag_g = %.3lf\n",gmag);
    
    //READ GALACTIC COORDINATES
    fscanf(fc,"%lf",&lep);
    fscanf(fc,"%lf",&bep);
    VPRINT(stdout,"\tl = %.17lf, b = %.17lf\n",lep,bep);
    
    //DISTANCE
    fscanf(fc,"%lf",&d);
    dd=d*dpar/par;
    VPRINT(stdout,"\td(pc) = %.17lf +/- %.3lf \n",d,dd);

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //SECONDARY
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    VPRINT(stdout,"\tSecondary properties:\n",gMag);
    
    //ABSOLUTE MAGNITUDE
    gMag=gmag-5*log10(d/10);
    VPRINT(stdout,"\tAbsolute magnitude = %.3lf\n",gMag);

    //POSITION TO STAR RESPECT TO J2000
    radrec_c(d,ra*DEG,dec*DEG,postar);
    VPRINT(stdout,"\tPosition in true epoch = %s\n",vec2str(postar,"%.17e "));

    /*
    //TRANSFORM FROM EPOCH TO J2000
    pxform_c("EARTHTRUEEPOCH","J2000",tstar,M_Epoch_J2000);
    mxv_c(M_Epoch_J2000,postar,postar);
    VPRINT(stdout,"\tPosition in J2000 = %s\n",vec2str(postar,"%.17e "));

    //J2000 COORDINATES
    recrad_c(postar,&tmp,&ra,&dec);
    VPRINT(stdout,"\tRA = %s, DEC = %s\n",dec2sex(ra*RAD/15.0),dec2sex(dec*RAD));
    */

    //TRANSFORM TO GALACTIC TO CHECK
    str2et_c("01/01/2000 00:00:00.000",&tstar);
    deltet_c(tstar,"et",&dt);
    tstar-=dt;
    pxform_c("J2000","GALACTIC",tstar,M_J2000_Galactic);
    mxv_c(M_J2000_Galactic,postar,postar);
    VPRINT(stdout,"\tPosition in Galactic = %s\n",vec2str(postar,"%.17e "));

    //GALACTIC COORDINATES
    recrad_c(postar,&tmp,&l,&b);
    VPRINT(stdout,"\tl = %.17lf, b = %.17lf\n",l*RAD,b*RAD);

    if(n>-1) break;
    n++;
  }
  fclose(fc);
  
  Nstars=n+1;
  VPRINT(stdout,"Number of stars: %d\n",Nstars);
  return 0;
}
