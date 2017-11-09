/*
  Transform state vector from cartesian (x,y,z) to polar coordinates
  (R,q,z)

  Example: if rho is in parsec and vrho in km/s units=1e3/PARSEC
 */
int cart2polar(double x[6],double p[6],double units=1)
{
  double *r=x,*v=x+3;
  double *rho=p,*vrho=p+3;

  //State vector in cartesian
  reccyl_c(r,&rho[0],&rho[1],&rho[2]);
  double cosq=cos(rho[1]),sinq=sin(rho[1]);
  double IM[][3]={{cosq,sinq,0},{-sinq,cosq,0},{0,0,1}};
  mxv_c(IM,v,vrho);
  vrho[1]=vrho[1]/rho[0]*units;
  
  return 0;
}
/*
  Transform state vector from polar (R,q,z) to cartesian (x,y,z)

  Example: if rho is in parsec and vrho in km/s units=1e3/PARSEC
 */
int polar2cart(double p[6],double x[6],double units=1.0)
{
  double *r=x,*v=x+3;
  double *rho=p,*vrho=p+3;
  double cosq=cos(rho[1]),sinq=sin(rho[1]);
  double M[][3]={{cosq,-sinq,0},{sinq,cosq,0},{0,0,1}};

  cylrec_c(rho[0],rho[1],rho[2],r);

  vrho[1]=vrho[1]*rho[0]/units;
  mxv_c(M,vrho,v);
  vrho[1]=vrho[1]/rho[0]*units;

  return 0;
}

  params[P_MUD]=MDISK*MSUN/UM;params[P_AD]=ADISK*PARSEC/UL;params[P_BD]=BDISK*PARSEC/UL;
  params[P_MUB]=MBULGE*MSUN/UM;params[P_AB]=0.0;params[P_BB]=BBULGE*PARSEC/UL;
  params[P_MUH]=MHALO*MSUN/UM;params[P_AH]=0.0;params[P_BH]=BHALO*PARSEC/UL;




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
    /*
    SpiceDouble elements[]={q*AU/1e3,e,inc*DEG,W*DEG,w*DEG,Mo*DEG,to,mu};
    */

    0:i
    1:tdb (terminal)
    2:tdb (future)
    3-8:Position Ecliptic J2000
    9-14:Position J2000
    15-20:Position Galactic J2000
    21:RA(h) (terminal)
    22:DEC(deg)
    23:l(deg)
    24:b(deg)
    25:d(AU)
    26-33:Asymptotic elements, q,e,i,W,w,Mo,to,mu
    34-39:Future Position Ecliptic J2000
    40-45:Future Position Galactic
    46:RA(h) (future)
    47:DEC(deg)
    48:l(deg)
    49:b(deg)
    50:d(pc)
    51-58:Initial elements, q,e,i,W,w,Mo,to,mu
