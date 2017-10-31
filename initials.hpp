/*
# obs. used (total)   	  86
   data-arc span   	  15 days
   first obs. used   	  2017-10-14
   last obs. used   	  2017-10-29
   planetary ephem.   	  DE431
   SB-pert. ephem.   	  SB431-N16
   fit RMS   	  .35829
   data source   	  ORB
   producer   	  Davide Farnocchia
   solution date   	  2017-Oct-29 10:45:18

*/
double ini_to_jed=2.45805250000000000e+06;
SpiceChar ini_to_date[]="2017-Oct-26.0";
double ini_e=1.19718870899035101e+00,ini_de=1.69079999999999990e-03;
double ini_a=-1.29042240119890494e+00,ini_da=7.76880000000000018e-03;//AU
double ini_q=2.54456727344640810e-01,ini_dq=6.49900000000000022e-04;//AU
double ini_i=1.22604016492382996e+02,ini_di=6.09149999999999969e-02;//deg
double ini_W=2.46029592565979804e+01,ini_dW=2.73820000000000010e-03;//deg
double ini_w=2.41542923826492512e+02,ini_dw=1.17150000000000004e-01;//deg
double ini_M=3.12832803597705613e+01,ini_dM=2.73990000000000011e-01;//deg
double ini_tp=2.45800597289258242e+06;//JD
SpiceChar ini_tp_date[]="2017-Sep-09.47289258";//JED
double ini_tp_jed=1.26619999999999996e-02;//JED
double ini_n=6.72366757707728468e-01,ini_dn=6.07179999999999959e-03;//deg/d
double ini_cov[][6]={
2.858709169167452E-6,1.098820139532213E-6,2.140740994127999E-5,-4.629574614074441E-6,.0001980724106465366,.0001029927307342494,
1.098820139532213E-6,4.223650568138116E-7,8.228257068002674E-6,-1.779505280075431E-6,7.613474148207401E-5,3.958804146112113E-5,
2.140740994127999E-5,8.228257068002674E-6,.0001603227078400257,-3.466805081263181E-5,.001483242149313819,.0007712542878842905,
-4.629574614074441E-6,-1.779505280075431E-6,-3.466805081263181E-5,7.497524956627533E-6,-.0003207714690047103,-.000166792920108298,
.0001980724106465366,7.613474148207401E-5,.001483242149313819,-.0003207714690047103,.01372394843890766,.007136100317461406,
.0001029927307342494,3.958804146112113E-5,.0007712542878842905,-.000166792920108298,.007136100317461406,.003710639376110803,
};
