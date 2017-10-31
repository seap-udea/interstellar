import numpy as np
import re

f=open("initials.conf")
qsolution=0
qelements=0
qcovariance=0

solution=""
elements=""
covariance=""

isolution=-1
ielements=-1
icovariance=-1
for line in f:
    if re.match("^\s*$",line):continue
    line=line.replace("\n","")
    #HEADER
    if "#*" in line:
        if "SOLUTION" in line:
            qcovariance=0
            qelements=0
            qsolution=1
        elif "ELEMENTS" in line:
            qcovariance=0
            qelements=1
            qsolution=0
        else:
            qcovariance=1
            qelements=0
            qsolution=0
        continue

    if qsolution:
        isolution+=1
        solution+="%s\n"%line
    if qelements:
        ielements+=1
        parts=line.split("\t")
        if ielements==0:
            parts=line.split(" ")
            elements+="double ini_to_jed=%.17e;\n"%(float(parts[4]))
            elements+="SpiceChar ini_to_date[]=\"%s\";\n"%(parts[5].replace("(","").replace(")",""))
        if ielements==1:
            elements+="double ini_e=%.17e,ini_de=%.17e;\n"%(float(parts[1]),
                                                    float(parts[2]))
        elif ielements==2:
            elements+="double ini_a=%.17e,ini_da=%.17e;//AU\n"%(float(parts[1]),
                                                        float(parts[2]))            
        elif ielements==3:
            elements+="double ini_q=%.17e,ini_dq=%.17e;//AU\n"%(float(parts[1]),
                                                        float(parts[2]))            
        elif ielements==4:
            elements+="double ini_i=%.17e,ini_di=%.17e;//deg\n"%(float(parts[1]),
                                                         float(parts[2]))
        elif ielements==5:
            elements+="double ini_W=%.17e,ini_dW=%.17e;//deg\n"%(float(parts[1]),
                                                         float(parts[2]))
        elif ielements==6:
            elements+="double ini_w=%.17e,ini_dw=%.17e;//deg\n"%(float(parts[1]),
                                                         float(parts[2]))
        elif ielements==7:
            elements+="double ini_M=%.17e,ini_dM=%.17e;//deg\n"%(float(parts[1]),
                                                         float(parts[2]))
        elif ielements==8:
            elements+="double ini_tp=%.17e;//JD\n"%(float(parts[1]))
        elif ielements==9:
            elements+="SpiceChar ini_tp_date[]=\"%s\";//JED\n"%(parts[0].replace("(","").replace(")",""))
            elements+="double ini_tp_jed=%.17e;//JED\n"%(float(parts[1]))
        elif ielements==14:
            elements+="double ini_n=%.17e,ini_dn=%.17e;//deg/d\n"%(float(parts[1]),
                                                           float(parts[2]))
    if qcovariance:
        icovariance+=1
        if icovariance==0:
            covariance+="double ini_cov[][6]={\n";
            continue
        if icovariance>0:
            parts=line.split("\t")
            linea=",".join(parts[1:])+",\n"
            covariance+=linea

covariance=covariance.strip(",")+"};\n";
f.close()

#STORE INFORMATION
f=open("initials.hpp","w")
f.write("/*\n%s\n*/\n"%solution)
f.write(elements)
f.write(covariance)
f.close()
