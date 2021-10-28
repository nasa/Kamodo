#import t89
import numpy as np
from geopack import geopack
#from geopack import t89,t96,t01,t04
import os
import datetime
from datetime import datetime,timezone,timedelta

def seconds_from_19700101(Time): #=datetime(2000,1,1,0,0,0,tzinfo=timezone.utc)):
    dt=Time-datetime(1970,1,1,0,0,0,tzinfo=timezone.utc)
    return dt.total_seconds();
                              
def transform(Xvec=np.array([[0.,0.,1.]]),
              Time=np.array([datetime(2000,1,1,0,0,0,tzinfo=timezone.utc)]),
              coord_in='GSE',
              coord_out='GEO',
              debug=False,
              **kwargs):
# no operation for same input and output coordinates
    if coord_in == coord_out:
        return Xvec
    Xvec=np.array(Xvec)
    Xvec_out=[]
    #        print(Xvec.shape)

# show whether we are using the correct time
    if debug:
        print(Time[0])
        print(seconds_from_19700101(Time[0]))
        psi=geopack.recalc(seconds_from_19700101(Time[0]))
        print('psi=',psi)

    if coord_in == "GSE":
        if coord_out == "GSM":
            for i in range(len(Xvec[:,0])):
                if i < len(Time):
                    psi=geopack.recalc(seconds_from_19700101(Time[i]))
                x_=Xvec[i,0]
                y_=Xvec[i,1]
                z_=Xvec[i,2]
                (x_out,y_out,z_out)=geopack.gsmgse(x_,y_,z_,-1)
                Xvec_out.append([x_out,y_out,z_out])
            return Xvec_out
                
        if coord_out == "SM":
            for i in range(len(Xvec[:,0])):
                if i < len(Time):
                    psi=geopack.recalc(seconds_from_19700101(Time[i]))
                x_=Xvec[i,0]
                y_=Xvec[i,1]
                z_=Xvec[i,2]
                (x_out,y_out,z_out)=geopack.gsmgse(x_,y_,z_,-1)
                (x_out2,y_out2,z_out2)=geopack.smgsm(x_out,y_out,z_out,-1)
                Xvec_out.append([x_out2,y_out2,z_out2])
            return Xvec_out            

        if coord_out == "GEO":
            for i in range(len(Xvec[:,0])):
                if i < len(Time):
                    psi=geopack.recalc(seconds_from_19700101(Time[i]))
                x_=Xvec[i,0]
                y_=Xvec[i,1]
                z_=Xvec[i,2]
                (x_out,y_out,z_out)=geopack.gsmgse(x_,y_,z_,-1)
                (x_out2,y_out2,z_out2)=geopack.geogsm(x_out,y_out,z_out,-1)
                Xvec_out.append([x_out2,y_out2,z_out2])
            return Xvec_out 

        if coord_out == "GEI":
            for i in range(len(Xvec[:,0])):
                if i < len(Time):
                    psi=geopack.recalc(seconds_from_19700101(Time[i]))
                x_=Xvec[i,0]
                y_=Xvec[i,1]
                z_=Xvec[i,2]
                (x_out,y_out,z_out)=geopack.gsmgse(x_,y_,z_,-1)
                (x_out2,y_out2,z_out2)=geopack.geogsm(x_out,y_out,z_out,-1)
                (x_out3,y_out3,z_out3)=geopack.geigeo(x_out2,y_out2,z_out2,-1)
                Xvec_out.append([x_out3,y_out3,z_out3])
            return Xvec_out

        if coord_out == "MAG":
            for i in range(len(Xvec[:,0])):
                if i < len(Time):
                    psi=geopack.recalc(seconds_from_19700101(Time[i]))
                x_=Xvec[i,0]
                y_=Xvec[i,1]
                z_=Xvec[i,2]
                (x_out,y_out,z_out)=geopack.gsmgse(x_,y_,z_,-1)
                (x_out2,y_out2,z_out2)=geopack.geogsm(x_out,y_out,z_out,-1)
                (x_out3,y_out3,z_out3)=geopack.geomag(x_out2,y_out2,z_out2,1)
                Xvec_out.append([x_out3,y_out3,z_out3])
            return Xvec_out

        print("target coordinate system not supported")
        return False

            
    if coord_in == "GSM":
        if coord_out == "GSE":
            for i in range(len(Xvec[:,0])):
                if i < len(Time):
                    psi=geopack.recalc(seconds_from_19700101(Time[i]))
                x_=Xvec[i,0]
                y_=Xvec[i,1]
                z_=Xvec[i,2]
                (x_out,y_out,z_out)=geopack.gsmgse(x_,y_,z_,1)
                Xvec_out.append([x_out,y_out,z_out])
            return Xvec_out
            
        if coord_out == "SM":
            for i in range(len(Xvec[:,0])):
                if i < len(Time):
                    psi=geopack.recalc(seconds_from_19700101(Time[i]))
                x_=Xvec[i,0]
                y_=Xvec[i,1]
                z_=Xvec[i,2]
                (x_out,y_out,z_out)=geopack.smgsm(x_,y_,z_,-1)
                Xvec_out.append([x_out,y_out,z_out])
            return Xvec_out

        if coord_out == "GEO":
            for i in range(len(Xvec[:,0])):
                if i < len(Time):
                    psi=geopack.recalc(seconds_from_19700101(Time[i]))
                x_=Xvec[i,0]
                y_=Xvec[i,1]
                z_=Xvec[i,2]
                (x_out,y_out,z_out)=geopack.geogsm(x_,y_,z_,-1)
                Xvec_out.append([x_out,y_out,z_out])
            return Xvec_out

        if coord_out == "GEI":
            for i in range(len(Xvec[:,0])):
                if i < len(Time):
                    psi=geopack.recalc(seconds_from_19700101(Time[i]))
                x_=Xvec[i,0]
                y_=Xvec[i,1]
                z_=Xvec[i,2]
                (x_out,y_out,z_out)=geopack.geogsm(x_,y_,z_,-1)
                (x_out2,y_out2,z_out2)=geopack.geigeo(x_out,y_out,z_out,-1)
                Xvec_out.append([x_out2,y_out2,z_out2])
            return Xvec_out
            
        if coord_out == "MAG":
            for i in range(len(Xvec[:,0])):
                if i < len(Time):
                    psi=geopack.recalc(seconds_from_19700101(Time[i]))
                x_=Xvec[i,0]
                y_=Xvec[i,1]
                z_=Xvec[i,2]
                (x_out,y_out,z_out)=geopack.gsmgse(x_,y_,z_,1)
                (x_out2,y_out2,z_out2)=geopack.geomag(x_out,y_out,z_out,1)
                Xvec_out.append([x_out2,y_out2,z_out2])
            return Xvec_out

        print("target coordinate system not supported")
        return False
           
    if coord_in == "SM":
        if coord_out == "GSM":
            for i in range(len(Xvec[:,0])):
                if i < len(Time):
                    psi=geopack.recalc(seconds_from_19700101(Time[i]))
                x_=Xvec[i,0]
                y_=Xvec[i,1]
                z_=Xvec[i,2]
                (x_out,y_out,z_out)=geopack.smgsm(x_,y_,z_,1)
                Xvec_out.append([x_out,y_out,z_out])
            return Xvec_out

        if coord_out == "GSE":
            for i in range(len(Xvec[:,0])):
                if i < len(Time):
                    psi=geopack.recalc(seconds_from_19700101(Time[i]))
                x_=Xvec[i,0]
                y_=Xvec[i,1]
                z_=Xvec[i,2]
                (x_out,y_out,z_out)=geopack.smgsm(x_,y_,z_,1)
                (x_out2,y_out2,z_out2)=geopack.gsmgse(x_out,y_out,z_out,1)
                Xvec_out.append([x_out2,y_out2,z_out2])
            return Xvec_out

        if coord_out == "MAG":
            for i in range(len(Xvec[:,0])):
                if i < len(Time):
                    psi=geopack.recalc(seconds_from_19700101(Time[i]))
                x_=Xvec[i,0]
                y_=Xvec[i,1]
                z_=Xvec[i,2]
                (x_out,y_out,z_out)=geopack.magsm(x_,y_,z_,-1)
                Xvec_out.append([x_out,y_out,z_out])
            return Xvec_out
            
        if coord_out == "GEO":
            for i in range(len(Xvec[:,0])):
                if i < len(Time):
                    psi=geopack.recalc(seconds_from_19700101(Time[i]))
                x_=Xvec[i,0]
                y_=Xvec[i,1]
                z_=Xvec[i,2]
                (x_out,y_out,z_out)=geopack.magsm(x_,y_,z_,-1)
                (x_out2,y_out2,z_out2)=geopack.geomag(x_out,y_out,z_out,-1)
                Xvec_out.append([x_out2,y_out2,z_out2])
            return Xvec_out
        
        if coord_out == "GEI":
            for i in range(len(Xvec[:,0])):
                if i < len(Time):
                    psi=geopack.recalc(seconds_from_19700101(Time[i]))
                x_=Xvec[i,0]
                y_=Xvec[i,1]
                z_=Xvec[i,2]
                (x_out,y_out,z_out)=geopack.magsm(x_,y_,z_,-1)
                (x_out2,y_out2,z_out2)=geopack.geomag(x_out,y_out,z_out,-1)
                (x_out3,y_out3,z_out3)=geopack.geigeo(x_out2,y_out2,z_out2,-1)
                Xvec_out.append([x_out3,y_out3,z_out3])
            return Xvec_out

        print("target coordinate system not supported")
        return False
                        
    if coord_in == "MAG":
        if coord_out == "SM":
            for i in range(len(Xvec[:,0])):
                if i < len(Time):
                    psi=geopack.recalc(seconds_from_19700101(Time[i]))
                x_=Xvec[i,0]
                y_=Xvec[i,1]
                z_=Xvec[i,2]
                (x_out,y_out,z_out)=geopack.magsm(x_,y_,z_,-1)
                Xvec_out.append([x_out,y_out,z_out])
            return Xvec_out
                        
        if coord_out == "GEO":
            for i in range(len(Xvec[:,0])):
                if i < len(Time):
                    psi=geopack.recalc(seconds_from_19700101(Time[i]))
                x_=Xvec[i,0]
                y_=Xvec[i,1]
                z_=Xvec[i,2]
                (x_out,y_out,z_out)=geopack.geomag(x_,y_,z_,-1)
                Xvec_out.append([x_out,y_out,z_out])
            return Xvec_out
            
        if coord_out == "GSM":
            for i in range(len(Xvec[:,0])):
                if i < len(Time):
                    psi=geopack.recalc(seconds_from_19700101(Time[i]))
                x_=Xvec[i,0]
                y_=Xvec[i,1]
                z_=Xvec[i,2]
                (x_out,y_out,z_out)=geopack.geomag(x_,y_,z_,-1)
                (x_out2,y_out2,z_out2)=geopack.geogsm(x_out,y_out,z_out,1)
                Xvec_out.append([x_out2,y_out2,z_out2])
            return Xvec_out

        if coord_out == "GSE":
            for i in range(len(Xvec[:,0])):
                if i < len(Time):
                    psi=geopack.recalc(seconds_from_19700101(Time[i]))
                x_=Xvec[i,0]
                y_=Xvec[i,1]
                z_=Xvec[i,2]
                (x_out,y_out,z_out)=geopack.geomag(x_,y_,z_,-1)
                (x_out2,y_out2,z_out2)=geopack.geogsm(x_out,y_out,z_out,1)
                (x_out3,y_out3,z_out3)=geopack.gsmgse(x_out2,y_out2,z_out2,1)
                Xvec_out.append([x_out3,y_out3,z_out3])
            return Xvec_out
            
        if coord_out == "GEI":
            for i in range(len(Xvec[:,0])):
                if i < len(Time):
                    psi=geopack.recalc(seconds_from_19700101(Time[i]))
                x_=Xvec[i,0]
                y_=Xvec[i,1]
                z_=Xvec[i,2]
                (x_out,y_out,z_out)=geopack.geomag(x_,y_,z_,-1)
                (x_out2,y_out2,z_out2)=geopack.geigeo(x_out,y_out,z_out,-1)
                Xvec_out.append([x_out2,y_out2,z_out2])
            return Xvec_out

        print("target coordinate system not supported")
        return False
            
    if coord_in == "GEO":
        if coord_out == "GEI":
            for i in range(len(Xvec[:,0])):
                if i < len(Time):
                    psi=geopack.recalc(seconds_from_19700101(Time[i]))
                x_=Xvec[i,0]
                y_=Xvec[i,1]
                z_=Xvec[i,2]
                (x_out,y_out,z_out)=geopack.geigeo(x_,y_,z_,-1)
                Xvec_out.append([x_out,y_out,z_out])
            return Xvec_out

        if coord_out == "MAG":
            for i in range(len(Xvec[:,0])):
                if i < len(Time):
                    psi=geopack.recalc(seconds_from_19700101(Time[i]))
                x_=Xvec[i,0]
                y_=Xvec[i,1]
                z_=Xvec[i,2]
                (x_out,y_out,z_out)=geopack.geomag(x_,y_,z_,1)
                Xvec_out.append([x_out,y_out,z_out])
            return Xvec_out

        if coord_out == "SM":
            for i in range(len(Xvec[:,0])):
                if i < len(Time):
                    psi=geopack.recalc(seconds_from_19700101(Time[i]))
                x_=Xvec[i,0]
                y_=Xvec[i,1]
                z_=Xvec[i,2]
                (x_out,y_out,z_out)=geopack.geomag(x_,y_,z_,1)
                (x_out2,y_out2,z_out2)=geopack.magsm(x_out,y_out,z_out,1)
                Xvec_out.append([x_out2,y_out2,z_out2])
            return Xvec_out

        if coord_out == "GSM":
            for i in range(len(Xvec[:,0])):
                if i < len(Time):
                    psi=geopack.recalc(seconds_from_19700101(Time[i]))
                x_=Xvec[i,0]
                y_=Xvec[i,1]
                z_=Xvec[i,2]
                (x_out,y_out,z_out)=geopack.geogsm(x_,y_,z_,1)
                Xvec_out.append([x_out,y_out,z_out])
            return Xvec_out

        if coord_out == "GSE":
            for i in range(len(Xvec[:,0])):
                if i < len(Time):
                    psi=geopack.recalc(seconds_from_19700101(Time[i]))
                x_=Xvec[i,0]
                y_=Xvec[i,1]
                z_=Xvec[i,2]
                (x_out,y_out,z_out)=geopack.geogsm(x_,y_,z_,1)
                (x_out2,y_out2,z_out2)=geopack.gsmgse(x_out,y_out,z_out,1)
                Xvec_out.append([x_out2,y_out2,z_out2])
            return Xvec_out

        print("target coordinate system not supported")
        return False

            
    if coord_in == "GEI":
        if coord_out == "GEO":
            for i in range(len(Xvec[:,0])):
                if i < len(Time):
                    psi=geopack.recalc(seconds_from_19700101(Time[i]))
                x_=Xvec[i,0]
                y_=Xvec[i,1]
                z_=Xvec[i,2]
                (x_out,y_out,z_out)=geopack.geigeo(x_,y_,z_,1)
                Xvec_out.append([x_out,y_out,z_out])
            return Xvec_out

        if coord_out == "MAG":
            for i in range(len(Xvec[:,0])):
                if i < len(Time):
                    psi=geopack.recalc(seconds_from_19700101(Time[i]))
                x_=Xvec[i,0]
                y_=Xvec[i,1]
                z_=Xvec[i,2]
                (x_out,y_out,z_out)=geopack.geigeo(x_,y_,z_,1)
                (x_out2,y_out2,z_out2)=geopack.geomag(x_out,y_out,z_out,1)
                Xvec_out.append([x_out,y_out,z_out])
            return Xvec_out

        if coord_out == "SM":
            for i in range(len(Xvec[:,0])):
                if i < len(Time):
                    psi=geopack.recalc(seconds_from_19700101(Time[i]))
                x_=Xvec[i,0]
                y_=Xvec[i,1]
                z_=Xvec[i,2]
                (x_out,y_out,z_out)=geopack.geigeo(x_,y_,z_,1)
                (x_out2,y_out2,z_out2)=geopack.geomag(x_,y_,z_,1)
                (x_out3,y_out3,z_out3)=geopack.magsm(x_out2,y_out2,z_out2,1)
                Xvec_out.append([x_out3,y_out3,z_out3])
            return Xvec_out

        if coord_out == "GSM":
            for i in range(len(Xvec[:,0])):
                if i < len(Time):
                    psi=geopack.recalc(seconds_from_19700101(Time[i]))
                x_=Xvec[i,0]
                y_=Xvec[i,1]
                z_=Xvec[i,2]
                (x_out,y_out,z_out)=geopack.geigeo(x_,y_,z_,1)
                (x_out2,y_out2,z_out2)=geopack.geogsm(x_out,y_out,z_out,1)
                Xvec_out.append([x_out2,y_out2,z_out2])
            return Xvec_out
        
        if coord_out == "GSE":
            for i in range(len(Xvec[:,0])):
                if i < len(Time):
                    psi=geopack.recalc(seconds_from_19700101(Time[i]))
                x_=Xvec[i,0]
                y_=Xvec[i,1]
                z_=Xvec[i,2]
                (x_out,y_out,z_out)=geopack.geigeo(x_,y_,z_,1)
                (x_out2,y_out2,z_out2)=geopack.geogsm(x_out,y_out,z_out,1)
                (x_out3,y_out3,z_out3)=geopack.gsmgse(x_out2,y_out2,z_out2,1)
                Xvec_out.append([x_out3,y_out3,z_out3])
            return Xvec_out
        
        print("target coordinate system not supported")
        return False
        
    print("origin coordinate system not supported")
    return False

