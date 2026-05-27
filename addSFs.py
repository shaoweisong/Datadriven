#/cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc11-opt/setup.sh
import awkward as ak
import numpy as np
import pandas as pd
import sys
sys.path.append("/eos/user/s/shsong/pkgs_condor/parquet_to_root-0.3.0/")
from parquet_to_root import parquet_to_root
from math import *

def load_parquet(fname):
    print("loading events from %s" % fname)
    events=ak.from_parquet(fname)
    events=events[events.Diphoton_minID_modified>-0.7]
    return events
def add_sale_factor(event,sclae_factor):
    event['weight_central']=sclae_factor*event.weight_central
    return event
######################################################################################
# ppsf=[1.41931,1.17995,7.55792,1.30833] #2017                                       #          
# ddsf=[0.990377,0.820495,0.65972,1.01199] #2017                                     #            
# ppsf=[1.70958,1.11875,7.46335,1.32682] #2018                                       #          
# ddsf=[0.752818,0.759332,0.830304,1.014] #2018                                      #           
# ppsf= [0.767877,0.629229,3.27182,0.778788] #2016pre datadriven file use 2017       #                                          
# ddsf= [0.506172,0.60107,0.555076,0.695838] #2016pre                                # HH                
# ppsf=[0.814871,0.692412,4.98402,0.913607] #2016post                                #                 
# ddsf=[0.785703,0.866525,0.504,1.02589] #2016post  cat1 use2017dd cat3 2018         #                                        
# ppsf=[5.17817,3.45248,21.2556,4.58472] #2016                                       #          
# ddsf=[0.902465,0.908307,1.09993,1.04744] #2016                                     #            
######################################################################################


######################################################################################
######################################################################################
#2016pre                                                                             #  
ppsf_16pre_highmY =[1.32499,0.985829,5.83325,1.2458]                                                          #          
ddsf_16pre_highmY =[0.494714,0.855295,0.592225,0.981137]
#2016post                                                                            #  
ppsf_16post_highmY =[1.45603,1.01543,7.41129,1.32791]
ddsf_16post_highmY =[0.976013,0.86819,1.10732,1.02984]
#2017                                                                                # YH
ppsf_17_highmY =[1.42568,1.14754,7.49827,1.35126]                                                             #          
ddsf_17_highmY =[0.70013,0.841792,0.770761,1.00546]                                                           #           
#2018                                                                                #
ppsf_18_highmY =[1.85362,1.20059,7.25905,1.4003]
ddsf_18_highmY =[0.65432,0.813623,0.979891,1.01301]
######################################################################################
years = ['16pre','16post','17','18']
# years = ['17']

for i in range(3,5):
    for year in years:
        if year=='16pre':
            ppsf=ppsf_16pre_highmY
            ddsf=ddsf_16pre_highmY
        elif year=='16post':
            ppsf=ppsf_16post_highmY
            ddsf=ddsf_16post_highmY
        elif year=='17':
            ppsf=ppsf_17_highmY
            ddsf=ddsf_17_highmY
        elif year=='18':
            ppsf=ppsf_18_highmY
            ddsf=ddsf_18_highmY
        # ppfile="/eos/cms/store/group/phys_b2g/shsong/HiggsDNA/PP"+year+"/DiPhotonJetsBox_cat"+str(i)+".parquet"
        # ddfile="/eos/user/s/shsong/HiggsDNA/dd/DatadrivenQCD_20"+year+"_cat"+str(i)+".parquet"
        ppfile="/eos/cms/store/group/phys_b2g/shsong/HiggsDNA/PP"+year+"_highmY/DiPhotonJetsBox_cat"+str(i)+".parquet"
        ddfile="/eos/user/s/shsong/HiggsDNA/dd/DatadrivenQCD_20"+year+"_cat"+str(i)+"_highmY.parquet"
        pp=load_parquet(ppfile)
        dd=load_parquet(ddfile)
        pp_reweighted=add_sale_factor(pp,ppsf[i-1])
        dd_reweighted=add_sale_factor(dd,ddsf[i-1])
        # ak.to_parquet(pp_reweighted,"/eos/user/s/shsong/HiggsDNA/dd/DiPhotonJetsBox_20"+year+"_cat"+str(i)+"_reweighted.parquet")
        # ak.to_parquet(dd_reweighted,"/eos/user/s/shsong/HiggsDNA/dd/DatadrivenQCD_20"+year+"_cat"+str(i)+"_reweighted.parquet")
        ak.to_parquet(pp_reweighted,"/eos/user/s/shsong/HiggsDNA/dd/DiPhotonJetsBox_20"+year+"_cat"+str(i)+"_highmY_reweighted.parquet")
        ak.to_parquet(dd_reweighted,"/eos/user/s/shsong/HiggsDNA/dd/DatadrivenQCD_20"+year+"_cat"+str(i)+"_highmY_reweighted.parquet")
    for year in years:
        if year=='16pre':
            ppsf=ppsf_16pre
            ddsf=ddsf_16pre
        elif year=='16post':
            ppsf=ppsf_16post
            ddsf=ddsf_16post
        elif year=='17':
            ppsf=ppsf_17
            ddsf=ddsf_17
        elif year=='18':
            ppsf=ppsf_18
            ddsf=ddsf_18
        # ppfile="/eos/cms/store/group/phys_b2g/shsong/HiggsDNA/PP"+year+"/DiPhotonJetsBox_cat"+str(i)+".parquet"
        # ddfile="/eos/user/s/shsong/HiggsDNA/dd/DatadrivenQCD_20"+year+"_cat"+str(i)+".parquet"
        ppfile="/eos/cms/store/group/phys_b2g/shsong/HiggsDNA/PP"+year+"/DiPhotonJetsBox_cat"+str(i)+".parquet"
        ddfile="/eos/user/s/shsong/HiggsDNA/dd/DatadrivenQCD_20"+year+"_cat"+str(i)+".parquet"
        pp=load_parquet(ppfile)
        dd=load_parquet(ddfile)
        pp_reweighted=add_sale_factor(pp,ppsf[i-1])
        dd_reweighted=add_sale_factor(dd,ddsf[i-1])
        # ak.to_parquet(pp_reweighted,"/eos/user/s/shsong/HiggsDNA/dd/DiPhotonJetsBox_20"+year+"_cat"+str(i)+"_reweighted.parquet")
        # ak.to_parquet(dd_reweighted,"/eos/user/s/shsong/HiggsDNA/dd/DatadrivenQCD_20"+year+"_cat"+str(i)+"_reweighted.parquet")
        ak.to_parquet(pp_reweighted,"/eos/user/s/shsong/HiggsDNA/dd/DiPhotonJetsBox_20"+year+"_cat"+str(i)+"_reweighted.parquet")
        ak.to_parquet(dd_reweighted,"/eos/user/s/shsong/HiggsDNA/dd/DatadrivenQCD_20"+year+"_cat"+str(i)+"_reweighted.parquet")
