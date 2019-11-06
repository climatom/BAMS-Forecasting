#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Compare the observations at South Col and Balcony with forecasts from the GFS
and from the ECMWF
"""

import datetime, numpy as np, pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
plt.rcParams.update({'font.size': 8})

# Function definitions
def mae(x,y):
    try:
        out=np.mean(np.abs(x.values-y.values))
    except:
        out=np.mean(np.abs(x-y))
    return out

def mse(x,y):
    
    try:
        x=x.values[:]; y = y.values
    except:
        x = x; y = y
    assert np.sum(np.isnan(x)) == 0 and np.sum(np.isnan(y)) == 0
    
    out=np.mean((x-y)**2)

    return out

def SS(obs,mod1,mod2):
    
    mse2=mse(obs,mod2); 
    mse1=mse(obs,mod1); 
    out=1-mse2/mse1
    
    # Components (see Wilks pg. 329)
    r_mod1=np.corrcoef(obs,mod1)[0,1]
    r_mod2=np.corrcoef(obs,mod2)[0,1]
    sobs=np.std(obs)
    s1=np.std(mod1)
    s2=np.std(mod2)
    ss1=[r_mod1**2,(r_mod1-s1/sobs)**2,((np.mean(mod1)-np.mean(obs))/sobs)**2]
    ss2=[r_mod2**2,(r_mod2-s2/sobs)**2,((np.mean(mod2)-np.mean(obs))/sobs)**2]
    return out, ss1, ss2
    

# Name of obs file(s)
obs_bal_f="/home/lunet/gytm3/Everest2019/AWS/Logging/balcony.csv"
obs_col_f="/home/lunet/gytm3/Everest2019/AWS/Logging/south_col.csv"

# Name of fc files
fc_bal_gfs_f=\
"/home/lunet/gytm3/Everest2019/Forecasting/NG/Logging/fclog_new_Balcony.txt"
fc_col_gfs_f=\
"/home/lunet/gytm3/Everest2019/Forecasting/NG/Logging/fclog_new_South_Col.txt"
fc_bal_ec_f=\
"/home/lunet/gytm3/Everest2019/Forecasting/ECMWF/Logging/fclog_Balcony.txt"
fc_col_ec_f=\
"/home/lunet/gytm3/Everest2019/Forecasting/ECMWF/Logging/fclog_SouthCol.txt"

# Number of hours we want to use to correct forecast.
# Also set the _hr parameters to suit the ECMWF 
k=48
min_hr=6.
max_hr=240.
step_hr=6.
force_end=datetime.datetime(year=2019,month=11,day=1)

# Read in the files
names_ec=["Valid","T","U","V","WS"]
names_gfs=["Valid","T","U","V","WS","P"]
obs_bal=pd.read_csv(obs_bal_f,parse_dates=True,index_col=0)
obs_col=pd.read_csv(obs_col_f,parse_dates=True,index_col=0)

# Fix columns so that ws is named similarly between locations
obs_col.columns=["WS_AVG_1" if x=="WS_AVG" else x for x in obs_col.columns]

fc_bal_gfs=pd.read_csv(fc_bal_gfs_f,parse_dates=True,index_col=0,delimiter="\t",\
                   names=names_gfs,skiprows=1)
fc_col_gfs=pd.read_csv(fc_col_gfs_f,parse_dates=True,index_col=0,delimiter="\t",\
                   names=names_gfs,skiprows=1)
fc_bal_ec=pd.read_csv(fc_bal_ec_f,parse_dates=True,index_col=0,delimiter="\t",\
                   names=names_ec,skiprows=1)
fc_col_ec=pd.read_csv(fc_col_ec_f,parse_dates=True,index_col=0,delimiter="\t",\
                   names=names_ec,skiprows=1)

# Bundle into list
fcst_in=[fc_bal_gfs,fc_col_gfs,fc_bal_ec,fc_col_ec]

# Get start and stop dates for the col and balcony
st_bal=obs_bal.index[0]; stop_bal=min([force_end,obs_bal.index[-1]])
st_col=obs_col.index[0]; stop_col=min([force_end,obs_col.index[-1]])

# Update forecast fields
fcst_out=[]
# Note that ECMWF starts later and finished earlier, so we will use these 
# dates to clip the gfs
st=fc_bal_ec.index[0]
end=fc_bal_ec.index[-1]
for i in fcst_in:
    i["Valid"]=pd.to_datetime(i["Valid"])
    i["dT"]=(i["Valid"]-i.index).dt.total_seconds()/3600.
    
    # Cut out values beyond end of the obs period (as they can't be checked). 
    # Also - clip to the ecmwf dates...
    fcst_out.append(i.loc[np.logical_and(np.logical_and(\
                        i["Valid"]<=np.max(obs_bal.index),i.index>=st),\
                        i.index<=end)])
        
# Now we can iterate over the dTs (in ECMWF) and compute correlations
dTs=np.unique(fcst_out[-1]["dT"]); 
#idx=np.logical_and(np.logical_and(dTs%step_hr<1,dTs>=min_hr),dTs<=max_hr);
#dTs=dTs[idx]; # These are the actual dTs  assert 1==2
dTs=np.arange(24,264,24); # Here are 24-hour dTs 
ndt=len(dTs)

# Preallocate
rt=np.zeros((ndt,5)); rw=np.zeros(rt.shape)
errt=np.zeros((rt.shape)); errw=np.zeros((rt.shape))
errw_corr=np.zeros((ndt,2)); me=np.zeros((ndt,2)) # latter is just the diff in 
skill_scores=np.zeros((ndt,2)) # MOS improvement
skill_scores_ref=np.zeros((ndt,2)) # Raw improvement (over persistence)
skill_components=np.zeros((ndt,12))
resid_pc=np.zeros((ndt,2))
names=["GFS","ECMWF"]
mod_col={names[0]:{},names[1]:{}}; 

# Ref MAEs -- persistence forecast 
ref_maes = [mae(obs_col["WS_AVG_1"].values[ii:],\
                obs_col["WS_AVG_1"].values[:-ii]) for ii in dTs]

# means
count = 0; 
for t in dTs:
    colcount=0
    for j in range(1,5,2): # note that every j is the col; [j-1] is the
        # Balcony; 1-2 is gfs; 3-4 is ecmwf
        
        #---------#
        # Balcony
        #---------#
        # subset
        sub_fc_bal=fcst_out[j-1] 
        idx=np.logical_and(sub_fc_bal["dT"]>t-0.1,sub_fc_bal["dT"]<t+0.1); #
        idx2=np.logical_and(sub_fc_bal["Valid"]>=st_bal,sub_fc_bal["Valid"]<=stop_bal)
        idx=np.logical_and(idx,idx2)
        sub_fc_bal=sub_fc_bal.loc[idx]
        sub_obs_bal=obs_bal.loc[obs_bal.index.isin(sub_fc_bal['Valid'])]
        assert (sub_obs_bal.index==sub_fc_bal["Valid"]).all(), "Dates mismatch"
        # Metrics
        rt_bal=np.corrcoef(sub_fc_bal["T"],sub_obs_bal["T_HMP"])[0,1]; 
        rw_bal=np.corrcoef(sub_fc_bal["WS"],sub_obs_bal["WS_AVG_1"])[0,1]; 
        errt_bal=mae(sub_fc_bal["T"],sub_obs_bal["T_HMP"])  
        errw_bal=mae(sub_fc_bal["WS"],sub_obs_bal["WS_AVG_1"])  
        
        #-----------#
        # South Col
        #-----------#
        # subset
        sub_fc_col=fcst_out[j] 
        idx=np.logical_and(sub_fc_col["dT"]>t-0.1,sub_fc_col["dT"]<t+0.1)
        idx2=np.logical_and(sub_fc_col["Valid"]>=st_col,sub_fc_col["Valid"]<=stop_col)
        idx=np.logical_and(idx,idx2)
        sub_fc_col=sub_fc_col.loc[idx]       
        sub_obs_col=obs_col.loc[obs_col.index.isin(sub_fc_col['Valid'])]
        assert (sub_obs_col.index==sub_fc_col["Valid"]).all(), "Dates mismatch"
        # metrics
        rt_col=np.corrcoef(sub_fc_col["T"],sub_obs_col["T_HMP"])[0,1]; 
        rw_col=np.corrcoef(sub_fc_col["WS"],sub_obs_col["WS_AVG_1"])[0,1]; 
        errt_col=mae(sub_fc_col["T"],sub_obs_col["T_HMP"])  
        errw_col=mae(sub_fc_col["WS"],sub_obs_col["WS_AVG_1"])  
        
        # Fit wind regression model and recompute the mae. Note: only on 
        # the Col
        ps=np.polyfit(sub_fc_col["WS"],sub_obs_col["WS_AVG_1"],1)
        mod=np.polyval(ps,sub_fc_col["WS"])
        # Store. Key is dT
        mod_col[names[colcount]][t]=pd.DataFrame(data=mod,\
               index=sub_fc_col["Valid"]) # <-- this is the MOS-corrected vers.
        
        err2=mae(mod,sub_obs_col["WS_AVG_1"]); 
        errw_corr[count,colcount]=err2
        me[count,colcount]=\
                    np.mean(sub_fc_col["WS"])-np.mean(sub_obs_col["WS_AVG_1"])
        _ss=SS(sub_obs_col["WS_AVG_1"],\
                    sub_fc_col["WS"],mod)
        skill_scores[count,colcount]=_ss[0]
        
        # Store jth percentile of the residual
        resid_pc[count,colcount]=np.percentile(sub_obs_col["WS_AVG_1"]-\
                mod,95)
         
        # Store
        rt[count,0]=t
        rt[count,j]=rt_bal; rt[count,j+1]=rt_col
        rw[count,0]=t
        rw[count,j]=rw_bal; rw[count,j+1]=rw_col
        errt[count,0]=t
        errt[count,j]=errt_bal; errt[count,j+1]=errt_col
        errw[count,0]=t
        errw[count,j]=errw_bal; errw[count,j+1]=errw_col
        
        colcount+=1
        
    count +=1

# Paste together errors from the col
# Format is before GFS, before ECMWF, after GFS, after ECMWF
errors=np.column_stack((errw[:,2],errw[:,4],errw_corr))
improv=(1-errors[:,2:]/errors[:,0:2])*100 # % increase 

# Compute the skill score against the persistance forecast here
skill_scores_ref=np.column_stack((1-errw[:,2]/ref_maes,1-errw[:,4]/ref_maes))

# Draw figure of Col
fig,ax=plt.subplots(2,2)
fig.set_size_inches(9,8)
# MAE
ax.flat[0].plot(dTs,errw[:,2],color="red",marker=".",markersize=20)
ax.flat[0].plot(dTs,errw[:,4],color="blue",marker=".",markersize=20)
ax.flat[0].set_xlim(22,242)
ax.flat[0].set_ylabel("MAE (m/s)")
ax.flat[0].set_xticks(dTs)
ax.flat[0].set_xticklabels([])
ax2=ax.flat[0].twinx()
ax2.scatter(dTs,skill_scores_ref[:,0],color="red",marker="*")
ax2.scatter(dTs,skill_scores_ref[:,1],color="blue",marker="*")
ax2.set_ylabel("Skill Score")

# Scatter at dT==24
gfs=fcst_out[1]; ec=fcst_out[-1]
gfs=gfs.loc[gfs["dT"]==24]; ec=ec.loc[ec["dT"]==24]
gfs=gfs.loc[np.logical_and(gfs["Valid"]>=st_col,gfs["Valid"]<=stop_col)]
ec=ec.loc[np.logical_and(ec["Valid"]>=st_col,ec["Valid"]<=stop_col)]
obs_gfs=obs_col.loc[obs_col.index.isin(gfs['Valid'])]["WS_AVG_1"]
obs_ec=obs_col.loc[obs_col.index.isin(ec['Valid'])]["WS_AVG_1"]
assert (obs_gfs.index==gfs["Valid"]).all(), "Dates mismatch"
assert (obs_ec.index==ec["Valid"]).all(), "Dates mismatch"
ps_gfs=np.polyfit(obs_gfs,gfs["WS"],1); rg=np.corrcoef(obs_gfs,gfs["WS"])[0,1]
print "GFS 24-r: %.2f"%rg
ps_ec=np.polyfit(obs_ec,ec["WS"],1); re=np.corrcoef(obs_ec,ec["WS"])[0,1]
print "ECMWF 24-r: %.2f"%re
ax.flat[1].plot(dTs,rw[:,2],color="red",marker=".",markersize=20)
ax.flat[1].plot(dTs,rw[:,4],color="blue",marker=".",markersize=20)
ax.flat[1].set_xlabel("Forecast Lead Time (hours)")
ax.flat[1].set_ylabel("Correlation")
ax.flat[1].yaxis.tick_right()
ax.flat[1].yaxis.set_label_position("right")

# Plot the MOS errors
ax.flat[2].plot(dTs,errw_corr[:,0],color="red",marker=".",markersize=20)
ax.flat[2].plot(dTs,errw_corr[:,1],color="blue",marker=".",markersize=20)
ax.flat[2].set_ylabel("MAE (m/s)")
ax.flat[2].axvline(48,linestyle="--",color="black")
ax.flat[2].set_xticks(dTs)
ax.flat[2].set_xlim(22,242)
ax.flat[2].set_xlabel("Forecast Lead Time (hours)")
#ax.flat[2].fill_between(dTs,errw_corr[:,1],resid_pc[:,1],alpha=0.2,color="blue")
#ax.flat[2].fill_between(dTs,errw_corr[:,0],resid_pc[:,0],alpha=0.2,color="red")
ax2=ax.flat[2].twinx()
ax2.scatter(dTs,skill_scores[:,0],color="red",marker="*")
ax2.scatter(dTs,skill_scores[:,1],color="blue",marker="*")
ax2.set_ylabel("Skill Score")

dates=pd.date_range(start='2019-06-07',end='2019-07-01',freq="120H")
dates_str=["%02d/%02d/%02d" % (ii.day,ii.month,ii.year) for ii in dates]
ax.flat[3].set_xticks(dates)
idx_obs=np.logical_and(obs_col.index.month == 6,obs_col.index.day>=6)
obs=obs_col.loc[idx_obs]["WS_AVG_1"]
mod_ec=mod_col["ECMWF"][k]; mod_ec = mod_ec.loc[mod_ec.index.month==6]; 
mod_ec.columns=["ECMWF"]
mod_gfs=mod_col["GFS"][k]; mod_gfs = mod_gfs.loc[mod_gfs.index.month==6]; 
mod_gfs.columns=["GFS"]
# Get indices
idx_obs=np.logical_and(obs.index.isin(mod_ec.index),obs.index.isin(mod_gfs.index))
idx_ec=np.logical_and(mod_ec.index.isin(mod_gfs.index),mod_ec.index.isin(obs.index))
idx_gfs=np.logical_and(mod_gfs.index.isin(mod_ec.index),mod_gfs.index.isin(obs.index))
# Plot
obs[idx_obs].plot(ax=ax.flat[3],color="black")
ax.flat[3].set_xticks([])
mod_ec[idx_ec].plot(ax=ax.flat[3],color="blue",label="ECMWF",linewidth=2)
ax.flat[3].set_xticks([], [])
mod_gfs.plot(ax=ax.flat[3],color="red",linewidth=2)
ax.flat[3].set_xticks([])
# Compute error distribution
err_ec=np.nanpercentile(np.abs(mod_ec["ECMWF"]-obs[idx_obs]),95)
err_gfs=np.nanpercentile(np.abs(mod_gfs["GFS"]-obs[idx_obs]),95)
print "Mu (obs) over same period = %.2f" % np.nanmean(obs[idx_obs])
print "95th %%-tile of ECMWF error = %.2f m/s" % err_ec
print "95th %%-tile of GFS error = %.2f m/s" % err_gfs
ax.flat[3].set_xticks([])
ax.flat[3].set_ylim(0,20)
ax.flat[3].yaxis.tick_right() 
ax.flat[3].set_ylabel("Wind Speed (m/s)")
ax.flat[3].yaxis.set_label_position("right")
ax.flat[3].set_xlabel("")
ax.flat[3].set_xticks([],minor=True)
ax.flat[3].set_xticks(dates)
ax.flat[3].fmt_xdata = mdates.DateFormatter('%Y-%m-%d')
ax.flat[3].set_xticklabels(dates_str,rotation=45)
ax.flat[3].fmt_xdata = mdates.DateFormatter('%Y-%m-%d')
plt.tight_layout()
fig.savefig("Forecast_Performance.png",dpi=300)



