#!/usr/bin/env python
# coding: utf-8

# In[3]:


from obspy import UTCDateTime
from obspy.clients.fdsn import Client
import matplotlib.pyplot as plt
from obspy.geodetics.base import gps2dist_azimuth,kilometer2degrees
from obspy.taup import TauPyModel
from obspy.signal.rotate import rotate_ne_rt
from matplotlib.ticker import MultipleLocator
from obspy import read, read_inventory, Stream, Trace
import numpy as np
from obspy.io.sac.sacpz import attach_paz
from obspy.signal.filter import envelope
import pandas as pd
import os,glob
import warnings
warnings.filterwarnings("ignore")
model = TauPyModel(model="iasp91")
client = Client("IRIS")
plt.rcParams['font.sans-serif']='Times New Roman'


# In[4]:


DATA_DIR = '/Volumes/home/Research/DataBase/00_'
PZ_DIR  =  '/Volumes/home/Research/DataBase/00_PZs'
INFO_DIR = '/Volumes/home/Research/DataBase/Armenia'
freqmin = 0.04
freqmax = 0.125
FIG_DIR = f'/Volumes/home/Research/Progress/01_AMTG_record_plot_SKS_{freqmin}-{freqmax}'
if not os.path.isdir(FIG_DIR):
    os.mkdir(FIG_DIR)

Rphase = ['S','Sdiff', 'SS' ,'SSS' ,'SSSS' ,'SKS' ,'SKKS' ,'SKKKS' ,'ScP' ,'SP' ,'PS' ,'PcS','SKP','PKS' ]
Tphase = ['S','Sdiff','SS','SSS','SSSS','ScS']
# phaselist = set( Rphase + Tphase)
phaselist = ['SKS','SKKS']
phasecolor = {'SKS':'lightcoral', 'SKKS':'lightblue','SKKKS':'lightgreen'}
network= ['AM','TG']


exg = 2
arr_size=20
mmm = 10**-4


# In[5]:


def SNRwindow(arr_t):
    signalbegin = arr_t -5
    signalend   = arr_t +25
    noiseend    = arr_t -20 
    noisebegin  = arr_t -20-60
    return signalbegin,signalend,noiseend,noisebegin
def checkday(num):
    if len(str(num)) == 1 :
        num = str(0)+str(num)
    return str(num)


# In[6]:


starttime = UTCDateTime("2010-10-01")
endtime = UTCDateTime("2015-10-01")
cat = client.get_events(starttime=starttime, endtime=endtime,
                        minmagnitude=6,latitude =41.115,longitude=43.8036,
                        minradius=85,maxradius=140)
# print(cat.__str__(print_all=True))


# In[16]:


Earthquake = []
MAGnitude=[]
NETWORK = []
STATION = []
DEGLIST = []
SKS_SNR_LIST = []
SKKS_SNR_LIST = []


for cata in cat[0:1]:
    eq_time = cata.origins[0].time
    print(eq_time)
    eq_lon = cata.origins[0].longitude
    eq_lat = cata.origins[0].latitude
    depth  = cata.origins[0].depth/1000
    mag    = cata.magnitudes[0].mag
    mag_type = cata.magnitudes[0].magnitude_type
    
    inv = client.get_stations(network="GO", station='*', channel="*",
                                    starttime=eq_time,endtime=eq_time+30*60)
    inv1 = client.get_stations(network="IU", station='GNI', channel="*",
                                    starttime=eq_time,endtime=eq_time+30*60)
    inv2 = client.get_stations(network="II", station='KIV', channel="*",
                                    starttime=eq_time,endtime=eq_time+30*60)
    inv = inv + inv1 + inv2   

##========================================================
    i = 0 
    for net_i in inv:
        NET=net_i.code
        if NET == 'GO':
            preflit = [0.001,0.005,45,50]
        elif NET == 'II' or NET == 'IU':
            preflit = [0.001,0.005,9,10]
        for sta_i in net_i:
            STA = sta_i.code
            st_lat = sta_i.latitude
            st_lon = sta_i.longitude
            dist,azi,baz = gps2dist_azimuth(eq_lat,eq_lon,st_lat,st_lon)
            dist = dist/1000
            deg = kilometer2degrees(dist)
            
            arrivals = model.get_travel_times(source_depth_in_km=depth,distance_in_degree=deg,phase_list=phaselist)
            arr_SKS = arrivals[0].time
            arr_SKKS = arrivals[1].time            
            
            i+=1
            try:
                ori_st = client.get_waveforms(NET, STA,'*','*',starttime=eq_time+arr_SKS-120, endtime=eq_time+arr_SKKS+100,attach_response=True)
                st = ori_st.copy()
                st.remove_response(pre_filt = preflit, output="DISP")
                st.detrend('linear')
                st.detrend('demean')
                st.taper(0.05,type='cosine')
                st.filter('bandpass',freqmin=freqmin,freqmax=freqmax,corners=4,zerophase=True)
                dt = 1 / st[0].stats.sampling_rate
                if NET == 'GO':
                    HHE = st.select(component='E')[0].data
                    HHN = st.select(component='N')[0].data
                    HHR,HHT = rotate_ne_rt(HHN,HHE,baz)
                elif NET == 'IU' or NET == 'II' :
                    HH2 = st.select(location='00', channel='BH2')[0].data
                    HH1 = st.select(location='00', channel='BH1')[0].data
                    HHR,HHT = rotate_ne_rt(HH1,HH2,baz)

                Earthquake.append(eq_time)
                MAGnitude.append(mag)
                NETWORK.append(NET)
                STATION.append(STA)
                DEGLIST.append(deg)
    ##==========================calculate SNR =========================       
                SKS_signalbegin,SKS_signalend,SKS_noiseend,SKS_noisebegin = SNRwindow(arr_SKS)
                temp_tr = Trace(data=HHR)
                temp_tr.stats.delta = dt
                temp_tr.stats.starttime = st[0].stats.starttime

                SKS_noise = temp_tr.slice(starttime=eq_time+SKS_noisebegin,endtime=eq_time+SKS_noiseend)
                SKS_signal = temp_tr.slice(starttime=eq_time+SKS_signalbegin ,endtime = eq_time+SKS_signalend)
                SKS_signal_envelope = envelope(SKS_signal.data)
                SKS_noise_envelope = envelope(SKS_noise.data)

                SKS_cal_signal = sum(SKS_signal_envelope**2)
                SKS_cal_noise = sum(SKS_noise_envelope**2)

                SKS_SNR = int(SKS_cal_signal * 2 / SKS_cal_noise)
                if SKS_signalbegin < arr_SKKS < SKS_signalend: 
#                     print('SKKS nan')
                    SKKS_SNR = np.nan
                else : 
                    SKKS_signalbegin,SKKS_signalend,SKKS_noiseend,SKKS_noisebegin = SNRwindow(arr_SKKS)
                    SKKS_signal = temp_tr.slice(starttime=eq_time+SKKS_signalbegin ,endtime = eq_time+SKKS_signalend)
                    SKKS_signal_envelope = envelope(SKKS_signal.data)
                    SKKS_cal_signal = sum(SKKS_signal_envelope**2)
                    SKKS_SNR = int(SKKS_cal_signal * 2 / SKS_cal_noise)
                SKS_SNR_LIST.append(SKS_SNR)
                SKKS_SNR_LIST.append(SKKS_SNR)
            except: pass
#                 print(f'{NET} {STA} byebye') 
#                 SKS_SNR = np.nan
#                 SKKS_SNR = np.nan


######==============變數刪除================    
        del locals()['arr_SKS']
        del locals()['arr_SKKS']  

    
SNRdf  = pd.DataFrame({'UTCDateTime':Earthquake,
                       'Magnitude':MAGnitude,
                      'Network':NETWORK,
                      'Station':STATION,
                      'Dist':DEGLIST,
                      'SKS_SNR':SKS_SNR_LIST,
                      'SKKS_SNR':SKKS_SNR_LIST})  


SNRdf.to_csv(f'/Volumes/home/Research/Progress/01_GO_phase_SNR_{freqmin}-{freqmax}.csv',index=False)


# In[19]:





# In[ ]:




