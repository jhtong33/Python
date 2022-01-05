#!/usr/bin/env python
# coding: utf-8

# In[1]:


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


# In[8]:


DATA_DIR = '/Volumes/home-2/Research/DataBase/00_'
PZ_DIR  =  '/Volumes/home-2/Research/DataBase/00_PZs'
INFO_DIR = '/Volumes/home-2/Research/DataBase/Armenia'
freqmin = 0.04
freqmax = 0.125
FIG_DIR = f'/Volumes/home-2/Research/Progress/01_AMTG_record_plot_SKS_{freqmin}-{freqmax}'
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

# Badstation = ['DDFL','DGRG','LICH','LGD','NAVR','BATM','CANZ','BAUR','GANZ','BKRG']


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


t1 = UTCDateTime("2015-10-01T00:00:00")
t2 = UTCDateTime("2020-12-31T23:59:59")
Cata= client.get_events(starttime=t1, endtime=t2, minmagnitude=6,latitude =41.115,longitude=43.8036,
                        minradius=80,maxradius=140,orderby='time-asc')


# In[7]:


df = pd.read_csv(INFO_DIR+'/Station_info.csv')


# In[21]:


Earthquake = []
MAGnitude=[]
NETWORK = []
STATION = []
DEGLIST = []
SKS_SNR_LIST = []
SKKS_SNR_LIST = []

for cata in Cata:
    eq_time = cata.origins[0].time
    print(eq_time)
    eq_lon = cata.origins[0].longitude
    eq_lat = cata.origins[0].latitude
    depth  = cata.origins[0].depth/1000
    mag    = cata.magnitudes[0].mag
    mag_type = cata.magnitudes[0].magnitude_type
    
    yyyy = eq_time.year
    mm = checkday(eq_time.month)
    dd = checkday(eq_time.day)
    hh = checkday(eq_time.hour)
    minn = checkday(eq_time.minute)
    

##========================================================
    i = 0 
    for net in network:
        NET_DIR = f'{DATA_DIR}{net}'
        NET_PZs = f'{PZ_DIR}/{net}'
        eq_DIR =  f'{NET_DIR}/{yyyy}{mm}{dd}{hh}{minn}'

        for path in sorted(glob.glob(f'{eq_DIR}/*Z')):
            
            STA = path.rsplit('.',2)[1]
            # print(STA)
            try : 
                st_lat = (df['lat'][df['station'] == STA]).values[0]
                st_lon = (df['lon'][df['station'] == STA]).values[0]
            except: 
                st_lat = (df['lat'][df['station'] == STA][0]).item()
                st_lon = (df['lon'][df['station'] == STA][0]).item()                

            dist,azi,baz = gps2dist_azimuth(eq_lat,eq_lon,st_lat,st_lon)
            dist = dist/1000
            deg = kilometer2degrees(dist)
            
            arrivals = model.get_travel_times(source_depth_in_km=depth,distance_in_degree=deg,phase_list=phaselist)
            arr_SKS = arrivals[0].time
            arr_SKKS = arrivals[1].time
            
            Earthquake.append(eq_time)
            MAGnitude.append(mag)
            NETWORK.append(net)
            STATION.append(STA)
            DEGLIST.append(deg)
            i+=1
            ori_st = Stream()
            for datapath in glob.glob(f'{eq_DIR}/*{STA}.HH?'):
                channel = datapath.rsplit('.',1)[-1]
                tr4pz = Trace()
                PZs = glob.glob(f'{NET_PZs}/{STA}/*{STA}_{channel}.txt')
                if PZs == [] and net == 'TG' :
                    PZs = glob.glob(f'{NET_PZs}/ABST/*ABST_{channel}.txt')
                elif PZs == [] and net == 'AM':
                    PZs = glob.glob(f'{NET_PZs}/ARZA/*ARZA_{channel}.txt')
                attach_paz(tr4pz,PZs[0])
                paz = dict(tr4pz.stats.paz)
                tr = read(datapath,starttime=eq_time+arr_SKS-120, endtime=eq_time+arr_SKKS+100)
                tr.simulate(paz_remove=paz,pre_filt=(0.033, 0.034, 45, 50))
                ori_st += tr
                
            try: 
                st = ori_st.copy()
                st.merge(fill_value=0)
                st.detrend('linear')
                st.detrend('demean')
                st.taper(0.05,type='cosine')
                st.filter('bandpass',freqmin=freqmin,freqmax=freqmax,corners=4,zerophase=True)
                dt = 1 / st[0].stats.sampling_rate
                HHE = st.select(component='E')[0].data
                HHN = st.select(component='N')[0].data
                HHZ = st.select(component='Z')[0].data
                HHR,HHT = rotate_ne_rt(HHN,HHE,baz)

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
                    # print('SKKS nan')
                    SKKS_SNR = np.nan
                else : 
                    SKKS_signalbegin,SKKS_signalend,SKKS_noiseend,SKKS_noisebegin = SNRwindow(arr_SKKS)
                    SKKS_signal = temp_tr.slice(starttime=eq_time+SKKS_signalbegin ,endtime = eq_time+SKKS_signalend)
                    SKKS_signal_envelope = envelope(SKKS_signal.data)
                    SKKS_cal_signal = sum(SKKS_signal_envelope**2)
                    SKKS_SNR = int(SKKS_cal_signal * 2 / SKS_cal_noise)
            except: 
                SKS_SNR = np.nan
                SKKS_SNR = np.nan
                # print(f'{net} {STA} byebye')
            SKS_SNR_LIST.append(SKS_SNR)
            SKKS_SNR_LIST.append(SKKS_SNR)
###==============變數刪除================
    del locals()['arr_SKS']
    del locals()['arr_SKKS']      
    



SNRdf  = pd.DataFrame({'UTCDateTime':Earthquake,
                       'Magnitude':MAGnitude,
                      'Network':NETWORK,
                      'Station':STATION,
                      'Dist':DEGLIST,
                      'SKS_SNR':SKS_SNR_LIST,
                      'SKKS_SNR':SKKS_SNR_LIST})  



SNRdf.to_csv(f'/Volumes/home-2/Research/Progress/01_AMTG_phase_SNR_{freqmin}-{freqmax}_2015-2021.csv',index=False)






