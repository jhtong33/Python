#!/usr/bin/env python
# coding: utf-8

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
phaselist = set( Rphase + Tphase)
phasecolor = {'SKS':'lightcoral', 'SKKS':'lightblue','SKKKS':'lightgreen'}
network= ['AM','TG']


exg = 2
arr_size=20
mmm = 10**-4

# Badstation = ['DDFL','DGRG','LICH','LGD','NAVR','BATM','CANZ','BAUR','GANZ','BKRG']



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


# t1 = UTCDateTime("2015-10-01T00:00:00")
# t1 = UTCDateTime("2018-02-02T00:00:00")
t2 = UTCDateTime("2020-12-31T23:59:59")
Cata= client.get_events(starttime=t1, endtime=t2, minmagnitude=6,latitude =41.115,longitude=43.8036,
                        minradius=80,maxradius=140,orderby='time-asc')

# print(cat)




df = pd.read_csv(INFO_DIR+'/Station_info.csv')





SNRlist=[]; old_SNRlist=[]
for cata in Cata[1:]:
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
    
    deglist = [];STAlist=[]; BAZlist = []
    for net in network:
        NET_DIR = f'{DATA_DIR}{net}'
        NET_PZs = f'{PZ_DIR}/{net}'
        eq_DIR =  f'{NET_DIR}/{yyyy}{mm}{dd}{hh}{minn}'
            
        for path in sorted(glob.glob(f'{eq_DIR}/*Z')):
            # print(path)
            STA = path.rsplit('.',2)[1]
#             print(STA)

            try : 
                st_lat = (df['lat'][df['station'] == STA]).values[0]
                st_lon = (df['lon'][df['station'] == STA]).values[0]
            except: 
                st_lat = (df['lat'][df['station'] == STA][0]).item()
                st_lon = (df['lon'][df['station'] == STA][0]).item() 

            dist,azi,baz = gps2dist_azimuth(eq_lat,eq_lon,st_lat,st_lon)
            dist = dist/1000
            deg = kilometer2degrees(dist)
            # print(deg)
            STAlist.append(STA)
            deglist.append(deg)
            BAZlist.append(baz)
##=============calculate traval time=====================
    min_arrivals = model.get_travel_times(source_depth_in_km=depth,distance_in_degree=min(deglist),phase_list=phaselist)
    for arr in min_arrivals:
        name=arr.phase.name
        if arr.time <1800:
            locals()['min_arr_%s' % (name)]= arr.time
    
    max_arrivals = model.get_travel_times(source_depth_in_km=depth,distance_in_degree=max(deglist),phase_list=phaselist)
    for arr in max_arrivals:
        name=arr.phase.name
        if arr.time <1800:       
            locals()['max_arr_%s' % (name)]= arr.time
##=============plot traval time curves====================  
    plt.figure(1,figsize=(10,12))
    ax1 = plt.subplot(121)
    ax2 = plt.subplot(122)
    for phase in Rphase:
        if phase == 'SKS' or phase == 'SKKS' or phase == 'SKKKS':
            ph_color = phasecolor[phase]
        else: ph_color = 'grey'
            
        try:
            min_arr_t= locals()['min_arr_%s' % (phase)]
            max_arr_t= locals()['max_arr_%s' % (phase)]
            if  min_arr_SKS-20 < min_arr_t < max_arr_SKKS+20:
                ax1.plot([min_arr_t,max_arr_t],[min(deglist)-0.1,max(deglist)+0.1],color=ph_color,ls='--') 
                ax1.text(min_arr_t,max(deglist)+0.15,phase,c='k',fontsize=8,ma='left',
                         bbox=dict(boxstyle="round",ec=ph_color,fc='white'))
        except KeyError: continue
    for phase in Tphase:
        try:
            min_arr_t= locals()['min_arr_%s' % (phase)]
            max_arr_t= locals()['max_arr_%s' % (phase)]
            if min_arr_SKS -20 <min_arr_t < max_arr_SKKS+20:
                ax2.plot([min_arr_t,max_arr_t],[min(deglist)-0.1,max(deglist)+0.1],ls='--',alpha=0.4) 
                ax2.text(min_arr_t,max(deglist)+0.15,phase,c='k',fontsize=8,ma='left',
                         bbox=dict(boxstyle="round",ec=ph_color,fc='white'))
        except KeyError: continue
##========================================================
    i = 0 
    for net in network:
        NET_DIR = f'{DATA_DIR}{net}'
        NET_PZs = f'{PZ_DIR}/{net}'
        eq_DIR =  f'{NET_DIR}/{yyyy}{mm}{dd}{hh}{minn}'

        for path in sorted(glob.glob(f'{eq_DIR}/*Z')):
            
            STA = path.rsplit('.',2)[1]
            if STA not in STAlist: 
                continue
            else : color = 'k'
                
                
            deg = deglist[i]
            baz = BAZlist[i]
            # print(net, STA)
            i+=1
            ori_st = Stream()
            for datapath in glob.glob(f'{eq_DIR}/*{STA}.???'):
                channel = datapath.rsplit('.',1)[-1]
                tr4pz = Trace()
                PZs = glob.glob(f'{NET_PZs}/{STA}/*{STA}_{channel}.txt')
                if PZs == [] and net == 'TG' :
                    PZs = glob.glob(f'{NET_PZs}/ABST/*ABST_{channel}.txt')
                elif PZs == [] and net == 'AM':
                    PZs = glob.glob(f'{NET_PZs}/ARZA/*ARZA_{channel}.txt')
                attach_paz(tr4pz,PZs[0])
                paz = dict(tr4pz.stats.paz)
                tr = read(datapath,starttime=eq_time+min_arr_SKS-100, endtime=eq_time+max_arr_SKKS+50)
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
                arrivals = model.get_travel_times(source_depth_in_km=depth,distance_in_degree=deg,phase_list=['SKS','SKKS'])
                arr_SKS = arrivals[0].time
                arr_SKKS = arrivals[1].time

                SKS_signalbegin,SKS_signalend,SKS_noiseend,SKS_noisebegin = SNRwindow(arr_SKS)
                ax1.text(SKS_signalbegin,deg,'|',color='r',fontsize=15,fontweight='bold')
                ax1.text(SKS_signalend,deg,'|',color='r',fontsize=15,fontweight='bold')
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
                
                SKKS_signalbegin,SKKS_signalend,SKKS_noiseend,SKKS_noisebegin = SNRwindow(arr_SKKS)

                SKKS_signal = temp_tr.slice(starttime=eq_time+SKKS_signalbegin ,endtime = eq_time+SKKS_signalend)
                ax1.text(SKKS_signalbegin,deg,'|',color='r',fontsize=15,fontweight='bold')
                ax1.text(SKKS_signalend,deg,'|',color='r',fontsize=15,fontweight='bold')
                SKKS_signal_envelope = envelope(SKKS_signal.data)
                SKKS_cal_signal = sum(SKKS_signal_envelope**2)
                SKKS_SNR = int(SKKS_cal_signal * 2 / SKS_cal_noise)
        ##================================================================= 
                times = st[0].times(reftime=eq_time)
                plot_HHR = HHR/ mmm
                plot_HHT = HHT/ mmm

                ax1.plot(times, plot_HHR*exg+deg,c=color,lw=1.5)
                ax1.set_title('R component')
                ax1.set_xlim(times[0],times[-1])
                ax1.set_ylim(min(deglist)-0.2,max(deglist)+0.2)
                ax1.set_xlabel('Time after origin time (s)',fontsize=12)
                ax1.set_ylabel('Distence (degree)',fontsize=12)
                ax1.text(SKS_noisebegin,deg,'|',color='grey',fontsize=15,fontweight='bold')
                ax1.text(SKS_noiseend,deg,'|',color='grey',fontsize=15,fontweight='bold')
                ax1.text(times[-1]+3,deg,f'{SKS_SNR},{SKKS_SNR}',fontsize=11,color=color)
                ax2.plot(times, plot_HHT*exg+deg, c=color,lw=1.5)
                ax2.set_title('T component')
                ax2.set_xlim(times[0],times[-1])
                ax2.set_ylim(min(deglist)-0.2,max(deglist)+0.2)
                ax2.set_xlabel('Time after origin time (s)',fontsize=12)

                ax2.text(times[-1]+5,deg,f'{STA}',fontsize=12,color=color)
                plt.suptitle(f'{eq_time}\n lat:{eq_lat} lon:{eq_lon} dep:{depth}km  {mag}{mag_type}\nbp:{freqmin}-{freqmax}Hz', 
                         fontsize=20)
            except: print(f'{net} {STA} byebye')
######=============== for GNI ================
    try: 
        inv1 = client.get_stations(network="IU", station='GNI', channel="*",starttime=eq_time+min_arr_SKS-100, endtime=eq_time+max_arr_SKKS+25)
        print('IU GNI')
        st_lat = inv1[0][0].latitude
        st_lon = inv1[0][0].longitude
        dist,azi,baz = gps2dist_azimuth(eq_lat,eq_lon,st_lat,st_lon)
        dist = dist/1000
        deg = kilometer2degrees(dist)
        
        gni_st = client.get_waveforms('IU', "GNI",'00','*',starttime=eq_time+min_arr_SKS-100, endtime=eq_time+max_arr_SKKS+25,attach_response=True)
        gni_st.remove_response(pre_filt = [0.001,0.005,9,10], output="DISP")
        gni_st.detrend('linear')
        gni_st.detrend('demean')
        gni_st.taper(0.05,type='cosine')
        gni_st.filter('bandpass',freqmin=0.04,freqmax=0.1,corners=4,zerophase=True)
        dt = 1 / gni_st[0].stats.sampling_rate
        HH2 = gni_st.select(location='00', channel='BH2')[0].data
        HH1 = gni_st.select(location='00', channel='BH1')[0].data
        HHZ = gni_st.select(location='00', channel='BHZ')[0].data
        HHR,HHT = rotate_ne_rt(HH1,HH2,baz)
# ##==========================calculate SNR =========================       
        arrivals = model.get_travel_times(source_depth_in_km=depth,distance_in_degree=deg,phase_list=['SKS','SKKS'])
        arr_SKS = arrivals[0].time
        arr_SKKS = arrivals[1].time
        SKS_signalbegin,SKS_signalend,SKS_noiseend,SKS_noisebegin = SNRwindow(arr_SKS)
        ax1.text(SKS_signalbegin,deg,'|',color='r',fontsize=15,fontweight='bold')
        ax1.text(SKS_signalend,deg,'|',color='r',fontsize=15,fontweight='bold')
        temp_tr = Trace(data=HHR)
        temp_tr.stats.delta = dt
        temp_tr.stats.starttime = gni_st[0].stats.starttime
        SKS_noise = temp_tr.slice(starttime=eq_time+SKS_noisebegin,endtime=eq_time+SKS_noiseend)
        SKS_signal = temp_tr.slice(starttime=eq_time+SKS_signalbegin ,endtime = eq_time+SKS_signalend)
        SKS_signal_envelope = envelope(SKS_signal.data)
        SKS_noise_envelope = envelope(SKS_noise.data)

        SKS_cal_signal = sum(SKS_signal_envelope**2)
        SKS_cal_noise = sum(SKS_noise_envelope**2)

        SKS_SNR = int(SKS_cal_signal * 2 / SKS_cal_noise)
        
        SKKS_signalbegin,SKKS_signalend,SKKS_noiseend,SKKS_noisebegin = SNRwindow(arr_SKKS)
        ax1.text(SKKS_signalbegin,deg,'|',color='r',fontsize=15,fontweight='bold')
        ax1.text(SKKS_signalend,deg,'|',color='r',fontsize=15,fontweight='bold')
        SKKS_signal = temp_tr.slice(starttime=eq_time+SKKS_signalbegin ,endtime = eq_time+SKKS_signalend)
        SKKS_signal_envelope = envelope(SKKS_signal.data)
        SKKS_cal_signal = sum(SKKS_signal_envelope**2)
        SKKS_SNR = int(SKKS_cal_signal * 2 / SKS_cal_noise)
##==================================================================
        times= gni_st[0].times(reftime=eq_time)
        plot_HHR = HHR / mmm
        plot_HHT = HHT / mmm
        ax1.plot(times,plot_HHR*exg+deg,lw=1.5,color='r')
        ax1.text(SKS_noisebegin,deg,'|',color='grey',fontsize=15,fontweight='bold')
        ax1.text(SKS_noiseend,deg,'|',color='grey',fontsize=15,fontweight='bold')
        ax1.text(times[-1]+3,deg,f'{SKS_SNR},{SKKS_SNR}',fontsize=11,color=color)
        ax2.plot(times,plot_HHT*exg+deg,lw=1.5,color='r')
        ax2.text(times[-1]+5,deg,f'GNI',fontsize=12,color='r')
    except: print('GNI byebye QQ')            
    ax1.text(times[0]+5,max(deglist)+0.15,f'Amp: {mmm}',fontsize=8)
    ax2.text(times[0]+5,max(deglist)+0.15,f'Amp: {mmm}',fontsize=8)
    plt.savefig(f'{FIG_DIR}/{eq_time}.png',dpi=200,facecolor='white')
    plt.close()    
###==============變數刪除================    
    for arr in min_arrivals:
        name=arr.phase.name
        var = 'min_arr_'+name
        if var in locals():
            del locals()['min_arr_%s' % (name)]
    for arr in max_arrivals:
        name=arr.phase.name
        var = 'max_arr_'+name
        if var in locals():
            del locals()['max_arr_%s' % (name)]        
