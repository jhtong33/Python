#!/usr/bin/env python
# coding: utf-8


from obspy import UTCDateTime
from obspy.clients.fdsn import Client
import matplotlib.pyplot as plt
from obspy.geodetics.base import gps2dist_azimuth,kilometer2degrees
from obspy.taup import TauPyModel
from obspy.signal.rotate import rotate_ne_rt
from matplotlib.ticker import MultipleLocator
from obspy.signal.filter import envelope
from obspy import read, read_inventory, Trace
import numpy as np
import os
import warnings
warnings.filterwarnings("ignore")
model = TauPyModel(model="iasp91")
client = Client("IRIS")
plt.rcParams['font.sans-serif']='Times New Roman'



# 2021-04-29T06:50:26.915000Z
starttime = UTCDateTime("2015-10-01")
# endtime = UTCDateTime("2021-06-30")
endtime = UTCDateTime("2021-04-30")
cat = client.get_events(starttime=starttime, endtime=endtime,
                        minmagnitude=6,latitude =41.115,longitude=43.8036,
                        minradius=85,maxradius=140)#,mindepth=100)
print(cat)



freqmin=0.04
freqmax=0.125
FIG_DIR=f'/Volumes/home/Research/Progress/01_GO_record_plot_SKS_{freqmin}-{freqmax}'
if not os.path.isdir(FIG_DIR):
    os.mkdir(FIG_DIR)
    
    
Rphase = ['S','Sdiff', 'SS' ,'SSS' ,'SSSS' ,'SKS' ,'SKKS' ,'SKKKS' ,'ScP' ,'SP' ,'PS' ,'PcS','SKP','PKS' ]
Tphase = ['S','Sdiff','SS','SSS','SSSS','ScS']
phaselist = set(Rphase + Tphase)
phasecolor = {'SKS':'lightcoral', 'SKKS':'lightblue','SKKKS':'lightgreen'}
EVILSTA=['SEAG','KZRT','BATM','BGD','DGRG','TRLG','CHVG','LGD','DDFL']
GOODSTA=['ONI','AKH']

mmm = 10**-4
exg = 2
arr_size=25



def SNRwindow(arr_t):
    signalbegin = arr_t -5
    signalend   = arr_t +25
    noiseend    = arr_t -20 
    noisebegin  = arr_t -20-60
    return signalbegin,signalend,noiseend,noisebegin



for cata in cat:
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

    
    deglist = [] ; STAlist = []; BAZlist = []
    for net_i in inv:
        NET=net_i.code   
        for sta_i in net_i:
            STA = sta_i.code
            if STA not in EVILSTA:
                st_lat = sta_i.latitude
                st_lon = sta_i.longitude

                dist,azi,baz = gps2dist_azimuth(eq_lat,eq_lon,st_lat,st_lon)
                dist = dist/1000
                deg = kilometer2degrees(dist)

                deglist.append(deg)
                STAlist.append(STA)
                BAZlist.append(baz)
##=============calculate traval time=====================
    min_arrivals = model.get_travel_times(source_depth_in_km=depth,distance_in_degree=min(deglist),phase_list=phaselist)
    for arr in min_arrivals:
        name=arr.phase.name
        if arr.time <1830:
            locals()['min_arr_%s' % (name)]= arr.time
            
    max_arrivals = model.get_travel_times(source_depth_in_km=depth,distance_in_degree=max(deglist),phase_list=phaselist)
    for arr in max_arrivals:
        name=arr.phase.name
        if arr.time <1830:       
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
            if min(deglist) <=140:
                if  min_arr_SKS-20 < min_arr_t < max_arr_SKKKS+20:
                    ax1.plot([min_arr_t,max_arr_t],[min(deglist)-0.1,max(deglist)+0.1],ls='--',color=ph_color) 
                    ax1.text(min_arr_t,max(deglist)+0.15,phase,c='k',fontsize=8,ma='left',
                             bbox=dict(boxstyle="round",ec=ph_color,fc='white'))
            else: 
                if min_arr_SKKS-20 < min_arr_t < max_arr_SKKS+20:
                    ax1.plot([min_arr_t,max_arr_t],[min(deglist)-0.1,max(deglist)+0.1],ls='--',color=ph_color) 
                    ax1.text(min_arr_t,max(deglist)+0.15,phase,c='k',fontsize=8,ma='left',
                             bbox=dict(boxstyle="round",ec=ph_color,fc='white'))                    
        except KeyError: continue
    for phase in Tphase:
        try:
            min_arr_t= locals()['min_arr_%s' % (phase)]
            max_arr_t= locals()['max_arr_%s' % (phase)]
            if min(deglist) <=140:
                if min_arr_SKS -20 <min_arr_t < max_arr_SKKKS+20:
                    ax2.plot([min_arr_t,max_arr_t],[min(deglist)-0.1,max(deglist)+0.1],ls='--',alpha=0.4) 
                    ax2.text(min_arr_t,max(deglist)+0.15,phase,c='k',fontsize=8,ma='left',
                             bbox=dict(boxstyle="round",ec=(1., 0.5, 0.5),fc='white'))
            else: 
                if min_arr_SKKS-20 < min_arr_t < max_arr_SKKS+20:
                    ax2.plot([min_arr_t,max_arr_t],[min(deglist)-0.1,max(deglist)+0.1],ls='--',color=ph_color) 
                    ax2.text(min_arr_t,max(deglist)+0.15,phase,c='k',fontsize=8,ma='left',
                             bbox=dict(boxstyle="round",ec=ph_color,fc='white'))                 
        except KeyError: continue
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
            if STA in STAlist:
                print(i,STA)
                baz = BAZlist[i]
                deg = deglist[i]
                i+=1
                try:
                    ori_st = client.get_waveforms(NET, STA,'*','*',eq_time+min_arr_SKS-100,eq_time+max_arr_SKKKS+25,attach_response=True)
                    st = ori_st.copy()
                    st.merge(fill_value=0)
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

                    if STA in GOODSTA:
                        color = 'k'
                    else : color='k'

                    
                    ax1.set_title('R component')
                    ax1.set_xlim(times[0],times[-1])
                    ax1.set_ylim(min(deglist)-0.2,max(deglist)+0.2)
                    ax1.set_xlabel('Time after origin time (s)',fontsize=12)
                    ax1.set_ylabel('Distence (degree)',fontsize=12)
                    ax1.text(SKS_noisebegin,deg,'|',color='grey',fontsize=15,fontweight='bold')
                    ax1.text(SKS_noiseend,deg,'|',color='grey',fontsize=15,fontweight='bold')
                    ax1.text(times[-1]+3,deg,f'{SKS_SNR},{SKKS_SNR}',fontsize=11,color=color)
                    ax1.plot(times, plot_HHR*exg+deg,c=color,lw=1.5)
                    
                    
                    ax2.set_title('T component')
                    ax2.set_xlim(times[0],times[-1])
                    ax2.set_ylim(min(deglist)-0.2,max(deglist)+0.2)
                    ax2.set_xlabel('Time after origin time (s)',fontsize=12)
                    ax2.text(times[-1]+5,deg,f'{STA}',fontsize=12,color=color)
                    ax2.plot(times, plot_HHT*exg+deg, c=color,lw=1.5)
                    plt.suptitle(f'{eq_time}\n lat:{eq_lat} lon:{eq_lon} dep:{depth}km  {mag}{mag_type}\nbp:{freqmin}-{freqmax}Hz', 
                             fontsize=20)
                except Exception as e  : print(f'{NET} {STA} no data'); continue
    ax1.text(times[0]+5,max(deglist)+0.15,f'Amp: {mmm}',fontsize=8)
    ax2.text(times[0]+5,max(deglist)+0.15,f'Amp: {mmm}',fontsize=8)
    plt.savefig(f'{FIG_DIR}/{eq_time}.png',dpi=200,facecolor='white')
    plt.close()

######==============變數刪除================    
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

