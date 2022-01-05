#!/usr/bin/env python
# coding: utf-8


from obspy import UTCDateTime
from obspy.clients.fdsn import Client
import matplotlib.pyplot as plt
from obspy.geodetics.base import gps2dist_azimuth,kilometer2degrees
from obspy.taup import TauPyModel
from obspy.signal.rotate import rotate_ne_rt
from matplotlib.ticker import MultipleLocator
from obspy import read, read_inventory
import numpy as np
import os
import warnings
warnings.filterwarnings("ignore")
model = TauPyModel(model="iasp91")
client = Client("IRIS")




starttime = UTCDateTime("2013-01-01")
endtime = UTCDateTime("2015-05-01")
cat = client.get_events(starttime=starttime, endtime=endtime,
                        minmagnitude=6.5,latitude =41.115,longitude=43.8036,
                        minradius=100,maxradius=155,mindepth=100)


FIG_DIR = '/Volumes/home/Research/Progress/00_GO_record_plot_ZRT_0.04-0.25'
freqmin = 0.04
freqmax = 0.25
phaselist = ["Pdiff","SKS",'SKKS','Sdiff']
phasecolor= {"Pdiff":'r',
             "SKS":'b',
             "SKKS":'g',
             "Sdiff":'purple'}




for cata in cat:
    plt.figure(1,figsize=(21,10))
    
    plt.rcParams['font.sans-serif']='Times New Roman'
    eq_time = cata.origins[0].time
    print(eq_time)
    eq_lon = cata.origins[0].longitude
    eq_lat = cata.origins[0].latitude
    depth  = cata.origins[0].depth/1000
    mag    = cata.magnitudes[0].mag
    mag_type = cata.magnitudes[0].magnitude_type
    
    inv  = client.get_stations(network="GO", station='*', channel="*",
                                    starttime=eq_time,endtime=eq_time+30*60)
    inv1 = client.get_stations(network="IU", station='GNI', channel="*",
                                    starttime=eq_time,endtime=eq_time+30*60)
    inv2 = client.get_stations(network="II", station='KIV', channel="*",
                                    starttime=eq_time,endtime=eq_time+30*60)
    inv = inv + inv1 + inv2
    deglist = []
    for net_i in inv:
        NET = net_i.code
        for sta_i in net_i:
            STA = sta_i.code
            print(NET,STA)
            st_lat = sta_i.latitude
            st_lon = sta_i.longitude
            
            dist,azi,baz = gps2dist_azimuth(eq_lat,eq_lon,st_lat,st_lon)
            dist = dist/1000
            deg = kilometer2degrees(dist)
            deglist.append(deg)
            arrivals = model.get_travel_times(source_depth_in_km=depth,
                                  distance_in_degree=deg,
                                  phase_list=phaselist)
            for arr in arrivals:
                name=arr.phase.name
                locals()['arr_%s' % (name)]= arr.time
            try:
                st = client.get_waveforms(NET, STA,'*','*',eq_time,eq_time+30*60,attach_response=True)
                gap=st.get_gaps()
                EVILSTA=['SEAG','KZRT','BATM']
                if gap != []:
                    for g in gap:
                        evil_sta= g[1]
                        EVILSTA.append(evil_sta)
                st.merge(fill_value=0)
                st.remove_response(pre_filt = [0.001, 0.005, 45, 50], output="VEL")
                st.detrend('linear')
                st.detrend('demean')
                st.taper(0.05,type='cosine')
                st.filter('bandpass',freqmin=freqmin,freqmax=freqmax,corners=4,zerophase=True)
                if NET == 'GO':
                    HHE = st.select(component='E')[0].data
                    HHN = st.select(component='N')[0].data
                    HHZ = st.select(component='Z')[0].data
                    HHR,HHT = rotate_ne_rt(HHN,HHE,baz)
                elif NET == 'IU' or NET == 'II' :
                    HH2 = st.select(location='00', channel='BH2')[0].data
                    HH1 = st.select(location='00', channel='BH1')[0].data
                    HHZ = st.select(location='00', channel='BHZ')[0].data
                    HHR,HHT = rotate_ne_rt(HH1,HH2,baz)
                
                times= st[0].times()
                HHZ=HHZ/(10**-5)
                HHR=HHR/(10**-5)
                HHT=HHT/(10**-5)
            

                exg = 0.7
                arr_size=25
                if STA in EVILSTA:
                    color='lightgrey'
                else : color='k'


                ax1 = plt.subplot(131)
                ax1.set_title('Z component',fontsize=12)
                for phase in phaselist:
                    arr_t= locals()['arr_%s' % (phase)]
                    if arr_t < 1800 : 
                        pha_c= phasecolor[phase]
                        ax1.text(times[0]+arr_t,deg,'|',c=pha_c,fontsize=arr_size,alpha=0.5,fontweight='bold')
                ax1.plot(times,HHZ*exg+deg,lw=1,color=color)
                ax1.set_xlabel('Time after origin time (s)',fontsize=12)
                ax1.set_ylabel('Distence (degree)',fontsize=12)
                ax1.yaxis.set_minor_locator(MultipleLocator(0.5)) ### 0.1為間隔
                ax2 = plt.subplot(132)
                ax2.set_title('R component',fontsize=12)
                for phase in phaselist:
                    arr_t= locals()['arr_%s' % (phase)]
                    if arr_t < 1800 : 
                        pha_c= phasecolor[phase]
                        ax2.text(times[0]+arr_t,deg,'|',c=pha_c,fontsize=arr_size,alpha=0.5,fontweight='bold')
                ax2.plot(times,HHR*exg+deg,lw=1,color=color)
                ax2.set_xlabel('Time after origin time (s)',fontsize=12)
                ax2.yaxis.set_minor_locator(MultipleLocator(0.5)) ### 0.1為間隔
                ax3 = plt.subplot(133)
                ax3.set_title('T component',fontsize=12)
                for phase in phaselist:
                    arr_t= locals()['arr_%s' % (phase)]
                    if arr_t < 1800 : 
                        pha_c= phasecolor[phase]
                        ax3.text(times[0]+arr_t,deg,'|',c=pha_c,fontsize=arr_size,alpha=0.5,fontweight='bold')
                ax3.plot(times,HHT*exg+deg,lw=1,color=color)
                ax3.set_xlabel('Time after origin time (s)',fontsize=12)
                ax3.yaxis.set_minor_locator(MultipleLocator(0.5)) ### 0.1為間隔

                ax1.text(times[-1]-150,max(deglist)+0.4,'Amp: 10**-5',fontsize=8)
                ax2.text(times[-1]-150,max(deglist)+0.4,'Amp: 10**-5',fontsize=8)
                ax3.text(times[-1]-150,max(deglist)+0.4,'Amp: 10**-5',fontsize=8)
                ax1.set_ylim(min(deglist)-0.5,max(deglist)+0.5)
                ax2.set_ylim(min(deglist)-0.5,max(deglist)+0.5)
                ax3.set_ylim(min(deglist)-0.5,max(deglist)+0.5)
                ax1.set_xlim(times[0]+arr_Pdiff-100,times[-1])
                ax2.set_xlim(times[0]+arr_Pdiff-100,times[-1])
                ax3.set_xlim(times[0]+arr_Pdiff-100,times[-1])
                ax3.text(times[-1]+50,deg,STA,fontsize=12,color=color)
                
                plt.suptitle(f'{eq_time}\n lat:{eq_lat} lon:{eq_lon} dep:{depth}km  {mag}{mag_type}\nbp:{freqmin}-{freqmax}Hz', 
                         fontsize=20)

            except: 
                print(f' {NET} {STA} pass')
                continue

    plt.savefig(f'{FIG_DIR}/{eq_time}.png',dpi=200,facecolor='white')
    plt.close()

