import numpy as np
from stdb import load_db 
import matplotlib.pyplot as plt
import pandas as pd
import pickle
import glob, os
import warnings  
from obspy.geodetics.base import gps2dist_azimuth
from obspy.core.event.catalog import read_events
warnings.filterwarnings("ignore")  

def REstart():
    RCcluster = {'-90-80':[],'-80-70':[],'-70-60':[],'-60-50':[],'-50-40':[],'-40-30':[],'-30-20':[],'-20-10':[],'-1000':[],
      '00-10':[],'10-20':[],'20-30':[],'30-40':[],'40-50':[],'50-60':[],'60-70':[],'70-80':[],'80-90':[]}
    SCcluster = {'-90-80':[],'-80-70':[],'-70-60':[],'-60-50':[],'-50-40':[],'-40-30':[],'-30-20':[],'-20-10':[],'-1000':[],
      '00-10':[],'10-20':[],'20-30':[],'30-40':[],'40-50':[],'50-60':[],'60-70':[],'70-80':[],'80-90':[]}
    return RCcluster,SCcluster
def Calc_rho(RCdt, SCdt):
    rho = RCdt/SCdt 
    return rho
def Calc_Phi(RCPhi, SCPhi):
    Phi = max (abs(RCPhi-SCPhi), abs(SCPhi-RCPhi))
    if Phi > 90: Phi = 180 - Phi
    return Phi
def PolarPlot(ax,title):
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    ax.set_thetagrids( range(-90,90+1, 30) )
    ax.set_rticks(range(1,6+1))
    ax.set_rmax(6)
    ax.set_rlabel_position(90)
    ax.yaxis.set_tick_params(labelsize=12)
    ax.xaxis.set_tick_params(labelsize=12)
    ax.set_title(title,fontsize=18,fontweight='bold')
    ax.set_thetamin(-90) # set the limits
    ax.set_thetamax(90)
def dtPlot(ax,dtmean=None,dterror=None):
    ax.set_ylabel('dt (s)',fontsize=15)
    ax.set_xlabel('Backazimuth',fontsize=15)
    ax.set_xlim(0,360)
    ax.set_ylim(0,4.2)
    ax.grid(axis='y',color='lightgrey',lw=1)
    ax.set_yticks(np.arange(0,4.1,0.5))
    ax.set_yticklabels(['0','0.5','1','1.5','2','2.5','3','3.5','Null'])
    ax.set_xticks(range(0,361,30))
    ax.vlines(x=[90,180,270],ymin=0,ymax=4.2,color='lightgrey',lw=1)
    if dtmean != None and dterror != None:
        ax.set_title(f'{round(dtmean,2)}$\pm${round(dterror,2)}(s)',fontsize=15)
def phiPlot(ax):
    ax.set_ylabel('Fast Direction',fontsize=15)
    ax.set_xlabel('Backazimuth',fontsize=15)
    ax.set_xlim(0,360)
    ax.set_ylim(-90,90)
    ax.grid(axis='y',color='lightgrey',lw=1)
    ax.set_xticks(range(0,361,30))
    ax.vlines(x=[90,180,270],ymin=-90,ymax=90,color='lightgrey',lw=1)
def Cluster(Phi):
    phii = Phi//10
    if phii >=0:
        start = f'{int(phii)}0'
        end = f'{int(phii+1)}0'
        strr = f'{start}-{end}'
    else: 
        start = f'{int(phii)}0'
        end = f'{int(phii+1)}0' 
        strr = f'{start}{end}'
    return strr
def angle_mean(dt, phi, ddt, dphi):
    dt=np.array(dt)
    phi=np.array(phi)
    ddt=np.array(ddt)
    dphi=np.array(dphi)
    x = dt*np.cos(2*phi*np.pi/180.0)
    y = dt*np.sin(2*phi*np.pi/180.0)
    c = x + 1j*y
    m = np.mean(c)

    phase = np.angle(m, deg=True)/2.
    radius = np.abs(m)
    dphase = np.sqrt(np.sum(dphi**2))/len(x)
    dradius = np.sqrt(np.sum(ddt**2))/len(x)

    return phase, dphase, radius, dradius

#=======================================
width = (2*np.pi) / 36
def plot_nonnull_average(NET, STA, df, stlat, stlon, SAVEpath, YearRange, ver):
    RCcluster,SCcluster = REstart()
    phiRC = []; DphiRC = []; dtRC = []; DdtRC = []
    phiSC = []; DphiSC = []; dtSC = []; DdtSC = []
    BAZ = []; 
    NULLBAZ_G = [] ; NULLdt_G = []; NULLBAZ_F = [] ; NULLdt_F = []
    plt.rc('font',family='Times New Roman')
    plt.figure(figsize=(10,8))
    ax1 = plt.subplot(221, polar=True)
    ax2 = plt.subplot(222, polar=True)
    ax3 = plt.subplot(223)
    ax4 = plt.subplot(224)
    # qqq=0
    for i in range(len(df)):
        evlat = df['Ev_lat'].values[i]
        evlon = df['Ev_lon'].values[i]
            
        RCPhi = df['RCPhi'].values[i]
        RCDPhi  = df['RCPhi_std'].values[i]
        RCdt  = df['RCdt'].values[i]
        RCDdt  = df['RCdt_std'].values[i]
        SCPhi = df['SCPhi'].values[i]
        SCDPhi  = df['SCPhi_std'].values[i]
        SCdt  = df['SCdt'].values[i]
        SCDdt  = df['SCdt_std'].values[i]

        rho = Calc_rho(RCdt, SCdt)
        phi = Calc_Phi(RCPhi, SCPhi)
        dist, az, baz = gps2dist_azimuth(evlat,evlon,stlat,stlon)
        
        if 25 < phi < 68 or df["CpH"].values[i] > 0.76 :
            NULLBAZ_G.append(baz)
            NULLdt_G.append(4)
        elif 0.8 < rho < 1.1 and phi < 8:## non-null good condition 
            print(df['Event'].values[i])
            phiRC.append(float(RCPhi));DphiRC.append(float(RCDPhi));dtRC.append(float(RCdt));DdtRC.append(float(RCDdt))
            phiSC.append(float(SCPhi));DphiSC.append(float(SCDPhi));dtSC.append(float(SCdt));DdtSC.append(float(SCDdt))
            RCstr = Cluster(RCPhi)
            SCstr = Cluster(SCPhi)
            RCcluster[RCstr].append(RCPhi)
            SCcluster[SCstr].append(SCPhi)
            BAZ.append(baz)
        elif 0.7 <= rho < 1.2 and phi <= 25: ## non-null FAIR condition 
            print(df['Event'].values[i])
            phiRC.append(float(RCPhi));DphiRC.append(float(RCDPhi));dtRC.append(float(RCdt));DdtRC.append(float(RCDdt))
            phiSC.append(float(SCPhi));DphiSC.append(float(SCDPhi));dtSC.append(float(SCdt));DdtSC.append(float(SCDdt))
            RCstr = Cluster(RCPhi)
            SCstr = Cluster(SCPhi)
            RCcluster[RCstr].append(RCPhi)
            SCcluster[SCstr].append(SCPhi)
            BAZ.append(baz)
        else: 
            NULLBAZ_G.append(baz)
            NULLdt_G.append(4)

    for key in SCcluster.keys():
        length = len(SCcluster[key])
        if length >0:
            if key != '-1000': num = key.rsplit('-',1)[0]
            else: num = '-10'
            PLOTPhi = int(num)+5 
            ax2.bar(PLOTPhi*np.pi/180, length , width=width, color='red')
    for key in RCcluster.keys():
        length = len(RCcluster[key])
        if length >0:
            if key != '-1000': num = key.rsplit('-',1)[0]
            else: num = '-10'
            PLOTPhi = int(num)+5 
            ax1.bar(PLOTPhi*np.pi/180, length , width=width, color='dodgerblue')
    #####============non-Null 
    if len(dtRC) > 0:
        meanphiRC, stdphiRC, meandtRC, stddtRC = angle_mean(dtRC, phiRC, DdtRC, DphiRC)
        meanphiSC, stdphiSC, meandtSC, stddtSC = angle_mean(dtSC, phiSC, DdtSC, DphiSC)
        ###============Fast direction        
        ax1.arrow(meanphiRC*np.pi/180,0,0,5.3, head_width=0.1,head_length=0.6,fc='dodgerblue',ec='k')
        ax1.text(meanphiRC*np.pi/180,7,f'{round(meanphiRC)}°',fontsize=13,color='dodgerblue',fontweight='bold')
        ax2.arrow(meanphiSC*np.pi/180,0,0,5.3, head_width=0.1,head_length=0.6,fc='darkred',ec='k')
        ax2.text(meanphiSC*np.pi/180,7,f'{round(meanphiSC)}°',fontsize=13,color='darkred',fontweight='bold')
        PolarPlot(ax1,title = f'RC: {round(meanphiRC)}°$\pm${round(stdphiRC)}°')
        PolarPlot(ax2,title = f'SC: {round(meanphiSC)}°$\pm${round(stdphiSC)}°')
        ###============dt
        dtPlot(ax3,meandtRC,stddtRC)
        dtPlot(ax4,meandtSC,stddtSC)
        ax3.errorbar(BAZ,dtRC,yerr=DdtRC,fmt='o',color ='dodgerblue' ,ecolor='k',capsize=3)
        plt.subplot(223)
        plt.scatter(NULLBAZ_G,NULLdt_G,marker='o',edgecolors='dodgerblue',linewidths=2,alpha=0.8,s=50,c='1')
        # plt.scatter(NULLBAZ_F,NULLdt_F,marker='d',edgecolors='dodgerblue',linewidths=2,alpha=0.4,s=50,c='1')
        ax4.errorbar(BAZ,dtSC,yerr=DdtSC,fmt='o',color ='red' ,ecolor='k',capsize=3)
        plt.subplot(224)
        plt.scatter(NULLBAZ_G,NULLdt_G,marker='o',edgecolors='red',linewidths=2,alpha=0.8,s=50,c='1')
        plt.scatter(NULLBAZ_F,NULLdt_F,marker='d',edgecolors='red',linewidths=2,alpha=0.4,s=50,c='1')
    #####============Null  
    else:
        meanphiRC, stdphiRC, meandtRC, stddtRC = np.NaN, np.NaN, np.NaN, np.NaN
        meanphiSC, stdphiSC, meandtSC, stddtSC = np.NaN, np.NaN, np.NaN, np.NaN
        dtPlot(ax3)
        dtPlot(ax4)
        plt.subplot(223)
        plt.scatter(NULLBAZ_G,NULLdt_G,marker='o',edgecolors='dodgerblue',linewidths=2,alpha=0.8,s=50,c='1')
        plt.subplot(224)
        plt.scatter(NULLBAZ_G,NULLdt_G,marker='o',edgecolors='red',linewidths=2,alpha=0.8,s=50,c='1')

        PolarPlot(ax1,title = f'RC')
        PolarPlot(ax2,title = f'SC')       
    plt.suptitle(f'{NET}.{STA}\nNon-Null: {len(phiRC)}\nNull: {len(NULLBAZ_F)+len(NULLBAZ_G)}',fontsize=18,fontweight='bold')
    noNull    = len(NULLBAZ_G)
    noNonNull = len(phiRC)
    plt.savefig(f'{SAVEpath}/{NET}.{STA}/{YearRange}/{NET}.{STA}_{YearRange}_meanphidt_v{ver}.png',facecolor='white')
    plt.savefig(f'{SAVEpath}/{NET}.{STA}/{YearRange}/{NET}.{STA}_{YearRange}_meanphidt_v{ver}.pdf')
    plt.close()
    return meanphiRC, stdphiRC, meandtRC, stddtRC, meanphiSC, stdphiSC, meandtSC, stddtSC, noNull, noNonNull
    
    
def plot_null_PCA09(NET, STA, df, stlat, stlon, SAVEpath, YearRange, meanphiSC, meandtSC, ver):
    import pygmt
    # define parameters for plotting
    pygmt.config(COLOR_BACKGROUND = 'white',
             MAP_GRID_PEN_PRIMARY = '0.3p,dimgrey',
             MAP_ANNOT_OFFSET_PRIMARY = '5p',
             MAP_ANNOT_OFFSET_SECONDARY = '5p', 
             MAP_ANNOT_OBLIQUE = '30',
             FONT_ANNOT_PRIMARY = '8p,Times-Roman', 
             FONT_LABEL = '8p',
             MAP_FRAME_WIDTH = '2p',
             MAP_FRAME_PEN = '1.2p',
             MAP_FRAME_TYPE = 'fancy',
             MAP_TICK_LENGTH_PRIMARY = '12p',
             MAP_LABEL_OFFSET = '5.5p',
             FORMAT_GEO_MAP = 'F')
    BAZnumber={'0-90':[],'90-180':[],'180-270':[],'270-360':[]}
    def Classify(BAZnumber,baz):
        if 0<=baz<90:      BAZnumber['0-90'].append(baz)
        elif 90<=baz<180:  BAZnumber['90-180'].append(baz)
        elif 180<=baz<270: BAZnumber['180-270'].append(baz)
        elif 270<=baz<360: BAZnumber['270-360'].append(baz)
    def BeginPygmt(fig):
        colall = '210/210/210'
        # first plot only basics
        fig.coast(region='g', 
              projection = f'E44/42/180/5c', 
              resolution = 'c', 
              land = colall, 
              shorelines ='1/0.1p,' + colall, 
              C = colall, 
              frame =True)
        fig.plot(data = 'PB2002_boundaries.dig.txt', pen = '0.35p,245.76/204.8/204.8')
        distlims = [90,120,150]
        for dists in distlims:
            fig.plot(x = 44, y = 42, style ='E-' + str(2 * dists) + 'd', 
                 pen ='0.3p,black,3_1:0p',t = '60')

        fig.text(x = 44, y = -48,  text = '90@.', font='4p')
        fig.text(x = 44, y = -79,  text = '120@.', font='4p')
        fig.text(x = -137, y = -71,   text = '150@.', font='4p')
    fig = pygmt.Figure()
    BeginPygmt(fig)    
    ppp = 0 
    
    for i in range(len(df)):
        evlat = df['Ev_lat'].values[i]
        evlon = df['Ev_lon'].values[i]
            
        dist, az, baz = gps2dist_azimuth(evlat,evlon,stlat,stlon)
            
        if df["CpH"].values[i] >= 0.9: 
            ppp += 1
            Classify(BAZnumber,baz)
            fig.plot(x=[stlon,evlon],y=[stlat,evlat],pen = '1,black@80')
    SCdata = [[stlon, stlat, meanphiSC-90, 300, meandtSC*5000]]
    fig.plot(data=SCdata, style="J", color="red1", pen="0.1p,black") 
     
    fig.plot(x = stlon, y = stlat,style ='c0.1c',pen = '0.3p,black', color = 'black',transparency = '5')
    fig.text(x = -135, y = -25, text='CpH>= 0.9', font = '5p,Times-Roman',justify='CM')
    fig.text(x = stlon, y = stlat-20, text= f'{NET}.{STA} ({ppp})', font = '5p,Times-Roman')
    textt = f"[{len(BAZnumber['0-90'])}/ {len(BAZnumber['90-180'])}/ {len(BAZnumber['180-270'])}/ {len(BAZnumber['270-360'])}]"
    fig.text(x = stlon, y = stlat-35, text= textt, font = '4p,Times-Roman,50')
    # fig.show()
    fig.savefig(f'{SAVEpath}/{NET}.{STA}/{YearRange}/{NET}.{STA}_{YearRange}_PCA09map_v{ver}.png')
    fig.savefig(f'{SAVEpath}/{NET}.{STA}/{YearRange}/{NET}.{STA}_{YearRange}_PCA09map_v{ver}.pdf')


def plot_null_PCA_angle_small(NET, STA, df, stlat, stlon, SAVEpath, YearRange, ver):
    RCcluster,SCcluster = REstart()
    phiRC = []; DphiRC = []; dtRC = []; DdtRC = []
    phiSC = []; DphiSC = []; dtSC = []; DdtSC = []
    BAZ = []; 
    NULLBAZ_G = [] ; NULLdt_G = []; NULLBAZ_F = [] ; NULLdt_F = []
    plt.rc('font',family='Times New Roman')
    plt.figure(figsize=(10,8))
    ax1 = plt.subplot(221, polar=True)
    ax2 = plt.subplot(222, polar=True)
    ax3 = plt.subplot(223)
    ax4 = plt.subplot(224)
    
    if os.path.isfile(f'{SAVEpath}/{NET}.{STA}/{YearRange}/plot_null_PCA_angle_small_%(YearRange)s_v{ver}.log'):
        # print('Remove log file')
        os.remove(f'{SAVEpath}/{NET}.{STA}/{YearRange}/plot_null_PCA_angle_small_%(YearRange)s_v{ver}.log')
    
    for i in range(len(df)):
        evlat = df['Ev_lat'].values[i]
        evlon = df['Ev_lon'].values[i]
            
        RCPhi = df['RCPhi'].values[i]
        RCDPhi  = df['RCPhi_std'].values[i]
        RCdt  = df['RCdt'].values[i]
        RCDdt  = df['RCdt_std'].values[i]
        SCPhi = df['SCPhi'].values[i]
        SCDPhi  = df['SCPhi_std'].values[i]
        SCdt  = df['SCdt'].values[i]
        SCDdt  = df['SCdt_std'].values[i]

        rho = Calc_rho(RCdt, SCdt)
        phi = Calc_Phi(RCPhi, SCPhi)
            
        dist, az, baz = gps2dist_azimuth(evlat,evlon,stlat,stlon)
        
        if phi <= 25 and 0.9 > df["CpH"].values[i] > 0.76: 
            phiRC.append(float(RCPhi));DphiRC.append(float(RCDPhi));dtRC.append(float(RCdt));DdtRC.append(float(RCDdt))
            phiSC.append(float(SCPhi));DphiSC.append(float(SCDPhi));dtSC.append(float(SCdt));DdtSC.append(float(SCDdt))
            RCstr = Cluster(RCPhi)
            SCstr = Cluster(SCPhi)
            RCcluster[RCstr].append(RCPhi)
            SCcluster[SCstr].append(SCPhi)
            BAZ.append(baz)
            
            evtlog = df['Event'].values[i]
            cmd = "echo %(evtlog)s >> %(SAVEpath)s/%(NET)s.%(STA)s/%(YearRange)s/plot_null_PCA_angle_small_%(YearRange)s_v{ver}.log" % locals()
            os.system(cmd )            
    for key in SCcluster.keys():
        length = len(SCcluster[key])
        if length >0:
            if key != '-1000': num = key.rsplit('-',1)[0]
            else: num = '-10'
            PLOTPhi = int(num)+5 
            ax2.bar(PLOTPhi*np.pi/180, length , width=width, color='red')
    for key in RCcluster.keys():
        length = len(RCcluster[key])
        if length >0:
            if key != '-1000': num = key.rsplit('-',1)[0]
            else: num = '-10'
            PLOTPhi = int(num)+5 
            ax1.bar(PLOTPhi*np.pi/180, length , width=width, color='dodgerblue')
    #####============non-Null 
    if len(dtRC) > 0:              
        dtPlot(ax3)
        dtPlot(ax4)
        ax3.errorbar(BAZ,dtRC,yerr=DdtRC,fmt='o',color ='dodgerblue' ,ecolor='k',capsize=3)
        ax4.errorbar(BAZ,dtSC,yerr=DdtSC,fmt='o',color ='red' ,ecolor='k',capsize=3)
    #####============Null     
    PolarPlot(ax1,title = f'RC')
    PolarPlot(ax2,title = f'SC')        
    plt.suptitle(f'{NET}.{STA}\n0.9>CpH>0.76 & $\Delta\Theta <25^\degree$: {len(dtRC)}',fontsize=18,fontweight='bold')
    plt.savefig(f'{SAVEpath}/{NET}.{STA}/{YearRange}/{NET}.{STA}_{YearRange}_PCA_angle_small_v{ver}.png',facecolor='white')
    plt.savefig(f'{SAVEpath}/{NET}.{STA}/{YearRange}/{NET}.{STA}_{YearRange}_PCA_angle_small_v{ver}.pdf')
    plt.close()
def plot_null_PCA_angle_large(NET, STA, df, stlat, stlon, SAVEpath, YearRange, ver):
    RCcluster,SCcluster = REstart()
    phiRC = []; DphiRC = []; dtRC = []; DdtRC = []
    phiSC = []; DphiSC = []; dtSC = []; DdtSC = []
    BAZ = []; 
    plt.rc('font',family='Times New Roman')
    plt.figure(figsize=(10,8))
    ax1 = plt.subplot(221, polar=True)
    ax2 = plt.subplot(222, polar=True)
    ax3 = plt.subplot(223)
    ax4 = plt.subplot(224)
    
    
    for i in range(len(df)):
        evlat = df['Ev_lat'].values[i]
        evlon = df['Ev_lon'].values[i]
            
        RCPhi = df['RCPhi'].values[i]
        RCDPhi  = df['RCPhi_std'].values[i]
        RCdt  = df['RCdt'].values[i]
        RCDdt  = df['RCdt_std'].values[i]
        SCPhi = df['SCPhi'].values[i]
        SCDPhi  = df['SCPhi_std'].values[i]
        SCdt  = df['SCdt'].values[i]
        SCDdt  = df['SCdt_std'].values[i]

        rho = Calc_rho(RCdt, SCdt)
        phi = Calc_Phi(RCPhi, SCPhi)
            
        dist, az, baz = gps2dist_azimuth(evlat,evlon,stlat,stlon)
        
        if phi > 25 and 0.9 > df["CpH"].values[i] > 0.76: 
            phiRC.append(float(RCPhi));DphiRC.append(float(RCDPhi));dtRC.append(float(RCdt));DdtRC.append(float(RCDdt))
            phiSC.append(float(SCPhi));DphiSC.append(float(SCDPhi));dtSC.append(float(SCdt));DdtSC.append(float(SCDdt))
            RCstr = Cluster(RCPhi)
            SCstr = Cluster(SCPhi)
            RCcluster[RCstr].append(RCPhi)
            SCcluster[SCstr].append(SCPhi)
            BAZ.append(baz)
            
    maxlenSC = 0 ;  maxlenRC = 0                   
    for key in SCcluster.keys():
        length = len(SCcluster[key])
        if length > maxlenSC : maxlenSC = length
        if length >0:
            if key != '-1000': num = key.rsplit('-',1)[0]
            else: num = '-10'
            PLOTPhi = int(num)+5 
            ax2.bar(PLOTPhi*np.pi/180, length , width=width, color='red')
    for key in RCcluster.keys():
        length = len(RCcluster[key])
        if length > maxlenRC : maxlenRC = length
        if length >0:
            if key != '-1000': num = key.rsplit('-',1)[0]
            else: num = '-10'
            PLOTPhi = int(num)+5 
            ax1.bar(PLOTPhi*np.pi/180, length , width=width, color='dodgerblue')
    #####============non-Null 
    if len(dtRC) > 0:              
        phiPlot(ax3)
        phiPlot(ax4)
        ax3.errorbar(BAZ,phiRC,yerr=DphiRC,fmt='o',color ='dodgerblue' ,ecolor='k',capsize=3)
        ax4.errorbar(BAZ,phiSC,yerr=DphiSC,fmt='o',color ='red' ,ecolor='k',capsize=3)
    #####============Null     
    PolarPlot(ax1,title = f'RC')
    PolarPlot(ax2,title = f'SC')  
    pmax = max(maxlenRC,maxlenSC)     
    ax1.set_rticks(range(1,pmax+1,2)); ax1.set_rmax(pmax)
    ax2.set_rticks(range(1,pmax+1,2)); ax2.set_rmax(pmax) 
    plt.suptitle(f'{NET}.{STA}\n0.9>CpH>0.76 & $\Delta\Theta >25^\degree$: {len(dtRC)}',fontsize=18,fontweight='bold')
    plt.savefig(f'{SAVEpath}/{NET}.{STA}/{YearRange}/{NET}.{STA}_{YearRange}_PCA_angle_large_v{ver}.png',facecolor='white')
    plt.savefig(f'{SAVEpath}/{NET}.{STA}/{YearRange}/{NET}.{STA}_{YearRange}_PCA_angle_large_v{ver}.pdf')
    plt.close()
def plot_null_angle(NET, STA, df, stlat, stlon, SAVEpath, YearRange, ver):
    RCcluster,SCcluster = REstart()
    phiRC = []; DphiRC = []; dtRC = []; DdtRC = []
    phiSC = []; DphiSC = []; dtSC = []; DdtSC = []
    BAZ = []; 
    plt.rc('font',family='Times New Roman')
    plt.figure(figsize=(10,8))
    ax1 = plt.subplot(221, polar=True)
    ax2 = plt.subplot(222, polar=True)
    ax3 = plt.subplot(223)
    ax4 = plt.subplot(224)
    
    for i in range(len(df)):
        evlat = df['Ev_lat'].values[i]
        evlon = df['Ev_lon'].values[i]
            
        RCPhi = df['RCPhi'].values[i]
        RCDPhi  = df['RCPhi_std'].values[i]
        RCdt  = df['RCdt'].values[i]
        RCDdt  = df['RCdt_std'].values[i]
        SCPhi = df['SCPhi'].values[i]
        SCDPhi  = df['SCPhi_std'].values[i]
        SCdt  = df['SCdt'].values[i]
        SCDdt  = df['SCdt_std'].values[i]

        rho = Calc_rho(RCdt, SCdt)
        phi = Calc_Phi(RCPhi, SCPhi)
            
        dist, az, baz = gps2dist_azimuth(evlat,evlon,stlat,stlon)
        if df["CpH"].values[i]>0.76:
            pass
        elif df["CpH"].values[i]<=0.76 and phi>25: #68>phi>25:
            phiRC.append(float(RCPhi));DphiRC.append(float(RCDPhi));dtRC.append(float(RCdt));DdtRC.append(float(RCDdt))
            phiSC.append(float(SCPhi));DphiSC.append(float(SCDPhi));dtSC.append(float(SCdt));DdtSC.append(float(SCDdt))
            RCstr = Cluster(RCPhi)
            SCstr = Cluster(SCPhi)
            RCcluster[RCstr].append(RCPhi)
            SCcluster[SCstr].append(SCPhi)
            BAZ.append(baz)
        elif 0.8 < rho < 1.1 and phi < 8:## non-null good condition 
            pass
        elif 0.7 <= rho < 1.2 and phi <= 25: ## non-null FAIR condition 
            pass
        else: 
            phiRC.append(float(RCPhi));DphiRC.append(float(RCDPhi));dtRC.append(float(RCdt));DdtRC.append(float(RCDdt))
            phiSC.append(float(SCPhi));DphiSC.append(float(SCDPhi));dtSC.append(float(SCdt));DdtSC.append(float(SCDdt))
            RCstr = Cluster(RCPhi)
            SCstr = Cluster(SCPhi)
            RCcluster[RCstr].append(RCPhi)
            SCcluster[SCstr].append(SCPhi)
            BAZ.append(baz)     
            

    for key in SCcluster.keys():
        length = len(SCcluster[key])
        if length >0:
            if key != '-1000': num = key.rsplit('-',1)[0]
            else: num = '-10'
            PLOTPhi = int(num)+5 
            ax2.bar(PLOTPhi*np.pi/180, length , width=width, color='red')
    for key in RCcluster.keys():
        length = len(RCcluster[key])
        if length >0:
            if key != '-1000': num = key.rsplit('-',1)[0]
            else: num = '-10'
#             num = key.rsplit('-',1)[0]
            PLOTPhi = int(num)+5 
            ax1.bar(PLOTPhi*np.pi/180, length , width=width, color='dodgerblue')
#####============non-Null 
    if len(dtRC) > 0:              
        dtPlot(ax3)
        dtPlot(ax4)
        ax3.errorbar(BAZ,dtRC,yerr=DdtRC,fmt='o',color ='dodgerblue' ,ecolor='k',capsize=3)
        ax4.errorbar(BAZ,dtSC,yerr=DdtSC,fmt='o',color ='red' ,ecolor='k',capsize=3)
#####============Null     
        PolarPlot(ax1,title = f'RC')
        PolarPlot(ax2,title = f'SC')        
    plt.suptitle(f'{NET}.{STA}\nPCA<0.76 & $\Delta \Theta$>25°: {len(dtRC)}',fontsize=18,fontweight='bold')
    plt.savefig(f'{SAVEpath}/{NET}.{STA}/{YearRange}/{NET}.{STA}_{YearRange}_null_v{ver}.png',facecolor='white')
    plt.savefig(f'{SAVEpath}/{NET}.{STA}/{YearRange}/{NET}.{STA}_{YearRange}_null_v{ver}.pdf')
    plt.close()
    