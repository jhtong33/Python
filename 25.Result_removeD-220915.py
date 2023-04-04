import pandas as pd
import glob, os

resultpath = '/Users/tong/Documents/99_temp/Research/STEP/02_Station_result_csv'
removepath = '/Users/tong/Documents/99_temp/Research/STEP/14_result_removeD-220913/checked'
savepath   = '/Users/tong/Documents/99_temp/Research/STEP/14_result_removeD-220913'

def calc_phi_rho (SCphi, RCphi, SCdt, RCdt):
    delta = abs(RCphi-SCphi)
    if delta >90: 
        delta = 180-delta
    rho = RCdt/SCdt
    return delta, rho
    
for path in glob.glob(f'{resultpath}/*/*csv'):
    print(path)
    dirr = path.rsplit('/')[-2]
    netsta = (path.rsplit('/')[-1]).rsplit('_')[0]
    if netsta not in ['BI.MAKU', 'XG.CMCY', 'XG.DGRL','TG.BGLV']:
        removecsv = glob.glob(f'{removepath}/{netsta}_*remove_event.csv')[0]
    
        dremove = pd.read_csv(removecsv)
        remove_event = dremove['Event'].to_list()
        df = pd.read_csv(path)
    #============ just remove piercing point on D''
        removeindex = []
        for i in range(len(df)):
            Event = df['Event'].values[i]
            if Event in remove_event:
                removeindex.append(i)
        dff = df.drop(removeindex)
    #=========== clean for non-null & null
        for j in range(len(dff)):
            Event = dff['Event'].values[j]
            cph = dff['CpH'].values[j]
            pick = dff['Pick'].values[j]
            if pick == True:
                if 0.76<= cph <0.8 :
                    dff['Null'].values[j] = False
                    SCphi = dff['SCPhi'].values[j]
                    SCdt  = dff['SCdt'].values[j]
                    RCphi = dff['RCPhi'].values[j]
                    RCdt  = dff['RCdt'].values[j]
                    del_phi,rho = calc_phi_rho(SCphi, RCphi,SCdt,RCdt)
                    if del_phi >25 or rho <0.3:
                        dff['Quality'].values[j] = 'Poor'
                        cmd = '''echo %(Event)s ====== change to Poor >> %(savepath)s/%(dirr)s/%(netsta)s_log.txt''' %locals()
                        os.system(cmd)
                        print(f'{Event}======change to Poor ')
                    elif del_phi <=25 and 0.7<rho <1.2:
                        dff['Quality'].values[j] = 'Fair'
                        cmd = '''echo %(Event)s ====== change to Fair >> %(savepath)s/%(dirr)s/%(netsta)s_log.txt''' %locals()
                        os.system(cmd)
                        print(f'{Event}======change to Fair ')
                    elif del_phi <=8 or 0.8<rho <1.1: 
                        dff['Quality'].values[j] = 'Good'  
                        cmd = '''echo %(Event)s ====== change to Good >> %(savepath)s/%(dirr)s/%(netsta)s_log.txt''' %locals()
                        os.system(cmd)
                        print(f'{Event}======change to Good ')
                else:
                    pass
            else:
                pass
            
        dff.to_csv(f'{savepath}/{dirr}/{netsta}_split_result_v2.csv',index=False)
    
    