import math
import pandas as pd
import numpy as np
from tanks_helper import *

df_chem = pd.read_csv('chemical_db.csv')
df_met = pd.read_csv('met_db.csv')

INPUT_CITY = 'Denver, Colorado'                     # From User Data
INPUT_TANK = [12, 8, 6]                             # tkshellht, skliqht, diameter
INPUT_CONTENTS = [8450, 'other stocks', 11.5, 4.5]  # throughput, productfactor, hlx, hln
CHEM_LIST = ['Cyclohexane', 'Benzene', 'Toluene']   # From User Data
ANNUAL_QUANTITY = [101, 2812, 258]                  # From User Data

df_met_sub = df_met[['ATMOS_PRS', 'INSOL_ANN', 'CTYST', 'CITY', 'STATE', 'TAX_ANN', 'TAN_ANN']]
df_met_filter = df_met_sub[df_met_sub['CTYST'] == INPUT_CITY]

MET_LIST = df_met_filter.values.tolist()

tank = VerticalFixedRoofTank(tkshellht=INPUT_TANK[0],           # From User Data
                             skliqht=INPUT_TANK[1],             # From User Data
                             tkrfslope=0.0625,                  # Default
                             diameter=INPUT_TANK[2],            # From User Data
                             ins=1491,                          # Default
                             solarabs=0.25,                     # From User Data
                             tax=MET_LIST[0][5],                # From Met Table
                             tan=MET_LIST[0][6],                # From Met Table
                             atmplocal=MET_LIST[0][0],          # From Met Table
                             throughput=INPUT_CONTENTS[0],      # From User Data
                             productfactor=INPUT_CONTENTS[1],   # From User Data
                             hlx=INPUT_CONTENTS[2],             # From User Data
                             hln=INPUT_CONTENTS[3],             # From User Data
                             ventsetting=1)                     # Default


name_ = np.array([df_chem[df_chem['NAME'].isin(CHEM_LIST)]['NAME'].tolist()])
cas_ = np.array([df_chem[df_chem['NAME'].isin(CHEM_LIST)]['CAS'].tolist()])
mw_ = np.array([df_chem[df_chem['NAME'].isin(CHEM_LIST)]['MOLWT'].tolist()])
vp_a = np.array([df_chem[df_chem['NAME'].isin(CHEM_LIST)]['VP_COEF_A'].tolist()])
vp_b = np.array([df_chem[df_chem['NAME'].isin(CHEM_LIST)]['VP_COEF_B'].tolist()])
vp_c = np.array([df_chem[df_chem['NAME'].isin(CHEM_LIST)]['VP_COEF_C'].tolist()])
arr = np.concatenate((name_, cas_, mw_, vp_a, vp_b, vp_c), axis=0).T

df1 = pd.DataFrame(data=arr, columns=['component',
                                      'cas_no',
                                      'mw',
                                      'antoine_coef_a',
                                      'antoine_coef_b',
                                      'antoine_coef_c'])

df1['comp_amt'] = ANNUAL_QUANTITY

df1['comp_vp'] = 10**(df1['antoine_coef_a'].astype(float) -
                      ((df1['antoine_coef_b'].astype(float)) / (tank.tla_c() +
                                                  (df1['antoine_coef_c'].astype(float))))) / 51.715

df1['comp_vp_tlx'] = 10**(df1['antoine_coef_a'].astype(float) -
                          ((df1['antoine_coef_b'].astype(float)) / (tank.tlx_c() +
                                                      (df1['antoine_coef_c'].astype(float))))) / 51.715

df1['comp_vp_tln'] = 10**(df1['antoine_coef_a'].astype(float) - ((df1['antoine_coef_b'].astype(float)) / (tank.tln_c() + (df1['antoine_coef_c'].astype(float))))) / 51.715

df1['comp_mole'] = df1['comp_amt'].astype(float) / df1['mw'].astype(float)
tot_moles = np.sum(df1['comp_mole'].tolist())

df1['comp_mole_xi'] = df1['comp_mole'].astype(float) / tot_moles
df1['comp_partial'] = df1['comp_mole_xi'].astype(float) * df1['comp_vp']
vp_mixture = np.sum(df1['comp_partial'].tolist())

df1['comp_vap_mole_frac'] = df1['comp_partial'].astype(float) / vp_mixture
df1['comp_vapor_mw_xi'] = df1['mw'].astype(float) * df1['comp_vap_mole_frac'].astype(float)
vapor_mw = np.sum(df1['comp_vapor_mw_xi'].tolist())

df1['comp_vp_tln'] = df1['comp_vp_tln'].tolist()
df1['comp_vp_tlx'] = df1['comp_vp_tlx'].tolist()

df1['comp_partial_tln'] = df1['comp_mole_xi'].astype(float) * df1['comp_vp_tln'].astype(float)
df1['comp_partial_tlx'] = df1['comp_mole_xi'].astype(float) * df1['comp_vp_tlx'].astype(float)

df1['vap_mole_xi'] = df1['comp_vap_mole_frac'].astype(float) * df1['mw'].astype(float) * 100
df1['vap_wt_xi'] = df1['vap_mole_xi'] / np.sum(df1['vap_mole_xi'].tolist())

tot_vp_tln = np.sum(df1['comp_partial_tln'].tolist())
tot_vp_tlx = np.sum(df1['comp_partial_tlx'].tolist())

tv = tank.tv()                  # Calculated Field
deltv = tank.deltv()            # Calculated Field
tla = tank.tla()                # Calculated Field
delbpv = tank.bventpress()      # Calculated Field
atmp = tank.atmp()              # From Met Table
hvo = tank.hvo                  # Calculated Field
vq = tank.vq()                  # Calculated Field
kn = tank.kn()                  # Calculated Field
kp = tank.kp()                  # Calculated Field
kb = tank.kb()                  # Calculated Field
vv = tank.vv()                  # Calculated Field

calc = EmissionCalculations(mv = vapor_mw, 
                            pva = vp_mixture, 
                            tv = tv, 
                            plx = tot_vp_tlx,
                            pln = tot_vp_tln,
                            deltv = deltv, 
                            tla = tla,
                            delbpv = delbpv, 
                            atmp = atmp,
                            hvo = hvo, 
                            vq = vq, 
                            kn = kn, 
                            kp = kp,
                            kb = kb, 
                            vv = vv)

df1['stand_loss_xi'] = df1['vap_wt_xi'] * calc.standingLosses()
df1['work_loss_xi'] = df1['vap_wt_xi'] * calc.workingLosses()
df1['total_loss_xi'] = df1['vap_wt_xi'] * calc.totalLosses()

print('Total Losses from Tank: {:.4f} lbs/year'.format(calc.totalLosses()))
print('Total Working Losses from Tank: {:.4f} lbs/year'.format(calc.workingLosses()))
print('Total Standing Losses from Tank: {:.4f} lbs/year'.format(calc.standingLosses()))
