import math
import pandas as pd
import numpy as np
import pandas as pd
import numpy as np

df_chem = pd.read_csv('chemical_db.csv')
df_met = pd.read_csv('met_db.csv')


class VerticalFixedRoofTank:

    def __init__(self,
                 tkshellht,
                 skliqht,
                 tkrfslope,
                 diameter,
                 ins,
                 solarabs,
                 tax,
                 tan,
                 atmplocal,
                 throughput,
                 productfactor,
                 hlx,
                 hln,
                 ventsetting):
        '''
        Parameters
        ----------

        tkshellht : Shell height of tank, in feet.
        skliqht :
        tkrfslope :
        diameter : Shell diamter of tank, in feet.

        ins :  Average daily total insulation, defined
               to be 1,491 (btu/ft^2*d) (from Table 7.1-7)

        solarabs :  Tank surface solar absorptance,
                    user defined, depends on tank and
                    condition color.

        tax :   Average Daily Maximum Ambient
                Temperature, user defined, in deg F.

        tan :   Average Daily Minimum Ambient
                Temperature, user defined, in deg F.

        atmplocal :  Atmospheric pressure at facility.

        throughput :   Annual tank throughput, in gallons.

        productfactor :
        hlx :
        hln :
        ventsetting :
        '''
        self.tkshellht = tkshellht
        self.skliqht = skliqht
        self.tkrfslope = tkrfslope
        self.diameter = diameter
        self.ins = ins
        self.solarabs = solarabs
        self.throughput = throughput
        self.productfactor = productfactor
        self.hlx = hlx
        self.hln = hln
        self.ventsetting = 1

        # Values are in deg. F
        self.tax = tax
        self.tan = tan

        self.rad_ = (1 / 2.) * self.diameter
        self.hro_ = (1 / 3.) * self.rad_ * self.tkrfslope

        self.hvo = self.tkshellht - self.skliqht + self.hro_

        # Breather vent pressure setting (p_bp: pressure setting; p_bv: vacuum setting)
        self.p_bp = 0.03
        self.p_bv = -0.03

        # Returns values in degrees Rankine
        self.tax_r = self.tax + 459.7
        self.tan_r = self.tan + 459.7

        self.delta_ta_r = self.tax_r - self.tan_r

        self.taa = (self.tax_r + self.tan_r) * (1 / 2.)
        self.tb = self.taa + (0.003 * self.solarabs * self.ins)

        # Atmospheric pressure at facility, user defined from Table 7.1-7.
        self.atmplocal = atmplocal

    def atmp(self):
        return self.atmplocal

    def vq(self):
        '''
        Returns net working loss throughput, in BBL,
        as calculated by Eqn: 1-35.
        '''
        return 5.614 * (self.throughput) * (1 / 42.)

    def kb(self):
        '''
        Returns the vent setting of the storage tank.
        '''
        return self.ventsetting

    def kn(self):
        '''
        Returns the turnover factor, as calculated
        in Eqn: 1-36 and 1-37.
        '''
        delta_liquid_height = (self.hlx - self.hln)
        turnover_factor_n = self.vq()
        turnover_factor_d = ((math.pi / 4.) * (self.diameter)**2)
        turnover_factor = turnover_factor_n / turnover_factor_d

        n = turnover_factor / delta_liquid_height

        if n <= 36:
            return 1
        else:
            return n

    def kp(self):
        '''
        Returns product factor, dimensionless constant,
        based on product type.
        '''
        if self.productfactor == 'crude oils':
            return 0.75
        elif self.productfactor == 'other stocks':
            return 1.0
        else:
            raise ValueError('Incorrect product type, \
            must be either crude oils or other stocks')

    '''
    See note from Eqn: 1-35, throughput is in gal and
    converted to BBL (1 BBL = 42 Gal)
    '''

    def hvo(self):
        '''
        Returns the Vapor Space Outage in units of
        Feet as calculated from Eqn: 1-16
        '''
        return self.hvo

    def vv(self):
        '''
        Returns the Vapor Space Volume of the storage tank.
        '''
        return (math.pi) * (1 / 4.) * (self.diameter**2) * self.hvo

    def tla(self):
        '''
        Returns the Daily Average Liquid Surface
        Temperature, as caluclated from Eqn: 1-28
        '''
        return (0.4 * self.taa) + (0.6 * self.tb) + ((0.005 * self.solarabs) * (self.ins))

    def tv(self):
        '''
        Returns Ave Vapor Temp.,
        in deg R. From Eqn: 1-33
        '''
        return (0.7 * self.taa) + (0.3 * self.tb) + (0.009 * self.solarabs * self.ins)

    def deltv(self):
        '''
        Returns the Daily Average Temperature Range, in Rankine as calculated from Eqn: 1-7
        '''
        return (0.7 * self.delta_ta_r) + (0.02 * self.solarabs * self.ins)

    def tlx_r(self):
        '''
        Returns the Maximum Liquid Temperature, in Rankine as calculated from Figure 7.1-17
        '''
        return self.tla() + (0.25 * (self.deltv()))

    def tln_r(self):
        '''
        Returns the Minimum Liquid Temperature, in Rankine as calculated from Figure 7.1-17
        '''
        return self.tla() - (0.25 * (self.deltv()))

    def tlx_f(self):
        '''
        Returns the Maximum Liquid Temperature, in F as calculated from Figure 7.1-17
        Assume: 1 R − 459.67 = -458.7 F
        '''
        return self.tlx_r() - 458.67

    def tln_f(self):
        '''
        Returns the Minimum Liquid Temperature, in F as calculated from Figure 7.1-17
        Assume: 1 R − 459.67 = -458.7 F
        '''
        return self.tln_r() - 458.67

    def bventpress(self):
        '''
        Returns Breather Vent Pressure, delta_pb, calculated from Eqn: 1-10
        '''
        return self.p_bp - self.p_bv


def stockDensity(mv,
                 pva,
                 tv):
    '''
    Uses R = 10.731 psia*ft3 / lb-mole* deg R
    Returns stock density in units of lbs/ft3
    '''
    return ((mv * pva) / (10.731 * tv))


def vapPressureRange(plx,
                     pln):
    '''
    Returns the vapor pressure range, in PSIA, as calculated from Eqn: 1-10.
    '''
    return plx - pln


def vaporSpaceExpansionFactor(deltv,
                              tla,
                              delpv,
                              delbpv,
                              atmp,
                              pva):
    '''
    Returns the Vapor Space Expansion Factor (ke), as calculated from Eqn: 1-5.
    '''
    return (deltv / tla) + ((delpv - delbpv) / (atmp - pva))


def ventedVaporSpaceSatFactor(pva,
                              hvo):
    '''
    Parameters
    ----------
    pva:

    hvo:

    Returns
    -------
    The Vented Vapor Space Saturation Factor.
    '''
    return 1 / (1 + (0.053 * pva * hvo))


def calculateLosses(vq,
                    kn,
                    kp,
                    wv,
                    kb,
                    ke,
                    ks,
                    vv):

    standing = 365 * vv * wv * ke * ks
    working = vq * kn * kp * wv * kb
    total = standing + working

    return [standing, working, total]


def rankineToCelsius(r):
    return (r - 491.7) * (5 / 9.)

###########################################################################


INPUT_CITY = 'Denver, Colorado'
INPUT_TANK = [12, 8, 6]  # tkshellht, skliqht, diameter
INPUT_CONTENTS = [8450, 'other stocks', 11.5, 4.5]  # throughput, productfactor, hlx, hln

df_met_sub = df_met[['ATMOS_PRS', 'INSOL_ANN', 'CTYST', 'CITY', 'STATE', 'TAX_ANN', 'TAN_ANN']]
df_met_filter = df_met_sub[df_met_sub['CTYST'] == INPUT_CITY]

MET_LIST = df_met_filter.values.tolist()

tank = VerticalFixedRoofTank(tkshellht=INPUT_TANK[0],  # From User Data
                             skliqht=INPUT_TANK[1],  # From User Data
                             tkrfslope=0.0625,  # Default
                             diameter=INPUT_TANK[2],  # From User Data
                             ins=1491,  # Default
                             solarabs=0.25,  # From User Data
                             tax=MET_LIST[0][5],  # From Met Table
                             tan=MET_LIST[0][6],  # From Met Table
                             atmplocal=MET_LIST[0][0],  # From Met Table
                             throughput=INPUT_CONTENTS[0],  # From User Data
                             productfactor=INPUT_CONTENTS[1],  # From User Data
                             hlx=INPUT_CONTENTS[2],  # From User Data
                             hln=INPUT_CONTENTS[3],  # From User Data
                             ventsetting=1)  # Default

# Converts from Rankine to deg C.
tla_c = rankineToCelsius(tank.tla())
tlx_c = rankineToCelsius(tank.tlx_r())
tln_c = rankineToCelsius(tank.tln_r())

###########################################################################
# TODO: POPULATE MATRIX OF CHEMICAL PROPERTIES AND QUERY DF

# Comes from SQL Database and is needed to build mixture.
# User would define the components in mixture.
cyH = ['Cyclohexane', '110-82-7', 84.16, 6.845, 1203.5, 222.86]
phMe = ['Toluene', '108-88-3', 92.14, 7.017, 1377.6, 222.64]
phH = ['Benzene', '71-43-2', 78.11, 6.906, 1211.0, 220.79]
header = ['component', 'cas_no', 'mw', 'antoine_coef_a', 'antoine_coef_b', 'antoine_coef_c']

df1 = pd.DataFrame([phH, phMe, cyH], columns=header)

df1['comp_vp'] = 10**(df1['antoine_coef_a'] -
                      ((df1['antoine_coef_b']) / (tla_c +
                                                  (df1['antoine_coef_c'])))) / 51.715

df1['comp_vp_tlx'] = 10**(df1['antoine_coef_a'] -
                          ((df1['antoine_coef_b']) / (tlx_c +
                                                      (df1['antoine_coef_c'])))) / 51.715

df1['comp_vp_tln'] = 10**(df1['antoine_coef_a'] - ((df1['antoine_coef_b']) / (tln_c + (df1['antoine_coef_c'])))) / 51.715


###########################################################################
# TODO: POPULATE MATRIX OF CHEMICAL QUANTITIES AND QUERY DF

# INTEGRATE THIS CODE TO BUILD THE MIXTURE PROFILE
CHEM_LIST = ['Cyclohexane', 'Benzene', 'Toluene']
ANNUAL_QUANTITY = [101, 2812, 258]

name_ = np.array([df_chem[df_chem['NAME'].isin(CHEM_LIST)]['NAME'].tolist()])
cas_ = np.array([df_chem[df_chem['NAME'].isin(CHEM_LIST)]['CAS'].tolist()])
mw_ = np.array([df_chem[df_chem['NAME'].isin(CHEM_LIST)]['MOLWT'].tolist()])
vp_a = np.array([df_chem[df_chem['NAME'].isin(CHEM_LIST)]['VP_COEF_A'].tolist()])
vp_b = np.array([df_chem[df_chem['NAME'].isin(CHEM_LIST)]['VP_COEF_B'].tolist()])
vp_c = np.array([df_chem[df_chem['NAME'].isin(CHEM_LIST)]['VP_COEF_C'].tolist()])
den_ = np.array([df_chem[df_chem['NAME'].isin(CHEM_LIST)]['L_DENS'].tolist()])
cat_ = np.array([df_chem[df_chem['NAME'].isin(CHEM_LIST)]['CATEGORY'].tolist()])
arr = np.concatenate((name_, cat_, cas_, mw_, den_, vp_a, vp_b, vp_c), axis=0).T

df_chem_filter = pd.DataFrame(data=arr, columns=['component',
                                                 'category',
                                                 'cas_no',
                                                 'mw',
                                                 'density',
                                                 'antoine_coef_a',
                                                 'antoine_coef_b',
                                                 'antoine_coef_c'])

CHEM_LIST = ['Cyclohexane', 'Benzene', 'Toluene']
ANNUAL_QUANTITY = [101, 2812, 258]
CAS = df_chem_filter['cas_no'].values.tolist()
MIXTURE_JSON = [{'component': chem,
                 'cas': cas,
                 'amount': amt} for chem, cas, amt in zip(CHEM_LIST, CAS, ANNUAL_QUANTITY)]
###########################################################################

phH1 = ['Benzene', '71-43-2', 78.11, 2812]
phMe1 = ['Toluene', '108-88-3', 92.14, 258]
cyH1 = ['Cyclohexane', '110-82-7', 84.16, 101]
data_list = [phH1, phMe1, cyH1]
header = ['component', 'cas_no', 'mw', 'comp_amt']
df2 = pd.DataFrame(data_list, columns=header)

df = pd.DataFrame({'component': df2['component'].tolist(),
                   'cas_no': df2['cas_no'].tolist(),
                   'amount_lb': df2['comp_amt'].tolist(),
                   'mw_lbs_mol': df2['mw'].tolist(),
                   'comp_vp': df1['comp_vp'].tolist()})

df['comp_mole'] = df['amount_lb'] / df['mw_lbs_mol']
tot_moles = np.sum(df['comp_mole'].tolist())

df['comp_mole_xi'] = df['comp_mole'] / tot_moles
df['comp_partial'] = df['comp_mole_xi'] * df['comp_vp']
vp_mixture = np.sum(df['comp_partial'].tolist())

df['comp_vap_mole_frac'] = df['comp_partial'] / vp_mixture
df['comp_vapor_mw_xi'] = df['mw_lbs_mol'] * df['comp_vap_mole_frac']
vapor_mw = np.sum(df['comp_vapor_mw_xi'].tolist())

df['comp_vp_tln'] = df1['comp_vp_tln'].tolist()
df['comp_vp_tlx'] = df1['comp_vp_tlx'].tolist()

df['comp_partial_tln'] = df['comp_mole_xi'] * df['comp_vp_tln']
df['comp_partial_tlx'] = df['comp_mole_xi'] * df['comp_vp_tlx']

tot_vp_tln = np.sum(df['comp_partial_tln'].tolist())
tot_vp_tlx = np.sum(df['comp_partial_tlx'].tolist())
###########################################################################

tv = tank.tv()  # Calculated Field
deltv = tank.deltv()  # Calculated Field
tla = tank.tla()  # Calculated Field
delbpv = tank.bventpress()  # Calculated Field
atmp = tank.atmp()  # From Met Table
hvo = tank.hvo  # Calculated Field
vq = tank.vq()  # Calculated Field
kn = tank.kn()  # Calculated Field
kp = tank.kp()  # Calculated Field
kb = tank.kb()  # Calculated Field
vv = tank.vv()  # Calculated Field

wv = stockDensity(mv=vapor_mw,
                  pva=vp_mixture,
                  tv=tv)

delpv = vapPressureRange(plx=tot_vp_tlx,
                         pln=tot_vp_tln)

ke = vaporSpaceExpansionFactor(deltv=deltv,
                               tla=tla,
                               delpv=delpv,
                               delbpv=delbpv,
                               atmp=atmp,
                               pva=vp_mixture)

ks = ventedVaporSpaceSatFactor(pva=vp_mixture, hvo=hvo)

losses = calculateLosses(vq=vq,
                         kn=kn,
                         kp=kp,
                         wv=wv,
                         kb=kb,
                         vv=vv,
                         ke=ke,
                         ks=ks)

print('Total Losses from Tank: {:.4f} lbs/year'.format(losses[2]))
print('Total Working Losses from Tank: {:.4f} lbs/year'.format(losses[1]))
print('Total Standing Losses from Tank: {:.4f} lbs/year'.format(losses[0]))
