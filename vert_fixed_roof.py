import math
import pandas as pd
import numpy as np

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
                 atm_p_local,
                 throughput,
                 product_factor,
                 hlx,
                 hln,
                 ventsetting):

        '''
        tax: Average Daily Maximum Ambient Temperature, user defined, in deg F.
        tan: Average Daily Minimum Ambient Temperature, user defined, in deg F.
        solarabs: Tank surface solar absorptance, user defined, depends on tank color.
        ins: Average daily total insulation, defined to be 1,491 (btu/ft^2*d) (from Table 7.1-7)
        '''
        self.tkshellht = tkshellht
        self.skliqht = skliqht
        self.tkrfslope = tkrfslope
        self.diameter = diameter
        self.ins = ins
        self.solarabs = solarabs
        self.throughput = throughput
        self.product_factor = product_factor
        self.hlx = hlx
        self.hln = hln
        self.ventsetting = 1

        #Values are in deg. F
        self.tax = tax
        self.tan = tan

        self.rad_ = (1/2.)*self.diameter
        self.hro_ = (1/3.)*self.rad_*self.tkrfslope

        self.hvo = self.tkshellht - self.skliqht + self.hro_

        #Breather vent pressure setting (p_bp: pressure setting; p_bv: vacuum setting)
        self.p_bp = 0.03
        self.p_bv = -0.03

        #Returns values in degrees Rankine
        self.tax_r = self.tax + 459.7
        self.tan_r = self.tan + 459.7

        self.delta_ta_r = self.tax_r - self.tan_r

        self.taa = (self.tax_r + self.tan_r) * (1/2.)
        self.tb = self.taa + (0.003*self.solarabs*self.ins)

        #Atmospheric pressure at facility, user defined from Table 7.1-7.
        self.atm_p_local = atm_p_local

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
        turnover_factor_n = (5.614*(self.throughput)*(1/42.))
        turnover_factor_d = ((math.pi/4.)*(self.diameter)**2)
        turnover_factor = turnover_factor_n/turnover_factor_d

        n = turnover_factor/delta_liquid_height

        if n <= 36:
            return 1
        else:
            return n

    def kp(self):
        '''
        Returns product factor, dimensionless constant,
        based on product type.
        '''
        if self.product_factor == 'crude oils':
            return 0.75
        elif self.product_factor == 'other stocks':
            return 1.0
        else:
            raise ValueError('Incorrect product type, \
            must be either crude oils or other stocks')

    '''
    See note from Eqn: 1-35, throughput is in gal and
    converted to BBL (1 BBL = 42 Gal)
    '''
    def vq(self):
        '''
        Returns net working loss throughput, in BBL,
        as calculated by Eqn: 1-35.
        '''
        return 5.614*(self.throughput)*(1/42.)

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
        return (math.pi)*(1/4.)*(self.diameter**2)*self.hvo

    def tla(self):
        '''
        Returns the Daily Average Liquid Surface Temperature, as caluclated from Eqn: 1-28
        '''
        return (0.4*self.taa) + (0.6*self.tb) + ((0.005*self.solarabs)*(self.ins))

    def tv(self):
        '''
        Returns Average Vapor Temperature, in Rankine as caluclated from Eqn: 1-33
        '''
        return (0.7*self.taa) + (0.3*self.tb) + (0.009*self.solarabs*self.ins)

    def deltv(self):
        '''
        Returns the Daily Average Temperature Range, in Rankine as calculated from Eqn: 1-7
        '''
        return (0.7*self.delta_ta_r) + (0.02*self.solarabs*self.ins)

    def tlx_r(self):
        '''
        Returns the Maximum Liquid Temperature, in Rankine as calculated from Figure 7.1-17
        '''
        return ((0.4*self.taa) + (0.6*self.tb) + ((0.005*self.solarabs)*(self.ins))) + (0.25*((0.7*self.delta_ta_r) + (0.02*self.solarabs*self.ins)))

    def tln_r(self):
        '''
        Returns the Minimum Liquid Temperature, in Rankine as calculated from Figure 7.1-17
        '''
        return ((0.4*self.taa) + (0.6*self.tb) + ((0.005*self.solarabs)*(self.ins))) - (0.25*((0.7*self.delta_ta_r) + (0.02*self.solarabs*self.ins)))

    def tlx_f(self):
        '''
        Returns the Maximum Liquid Temperature, in F as calculated from Figure 7.1-17
        Assume: 1 R − 459.67 = -458.7 F
        '''
        return (((0.4*self.taa) + (0.6*self.tb) + ((0.005*self.solarabs)*(self.ins))) + (0.25*((0.7*self.delta_ta_r) + (0.02*self.solarabs*self.ins)))) - 458.67

    def tln_f(self):
        '''
        Returns the Minimum Liquid Temperature, in F as calculated from Figure 7.1-17
        Assume: 1 R − 459.67 = -458.7 F
        '''
        return (((0.4*self.taa) + (0.6*self.tb) + ((0.005*self.solarabs)*(self.ins))) - (0.25*((0.7*self.delta_ta_r) + (0.02*self.solarabs*self.ins)))) - 458.67

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
    return ((mv*pva) / (10.731 *tv))

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
    return (deltv/tla) + ((delpv - delbpv) / (atmp - pva))

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
    return 1/(1 + (0.053*pva*hvo))

def calculateStandingLosses(vv,
                            wv,
                            ke,
                            ks):
    '''
    Parameters
    ----------
    vv:

    wv:

    ke:

    ks:

    Returns
    -------
    Standing losses from the storage tank in lbs/year.
    '''
    return 365*vv*wv*ke*ks

def calculateWorkingLosses(vq,
                           kn,
                           kp,
                           wv,
                           kb):
    '''
    Parameters
    ----------
    vq:

    kn:

    kp:

    wv:

    kb:

    Returns
    -------
    Working losses from the storage tank in lbs/year.
    '''
    return vq*kn*kp*wv*kb

tank = VerticalFixedRoofTank(tkshellht=12,
                             skliqht=8,
                             tkrfslope=0.0625,
                             diameter=6,
                             ins=1491,
                             solarabs=0.25,
                             tax=63.5,
                             tan=37.9,
                             atm_p_local=12.08,
                             throughput=8450,
                             product_factor='other stocks',
                             hlx=11.5,
                             hln=4.5,
                             ventsetting=1)

#Converts from Rankine to deg C.
tla_c = (tank.tla() - 491.7) * (5/9.)

#Converts from Rankine to deg C.
tlx_c = (tank.tlx_r() - 491.7) * (5/9.)
tln_c = (tank.tln_r() - 491.7) * (5/9.)

'''
Comes from SQL Database and is needed to build mixture.
User would define the components in mixture.
'''
cyH = ['Cyclohexane', '110-82-7', 84.16, 6.845, 1203.5, 222.86]
phMe = ['Toluene', '108-88-3', 92.14, 7.017, 1377.6, 222.64]
phH = ['Benzene', '71-43-2', 78.11, 6.906, 1211.0, 220.79]
header = ['component', 'cas_no', 'mw', 'antoine_coef_a', 'antoine_coef_b', 'antoine_coef_c']

df1 = pd.DataFrame([phH, phMe, cyH], columns=header)

df1['comp_vp'] = 10**(df1['antoine_coef_a'] - \
                 ((df1['antoine_coef_b']) / (tla_c + \
                 (df1['antoine_coef_c'])))) / 51.715

df1['comp_vp_tlx'] = 10**(df1['antoine_coef_a'] - \
                         ((df1['antoine_coef_b']) / (tlx_c + \
                         (df1['antoine_coef_c'])))) / 51.715

df1['comp_vp_tln'] = 10**(df1['antoine_coef_a'] - ((df1['antoine_coef_b']) / (tln_c + (df1['antoine_coef_c']))))/ 51.715

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

wv = stockDensity(mv=vapor_mw,
                  pva=vp_mixture,
                  tv=tank.tv())

delpv = vapPressureRange(plx=tot_vp_tlx,
                         pln=tot_vp_tln)

ke = vaporSpaceExpansionFactor(deltv=tank.deltv(),
                          tla=tank.tla(),
                          delpv=delpv,
                          delbpv=tank.bventpress(),
                          atmp=12.08,
                          pva=vp_mixture)

ks = ventedVaporSpaceSatFactor(pva=vp_mixture, hvo=tank.hvo)

standing = calculateStandingLosses(vv=tank.vv(),
                                   wv=wv,
                                   ke=ke,
                                   ks=ks)

working = calculateWorkingLosses(vq=tank.vq(),
                                 kn=tank.kn(),
                                 kp=tank.kp(),
                                 wv=wv,
                                 kb=tank.kb())

total_losses = standing + working
print('Total Losses from Tank: {:.3f}'.format(total_losses))
