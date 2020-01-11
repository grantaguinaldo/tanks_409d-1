import math
import pandas as pd
import numpy as np
import os


def filterMetList(df, input_city):
    df_sub = df[['ATMOS_PRS', 'INSOL_ANN',
                 'CTYST', 'CITY', 'STATE',
                 'TAX_ANN', 'TAN_ANN']]
    df_sub_filter = df_sub[df_sub['CTYST'] == input_city]

    return df_sub_filter.values.tolist()


def solarabsLookUp(df, color, shade, condition):
    if len(df[(df['Surface_Color'] == color) & (df['Shade'] == shade)][condition]) != 0:
        return float(df[(df['Surface_Color'] == color) & (df['Shade'] == shade)][condition])
    else:
        raise ValueError('Error: Color, Shade, Condition combination is not valid.')


def calculation(df, chem_list, annual_qty, tank, file_name):

    name_ = np.array([df[df['NAME'].isin(chem_list)]['NAME'].tolist()])
    cas_ = np.array([df[df['NAME'].isin(chem_list)]['CAS'].tolist()])
    mw_ = np.array([df[df['NAME'].isin(chem_list)]['MOLWT'].tolist()])
    vp_a = np.array([df[df['NAME'].isin(chem_list)]['VP_COEF_A'].tolist()])
    vp_b = np.array([df[df['NAME'].isin(chem_list)]['VP_COEF_B'].tolist()])
    vp_c = np.array([df[df['NAME'].isin(chem_list)]['VP_COEF_C'].tolist()])
    arr = np.concatenate((name_, cas_, mw_, vp_a, vp_b, vp_c), axis=0).T

    df1 = pd.DataFrame(data=arr, columns=['component',
                                          'cas_no',
                                          'mw',
                                          'antoine_coef_a',
                                          'antoine_coef_b',
                                          'antoine_coef_c'])

    df1['comp_amt'] = annual_qty

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

    calc = EmissionCalculations(mv=vapor_mw,
                                pva=vp_mixture,
                                tv=tv,
                                plx=tot_vp_tlx,
                                pln=tot_vp_tln,
                                deltv=deltv,
                                tla=tla,
                                delbpv=delbpv,
                                atmp=atmp,
                                hvo=hvo,
                                vq=vq,
                                kn=kn,
                                kp=kp,
                                kb=kb,
                                vv=vv)

    df1['stand_loss_xi'] = df1['vap_wt_xi'] * calc.standingLosses()
    df1['work_loss_xi'] = df1['vap_wt_xi'] * calc.workingLosses()
    df1['total_loss_xi'] = df1['vap_wt_xi'] * calc.totalLosses()

    print('Exporting dataframe ...')

    return df1.to_html(os.path.join('templates', file_name))


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

    def tla_c(self):

        return (self.tla() - 491.7) * (5. / 9.)

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

    def tlx_c(self):

        return (self.tlx_r() - 491.7) * (5. / 9.)

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

    def tln_c(self):
        return (self.tln_r() - 491.7) * (5. / 9.)

    def bventpress(self):
        '''
        Returns Breather Vent Pressure, delta_pb, calculated from Eqn: 1-10
        '''
        return self.p_bp - self.p_bv


class EmissionCalculations:
    def __init__(self,
                 mv,
                 pva,
                 tv,
                 plx,
                 pln,
                 deltv,
                 tla,
                 delbpv,
                 atmp,
                 hvo,
                 vq,
                 kn,
                 kp,
                 kb,
                 vv):

        self.mv = mv
        self.pva = pva
        self.tv = tv
        self.plx = plx
        self.pln = pln
        self.deltv = deltv
        self.tla = tla
        self.delbpv = delbpv
        self.atmp = atmp
        self.hvo = hvo
        self.vq = vq
        self.kn = kn
        self.kp = kp
        self.kb = kb
        self.vv = vv

    def vapPressureRange(self):
        return self.plx - self.pln

    def stockDensity(self):
        return ((self.mv * self.pva) / (10.731 * self.tv))

    def vaporSpaceExpansionFactor(self):
        return (self.deltv / self.tla) + ((self.vapPressureRange() - self.delbpv) / (self.atmp - self.pva))

    def ventedVaporSpaceSatFactor(self):
        return 1 / (1 + (0.053 * self.pva * self.hvo))

    def vapPressureRange(self):

        return self.plx - self.pln

    def standingLosses(self):

        return 365 * self.vv * self.stockDensity() * self.vaporSpaceExpansionFactor() * self.ventedVaporSpaceSatFactor()

    def workingLosses(self):

        return self.vq * self.kn * self.kp * self.stockDensity() * self.kb

    def totalLosses(self):

        return self.standingLosses() + self.workingLosses()
