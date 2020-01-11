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

        return (self.tla() - 491.7) * (5./9.)

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

        return (self.tlx_r() - 491.7) * (5./9.)

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
        return (self.tln_r() - 491.7) * (5./9.)

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