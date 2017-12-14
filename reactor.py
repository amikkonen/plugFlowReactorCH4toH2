#!/usr/bin/env python3
# -*- coding: utf-8 -*-
try:
    from IPython import get_ipython
    get_ipython().magic('reset -sf')
except:
    pass
"""
Created on Thu Nov 30 12:11:40 2017

Based on Tiina Keipi's PhD thesis Technology Development and Techno-Economic 
Analysis of Hydrogen Production by Thermal Decomposition of Methane.

The nesessary information is also available in Keipi's paper

    Methane thermal decomposition in regenerative heat exchanger reactor: 
    Experimental and modeling study. / Keipi, Tiina; Li, Tian; Løvås , 
    Terese; Tolvanen, Henrik; Konttinen, Jukka. Energy, Vuosikerta 135, 
    15.09.2017, s. 823-832.

Useful "switches":
    In the __main__ section at the end of the code you can choose between basic
    run, optimization with pso, and optimization results post processing.

        reactor = main()
    #    optimize()
    #    plot_optimization()

    You can also choose the number od cases studies in the main with the 
    
        nros = range(1, len(cases_measured_data))
   #    nros = range(46, 47) # One case

    lines. 
    
    If you study a lot of cases you'll propaply want to leave the plotting
    function commented as shown below.

        cases = []
        for k in nros:
            case = cases_measured_data[k]
            case.update(common)
        
            ##################
            reactor = Reactor(case)
            reactor.solve()
            reactor.post()
    #        reactor.plot_all()
            print(reactor)

    If on the other hand you want detailed information same indivudual case,
    choose the case with nros and uncommnet the plotting.

@author: Antti Mikkonen, a.mikkonen@iki.fi
"""

import time
import scipy as sp
from matplotlib import pyplot as plt
from scipy import interpolate
from scipy import optimize
from xlrd import open_workbook
from CoolProp.CoolProp import PropsSI # Only used for viscosity, you can easily replace the dependacy


###############################################################################
def dm3minTom3s(dm3min):
    return dm3min*1.66666667e-5 #1e-3/60

def rhoIdeal(T,p, M):
    R_u = 8.314459848 # J/mol/K
    return (p*M) / (R_u*T)

def rhoN2(T, p):
    M = 28.02e-3 # kg/mol
    return rhoIdeal(T,p, M)

def rhoCH4(T, p):
    M = 16.04e-3 # kg/mol
    return rhoIdeal(T,p, M)

def rhoH2(T, p):
    M = 2.0158814e-3 # kg/mol
    return rhoIdeal(T,p, M)

###############################################################################

class Reactor(object):
    def __init__(self, parms):
        self.label                          = parms["label"]
        self.T_C_m                          = parms["T_C"]
        self.x_m                            = parms["x"]
        self.Q_CH4_in                       = dm3minTom3s(parms["Q_CH4_in"])
        self.Q_N2_in                        = dm3minTom3s(parms["Q_N2_in"])
        self.CH4_conversion_rate_measured   = parms["CH4_conversion_rate_measured"]
        self.d                              = parms["d"]
        self.L_tot                          = parms["L_tot"]
        self.eps                            = parms["eps"]
        self.p_in                           = parms["p_in"]
        self.Af                             = parms["Af"]
        self.Ab                             = parms["Ab"]
        self.Eaf                            = parms["Eaf"]
        self.Eab                            = parms["Eab"]
        self.nf                             = parms["nf"]
        self.mb                             = parms["mb"]
        self.n                              = parms["n"]
        self.d_p                            = parms["d_p"]   
        
        # Constants
        self.R_u    = 8.314459848   # J/mol/K
        self.M_CH4  = 16.04e-3      # kg/mol
        self.M_N2   = 28.02e-3      # kg/mol
        self.M_H2   = 2.0158814e-3  # kg/mol
        self.M_Cs   = 12.0107e-3    # kg/mol
        
        T_C_in = 20 # Arbitrary
        dT_end = -30

        self.T_C_m = sp.array([T_C_in] + list(self.T_C_m) + [self.T_C_m[-1] + dT_end])
        self.x_m   = sp.array([0] + list(self.x_m) + [self.L_tot])
        self.T_m = self.T_C_m + 273.15

        # Loop variable
        self.k = 0

        # Interpolate values
        self.x = sp.linspace(0, self.L_tot, self.n)
        self.T =  interpolate.interp1d(self.x_m, self.T_m, kind='linear')(self.x)  

        # Calculate values
        self.Ac = 0.25*sp.pi*self.d**2
        self.dx = self.L_tot / self.n
        
        self.Lc = self.d_p*(self.eps / (1-self.eps))

        self.mu = sp.array([PropsSI("V","T", T, "P", self.p_in, "N2") 
                            for T in self.T])
        
        # Units cm3, J, s
        self.kf     = self.Arrheniues_forward()
        self.kb     = self.Arrheniues_backward()

        # Storages for components run time components
        self.x_CH4  = sp.zeros_like(self.x, dtype="float")
        self.x_N2   = sp.zeros_like(self.x, dtype="float")
        self.x_H2   = sp.zeros_like(self.x, dtype="float")
#        self.x_Cs   = sp.zeros_like(self.x, dtype="float")
        
        self.y_CH4  = sp.zeros_like(self.x, dtype="float")
        self.y_N2   = sp.zeros_like(self.x, dtype="float")
        self.y_H2   = sp.zeros_like(self.x, dtype="float")

        self.c_CH4  = sp.zeros_like(self.x, dtype="float")
        self.c_N2   = sp.zeros_like(self.x, dtype="float")
        self.c_H2   = sp.zeros_like(self.x, dtype="float")
        self.c_Cs   = sp.zeros_like(self.x, dtype="float")        


        # Storages for mixture        
        self.M      = sp.zeros_like(self.x, dtype="float") # Mean molar mass kg/m3
        self.rho    = sp.zeros_like(self.x, dtype="float")
        self.Vc     = sp.zeros_like(self.x, dtype="float")
        self.Q      = sp.zeros_like(self.x, dtype="float")
        self.N_g      = sp.zeros_like(self.x, dtype="float")
        self.m_g    = sp.zeros_like(self.x, dtype="float")
        
        self.r      = sp.zeros_like(self.x, dtype="float") # Reaction rate, mol/s
        self.dt     = sp.zeros_like(self.x, dtype="float") # Pass through time to next (s)

        self.p      = sp.zeros_like(self.x, dtype="float") 
        
        self.rhoCH4 = sp.zeros_like(self.x, dtype="float") 
        self.rhoN2  = sp.zeros_like(self.x, dtype="float") 
        self.rhoH2  = sp.zeros_like(self.x, dtype="float") 

        # Inlet
        self.p[0]      = self.p_in
        self.rhoCH4[0] = rhoCH4(self.T[0], self.p[0]) 
        self.rhoN2[0]  = rhoN2(self.T[0], self.p[0])
        self.rhoH2[0]  = rhoH2(self.T[0], self.p[0])

        m_CH4_in    = self.Q_CH4_in * self.rhoCH4[0]
        m_N2_in     = self.Q_N2_in  * self.rhoN2[0]
        N_CH4_in    = m_CH4_in / self.M_CH4
        self.N_N2   = m_N2_in  / self.M_N2
        N_in        = N_CH4_in + self.N_N2
        
        self.m_N2   = self.N_N2 * self.M_N2
        
        self.m_tot  = m_CH4_in + m_N2_in
        y_CH4_in    = m_CH4_in / self.m_tot
        y_N2_in     = m_N2_in / self.m_tot
        
        self.x_CH4[0]  = N_CH4_in / N_in
        self.x_N2[0]   = self.N_N2 / N_in

        self.M[0]      = self.x_CH4[0]*self.M_CH4 + self.x_N2[0]*self.M_N2
        self.rho[0]    = self.p[0] * self.M[0] / (self.R_u * self.T[0])
    
    
        self.Vc[0]     = self.m_tot / (self.rho[0]*self.eps*self.Ac)

        self.dt[0]     = self.dx / self.Vc[0]
    
    
        self.c_CH4[0] = y_CH4_in*self.rho[0] / self.M_CH4
        self.c_N2 [0] = y_N2_in*self.rho[0]  / self.M_N2
        self.c_H2[self.k] = 0
#        self.c_Cs[self.k] = 0
        
        self.m_g[0]    = self.m_tot
        
        
        
        # TEST
        self.N_CH4 = sp.zeros_like(self.x, dtype="float")
        self.N_H2  = sp.zeros_like(self.x, dtype="float")
        self.N_Cs  = sp.zeros_like(self.x, dtype="float")
        
        self.m_CH4 = sp.zeros_like(self.x, dtype="float")
        self.m_H2  = sp.zeros_like(self.x, dtype="float")
        self.m_Cs  = sp.zeros_like(self.x, dtype="float")    
        
        self.m_CH4[0] = m_CH4_in
        self.N_CH4[0] = N_CH4_in
        self.N_g[0]   = N_in
        self.Q[0]     = self.Q_CH4_in + self.Q_N2_in


    def __str__(self):
        out = \
"""        
Label %i
Rate (measured)   %.4f
Rate (calculated) %.4f
Difference        %.2f %s

""" % (self.label, self.CH4_conversion_rate_measured, self.CH4_conversion_rate,
       100*self.diff_in_rates,
       "%")
        
        return out

    
    def post(self):
        self.CH4_conversion_rate = (self.N_CH4[0]-self.N_CH4[-1])/self.N_CH4[0]
        self.diff_in_rates = (self.CH4_conversion_rate \
                              - self.CH4_conversion_rate_measured)\
                             / self.CH4_conversion_rate_measured
    
    
    def solve(self):

        while self.k < self.n-1:
            self.progress_one_step()
            self.k += 1


    def progress_one_step(self):
        
        k  = self.k
        kn = self.k + 1
        self.r[k] = -(self.kf[k] * (self.c_CH4[k]*1e-6)**self.nf 
                     - self.kb[k] * (self.c_H2[k]*1e-6)**self.mb
                    )*1e6
                    
        dcCH4 = self.r[k] * self.dt[k]
        dcCs  = -dcCH4
        dcH2  = -2*dcCH4             
        
        self.N_CH4[kn] = (self.c_CH4[k]+dcCH4)*self.Q[k]
        self.N_H2[kn]  = (self.c_H2[k] +dcH2) *self.Q[k]
        self.N_Cs[kn]  = (self.c_Cs[k] +dcCs) *self.Q[k]
        self.N_g[kn]   = self.N_CH4[kn] + self.N_H2[kn] + self.N_N2

        self.x_CH4[kn] = self.N_CH4[kn] / self.N_g[kn]
        self.x_H2[kn]  = self.N_H2[kn] / self.N_g[kn]
        self.x_N2[kn]  = self.N_N2 / self.N_g[kn]

        assert sp.isclose(self.x_CH4[kn] + self.x_H2[kn] + self.x_N2[kn], 1)
        
        self.M[kn] = self.x_CH4[kn]*self.M_CH4 + self.x_N2[kn]*self.M_N2 \
                   + self.x_H2[kn]*self.M_H2 

        self.p[kn] = self.p[k] - self.dx * \
                                 (150*self.mu[k]*self.Vc[k]/self.Lc**2 
                                  + 1.75*self.rho[k]*self.Vc[k]**2/self.Lc
                                  )

        self.rho[kn] = self.p[kn] * self.M[kn] / (self.R_u * self.T[kn])

        self.m_CH4[kn] = self.N_CH4[kn] * self.M_CH4
        self.m_H2[kn]  = self.N_H2[kn] * self.M_H2
        self.m_Cs[kn]  = self.N_Cs[kn] * self.M_Cs
        self.m_g[kn]   = self.m_CH4[kn] + self.m_H2[kn] + self.m_N2

        self.rhoCH4[kn] = rhoCH4(self.T[kn],self.p[kn])
        self.rhoH2[kn]  = rhoH2(self.T[kn],self.p[kn])
        self.rhoN2[kn]  = rhoN2(self.T[kn],self.p[kn])

        self.Q[kn] = self.m_CH4[kn] / self.rhoCH4[kn] \
                   + self.m_H2[kn]  / self.rhoH2[kn] \
                   + self.m_N2      / self.rhoN2[kn]

        self.Vc[kn]  = self.Q[kn] / (self.eps*self.Ac)# / 0.34 #/ 0.37
        self.dt[kn]  = self.dx / self.Vc[kn]

        self.c_CH4[kn] = self.N_CH4[kn] / self.Q[kn]
        self.c_H2[kn]  = self.N_H2[kn] / self.Q[kn]
        self.c_N2[kn]  = self.N_N2 / self.Q[kn]
        self.c_Cs[kn]  = self.N_Cs[kn] / self.Q[kn]

    def Arrheniues_forward(self):
        return self.Af*sp.exp(-self.Eaf/(self.R_u*self.T))
    
    def Arrheniues_backward(self):
        return self.Ab*sp.exp(-self.Eab/(self.R_u*self.T))    

    def plot_T(self):#x, T, caseNro=""):
        fig, ax = plt.subplots()
        ax.plot(self.x_m,self.T_m, "k-d")
        ax.plot([self.x[0], self.x_m[0]],[self.T[0], self.T_m[0]], "b--")
        ax.plot([self.x_m[-1], self.x[-1]],[self.T_m[-1], self.T[-1]], "b--")
        ax.plot(self.x[-1],self.T[-1], "bd")
        ax.plot(self.x[0],self.T[0], "bd")
        ax.set_title("Case " + str(self.label))
        ax.set_xlabel("x (m)")
        ax.set_ylabel("T (K)")
        ax.set_xticks([0]+list(self.x_m)+[self.L_tot])
        ax.set_xlim(0, self.L_tot)
        for k in range(len(self.x_m)):
            ax.text(self.x_m[k], self.T_m[k]+20, self.T_m[k])
            
            
        ax.set_ylim(None, self.T_m.max()+100)    
        ax.grid()
        fig.savefig("temperatureDistribution.pdf")
    

    def plot_ks(self):#x, T, caseNro=""):
        fig, axes = plt.subplots(2)
        axes[0].plot(self.x, self.kf, 'b-')
        axes[1].plot(self.x, self.kb, 'y-')
        axes[0].set_title("Case " + str(self.label))
        axes[1].set_xlabel("x (m)")
        axes[0].set_ylabel("kf")        
        axes[1].set_ylabel("kb")        
        for ax in axes:
            ax.grid()
        fig.savefig("ks.pdf")    
      
    def plot_concentrations(self):
        fig, axs = plt.subplots(4)
        axs[0].plot(self.x, self.c_CH4, label="c CH4")
        axs[1].plot(self.x, self.c_H2, label="c H2")
        axs[2].plot(self.x, self.c_Cs, label="c Cs")        
        axs[3].plot(self.x, self.c_N2, label="c N2")        
        
        
        for ax in axs:
            ax.legend()
            ax.set_xlim(0, self.L_tot)
        fig.savefig("consentrations.pdf")
    

    def plot_massflow(self):
        fig, axs = plt.subplots(4)
        axs[0].plot(self.x, self.m_CH4, label="m CH4")
        axs[1].plot(self.x, self.m_H2, label="m H2")
        axs[2].plot(self.x, self.m_Cs, label="m Cs")        
        axs[3].plot(self.x, sp.ones_like(self.x)*self.m_N2, label="m N2")        
        
        
        for ax in axs:
            ax.legend()
            ax.set_xlim(0, self.L_tot)
        
        fig.savefig("massFlow.pdf")    
  
    def plot_rate(self):
        fig, axs = plt.subplots(2)
        
        ax = axs[0]
        ax.plot(self.x_m,self.T_m, "k-d")
        ax.plot([self.x[0], self.x_m[0]],[self.T[0], self.T_m[0]], "b--")
        ax.plot([self.x_m[-1], self.x[-1]],[self.T_m[-1], self.T[-1]], "b--")
        ax.plot(self.x[-1],self.T[-1], "bd")
        ax.plot(self.x[0],self.T[0], "bd")
        ax.set_title("Case " + str(self.label))
        ax.set_xlabel("x (m)")
        ax.set_ylabel("T (K)")
        ax.set_xticks([0]+list(self.x_m)+[self.L_tot])
        ax.set_xlim(0, self.L_tot)
        for k in range(len(self.x_m)):
            ax.text(self.x_m[k], self.T_m[k]+20, self.T_m[k])
        ax.set_ylim(None, self.T_m.max()+100)    
        
        
        ax = axs[1]
        ax.plot(self.x, self.r, label="r")
        ax.axhline(0, color='k')
        ax.legend()
        
        fig.savefig("rate.pdf")    
    
    def plot_all(self):
        self.plot_rate()
        self.plot_concentrations()
        self.plot_massflow()
        
        
###############################################################################

def read_experimental():
    wb = open_workbook('keipiData.xls')
    sheet = wb.sheets()[0]

    first_row = 1
    last_row = 49
    cases = []
    
#    for col in range(0,12):
#        print "|", sheet.cell(first_row,col).value, 
    
    for row in range(first_row, last_row):
        case = {}
        case["label"] = int(sheet.cell(row,0).value)
        case["Q_CH4_in"] = sheet.cell(row,1).value
        case["Q_N2_in"] = sheet.cell(row,2).value
        T_C = []
        for k in range(3,11):
            T_C.append(sheet.cell(row,k).value)
        case["T_C"] = T_C
        case["CH4_conversion_rate_measured"] = sheet.cell(row,11).value
#        print case
        cases.append(case)
    return cases


    
def optimize():
    from pyswarm import pso
    
    def to_min(parms):
        Af  = parms[0] # 8.5708e12 # forward
        Ab  = parms[1] # 1.1190e7  # backward
        Eaf = parms[2] # 337.12e3 # forward J/mol
        Eab = parms[3] # 243.16e3 # backward J/mol
        nf  = parms[4] # 1.123     # forward
        mb  = parms[5] # 0.9296    # backward
        
        x = sp.array([0.4, 0.6, 0.9, 1.2, 1.4, 1.7, 2, 2.4]) # m
        d = 73e-3 # m
        eps = 0.4 # Void fraction, empty/total
        L_tot = 2.5
        p_in = 1e5
        n = 20
    
        common = {
                "x"        : x,
                "d"        : d,
                "L_tot"    : L_tot,
                "eps"      : eps,
                "p_in"     : p_in, 
                "n"        : n,
                "Af"       : Af,
                "Ab"       : Ab,
                "Eaf"      : Eaf,
                "Eab"      : Eab,
                "nf"       : nf,
                "mb"       : mb
                }
    
        cases_measured_data = read_experimental()
        
        nros = range(1, len(cases_measured_data))
    
    
        diffs = sp.zeros_like(nros)
        for k in nros:
            case = cases_measured_data[k]
            case.update(common)
        
            try:
                ##################
                reactor = Reactor(case)
                reactor.solve()
                reactor.post()
                
                ##################
                
                diffs[k-1] = reactor.diff_in_rates * 100
            except AssertionError:
                diffs[k-1] = max(10000, diffs.max())
        
            
        rms = sp.sqrt(sp.mean(sp.square(diffs)))
        
        with open("optimization.txt", "a") as ofile:
            out = ""
            for parm in parms:
                out += str(parm) + "\t" 
            out += str(rms) + "\n"
        
            ofile.write(out)
        
        return rms

    parms    = sp.zeros(6)
    parms[0] = 8.5708e12 # forward
    parms[1] = 1.1190e7  # backward
    parms[2] = 337.12e3 # forward J/mol
    parms[3] = 243.16e3 # backward J/mol
    parms[4] = 1.123     # forward
    parms[5] = 0.9296    # backward
    print( "er", to_min(parms))
    
    lb = parms*0.5
    ub = parms*2

    xopt, fopt = pso(to_min, lb, ub)
    print(xopt)
    print(fopt)

def plot_optimization():
    parms    = sp.zeros(6)
    parms[0] = 8.5708e12 # forward
    parms[1] = 1.1190e7  # backward
    parms[2] = 337.12e3 # forward J/mol
    parms[3] = 243.16e3 # backward J/mol
    parms[4] = 1.123     # forward
    parms[5] = 0.9296    # backward
    
    with open("optimization.txt", "r") as ifile:
        lines = ifile.readlines()
    
    ers = sp.zeros(len(lines))
    for k, line in enumerate(lines):
        parts = [float(part) for part in line.split()]
        
        ers[k] = parts[-1]
    
    fig, axs = plt.subplots(2)
    markersize=0.1
    axs[0].semilogy(ers, 'ko',markersize=markersize)
    axs[1].plot(ers, 'ko',markersize=markersize)
    axs[1].set_ylim(0, 50)
    axs[0].grid()
    axs[1].grid()
    fig.savefig("optimization.pdf")
    
def main():
    x = sp.array([0.4, 0.6, 0.9, 1.2, 1.4, 1.7, 2, 2.4]) # m
    d = 73e-3 # m
    eps = 0.4 # Void fraction, empty/total
    L_tot = 2.5
    p_in = 2e5
    Af = 8.5708e12 # forward
    Ab = 1.1190e7  # backward
    Eaf = 337.12e3 # forward J/mol
    Eab = 243.16e3 # backward J/mol
    nf = 1.123     # forward
    mb = 0.9296    # backward
    d_p = 0.01


    n = 100

    common = {
            "x"        : x,
            "d"        : d,
            "L_tot"    : L_tot,
            "eps"      : eps,
            "p_in"     : p_in, 
            "n"        : n,
            "Af"       : Af,
            "Ab"       : Ab,
            "Eaf"      : Eaf,
            "Eab"      : Eab,
            "nf"       : nf,
            "mb"       : mb,
            "d_p"      : d_p  
            }

    cases_measured_data = read_experimental()
    
    nros = range(1, len(cases_measured_data))
#    nros = range(46, 47) # One case

    cases = []
    for k in nros:
        case = cases_measured_data[k]
        case.update(common)
    
        ##################
        reactor = Reactor(case)
        reactor.solve()
        reactor.post()
#        reactor.plot_all()
        print(reactor)
        
        ##################
        
        cases.append(reactor)
    
    diffs = sp.zeros(len(cases))
    calc = sp.zeros(len(cases))
    meas = sp.zeros(len(cases))
    for k,case in enumerate(cases):
        diffs[k] = case.diff_in_rates * 100
        meas[k]  = case.CH4_conversion_rate_measured
        calc[k]  = case.CH4_conversion_rate
        
        
    fig, ax = plt.subplots()
    ax.plot(nros, calc, 'o',label="Calculated")
    ax.plot(nros, meas, 'd', label="Measured")
    ax.legend()
    fig.show()
    fig.savefig("comparison.pdf")

###############################################################################    

if __name__ == "__main__":
    print("START")
    start = time.time()
    reactor = main()
#    optimize()
#    plot_optimization()
    
    print("END", time.time() - start, "s")
