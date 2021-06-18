#!/usr/bin/env python

import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LinearSegmentedColormap
from textwrap import wrap
from lmfit.models import PolynomialModel, GaussianModel, PseudoVoigtModel, ExponentialModel, Model
import scipy as scipy

#General Functions

def approveChoice(InputQuestionString, AcceptableAnswerList):
    works = False
    while works == False:
        answer = input(InputQuestionString)
        if answer in AcceptableAnswerList:
            works = True
        else:
            print('Not a valid input')
            works = False
    return answer

def WLtoIndex(wl, w):
    details = min(enumerate(wl.iloc[:,0]), key=lambda x: abs(x[1]-w))
    return details[0]


def polyfit(SSorTR, df, wl, order, wlLow, wlHigh):
    k = WLtoIndex(wl, wlLow)
    j = WLtoIndex(wl, wlHigh)
    if SSorTR == 'SS':
        return ss_polyfit(df,wl,order,j,k)
    elif SSorTR == 'TR':
        return tr_polyfit(df,wl,order,j,k)
    else:
        return 'SSorTR'

def tr_polyfit(df, wl, order, j, k):
    wl_full = wl.values.T[0]
    h = 0
    l = df.shape[1]
    wl = wl_full[h:l]
    wvn1 = wl_full[h:j]
    wvn2 = wl_full[k:l]
    wvn = np.append(wvn1,wvn2)
    new_list = []
    for item in wl:
        new_list = np.append(new_list, str(item))
    
    new_df = pd.DataFrame()
    
    for m in np.arange(len(df.index)):
        data2 = df.iloc[m, h:l].to_numpy()
        dataA = df.iloc[m, h:j]
        dataB = df.iloc[m, k:l]
        data3 = np.append(dataA,dataB)

        coeff = np.polyfit(wvn, data3, order)
        correction = np.polyval(coeff, wl_full[h:l])
        data2c = data2 - correction
        new_df = new_df.append(pd.DataFrame(data2c).T)
    return new_df

def ss_polyfit(df, wl, order, j, k):
    h = 0
    l = len(df)
    wl = wl.iloc[h:l, 0]
    wvn1 = wl.iloc[h:j]
    wvn2 = wl.iloc[k:l]
    wvn = np.append(wvn1,wvn2)
    data2 = df.iloc[h:l]
    dataA = df.iloc[h:j]
    dataB = df.iloc[k:l]
    data3 = np.append(dataA,dataB)
    
    coeff = np.polyfit(wvn, data3, order)
    correction = np.polyval(coeff, wl[h:l])
    data2c = data2 - correction
    return data2c, correction


def fit_exponential_model(Peaks, time, with_plot, with_report, no_of_terms, **kwargs):
    mod = ExponentialModel()
    
    x = time
    y = Peaks

    def test_func(x, A1, A2, tau, **kwargs):
        if 'A3' in kwargs and 'tau2' in kwargs:
            return A2*np.exp(-x/tau) + kwargs.get('A2')*np.exp(-x/kwargs.get('tau2'))+A1
        else:
            return A2*np.exp(-x/tau) +A1

    if no_of_terms == 1: 
        mod = Model(test_func)
        out = mod.fit(y, x=x, A1=2221, tau=10*60, A2=7.5)
        params = [out.values['A1'], out.values['A2'], out.values['tau']]

    elif no_of_terms == 2: 
        print('terms = 2')
        mod = Model(test_func)
        out = mod.fit(y, x=x, A1=2221, tau=3*60, A2 = 1.5, A3=3.5, tau2=(21.06*60))
        params = [out.values['A1'], out.values['A2'], out.values['tau'], 
                  out.values['A3'], out.values['tau2']]

    if with_report == True:
        print(out.fit_report(min_correl=0.25))
        out.plot_residuals('k.')
        
    if with_plot == True:
        plt.figure()
        fig, ax = plt.subplots(figsize = (5,3))
        plt.plot(time/60, out.best_fit, 'r-', label='best fit')
        #plt.xlim(0, 90)

    return params

"""
Returns a dictionary of the amplitude, center, and sigma for the beginning and end of a specific identifier.
"""
def getUniversalVariables(identifier, SteadyStateSpectra):
    paramsDict = {}
    for key in SteadyStateSpectra[identifier].keys():
        spectrum = SteadyStateSpectra[identifier][key]
        y = np.zeros(len(spectrum))
        x = np.zeros(len(spectrum))
        for i in range(len(spectrum)):
            y[i] = spectrum.iloc[i,1]
            x[i] = spectrum.iloc[i,0]
        data = np.column_stack([x,y])
        dat = data
        x = dat[:, 0]
        y = dat[:, 1]

        mod = PseudoVoigtModel()
        pars = mod.guess(y, x=x)

        init = mod.eval(pars, x=x)
        out = mod.fit(y, pars, x=x)
        paramsDict[key] = out.best_values
    return paramsDict #returns a dictionary of amplitude, center, sigma


def GetStats(dictionary):
    List = []
    for key in dictionary:
        List.append(dictionary[key])
    mean = np.mean(List)
    std = np.std(List)
    return mean, std

#Data Handling Functions

"""
With user input, retrieves spectra to be used with Steady State functions.
"""
def getAllSteadyStateData():
    works = False
    while works == False: 
        try:
            num = eval(input('How many conditions are you working with?: '))
            works = True
        except NameError:
            print("Must be an integer")
            works = False
    SteadyStateSpectra = {}
    for i in range(num):
        conditionSpectraDict = {}
        condition = input('What is the condition #' + str(i+1) + '?: ')
        
        works = False
        while works == False: 
            try:
                numSpectra = eval(input('How many steady-state spectra are you entering for this condition?: '))
                works = True
            except NameError:
                print("Must be an integer")
                works = False
        for i in range(numSpectra):
            replicate = input('Enter name of replicate ' + str(i+1) + ': ')
            works = False
            while works == False: 
                try:
                    numSpectra2 = eval(input('How many spectra would you like to overlay for this replicate?: '))
                    works = True
                except NameError:
                    print("Must be an integer")
                    works = False
            for j in range(numSpectra2):

                works = False
                while works == False: 
                    try:
                        specLoc = input('Enter location of spectrum '+condition+' '+replicate+': ')
                        Spectrum = pd.read_table(specLoc, header=None)
                        conditionSpectraDict[j] = Spectrum
                        works = True
                    except FileNotFoundError:
                        print('Spectrum '+condition+' '+name+'not found. Try again.')
                        works = False
                
        SteadyStateSpectra[condition] = conditionSpectraDict
    return SteadyStateSpectra


"""
With user input, retrieves spectra to be used with Steady State functions.
"""
def getAllSteadyStateData_forSSonly():
    works = False
    while works == False: 
        try:
            num = eval(input('How many plots would you like to include? (up to 10): '))
            works = True
        except NameError:
            print("Must be an integer")
            works = False
    SteadyStateSpectra = {}
    for i in range(num):
        conditionSpectraDict = {}
        condition = input('What is the condition for plot #' + str(i+1) + '?: ')
        
        works = False
        while works == False: 
            try:
                numSpectra = eval(input('How many spectra would you like to overlay for plot #' + str(i+1) + '?: '))
                works = True
            except NameError:
                print("Must be an integer")
                works = False
        for j in range(numSpectra):
            replicate = input('Enter name of spectrum ' + str(j+1) + ': ')
            works = False
            while works == False: 
                try:
                    specLoc = input('Enter location of spectrum '+condition+' '+replicate+': ')
                    Spectrum = pd.read_table(specLoc, header=None)
                    conditionSpectraDict[replicate] = Spectrum
                    works = True
                except FileNotFoundError:
                    print('Spectrum '+condition+' '+replicate+'not found. Try again.')
                    works = False
                
        SteadyStateSpectra[condition] = conditionSpectraDict
    return SteadyStateSpectra


#SteadyState Functions
def plotSteadyState(steadyState, wlLow, wlHigh, pdf, title='Mean IR Band for Each Condition'):
    f = plt.figure(figsize=(7,20))
        
    index = 1
    for key in steadyState.keys():
        s = WLtoIndex(steadyState[key]['beg'], wlHigh)
        e = WLtoIndex(steadyState[key]['beg'], wlLow)
        ax = f.add_subplot(10,1,index)
        begData = steadyState[key]['beg']
        endData = steadyState[key]['end']
        begLabel = (str(key)+' Beg')
        endLabel = (str(key)+' End')
        ax.plot(begData.iloc[s:e,0], begData.iloc[s:e, 1]*1000, 'b', label=begLabel)
        ax.plot(endData.iloc[s:e,0], endData.iloc[s:e, 1]*1000, 'r', label=endLabel)
        ax.margins(x=0)
        ax.set_ylabel('Absorbance (m0D)')
        ax.set_xlabel('Wavenumber (cm\N{SUPERSCRIPT MINUS}\N{SUPERSCRIPT ONE})')
        #ax.spines['right'].set_visible(False)
        ax.tick_params(bottom='off', left='off')
        ax.legend(loc = 'best')
        index += 1
    f.subplots_adjust(hspace=0)
    wrapped_title = "\n".join(wrap(title, 40))
    f.suptitle(wrapped_title, fontsize=14, ha = 'right', y = 0.91)
    pdf.savefig(bbox_inches = "tight")

def plotSteadyStateOnly(steadyState, wlLow, wlHigh, pdf, title='Mean IR Band for Each Condition'):
    f = plt.figure(figsize=(7,20))

    index = 1
    for condition in steadyState.keys():
        ax = f.add_subplot(10,1,index)
        ax.margins(x=0)
        ax.set_ylabel('Absorbance (m0D)')
        ax.set_xlabel('Wavenumber (cm\N{SUPERSCRIPT MINUS}\N{SUPERSCRIPT ONE})')
        ax.tick_params(bottom='off', left='off')
        i = 0
        for replicate in steadyState[condition].keys():
            i += 1
            s = WLtoIndex(steadyState[condition][replicate], wlHigh)
            e = WLtoIndex(steadyState[condition][replicate], wlLow)
            Data = steadyState[condition][replicate]
            Label = (str(condition)+' '+str(replicate))
            ax.plot(Data.iloc[s:e,0], Data.iloc[s:e, 1]*1000, label=Label)
        index += 1
        ax.legend(loc = 'best')
    f.subplots_adjust(hspace=0)
    wrapped_title = "\n".join(wrap(title, 40))
    plt.title(wrapped_title, fontsize=14, y=1.02, loc = 'left')
    pdf.savefig(bbox_inches = "tight")


def main():
    ##Enter steady state data##
    print()
    print("\n" + "Now collecting steady state data...")
    
    steadyState = getAllSteadyStateData_forSSonly()

    ##Get dictionaries of amp, center, sigma at each condition##
    params = {}

    for key in steadyState.keys():
        params[key] = getUniversalVariables(key, steadyState)
        
    all_keys = ""
    for i in range(len(steadyState.keys())):
        if i < len(steadyState.keys()) - 1:
            all_keys = all_keys + str(list(steadyState.keys())[i])+ ", " 
        else:
            all_keys = all_keys + str(list(steadyState.keys())[i])
        

    print()

    #create pdf to save plots to
    pdfName = input('What do you want to name the pdf where plots will be saved?: ')
    pdf = PdfPages(pdfName+'.pdf')

    ##Plot##
    
    wl_start = float(input("Please input the lower bound for your steady-state plot x-axis "))
    wl_end = float(input("Please input the upper bound for your steady-state plot x-axis "))

    name = input('What would you like to title this plot? Type the word default, none, or input your own title. (default title = Mean IR Band for Each Condition) ')
    if name == 'default':
        plotSteadyStateOnly(steadyState, wl_start, wl_end, pdf)
    else:
        if name == 'none':
            name = ''
        plotSteadyStateOnly(steadyState, wl_start, wl_end, pdf, name)


    pdf.close()

if __name__ == '__main__':
    main()
