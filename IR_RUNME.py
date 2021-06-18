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
With user input, retrieves data to be used with Time Resolved functions.
"""
def getTimeResolvedRawData():
    works = False
    while works == False: 
        try:
            num = eval(input('How many conditions are you working with?: '))
            works = True
        except NameError:
            print("Must be an integer")
            works = False
            
   
    allRawSpectraDict = {}
    allRawData = {}

    for i in range(num):
        conditionSpectraDict = {}
        conditionDataDict = {}
        condition = input('What is the condition #' + str(i+1) + '?: ')
        works = False
        while works == False: 
            try:
                numSpectra = eval(input('How many time-resolved spectra are you entering for this condition?: '))
                works = True
            except NameError:
                print("Must be an integer")
                works = False
        for i in range(numSpectra):
            name = input('Enter name of replicate ' + str(i+1) + ': ')
            works = False
            while works == False:
                try:
                    loc = input('Enter location of replicate ' + str(i+1) + ': ')
                    spectrum, time, wl = ReadIn(loc)
                    data = {'spectrum':spectrum, 'time':time, 'wl':wl}
                    conditionSpectraDict[name] = spectrum
                    conditionDataDict[name] = data
                    works = True
                except FileNotFoundError:
                    print('Data not found. Try again.')
                    works = False
        allRawSpectraDict[condition] = conditionSpectraDict
        allRawData[condition] = conditionDataDict
    return allRawSpectraDict, allRawData


"""
With user input, retrieves spectra to be used with Steady State functions.
"""
def getAllSteadyStateDataForTR(dictionary):
    num = len(dictionary.keys())

    SteadyStateSpectra = {}
    for i in range(num):
        conditionSpectraDict = {}

        condition = list(dictionary.keys())[i]
        print('For the condition: '+str(condition))

        works = False
        while works == False: 
            try:
                begLoc = input('Enter location of beginning spectrum: ')
                begSpectrum = pd.read_table(begLoc, header=None)
                conditionSpectraDict['beg'] = begSpectrum
                works = True
            except FileNotFoundError:
                print('Beginning spectrum not found. Try again.')
                works = False
        works = False
        while works == False: 
            try:
                endLoc = input('Enter location of ending spectrum: ')
                endSpectrum = pd.read_table(endLoc, header=None)
                conditionSpectraDict['end'] = endSpectrum
                works = True
            except FileNotFoundError:
                print('Ending spectrum not found. Try again.')
                works = False
                
        SteadyStateSpectra[condition] = conditionSpectraDict
    return SteadyStateSpectra

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
                begLoc = input('Enter location of beginning spectrum: ')
                begSpectrum = pd.read_table(begLoc, header=None)
                conditionSpectraDict['beg'] = begSpectrum
                works = True
            except FileNotFoundError:
                print('Beginning spectrum not found. Try again.')
                works = False
        works = False
        while works == False: 
            try:
                endLoc = input('Enter location of ending spectrum: ')
                endSpectrum = pd.read_table(endLoc, header=None)
                conditionSpectraDict['end'] = endSpectrum
                works = True
            except FileNotFoundError:
                print('Ending spectrum not found. Try again.')
                works = False
                
        SteadyStateSpectra[condition] = conditionSpectraDict
    return SteadyStateSpectra


def getTRPVariables(allRawData):
    print("Beginning indicates the scan number at which the “background” preinjection section stops. Background 2 indicates which scan number will be the first data point included in analysis.")
    b = eval(input('Enter beginning: '))
    b2 = eval(input('Enter beginning 2: '))
    print("Please give bounds for the region of interest")
    left_bound = eval(input('Enter lower wavenumber bound: '))
    right_bound = eval(input('Enter upper wavenumber bound: '))
    print("The average number determines the number of spectra included in a rolling average to improve signal to noise in the data.")
    ave_n = eval(input('Enter average number: '))
    print("The poly number determines the order of the polynomial subtracted from the data to correct the baseline.")
    poly_n = eval(input('Enter poly number: '))
    TRP = {'b':b,'b2':b2,'left_bound':left_bound,'right_bound':right_bound,'ave_n':ave_n,'poly_n':poly_n}
    return TRP


"""
Returns a dictionary of prepped spectra and a separate dictionary of respective times and wavelengths.
"""
def prepareAllTimeResolvedRawData(allRawData, allRawSpectraDict, paramsDict, TRP):
    All_Prepped_Spectra = {}
    AllTimeAndWL = {}
    for key in allRawSpectraDict.keys():
        params = paramsDict[str(key)]
        preppedSpectra, timeAndWL, TRP = Dictionary_of_Prepped_Spectra(key, allRawSpectraDict, allRawData, TRP, params)
        All_Prepped_Spectra[key] = preppedSpectra
        AllTimeAndWL[key] = timeAndWL
    return All_Prepped_Spectra, AllTimeAndWL, TRP


"""
Retrieves the time and wavelength for a specific identifier.
"""
def getTimeAndWL(identifier, AllTimeAndWL):
    time = AllTimeAndWL[identifier]['time']
    wl = AllTimeAndWL[identifier]['wl']
    return time, wl


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

#TimeResolved Functions
def ReadIn (name):
    table = pd.read_csv(name + '.csv', header = None)
    time = pd.read_table(name + '.z.txt', header = None)
    wavelength = pd.read_table(name + '.x.txt', header = None)
    return (table, time, wavelength)

def Dictionary_of_Prepped_Spectra(identifier, allSpectraDict, allDataDict, TRP, paramsDict):
    Prepped_Spectra = {}
    timeAndWLDict ={}
    sumCenters = paramsDict['beg']['center']+paramsDict['end']['center']
    sumSigmas = paramsDict['beg']['sigma']+paramsDict['end']['sigma']

    peak_beg = sumCenters/2-(sumSigmas/2)*2.3548200
    peak_end = sumCenters/2+(sumSigmas/2)*2.3548200
    spectraDict = allSpectraDict[identifier]
    TRP[identifier] = {'peak_beg':peak_beg,'peak_end':peak_end}
    for key in spectraDict.keys():
        spectrum = spectraDict[key]
        time = allDataDict[identifier][key]['time']
        wl = allDataDict[identifier][key]['wl']
        Out_spectrum, time, wl = PrepareSpectra(spectrum, time, wl, 
                                                    TRP['b'],TRP['b2'], 
                                                    TRP['left_bound'], TRP['right_bound'], 
                                                    TRP['ave_n'], TRP['poly_n'], 
                                                    peak_beg, peak_end)

        Prepped_Spectra[key] = Out_spectrum
        timeAndWLDict = {'time':time, 'wl':wl}
    return Prepped_Spectra, timeAndWLDict, TRP

"""
Use this function to take a dictionary of prepped spectra (already background and polynomial corrected) and get an average spectrum
"""
def average_of_prepped_spectra(dictionary, identifier):
    spectra_keys = list(dictionary[identifier].keys())
    Sum = dictionary[identifier][spectra_keys[0]]-dictionary[identifier][spectra_keys[0]]
    for key in spectra_keys:
        Sum += dictionary[identifier][key]
    table_mean = Sum/len(spectra_keys)
    return table_mean

"""
This function takes a Time Resolved spectrum given as a dataframe (pd.read_csv), averages over a rolling average specified by n, then subtracts the mean of the "before injection" period, fits and subtracts a polynomial and outputs the converted spectra and relevant time file. Reasonable way to call this function is as follows: PrepareSpectra(test, wllDF, 21,21, 104, 256,100, 4)
"""
def PrepareSpectra (TRspec, time, wl, before, beginning, endpoint_wl2, endpoint_wl1,ave_n, poly_n, wlLow, wlHigh):
    beg = beginning
    end1 = WLtoIndex(wl, endpoint_wl1)
    end2 = WLtoIndex(wl, endpoint_wl2)
    S = averaging(TRspec.iloc[beg:, end1:end2], ave_n) - np.mean(TRspec.iloc[0:before, end1:end2].values, axis = 0)
    new_table = polyfit('TR', S, wl.iloc[end1:end2,:], poly_n, wlLow, wlHigh)
    
    AdjTime = averaging(time.iloc[0:len(time)-beg], ave_n)
    TimeArray = AdjTime.values.T[0]
    return (new_table, TimeArray, wl.iloc[end1:end2,:])

"""
This function takes a Time Resolved spectrum given as a dataframe (pd.read_csv) and averages over a rolling average specified by ave_n.
"""
def averaging(df,ave_n, **kwargs):  
    new_df = pd.DataFrame()
    reps = len(df)
    start = ave_n
    finish = len(df)+1
    for j in np.arange(start, finish):
        sector = df[(j-(ave_n)):j]
        aves = sector.mean(axis = 0)
        aves_T = pd.DataFrame(aves).T
        new_df = new_df.append(aves_T)
        
    return new_df

def PeakFromSingleGaussian(identifier, All_Prepped_Spectra, tAndWL, paramsDict):
    def _1gaussian(x, amp1,cen1,sigma1):
        return amp1*(1/(sigma1*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x-cen1)/sigma1)**2)))
    time, wl = getTimeAndWL(identifier, tAndWL)
    peaks = []
    ts = []
    begParams = paramsDict[identifier]['beg']
    endParams = paramsDict[identifier]['end']
    begAmp = begParams['amplitude']
    endAmp = endParams['amplitude']
    begSigma = begParams['sigma']
    endSigma = endParams['sigma']
    begCenter = begParams['center']
    endCenter = endParams['center']

    amp = np.mean([begAmp, endAmp])
    sigma = np.mean([begSigma, endSigma])
    center = np.mean([begCenter, endCenter])

    spectrum = average_of_prepped_spectra(All_Prepped_Spectra, identifier)

    for t in np.arange(0,np.shape(spectrum)[0],1):
        y = np.zeros(len(spectrum.iloc[0]))
        x = np.zeros(len(wl))
        for i in range(len(spectrum.iloc[0])):
            y[i] = spectrum.iloc[t, i]
            x[i] = wl.iloc[i]
        data = np.column_stack([x,y])

        dat = data
        x = dat[:, 0]
        y = dat[:, 1]

        popt_gauss, pcov_gauss = scipy.optimize.curve_fit(_1gaussian, x, y, p0=[amp, center, sigma])

        perr_gauss = np.sqrt(np.diag(pcov_gauss))
        peaks = np.append(peaks, popt_gauss[1])
        ts = np.append(ts, time[t])
    return(peaks, ts)

def PeakFromHalfIntegral(identifier, All_Prepped_Spectra, tAndWL, TRP):
    import warnings
    warnings.filterwarnings("ignore")
    
    spectrum = average_of_prepped_spectra(All_Prepped_Spectra, identifier)
    time, wl = getTimeAndWL(identifier, tAndWL)
    wlLow = TRP[identifier]['peak_beg']
    wlHigh = TRP[identifier]['peak_end']
    k = WLtoIndex(wl, wlLow)
    j = WLtoIndex(wl, wlHigh)
    wl_testSize = wl.iloc[j:k]
    
    peaks = []
    ts = []
    for t in np.arange(0,np.shape(spectrum)[0],1):

        time_slice = spectrum.iloc[t, j:k]
        total = np.trapz(time_slice)
        percents = []
        for i in range(len(time_slice)):
            integral = np.trapz(time_slice[0:i])/total
            percents=np.append(percents, integral)
        def test_func(w, m, b,c):
            return m/(b+(np.exp(-c*w)))
        params, params_covariance = scipy.optimize.curve_fit(test_func, np.arange(k-j), percents, p0=[-1, 0.5, 2])
        peakpos = (-1/params[2])*np.log((params[0]-0.5*params[1])/0.5) 
        peak = wlHigh - peakpos*(wlHigh-wlLow)/(k-j)
        peaks = np.append(peaks, peak)
        ts = np.append(ts, time[t])
    return peaks, ts


def findPeakWithSingleGaussianFit(All_Prepped_Spectra, tAndWL, identifier, description, paramsDict, pdf):
    peaks, time = PeakFromSingleGaussian(identifier, All_Prepped_Spectra, tAndWL, paramsDict)
    plot_any_peaks_over_time(time, peaks, description, with_plot = True, with_report = True, to_save = True)
    plt.title('$\\tau_g$ = 9.6 $\pm$ 3.4 min  ' , loc = 'Right', y = 0.85, fontsize = 14)
    pdf.savefig()



def findPeaksAtPointofHalfTotalIntegration(identifier, All_Prepped_Spectra, tAndWL, TRP, description, pdf):
    peaks, time = PeakFromHalfIntegral(identifier, All_Prepped_Spectra, tAndWL, TRP)
    plot_any_peaks_over_time(time, peaks, description, pdf, to_save = False,  with_plot = True, with_report = True)


    ##TimeResolved Plotting Functions###

def plot_shift(spectrum, wl, n, **kwargs):
    for i in np.arange(0,len(spectrum),n):
        if 'color' in kwargs:
            colors = kwargs.get('color')
        else:
            colors = plt.cm.RdYlBu(round(1 - i/700, 2))
        plt.plot(wl, spectrum.iloc[i,:].T, zorder=len(spectrum)-i, color=colors)
        plt.title(kwargs.get('title'))
        plt.xlabel(kwargs.get('xlabel'))
        plt.ylabel(kwargs.get('ylabel'))

"""
Takes a spectrum from PrepareSpectra (with generated time and wavelengths) and plots it with colors representing change over time
"""
def ShiftOverTime(spectrum, time, wl, string, Title, Title2):
    time = time/60
    
    fig, ax = plt.subplots(figsize = (10,3))
    
    plot_shift(spectrum*1000,wl, 10 )

    wrapped_title = "\n".join(wrap(Title, 40))
    plt.title(wrapped_title, fontsize=14, y = 1.02, loc = 'left')
    fig.suptitle(Title2, fontsize=14, y = 0.85, x = 0.20)
    
    plt.xlabel('Wavenumber (cm\N{SUPERSCRIPT MINUS}\N{SUPERSCRIPT ONE})')
    plt.ylabel('Absorbance (mOD)')

    N = 1000
    cmap = plt.get_cmap('RdYlBu_r',N)
    norm = mpl.colors.Normalize(vmin=min(time),vmax=max(time))
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ticks=np.arange(0,max(time), 10), boundaries=np.arange(-0.05,max(time),.0005))
    cbar.set_label('Time (min)')

def plot_any_peaks_over_time(time, peaks, description, pdf, **kwargs):
    with_plot = kwargs.get("with_plot")
    with_report = kwargs.get("with_report")
    to_save = kwargs.get("to_save")
    params = fit_exponential_model(peaks, time, with_plot, with_report, 1)
    if to_save == True: 
        time_label = '$\\tau_g$ = ' +str(round(params[2]/60,2)) +'min'
    else:
        time_label = None
        

    plt.plot(time/60, peaks, label= time_label)
    plt.xlabel('time (min)', fontsize = 14)
    plt.ylabel('wavenumber (cm\N{SUPERSCRIPT MINUS}\N{SUPERSCRIPT ONE})', fontsize = 14)
    #plt.ylim(2220,2229.5)
    
    if to_save == False:
        plt.legend(loc = 'best')
        wrapped_title = "\n".join(wrap(description, 40))
        plt.title(wrapped_title, fontsize=14, y = 1.02, loc = 'left')
        txt = 'TimeConstant: '+str(round(params[2]/60,2)) +'minutes'
        plt.text(100, 5, txt, ha='center', va='center', transform=None)

        pdf.savefig(bbox_inches = "tight")

def plotPeaksOverTime(All_Prepped_Spectra, AllTimeAndWL, pdf):
    for key in All_Prepped_Spectra.keys():
        ave = average_of_prepped_spectra(All_Prepped_Spectra, key)
        time, wl = getTimeAndWL(key, AllTimeAndWL)
        spectrum = ave
        time = time/60
        fig, ax = plt.subplots(figsize = (10,3))

        plot_shift(spectrum*1000,wl, 10 )
        ax.margins(x=0)

        ax.tick_params(axis="x", labelsize=14)
        ax.tick_params(axis="y", labelsize=14)
        #ax.title(title, fontsize=20, y = 1.02, loc = 'left')
        fig.suptitle(key, fontsize=14, y = 0.85, x = 0.20)

        plt.xlabel('wavenumber (cm\N{SUPERSCRIPT MINUS}\N{SUPERSCRIPT ONE})', fontsize = 14)
        plt.ylabel('absorbance (10\N{SUPERSCRIPT MINUS}\N{SUPERSCRIPT THREE})', fontsize = 14)

        N = 1000
        cmap = plt.get_cmap('RdYlBu_r',N)
        norm = mpl.colors.Normalize(vmin=min(time),vmax=max(time))
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = plt.colorbar(sm, ticks=np.arange(0,max(time), 10), boundaries=np.arange(-0.05,max(time),.0005), pad = 0.01);
        cbar.set_label('time (min)', fontsize = 14)
        #plt.savefig(key+'PeakShiftAverage.png', bbox_inches='tight', pad_inches=0.1)
        #plt.show()
        pdf.savefig(bbox_inches='tight', pad_inches=1.0)


def plotPeakShiftAverage(All_Prepped_Spectra, AllTimeAndWL, identifier, title, title2, pdf):
    ave = average_of_prepped_spectra(All_Prepped_Spectra, identifier)
    time, wl = getTimeAndWL(identifier, AllTimeAndWL)
    ShiftOverTime(ave, time, wl, identifier, title, title2)
    #plt.savefig(identifier+'PeakShiftAverage.png', bbox_inches='tight', pad_inches=0.1)
    #plt.show()
    pdf.savefig(bbox_inches='tight', pad_inches=1.0)


def plotUnaveragedInitFin(All_Prepped_Spectra, AllTimeAndWL, identifier, description, pdf):
    time, wl = getTimeAndWL(identifier, AllTimeAndWL)    
    
    for individual_spectrum in All_Prepped_Spectra[identifier].keys():
        if description == 'default':
            title = 'Initial and Final Spectrum From Time Resolved Data Set ' + identifier + ' '+ individual_spectrum
        elif description == 'none':
            title = ''
        else:
            title = description
        spectraDict = All_Prepped_Spectra[identifier]
        plt.figure()
        plt.plot(wl, spectraDict[individual_spectrum].iloc[0, :], label = 'initial spectrum')
        plt.plot(wl, spectraDict[individual_spectrum].iloc[len(spectraDict[individual_spectrum])-1, :], label = 'final spectrum')
        plt.legend(loc = 'best')


        wrapped_title = "\n".join(wrap(title, 40))
        plt.title(wrapped_title, fontsize=14, y=1.02, loc = 'left')
        
        pdf.savefig(bbox_inches='tight', pad_inches=0.1)

def combined_time_plots(All_Prepped_Spectra, tAndWL, paramsDict, TRP, identifier, name, method, pdf):
    time, wl = getTimeAndWL(identifier, tAndWL)
    time_constants = dict()
    if method == 'HalfIntegral':
        peaks, time = PeakFromHalfIntegral(identifier, All_Prepped_Spectra, tAndWL, TRP)
    if method == 'SingleGaussian':
        peaks, time = PeakFromSingleGaussian(identifier, All_Prepped_Spectra, tAndWL, paramsDict)
    plot_any_peaks_over_time(time, peaks, name, pdf, to_save = False,  with_plot = True, with_report = True)
    time_constants[identifier] = fit_exponential_model(peaks, time, 0,0,1)[2]/60
    stats = GetStats(time_constants)
    print('Mean Time Constant: '+ str(np.round(stats[0],2)), u"\u00B1"+str(np.round(stats[1],2))+' minutes')
    return time_constants





def main():
    ##Enter raw time-resolved data##
    print("\n" + "Now collecting time-resolved data...")
    print("Note: each condition will require time resolved data and steady state data.")
    allRawSpectraDict, allRawData = getTimeResolvedRawData()
    ##Enter steady state data##
    print()
    print("\n" + "Now collecting steady state data...")
    
    steadyState = getAllSteadyStateDataForTR(allRawSpectraDict)
    #steadyState = getAllSteadyStateData()

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

    #Get TRP variables
    TRP = getTRPVariables(allRawData)

    print("\nProcessing data and preparing to plot...\n")
    
    ##Prepare data for plotting##
    preppedSpectra, tAndWL, TRP = prepareAllTimeResolvedRawData(allRawData, allRawSpectraDict, params, TRP)

    #create pdf to save plots to
    pdfName = input('What do you want to name the pdf where plots will be saved?: ')
    pdf = PdfPages(pdfName+'.pdf')

    ##Plot##

    #PLOT 1
    numPlots = approveChoice("\nDo you want to plot the shifting peak averaged over all replicates for every condition or for specific conditions? (every/specific/none): ", ["every", "specific", "none"])

    if numPlots == "every":
        plotPeaksOverTime(preppedSpectra, tAndWL, pdf)
        
    elif numPlots == 'specific':
        done = False
        while done == False:
            print("As a reminder, here are the identifiers you entered for each condition: " + all_keys)
            identifier = approveChoice('What condition do you want to plot?: ', all_keys)
            
            title = input('What would you like to title the plot for condition ' + identifier + '? Type the word default, none, or input your own title. (default title = Peak Shifting over Time): ')
            if title == 'default':
                title = 'Peak Shifting over Time'
            if title == 'none':
                title = ''
            
            title2 = input('What would you like as a subtitle the plot for condition ' + identifier + '? Type the word default, none, or input your own title. (default subtitle = '+ identifier+ '): ')
            if title2 == 'default':
                title2 = identifier
            if title2 == 'none':
                title2 = ''
            
            plotPeakShiftAverage(preppedSpectra, tAndWL, identifier, title, title2, pdf)
            keepPlotting = approveChoice('\nDo you want to make another plot of this type? (yes/no): ', ["yes", "no"])
            if keepPlotting == 'yes':
                done = False
            if keepPlotting == 'no':
                done = True


    #PLOT 2

    numUnaveragedPlots = approveChoice('\nDo you want to overlay the unaveraged initial and final spectra for every condition or for specific conditions? (every/specific/none): ', ["every", "specific", "none"])
    
    if numUnaveragedPlots == "every":
        titleDefault = approveChoice('Would you like to use the default title? (default title = Peak Position over Time from Method for Condition) (yes/no): ', ["yes", "no"])
        
        for identifier in preppedSpectra.keys():
            if titleDefault == "yes":
                title = 'default'
            else:
                title = input('What would you like to title the plot for condition ' + identifier + '? Type the word default, none, or input your own title. (default title = Initial and Final Spectrum From Time Resolved Data Set): ')
            plotUnaveragedInitFin(preppedSpectra, tAndWL, identifier, title, pdf)
            
    elif numUnaveragedPlots == "specific":
        done = False
        while done == False:
            print('As a reminder, here are the identifiers you entered for each condition: ' + all_keys)
            identifier = approveChoice('What condition do you want to plot?: ', all_keys)
            title = input('What would you like to title the plot for condition ' + identifier + '? Type the word default, none, or input your own title. (default title = Initial and Final Spectrum From Time Resolved Data Set): ')
            plotUnaveragedInitFin(preppedSpectra, tAndWL, identifier, title, pdf)
            
            keepPlotting = approveChoice('\nDo you want to make another plot of this type? (yes/no): ', ["yes", "no"])
            if keepPlotting == 'yes':
                done = False
            if keepPlotting == 'no':
                done = True
        

    #PLOT 3

    numCombinedTimePlots = approveChoice('\nThe next set of plots will plot the peak position over time from data averaged over all replicates. Would you like to plot every condition or specific conditions? (every/specific/none): ', ["every", "specific", "none"])
    
    if numCombinedTimePlots == 'every':
        method = approveChoice('Which method would you like to use? (HalfIntegral/SingleGaussian): ', ["HalfIntegral", "SingleGaussian"])
        nameDefault = approveChoice('Would you like to use the default title? (default title = Peak Position over Time from Method for Condition) (yes/no): ', ["yes", "no"])

        for identifier in preppedSpectra.keys():
            if nameDefault == "yes":
                name = 'Peak Position over Time from '+method+' for Condition '+identifier
            else:
                name = input('What would you like to title the plot for condition ' + identifier + '? Type the word default, none, or input your own title. (default title = Peak Position over Time from '+method+' for Condition '+identifier+'): ')
                if name == 'default':
                     name = 'Peak Position over Time from '+method+' for Condition '+identifier
                if name == 'none':
                    name = ''
            print("\nCalculating peak positions...\n")
            combined_time_plots(preppedSpectra, tAndWL, params, TRP, identifier, name, method, pdf)
            
    elif numCombinedTimePlots == 'specific':
        done = False
        while done == False:
            print('As a reminder, here are the identifiers you entered for each condition: ' + all_keys)
            identifier = approveChoice('What condition do you want to plot?: ', all_keys)
            method = approveChoice('Which method would you like to use? (HalfIntegral/SingleGaussian): ', ["HalfIntegral", "SingleGaussian"])
            
            name = input('What would you like to title the plot for condition ' + identifier + '? Type the word default, none, or input your own title. (default title = Peak Position over Time from '+method+' for Condition '+identifier+'): ')
            if name == 'default':
                 name = 'Peak Position over Time from '+method+' for Condition '+identifier
            if name == 'none':
                name = ''
                
            print("\nCalculating peak positions...\n")
            combined_time_plots(preppedSpectra, tAndWL, params, TRP, identifier, name, method, pdf)

            keepPlotting = approveChoice('\nDo you want to make another plot of this type? (yes/no): ', ["yes", "no"])
            if keepPlotting == 'yes':
                done = False
            if keepPlotting == 'no':
                done = True


    #PLOT 4

    SteadyStatePlotting = approveChoice('\nWould you like to plot the Steady State data? This step produces a plot of the mean IR band for each condition. (yes/no)', ["yes", "no"])
    wl_start = TRP['left_bound']#float(input("Please input the lower bound for your steady-state plot x-axis"))
    wl_end = TRP['right_bound']#float(input("Please input the upper bound for your steady-state plot x-axis"))
    if SteadyStatePlotting == 'yes':
        name = input('What would you like to title this plot? Type the word default, none, or input your own title. (default title = Mean IR Band for Each Condition) ')
        if name == 'default':
            plotSteadyState(steadyState, wl_start, wl_end, pdf)
        else:
            if name == 'none':
                name = ''
            plotSteadyState(steadyState, wl_start, wl_end, pdf, name)

    pdf.close()

if __name__ == '__main__':
    main()
