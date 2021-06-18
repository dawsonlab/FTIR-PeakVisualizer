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

def polyfit(SSorTR, df, wl, order, wlLow, wlHigh, wl_lowerbound, wl_upperbound):
    k = WLtoIndex(wl, wlLow)
    j = WLtoIndex(wl, wlHigh)
    l = WLtoIndex(wl, wl_lowerbound)
    h = WLtoIndex(wl, wl_upperbound)

    if SSorTR == 'SS':
        return ss_polyfit(df,wl,order,j,k, h, l)

def ss_polyfit(df, wl, order, j, k, h, l):

    wl = wl.iloc[:, 0]
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
    wl = wl.iloc[h:l]
    return data2c, correction, wl


def main():

    works = False
    while works == False:
        try:
            loc = input('Enter location of spectrum: ')
            works = True
        except FileNotFoundError:
            print('Data not found. Try again.')
            works = False
    
    full_data = pd.read_table(loc, header = None)
    wl = full_data
    data = full_data.iloc[:, 1]
    
    works = False
    while works == False: 
        try:
            order = eval(input('Enter order of polynomial to be subtracted from data: '))
            wl_lowerbound = float(input("Please input the lower bound for your data set "))
            wl_upperbound = float(input("Please input the upper bound for your data set "))
            wlLow = float(input("Please input the lower bound for your expected peak position "))
            wlHigh = float(input("Please input the upper bound for your expected peak position "))

            works = True
        except NameError:
            print("Must be an integer")
            works = False

    CorrectedData, correction, WL = polyfit('SS', data, wl, order, wlLow, wlHigh, wl_lowerbound, wl_upperbound)
    WL = pd.DataFrame(WL)
    CD = pd.DataFrame(CorrectedData)
    df = pd.concat([WL, CD], join = 'outer', axis = 1)
    name = loc.rsplit('.',1)[0]
    np.savetxt(name+'PolyCorrected.txt', df.values, fmt='%E', delimiter = '\t')

if __name__ == '__main__':
    main()
