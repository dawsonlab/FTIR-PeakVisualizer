# FTIR-PeakVisualizer

==========
PAPER TITLE: Direct Kinetic Observation of MAX1 Peptide Self-Assembly via Site-Specific Carbon-Deuterium Bonds as Infrared Probes

DATE: 2021-06-18

#################################################################
# If you used these scripts in your work, please cite:          #
#                                                               #
# Z.C. Adams, E.J. Olson, Z. Lian, A.Y. Kim, M. Holcomb,        #
# R. Adhikary, J. Zimmermann, and P.E. Dawson                   #
# Direct Observation of Peptide Hydrogel Self-Assembly          #
#                                                               #
#                                                               #
# DOI: tbd                                                      #
#                                                               #
#################################################################

HOW TO START THIS PROJECT: This project requires pandas, numpy, matplotlib, lmfit, and scipy. All necessary packages are compatible with pip install. For more information on the lmfit package please see https://lmfit.github.io/lmfit-py/.


USER INSTRUCTIONS: In this folder:

Jupyter Notebook Files:
IR Processing Notebook.ipynb
GeneralFunctions.ipynb
TimeResolvedFunctions.ipynb
SteadyStateFunctions.ipynb

Python Files:
IR_RUNME.py
IR_RUNME_SteadyState.py
PolyFitSS.py. 

Example Data Folders:
SteadyState
TimeResolved

To run the time-resolved spectra processing program from the command line, simply run IR_RUNME.py and follow the prompts. 

We have also provided a file specific for overlaying multiple conditions of steady-state data. To run this program, simply run IR_RUNME_SteadyState.py, and follow the prompts.

All figures will be saved to a pdf file named by the user. 

The user may encounter a situation where the .py files provided do not fit the use case. Although we tried to predict the most useful options, the functions provided are not limited to the options given in the main files. Jupyter notebook files have been provided for ease of access to these functions to fit a user's specific needs.

DATA TYPES: 

Import all IR data files into a folder within the program folder for ease of access.

IR_RUNME.py: Time Resolved data should contain three files for each condition, matrix of absorbance values, wavelength axis, time axis, named as follows: [NAME].csv, [NAME].x.txt, [NAME].y.txt, respectively. Time Resolved data should be Fourier-transformed, and saved as a .csv. If the user is running SF-FTIR, background can be subtracted from the data using the data collected pre-injection, defined when the code asks for Beginning and Beginning 2 values. If the user chooses to background correct Time Resolved data with a separately collected background spectrum, Beginning and Beginning 2 may be set to 0. Steady state data should consist of a single .txt file for each condition, which contains data that has already been background subtracted, polynomial corrected, and should be an average of all available steady state replicates. PolyFitSS.py has been provided to help with polynomial correction if necessary. 

IR_RUNME_SteadyState.py: Steady state data should consist of .txt files containing data that has already been background subtracted and polynomial corrected. PolyFitSS.py has been provided to help with polynomial correction if necessary. 

NOTE: When time-resolved data is imported, all three files will be necessary, so the name given should be "[NAME]" alone. However, when steady-state data is imported, the name given should be "[NAME].txt". 

EXAMPLES:

EXAMPLE1: For your convenience, sample data has been included. To run the time resolved analysis on the test data provided please follow the prompts using the answers shown below: 

python3 IR_RUNME.py 

Now collecting time-resolved data...
Note: each condition will require time resolved data and steady state data.
How many conditions are you working with?: 1
What is the condition #1?: test
How many time-resolved spectra are you entering for this condition?: 1
Enter name of replicate 1: a
Enter location of replicate 1: TimeResolved/TR_example


Now collecting steady state data...
For the condition: test
Enter location of beginning spectrum: SteadyState/SS_early.txt
Enter location of ending spectrum: SteadyState/SS_late.txt

Beginning indicates the scan number at which the “background” preinjection section stops. Background 2 indicates which scan number will be the first data point included in analysis.
Enter beginning: 21
Enter beginning 2: 21
Please give bounds for the region of interest
Enter lower wavenumber bound: 2180
Enter upper wavenumber bound: 2270
The average number determines the number of spectra included in a rolling average to improve signal to noise in the data.
Enter average number: 10
The poly number determines the order of the polynomial subtracted from the data to correct the baseline.
Enter poly number: 4

EXAMPLE2: For your convenience, sample data has been included. To run the steady state analysis on the test data provided please follow the prompts using the answers shown below: 

python3 IR_RUNME_SteadyState.py


Now collecting steady state data...
How many plots would you like to include? (up to 10): 1
What is the condtion for plot #1?: test
How many spectra would you like to overlay for plot #1?: 2
Enter name of spectrum 1: early
Enter location of spectrum test early: SteadyState/SS_early.txt
Enter name of spectrum 2: late 
Enter location of spectrum test late: SteadyState/SS_late.txt

What do you want to name the pdf where plots will be saved?: test

Please input the lower bound for your steady-state plot x-axis2180
Please input the upper bound for your steady-state plot x-axis2270
What would you like to title this plot? Type the word default, none, or input your own title. (default title = Mean IR Band for Each Condition)default



OVERVIEW OF CODE: The code ultimately performs three major functions: collects raw data, prepares said data, and then plots the prepared data in a variety of ways. This is done using a number of functions that are contained within the three files. Here is a walk-through of functions included in IR_RUNME.py.

DATA COLLECTION:
getTimeResolvedRawData() is run: This collects raw data in the form of three separate csv files (spectra matrix, time, and wavelength) for a user-specified number of conditions. Time's units are seconds, and wavelength's units are od. The function returns two dictionaries, allRawSpectraDict and allRawData. allRawSpectraDict, as the name suggests, contains only the spectra of the conditions. For example, the dictionary could look something like allRawSpectraDict = {'a':spectraA, 'b':spectraB, 'c':spectraC}, with 'a', 'b', and 'c' being the user-entered names of the conditions and spectraA, spectraB, and spectraC being the spectra for each condition. On the other hand, allRawData contains the spectrum, time, and wavelength for each condition. This dictionary is structured as such: allRawData = {'a':{'time':t, 'wl':wavelength, 'spectrum':spectrumA}, 'b':{'time':t, 'wl':wavelength, 'spectrum':spectrumB}, 'c':{'time':t, 'wl':wavelength, 'spectrum':spectrumC}). Thus, allRawData is a dictionary that contains an additional dictionary in it for each condition. 

getAllSteadyStateDataforTR() is run: This function collects raw Steady State data (only spectra). For each condition specified, a beginning and ending spectrum is collected and organized into a dictionary. At the end, the dictionary SteadyStateSpectra is returned, with a structure as such: SteadyStateSpectra = {'a':{'beg':begSpectra, 'end':endSpectra}, 'b':{'beg':begSpectra, 'end':endSpectra}, 'c':{'beg':begSpectra, 'end':endSpectra}}. 

NOTE: The names of the conditions will be consistent when entering both the time resolved data and the steady state data. For example, the user cannot enter "1" as the condition's name for the time resolved data and "condition1" for the condition's name for the steady state data. The condition must be inputted as "1" or "condition1" in both scenarios. 

getUniversalVariables(identifier, SteadyStateSpectra) is run: This function returns a dictionary of the amplitude, center, and sigma for the beginning and end of a specific condition send in through the parameter "identifier". In order to find the amplitude, center, and sigma, each spectrum is fit to the Pseudo Voigt Model from the lmfit package. lmfit returns a dictionary of the best-fit amplitude, center, and sigma using best_values. getUniversalVariables returns the dictionary paramsDict that contains the best-fit variables for each beginning and ending spectrum of each condition. paramsDict's structure looks like this: paramsDict = {'a':{'beg':{'amplitude':amplitude, 'center':center, 'sigma':sigma}, 'end':{'amplitude':amplitude, 'center':center, 'sigma':sigma}}, 'b':{'beg':{'amplitude':amplitude, 'center':center, 'sigma':sigma}, 'end':{'amplitude':amplitude, 'center':center, 'sigma':sigma}}, 'c':{'beg':{'amplitude':amplitude, 'center':center, 'sigma':sigma}, 'end':{'amplitude':amplitude, 'center':center, 'sigma':sigma}}}. 

getTRPVariables() is run: This function collects a number of variables from the user. These variables include beginning, beginning2, left-bound, right-bound, average number, and poly number. This information is all returned in the dictionary TRP. Beginning indicates the scan number at which the “background” pre-injection section stops. Background 2 indicates which scan number will be the first data point included in analysis.The left- and right- bounds are bounds for the region of interest. The average number determines the number of spectra included in a rolling average to improve signal to noise in the data. The poly number determines the order of the polynomial subtracted from the data to correct the baseline.

DATA PREPARATION:
prepareAllTimeResolvedRawData(allRawData, spectraDict, paramsDict, TRP) is run: This function prepares all of the Time Resolved data that we have previously retrieved through the user. For each condition, the parameters found through getUniversalVariables(), spectra contained in allRawSpectraDict, and times and wavelengths contained in allRawData are prepared. The preparation is done by the function Dictionary_of_Prepped_Spectra(), which is found in the file TimeResolvedFunctions. The region of interest for the peak location is determined from provided Steady State spectra, using the average of the beginning and ending peak positions +/- the average fwhm. Each spectrum is then background subtracted, averaged with a rolling average across the time axis, and polynomial subtracted to flatten the baseline. 

PLOTTING:
Now that the data is all prepared, plotting can begin. First, a pdf file is created to save plots to. This pdf is created using PdfPages from matplotlib.backends.backend_pdf. The user will be asked what they would like to name said pdf, and the pdf will automatically be saved to the same folder from which the user is running this program. 

Next, the user gets to choose whether to plot all of the peaks over time or a specific peak over time. If the user decides to plot all of their data, plotPeaksOverTime(All_Prepped_Spectra, AllTimeAndWL, pdf) is run. If the user decides to plot specific peaks, plotPeakShiftAverage(All_Prepped_Spectra, AllTimeAndWL, identifier, title, title2, pdf) is run. Note that the main difference between plotPeaksOverTime (for all data) and plotPeakShiftAverage (for specific data) is that plotPeakShiftAverage takes in an identifier, which will specify which dataset to plot. plotPeakShiftAverage will also ask the user for titles. These plots will all be saved to the named pdf created above.

The next plots produced overlay the unaveraged initial and final spectra for either every condition or for specific conditions. The user is able to title each plot. 

The next function run is combined_time_plots(All_Prepped_Spectra, allTimeAndWL, identifier, description, paramsDict, pdf). This function both plots and prints a report containing important information, most notably the Mean Time Constant. The user is able to make the choice between running combined_time_plots() for all of the spectra or for specified identifiers. The user is also able to choose whether or not they want to fit the spectra to a Single Gaussian plot or to a Half Integral plot. If the user specifies "SingleGaussian," PeakFromSingleGaussian(identifier, All_Prepped_Spectra, allTimeAndWL, paramsDict) is run, and if the user specifies "HalfIntegral," PeakFromHalfIntegral(identifier, All_Prepped_Spectra, allTimeAndWL, TRP) is run. Within each method, an average spectrum is found using average_of_prepped_spectra(dictionary, identifier). average_of_prepped_spectra takes all of the various spectra entered for each condition and produces an average matrix. Then, this average matrix is plotted, and a time constant is produced. This time constant is printed in each plot and is also printed in the Terminal. 

Finally, the next function plots out the Steady State data, producing a plot of the mean Carbon-Deuterium band for each condition. 

After the user finishes plotting, the pdf is closed, and the user can now access their plots. 

Here is a run through of functions included in IR_RUNME_SteadyState.py:

DATA COLLECTION:

getAllSteadyStateData_forSSonly() is run: This function collects raw Steady State data (only spectra). For each condition specified (up to 10), a series of spectra to be overlaid are collected and organized into a dictionary. At the end, the dictionary SteadyStateSpectra is returned, with a structure as such: SteadyStateSpectra = {'a':{'1':spectrum, '2':spectrum, '3':spectrum}, 'b':{'1':spectrum, '2':spectrum, '3':spectrum}}. 

plotSteadyStateData_forSSonly is run: This function takes the data imported by getAllSteadyStateData_forSSonly() and makes up to 10 consecutive plots, each with the defined number of overlaid spectra. 

Here is a run through of functions included in PolyFitSS.py:

Data is imported as a .txt file with a wavenumber column and an absorbance column. Bounds for the region of interest and the expected location of the peak are specified by the user. The function polyfit() is run on the data provided, which removes the region the peak is expected to appear in, fits a polynomial to the remaining data, then subtracts the polynomial from the original data, outputting the corrected data, the correction, and the x-axis. The corrected data and x-axis values are then saved into a .txt file with the phrase "PolyCorrected" appended to the original name of the data set. 


**************************************************************************
**  CITING THESE TOOLS IN PUBLICATIONS
**************************************************************************

We ask that you acknowledge using these tools in any publications arising from 
their use.

**************************************************************************
** LICENSES
**************************************************************************

You can refer to the LICENSES.txt file included in the distribution

**************************************************************************
** HOW TO CONTACT US:
**************************************************************************

zadams@scripps.edu
dawson@scripps.edu
