#######
####### THIS PROGRAM HAS ALL ANALYSIS REALISED INTO MY THESIS
####### ----> Lateral distribution
####### ----> Energy distribution
####### ----> All energy deposited for the showers
####### -------------------ETC ETC -------------------

####* PARTICLES WILL BE TREATED AS OBJECT OF A CLASS CALLED "PARTICLES"

#IMPORTING LIBRARIES TO USING

import pandas as pd
import numpy as np
import math 
import scipy as sp
import fnmatch
from natsort import natsorted
import re 
import os

#GENERAL VARIABLES

eR = np.concatenate((np.arange(0.1, 10, 0.3),np.power(10, np.arange(1, 5+0.025, 0.025)))) ##energy range, this variable can be selected using any method of python or python modules
seq = "Data*.txt" ##Secuence of the name, this variable is an string to represent the names of the output files of readPartExample program
NSHOWS = 1000

#GENERAL FUNCTIONS

def filesFilter(height, typeP, seq):
    allFiles = os.listdir("./{}/{}".format(height, typeP))
    dFilesSelect = np.array([])
    for i in allFiles:
        if fnmatch.fnmatch(i, seq):
            dFilesSelect = np.append(dFilesSelect, i)
    return np.array(natsorted(dFilesSelect, key = lambda y: y.lower()))

def closeIt(listNumbers, value): #* pass it any value and it returns the more approximate value of the list 
    return(min(listNumbers, key = lambda y : abs(y - value)))

def selectEnergy(listNumbers, valueMin, valueMax):
    energiesSelected = np.array([])
    for i in listNumbers:
        if i >= closeIt(listNumbers, valueMin) and i <= closeIt(listNumbers, valueMax):
            energiesSelected = np.append(energiesSelected, i)
    return energiesSelected

def specificEnergy(height, energyS, typeP, seq):
    filterFilesV = filesFilter(height, typeP, seq)
    energySelected = None
    for i in filterFilesV:
        if np.where(filterFilesV == i)[0][0] == np.where(eR == closeIt(eR, energyS)):
            energySelected = i
    return energySelected 

def workFile(height, typeP, seq, minE, maxE): #*Files to work with, IF maxE == None => all energies since minE will be select
    filterFilesV = filesFilter(height, typeP, seq)
    selectIds = np.array([])
    for i in filterFilesV:
        if maxE == None:
            if np.where(filterFilesV == i)[0][0] >= np.where(eR == closeIt(eR, minE))[0][0]: #-1
                selectIds = np.append(selectIds, i)
        else:
            if np.where(filterFilesV == i)[0][0] >= np.where(eR == closeIt(eR, minE))[0][0] and np.where(filterFilesV == i)[0][0] <= np.where(eR == closeIt(eR, maxE))[0][0]: #-1
                selectIds = np.append(selectIds, i)
    return selectIds

class particleAnalysis():
    description = "This class return an DataFrame object"
    def __init__(self, height, pParticle, component, eMin, eMax, seq) -> description:
        self.__height = height
        self.__eMin = eMin
        self.__eMax = eMax
        self.__pParticle = pParticle
        self.seq = seq
        self.__component = component

    def separateComponent(self, dataEach):
        if self.__component == "muonic":
            id = np.array([5, 6, 75, 76])
            muons = pd.DataFrame(dataEach.loc[dataEach.index.intersection(id), : ]).drop_duplicates()
            return muons
        elif self.__component == "electromagnetic":
            id = np.array([1, 2, 3])
            electromagnetic = pd.DataFrame(dataEach.loc[dataEach.index.intersection(id), : ]).drop_duplicates()
            return electromagnetic
        elif self.__component == "pionic":
            id = np.array([7, 8, 9])
            pionic = pd.DataFrame(dataEach.loc[dataEach.index.intersection(id), : ]).drop_duplicates()
            return pionic
        elif self.__component == "muonicelectronic":
            id = np.array([5, 6, 75, 76, 2, 3])
            muonicelectronic = pd.DataFrame(dataEach.loc[dataEach.index.intersection(id), : ]).drop_duplicates()
            return muonicelectronic
        else: 
            print(":C I can't study this type of particle, you say to my programer that implement this analysis!!!")

    def totalEnergyDist(self): #This method return a graph of secondary average energy particle 
        filesToStudy = workFile(self.__height, self.__pParticle, self.seq, self.__eMin, self.__eMax)
        exportData = {"Energy_Mean" : np.array([]), "Std_Mean" : np.array([])}
        for i in filesToStudy:
            data = pd.read_csv("./{}/{}/{}".format(self.__height, self.__pParticle, i), header=0, delim_whitespace=True).set_index("NEvent")
            datasubShowerMean = pd.DataFrame()
            for j in np.arange(1, NSHOWS + 1, 1):
                dataEach = pd.DataFrame(data.loc[data.index.intersection([j]), : ]).drop_duplicates().set_index("Type")
                typeParData = self.separateComponent(dataEach) ###Data of muonic or electromagnetic component of each subshower
                eMean = np.sum(typeParData["Energy"]*typeParData["Weight"])/np.sum(typeParData["Weight"])
                datasubShowerMean["Shower{}".format(j)] = np.array([eMean])
            exportData["Energy_Mean"] = np.append(exportData["Energy_Mean"], np.nanmean(datasubShowerMean.loc[0, : ]))
            exportData["Std_Mean"] = np.append(exportData["Std_Mean"], np.nanstd(datasubShowerMean.loc[0, : ]))
        res = pd.DataFrame(exportData)
        res["Primary_Energy"] = selectEnergy(eR, self.__eMin) ###FINAL DATAFRAME
        res.to_csv("{}_{}_{}.txt".format(self.__height, self.__pParticle, self.__component), sep = "\t", index = False)
        return res

    def lateralparticleDist(self, dataFrame, option, rmin, rmax, bin):
        if type(dataFrame) == "str":
            data = pd.read_csv("./{}/{}/{}".format(self.__height, self.__pParticle, dataFrame), header = 0, delim_whitespace=True).set_index("NEvent")
        else:
            file = specificEnergy(self.__height, dataFrame, self.__pParticle, seq)
            data = pd.read_csv("./{}/{}/{}".format(self.__height, self.__pParticle, file), header = 0, delim_whitespace=True).set_index("NEvent")
        bins = pd.interval_range(rmin, rmax, freq=bin)
        rMean = np.array([(i + j)/2 for i, j in zip(bins.left, bins.right)])
        dataFinalShowers = pd.DataFrame()
        for i in np.arange(1, NSHOWS + 1, 1):
            dataEach = pd.DataFrame(data.loc[data.index.intersection([i]), : ]).drop_duplicates().set_index("Type")
            dataSeparated = self.separateComponent(dataEach) ###MUONS FOR EXAMPLE
            if option == "perpenticular":
                if self.__component == "electromagnetic":
                    dataSeparated["Total_Energy"] = dataSeparated["Weight"]*dataSeparated["Energy"]
                    sum = dataSeparated.groupby(pd.cut(dataSeparated.R, bins = bins))["Total_Energy"].sum()
                if self.__component == "muonic" or self.__component == "pionic":
                    sum = dataSeparated.groupby(pd.cut(dataSeparated.R, bins = bins))["Weight"].sum()
                dataFinalShowers["Shower{}".format(i)] = np.array(sum.values)
            elif option == "inclined":
                dataSeparated.R = dataSeparated.R*np.sqrt((np.sin(dataSeparated.primPhi-np.arctan(dataSeparated.Y/dataSeparated.X)))**2+(np.cos(dataSeparated.primPhi-np.arctan(dataSeparated.Y/dataSeparated.X)))**2*(np.cos(dataSeparated.primTheta))**2)
                if self.__component == "electromagnetic":
                    dataSeparated["Total_Energy"] = dataSeparated["Weight"]*dataSeparated["Energy"]
                    sum = dataSeparated.groupby(pd.cut(dataSeparated.R, bins = bins))["Total_Energy"].sum()
                if self.__component == "muonic" or self.__component == "pionic":
                    sum = dataSeparated.groupby(pd.cut(dataSeparated.R, bins = bins))["Weight"].sum()
                dataFinalShowers["Shower{}".format(i)] = np.array(sum.values)
            else:
                print ("option selected doesn't existe")
                return
        ####* FOR THE AVERAGES
        FINAL = {"mean" : np.array([]), "std" : np.array([])}
        for i in np.arange(0, rMean.size):
            FINAL["mean"] = np.append(FINAL["mean"], np.nanmean(dataFinalShowers.loc[i, :]))
            FINAL["std"] = np.append(FINAL["std"], np.nanstd(dataFinalShowers.loc[i, :]))
        if self.__component == "electromagnetic":
            titleDataFinal = "Total_Energy"
        if self.__component == "muonic" or self.__component == "pionic":
            titleDataFinal = "Total_Weight"
        res = pd.DataFrame({'Internal_radius': bins.left, 'External_radius': bins.right, 'R':rMean, '{}'.format(titleDataFinal): FINAL["mean"], "Standar_deviation" : FINAL["std"]})
        res.to_csv('LateralDistribution{}_{}_{}.txt'.format(self.__height, self.__pParticle, self.__component), sep='\t', index=False)
        return res

    def lateralEnergyDist(self, dataFrame, emax, bin):
        if type(dataFrame) == str:
            data = pd.read_csv("./{}/{}/{}".format(self.__height, self.__pParticle, dataFrame), header = 0, delim_whitespace=True).set_index("NEvent")
        else:
            file = specificEnergy(self.__height, dataFrame, self.__pParticle, seq)
            data = pd.read_csv("./{}/{}/{}".format(self.__height, self.__pParticle, file), header = 0, delim_whitespace=True).set_index("NEvent")
        data["Energy"] = data["Energy"]/10**9
        bins = pd.interval_range(0, emax + bin, freq=bin)
        eMean = np.array([(i+j)/2 for i,j in zip(bins.left, bins.right)])
        dataFinalShowers = pd.DataFrame()
        for i in np.arange(1, NSHOWS + 1, 1):
            dataEach = pd.DataFrame(data.loc[data.index.intersection([i]), : ]).drop_duplicates().set_index("Type")
            dataSeparated = self.separateComponent(dataEach) ###MUONS FOR EXAMPLE
            sumParticles = dataSeparated.groupby(pd.cut(dataSeparated.Energy, bins = bins))["Weight"].sum()
            dataFinalShowers["Shower{}".format(i)] = np.array(sumParticles.values)
        ####* FOR THE AVERAGES
        FINAL = {"mean" : np.array([]), "std" : np.array([])}
        for i in np.arange(0, eMean.size):
            FINAL["mean"] = np.append(FINAL["mean"], np.nanmean(dataFinalShowers.loc[i, :]))
            FINAL["std"] = np.append(FINAL["std"], np.nanstd(dataFinalShowers.loc[i, :]))
        res = pd.DataFrame({'Internal_Energy': bins.left, 'External_Energy': bins.right, 'eMean': eMean, "Total_Weight" : FINAL["mean"], "Standar_deviation" : FINAL["std"]})
        res.to_csv('LateralDistribution{}_{}_{}.txt'.format(self.__height, self.__pParticle, self.__component), sep='\t', index=False)
        return res

    def e0Predict(self):
        filesToStudy = workFile(self.__height, self.__pParticle, self.seq, self.__eMin, self.__eMax)
        exportData = {"Electronic_Num" : np.array([]), "Std_Mean_Electronic" : np.array([]), "Muonic_Num" : np.array([]), "Std_Mean_Muonic" : np.array([])}
        for i in filesToStudy:
            data = pd.read_csv("./{}/{}/{}".format(self.__height, self.__pParticle, i), header=0, delim_whitespace=True).set_index("NEvent")
            dataMuonicMean = pd.DataFrame()
            dataElectronicMean = pd.DataFrame()
            for j in np.arange(1, NSHOWS + 1, 1):
                dataEach = pd.DataFrame(data.loc[data.index.intersection([j]), : ]).drop_duplicates().set_index("Type")
                typeParData = self.separateComponent(dataEach) ###Data of muonic or electromagnetic component of each subshower
                muonic = pd.DataFrame(typeParData.loc[typeParData.index.intersection([5, 6, 75, 76]), : ])
                electronic = pd.DataFrame(typeParData.loc[typeParData.index.intersection([2, 3]), : ])
                nMuonic = np.nansum(muonic["Weight"])
                nElectronic = np.nansum(electronic["Weight"]) 
                dataMuonicMean["Shower{}".format(j)] = np.array([nMuonic])
                dataElectronicMean["Shower{}".format(j)] = np.array([nElectronic])
            exportData["Electronic_Num"] = np.append(exportData["Electronic_Num"], np.nanmean(dataElectronicMean.loc[0, : ]))
            exportData["Std_Mean_Electronic"] = np.append(exportData["Std_Mean_Electronic"], np.nanstd(dataElectronicMean.loc[0, : ]))
            exportData["Muonic_Num"] = np.append(exportData["Muonic_Num"], np.nanmean(dataMuonicMean.loc[0, : ]))
            exportData["Std_Mean_Muonic"] = np.append(exportData["Std_Mean_Muonic"], np.nanstd(dataMuonicMean.loc[0, : ]))
        res = pd.DataFrame(exportData)
        res["Primary_Energy"] = selectEnergy(eR, self.__eMin, self.__eMax) ###FINAL DATAFRAME
        res.to_csv("{}_{}_{}.txt".format(self.__height, self.__pParticle, self.__component), sep = "\t", index = False)
        return res

#def __init__(self, height, pParticle, component, eMin, eMax, seq) -> description:

sealeveliron10TeV = particleAnalysis("iron", 0, "muonic", None, None, 50)
sealeveliron10TeV.lateralEnergyDist("Data.txt", 10, 1)
yanqueiron10TeV = particleAnalysis("iron", 4850, "muonic", None, None, 50)
yanqueiron10TeV.lateralEnergyDist("Data.txt", 10, 1)
chachaniiron10TeV = particleAnalysis("iron", 6500, "muonic", None, None, 50)
chachaniiron10TeV.lateralEnergyDist("Data.txt", 10, 1)
overeverestiron10TeV = particleAnalysis("iron", 10000, "muonic", None, None, 50)
overeverestiron10TeV.lateralEnergyDist("Data.txt", 10, 1)

# sealevelproton10TeV = particleAnalysis("proton", 0, "muonic", None, None, 50)
# sealevelproton10TeV.lateralEnergyDist("Data.txt", 10, 1)
# yanqueproton10TeV = particleAnalysis("proton", 4850, "muonic", None, None, 50)
# yanqueproton10TeV.lateralEnergyDist("Data.txt", 10, 1)
# chachaniproton10TeV = particleAnalysis("proton", 6500, "muonic", None, None, 50)
# chachaniproton10TeV.lateralEnergyDist("Data.txt", 10, 1)
# overeverestproton10TeV = particleAnalysis("proton", 10000, "muonic", None, None, 50)
# overeverestproton10TeV.lateralEnergyDist("Data.txt", 10, 1)

        

        
