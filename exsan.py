"""
Title:      EXSAN: An (E)NDF (C)ross (S)ection (An)alysis

LA-CC ID:   Los Alamos National Laboratory Computer Code LA-CC-19-002

Author:     H. Omar Wooten, PhD, DABR
            Los Alamos National Laboratory
            hasani@lanl.gov

Abstract:   This application provides an intutive  graphical interface for:
            1) Downloading complete ENDF libraries from the National Nuclear Data Center (NNDC)

            2) Automated batch processing of ENDF libraries with the Los Alamos code NJOY

            3) Plotting continuous energy and multigroup neutron reaction and photoatomic cross sections

            4) Plotting cross sections for isotopic mixes (multigroup only)

            5) Automated plotting and saving of cross section figures with an input file.

            6) Automated analysis of ENDF/B-VII and ENDF/B-VIII.0 radioactive decay libraries and plotting of the decay particles.

References:
            Wooten, H.O., "An application for streamlined and automated ENDF Cross Section Analysis and visualization," Ann Nuc Energy, 129, 482-486, 2019

            Wooten, H.O., A Python Based ENDF Cross Section Explorer,
            LA-UR-18-25324, Los Alamos National Laboratory, 2018.
"""
from Tkinter import *
from Tkinter import _setit as setit
import ttk
import tkFileDialog
import tkFont
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as pl
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from matplotlib.font_manager import FontProperties
import numpy as np
import sys, os, operator
from glob import glob
from datadic import *
import urllib, urllib2
import re
from scipy.interpolate import CubicSpline as cs
import time, datetime
from copy import copy
from pdb import set_trace as st
from pylab import connect, draw
import subprocess as sp
from shutil import move as mv
from collections import OrderedDict as od
from argparse import ArgumentParser


opt = ArgumentParser(usage='python exsan.py [options]')
opt.add_argument('-n',  '--njoy', dest='njoyPath',  help='Absolute path to NJOY2016')
opt.add_argument('-b',  '--batch', dest='batchFile',  help='batch input file')
opt.add_argument('-v',   '--verbose', action='store_true', help='verbosity')
opt.add_argument('-d', '--demo', action='store_true', help='demo mode')
opt.add_argument('-a', '--auto', action='store_true', help='minimal automation when starting up')
opt.add_argument('-l', '--log', action='store_true', help='write log file exsan.log')
options = opt.parse_args()

demoIsotopes = ['H-1 total','Li-6 total']
os.system('clear')

#===========================================================
# set up a class for global variables.
# Consise description and reasoning given here:
# https://pythonconquerstheuniverse.wordpress.com/2010/10/20/a-globals-class-pattern-for-python/
# All global variables are in the 'mem' namespace. This includes Tkinter widgets
#===========================================================
class mem: pass

#===========================================================
# set plot and font sizes
pl.rcParams.update({'font.size': 10}) # set plot fonts
#===========================================================

class makeAlphabetButton(object):
    '''
    This class makes a cascading Menubutton of Letter > Elements > Isotopes.
    '''
    def __init__(self, frame, elementsIsotopesDict,row, selection, isotopeList, labelVar, weightVar):
        self.frame = frame
        self.elementsIsotopesDict = elementsIsotopesDict
        self.periodicTableDict = periodicTableDict
        self.row = row
        self.selection = selection
        self.isotopeList = isotopeList
        self.labelVar = labelVar
        self.weightVar = weightVar
        self.makeAZdict()
        self.makeAZbutton()

    def makeAZdict(self):
        azList, elements = ([] for i in range(2))
        for iso in self.isotopeList:
            azList.append(iso[0])
            elements.append(iso.split('-')[0])
        azList = np.unique(azList)
        elements = np.unique(elements)
        self.azDict = {}
        for letter in azList:
            self.azDict[letter]=[]
            for el in elements:
                if el[0]==letter:
                    self.azDict[letter].append(el)

    def makeAZbutton(self):
        az_btn = Menubutton(self.frame,text='',font=mem.MENUBUTTON)
        az_btn.menu = Menu(az_btn)
        az_btn["menu"] = az_btn.menu

        for letter in sorted(self.azDict.keys()):
            az_menu = Menu(az_btn.menu)
            for el in self.azDict[letter]:
                el_menu = Menu(az_menu)
                for iso in self.elementsIsotopesDict[el]:
                    el_menu.add_command(label=iso,command=lambda iso=iso: self.setVars(iso))
                az_menu.add_cascade(label=el, menu=el_menu)
            az_btn.menu.add_cascade(label=letter, menu=az_menu)
        az_btn.configure(width=10)
        az_btn.grid(column=2, row=self.row, padx=3, pady=3)

    def setVars(self, iso):
        self.selection.set(iso)
        self.labelVar.set(iso)
        if mem.percentType_var.get() == 0:
            self.weightVar.set('1.0')
        else:
            self.weightVar.set('1')

#===========================================================
# A class for complex multilevel button in the periodic table
#===========================================================
class makePeriodicTableButton(object):
    '''
    This class makes a cascading Menubutton of element > isotopes > MT cross sections.
    '''
    def __init__(self, el, frame, elementsIsotopesDict, periodicTableDict, isotopesMTdict, decayFlag, fileDir, flag_pOrM=False ):
        for k,v in locals().iteritems():
            if k != 'self':
                setattr(self, k, v)

        self.allFlag =[False, True]

        self.makeComplexButton()

    def makeComplexButton(self):
        allListMaster = {}
        allListMasterMaster = {}
        allListMasterMaster[self.el] = {}

        el_btn = Menubutton(self.frame, text=self.el+'_'+str(self.periodicTableDict[self.el][0]), font=mem.MENUBUTTON)
        el_btn.menu  =  Menu(el_btn)
        el_btn["menu"]  =  el_btn.menu

        for iso in sorted(self.elementsIsotopesDict[self.el]):
            allListMaster[iso]=[]
            allListMasterMaster[self.el][iso]=[]
            allList = []
            iso_menu = Menu(el_btn.menu)
            for mtItem in self.isotopesMTdict[iso]:
                oneList = []
                oneList.append(iso+' '+mtItem)
                if self.decayFlag:
                    iso_menu.add_command(label=mtItem, command=lambda iso=iso, fileDir=self.fileDir: get_decay_data(iso, fileDir))
                else:
                    iso_menu.add_command(label=mtItem,command=lambda oneList=oneList:getInfo2(oneList, self.flag_pOrM, self.allFlag[0]))
                    allList.append(iso+' '+mtItem)
                    allListMaster[iso].append(iso+' '+mtItem)
                    allListMasterMaster[self.el][iso].append(iso+' '+mtItem)
            if not self.decayFlag:
                iso_menu.add_command(label='All',command=lambda allList=allList:getInfo2(allList, self.flag_pOrM, self.allFlag[1]))
            el_btn.menu.add_cascade(label=iso, menu=iso_menu)
        el_btn.grid(row=self.periodicTableDict[self.el][2], column=self.periodicTableDict[self.el][1])
        for x in range(19):
            Grid.columnconfigure(self.frame, x, weight=1)

        for y in range(11):
            Grid.rowconfigure(self.frame, y, weight=1)

        for k,v in allListMaster.iteritems():
            mem.allListMaster[k] = v

        for k,v in allListMasterMaster.iteritems():
            mem.allListMasterMaster[k] = v

#===========================================================
#  A base class that makes a dictionary of MF/MT combination
#  from ENDF file
#===========================================================
class mfmtDict(object):
    def __init__(self, lines):
        self.dic = {}
        self.dicNest = {}
        self.lines = lines
        self.get_section_data()

    def get_section_data(self):
        MF = [l[70:72] for l in self.lines]
        MT = [l[72:75] for l in self.lines]
        v = [l[70:75] for l in self.lines]
        for i, line in enumerate(self.lines):
            mfi = int(MF[i].strip())
            mti = int(MT[i].strip())
            if mfi >=3:
                if (not mfi in self.dic.keys()):
                    self.dic[mfi]=[]
                else:
                    if mti >0:
                        self.dic[mfi].append(mti)
        for key in self.dic.keys():
            self.dic[key] = np.unique(self.dic[key])

        for key in self.dic.keys():
            self.dicNest[key]={}
            for value in self.dic[key]:
                mfmt = "%2s%3s"%(key, value)
                n = len(v)
                iFirst=0
                iFirst = v.index(mfmt)    # first occurrence of MT/MF section
                iLast = n - v[::-1].index(mfmt)
                self.dicNest[key][value]=[iFirst, iLast]

#===========================================================
# This class inherits from mtmfdict class and creates an
# instance of multigroup cross sections for plotting
#===========================================================
class groupr(mfmtDict):
    def __init__(self, lines, MF, MT):
        mfmtDict.__init__(self, lines)
        self.n = False
        self.g = False
        self.linesOrig = self.lines
        self.lines = [line[:66] for line in self.lines]
        self.lines = filter(None,[line.replace('+','e+').replace('-','e-') for line in self.lines[2:]])
        self.temp, = [float(j) for j in self.lines[0].split()[:1]]
        self.nGroups, self.gGroups = [int(j) for j in self.lines[0].split()[2:4]]
        self.getInfoLines
        self.headerInfo = []
        self.getInfoLines()
        self.eBins_n = self.headerInfo[2:self.nGroups+3]
        self.eBins_g = self.headerInfo[self.nGroups+3:]
        self.MF = MF
        self.MT = MT
        self.getFluxAndXS()
        self.xs[0] = self.xs[1]

    def getInfoLines(self):
        v = [l[70:75] for l in self.linesOrig]
        n = len(v)
        iFirst=0
        iFirst = v.index(' 1451')     # first occurrence of MT/MF section
        iLast = n - v[::-1].index(' 1451')
        for line in self.lines[iFirst:iLast-2]:
            self.headerInfo.append(line[:72].split())
        self.headerInfo = [float(item) for items in self.headerInfo for item in items]

    def getFluxAndXS(self):
        data = self.lines[ self.dicNest[self.MF][self.MT][0]-2 : self.dicNest[self.MF][self.MT][1]-2:2][1:]
        nGroups = self.lines[ self.dicNest[self.MF][self.MT][0]-1: self.dicNest[self.MF][self.MT][1]-2:2]

        self.flux = []
        xsTmp = []
        self.xs =[1.e-10]*len(self.eBins_n)
        g = []
        for n in nGroups:
            g.append(int(n.split()[5]))

        for i, d in enumerate(data):
            flux, xs = [float(item) for item in d.split()]
            self.flux.append(flux)
            xsTmp.append(xs)
        for i, ig in enumerate(g):
            self.xs[ig] = xsTmp[i]

#===========================================================
# This class allows for multicursor snap to data
#===========================================================
class SnapToCursor(object):
    """
    This is based on the SnaptoCursor class from https://goo.gl/U62eNR.
    I modified it to be able to snap to multiple datasets.
    H. Omar Wooten
    7/4/18
    """
    def __init__(self, ax, data, fig, titles, textXY, textLabel):
        self.ax = ax
        self.lx = ax.axhline(color='k')  # the horiz line
        self.ly = ax.axvline(color='k')  # the vert line
        self.data = data # a list of tuples of x,y datasets. [(xData1,yData1),(xData2,yData2)]
        self.x = []
        for d in self.data:
            self.x.append(d[0])
        self.textXY = textXY
        self.textLabel = textLabel
        self.txt = ax.text(self.textXY[0], self.textXY[1], '', transform=ax.transAxes)
        self.txt.set_bbox(dict(facecolor='white', alpha=0.9, edgecolor='black'))
        self.fig = fig
        self.titles = titles

    def mouse_move(self, event):
        if not event.inaxes:
            return
        if len(self.data) == 0:
            return
        x, y = event.xdata, event.ydata
        indx = []
        for xData in self.x:
            indx.append(min(range(len(xData)), key=lambda i: abs(xData[i]-x)))
        tmp = []
        for i in range(len(self.data)):
            tmp.append(np.abs(y-self.data[i][1][indx[i]]))
        indy = np.argmin(tmp)
        x = self.data[indy][0][indx[indy]]
        y = self.data[indy][1][indx[indy]]
        # update the line positions
        self.lx.set_ydata(y)
        self.ly.set_xdata(x)
        self.txt.set_text('%s: \neV: %.3e\n%s: %.3e'%(self.titles[indy], x, self.textLabel, y))
        self.fig.draw()

#===========================================================
# This class allows standard output to be written both to
# a log file and to the terminal.
#===========================================================
class Log(object):
    '''
    This class allows standard output to be written both to a log file and to the terminal.
    '''
    def __init__(self, message=None):
        self.terminal = sys.stdout
        self.log = open('exsan.log', 'a')
        self.message = message
        if not self.message == None:
            self.write(message)

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)



#===========================================================
#  This function serves as a gate for reading NJOY-processed
#  files based on whether they are a GROUPR output or not
#===========================================================
def makeMultiGroup(lines, MF, MT):
    if lines[1].split()[4] == '-1':
        g = groupr(lines, int(MF), int(MT))
        return g
    else:
        g = mfmtDict(lines)
        return g

#===========================================================
#    Return data lines for specific MF & MT as x,y pairs
#===========================================================
def get_section_data(MF=3, MT=1):
    v = [l[70:75] for l in mem.lines]
    n = len(v)
    mfmtstr = '%2s%3s' % (MF, MT)       # search string
    mem.iFirst=0
    mem.iFirst = v.index(mfmtstr)           # first occurrence of MT/MF section
    mem.iLast = n - v[::-1].index(mfmtstr)  # last occurrence of MT/MF section

    mem.data = []
    checkPos = [i*11 for i in range(0,6)]

    for l in mem.lines[mem.iFirst+3:mem.iLast]:
        if l.count('E') > 0:
            a = filter(None, [l[i*11:i*11+11] for i in range(0,6)])
            b = filter(None,[ai.replace(' ','') for ai in a])
            b = [float(bi) for bi in b]
        else:
            lTmp = list(l)
            for pos in checkPos:
                lTmp[pos]=' '
            l = ''.join(lTmp)
            a = filter(None, [l[i*11:i*11+11].replace(' ','').replace('-','e-').replace('+','e+') for i in range(0,6)])
            b = [float(ai) for ai in a]

        mem.data.append(b)

    flatData = [item for sublist in mem.data for item in sublist]

    energies, xsecs = flatData[::2],flatData[1::2]

    xNoZero = []
    yNoZero = []
    zeroTally=0
    for i in range(len(energies)):
        if energies[i] == 0 or xsecs[i] == 0:
            zeroTally+=1
            continue
        else:
            xNoZero.append(energies[i])
            yNoZero.append(xsecs[i])

    [x for (y,x) in sorted(zip(yNoZero,xNoZero), key=lambda pair: pair[0])]
#     if zeroTally > 0:
#         print str(zeroTally)+' points = 0 were removed for to allow log-log plotting'
    return xNoZero, yNoZero


#===========================================================
#    Get radioactive decay data
#===========================================================
def get_decay_data(iso, dir, MF=8, MT=457, miniParse=False):
    '''
    Parse ENDF decay files.

    Reference:
    A. Trkov, D. A. Brown, "ENDF-6 Formats Manual: Data Formats and Procedures for the Evaluated Nuclear Data Files ENDF/B-VI, ENDF/B-VII and ENDF/B-VIII," BNL-203218-2018-INRE, Brookhaven National Laboratory, 2018.
    '''

    def parseLine(line):
        '''
        Parse a line of 6*11 characters, [:66], of an ENDF tape converting words with +/- to floats, and everything else as integers. Return the list of 6 items.
        '''
        line = line[:66]
        # If a number is negative, Don't replace its '-' sign
        if ' -' in line:
            line = line.replace(' -','zz')
        line = ''.join(list(line.replace('-','e-').replace('+','e+').replace('zz',' -'))).split()
        return [float(h) if '+' in h or '-' in h else int(h) for h in line]

    def parseLineIndices(line, indices=None):
        '''
        Parse a line of 6*11 characters, [:66], and return specific indices of the 6-item list.
        '''
        if indices == None:
            indices = [i for i in range(6)]
        return [parseLine(line)[x] for x in indices]

    def resetDecay():
        '''
        Reset dictionary keys to blank lists
        '''
        decay['rtyp2'],  \
        decay['er'],     \
        decay['d_er'],   \
        decay['type'],   \
        decay['ri'],     \
        decay['d_ri'],   \
        decay['ris'],    \
        decay['d_ris'],  \
        decay['ricc'],   \
        decay['d_ricc'], \
        decay['rick'],   \
        decay['d_rick'], \
        decay['ricl'],   \
        decay['d_ricl'], \
         = ([] for i in range(14))

    decayModes = copy(decayTypes)
    decayModes[3.]= 'IT'

    fileName = '%s/%s.txt'%(dir, ''.join(iso.split('-')))

    with open(fileName, 'r') as f:
        lines = f.readlines()


    v = [l[70:75] for l in lines]
    n = len(v)
    mfmtstr = '%2s%3s' % (MF, MT)       # search string
    iFirst=0
    iFirst = v.index(mfmtstr)           # first occurrence of MT/MF section
    iLast = n - v[::-1].index(mfmtstr)  # last occurrence of MT/MF section

    data = []
    decayAll = {}
    decay = {}
    decaySummary = {}

    # there's probably a better way to to this....
    decayAll['avgDecayEnergy'], \
    decayAll['rtyp'],  \
    decayAll['rfs'],    \
    decayAll['Q'],      \
    decayAll['dQ'],     \
    decayAll['br'],     \
    decayAll['d_br'],   \
    decay['rtyp2'],  \
    decay['er'],     \
    decay['d_er'],   \
    decay['type'],   \
    decay['ri'],     \
    decay['d_ri'],   \
    decay['ris'],    \
    decay['d_ris'],  \
    decay['ricc'],   \
    decay['d_ricc'], \
    decay['rick'],   \
    decay['d_rick'], \
    decay['ricl'],   \
    decay['d_ricl'], \
     = ([] for i in range(21))

    decayAll['za'], decayAll['awr'], decayAll['nsp'] = parseLineIndices(lines[iFirst], [0,1,5])

    decayAll['tHalf'], decayAll['nc'] = parseLineIndices(lines[iFirst+1], [0,4])
    decayAll['nc'] = decayAll['nc']/2

    if decayAll['nsp']==0. and decayAll['tHalf']==0:
        mem.resultsText.delete('1.0', END)
        mem.resultsText.insert(INSERT, '%s is stable.\n'%(iso))
        return

    if decayAll['nc'] == 3:
        iStart, iEnd = iFirst+2, iFirst+3
    else:
        iStart, iEnd = iFirst+2, iFirst+12

    for line in lines[iStart:iEnd]:
        decayAll['avgDecayEnergy'] += parseLine(line)

    decayAll['spin'], decayAll['par'], decayAll['nDecayModes'] = [parseLine(lines[iEnd])[x] for x in [0, 1, 5]]

    # Loop over decay modes (NDK)
    data = []
    iStart = iEnd+1
    for ndki in range(decayAll['nDecayModes']):
        data += parseLine(lines[iStart+ndki])

    iStart += ndki+1
    decayAll['rtyp'] = [decayModes[x] for x in data[::6]]
    decayAll['rfs'] = data[1::6]
    decayAll['Q'] = data[2::6]
    decayAll['d_Q'] = data[3::6]
    decayAll['br'] = data[4::6]
    decayAll['d_br'] = data[5::6]

    decaySummary[u'T-\u00BD'] = [decayAll['tHalf']]
    decaySummary['Decay Modes'] = zip(decayAll['rtyp'], decayAll['br'])


    # Loop over spectra (NSP)
    for nspi in range(decayAll['nsp']):
        # resetDecay()
        decay['styp'], \
        decay['lcon'], \
        decay['lcov'], \
        decay['ner'] = parseLineIndices(lines[iStart], [1, 2, 3, 5])
        decay['styp'] = decayTypes[decay['styp']]
        iStart += 1

        decay['fd'],     \
        decay['d_fd'],   \
        decay['erAv'],   \
        decay['d_erAv'], \
        decay['fc'],     \
        decay['d_fc'] = parseLineIndices(lines[iStart])
        iStart+=1

        # Loop over discrete energies for a given spectral type (NER)
        for neri in range(decay['ner']):
            er, d_er, nt = parseLineIndices(lines[iStart], [0, 1, 4])
            iStart += 1

            rtyp2, type, ri, d_ri, ris, d_ris = parseLineIndices(lines[iStart])
            # is there any way to do this in a loop? Like this?
            # for i in [er, d_er]:
            #     varName = [k for k,v in locals().iteritems() if v==i[0]]
            #     mem.decay[varName].append(i)
            decay['er'].append(er)
            decay['d_er'].append(d_er)
            decay['rtyp2'].append(rtyp2)
            decay['type'].append(type)
            decay['ri'].append(ri)
            decay['d_ri'].append(d_ri)
            decay['ris'].append(ris)
            decay['d_ris'].append(d_ris)
            iStart += 1

            if nt > 6:
                ricc, d_ricc, rick, d_rick, ricl, d_ricl = parseLineIndices(lines[iStart])
                decay['ricc'].append(ricc)
                decay['d_ricc'].append(d_ricc)
                decay['rick'].append(ricc)
                decay['d_rick'].append(ricc)
                decay['ricl'].append(ricc)
                decay['d_ricl'].append(ricc)
                iStart += 1

        radSummary = zip(decay['er'], decay['ri'])
        radSummary.sort(key=lambda x: x[1])
        decaySummary[decay['styp']] = radSummary

    if miniParse:
        mem.decaySummary = copy(decaySummary)
        mem.decaySummary['isotope']=iso
        return

    mem.resultsText.delete('1.0', END)
    mem.resultsText.insert(INSERT, '     -- Radioactive Decay Summary: %s --\n\n'%(iso))

    value, units = getDecayUnits(decaySummary[u'T-\u00BD'][0])
    if value>1000:
        mem.resultsText.insert(INSERT, u'%18s %10.3e %s\n\n'%(u'T-\u00BD:', value, units))
    else:
        mem.resultsText.insert(INSERT, u'%18s %10.3f %s\n\n'%(u'T-\u00BD:', value, units))

    mem.resultsText.insert(INSERT, '='*50+'\n\n')

    for k,v in decaySummary.iteritems():
        # if k==u'T-\u00BD': # already printing half life at the top of the summary
        #     continue
        mem.resultsText.insert(INSERT, '%s\n\n'%(k))
        if k == 'Decay Modes':
            for vi in v:
                mem.resultsText.insert(INSERT,'%20s %14.4e \n'%(vi[0], vi[1]))
            mem.resultsText.insert(INSERT, '\n'+'='*50+'\n\n')
            continue


        # convert half life from seconds to reasonable units
        if k == u'T-\u00BD': # half life
            for kk,vv in mem.halfLifeScale.iteritems():
                if decaySummary[k][0] >= vv[0] and decaySummary[k][0] < vv[1]:
                    if decaySummary[k][0] >3.15e11:
                        mem.resultsText.insert(INSERT,'%20.4e %s\n\n'%(decaySummary[k][0]/vv[0], kk))
                    else:
                        if kk=='sec':
                            mem.resultsText.insert(INSERT,'%20.4f %s\n\n'%(decaySummary[k][0], kk))
                        else:
                            mem.resultsText.insert(INSERT,'%20.4f %s\n\n'%(decaySummary[k][0]/vv[0], kk))
                    mem.resultsText.insert(INSERT, '='*50+'\n\n')
            continue

        elif k in decayTypes.values():
            mem.resultsText.insert(INSERT, '%20s %16s\n\n'%('E (eV)','Fraction'))

        for vi in v:
            try:
                mem.resultsText.insert(INSERT,'%21.4e %16.4e \n'%(vi[0], vi[1]))
            except:
                mem.resultsText.insert(INSERT, '\t%s\n'%(str(vi)))
        mem.resultsText.insert(INSERT, '\n'+'='*50+'\n')

        mem.decaySummary = copy(decaySummary)
        mem.decaySummary['isotope']=iso
        mem.decayPlotButton['state'] = 'normal'

#===========================================================
#   Create a dictionary of MT/MF combinations from lines
#===========================================================
def make_dict(fileDir):
    v = [l[70:75] for l in mem.lines]
    MF = set([l[70:72].strip() for l in mem.lines])
    MT = [l[72:75].strip() for l in mem.lines]
    MTMF_dict = dict((mf,[]) for mf in MF)
    try:
        MTMF_dict.pop('0')
    except:
        pass
    mtTmp =''
    for vi in v:
        mf = vi[0:2].strip()
        mt = vi[2:5].strip()
        if mt=='0':
            continue
        if not  mt == mtTmp:
            MTMF_dict[mf].append(mt)
        mtTmp = mt

    if mem.particle.get()==1:
        range = '23'
    else:
        if 'decay' in fileDir:
            range = '8'
        else:
            range = '3'

    MT_verbose = []
    for item in MTMF_dict[range]:
        if item in mtdic2.keys():
            MT_verbose.append(mtdic2[item])
    return MT_verbose

#===========================================================
#  Incorporate progress bars as a visual aid
#===========================================================
def progressBar(lenFiles):
    popup = tk.Toplevel()
    popup_lab = Label(popup, text="Reading files ...")
    popup_lab.grid(row=0, column=0)

    status = 0
    status_var = DoubleVar()
    progress_bar = Progressbar(popup, variable=status_var, maximum=100)
    progress_bar.grid(row=1, column=0)

    statusInc = 100.0/lenFiles
    for i in range(lenFiles):
        popup.update()
        status += statusInc
        status_var.set(status)
    popup.destroy()

#===========================================================
#   Update file list options
#===========================================================
def update_files(event=None):
    if sorted(dirDict2.keys())[mem.verSelect.get()] in ['ENDF/B-VII.1','ENDF/B-VIII.0']:
        fileDirectory = dirDict2[np.sort(dirDict2.keys())[mem.verSelect.get()]][mem.particle.get()]
    else:
        mem.particle.set(0)
        fileDirectory = dirDict2[np.sort(dirDict2.keys())[mem.verSelect.get()]][0]

    particles = []
    for vi in dirDict2[sorted(dirDict2.keys())[mem.verSelect.get()]]:
        if 'parsed' in vi or 'neutrons' in vi:
            particles.append('n')
        if 'photo' in vi:
            particles.append('x')

    if 'n' in particles:
        mem.particleSelect['neutrons'].configure(state='normal')
    else:
        mem.particleSelect['neutrons'].configure(state='disabled')

    if 'x' in particles:
        mem.particleSelect['x-rays'].configure(state='normal')
    else:
        mem.particleSelect['x-rays'].configure(state='disabled')

    fileDirs = {}
    fileDirs[fileDirectory] = mem.tab11
    if mem.particle.get() == 0:
        fileDirs[fileDirectory+'/multigroup'] = mem.tab12
        fileDirs['/'.join(fileDirectory.split('/')[:-1])+'/decay'] = mem.tab13
    # mem.allListMaster = {}
    # mem.allListMasterMaster = {}
    mem.masterTracker = {}

    for mem.fileDir, tab in fileDirs.iteritems():
        mem.allListMaster = {}
        mem.allListMasterMaster = {}

        decayFlag = True if 'decay' in mem.fileDir else False

        mem.files = glob('./'+mem.fileDir+'/*.txt')
        isotopesTmp = [f.split('/')[-1].split('.txt')[0] for f in mem.files]
        mem.elements = []
        mem.isotopes = []

        if len(mem.files) > 0:
            logTxtAndPrint('%i isotopes in %s directory\n'%(len(mem.files), mem.fileDir))
        else:
            logTxtAndPrint('%s directory is empty.\n'%(mem.fileDir))

        # initialize progress bar
        popup = Toplevel()
        popup_lab = Label(popup, text='Reading files:  /'+mem.fileDir, width=50)
        popup_lab.grid(row=0, column=0)
        status = 0
        status_var = DoubleVar()
        progress_bar = ttk.Progressbar(popup, variable=status_var, maximum=100)
        try:
            statusInc = 100.0/len(mem.files)
            progress_bar.grid(row=1, column=0, sticky=EW)
        except:
            pass

        trackerA = 0
        trackerB = 0

        reactionsDict = {}

        deleteFiles = []
        skipFiles = ['CHANGELOG.txt','README.txt','n1.txt']
        for file in mem.files:
            for sf in skipFiles:
                if sf in file:
                    deleteFiles.append(file)
        for file in deleteFiles:
            mem.files.remove(file)

        for file in mem.files:
            if options.verbose:
                print 'Reading %s'%file
            popup.update()
            a = file.split('/')[-1].split('.')[0]
            b = re.split('(\d+)',a)
            if len(b) > 1:
                mem.elements.append(b[0])
                if len(b) > 2:
                    c = b[0]+'-'+''.join(b[1:-1])
                else:
                    c = b[0]+'-'+b[1]
            else:
                mem.elements.append(a.split('nat')[0])
                c = b[0].split('nat')[0]+'-nat'

            mem.isotopes.append(c)
            mem.fileDict_all[mem.note1.index(tab)][file] = c
            f=open(file)
            mem.lines = f.readlines()

            # if there is an error with an NJOY processed file, keep going
            try:
                reactions = make_dict(mem.fileDir)
                mem.isotopesMTdict[c] = reactions
            except:
                continue

            if not (mem.fileDir == 'endfvii-atomic' or 'decay' in mem.fileDir):
                mem.isotopesMTdict[c].append('absorption')

            # create local dictionary of reactions and isotopes containing that reaction. For quicker rank-order analysis
            for rx in reactions:
                if not rx in reactionsDict.keys():
                    reactionsDict[rx]=[[c,file]] # why the extra []? not sure.
                else:
                    reactionsDict[rx].append([c,file])

            status += statusInc
            status_var.set(status)

        # create master dictionary of reactions with tabs as keys
        mem.reactionsDict[tab]=reactionsDict

        mem.elements = np.unique(mem.elements)

        mem.elementsIsotopesDict = {}
        tmpList=[]

        for e in mem.elements:
            for i in mem.isotopes:
                if i.split('-')[0] == e:
                    tmpList.append(i)
            mem.elementsIsotopesDict[e] = tmpList
            tmpList=[]
        for widget in fileDirs[mem.fileDir].winfo_children():
            widget.destroy()
        elIndex=-1
        if 'multigroup' in mem.fileDir:
            flag_pOrM = True
        else:
            flag_pOrM = False

        for element in mem.elements:
            # if there's an error while making a periodic table button, keep going
            try:
                btn = makePeriodicTableButton(element, tab, \
                  mem.elementsIsotopesDict, periodicTableDict, mem.isotopesMTdict, decayFlag, mem.fileDir, flag_pOrM)
            except:
                continue
        popup.destroy()

        az = [chr(i) for i in range(65,91)]
        mem.azDict = {}
        for letter in az: mem.azDict[letter] = []
        for el in mem.elements: mem.azDict[el[0]].append(el)

        if mem.fileDir.split('/')[-1] == 'neutrons':
            mem.masterTracker['pointwise'] = mem.allListMasterMaster
        else:
            mem.masterTracker[mem.fileDir.split('/')[-1]] = mem.allListMasterMaster

#===========================================================
#   Get everything ready to plot
#===========================================================
def getInfo2(isotopes, flag_pOrM, allFlag):
    multiLibFlag = False

    if allFlag and len(mem.plotTitles)>0:
        mem.plotX = {}
        mem.plotY= {}
        for title in mem.plotTitles:
            logTxtAndPrint('%s removed from plot. \n'%(title))
        mem.plotTitles = []
        mem.plotTherm = []
        mem.plotFiss = []
        mem.plotFus = []
        mem.plotUser = []

    else:
        for isotope in isotopes:
            iso, mtDescript = isotope.split()
            mem.element=iso.split('-')[0]
            file = mem.fileDict_all[mem.note1.index("current")].keys()[mem.fileDict_all[mem.note1.index("current")].values().index(iso)]
            mtID = mtdic2.keys()[mtdic2.values().index(mtDescript)]
            with open(file, 'r') as f:
                mem.lines = f.readlines()
            mem.Title_var.set(iso)
            if mem.particle.get() == 1:  # For X-ray files, MF==23
                mem.MF_var.set('23')
            else:
                mem.MF_var.set('3')
            mem.MT_var.set(str(mtID))
            mem.MT_descript.set(mtDescript)
            for k,v in mem.downloadList.iteritems():
                if v.get() != 0:
                    type = 'm' if 'multigroup' in file else 'p'
                    mem.multiLibDic[k].append(
                        [isotope,
                        flag_pOrM,
                        mem.MF_var.get(),
                        mem.MT_var.get(),
                        mem.MT_descript.get(),
                        file])
                    multiLibFlag = True

            if multiLibFlag:
                continue
            else:
                addPlot(flag_pOrM, allFlag)
        if multiLibFlag:
            makeMultiLib()


def makeMultiLib():
    for k,v in mem.multiLibDic.iteritems():
        if len(v)>0:
            for vi in v:
                mem.Title_var.set('%s_%s'%(k, vi[0]))
                flag_pOrM = vi[1]
                mem.MF_var.set(vi[2])
                mem.MT_var.set(vi[3])
                mem.MT_descript.set(vi[4])
                with open(vi[5], 'r') as f:
                    mem.lines = f.readlines()
                addPlot(flag_pOrM, False)




#===========================================================
#   A tiny function to find index of nearest value of array
#===========================================================
def find_nearest(array, value):
    return (np.abs(array-value)).argmin()

#===========================================================
#   Create data for plots
#===========================================================
def addPlot(flag_pOrM, allFlag, mixFlag=False):#event=None, flag_pOrM=False):
    Av = 6.022e23 # Avogado's number
    if mixFlag:
        x = mem.mixDict['eBins']
        y = mem.mixXS
        title = mem.mixTitle

    elif flag_pOrM:
        if mem.MT_var.get() =='999':
            gTot = makeMultiGroup(mem.lines, mem.MF_var.get(),1)
            gSca = makeMultiGroup(mem.lines, mem.MF_var.get(),2)
            x = np.array(gTot.eBins_n)
            y = np.array(gTot.xs) - np.array(gSca.xs)

        else:
            g = makeMultiGroup(mem.lines, mem.MF_var.get(), mem.MT_var.get())
            x = np.array(g.eBins_n)
            y = np.array(g.xs)

        title = mem.Title_var.get()+' '+mem.MT_descript.get()+'_mg'

    else:
        if mem.MT_var.get() == '999':    # compute absorption as total - scatter
            if mem.particle.get()==0:
                x1,yTot = get_section_data(mem.MF_var.get(),1)
                x1,ySca = get_section_data(mem.MF_var.get(),2)
            else:
                x1,yTot = get_section_data(mem.MF_var.get(),501)
                x2,ySca_coh = get_section_data(mem.MF_var.get(),502)
                x3,ySca_inc = get_section_data(mem.MF_var.get(),504)
                # scattered data are on a coarser resolution than total
                # fix this by fitting scattered data to a cubic spline with
                # the same energy resolution as total
                ySca_coh = cs(x2,ySca_coh)(x1)
                ySca_inc = cs(x3,ySca_inc)(x1)
                ySca = np.array(ySca_inc)

            x = np.array(x1)
            y = np.array(yTot)-np.array(ySca)

        else:
            x,y = get_section_data(mem.MF_var.get(),mem.MT_var.get())
            x = np.array(x)
            y = np.array(y)

        title = mem.Title_var.get()+' '+mem.MT_descript.get()

    N = periodicTableDict[mem.element][4]/periodicTableDict[mem.element][3]*Av

    if mem.microMacro.get() == 1 and len(periodicTableDict[mem.element])== 5:
        y = N*y*1.e-24

    elif mem.microMacro.get() == 2 and len(periodicTableDict[mem.element])== 5:
        y = 1/(N*y*1.e-24)

    elif mem.microMacro.get() == 3 and len(periodicTableDict[mem.element])== 5:
        y = N*y*1.e-24/periodicTableDict[mem.element][4]

    logTxtAndPrint('\nNumber density for %s = %.3e\n'%(mem.element, N))

    e_th    = 0.025 # thermal neutron energy
    e_fiss  = 1.0e6 # fission neutron energy
    e_fus   = 14.e6 # fusion neutron energy
    e_user  = mem.Select_E1.get()

    if x.min() <= e_th <= x.max():
        xs_th   = y[find_nearest(x,e_th)]
    else:
#         xs_th = 'Not tabulated'
        xs_th = 0.0

    if x.min() <= e_fiss <= x.max():
        try:
            xs_fiss   = y[find_nearest(x,e_fiss)]
        except:
            xs_fiss = x[-1]
    else:
#         xs_fiss = 'Not tabulated'
        xs_fiss = 0.0

    if x.min() <= e_fus <= x.max():
        try:
            xs_fus   = y[find_nearest(x,e_fus)]
        except:
            xs_fus = x[0]
    else:
#         xs_fus = 'Not tabulated'
        xs_fus = 0.0

    if x.min() <= e_user <= x.max():
        xs_user   = y[find_nearest(x,e_user)]
    else:
#          xs_fus = 'Not tabulated'
        xs_user = 0.0

    if not (len(x)==len(y)):
        x = x[-len(y):]

    mem.plotX[title] = x
    mem.plotY[title] = y
    mem.plotTitles.append(title)
    mem.plotTherm.append(xs_th)
    mem.plotFiss.append(xs_fiss)
    mem.plotFus.append(xs_fus)
    mem.plotUser.append(xs_user)

    logStr =  '%s added to plot.\n'%(title)
    logStr += 'Atomic mass: %.3f\n'%(periodicTableDict[mem.element][3])
    logStr += u'Density (g/cm\u00b3): %.3f\n'%(periodicTableDict[mem.element][4])
    logStr += '-'*25+'\n'
    logTxtAndPrint('%s\n'%(logStr))

#===========================================================
#   A function to remove plots from the list
#===========================================================
def delPlot(event):
    mem.plotTherm = mem.plotTherm[:-1]
    mem.plotFiss = mem.plotFiss[:-1]
    mem.plotFus = mem.plotFus[:-1]
    mem.plotUser = mem.plotUser[:-1]

    title = mem.plotTitles[-1]
    mem.plotX.pop(mem.plotTitles[-1])
    mem.plotY.pop(mem.plotTitles[-1])
    mem.plotTitles = mem.plotTitles[:-1]
    logTxtAndPrint('%s removed from plot.\n'%(title))
    plotMe2()

#===========================================================
#   A  function to save plot data to file
#===========================================================
def savePlotData():
    print 'saving plot data...'
    f = open('zz.txt','w')

    t = '*%45s'%(' ')
    for i in mem.columnText[1:]:
        t += '%20s'%(i)
    f.write(t+'\n')

    for j, row in enumerate(mem.cell_text):
        t = '%45s'%(mem.plotTitles[j])
        rt = [float(i) for i in row[1:]]
        for i in rt:
            t += '%20.3e'%i
        f.write(t+'\n')
    f.close()

#===========================================================
#   Create cross section plots, tabular data, and pie charts
#===========================================================
def plotMe2(event=None, allFlag=False, saveFlag=False, saveDir='figs'):
    # same as findNearest function... why repeat it?
    def fn(array,value):
        return (np.abs(array-value)).argmin()

    # Average sigma, integrated
    def AvgSigmaInt(E, sig, phi , Elow, Ehigh):
        Elow = fn(E, Elow)
        Ehigh = fn(E, Ehigh)
        E = E[Elow:Ehigh]
        sig = sig[Elow:Ehigh]
        phi = phi[Elow:Ehigh]
        intTop = np.trapz(sig*phi, E)
        intBtm = np.trapz(phi, E)
        if intBtm == 0:
            return 0.
        else:
            return intTop/intBtm

    # resonance integral
    def ResInt(E, sig, Elow, Ehigh):
        Elow = fn(E, Elow)
        Ehigh = fn(E, Ehigh)
        E = E[Elow:Ehigh]
        sig = sig[Elow:Ehigh]
        intTop = np.trapz(sig/E, E)
        return intTop

    # set allFlag if only 1 isotope is selected
    # mixFlag = [True if '$' in i else False for i in mem.plotTitles]
    # mixFlag = True if True in mixFlag else False

    lineStyles = ['-.' if ('$' or 'mg') in i else '-' for i in mem.plotTitles]
    lineStyles = []
    mgFlag = False
    for i in mem.plotTitles:
        if ('$' or 'mg') in i:
            mgFlag = True
        lineStyles.append('-')

    # if not '-.' in lineStyles:
    if not mgFlag:
        a = [title.split()[0] for title in mem.plotTitles]
        if len(np.unique(a))==1 and len(a)==len(mem.allListMaster[title.split()[0]]):
            allFlag=True

    pl.cla()
    if allFlag:
        fig = pl.figure(1)
        figImg = fig.set_size_inches(15., 10.)
        ax1 = pl.subplot2grid((2, 4), (0, 0), colspan=2)

    else:
        fig = pl.figure(1)
        figImg = fig.set_size_inches(15., 10)
        ax1 = pl.subplot2grid((1,2), (0, 0))

    # fig1, ax1 = pl.subplots()
    e_th    = 0.025 # thermal neutron energy
    e_fiss  = 1.0e6 # fission neutron energy
    e_fus   = 14.e6 # fusion neutron energy
    e_user = mem.Select_E1.get()
    en = [e_th,  e_fiss, e_fus, e_user] # thermal, fission, fusion neutron energies

    mem.minimumX = []
    mem.maximumX = []
    mem.minimumY = []
    mem.maximumY = []

    mem.cell_text = []
    mem.rowText = []
    plotOrder = []
    unitsDict2 = {0:"barns",
                 1:"1/cm",
                 2:"cm",
                 3:"cm^2/g"}

    eCust = 'N/A' if e_user == 0. else '%3.2e\n(eV)'%(e_user)
    mem.columnText = ['',
        r'$E_{th}$',
        r'$E_{fiss}$',
        r'$E_{fus}$',
        # r'$E_{cust}$',
        # '%3.2e eV'%(mem.Select_E1.get()),
        eCust,
        r'$\overline{\sigma}$',
        r'$\int{\sigma_{res}}$',
        r'$\overline{\phi_{fiss}(E)}$']


    lineWidth = 2
    legendFontSize = 14 # useful for journal article plots
    if allFlag:
        legendFontSize = 10

    yTmp = np.array([0])
    snapData = []
    for i in range(0,len(mem.plotX.keys())):
        plotOrder.append(mem.plotTitles[i])
        x = mem.plotX[mem.plotTitles[i]]
        y = mem.plotY[mem.plotTitles[i]]
        xMin = x.min()
        xMax = x.max()
        yMin = y.min()
        yMax = y.max()
        mem.minimumX.append(xMin)
        mem.maximumX.append(xMax)
        mem.minimumY.append(yMin)
        mem.maximumY.append(yMax)
        title = mem.plotTitles[i]#.split(' ')[2:]
        therm = mem.plotTherm[i]
        fiss = mem.plotFiss[i]
        fus = mem.plotFus[i]
        user = mem.plotUser[i]

        # IAEA spectrum-averaged and resonance integral cross sections
        E0  = 2.53e-2
        E1  = 1.e-5
        E2  = 10.
        E3  = 0.5
        E4  = 1.e5
        Ef  = 1.35e6
        Ef1 = 1.e3
        Ef2 = 2.e7
        E14 = 1.4e7
        Phi_m = (x/E0**2)*np.exp(-x/E0)
        Phi_f = np.sqrt(x/Ef)/Ef*np.exp(-x/Ef)
        avgSigma = 2./np.sqrt(np.pi)*AvgSigmaInt(x, y, Phi_m, E1, E2)
        resInt = ResInt(x, y, E3, E4)
        fissSpecAvg = AvgSigmaInt(x, y, Phi_f, Ef1, Ef2)

        colorIdx = i%len(mem.c) # loop over mem.c colors
        for e in en:
            ax1.axvline(x=e,color='darkcyan',linestyle='--')

        l = lineStyles[i]
        if yTmp.tolist()==y.tolist():
            l = '--'
        yTmp = y

        if 'mg' in title or '$' in title:
            ax1.loglog(x,y, linestyle=l, c=mem.c[colorIdx],lw=lineWidth, drawstyle='steps', label=str(i+1)+'. '+title)
            ax1.set_xscale('log')
            ax1.set_yscale('log')
            mixFlag = False

        else:
            ax1.loglog(x,y, linestyle=l, c=mem.c[colorIdx],lw=lineWidth, label=str(i+1)+'. '+title)
        mem.cell_text.append(['%-2s'%(i+1),'%5.3e'%(therm),'%5.3e'%(fiss),'%5.3e'%(fus),'%5.3e'%(user), '%5.3e'%(avgSigma),'%5.3e'%(resInt), '%5.3e'%(fissSpecAvg)])
        mem.rowText.append(title)

        # cursor snap
        snapData.append((x,y))

    unitsDict = {0:"Cross Section (barns)",
                 1:"Cross Section (1/cm)",
                 2:"Mean free path (cm)",
                 3:"Mean free path (cm^2/g)"}


    tableFormat = [1.05, 0.0, 1.2, 1.0] if allFlag else [1.05, 1.-0.06*(len(mem.plotTitles)+1), 1.35, 0.06*(len(mem.plotTitles)+1)]
    the_table = ax1.table(cellText=mem.cell_text,
                        loc='right',
                        colWidths=[0.15,0.6,0.5,0.5,0.6,0.5,0.5,0.6],
                        colLabels=mem.columnText,
                        colLoc='center',
                        rowLoc='center',
                        cellLoc='center',
                        bbox=tableFormat)

    for (row, col), cell in the_table.get_celld().items():
        cell.set_linewidth(0)
        if (col == 0):
            cell.set_text_props(fontproperties=FontProperties(weight='bold'))
        if row == 0:
            cell._text.set_fontsize(14)
            cell.set_facecolor('lightgrey')
            cell.set_height(.08)
        if col==0 and row>=1:
            colorIdx = row%len(mem.c)
            cell._text.set_color(mem.c[row-1])


    the_table.auto_set_column_width([-1,-1,-1,-1,-1,-1,-1,-1])
    the_table.auto_set_font_size(True)
    # the_table.set_fontsize(18)
    the_table.scale(1, 0.25)
    ax1.set_xlim([0.6*np.array(mem.minimumX).min(), 1.4*np.array(mem.maximumX).max()])
    ax1.set_ylim([max(0.6*np.array(mem.minimumY).min(), 1.e-7), 1.4*np.array(mem.maximumY).max()])
    ax1.grid(which='both', linestyle='--',color='0.8')


    ax1.set_ylabel(unitsDict[mem.microMacro.get()], fontsize=legendFontSize)
    ax1.set_xlabel("Energy (eV)", fontsize=legendFontSize)
    ax1.set_title(unitsDict[mem.microMacro.get()], fontsize=legendFontSize+2)
    ax1.tick_params(labelsize=legendFontSize-2)
    if len(mem.plotTitles) >5:
        ax1.legend(loc=3, prop={'size':legendFontSize/1.5})
    else:
        ax1.legend(loc=3, prop={'size':legendFontSize})

    # Create pie charts of cross section values at each energy only if plotting all x-sects (allFlag=True)
    pieList = [mem.plotTherm, mem.plotFiss, mem.plotFus, mem.plotUser]
    titleList = mem.columnText[1:]

    # no need for thermal energy for photatomic data
    if mem.particle.get()==1:
        pieList = [mem.plotFiss, mem.plotFus, mem.plotUser]
        titleList = mem.columnText[2:]

    if allFlag:
        for i, eGroup in enumerate(pieList):
            if not sum(eGroup[1:])==0.:
                # for x-rays, absorption is not included, so ok to loop over everything
                if mem.particle.get()==1:
                    sizes = np.array(eGroup[1:])/sum(eGroup[1:])*100
                    labels = plotOrder[1:]
                    ax = pl.subplot2grid((2, 4), (1,i))
                # for neutrons, exclude the absorption
                else:
                    sizes = np.array(eGroup[1:-1])/sum(eGroup[1:-1])*100
                    labels = plotOrder[1:-1]
                    ax = pl.subplot2grid((2, 4), (1,i))
                sl = sorted(zip(sizes,labels))
                sizes = [si for si in sizes if si>=0.05]
                sl = sorted([(si,li) for si,li in sl if si>=0.05], reverse=True)
                ax.pie(sorted(sizes,reverse=True))
                ax.axis('equal')
                ax.legend(prop={'size':legendFontSize}, bbox_to_anchor=[0.75+i*0.1, 0.15], labels=['%1.1f %%, %s' % (s,l) for s,l in sl])
                ax.set_title(titleList[i])

    isotopes = [i.split()[0] for i in mem.plotTitles]

    if len(np.unique(isotopes))>1:
        pl.suptitle('Various Isotopes')
    else:
        pl.suptitle(mem.plotTitles[0].split()[0])
    logTxtAndPrint('Processing %s\n'%(mem.plotTitles[0].split()[0]))

    # snap
    if mem.Cursor_var.get() and not allFlag:
        cursor = SnapToCursor(ax1, snapData, fig.canvas, plotOrder, (0.65, 0.9), 'xs')
        pl.connect('motion_notify_event', cursor.mouse_move)


    if saveFlag==True:
        if not os.path.exists(saveDir):
            os.mkdir(saveDir)
        name = []
        for p in mem.plotTitles:
            name.append(p.replace(' ','_'))
        name = '__'.join(name)
        if len(np.unique(isotopes)) == 1 and len(isotopes) >1:
            isotope = np.unique(isotopes)[0]
            el = isotope.split('-')[0]
            if mem.allListMasterMaster[el][isotope]==mem.plotTitles:
                name = '%s_All'%(isotope)
        pdf = mem.savePDF('%s/%s.pdf'%(saveDir, name))
        pdf.savefig(figImg)
        pdf.close()
    else:
        pl.show()


def plotDecay(saveFlag=False, saveDir='figs'):
    '''
    This function will plot decay data
    '''
    ps = 3. # scale plot markers by this value
    pDict = {
        'Gamma': ['*', 150*ps, 0.5, 'r'],
        'X-rays': ['x', 40*ps, 0.8, 'k'],
        'Beta-': ['<', 90*ps, 0.5, 'b'],
        'EC, Beta+': ['o', 60*ps, 0.5, 'y'],
        'Alpha':  ['o', 100*ps, 0.5, 'g'],
        'Neutrons': ['>', 100*ps, 0.5, 'm'],
        'SF' : ['1', 100*ps*ps, 0.5, 'orange'],
        'Protons' : ['+', 100, 0.5, 'navy'],
        'e- (Auger, conversion)': ['o', 20*ps, 0.5, 'orange'],
        'e- anti-neutrinos': ['|', 60*ps, 0.5, 'violet'],
        'e- neutrinos': ['-', 60*ps, 0.5, 'slategray']}

    fig = pl.figure('Radioactive Decay', figsize=(15,8))
    # ax = fig.add_subplot(111)
    ax = pl.subplot2grid((1,2), (0, 0))
    ax.set_axisbelow(True)

    i = 0
    pl.xscale('log')
    pl.yscale('log')
    ms = 70
    snapData = []
    titles = []
    minX, maxX, minY, maxY, cellText, tableColors = ([] for i in range(6))
    avgE, maxE = ({} for i in range(2))
    xData = None
    for k,v in mem.decaySummary.iteritems():
        if k in ['Half Life', 'Decay Modes', u'T-\xbd','isotope']:
            continue
        xData = [vi[0] for vi in v]
        yData = [vi[1] for vi in v]
        totFrac = np.sum([j[1] for j in mem.decaySummary[k]])
        avgE[k] = np.sum([j[0]*j[1]/totFrac for j in mem.decaySummary[k]])
        maxE[k] = max(i[0] for i in mem.decaySummary[k])
        cellText.append(['%s'%(k), '%.3e'%(avgE[k]), '%.3e'%(maxE[k])])
        minX.append(min(xData))
        minY.append(min(yData))
        maxX.append(max(xData))
        maxY.append(max(yData))
        snapData.append((xData, yData))
        titles.append(k)
        colorIdx = i%len(mem.c)
        c = mem.c[colorIdx]
        tableColors.append(pDict[k][3])
        pl.scatter(xData, yData, marker=pDict[k][0], s=pDict[k][1], alpha=pDict[k][2], facecolors=pDict[k][3], edgecolors='k', label=k)

        i += 1

    if not xData:
        return None
    # Table
    columnText = ['Particle', r'$E_{avg}\ (eV)$', r'$E_{max}\ (eV)$']
    tableFormat = [1.05, 1-0.06*(len(titles)+1), 1.2, 0.06*(len(titles)+1)]
    the_table = ax.table(cellText=cellText,
                        loc = 'right',
                        colWidths = [0.25, 0.25, 0.25],
                        colLabels = columnText,
                        colLoc = 'center',
                        rowLoc = 'center',
                        cellLoc = 'center',
                        bbox = tableFormat)

    for (row, col), cell in the_table.get_celld().items():
        cell.set_linewidth(0)
        if (row == 0) or (col == 0):
            cell.set_text_props(color='k')
        if row==0:
            cell.set_facecolor('lightgrey')
        if col == 0 and row >= 1:
            colorIdx = i%len(mem.c)
            cell._text.set_color(tableColors[row-1])
            cell.set_text_props(fontweight='bold')

    the_table.auto_set_column_width([-1, -1, -1])
    the_table.auto_set_font_size(True)
    the_table.scale(4,2)

    val, unit = getDecayUnits(mem.decaySummary[u'T-\xbd'][0])
    fs = 20
    fw = 'bold'
    pl.legend(loc=3, fontsize=fs/2., markerscale=1)
    if val>1000:
        titleString = u'%s\n$T_{\xbd} = %.3e\ %s$'%(mem.decaySummary['isotope'], val, unit)
    else:
        titleString = u'%s\n$T_{\xbd} = %.3f\ %s$'%(mem.decaySummary['isotope'], val, unit)
    pl.suptitle(titleString, fontsize=fs, fontweight=fw)
    pl.xlabel('Energy (eV)', fontsize=fs-4)
    pl.ylabel('Probability per decay', fontsize=fs-4)
    ax.set_xlim([0.1*min(minX), 6*max(maxX)])
    ax.set_ylim([0.1*min(minY), 6*max(maxY)])
    pl.grid(which='both', linestyle='--',color='0.9')
    pl.subplots_adjust(wspace=0, hspace=0)

    if mem.Cursor_var.get():
        cursor = SnapToCursor(ax, snapData, fig.canvas, titles, (0.1, 0.9), 'prob')
        pl.connect('motion_notify_event', cursor.mouse_move)

    if saveFlag==True:
        if not os.path.exists(saveDir):
            os.mkdir(saveDir)
        name = '%s_decay'%(mem.decaySummary['isotope'])
        pdf = mem.savePDF('%s/%s.pdf'%(saveDir,name))
        pdf.savefig(fig)
        pdf.close()
    else:
        pl.show()

#===========================================================
#  This function unzips NNDCD ENDF files, and extracts each
#  isotope as an individual file from the tape files.
#  NNDC tar/gz files are messy. Very. Messy.
#===========================================================
def nndcParse(ver):

    def getIsotopeFromID(g):
        # this function converts an ENDF isotope ID into human-readable isotopes
        # e.g.  6.61640+ 4 --> Dy-164
        #       4.009000+4 --> Zr-90
        #       5.01100+ 3 --> B-11
        b = int(g.strip().split('+')[1])
        end = (2,1) if b==3 else (1,2)
        id = g.strip().split('+')[0].replace('.','')[:4] if b==3 else g.strip().split('+')[0].replace('.','')[:5]
        el = newDict.keys()[newDict.values().index(int(id[:end[1]]))]
        a = id[end[1]:].lstrip('0') if not len(id[end[1]:].lstrip('0'))==0 else 'nat'
        return '%s%s'%(el,a)

    def renameEndf(name):
        subdir = '/'.join(name.split('/')[:3])
        newName = name.split('/')[-1]
        el = newName.split('_')[1]
        iso = newName.split('_')[-1].split('.')[0].lstrip('0')
        if 'neutrons' or 'decay' in name:
            os.system('mv %s %s/%s%s.txt'%(name,subdir,el,iso))
        else:
            os.system('mv %s %s/%snat.txt'%(name,subdir,el))
        return

    # slimmed down version of periodTableDict
    newDict = {}
    for k,v in periodicTableDict.iteritems():
        newDict[k]=v[0]

    if not ver in nndcDict.keys():
        sys.exit()

    if ver in ['ENDF/B-VII.1', 'ENDF/B-VIII.0']:
        for subdir in nndcDict[ver][1:]:
            files = glob(subdir+'/*.endf')
            for f in files:
                renameEndf(f)
    else:
        files = glob(nndcDict[ver][1])
        tally=0
        isotopesAll = []
        crazy = []
        isotopeID = []
        for file in files:
            if '307' in file: # skip 307 for ENDF/B-III...what about others?
                continue
            f = open(file,'r')
            lines = f.readlines()
            s = '1451    1'
            firstLines = [i for i,l in enumerate(lines) if s in l]
            for j,i in enumerate(firstLines):
                try:
                    isotope = getIsotopeFromID(lines[i][:11].lstrip())
                except:
                    continue
                crazy.append([file,i,isotope])
                if len(isotope)>0:
                    isotopesAll.append(isotope)
                    tally+=1

                    if options.verbose: print file, isotope

                    verDir = ver.replace('/','_').replace('-','_')
                    if not os.path.isdir(verDir+'/parsed'):
                        os.mkdir(verDir+'/parsed')
                    outFileName = verDir+'/parsed/'+isotope+'.txt'
                    if os.path.isfile(outFileName):
                        for letter in ['a','b','c','d','e','f','g']:
                            if os.path.isfile(verDir+'/parsed/'+isotope+'_'+letter+'.txt'):
                                print 'ver',letter,'exists'
                                continue
                            else:
                                outFileName = verDir+'/parsed/'+isotope+'_'+letter+'.txt'
                                break
                    fOut = open(outFileName,'w')
                    # fOut.write(lines[0]) # NJOY requires this to be first line of the deck
                    fOut.write('ENDF TAPE HEADER FOR NJOY                                          300 0  0    0\n')
                    if not i==firstLines[-1]:
                        for l in range(firstLines[j], firstLines[j+1]):
                            fOut.write(lines[l])
                            tmp = lines[l]
                        tmp = tmp[:-13]+'-1'+tmp[-11:]
                        fOut.write(tmp) # NJOY requires this to be last line of the deck
                    else:
                        for l in range(firstLines[j], len(lines)):
                            fOut.write(lines[l])
                            tmp = lines[l]
                        tmp = tmp[:-13]+'-1'+tmp[-11:]
                        # fOut.write(tmp) # NJOY requires this to be last line of the deck
                        # fOut.write(' 0.000000+0 0.000000+0          0          0          0          0  -1 0  0    0\n')
                    fOut.close()
                f.close()
        print 'total # isostopes',tally

#===========================================================
#   These functions control the download select buttons
#===========================================================
def addMe(root, processOnly=False):
    def download(url,destination):
        f = urllib2.urlopen(url)
        totalSize = int(f.info().getheader('Content-Length').strip())/1e6
        print 'Total File Size: %i MB'%(totalSize)
        increment = 10 # download in this many MB chunks
        CHUNK = increment * 1048576
        bytesSoFar = 0
        with open(destination,'wb') as o:
            while True:
                chunkStart = time.time()
                chunk = f.read(CHUNK)
                chunkStop = time.time()
                dt = chunkStop - chunkStart # time to download increment # of MB
                if not chunk:
                    break
                o.write(chunk)
                # inform user how much time remains
                bytesSoFar += len(chunk)
                bytesRemaining = int(totalSize-int(bytesSoFar/1e6))
                timeRemaining = float(bytesRemaining)/(increment/(dt/60))
                print '%i MB remaining (approx. %.1f minutes)'%(bytesRemaining, timeRemaining)

    def setReDownload(i):
        yesNo_var.set(i)
        priorDownload.destroy()

    yesNo_var = IntVar()
    yesNo_var.set(0)
    for k in sorted(nndcDict.keys()):
        if mem.downloadList[k].get() == 1:
            if urllib2.urlopen(nndcDict[k][0]).geturl() == nndcDict[k][0]:

                kDir = k.replace('/','_').replace('-','_')
                dest = nndcDict[k][0].split('/')[-1].replace('-','_')

                if os.path.isdir(kDir) and not processOnly:
                    priorDownload = Toplevel()
                    priorDownload.geometry('+400+200')
                    pdl_text1 = '%s has already been downloaded.\n Subdirectories include:\n'%(k)
                    for i, dir in enumerate(dirDict2[k], start=1):
                        if os.path.isdir(dir):
                            pdl_text1+= '\t%i) %s\n'%(i, dir)
                    pdl_text1+='Download and overwrite?'
                    pdl_label1 = Label(priorDownload, text=pdl_text1, justify='left')
                    pdl_label1.grid(row=0, column=0, columnspan=2)
                    yes_btn = Button(priorDownload, text='Yes, download and overwrite', command = lambda dlVar = 1: setReDownload(dlVar))
                    no_btn1 = Button(priorDownload, text='keep current downloaded library and clean up', command = lambda dlVar = 2: setReDownload(dlVar))
                    no_btn2 = Button(priorDownload, text='No, keep current downloaded library', command = lambda dlVar = 0: setReDownload(dlVar))
                    yes_btn.grid(row=1, column=0)
                    no_btn1.grid(row=1, column=1)
                    no_btn2.grid(row=1, column=2)
                    # wait for user feedback
                    root.wait_window(priorDownload)

                elif not os.path.isdir(kDir):
                    yesNo_var.set(1)

                if yesNo_var.get() == 2:
                    processOnly = True

                if processOnly:
                    yesNo_var.set(1)
                    mv('%s/%s'%(kDir, dest),'.')

                if yesNo_var.get() == 1:
                    if os.path.isdir(kDir):
                        sp.check_call(['rm','-rf', kDir])

                elif yesNo_var.get() == 0:
                    continue

                if not os.path.isdir(kDir):
                    sp.check_call(['mkdir', kDir])

                if not options.demo:
                    # download tar file from NNDC
                    print 'downloading %s data file from\n %s\n'%(k,nndcDict[k][0])
                    start = time.time()
                    # mem.urlFile.retrieve(nndcDict[k][0], '/'.join([kDir, dest]))
                    if processOnly:
                        mv(dest, kDir)
                    else:
                        download(nndcDict[k][0], '/'.join([kDir, dest]))
                    end = time.time()
                    print '%s downloaded in %7.3f sec'%(dest,(end-start))
                # untar the file in its respective fileDirectory
                sp.check_call(['tar','-xvf', '/'.join([kDir, dest]), '-C',kDir+'/'])
                print k,'has been dowloaded. Post-processing in progress...'
                logTxtAndPrint('\n%s has been dowloaded. Post-processing in progress...\n'%k)
                nndcParse(k)
                logTxtAndPrint('\n%s has been successfully post-processed.\n'%k)
                load = np.sort(nndcDict.keys()).tolist().index(k)
                mem.verSelect.set(load)
                update_files()

            else:
                print k,'cannot be downloaed from NNDC at this time.'
                logTxtAndPrint('\n%s hcannot be downloaed from NNDC at this time.\n'%k)


def fileChecker(root):
    numFilesFound=0
    for k,v in dirDict2.iteritems():
        for i, value in enumerate(v):
            checkDir = './'+value+'/*.txt'
            if os.path.isdir('./'+value):
                files = glob(checkDir)
                if len(files) > 0:
                    logTxtAndPrint('%i files detected in %s subdirectory.\n'%(len(files), value))
                    numFilesFound += len(files)
                    mem.verSelect.set(i)
                    mem.verSelectRad[k].configure(state='normal')

    if numFilesFound ==0:
        logTxtAndPrint('Must download ENDF files first.\n')
        root.bind('<Escape>', close)
        root.mainloop()

#===========================================================
#   Mass process all files in directory and return
#   arrays of total cross section at thermal, fiss, and fus
#   energies
#===========================================================
def preCheck(file,MF,MT):
    f=open(file)
    lines = f.readlines()
    v = [l[70:75] for l in lines]
    n = len(v)
    mfmtstr = '%2s%3s' % (MF, MT)  # search string

    tally=0
    for iv in v:
        if iv == mfmtstr:
            tally+=1
            f.close()
            return file

    if tally==0:
        lines=[]
        return False

#===========================================================
#   Create a dictionary of available MT and MF combinations within this ENDF file
#===========================================================
def master_list(tab):
    import operator
    Av = 6.022e23

    mem.isotopeList = []
    mem.thermalList = []
    mem.fissionList = []
    mem.fusionList = []
    mem.thermalSort = []
    mem.thermalDict = {}
    mem.fissionDict={}
    mem.fusionDict={}
    mem.eUserDict={}
    absFlag = False
    mem.mixDict = {}

    logTxtAndPrint('Searching for isotopes with %s reaction...'%(mem.Select_MT.get()))

    mem.preCheckedFiles=[]
    if mem.Select_MT.get() == 'absorption':
        absFlag = True

    if mem.Select_MT.get() == 'radioactive_decay':
        mem.note1.select(2)

    if mem.note0.index(mem.note0.select()) == 1:
        whichTab = mem.note1.index(1)

    else:
        whichTab = mem.note1.index("current")

    popup = Toplevel()
    popupText = StringVar()
    popupText.set('Analyzing isotopes with %s reaction'%mem.Select_MT.get())
    popup_lab2 = Label(popup, textvariable=popupText, width=50)
    status2 = 0
    status_var2 = DoubleVar()
    progress_bar2 = ttk.Progressbar(popup, variable=status_var2, maximum=100)
    try:
        statusInc1 = 100.0/len(mem.fileDict_all[whichTab].keys())
        popup_lab2.grid(row=0, column=0)
        progress_bar2.grid(row=1, column=0, sticky=EW)
    except:
        pass

    mem.preCheckedFiles = []
    for i in mem.reactionsDict[mem.periodicTableTabs[whichTab][0]][mem.Select_MT.get()]:
        mem.preCheckedFiles.append(i[1])

    logTxtAndPrint('%i isotopes found for %s.\n'%(len(mem.preCheckedFiles), mem.Select_MT.get()))

    status = 0
    tally = 0

    try:
        statusInc2 = 100.0/len(mem.preCheckedFiles)
    except:
        popup_lab2 = Label(popup, text='Select ENDF library (with multigroup data) first', width=50)
        popup_lab2.grid(row=0, column=0)
        popup.update()
        return

    for i, file in enumerate(mem.preCheckedFiles):
        start = time.time()
        popupText.set('Analyzing file %i/%i with %s'%(i+1, len(mem.preCheckedFiles),mem.Select_MT.get()))

        popup.update()
        isotope = mem.fileDict_all[whichTab][file]
        mem.element = isotope.split('-')[0]
        # This is probably a good place to call another function that assimilates all of the radioactive decay data into a large data structure
        if mem.Select_MT.get() == 'radioactive_decay':
            if mem.decayOpt['state'] == 'disabled':
                mem.decayOpt['state'] = 'normal'
                mem.decaySortName['state'] = 'normal'
                mem.decaySortEnergy['state'] = 'normal'
                mem.decaySortOccurrence['state'] = 'normal'
                mem.decay_analysis_opt.set(sorted(mem.decayTypeDict.keys())[0])
            dirTmp = '/'.join(file.split('/')[:-1])
            try:
                decaySummary = get_decay_data(isotope, dirTmp, miniParse=True)
                mem.tHalfDict[isotope] = decaySummary[u'T-\xbd'][0]
                mem.decayModeDict[isotope] = decaySummary['Decay Modes']
                del decaySummary['Decay Modes']
                del decaySummary[u'T-\xbd']

                # assemble all of the decay data
                for k,v in mem.decayTypeDict.iteritems():
                    try:
                        # each key (isotope) = most probable particle energy
                        v[isotope] = sorted(decaySummary[k], key=lambda x: x[1])[-1]
                    except:
                        pass
            except:
                pass
            status2 += statusInc2
            status_var2.set(status2)
            continue

        else:
            with open(file, 'r') as f:
                mem.lines = f.readlines()

            if mem.lines[1].split()[4] == '-1':
                if mem.Select_MT.get() =='absorption':
                     gTot = makeMultiGroup(mem.lines, 3,1)
                     gSca = makeMultiGroup(mem.lines, 3,2)
                     x = np.array(gTot.eBins_n)
                     y = np.array(gTot.xs) - np.array(gSca.xs)

                else:
                     g = makeMultiGroup(mem.lines, 3,mtdic2.keys()[mtdic2.values().index(mem.Select_MT.get())])
                     x = np.array(g.eBins_n)
                     y = np.array(g.xs)

            else:
                if absFlag:    # compute absorption as total-scatter
                    x1,yTot = get_section_data(3,1)
                    x1,ySca = get_section_data(3,2)
                    x=np.array(x1)
                    y=np.array(yTot)-np.array(ySca)

                elif mem.particle.get()==1:
                    x,y = get_section_data(23,mtdic2.keys()[mtdic2.values().index(mem.Select_MT.get())])
                else:
                    x,y = get_section_data(3,mtdic2.keys()[mtdic2.values().index(mem.Select_MT.get())])

        x = np.array(x)
        y = np.array(y)
        if not 'eBins' in mem.mixDict.keys():
            mem.mixDict['eBins'] = x
        mem.mixDict[isotope] = y

        if mem.microMacro.get() == 1 and len(periodicTableDict[mem.element])==5:
            N = periodicTableDict[mem.element][4]/periodicTableDict[mem.element][3]*Av
            y = N*y*1.e-24
        elif mem.microMacro.get() == 2 and len(periodicTableDict[mem.element])==5:
            N = periodicTableDict[mem.element][4]/periodicTableDict[mem.element][3]*Av
            yTmp = N*y*1.e-24
            y = 1./yTmp
        elif mem.microMacro.get() == 3 and len(periodicTableDict[mem.element])==5:
            N = periodicTableDict[mem.element][4]/periodicTableDict[mem.element][3]*Av
            y = N*y*1.e-24/(periodicTableDict[mem.element][4])

        e_th    = 0.025 # thermal neutron energy
        e_fiss  = 1.0e6 # fission neutron energy
        e_fus   = 14.e6 # fusion neutron energy
        e_user = mem.Select_E1.get()
        en = [e_th,  e_fiss,e_fus, e_user] # thermal, fission, fusion neutron energies


        if x.min() <= e_th <= x.max():
            xs_th   = y[find_nearest(x,e_th)]
        else:
            #         xs_th = 'Not tabulated'
            xs_th = 0.0

        if x.min() <= e_fiss <= x.max():
            try:
                xs_fiss   = y[find_nearest(x,e_fiss)]
            except:
                xs_fiss = x[-1]
        else:
            xs_fiss = 'Not tabulated'
            xs_fiss = 0.0

        if x.min() <= e_fus <= x.max():
            try:
                xs_fus   = y[find_nearest(x,e_fus)]
            except:
                xs_fus = x[-1]
        else:
            # xs_fus = 'Not tabulated'
            xs_fus = 0.0

        if x.min() <= e_user <= x.max():
            xs_user   = y[find_nearest(x,e_user)]
        else:
            # xs_fus = 'Not tabulated'
            xs_user = 0.0

        if not (len(x)==len(y)):
            x = x[-len(y):]

        mem.thermalDict[isotope] = xs_th
        mem.fissionDict[isotope] = xs_fiss
        mem.fusionDict[isotope] = xs_fus
        mem.eUserDict[isotope] = xs_user

        mem.isotopeList.append(isotope)
        mem.thermalList.append(xs_th)
        mem.fissionList.append(xs_fiss)
        mem.fusionList.append(xs_fus)
        mem.thermalSort.append(xs_user)

        status2 += statusInc2
        status_var2.set(status2)

    mem.decayTypeDict['Half Life'] = mem.tHalfDict
    popup.destroy()
    # mem.logText.insert(INSERT, 'Analysis complete.\n')
    logTxtAndPrint('Analysis complete.\n')

    analysisType = 'decay' if mem.Select_MT.get() == 'radioactive_decay' else 'xs'
    displayAnalysis()

def getDecayUnits(d):
    # convert seconds to hrs, days, months, years
    for k,v in mem.halfLifeScale.iteritems():
        if d >= v[0] and d < v[1]:
            if k == 'sec':
                return d, k
            else:
                return d/v[0], k


def displayAnalysis(event=None):
    def dictToArray(dict,num):
        if num == 2:
            outList = []
            a = sorted(dict.items(), key=operator.itemgetter(0))
            elements = np.unique([i[0].split('-')[0] for i in a]).tolist()
            dicElements = {}
            for el in elements:
                dicElements[el] = [i for i in a if i[0].split('-')[0]==el]
                dicElements[el] = sorted(dicElements[el], key = lambda x: x[1], reverse=True)
            a = [vi for k,v in sorted(dicElements.iteritems()) for vi in v]

        elif num == 1:
            a = sorted(dict.items(), key=operator.itemgetter(num), reverse=True)

        elif num == 0:
            a = sorted(dict.items(), key=operator.itemgetter(num))
        a = np.array(a)
        return a

    def dictToTuples(dict):
        tuples = []
        for k,v in dict.iteritems():
            if type(v) == float:
                tuples.append((k, v))
            elif len(v) == 2:
                tuples.append((k, v[0], v[1]))
        return tuples

    if mem.Select_MT.get() == 'radioactive_decay':
        if mem.decay_analysis_opt.get()=='Half Life':
            mem.decaySortEnergy['text'] = u'T-\xbd'
            mem.decaySortOccurrence.configure(state='disabled')
        else:
            mem.decaySortEnergy['text'] = 'Energy'
            mem.decaySortOccurrence.configure(state='normal')

        results = sorted(dictToTuples(mem.decayTypeDict[mem.decay_analysis_opt.get()]), key=lambda x: x[mem.decaySort.get()])
        mem.resultsText.delete('1.0', END)
        if len(results[0]) == 2:
            mem.resultsText.insert(INSERT, u'%15s %15s\n\n'%('Isotope',u'T-\xbd'))
        elif len(results[0]) == 3:
            mem.resultsText.insert(INSERT, '%15s %15s %15s\n\n'%('Isotope','E (eV)', 'Fraction'))
        for i, r in enumerate(results, start=1):
            if len(r) == 2:
                value, unit = getDecayUnits(r[1])
                if value >=1000. and unit=='years':
                    mem.resultsText.insert(INSERT, '%4i %10s %15.3e %8s\n'%(i, r[0], value, unit))
                else:
                    mem.resultsText.insert(INSERT, '%4i %10s %15.3f %8s\n'%(i, r[0], value, unit))
            if len(r) == 3:
                mem.resultsText.insert(INSERT, '%4i %10s %15.3e %15.3e\n'%(i, r[0], r[1], r[2]))


    else:# mem.Select_MT.get() == 'xs':
        xs_list = [mem.thermalDict, mem.fissionDict, mem.fusionDict, mem.eUserDict]
        xs_listS = ['thermalDictSorted', 'fissionDictSorted', 'fusionDictSorted', 'eUserDictSorted']
        e_list  = ['thermal', 'fiss', 'fus', 'user']

        mem.resultsToPrint = dictToArray(xs_list[mem.eRange.get()], mem.sortBy.get())

        mem.resultsText.delete('1.0', END)
        if mem.microMacro.get() == 1:
            mem.resultsText.insert(INSERT, u"   Isotope      \u03a3(1/cm)\n\n")
        elif mem.microMacro.get() == 2:
            mem.resultsText.insert(INSERT, u"   Isotope      Mean_Free_Path(cm)\n\n")
        elif mem.microMacro.get() == 3:
            mem.resultsText.insert(INSERT, u"   Isotope      \u03a3/\u03c1(cm^2/g)\n\n")
        else:
            mem.resultsText.insert(INSERT, u"   Isotope      \u03c3(barns)\n\n")

        for i in range(len(mem.resultsToPrint)):
            if isinstance(mem.resultsToPrint[i,1], basestring):
                mem.resultsText.insert(INSERT, '%3i %-10s %8.4e\n' %(i+1, mem.resultsToPrint[i,0],float(mem.resultsToPrint[i,1]))) #'%2s%3s' % (MF, MT)
            else:
                mem.resultsText.insert(INSERT, '%3i %-10s %8.4e\n' %(i+1, mem.resultsToPrint[i,0],float(mem.resultsToPrint[i,1]))) #'%2s%3s' % (MF, MT)
            if mem.sortBy.get() in [0, 2] and i< len(mem.resultsToPrint)-1:
                if mem.resultsToPrint[i+1,0].split('-')[0] !=mem.resultsToPrint[i,0].split('-')[0]:
                    mem.resultsText.insert(INSERT, '\n')

#===========================================================
#   Process all raw ENDF files with NJOY
#===========================================================
def run_njoy(root):

    def setReprocess(i):
        reprocess.set(i)
        priorProcess.destroy()

    def NJOY(file):
        with open(file, 'r') as f:
            line = f.readlines()[1]
        matID = line[66:70].rstrip()
        matName = file.split('/')[-1].split('.')[0]
        ver = file.split('/')[-2]
        tag = 'NJOY processed '+ver+' '+matName+' '+matID
        print tag
        os.system('mv '+file+' ./working/tape20')
        f=open('./working/njoy.in','w')
        f.write("reconr\n")
        f.write("20 21\n")
        f.write("'"+tag+"'/\n")
        f.write(matID+" 1/\n")
        f.write(".01 0.0/ 1% linearization\n") # processed at 300 deg K
        f.write("'"+tag+"'/\n")
        f.write("0/\n")
        f.write("broadr\n")
        f.write("20 21 22\n")
        f.write(matID+" 1/\n")
        f.write("0.001 /\n")
        f.write("%.1f /\n"%(mem.njoy_temp.get()))
        # f.write("300 /\n")
        f.write("0 /\n")
        f.write("groupr\n")
        f.write("20 22 0 23\n")
        # NJOY 2016
        f.write(matID+" %i 0 3 3 1 1 1 1\n"%(njoy_ign_dict[mem.njoy_ign.get()]))
        # f.write(matID+" 3 3 3 3 1 1 1\n")
        #f.write(matID+" 10 3 3 3 1 1 1\n") # 100 group
        f.write(matName+"\n")
        f.write("%.1f /\n"%(mem.njoy_temp.get()))
        # f.write("300 /\n")
        f.write("1.e10 /\n")
        f.write("3 /\n")
        f.write("3 455\n")
        f.write("5 455\n")
        f.write("0 /\n")
        f.write("0 /\n")
        f.write("stop\n")
        f.close()
        os.chdir('./working')
        os.system('./njoy < njoy.in | tee -a %s'%(logFile))
        os.system('mv tape22 '+file.split('/')[-1])
        os.system('mv tape23 multigroup/'+file.split('/')[-1])
        os.chdir('../')
        os.system('pwd')

    processList = []
    reprocess = IntVar()
    reprocess.set(0)

    for key,value in mem.downloadList.iteritems():
        if value.get()==1 and key!='ver_xRay':
            processList.append(key)
    if not os.path.isdir('./working/multigroup'):
        os.mkdir('./working/multigroup')

    start = time.time()

    for value in processList:
        newValue = dirDict2[value][mem.particle.get()]
        checkDir = './'+newValue+'/*.txt'
        myDir = './'+newValue+'/'
        logFile = 'njoy_%s.log'%(value.replace('/','_').replace('-','_'))
        if not os.path.isfile(logFile):
            os.system('touch ./working/%s'%(logFile))

        files = glob(checkDir)
        with open(files[0],'r') as f:
            line = f.readline()

        if 'NJOY processed' in line:
        # if os.path.isdir('./'+newValue+'/multigroup'):
            priorProcess = Toplevel()
            priorProcess.geometry('+300+200')
            pdl_text1 = '%s has already been processed with NJOY.\nProcess again and overwrite?'%(newValue)
            pdl_label1 = Label(priorProcess, text=pdl_text1, justify='left')
            pdl_label1.grid(row=0, column=0, columnspan=2)
            yes_btn = Button(priorProcess, text='Yes, process and overwrite', command = lambda dlVar = 1: setReprocess(dlVar))
            no_btn = Button(priorProcess, text='No, keep current NJOY-processed files', command = lambda dlVar = 0: setReprocess(dlVar))
            yes_btn.grid(row=1, column=0)
            no_btn.grid(row=1, column=1)
            # wait for user feedback
            root.wait_window(priorProcess)
            if reprocess.get()==0:
                continue


        if reprocess.get() ==1:
            addMe(root, processOnly=True)
            os.mkdir('./'+newValue+'/multigroup')
            files = glob(checkDir)
            # mem.logText.insert(INSERT, 'NJOY is starting to process %s\n'%(value))
            logTxtAndPrint('NJOY is starting to process %s.\n'%(value))
            for file in files:
                if options.verbose: print 'NJOY working on %s'%(file)
                NJOY(file)

            end = time.time()
            os.system('mv ./working/*.txt '+myDir)
            os.system('mv ./working/multigroup/*.txt '+myDir+'/multigroup')
            tNJOY = (end - start)/60.
            logTxtAndPrint('NJOY processing of %s complete.\n\n'%(value))

        elif os.path.isdir('./'+newValue):
            if not os.path.isdir('./'+newValue+'/multigroup'):
                os.mkdir('./'+newValue+'/multigroup')
            logTxtAndPrint('NJOY is starting to process %s.\n'%(value))
            files = glob(checkDir)
            for file in files:
                if options.verbose: print 'NJOY working on %s'%(file)
                NJOY(file)
            end = time.time()
            os.system('mv ./working/*.txt '+myDir)
            os.system('mv ./working/multigroup/*.txt '+myDir+'/multigroup')
            tNJOY = (end - start)/60.
            logTxtAndPrint('NJOY processing of %s complete\n\n'%(value))
        else:
            pass
    if value==processList[-1]:
        logTxtAndPrint('Loading the %s library...\n'%(value))
        load = np.sort(nndcDict.keys()).tolist().index(value)
        mem.verSelect.set(load)
        update_files()
    sys.stdout = sys.__stdout__


#===========================================================
#   Save results of rank order analysis to file
#===========================================================
def saveText(*args):
    def fileNameFix(a):
        return a.replace('(','').replace(')','').replace(',','_')

    units = ['barns','1/cm','cm','g/cm^2']
    unit = units[mem.microMacro.get()]

    mem.testMe = mem.resultsText.get("1.0",END)

    uHead = [unit+'(Th)', unit+'(Fis)', unit+'(Fus)', unit+'(%8.2e)'%(mem.Select_E1.get())]
    mem.header = '%-18s%-18s%-18s%-18s%-18s\n'%('Isotope', uHead[0], uHead[1], uHead[2], uHead[3])

    for i in np.arange(1,100,1):
        if not os.path.isfile('./out'+str(i)+'.txt'):
            outFileName = 'out'+str(i)+'.txt'
            break
    f = open(outFileName,'w')
    f.write(mem.header)
    f.close()

    for i in range(4):
        mem.eRange.set(i)
        displayAnalysis()
        mem.a = mem.resultsText.get("1.0",END)
        mem.a = mem.a.split()[2:]
        mem.name = mem.a[1::3]
        mem.data = mem.a[2::3]

        mem.name = [n.encode("utf-8") for n in mem.name]
        mem.data = [n.encode("utf-8") for n in mem.data]

        if mem.sortBy.get()==0:
            if i==0:
                mem.all = np.column_stack((mem.name,mem.data))
            else:
                mem.all = np.column_stack((mem.all,mem.data))
        else:
            if i==0:
                mem.all = np.column_stack((mem.name,mem.data))
            else:
                mem.all = np.column_stack((mem.all,mem.data))

    with open(outFileName,'a') as f:
        np.savetxt(f,mem.all, delimiter="", fmt="%-18s")
    f.close()
    logTxtAndPrint('Analysis results written to file %s.\n'%(outFileName))

def reset():
    for ver in nndcDict.keys():
        mem.downloadList[ver].set(0)

def selectAll():
    for ver in nndcDict.keys():
        mem.downloadList[ver].set(1)

#===========================================================
#   Close All
#===========================================================
def close(event=None):
    clearAll()
    sys.exit() # if you want to exit the entire thing

def clearAll(event=None):
    mem.plotX={}
    mem.plotY={}
    mem.plotTitles=[]
    mem.plotTherm=[]
    mem.plotFiss=[]
    mem.plotFus=[]
    mem.plotUser=[]
    pl.close('all')
    mem.resultsText.delete('1.0', END)
    mem.logText.delete('1.0', END)
    logTxtAndPrint('Plots cleared.\n')
    mem.batchFile = None
    for k,v in mem.multiLibDic.iteritems():
        mem.multiLibDic[k] = []

def makeMixList(*args):
    master_list(mem.tab02)
    for i in range(10):
        row = i+3
        mem.weightEntryVar = StringVar()
        mem.weightEntry = Entry(mem.tab2_frameL, textvariable=mem.weightEntryVar, bg='aliceblue')
        mem.weightEntry.configure(width=6)
        mem.Select_isotope = StringVar()

        mem.isoLabelVar = StringVar()
        mem.isoLabel = Label(mem.tab2_frameL, textvariable=mem.isoLabelVar,bg='aliceBlue')

        mem.Select_isotope_opt = makeAlphabetButton(mem.tab2_frameL,mem.elementsIsotopesDict, row, mem.Select_isotope,mem.isotopeList,mem.isoLabelVar,mem.weightEntryVar)

        mem.weightEntry.grid(column=0, row=row, padx=3, pady=3)

        mem.isoLabel.configure(width=10)
        mem.isoLabel.grid(column=1,row=row,padx=3,pady=3)
        # keep track of widgets as per: https://goo.gl/aLSkKF
        mem.entryVars[i] = mem.weightEntryVar
        mem.selectVars[i] = mem.Select_isotope
        mem.entry[i] = mem.weightEntry
        mem.options[i] =  mem.Select_isotope_opt
        mem.labelVars[i] = mem.isoLabelVar
    mem.plotMix_btn = Button(mem.tab2_frameL, text='Plot mix', command=plotMix)
    mem.plotMix_btn.grid(row=i+4, column=0)

def clearMix():
    mem.plotMixCount = 0
    for i in range(10):
        mem.entryVars[i].set('')
        mem.selectVars[i].set('')
        mem.labelVars[i].set('')

def clearPlotMix():
    clearMix()
    mem.ax.clear()
    mem.canvas_tab02.draw()

def plotMix():
    def normalizeList(a):
        a = np.array(a)
        return a/sum(a)

    def saveData(x,y,title):
        np.savetxt('zz_'+mem.Select_MT.get()+'_'+title+'.txt', np.vstack((x,y)).T, '%1.4e')

    mixLength = 10
    mem.mixXS = np.zeros(len(mem.mixDict['eBins']))
    title = ''
    weights = []
    weightsTmp = []

    if mem.percentType_var.get() == 1:
        massSum = 0.
        for i in range(mixLength):
            if not (mem.entry[i].get() == ''):
                el = mem.selectVars[i].get().split('-')[0]
                elMass = periodicTableDict[el][3]
                massSum += float(mem.entry[i].get())*elMass
                weightsTmp.append(float(mem.entry[i].get())*elMass)
        print massSum
        for i, w in enumerate(weightsTmp):
            weights.append(w/massSum)

    else:
        for i in range(mixLength):
            if not (mem.entry[i].get() == ''):
                weights.append(float(mem.entry[i].get())/100.)
        weights = normalizeList(weights)


    for i in range(mixLength):
        if not (mem.entry[i].get() == ''):
            mem.mixXS += weights[i]*np.array(mem.mixDict[mem.selectVars[i].get()])
            if len(title) == 0:
                title += r'%5.3f$\times$%s '%(weights[i], mem.selectVars[i].get())
            else:
                title += r'%5.3f$\times$%s '%(weights[i], mem.selectVars[i].get())

    title += mem.Select_MT.get()
    mem.mixTitle = title
    e_th    = 0.025 # thermal neutron energy
    e_fiss  = 1.0e6 # fission neutron energy
    e_fus   = 14.e6 # fusion neutron energy
    e_user = mem.Select_E1.get()
    en = [e_th,  e_fiss, e_fus, e_user] # thermal, fission, fusion neutron

    addPlot(False,False,True)

    mem.plotMixSave_btn = Button(mem.tab2_frameL, text='Save', command=lambda \
                            x=mem.mixDict['eBins'],
                            y=mem.mixXS,
                            title=title:saveData(x,y,title))
    mem.plotMixSave_btn.grid(row=i+4, column=1)
    mem.clearMix = Button(mem.tab2_frameL, text='Clear', command=clearMix)
    mem.plotClear_btn = Button(mem.tab2_frameL, text='Clear Plots', command=clearPlotMix)
    mem.clearMix.grid(row=i+5, column=0)
    mem.plotClear_btn.grid(row=i+6, column=0)
    mem.plotMixCount+=1


def initializeVariables():
    #===========================================================
    #   Main code
    #===========================================================
    mem.lines = []
    mem.isotopeTitle = ''
    mem.fileToGet=''
    mem.MTMF_dict={}
    mem.fileDict_p = {} # fileDict for pointwise files
    mem.fileDict_m = {} # fileDict for multigroup fileDirs
    mem.fileDict_decay = {} # fileDict for multigroup fileDirs
    # list of fileDicts corresponding to tab1 (point), tab2 (multigroup)
    mem.fileDict_all = [mem.fileDict_p, mem.fileDict_m, mem.fileDict_decay]
    mem.labels = ['MF','MT']
    mem.MTlist=['']
    mem.MFlist=['']
    mem.MTlist_v=[] # verbose version that includes actual ENDF interpretations
    mem.MFlist_v=[] # verbose version that includes actual ENDF interpretations
    mem.dlList=[]
    mem.dlTotal={}
    mem.linksDict={}
    mem.urlFile = urllib.URLopener()
    mem.el_isotope = []
    mem.elements = []
    mem.isotopes = []
    mem.elementsIsotopesDict={}
    mem.isotopesMTdict = {}
    mem.remoteFile=''
    mem.thermalDict={}
    mem.fissionDict={}
    mem.fusionDict={}
    mem.eUserDict={}
    mem.isotopeList=[]
    mem.thermalList=[]
    mem.fissionList=[]
    mem.fusionList=[]
    mem.userList=[]
    mem.preCheckedFiles=[]
    mem.plotX={}
    mem.plotY={}
    mem.plotTitles=[]
    mem.plotTherm=[]
    mem.plotFiss=[]
    mem.plotFus=[]
    mem.plotUser=[]
    mem.dlList1 = []
    mem.dlList2 = []
    mem.dlList3 = []
    mem.mixTitle=''
    mem.plotYtmp = []
    mem.mixSum = 0.0
    mem.plotMixCount = 0
    mem.c = ['b','r','g','y','m','orange','lightgreen','tan','grey','lightgrey','k']
    mem.reactionsDict={}

    # decay radiation dicts
    mem.tHalfDict = {}
    mem.alphaDict = {}
    mem.xrayDict = {}
    mem.gammaDict = {}
    mem.betaDict = {}
    mem.neutronDict = {}
    mem.sfDict = {}
    mem.ecDict = {}
    mem.protonDict = {}
    mem.eDict = {}
    mem.antiNeutrinoDict = {}
    mem.neutrinoDict = {}
    mem.decayModeDict = {}
    mem.tHalfList = []
    mem.alphaList = []
    mem.xrayList = []
    mem.gammaList = []
    mem.multiLibDic = {}

    mem.decayTypeDict = {
        'Gamma': mem.gammaDict,
        'Beta-': mem.betaDict,
        'EC, Beta+': mem.ecDict,
        'Alpha': mem.alphaDict,
        'Neutrons': mem.neutronDict,
        'SF' : mem.sfDict,
        'Protons' : mem.protonDict,
        'e- (Auger, conversion)': mem.eDict,
        'X-rays': mem.xrayDict,
        'e- anti-neutrinos': mem.antiNeutrinoDict,
        'e- neutrinos': mem.neutrinoDict}

    mem.halfLifeScale = {
        'sec':[0., 60.],
        'min':[60., 3600.],
        'hrs':[3600., 86400.],
        'days':[86400., 2.592e6],
        'months':[2.592e6, 3.1536e7],
        'years':[3.1536e7, 1.e99]}


def makeWidgets(root):
    initializeVariables()
    mem.TITLE  = tkFont.Font(family='Garamond',size=int(30*mem.scale_root),weight='bold')
    mem.AUTHOR = tkFont.Font(family='Garamond',size=int(18*mem.scale_root),weight='bold')
    mem.GREEK = tkFont.Font(family='Georgia',size=int(15*mem.scale_root))
    mem.HDG1 = tkFont.Font(family='Helvetica',size=int(14*mem.scale_root))
    mem.DATA = tkFont.Font(family='Courier',size=int(14*mem.scale_root))
    mem.MENUBUTTON = tkFont.Font(family='Helvetica',size=int(12*mem.scale_root))#,'bold')

    mem.note0 = ttk.Notebook(root, width=int(mem.rootWidth*0.6), height=mem.rootHeight)
    mem.tab01 = Frame(mem.note0)#, width=int(rootWidth*scale_root))
    mem.tab02 = Frame(mem.note0)#, width=int(rootWidth*scale_root))
    mem.tab03 = Frame(mem.note0)#, width=int(rootWidth*scale_root))
    Grid.rowconfigure(mem.tab01, 0, weight=1)
    Grid.columnconfigure(mem.tab02, 0, weight=1)
    mem.tab01.columnconfigure(0, weight=1)
    mem.tab02.columnconfigure(0, weight=1)
    mem.tab01.rowconfigure(0, weight=0)
    mem.tab01.rowconfigure(1, weight=1)
    mem.tab01.rowconfigure(2, weight=0)
    mem.tab02.rowconfigure(0, weight=1)
    mem.tab01.grid(row=0, column=0, sticky=NSEW)
    mem.tab02.grid(row=0, column=0, sticky=NSEW)
    mem.note0.add(mem.tab01,text='Download, Analyze, Plot')
    mem.note0.add(mem.tab02,text='Mix')
    mem.note0.pack(side='top', fill='both', expand=True)

    mem.topFrameScale1 = 1
    mem.topFrameScale2 = 1

    mem.top_frame = Frame(mem.tab01, width=mem.rootWidth, height=mem.rootHeight/mem.topFrameScale1, pady=0)#, relief=GROOVE, borderwidth=2)
    mem.download_frame = Frame(mem.top_frame, width=mem.rootWidth, height=mem.rootHeight/mem.topFrameScale1/mem.topFrameScale2, pady=1, relief=GROOVE, borderwidth=2)
    mem.version_frame = Frame(mem.top_frame, width=mem.rootWidth, height=mem.rootHeight/mem.topFrameScale1/mem.topFrameScale2, pady=1, relief=GROOVE, borderwidth=2)
    mem.microMacro_frame = Frame(mem.top_frame, width=mem.rootWidth, height=mem.rootHeight/mem.topFrameScale1/mem.topFrameScale2, pady=1, relief=GROOVE, borderwidth=2)

    mem.isotopes_frame = Frame(mem.tab01,width=int(mem.rootWidth*mem.scale_root), height=mem.rootHeight/mem.topFrameScale1, pady=1)#, relief=GROOVE, borderwidth=2)

    mem.note1 = ttk.Notebook(mem.isotopes_frame, width=int(mem.rootWidth*0.6), height=mem.rootHeight/mem.topFrameScale1)
    mem.tab11 = Frame(mem.note1)
    mem.tab12 = Frame(mem.note1)
    mem.tab13 = Frame(mem.note1)

    # keep track of these tabs for later
    mem.periodicTableTabs = [(mem.tab11,'Pointwise'), (mem.tab12,'Multigroup'), (mem.tab13,'Decay')]

    Grid.rowconfigure(mem.tab11, 0, weight=1)
    Grid.columnconfigure(mem.tab12, 0, weight=1)
    mem.tab11.columnconfigure(0, weight=1)
    mem.tab12.columnconfigure(0, weight=1)
    mem.tab13.columnconfigure(0, weight=1)
    mem.tab11.grid(row=0, column=0, sticky=NSEW)
    mem.tab12.grid(row=0, column=0, sticky=NSEW)
    mem.tab13.grid(row=0, column=0, sticky=NSEW)
    mem.isotopeGrid11 = Frame(mem.tab11)
    mem.isotopeGrid12 = Frame(mem.tab12)
    mem.isotopeGrid13 = Frame(mem.tab13)
    mem.isotopeGrid11.grid(column=0, row=7, sticky=NSEW)
    mem.isotopeGrid11.grid(column=0, row=7, sticky=NSEW)
    mem.note1.add(mem.tab11,text='Pointwise')
    mem.note1.add(mem.tab12,text='Multigroup')
    mem.note1.add(mem.tab13,text='Decay')
    mem.note1.pack(side='top', fill='both', expand=True)
    s = ttk.Style()
    s.configure('TNotebook.Tab', padding=(20, 8, 20, 0))

    mem.user_frame =  Frame(mem.tab01,width=mem.rootWidth, height=50, pady=1,padx=5)
    mem.user_frame.columnconfigure(0, weight=1)
    mem.userE_frame =  Frame(mem.user_frame,width=mem.rootWidth/2, height=50, pady=1)
    mem.userMT_frame =  Frame(mem.user_frame,width=mem.rootWidth/2, height=50, pady=1)

    mem.sortBy_frame = Frame(mem.tab01, width=mem.rootWidth, height=50, pady=7)#, relief=GROOVE, borderwidth=2)
    mem.sortBy_frame.columnconfigure(0, weight=1)
    mem.sortByXS_frame = LabelFrame(mem.sortBy_frame, width=mem.rootWidth/3.8, height=52, padx=20, text='Sort analysis')
    mem.sortByE_frame = LabelFrame(mem.sortBy_frame, width=mem.rootWidth/5, height=52, padx=20, text='Analysis energy')
    mem.sortDecay_frame = LabelFrame(mem.sortBy_frame, width=mem.rootWidth/2.5, height=52, padx=20, text='Decay Analysis')
    # mem.mix_frame = Frame(mem.sortBy_frame, width= 200, height=50, padx=20)

    mem.dataText_frame = Frame(mem.tab01, width=mem.rootWidth/2, height=50, relief=GROOVE, borderwidth=2)
    Grid.rowconfigure(mem.dataText_frame, 0, weight=1)
    Grid.columnconfigure(mem.dataText_frame, 0, weight=1)

    mem.dataTextData_frame = Frame(mem.dataText_frame,  padx=1)
    mem.dataTextData_frame.rowconfigure(0, weight=1)
    mem.dataTextData_frame.columnconfigure(0, weight=1)

    mem.dataTextLog_frame = Frame(mem.dataText_frame, padx=1,relief=GROOVE, borderwidth=2)
    mem.dataTextLog_frame.rowconfigure(0, weight=1)
    mem.dataTextLog_frame.columnconfigure(0, weight=1)

    mem.top_frame.columnconfigure(0, weight=1)
    mem.top_frame.rowconfigure(0, weight=1)
    mem.download_frame.columnconfigure(0, weight=1)
    mem.download_frame.rowconfigure(0, weight=1)
    mem.version_frame.columnconfigure(0, weight=1)
    mem.version_frame.rowconfigure(0, weight=1)
    mem.microMacro_frame.columnconfigure(0, weight=1)
    mem.microMacro_frame.rowconfigure(0, weight=1)
    mem.top_frame.grid(row=0, sticky=EW)
    mem.download_frame.grid(row= 0, column=0)#, sticky=EW)
    mem.version_frame.grid(row= 1, column=0)#, sticky=EW)
    mem.microMacro_frame.grid(row= 2, column=0)#, sticky=EW)

    mem.isotopes_frame.columnconfigure(0, weight=1)
    mem.isotopes_frame.rowconfigure(0, weight=1)
    mem.isotopes_frame.grid(row=1, sticky=NSEW)

    mem.user_frame.grid(row=2)
    mem.userE_frame.grid(row=0,column=0)
    mem.userMT_frame.grid(row=0,column=1)
    mem.sortBy_frame.grid(row=3, sticky=NS)

    mem.sortByXS_frame.grid(row=0, column=0, padx=3)
    mem.sortByXS_frame.grid_propagate(0)
    mem.sortByE_frame.grid(row=0, column=1, padx=3)
    mem.sortByE_frame.grid_propagate(0)
    mem.sortDecay_frame.grid(row=0, column=2, padx=3)
    mem.sortDecay_frame.grid_propagate(0)

    mem.dataText_frame.grid(row=4, column=0, sticky=EW)
    mem.dataTextData_frame.grid(row=0,column=0, sticky=EW)
    mem.dataTextLog_frame.grid(row=0,column=1, sticky=EW)

    mem.downloadList={}
    for i, ver in enumerate(np.sort(nndcDict.keys())):
            mem.downloadList[ver] = IntVar()
            cb = Checkbutton(mem.download_frame, variable=mem.downloadList[ver], text='%-15s'%(ver), font=mem.HDG1)
            cb.grid(row=0, column=i)
            mem.multiLibDic[ver] = []

    mem.selectAll_btn = Button(mem.download_frame,text='Select All',command=selectAll,font=mem.HDG1)
    mem.selectAll_btn.grid(row=1, column=1)
    mem.reset = Button(mem.download_frame,text='Deselect All',command=reset,font=mem.HDG1)
    mem.reset.grid(row=1, column=2)
    mem.get = Button(mem.download_frame,text='Download Now',command=lambda root=root: addMe(root),font=mem.HDG1)
    mem.get.grid(row=1, column=3)
    mem.njoy = Button(mem.download_frame,text='Process with NJOY',command=lambda root=root: run_njoy(root),font=mem.HDG1)
    mem.njoy.grid(row=1, column=4)
    mem.njoy_temp = DoubleVar()
    mem.njoy_temp.set(300.)
    mem.njoy_ign = StringVar()
    mem.njoy_ign.set(njoy_ign_dict.keys()[njoy_ign_dict.values().index(3)])
    mem.njoy_temp_ent = Entry(mem.download_frame, textvariable=mem.njoy_temp, bg='aliceblue',font=mem.HDG1)
    mem.njoy_temp_ent.configure(width=5)
    mem.njoy_temp_ent.grid(row=1, column=5)
    mem.njoy_temp_label = Label(mem.download_frame, text='Temp (K)')
    mem.njoy_temp_label.grid(row=1, column=6)
    mem.njoy_ign_opt = OptionMenu(mem.download_frame, mem.njoy_ign, *njoy_ign_dict.keys())
    mem.njoy_ign_opt.configure(width=15)
    mem.njoy_ign_opt.grid(row=1, column=7)
    mem.njoy_ign_label = Label(mem.download_frame, text='E-group')
    mem.njoy_ign_label.grid(row=1, column=8)


    mem.verSelect = IntVar()
    mem.verSelectRad = {}
    for i, ver in enumerate(np.sort(nndcDict.keys())):
        mem.verSelectRad[ver] = Radiobutton(mem.version_frame, variable=mem.verSelect,text=ver,
            value=i, command=update_files, font=mem.HDG1, state='disabled')
        mem.verSelectRad[ver].grid(row=1, column=i, padx=1, pady=1)
    mem.particle = IntVar()
    mem.particleSelect = {}
    particles = ['neutrons','x-rays']
    for j, ver in enumerate(particles):
        mem.particleSelect[ver] = Radiobutton(mem.version_frame,variable=mem.particle,text=ver,
           value=j, font=mem.HDG1, command=update_files, state='disabled')
        mem.particleSelect[ver].grid(row=1,column=i+1)
        i+=1

    mem.microMacro = IntVar()
    mem.micro = Radiobutton(mem.microMacro_frame, variable=mem.microMacro, text=u'\u03c3', font=mem.GREEK, value=0,padx=15, pady=1)
    mem.macro = Radiobutton(mem.microMacro_frame, variable=mem.microMacro, text=u'\u03a3', font=mem.GREEK, value=1,padx=15, pady=1)
    mem.mfp = Radiobutton(mem.microMacro_frame, variable=mem.microMacro, text=u'(1 / \u03a3)', font=mem.GREEK, value=2,padx=15, pady=1)
    mem.massAtten = Radiobutton(mem.microMacro_frame, variable=mem.microMacro, text=u'(\u03a3/\u03c1)', font=mem.GREEK, value=3,padx=15, pady=1)
    mem.micro.grid(row=0,column=0)
    mem.macro.grid(row=0,column=1)
    mem.mfp.grid(row=0,column=2)
    mem.massAtten.grid(row=0,column=3)
    mem.microMacro.set(0)

    mem.MF_var = StringVar()
    mem.MT_var = StringVar()
    mem.MT_descript = StringVar()
    mem.Title_var = StringVar()

    mem.verSelect.set(0)
    # update_files()

    mem.Select_E1 = DoubleVar()
    mem.Select_lab1 = Label(mem.userE_frame, text='Enter additional energy (in eV): ',font=mem.HDG1)
    mem.Select_E1_ent = Entry(mem.userE_frame, textvariable=mem.Select_E1, bg='aliceblue',font=mem.HDG1)
    mem.Select_E1_ent.configure(width=5)

    mem.Select_MT = StringVar()
    mem.Select_lab2 = Label(mem.userMT_frame, text=' for comparison: ',font=mem.HDG1)
    mem.tmpList = list(mtdic2.values())
    mem.tmpList.sort()
    mem.tmpList.append('')
    mem.Select_MT_opt = OptionMenu(mem.userMT_frame, mem.Select_MT, *mem.tmpList)
    mem.Select_MT_opt.configure(width=10)

    mem.ClearAll_btn = Button(mem.userE_frame, width=10, text='Clear All', command=clearAll,padx=2,font=mem.HDG1)
    mem.Analyze_btn = Button(mem.userMT_frame, width=10, text = 'Analyze All', command=makeMixList,padx=2,font=mem.HDG1)
    # mem.Analyze_btn = Button(mem.userMT_frame, width=10, text = 'Analyze All', command=lambda tab01=mem.tab01:master_list(mem.tab01),padx=2,font=mem.HDG1)
    mem.SaveAnalysis_Btn = Button(mem.userMT_frame, width=13, text='Save Analysis', command=saveText,padx=2,font=mem.HDG1)
    mem.SavePlotData_Btn = Button(mem.userMT_frame, width=15, text='Save Plot Data', command=savePlotData,padx=2,font=mem.HDG1)
    mem.RunBatch_Btn = Button(mem.userMT_frame, width=10, text='Batch', command=lambda root=root: runBatch(root), padx=2, font=mem.HDG1)
    mem.Plot_Btn = Button(mem.userMT_frame, width=10, text='Plot', command=plotMe2, padx=2, font=mem.HDG1)
    mem.Cursor_var = IntVar()
    mem.Cursor = Checkbutton(mem.userMT_frame, variable=mem.Cursor_var, text='Cursor', font=mem.HDG1)

    mem.Select_lab1.grid(row=0,column=0)
    mem.Select_E1_ent.grid(row=0,column=1)
    mem.ClearAll_btn.grid(row=0,column=2)
    mem.Analyze_btn.grid(row=0,column=2)
    mem.SaveAnalysis_Btn.grid(row=0,column=3)
    mem.SavePlotData_Btn.grid(row=0,column=4)
    mem.RunBatch_Btn.grid(row=0, column=5)
    mem.Select_lab2.grid(row=0,column=0)
    mem.Select_MT_opt.grid(row=0,column=1)
    mem.Plot_Btn.grid(row=0, column=6)
    mem.Cursor.grid(row=0, column=7)

    mem.sortBy = IntVar()
    mem.sortByName = Radiobutton(mem.sortByXS_frame, variable=mem.sortBy, value=0,text='by isotope',command=displayAnalysis,font=mem.HDG1)
    mem.sortByXS = Radiobutton(mem.sortByXS_frame, variable=mem.sortBy, value=1,text='by value',command=displayAnalysis,font=mem.HDG1)
    mem.sortByNameThenXS = Radiobutton(mem.sortByXS_frame, variable=mem.sortBy, value=2,text='by isotope, then value',command=displayAnalysis,font=mem.HDG1)
    mem.sortByName.grid(row=0,column=0)
    mem.sortByXS.grid(row=0,column=1)
    mem.sortByNameThenXS .grid(row=0, column=2)

    mem.eRange = IntVar()
    mem.eRangeTh = Radiobutton(mem.sortByE_frame,variable=mem.eRange,value=0,text='Th',command=displayAnalysis,font=mem.HDG1)
    mem.eRangeFiss = Radiobutton(mem.sortByE_frame,variable=mem.eRange,value=1,text='Fiss',command=displayAnalysis,font=mem.HDG1)
    mem.eRangeFus = Radiobutton(mem.sortByE_frame,variable=mem.eRange,value=2,text='Fus',command=displayAnalysis,font=mem.HDG1)
    mem.eRangeUser = Radiobutton(mem.sortByE_frame,variable=mem.eRange,value=3,text='E_User',command=displayAnalysis,font=mem.HDG1)
    mem.eRangeTh.grid(row=0,column=0)
    mem.eRangeFiss.grid(row=0,column=1)
    mem.eRangeFus.grid(row=0,column=2)
    mem.eRangeUser.grid(row=0,column=3)
    mem.sortByName.select()
    mem.eRangeTh.select()

    mem.decaySort = IntVar()
    mem.decay_analysis_opt = StringVar()

    mem.decayOpt = OptionMenu(mem.sortDecay_frame, mem.decay_analysis_opt, *sorted(mem.decayTypeDict.keys()+['Half Life']),command=displayAnalysis)

    mem.decaySortName = Radiobutton(mem.sortDecay_frame,variable=mem.decaySort,value=0,text='Isotope',command=displayAnalysis,font=mem.HDG1)

    mem.decaySortEnergy = Radiobutton(mem.sortDecay_frame,variable=mem.decaySort,value=1,text='Energy',command=displayAnalysis,font=mem.HDG1)

    mem.decaySortOccurrence =    Radiobutton(mem.sortDecay_frame,variable=mem.decaySort,value=2,text='Freq',command=displayAnalysis,font=mem.HDG1)

    mem.decayPlotButton = Button(mem.sortDecay_frame, text='Plot Decay', command=plotDecay)

    mem.decayOpt.configure(width=10, state='disabled')
    mem.decaySortName.configure(width=10, state='disabled')
    mem.decaySortEnergy.configure(width=10, state='disabled')
    mem.decaySortOccurrence.configure(width=10, state='disabled')
    mem.decayPlotButton.configure(width=10, state='disabled')
    mem.decayOpt.grid(row=0,column=0)
    mem.decaySortName.grid(row=0,column=1)
    mem.decaySortEnergy.grid(row=0,column=2)
    mem.decaySortOccurrence.grid(row=0,column=3)
    mem.decayPlotButton.grid(row=0, column=4)


    mem.mix = IntVar()
    mem.mixFrac = DoubleVar()
    mem.mixFrac.set(1.0)

    mem.resultsText = Text(mem.dataTextData_frame, height=20, bg='aliceblue', font=mem.DATA)
    mem.scrollbar1 = Scrollbar(mem.dataTextData_frame, command=mem.resultsText.yview)
    mem.resultsText.configure(yscrollcommand=mem.scrollbar1.set)
    mem.resultsText.grid(row=0, column=0,sticky=NSEW,padx=1)
    mem.scrollbar1.grid(row=0, column=1, sticky=NS)

    mem.logText = Text(mem.dataTextLog_frame, height=20, bg='aliceblue',font=mem.DATA)
    mem.scrollbar2 = Scrollbar(mem.dataTextLog_frame, command=mem.logText.yview)
    mem.logText.configure(yscrollcommand=mem.scrollbar2.set)
    mem.logText.grid(row=0, column=0,sticky=NSEW,padx=1)
    mem.scrollbar2.grid(row=0, column=1, sticky=NS)


    # tab 2 widgets
    mem.tab2_frameL = Frame(mem.tab02, width=500, height=1000, padx=3, relief=GROOVE, borderwidth=2)
    mem.tab2_frameR = Frame(mem.tab02, height=1000, padx=3, relief=GROOVE, borderwidth=5)
    mem.tab2_frameL.pack(side='left', fill='both', expand=True)#side=LEFT)
    mem.tab2_frameR.pack(side='left', fill='both', expand=True)


    mem.fig = Figure()
    mem.ax = mem.fig.add_subplot(111)
    mem.canvas_tab02 = FigureCanvasTkAgg(mem.fig,master=mem.tab2_frameR)
    mem.canvas_tab02.show()
    mem.canvas_tab02.get_tk_widget().pack(side='top', fill='both', expand=1)

    mem.entryVars = {}
    mem.selectVars = {}
    mem.labelVars = {}
    mem.entry = {}
    mem.label = {}
    mem.options = {}

    mem.MixReaction_lab = Label(mem.tab2_frameL, text='Select reaction for mixing: ')
    mem.MixReaction_opt = OptionMenu(mem.tab2_frameL, mem.Select_MT, *mem.tmpList, command=makeMixList)
    mem.MixReaction_opt.configure(width=15)
    mem.MixReaction_lab.grid(row=0, column=0)
    mem.MixReaction_opt.grid(row=1, column=0)
    mem.percentType_var = IntVar()
    mem.percentTypeStr = ['Weight-%', 'Atom Ratio']

    for i, mem.percentType in enumerate(mem.percentTypeStr):
        mem.typeRadio_btn = Radiobutton(mem.tab2_frameL, variable=mem.percentType_var, text=mem.percentType, value=i)
        mem.typeRadio_btn.grid(row=2, column=i)

    # Launch batch mode if input file is provided
    if options.batchFile:
        mem.batchFile = options.batchFile
        runBatch(root)
    else:
        mem.batchFile = None



def runBatch(root):
    from matplotlib.backends.backend_pdf import PdfPages as savePDF
    mem.savePDF = savePDF

    def writeTexIndexStyleFile():
        '''
        This makes a nicely-formatted index for LaTexFiles
        '''
        with open('index_style.ist', 'w') as f:
            f.write('headings_flag 1\n')
            f.write(r'heading_prefix "{\\large\\rmfamily\\bfseries "'+'\n')
            f.write(r'heading_suffix "}\\nopagebreak\n"'+'\n')
            f.write(r'delim_0 " \\dotfill "'+'\n')
            f.write(r'delim_1 " \\dotfill "'+'\n')
            f.write(r'delim_2 " \\dotfill "'+'\n')

    def allXS(plotType, dirs):

        writeTexIndexStyleFile()
        flag_pOrM_master = False
        saveFlag_master = True
        allFlag_master = True
        loopList = []
        log = []
        if plotType == 'decay':
            mem.note1.select(2)
            for k,v in sorted(mem.masterTracker[plotType].iteritems()):
                loopList.append(sorted(v.keys()))
            loopList = [item for sublist in loopList for item in sublist]
            dir = [d for d in dirs if 'decay' in d][0]
            for i, isotope in enumerate(loopList, start=1):
                try:
                    get_decay_data(isotope, dir, miniParse=True)
                    plotDecay(saveFlag=True, saveDir=saveDir)
                except:
                    log.append(isotope)
            for i, iLog in enumerate(log):
                print i, iLog


        elif plotType =='pointwise' or plotType =='multigroup':
            mem.note1.select(0)
            if plotType=='multigroup':
                flag_pOrM_master = True
                mem.note1.select(1)
            for k0, v0 in sorted(mem.masterTracker[plotType].iteritems()):
                for k1, v1 in sorted(mem.masterTracker[plotType][k0].iteritems()):
                    try:
                        getInfo2(v1, flag_pOrM_master, allFlag_master)
                        plotMe2(None, allFlag_master, saveFlag_master, saveDir=saveDir)
                    except:
                        pass
                    clearAll()


    def allAnalyze(lib, dataType, reportType):
        writeTexIndexStyleFile()

        def addBody(s, caption):
            b = r'''
                \begin{center}
                \begin{longtable}[p]{ccc}
                \caption{%s.}\\
                \hyperref[index]{Index} / \hyperlink{page.1}{TOC}\\
                \hline
                \index{%s}
                \textbf{Item} & \textbf{Isotope} & \textbf{$\sigma$ (barns)}\\
                \hline
                \endfirsthead
                \multicolumn{3}{c}
                {\tablename\ \thetable\ -- \textit{Continued from previous page}} \\
                \hline
                \textbf{Item} & \textbf{Isotope} & \textbf{$\sigma$ (barns)}\\
                \hline
                \endhead
                \hline \multicolumn{3}{c}{\textit{Continued on next page}} \\
                \endfoot
                \hline
                \endlastfoot
                '''%(caption[1], caption[1])

            c = r'''
            '''
            for i, si in enumerate(s.split('\n')):
                si = si.split()
                if len(si)>1:
                    c += \
                    r'''
                    %s & %s & %s \\'''%(si[0], si[1], si[2])
                else:
                    c += \
                    r'''\\'''

            c +=r'''
            \end{longtable}
            \end{center}
            \newpage
            '''
            return b+c


        sortByDic = {
            0: 'By isotope',
            1: 'By value',
            2: 'By isotope then by value'}

        eRangeDic = {
            0: 'Thermal',
            1: 'Fission',
            2: 'Fusion',
            3: 'User'}

        if type == 'multigroup':
            mem.note1.select(1)
        else:
            mem.note1.select(0)
        MTlist = sorted([int(mt) for mt in mtdic2.keys() if int(mt) <=107])
        MTlist.append('849')
        MTlist.append('999')
        MTlist = [str(mt) for mt in MTlist]

        top, tail = texTopAndTail(lib, dataType, reportType)

        body = r''''''

        energiesNum = 3 if float(mem.Select_E1.get()) == 0. else 4
        for i, MT in enumerate(MTlist):
            mem.Select_MT.set(mtdic2[MT])
            master_list(mem.tab01)
            for energy in range(energiesNum):
                mem.eRange.set(energy)
                for sortBy in range(3):
                    mem.sortBy.set(sortBy)
                    displayAnalysis()
                    caption = [i, '%s.%s.%s'%(mtdic2[MT].replace('_',' ').capitalize(), eRangeDic[energy], sortByDic[sortBy])]
                    body += addBody(mem.resultsText.get('3.0',END), caption)


        tex = top + body + tail
        outputFileName = 'texAnalysis_%s.tex'%(datetime.datetime.now().strftime("%Y%m%d%H%M"))

        with open(outputFileName, 'w') as f:
            f.write(tex)
        for iTex in range(2): # compile twice
            sp.check_call(['pdflatex', outputFileName])
        return None



    if mem.batchFile == None:
        root.update()
        batchFile = tkFileDialog.askopenfilename()
    else:
        batchFile = mem.batchFile

    with open(batchFile,'r') as f:
        tmp = f.readlines()
        lines = []
        for line in tmp:
            lines.append(line.strip())

    commentLines = [i for (i,l) in enumerate(lines) if l[0]=='#']
    for cl in commentLines[::-1]:
        lines.pop(cl)

    tmp = {}
    for i, (k,v) in enumerate(sorted(dirDict2.iteritems())):
        tmp[k]=i



    libDict = {}
    dataType = None
    reportType = None
    for line in lines:
        if line in dirDict2.keys():
            lib = line
            libDict[line] = {}
            libDict[lib]['save'] = []
        elif 'E=' in line or 'E =' in line:
            libDict[lib]['E'] = float(line.split('=')[-1])
            libDict[lib]['xs'] = []
        elif line == 'save':
            libDict[lib]['save'].append(libDict[lib]['xs'])
            libDict[lib]['xs'] = []
        elif 'all' in line.split()[0]:
            saveDir = 'figs_%s'%(datetime.datetime.now().strftime("%Y%m%d%H%M"))
            if not mem.verSelect.get()==tmp[lib]:
                mem.verSelect.set(tmp[lib])
                update_files()
            dataType = line.split()[-1]
            reportType = 'All Reactions'
            allXS(dataType, dirDict2[lib])
            for iTex in range(2): # compile 2x for hyperrefs
                texBatchPlot(saveDir, lib, dataType, reportType)
        elif 'analysis' in line.split()[0]:
            saveDir = 'figs_%s'%(datetime.datetime.now().strftime("%Y%m%d%H%M"))
            if not mem.verSelect.get()==tmp[lib]:
                mem.verSelect.set(tmp[lib])
                update_files()
            dataType = line.split()[-1]
            reportType = 'Analysis'
            allAnalyze(lib, dataType, reportType)
        else:
            libDict[lib]['xs'].append(line)

    for k,v in libDict.iteritems():
        if not mem.verSelect.get()==tmp[k]:
            mem.verSelect.set(tmp[k])
            update_files()

        logTxtAndPrint('Writing out all available isotopes and reactions for %s.'%(lib))
        with open('zz_all_%s.txt'%(lib.replace('/','_').replace('-','_')), 'w') as f:
            for k_el, v_iso in sorted(mem.allListMasterMaster.iteritems()):
                for iso, xs in sorted(v_iso.iteritems()):
                    for item in sorted(xs):
                        f.write('%s\n'%item)

        saveDir = 'figs_%s'%(datetime.datetime.now().strftime("%Y%m%d%H%M"))
        mem.Select_E1.set(libDict[k]['E'])
        for xsList in libDict[k]['save']:
            for xs in xsList:
                flag_pOrM = False
                allFlag = False
                saveFlag = True
                mem.note1.select(0)
                if 'm' in xs.split():
                    flag_pOrM = True
                    mem.note1.select(1)
                if 'all' in xs.split():
                    allFlag = True
                    el = xs.split('-')[0]
                    iso = xs.split()[0]
                    xs = mem.allListMasterMaster[el][iso]
                    logTxtAndPrint('batch analysis: %s %s'%(lib, xs))
                    getInfo2(xs, flag_pOrM, allFlag)
                    continue
                else:
                    xs = ' '.join(xs.split()[0:2])
                    logTxtAndPrint('batch analysis: %s %s'%(lib, xs))
                    getInfo2([xs], flag_pOrM, allFlag)

            # from matplotlib.backends.backend_pdf import PdfPages as savePDF
            # mem.savePDF = savePDF
            plotMe2(None, allFlag, saveFlag, saveDir=saveDir)
            clearAll()
            mem.Select_E1.set(libDict[k]['E'])
        for iTex in range(2): # compile 2x for hyperrefs
            texBatchPlot(saveDir, lib, dataType, reportType)


def texTopAndTail(endf, dataType, reportType):
        texTop = r'''
        \documentclass{article}
        \usepackage{graphicx}
        \usepackage{pdflscape}
        \usepackage[margin=0.8in]{geometry}
        \usepackage{longtable}
        \usepackage{fancyhdr}
        \usepackage{imakeidx}
        \makeindex[columns=2, title={Alphabetical Index\label{index}},options={-s index_style.ist}, intoc]
        \usepackage{hyperref}
        \hypersetup{
            colorlinks=true,
            linkcolor=blue,
            filecolor=magenta,
            urlcolor=cyan
        }
        \usepackage[utf8]{inputenc}
        \usepackage[T1]{fontenc}
        \usepackage[mmddyyyy]{datetime}
        \fancyhead{}
        \fancyfoot{}
        \fancyhead[CO,CE]{---    Compiled by EXSAN    ---}
        \fancyfoot[C]{EXSAN}
        \fancyfoot[R] {\thepage}
        \begin{document}
        \title{ENDF Cross Section \& Nuclear Data Analysis \\
        \large \url{https://github.com/lanl/EXSAN}}
        \maketitle
        \pagestyle{fancy}
        \thispagestyle{fancy}
        \tableofcontents
        \listoffigures
        \listoftables
        \newpage
        \section{Details}
        \begin{tabular}{ll}
        Date: & \today \\
        ENDF Library: & %s \\
        Data Type: & %s \\
        Report Type: & %s \\
        \end{tabular}\\ \\ \\
        \hyperref[index]{Index} / \hyperlink{page.1}{TOC}
        \newpage
        \section{Plots, Tables, and Data}
        '''%(endf, dataType, reportType)

        texTail = r'''
        \printindex
        \hyperlink{page.1}{TOC}
        \end{document}
        '''
        return texTop, texTail

#===========================================================
# Write LaTeX file containing plots from batch file
#===========================================================

def texBatchPlot(figDir, lib, dataType, reportType):
    '''
    This module writes a LaTeX file that assembles all of the figures
    from a batch run, and adds an alphabetical index with hyperlinks
    at the end.

    Helpful links:
        https://es.overleaf.com/learn/latex/hyperlinks
        https://www.overleaf.com/learn/latex/Indices
        https://texblog.org/2007/11/07/headerfooter-in-latex-with-fancyhdr/
    '''

    def addBody(s):
        '''
        Add lines for an individual figure including caption, index
        labels, and and a hyperlink to the Index. This is useful for
        very large files!
        '''
        contents = s.split('/')[-1].split('.')[0].split('__')
        label = '.'.join(contents)
        contentsCaption = r''''''
        indexList = []
        print s
        for c in contents:
            iso, xs = c.split('_')[0:2]
            if c.split('_') == 3:
                mg = c.split('_')[-1]
            el, mass = iso.split('-')
            contentsCaption += \
            r'''\textsuperscript{%s}%s %s, '''%(mass, el, xs)
            indexList.append(r'''%s %s'''%(iso, xs))

        b = \
        r'''
        \begin{landscape}
        \begin{centering}
        \begin{figure}
          \includegraphics[width=1.0\linewidth]{%s}
          \caption{%s}'''%(s, contentsCaption)

        c = r''''''
        for item in indexList:
            c += \
            r'''
            \index{%s}'''%(item)

        b += c
        b += \
        r'''
          \label{fig:%s}
        \hyperref[index]{Index} / \hyperlink{page.1}{TOC}
        \end{figure}
        \end{centering}
        \end{landscape}
        '''%(label)
        return b


    #===========================================================
    # list of figures output by EXSAN's batch analysis
    #===========================================================
    figs = glob('%s/*'%(figDir))
    outputFileName = 'tex_%s.tex'%(figDir.split('_')[-1])


    #===========================================================
    # Top and bottom of the tex file
    #===========================================================
    top, tail = texTopAndTail(lib, dataType, reportType)


    #===========================================================
    # Body of the tex file, one for each figure
    #===========================================================
    body = r''''''
    for fig in sorted(figs):
        body += addBody(fig)

    #===========================================================
    # Make and write the entire tex file and PDF
    #===========================================================
    tex = top + body + tail

    with open(outputFileName, 'w') as f:
        f.write(tex)
    sp.check_call(['pdflatex', outputFileName])
    return



#===========================================================
# Write LaTeX file for analysis
#===========================================================
def texAnalysis(figDir, lib, datatype, reportType):
    '''
    This module writes a LaTeX file that assembles all of the analyses
    and sort options.

    Helpful links:
        https://texblog.org/2011/05/15/multi-page-tables-using-longtable/
    '''

    def addBody(s):
        '''
        Add lines for an individual figure including caption, index
        labels, and and a hyperlink to the Index. This is useful for
        very large files!
        '''
        contents = s.split('/')[-1].split('.')[0].split('__')
        label = '.'.join(contents)
        contentsCaption = r''''''
        indexList = []
        print s
        for c in contents:
            iso, xs = c.split('_')[0:2]
            if c.split('_') == 3:
                mg = c.split('_')[-1]
            el, mass = iso.split('-')
            contentsCaption += \
            r'''\textsuperscript{%s}%s %s, '''%(mass, el, xs)
            indexList.append(r'''%s %s'''%(iso, xs))

        b = \
        r'''
        \begin{landscape}
        \begin{centering}
        \begin{figure}
          \includegraphics[width=1.0\linewidth]{%s}
          \caption{%s}'''%(s, contentsCaption)

        c = r''''''
        for item in indexList:
            c += \
            r'''
            \index{%s}'''%(item)

        b += c
        b += \
        r'''
          \label{fig:%s}
        \hyperref[index]{Index} / \hyperlink{page.1}{TOC}
        \end{figure}
        \end{centering}
        \end{landscape}
        '''%(label)
        return b


    #===========================================================
    # list of figures output by EXSAN's batch analysis
    #===========================================================
    figs = glob('%s/*'%(figDir))
    outputFileName = 'tex_%s.tex'%(figDir.split('_')[-1])


    #===========================================================
    # Top and bottom of the tex file
    #===========================================================
    top, tail = texTopAndTail(lib, datatype, reportType)


    #===========================================================
    # Body of the tex file, one for each figure
    #===========================================================
    body = r''''''
    for fig in sorted(figs):
        body += addBody(fig)

    #===========================================================
    # Make and write the entire tex file and PDF
    #===========================================================
    tex = top + body + tail

    with open(outputFileName, 'w') as f:
        f.write(tex)
    sp.check_call(['pdflatex', outputFileName])
    return





#===========================================================
# create a symlink for the NJOY2016 executable file
#===========================================================
def linkNjoy():
    if not os.path.isdir('working'):
        os.mkdir('working')
    else:
        pass
    if not os.path.isfile('./working/njoy') and not os.path.islink('./working/njoy'):
        if not options.njoyPath:
            print '\n***Note: NJOY2016 is required to process raw ENDF files'
            print 'If you have an NJOY2016 executable, you can supply it'
            print 'when launching EXSAN:\n\tpython exsan.py -n /Users/user/Desktop/NJOY2016/bin/njoy'
        else:
            os.symlink('%s'%(options.njoyPath), 'working/njoy')
            print 'NJOY soft link created'
        return

#===========================================================
# Update the log file and log Text area of the GUI
#===========================================================
def logTxtAndPrint(s):
    s = s.encode('utf-8')
    mem.logText.insert(END, s)
    # UNIX can't print unicode characters to the terminal? so try this:

    if options.log:
        try:
            Log(s)
        except:
            pass
    elif options.verbose:
        print s.strip()
    else:
        pass

def main():
    if options.log:
        if os.path.isfile('./exsan.log'):
            sp.check_call(['rm', './exsan.log'])
        # sys.stdout = Log()
    else:
        sys.stdout = sys.__stdout__
    banner = "                    _______  ______    _    _   _ \n"
    banner +="                   | ____\ \/ / ___|  / \  | \ | |\n"
    banner +="                   |  _|  \  /\___ \ / _ \ |  \| |\n"
    banner +="                   | |___ /  \ ___) / ___ \| |\  |\n"
    banner +="                   |_____/_/\_\____/_/   \_\_| \_|\n"
    banner +="\n-- The (E)NDF (C)ross (S)ection (An)alysis and visualization application --\n"
    banner +="                             Created by:\n"
    banner +="                    H. Omar Wooten, PhD, DABR\n"
    banner +="                          hasani@lanl.gov\n"
    banner +="            Copyright Los Alamos National Laboratory, 2018 \n"
    banner +="        Los Alamos National Laboratory Computer Code LA-CC-19-002\n"
    print banner
    root = Tk()
    root.title('EXSAN 1.0')
    mem.scale_root = 0.95
    # mem.rootHeight = int(root.winfo_screenheight()*mem.scale_root)
    # mem.rootWidth  = int(root.winfo_screenwidth()*mem.scale_root)
    mem.rootWidth  = 1300
    mem.rootHeight = 850
    root.geometry('{}x{}'.format(mem.rootWidth, mem.rootHeight))
    root.grid_columnconfigure(0, weight=1)
    root.grid_rowconfigure(0, weight=1)

    linkNjoy()
    makeWidgets(root)
    fileChecker(root)

    # Keyboard shortcuts.
    root.bind('<Control-Escape>', close)
    root.bind('a',addPlot)
    root.bind('d',delPlot)
    root.bind('p',plotMe2)
    root.bind('c',clearAll)

    # for testing, automatically load ENDF 8
    if options.auto:
        mem.verSelect.set(6)
        update_files()
        # getInfo2(demoIsotopes, False, False)
        # mem.Cursor_var.set(True)
        # plotMe2()

        #mem.Select_MT.set('(z,fission)')
        #makeMixList()

    root.mainloop()

if __name__ == '__main__':
    main()
