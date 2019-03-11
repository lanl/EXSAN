"""
Title:      EXSAN: An (E)NDF (C)ross (S)ection (An)alysis

Author:     H. Omar Wooten, PhD, DABR
            Los Alamos National Laboratory
            hasani@lanl.gov

Abstract:   This application provides an intutive  graphical interface for:
            1) Downloading complete ENDF libraries from the National Nuclear Data Center (NNDC)

            2) Automated batch processing of ENDF libraries with the Los Alamos code NJOY

            3) Plotting continuous energy and multigroup neutron reaction and photoatomic cross sections

            4) Plotting cross sections for isotopic mixes (multigroup only)

            5) Automated plotting and saving of cross section figures with an input file.

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
import numpy as np
import sys, os, operator
from glob import glob
from datadic import *
import urllib, urllib2, lxml.html
import re
from scipy.interpolate import CubicSpline as cs
import time
from pdb import set_trace as st
from pylab import connect, draw

# Global control variables
verbose = False     # invoke print statements for debugging
demo    = False      # for demo purposes (i.e. speed), do not download large tar files
auto    = True
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
                # az_menu.add_command(label=el,command=self.dummy)
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
    def __init__(self, el, frame, elementsIsotopesDict, periodicTableDict, isotopesMTdict, flag_pOrM=False ):
        self.el = el
        self.frame = frame
        self.elementsIsotopesDict = elementsIsotopesDict
        self.periodicTableDict = periodicTableDict
        self.isotopesMTdict = isotopesMTdict
        self.flag_pOrM = flag_pOrM # flag point or multigroup
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
                iso_menu.add_command(label=mtItem,command=lambda oneList=oneList:getInfo2(oneList, self.flag_pOrM, self.allFlag[0]))
                allList.append(iso+' '+mtItem)
                allListMaster[iso].append(iso+' '+mtItem)
                allListMasterMaster[self.el][iso].append(iso+' '+mtItem)
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
    def __init__(self, ax, data, fig, titles):
        self.ax = ax
        self.lx = ax.axhline(color='k')  # the horiz line
        self.ly = ax.axvline(color='k')  # the vert line
        self.data = data # a list of tuples of x,y datasets. [(xData1,yData1),(xData2,yData2)]
        self.x = []
        for d in self.data:
            self.x.append(d[0])
        self.txt = ax.text(0.65, 0.9, '', transform=ax.transAxes)
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
        self.txt.set_text('%s: \neV: %.3e\nxs: %.3e'%(self.titles[indy], x, y))
        self.fig.draw()

#===========================================================
#  This function serves as a gate for reading NJOY-processed
#  files based on whether they are a GROUPR output or not
#===========================================================
def makeMultiGroup(lines, MF, MT):
    # global groupr
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
#   Create a dictionary of MT/MF combinations from lines
#===========================================================
def make_dict():
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

    fileDirs = {}
    fileDirs[fileDirectory] = mem.tab11
    fileDirs[fileDirectory+'/multigroup'] = mem.tab12
    mem.allListMaster = {}
    mem.allListMasterMaster = {}

    for mem.fileDir, tab in fileDirs.iteritems():
        mem.files = glob('./'+mem.fileDir+'/*.txt')
        isotopesTmp = [f.split('/')[-1].split('.txt')[0] for f in mem.files]
        mem.elements = []
        mem.isotopes = []

        if len(mem.files) > 0:
            mem.logText.insert(INSERT, str(len(mem.files))+' isotopes in '+mem.fileDir+' directory. \n')
        else:
            mem.logText.insert(INSERT, mem.fileDir+' directory is empty. \n')

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
            if verbose:
                print 'Reading %s'%file
            popup.update()
            a = file.split('/')[-1].split('.')[0]
            b = re.split('(\d+)',a)
            if len(b) > 1:
                mem.elements.append(b[0])
                if len(b) > 2:
                    c = b[0]+'-'+b[1]+b[2]
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
                reactions = make_dict()
                mem.isotopesMTdict[c] = reactions
            except:
                continue

            if not mem.fileDir == 'endfvii-atomic':
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
                  mem.elementsIsotopesDict, periodicTableDict, mem.isotopesMTdict, flag_pOrM)
            except:
                continue
        popup.destroy()

        az = [chr(i) for i in range(65,91)]
        mem.azDict = {}
        for letter in az: mem.azDict[letter] = []
        for el in mem.elements: mem.azDict[el[0]].append(el)

#===========================================================
#   Get everything ready to plot
#===========================================================
def getInfo2(isotopes, flag_pOrM, allFlag):
    if allFlag and len(mem.plotTitles)>0:
        mem.plotX = {}
        mem.plotY= {}
        for title in mem.plotTitles:
            mem.logText.insert(INSERT, title+' removed from plot\n')
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
            f = open(file)
            mem.lines = f.readlines()
            mem.Title_var.set(iso)
            if mem.particle.get() == 1:  # For X-ray files, MF==23
                mem.MF_var.set('23')
            else:
                mem.MF_var.set('3')
            mem.MT_var.set(str(mtID))
            mem.MT_descript.set(mtDescript)
            addPlot(flag_pOrM, allFlag)

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

    if mem.microMacro.get() == 1 and len(periodicTableDict[mem.element])== 5:
        N = periodicTableDict[mem.element][4]/periodicTableDict[mem.element][3]*Av
        y = N*y*1.e-24
        mem.logText.insert(INSERT,'Number density for '+mem.element+' = '+str(N)+'\n')

    elif mem.microMacro.get() == 2 and len(periodicTableDict[mem.element])== 5:
        N = periodicTableDict[mem.element][4]/periodicTableDict[mem.element][3]*Av
        y = 1/(N*y*1.e-24)
        mem.logText.insert(INSERT,'Number density for '+mem.element+' = '+str(N)+'\n')

    elif mem.microMacro.get() == 3 and len(periodicTableDict[mem.element])== 5:
        N = periodicTableDict[mem.element][4]/periodicTableDict[mem.element][3]*Av
        y = N*y*1.e-24/periodicTableDict[mem.element][4]
        mem.logText.insert(INSERT,'Number density for '+mem.element+' = '+str(N)+'\n')

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
    mem.logText.insert(INSERT, title+' added to plot\n')
    mem.logText.insert(INSERT, 'Atomic mass  : '+str(periodicTableDict[mem.element][3])+'\n')
    mem.logText.insert(INSERT, 'Density(g/cc): '+str(periodicTableDict[mem.element][4])+'\n')
    mem.logText.insert(INSERT, '--------------------------\n\n')

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
    mem.logText.insert(INSERT, title+' removed from plot\n')
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
def plotMe2(event=None, allFlag=False, saveFlag=False):
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
        # return intTop/intBtm
        # return 1.0

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
        figImg = fig.set_size_inches(15., 15.)
        ax1 = pl.subplot2grid((2, 4), (0, 0), colspan=2)
    else:
        fig = pl.figure(1)
        figImg = fig.set_size_inches(15., 7.5)
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

    mem.columnText = ['','Th\n(0.025 eV)','Fiss\n(1 MeV)','Fus\n(14 MeV)',
                  'Custom\n('+'%3.2e'%(mem.Select_E1.get())+' eV)','AvgSigma','ResIntegral',
                  'FissSpecAvg']


    lineWidth = 2
    legendFontSize = 14 # useful for journal article plots
    if allFlag:
        legendFontSize = 10

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

        if 'mg' in title or '$' in title:
            ax1.loglog(x,y, linestyle=l, c=mem.c[colorIdx],lw=lineWidth, drawstyle='steps', label=str(i+1)+'. '+title)
            ax1.set_xscale('log')
            ax1.set_yscale('log')
            mixFlag = False

        else:
            ax1.loglog(x,y, linestyle=l, c=mem.c[colorIdx],lw=lineWidth, label=str(i+1)+'. '+title)
        mem.cell_text.append(['%-2s'%(i+1),'%5.3e'%(therm),'%5.3e'%(fiss),'%5.3e'%(fus),'%5.3e'%(user),
                          '%5.3e'%(avgSigma),'%5.3e'%(resInt), '%5.3e'%(fissSpecAvg)])
        mem.rowText.append(title)

        # cursor snap
        snapData.append((x,y))

    unitsDict = {0:"Cross Section (barns)",
                 1:"Cross Section (1/cm)",
                 2:"Mean free path (cm)",
                 3:"Mean free path (cm^2/g)"}


    tableFormat = [1.1, 0.0, 1.2, 1.0] if allFlag else [1.1, 0.0, 1.35, 0.75]
    the_table = ax1.table(cellText=mem.cell_text,
                        loc='right',
                        colWidths=[0.15,0.6,0.5,0.5,0.6,0.5,0.5,0.6],
                        colLabels=mem.columnText,
                        colLoc='center',
                        rowLoc='center',
                        cellLoc='center',
                        bbox=tableFormat)

    the_table.auto_set_column_width([-1,-1,-1,-1,-1,-1,-1,-1])
    the_table.auto_set_font_size(True)
    the_table.scale(4,2)
    ax1.set_xlim([0.6*np.array(mem.minimumX).min(),1.4*np.array(mem.maximumX).max()])
    ax1.set_ylim([0.6*np.array(mem.minimumY).min(),1.4*np.array(mem.maximumY).max()])
    ax1.grid(which='both', linestyle='--',color='0.8')


    ax1.set_ylabel(unitsDict[mem.microMacro.get()], fontsize=legendFontSize)
    ax1.set_xlabel("Energy (eV)", fontsize=legendFontSize)
    ax1.set_title(unitsDict[mem.microMacro.get()], fontsize=legendFontSize+2)
    ax1.legend(loc=3, prop={'size':legendFontSize})
    ax1.tick_params(labelsize=legendFontSize-2)

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
                box = ax.get_position()
                ax.legend(prop={'size':legendFontSize}, bbox_to_anchor=[0.75+i*0.1, 0.15], labels=['%1.1f %%, %s' % (s,l) for s,l in sl])
                ax.set_title(titleList[i])

    isotopes = [i.split()[0] for i in mem.plotTitles]

    if len(np.unique(isotopes))>1:
        pl.suptitle('Various Isotopes')
    else:
        pl.suptitle(mem.plotTitles[0].split()[0])
    mem.logText.insert(INSERT, 'Processing %s\n'%(mem.plotTitles[0].split()[0]))

    # snap
    if mem.Cursor_var.get() and not allFlag:
        cursor = SnapToCursor(ax1, snapData, fig.canvas, plotOrder)
        pl.connect('motion_notify_event', cursor.mouse_move)


    if saveFlag==True:
        if not os.path.exists('figs'):
            os.mkdir('figs')
        pdf = mem.savePDF('figs/'+mem.plotTitles[0][:7]+'.pdf')
        pdf.savefig(figImg)
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
        iso = newName.split('_')[-1][:3].lstrip('0')
        if 'neutrons' in name:
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

                    if verbose: print file, isotope

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
def addMe(*args):
    import subprocess as sp

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


    for k in sorted(nndcDict.keys()):
        if mem.downloadList[k].get() == 1:
            if urllib2.urlopen(nndcDict[k][0]).geturl() == nndcDict[k][0]:

                kDir = k.replace('/','_').replace('-','_')
                dest = nndcDict[k][0].split('/')[-1].replace('-','_')
                if not os.path.isdir(kDir):
                    sp.check_call(['mkdir', kDir])

                if not demo:
                    # download tar file from NNDC
                    print 'downloading %s data file from\n %s\n'%(k,nndcDict[k][0])
                    start = time.time()
                    # mem.urlFile.retrieve(nndcDict[k][0], '/'.join([kDir, dest]))
                    download(nndcDict[k][0], '/'.join([kDir, dest]))
                    end = time.time()
                    print '%s downloaded in %7.3f sec'%(dest,(end-start))
                # untar the file in its respective fileDirectory
                sp.check_call(['tar','-xvf', '/'.join([kDir, dest]), '-C',kDir+'/'])
                print k,'has been dowloaded. Post-processing in progress...'
                nndcParse(k)
                print k,'has been successfully post-processed.'
            else:
                print k,'cannot be downloaed from NNDC at this time.'



def fileChecker(root):
    numFilesFound=0
    for k,v in dirDict2.iteritems():
        for i, value in enumerate(v):
            checkDir = './'+value+'/*.txt'
            if os.path.isdir('./'+value):
                files = glob(checkDir)
                if len(files) > 0:
                    mem.logText.insert(INSERT, str(len(files))+' files detected in '+value+' folder  \n')
                    numFilesFound += len(files)
                    mem.verSelect.set(i)

    if numFilesFound ==0:
        mem.logText.insert(INSERT, 'Must download ENDF files first.\n')
        root.bind('<Escape>', close)
        root.mainloop()

#===========================================================
#   Attempt to mass process all files in directory and return
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



    print 'Searching for isotopes with '+mem.Select_MT.get()+' reaction...'
#     logText.insert(INSERT,'Searching for isotopes with '+Select_MT.get()+' reaction...\n')

    mem.preCheckedFiles=[]
    if mem.Select_MT.get() == 'absorption':
        absFlag = True
    if mem.note0.index(mem.note0.select()) == 1:
        whichTab = mem.note1.index(1)
    else:
        whichTab = mem.note1.index("current")

    popup = Toplevel()
    popup_lab2 = Label(popup, text='Analyzing isotopes with '+mem.Select_MT.get()+' reaction', width=50)
    status2 = 0
    status_var2 = DoubleVar()
    progress_bar2 = ttk.Progressbar(popup, variable=status_var2, maximum=100)
    try:
        statusInc1 = 100.0/len(mem.fileDict_all[whichTab].keys())
        popup_lab2.grid(row=2, column=0)
        progress_bar2.grid(row=3, column=0, sticky=EW)#.pack(fill=tk.X, expand=1, side=tk.BOTTOM)
    except:
        pass


    mem.preCheckedFiles = []
    for i in mem.reactionsDict[mem.periodicTableTabs[whichTab][0]][mem.Select_MT.get()]:
        mem.preCheckedFiles.append(i[1])

    mem.logText.insert(INSERT,'%s isotopes found\n'%(str(len(mem.preCheckedFiles))))

    status = 0
    tally = 0

    try:
        statusInc2 = 100.0/len(mem.preCheckedFiles)
    except:
        popup_lab2 = Label(popup, text='Select ENDF library (with multigroup data) first', width=50)
        popup_lab2.grid(row=0, column=0)
        # popup.update()
        return

    for i, file in enumerate(mem.preCheckedFiles):
        popup_lab2 = Label(popup, text='Analyzing file '+str(i+1)+'/'+str(len(mem.preCheckedFiles))+' with '+mem.Select_MT.get()+' reaction.', width=50)
        popup_lab2.grid(row=2, column=0)
        popup.update()
        tally+=1
        status+=1

        f=open(file)
        isotope = mem.fileDict_all[whichTab][file]
        mem.element = isotope.split('-')[0]
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

    popup.destroy()
    mem.logText.insert(INSERT, 'Analysis complete.\n')
    displayAnalysis()


def displayAnalysis():

    def dictToArray(dict,num):
        a = sorted(dict.items(), key=operator.itemgetter(num))
        a = np.array(a)
        return a


    xs_list = [mem.thermalDict, mem.fissionDict, mem.fusionDict, mem.eUserDict]
    xs_listS = ['thermalDictSorted', 'fissionDictSorted', 'fusionDictSorted', 'eUserDictSorted']
    e_list  = ['thermal', 'fiss', 'fus', 'user']

    mem.resultsToPrint = dictToArray(xs_list[mem.eRange.get()], mem.sortBy.get())

    mem.resultsText.delete('1.0', END)
    if mem.microMacro.get() == 1:
        mem.resultsText.insert(INSERT, u"Isotope         \u03a3(1/cm)\n")
    elif mem.microMacro.get() == 2:
        mem.resultsText.insert(INSERT, u"Isotope         Mean_Free_Path(cm)\n")
    elif mem.microMacro.get() == 3:
        mem.resultsText.insert(INSERT, u"Isotope         \u03a3/\u03c1(cm^2/g)\n")
    else:
        mem.resultsText.insert(INSERT, u"Isotope         \u03c3(barns)\n")

    for i in range(len(mem.resultsToPrint)):
        if isinstance(mem.resultsToPrint[i,1], basestring):
            mem.resultsText.insert(INSERT, '%3i %-10s %8.4e\n' %(i+1, mem.resultsToPrint[i,0],float(mem.resultsToPrint[i,1]))) #'%2s%3s' % (MF, MT)
        else:
            mem.resultsText.insert(INSERT, '%3i %-10s %8.4e\n' %(i+1, mem.resultsToPrint[i,0],float(mem.resultsToPrint[i,1]))) #'%2s%3s' % (MF, MT)

#===========================================================
#   Process all raw ENDF files with NJOY
#===========================================================
def run_njoy(*args):

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
        os.system('./njoy < njoy.in')
        os.system('mv tape22 '+file.split('/')[-1])
        os.system('mv tape23 multigroup/'+file.split('/')[-1])
        os.chdir('../')
        os.system('pwd')

    processList = []
    for key,value in mem.downloadList.iteritems():
        if value.get()==1 and key!='ver_xRay':
            processList.append(key)
    if not os.path.isdir('./working/multigroup'):
        os.mkdir('./working/multigroup')
    for value in processList:
        newValue = dirDict2[value][mem.particle.get()]
        checkDir = './'+newValue+'/*.txt'
        myDir = './'+newValue+'/'
        if os.path.isdir('./'+newValue):
            if not os.path.isdir('./'+newValue+'/multigroup'):
                os.mkdir('./'+newValue+'/multigroup')
            files = glob(checkDir)
            for file in files:
                if verbose: print 'NJOY working on %s'%(file)
                NJOY(file)

            os.system('mv ./working/*.txt '+myDir)
            os.system('mv ./working/multigroup/*.txt '+myDir+'/multigroup')
        else:
            pass

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
    mem.logText.insert(INSERT, 'Analysis results written to file '+outFileName+'\n')

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
    mem.logText.insert(INSERT, 'Plots cleared \n')
    mem.batchFile = None

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
    # list of fileDicts corresponding to tab1 (point), tab2 (multigroup)
    mem.fileDict_all = [mem.fileDict_p, mem.fileDict_m]
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

    # keep track of these tabs for later
    mem.periodicTableTabs = [(mem.tab11,'Pointwise'), (mem.tab12,'Multigroup')]

    Grid.rowconfigure(mem.tab11, 0, weight=1)
    Grid.columnconfigure(mem.tab12, 0, weight=1)
    mem.tab11.columnconfigure(0, weight=1)
    mem.tab12.columnconfigure(0, weight=1)
    mem.tab11.grid(row=0, column=0, sticky=NSEW)
    mem.tab12.grid(row=0, column=0, sticky=NSEW)
    mem.isotopeGrid11 = Frame(mem.tab11)
    mem.isotopeGrid12 = Frame(mem.tab12)
    mem.isotopeGrid11.grid(column=0, row=7, sticky=NSEW)
    mem.isotopeGrid11.grid(column=0, row=7, sticky=NSEW)
    mem.note1.add(mem.tab11,text='Pointwise')
    mem.note1.add(mem.tab12,text='Multigroup')
    mem.note1.pack(side='top', fill='both', expand=True)
    s = ttk.Style()
    s.configure('TNotebook.Tab', padding=(20, 8, 20, 0))

    mem.user_frame =  Frame(mem.tab01,width=mem.rootWidth, height=50, pady=1,padx=5)
    mem.user_frame.columnconfigure(0, weight=1)
    mem.userE_frame =  Frame(mem.user_frame,width=mem.rootWidth/2, height=50, pady=1)
    mem.userMT_frame =  Frame(mem.user_frame,width=mem.rootWidth/2, height=50, pady=1)

    mem.sortBy_frame = Frame(mem.tab01, width=mem.rootWidth, height=50, pady=7, relief=GROOVE, borderwidth=2)
    mem.sortBy_frame.columnconfigure(0, weight=1)
    mem.sortByXS_frame = Frame(mem.sortBy_frame, width=mem.rootWidth/2, height=50, padx=20)
    mem.sortByE_frame = Frame(mem.sortBy_frame, width=mem.rootWidth/2, height=50, padx=20)
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
    mem.sortBy_frame.grid(row=3)
    mem.sortByXS_frame.grid(row=0, column=0)
    mem.sortByE_frame.grid(row=0, column=1)
    mem.dataText_frame.grid(row=4, column=0, sticky=EW)
    mem.dataTextData_frame.grid(row=0,column=0, sticky=EW)
    mem.dataTextLog_frame.grid(row=0,column=1, sticky=EW)

    mem.downloadList={}
    for i, ver in enumerate(np.sort(nndcDict.keys())):
            mem.downloadList[ver] = IntVar()
            cb = Checkbutton(mem.download_frame, variable=mem.downloadList[ver], text='%-15s'%(ver), font=mem.HDG1)
            cb.grid(row=0, column=i)

    mem.selectAll_btn = Button(mem.download_frame,text='Select All',command=selectAll,font=mem.HDG1)
    mem.selectAll_btn.grid(row=1, column=1)
    mem.reset = Button(mem.download_frame,text='Deselect All',command=reset,font=mem.HDG1)
    mem.reset.grid(row=1, column=2)
    mem.get = Button(mem.download_frame,text='Download Now',command=addMe,font=mem.HDG1)
    mem.get.grid(row=1, column=3)
    mem.njoy = Button(mem.download_frame,text='Process with NJOY',command=run_njoy,font=mem.HDG1)
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
            value=i, command=update_files, font=mem.HDG1)
        mem.verSelectRad[ver].grid(row=1, column=i, padx=1, pady=1)
    mem.particle = IntVar()
    mem.particleSelect = {}
    particles = ['neutrons','x-rays']
    for j, ver in enumerate(particles):
        mem.particleSelect[ver] = Radiobutton(mem.version_frame,variable=mem.particle,text=ver,
           value=j, font=mem.HDG1)
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
    mem.RunBatch_Btn = Button(mem.userMT_frame, width=10, text='Batch', command=runBatch, padx=2, font=mem.HDG1)
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
    mem.sortByName = Radiobutton(mem.sortByXS_frame, variable=mem.sortBy, value=0,text='Sort by Isotope',command=displayAnalysis,font=mem.HDG1)
    mem.sortByXS = Radiobutton(mem.sortByXS_frame, variable=mem.sortBy, value=1,text='Sort by x-sec',command=displayAnalysis,font=mem.HDG1)
    mem.sortByName.grid(row=0,column=0)
    mem.sortByXS.grid(row=0,column=1)

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

    mem.mix = IntVar()
    mem.mixFrac = DoubleVar()
    mem.mixFrac.set(1.0)

    mem.scrollbar1 = Scrollbar(mem.dataTextData_frame)
    mem.resultsText = Text(mem.dataTextData_frame, height=20, bg='aliceblue',font=mem.DATA,yscrollcommand=mem.scrollbar1.set)
    mem.resultsText.grid(row=0, column=0,sticky=NSEW,padx=1)
    mem.scrollbar1.grid(row=0, column=1, sticky=NS)

    mem.scrollbar2 = Scrollbar(mem.dataTextLog_frame)
    mem.scrollbar2.grid(row=0, column=1, sticky=NS)
    mem.logText = Text(mem.dataTextLog_frame, height=20, bg='aliceblue',font=mem.DATA,yscrollcommand=mem.scrollbar2.set)
    mem.logText.grid(row=0, column=0,sticky=NSEW,padx=1)


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
    if len(sys.argv) > 1:
        mem.batchFile = sys.argv[1]
        runBatch()
    else:
        mem.batchFile = None

def runBatch():
    if mem.batchFile == None:
        batchFile = tkFileDialog.askopenfilename()
    else:
        batchFile = mem.batchFile
    # with open(tkFileDialog.askopenfilename()) as f:
    with open(batchFile,'r') as f:
        lines = f.readlines()
    f.close()
    commentLines = [i for (i,l) in enumerate(lines) if l[0]=='#']
    for cl in commentLines[::-1]:
        lines.pop(cl)
    # exec(lines[0])
    i=0
    tmp = {}
    for k,v in sorted(dirDict2.iteritems()):
        tmp[k]=i; i+=1

    if lines[0].strip() in sorted(dirDict2.keys()):
        if not mem.verSelect.get()==tmp[lines[0].strip()]:
            mem.verSelect.set(tmp[lines[0].strip()])
            update_files()

    lines.pop(0)
    mem.logText.delete('1.0', END)
    mem.logText.insert(INSERT, 'Batch processing: %s\n'%(batchFile))

    saveFlag_master = flag_pOrM_master = allFlag_master = False
    saveFlag = flag_pOrM = allFlag = False
    if 'E' in lines[0]:
        mem.Select_E1.set(lines[0].split('=')[1].split()[0])
        popLast = True
    if 'save' or 'all' in lines[0]:
        from matplotlib.backends.backend_pdf import PdfPages as savePDF
        mem.savePDF = savePDF
        saveFlag_master = True
        print 'save master'
        popLast = True
    if 'm' in lines[0]:
        flag_pOrM_master = True
        mem.note1.select(1)
        print 'pOrM master'
        popLast = True

    if 'all' in lines[0]:
        allFlag = True
        for key in sorted(mem.allListMasterMaster):
            for key1 in sorted(mem.allListMasterMaster[key]):
                print key1
                getInfo2(mem.allListMasterMaster[key][key1], flag_pOrM_master, allFlag_master)
                plotMe2(None, allFlag_master, saveFlag_master)
                clearAll()
    elif 'rank' in lines[0]:
        print 'rank mode'
        lines.pop(0)
        for line in lines:
            if line.split()[0] in list(mem.tmpList):
                mem.Select_MT.set(line.split()[0])
                master_list(mem.note1.index("current"))
                saveText()

    else:
        if popLast:
            lines.pop(0)
        for line in lines:
            if line.split()[0] in mem.isotopesMTdict.keys() and line.split()[1] in mem.isotopesMTdict[line.split()[0]]:
                if not flag_pOrM_master:
                    if 'm' in line.split():
                        flag_pOrM = True
                        mem.note1.select(1)
                    else:
                        flag_pOrM = False
                        mem.note1.select(0)
                else:
                    flag_pOrM = True

                if not saveFlag_master:
                    if 'save' in line.split():
                        saveFlag = True
                    else:
                        saveFlag = False
                else:
                    saveFlag = True

                if 'all' in line.split():
                    allFlag = True
                    line = mem.allListMaster[line.split()[0]]
                else:
                    allFlag = False
                    line = [' '.join(line.split()[:2])]
                getInfo2(line, flag_pOrM, allFlag)
                # plotMe2(None, allFlag, saveFlag)
                # clearAll()
            else:
                print line.split()[0], 'not found'
                continue
        plotMe2(None, allFlag, saveFlag)
        # clearAll()


def main():
    print "                    _______  ______    _    _   _ "
    print "                   | ____\ \/ / ___|  / \  | \ | |"
    print "                   |  _|  \  /\___ \ / _ \ |  \| |"
    print "                   | |___ /  \ ___) / ___ \| |\  |"
    print "                   |_____/_/\_\____/_/   \_\_| \_|\n"
    print "-- The (E)NDF (C)ross (S)ection (An)alysis and visualization appication --\n"
    print "                             Created by:"
    print "                    H. Omar Wooten, PhD, DABR"
    print "                          hasani@lanl.gov"
    print "            Copyright Los Alamos National Laboratory, 2018  "
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

    makeWidgets(root)
    fileChecker(root)

    # Keyboard shortcuts.
    root.bind('<Control-Escape>', close)
    root.bind('a',addPlot)
    root.bind('d',delPlot)
    root.bind('p',plotMe2)
    root.bind('c',clearAll)

    # for testing, automatically load ENDF 8
    if auto:
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