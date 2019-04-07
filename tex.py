from glob import glob
import subprocess as sp
# from pdb import set_trace as st
'''
This module writes a LaTeX file that includes all of the figures
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
        r'''\textsuperscript{%s}%s (%s), '''%(mass, el, xs)
        indexList.append(r'''%s (%s)'''%(iso, xs))

    b = \
    r'''
    \begin{landscape}
    \begin{centering}
    \begin{figure}
      \includegraphics[width=0.7\linewidth]{%s}
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
    \hyperref[index]{index}
    \end{figure}
    \end{centering}
    \end{landscape}
    '''%(label)
    return b


#######################################################
# list of figures output by EXSAN's batch analysis
#######################################################
figs = glob('./figs/*')


#######################################################
# Top of the tex file
#######################################################
top =\
r'''
\documentclass{article}
\usepackage{graphicx}
\usepackage{pdflscape}
\usepackage[margin=0.8in]{geometry}
\usepackage{fancyhdr}
\usepackage{imakeidx}
\usepackage{hyperref}
\hypersetup{
    colorlinks=true,
    linkcolor=blue,
    filecolor=magenta,
    urlcolor=cyan
}
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage{imakeidx}
\usepackage[mmddyyyy]{datetime}
\fancyhead{}
\fancyfoot{}
\fancyhead[CO,CE]{---    Compiled by EXSAN on \today    ---}
\fancyfoot[C]{EXSAN}
\fancyfoot[R] {\thepage}
\makeindex[columns=2, title={Alphabetical Index\label{index}}]
\begin{document}
\title{EXSAN Cross Sections}
\maketitle
\pagestyle{fancy}
\thispagestyle{fancy}
'''

#######################################################
# Bottom of the tex file
#######################################################
tail = \
r'''
\printindex
\end{document}
'''

#######################################################
# Body of the tex file, one for each figure
#######################################################
body = r''''''
for fig in figs:
    body += addBody(fig)

#######################################################
# Make and write the entire tex file and PDF
#######################################################
tex = top + body + tail

with open('tex.tex', 'w') as f:
    f.write(tex)
sp.check_call(['pdflatex', 'tex.tex'])
