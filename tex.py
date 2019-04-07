from glob import glob
import subprocess as sp
from pdb import set_trace as st


'''
Helpful links:
    https://es.overleaf.com/learn/latex/hyperlinks
    https://www.overleaf.com/learn/latex/Indices
'''

def addBody(s):
    isotope = s.split('/')[-1].split('_')[0]
    el, mass = isotope.split('-')
    b = \
    r'''
    \fancyfoot[R]{\thepage}
    \begin{landscape}
    \begin{centering}
    \begin{figure}
      \includegraphics[width=0.9\linewidth]{%s}
      \caption{\textsuperscript{%s}%s.}
      \index{%s.}
      \label{fig:%s}
      \hyperref[index]{index}
    \end{figure}
    \end{centering}
    \end{landscape}
    '''%(s, mass, el, isotope, isotope)
    return b



figs = glob('./figs/*')

'''
\documentclass{article}
\usepackage{graphicx}
\usepackage{pdflscape}
\usepackage[margin=0.8in]{geometry}
\usepackage{fancyhdr}
\usepackage{imakeidx}
\makeindex[columns=3, title=Alphabetical Index]
% Turn on the style
\fancyhead{}
\fancyfoot{}
% Set the right side of the footer to be the page number
\fancyfoot[R]{\thepage}
\pagestyle{fancy}

\begin{document}
\title{EXSAN Cross Sections}
\author{H. Omar Wooten, Ph.D.}
\maketitle

\begin{landscape}
\begin{figure}
  \includegraphics[width=\linewidth]{./figs/Li-6.pdf}
  \caption{\index{\textsuperscript{6}Li.}}
  \label{fig:Li-6}
\end{figure}
\end{landscape}


\printindex
\end{document}
'''

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
\makeindex[columns=4, title={Alphabetical Index\label{index}}]
\fancyhead{}
\fancyfoot{}
% Set the right side of the footer to be the page number
\fancyfoot[R]{\thepage}
% Turn on the style
\pagestyle{fancy}

\begin{document}
\title{EXSAN Cross Sections}
\author{H. Omar Wooten, Ph.D.}
\maketitle
\thispagestyle{fancy}
'''

tail = \
r'''
\printindex
\end{document}
'''

body = r''''''
for fig in figs:
    body += addBody(fig)



tex = top + body + tail
# tex = head + main + foot

with open('tex.tex', 'w') as f:
    f.write(tex)
sp.check_call(['pdflatex', 'tex.tex'])
