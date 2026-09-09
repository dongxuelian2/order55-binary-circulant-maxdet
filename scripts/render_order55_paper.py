"""Typeset the complete Markdown research note with local LaTeX."""
from pathlib import Path
import re,subprocess
ROOT=Path(__file__).resolve().parents[1]
work=ROOT/'tmp/pdfs';work.mkdir(parents=True,exist_ok=True)
out=ROOT/'output/pdf';out.mkdir(parents=True,exist_ok=True)
symbols={0x03bb:r'\lambda',0x03a3:r'\sum',0x03a6:r'\Phi',0x03b1:r'\alpha',0x03b2:r'\beta',0x03bc:r'\mu',0x03c0:r'\pi',0x03b6:r'\zeta',0x03c6:r'\varphi',0x00d7:r'\times',0x2212:'-',0x21a6:r'\mapsto',0x2265:r'\ge',0x2264:r'\le',0x2208:r'\in',0x00b7:r'\cdot',0x2261:r'\equiv',0x00b1:r'\pm'}
escapes={'_':r'\_','%':r'\%','&':r'\&','#':r'\#','$':r'\$','{':r'\{','}':r'\}','^':r'\textasciicircum{}','~':r'\textasciitilde{}','\\':r'\textbackslash{}'}
def mathtext(s):
    s=s.replace('D_circ',r'D_{\rm circ}')
    return ''.join((symbols[ord(c)]+' ') if ord(c) in symbols else ('^2' if c=='\u00b2' else c) for c in s)
def plain(s):
    tokens=re.split(r'(\S*[_^]\S*)',s);out=[]
    for t in tokens:
        if '_' in t or '^' in t:
            out.append(r'\ensuremath{'+mathtext(t)+'}')
        else:
            for c in t:
                if ord(c) in symbols:out.append(r'\ensuremath{'+symbols[ord(c)]+'}')
                elif c=='\u00b2':out.append(r'\textsuperscript{2}')
                elif c in ('\u2013','\u2014'):out.append('--')
                elif c=='\u0107':out.append(r"\'{c}")
                else:out.append(escapes.get(c,c))
    return ''.join(out)
def inline(s):
    chunks=re.split(r'(`[^`]+`|\*\*[^*]+\*\*|\*[^*]+\*|<https?://[^>]+>)',s);out=[]
    for t in chunks:
        if t.startswith('`'):out.append(r'\path|'+t[1:-1]+'|')
        elif t.startswith('**'):out.append(r'\textbf{'+plain(t[2:-2])+'}')
        elif t.startswith('*'):out.append(r'\emph{'+plain(t[1:-1])+'}')
        elif t.startswith('<http'):out.append(r'\url{'+t[1:-1]+'}')
        else:out.append(plain(t))
    return ''.join(out)
lines=(ROOT/'paper/order55_global/manuscript.md').read_text(encoding='utf-8-sig').splitlines();body=[];i=4
while i<len(lines):
    line=lines[i]
    if line.startswith('```'):
        block=[];i+=1
        while i<len(lines) and not lines[i].startswith('```'):block.append(lines[i]);i+=1
        body.extend([r'\begin{lstlisting}',*block,r'\end{lstlisting}'])
    elif line.strip()==r'\[':
        body.append(r'\[');i+=1
        while i<len(lines) and lines[i].strip()!=r'\]':body.append(lines[i]);i+=1
        body.append(r'\]')
    elif line.startswith('|'):
        rows=[]
        while i<len(lines) and lines[i].startswith('|'):
            cells=[x.strip() for x in lines[i].strip('|').split('|')]
            if not all(re.fullmatch(r'[-:]+',c) for c in cells):rows.append(cells)
            i+=1
        n=len(rows[0]);body.extend([r'\begin{center}',r'\footnotesize',r'\setlength{\tabcolsep}{4pt}',r'\begin{tabular}{'+'l'+'r'*(n-1)+'}',r'\toprule'])
        for j,row in enumerate(rows):
            if j==0 and n==5:row=['k','Allowed profiles','Enumerated','Canonical','Lift targets']
            body.append(' & '.join(inline(c) for c in row)+r' \\')
            if j==0:body.append(r'\midrule')
        body.extend([r'\bottomrule',r'\end{tabular}',r'\end{center}']);continue
    elif line.startswith('### '):body.append(r'\subsection*{'+inline(line[4:])+'}')
    elif line.startswith('## '):body.append(r'\section*{'+inline(line[3:])+'}')
    elif line.startswith('> '):body.append(inline(line[2:]))
    elif re.match(r'^\d+\. ',line):
        para=[line];i+=1
        while i<len(lines) and lines[i] and not re.match(r'^\d+\. ',lines[i]) and not lines[i].startswith('#'):
            para.append(lines[i].strip());i+=1
        body.append(r'\par\medskip\noindent '+inline(' '.join(para)));continue
    elif line and not line.startswith('#'):
        para=[line];i+=1
        while i<len(lines) and lines[i] and not lines[i].startswith(('#','```','|','>')) and lines[i].strip()!=r'\[':
            para.append(lines[i]);i+=1
        body.append(inline(' '.join(para)));continue
    else:body.append(inline(line))
    i+=1
preamble=r'''\documentclass[11pt,a4paper]{article}
\usepackage[T1]{fontenc}
\usepackage[utf8]{inputenc}
\usepackage{lmodern,amsmath,amssymb,booktabs,listings,hyperref,microtype}
\usepackage[margin=26mm]{geometry}
\hypersetup{colorlinks=true,urlcolor=blue,linkcolor=blue,pdftitle={The Maximum Determinant of a 55 by 55 Binary Circulant Matrix}}
\lstset{basicstyle=\ttfamily\footnotesize,breaklines=true,columns=fullflexible,keepspaces=true}
\setlength{\parindent}{0pt}
\setlength{\parskip}{5pt}
\setlength{\emergencystretch}{2em}
\title{The Maximum Determinant of a\\55 $\times$ 55 Binary Circulant Matrix}
\author{Reproducible computational research note}
\date{9 September 2026}
\begin{document}
\maketitle
'''
tex=ROOT/'paper/order55_global/manuscript.tex';tex.write_text(preamble+'\n'.join(body)+'\n\\end{document}\n',encoding='utf8')
for _ in range(2):
    result=subprocess.run(['pdflatex','-interaction=nonstopmode','-halt-on-error',f'-output-directory={work}',str(tex)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    if result.returncode:print(result.stdout.decode(errors='replace'));raise SystemExit(result.returncode)
(out/'order55_global.pdf').write_bytes((work/'manuscript.pdf').read_bytes());print(out/'order55_global.pdf')
