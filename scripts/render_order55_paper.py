"""Build the canonical joint LaTeX manuscript and retain compatible PDF paths."""
from pathlib import Path
import subprocess,shutil
ROOT=Path(__file__).resolve().parents[1];paper=ROOT/"paper/order55_global";source=paper/"arxiv/main.tex"
work=ROOT/"tmp/pdfs/joint";work.mkdir(parents=True,exist_ok=True)
out=ROOT/"output/pdf";out.mkdir(parents=True,exist_ok=True)
for stage in ("pdflatex","bibtex","pdflatex","pdflatex"):
 if stage=="bibtex":
  # Copy the unchanged BibTeX database beside the build output.
  shutil.copyfile(paper/"arxiv/references.bib",work/"references.bib")
  cmd=["bibtex","main"]
 else:cmd=["pdflatex","-interaction=nonstopmode","-halt-on-error",f"-output-directory={work}",str(source)]
 run=subprocess.run(cmd,cwd=work,capture_output=True)
 if run.returncode:print(run.stdout.decode(errors="replace"));raise SystemExit(run.returncode)
pdf=work/"main.pdf"
dist=ROOT/"dist";dist.mkdir(parents=True,exist_ok=True)
for path in (dist/"order55-circulant-maxdet-v1.0.0.pdf",out/"order55_joint.pdf",out/"order55_global.pdf",paper/"arxiv/main.pdf"):shutil.copyfile(pdf,path)
shutil.copyfile(source,paper/"manuscript.tex")
print(out/"order55_joint.pdf")
