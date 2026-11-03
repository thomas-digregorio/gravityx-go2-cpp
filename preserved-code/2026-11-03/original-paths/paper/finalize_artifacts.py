"""Check the compiled paper and package source files. Does not run a solver."""
from pathlib import Path
import json
import hashlib
import re
import shutil
import zipfile
from pypdf import PdfReader

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
BUILD = ROOT/'tmp/pdfs'
OUT = ROOT/'output/pdf'
for path in (HERE,BUILD,OUT):
    if 'onedrive' in str(path.resolve()).lower():
        raise RuntimeError('OneDrive output prohibited')
OUT.mkdir(parents=True,exist_ok=True)
data=json.loads((HERE/'data/evidence.json').read_text())
pdf=BUILD/'main.pdf'
reader=PdfReader(pdf)
texts=[p.extract_text() or '' for p in reader.pages]
compact=re.sub(r'\s+','', '\n'.join(texts))
assert len(texts)==16, len(texts)
assert 'Conclusion' in texts[4] and 'References' in texts[5]
assert 'Complete frozen-campaign results' in texts[6]
assert len(data['cases'])==37
for c in data['cases']:
    assert f"{c['objective']:,.2f}" in compact, (c['family'],c['scenario'],'objective')
    assert f"{c['seconds']:.2f}" in compact, (c['family'],c['scenario'],'seconds')
    for t in c['official_top5']:
        assert f"{t['objective']:,.2f}" in compact, (t['team'],c['scenario'])
        assert f"{t['workbook_score']:,.0f}" in compact
assert '[?]' not in compact and '??' not in compact
log=(BUILD/'main.log').read_text()
assert 'Overfull' not in log and 'undefined' not in log
assert all(c['base_commitment_matches_source'] for c in data['cases'])
assert sum(c['contingencies'] for c in data['cases'])==44517
assert sum(c['routing']['linearized_contingency_accepted_count'] for c in data['cases'])==150
assert sum(c['routing']['ipopt_fallback_count'] for c in data['cases'])==0
assert all(c['model_constants']['baseMVA']==100 for c in data['cases'])
for c in data['cases']:
    interval = 0.0833 if c['family'] in ('00617','02020') else 0.0833333333 if c['family']=='19402' else 0.083333
    ramp_window = 0.16667 if c['family'] in ('04200','08300','17700') else interval
    assert all(c['model_constants'][p]==interval for p in ('delta','deltar','deltactg'))
    assert c['model_constants']['deltarctg']==ramp_window
assert all(c['model_constants']['sm_vio_limit']==0.2 for c in data['cases'])
portable_pdf=BUILD/'standalone-source/main.pdf'
portable_verified=False
if portable_pdf.exists():
    portable_texts=[p.extract_text() or '' for p in PdfReader(portable_pdf).pages]
    assert portable_texts==texts, 'Standalone source compile changed document text'
    portable_verified=True

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

target=OUT/'verified_ac_dispatch.pdf'
shutil.copy2(pdf,target)
files=[(HERE/'main.tex','main.tex'),(HERE/'references.bib','references.bib'),(BUILD/'main.bbl','main.bbl')]
files += [(p,'generated/'+p.name) for p in sorted((HERE/'generated').iterdir()) if p.suffix in ('.tex','.pdf')]
bundle=OUT/'verified_ac_dispatch_latex.zip'
with zipfile.ZipFile(bundle,'w',compression=zipfile.ZIP_DEFLATED) as z:
    for p,name in files:
        z.write(p,name)
with zipfile.ZipFile(bundle) as z:
    assert z.testzip() is None
    for p,name in files:
        assert z.read(name)==p.read_bytes()

report={
    'pdf_pages':len(texts),'main_pages':5,'reference_pages':1,'appendix_pages':10,
    'cases':37,'source_contingencies':44517,'official_reference_entries':185,
    'all_case_objectives_and_times_found_in_pdf':True,
    'all_official_objectives_and_scores_found_in_pdf':True,
    'undefined_references':False,'overfull_boxes':False,
    'pdf_sha256':sha(target),'latex_zip_sha256':sha(bundle),
    'evidence_sha256':sha(HERE/'data/evidence.json'),
    'standalone_source_compile_verified':portable_verified,
    'algorithm_revision':data['algorithm_revision'],
    'optimizations_run_for_manuscript':0,
    'visual_review':'All 16 rendered pages inspected; final changed pages rechecked.',
    'release_status':'Local draft only; no submission or push.'
}
(HERE/'qa_report.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
