"""Read-only inspection for the V34 manuscript. Never invokes a solver."""
from pathlib import Path
import json
import hashlib
import zipfile
import math
import statistics
from collections import defaultdict, Counter
from xml.etree import ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
NS = '{http://schemas.openxmlformats.org/spreadsheetml/2006/main}'
WORKBOOK = Path(r'C:\Users\thoma\Documents\goc2-ac-score-check\final.xlsx')

def digest(path):
    with path.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()

def strings(z):
    out = []
    with z.open('xl/sharedStrings.xml') as f:
        for _, e in ET.iterparse(f, events=('end',)):
            if e.tag == NS+'si':
                out.append(''.join(n.text or '' for n in e.iter(NS+'t')))
                e.clear()
    return out

def value(c, shared):
    t = c.get('t')
    v = c.find(NS+'v')
    if t == 'inlineStr':
        return ''.join(n.text or '' for n in c.iter(NS+'t'))
    if v is None or v.text is None:
        return None
    if t == 's':
        return shared[int(v.text)]
    if t in ('str', 'e'):
        return v.text
    if t == 'b':
        return v.text == '1'
    try:
        return float(v.text)
    except ValueError:
        return v.text

def workbook_rows(z, path, shared):
    with z.open(path) as f:
        for _, e in ET.iterparse(f, events=('end',)):
            if e.tag == NS+'row':
                yield int(e.get('r')), {''.join(x for x in c.get('r') if x.isalpha()): value(c, shared) for c in e.findall(NS+'c')}
                e.clear()

def main():
    summary_path = ROOT/'docs/evidence/V34_ALL_SCENARIOS_SUMMARY_20260916.json'
    summary = json.loads(summary_path.read_text())
    wanted = {('C2FEN'+c['family'], int(c['scenario'])) for c in summary['cases']}
    official = defaultdict(list)
    workbook_hash = digest(WORKBOOK)
    assert workbook_hash == '56146934c9c90e5eff6477246243caa52d1483df4da24b8a5dc0dc0509ee63f0'
    with zipfile.ZipFile(WORKBOOK) as z:
        shared = strings(z)
        for number, row in workbook_rows(z, 'xl/worksheets/sheet1.xml', shared):
            key = (row.get('A'), row.get('B'))
            if key not in wanted or row.get('C') != 1:
                continue
            official[key].append(dict(excel_row=number, team=row.get('E'),
                objective=row.get('U'), workbook_score=row.get('Q'),
                mspp=row.get('I'), msgain=row.get('P'),
                infeasible=row.get('AE'), evaluator_infeasible=row.get('BL'),
                code1_seconds=row.get('Z'), code2_seconds=row.get('AA'),
                code2_limit=row.get('AB'), contingencies=row.get('AC'),
                code1_resources=row.get('X'),code2_resources=row.get('Y'),
                native_bus_count=row.get('IT'), evaluator_seconds=row.get('AF'),
                code1_timedout=row.get('AO'),code2_timedout=row.get('AP')))
    cases = []
    for case in summary['cases']:
        root = ROOT / case['run_root']
        status = json.loads((root/'run_status.json').read_text())
        run = json.loads((root/'run_summary.json').read_text())
        cert = run['official_evaluation_certificate']
        assert case['passed'] and status['success'] and status['stage']=='complete'
        assert run['git_revision']==summary['algorithm_revision']
        assert run['base_exact_executable_sha256']==run['fast_screen_executable_sha256']==summary['executable_sha256']
        assert status['completed_contingency_count']==case['contingencies']==run['contingency_count']
        assert cert['complete_label_set'] and cert['observed_unique_detail_count']==case['contingencies']+1
        assert run['official_evaluation']['infeas']==cert['reported_infeasibility']==0
        assert run['official_evaluation']['obj']==case['objective']
        assert max(run['max_independent_contingency_residual'],run['base']['base_validation']['max_residual'])==case['max_residual']<=1e-5
        assert run['total_wall_seconds']==case['seconds']<=300
        base_path = root/run['base']['artifact']
        assert digest(base_path)==run['base']['sha256']
        base = json.loads(base_path.read_text())
        family, scenario = case['family'],case['scenario']
        prefix = '617' if family=='00617' and scenario=='005' else family
        model_path = ROOT/'.data'/f'final_{family}'/f'scenario_{scenario}'/'model.json' if family in ('17700','19402') else ROOT/'.data'/f'final_{prefix}_{scenario}'/'model.json'
        model = json.loads(model_path.read_text())
        row = dict(case)
        row['summary_sha256']=digest(root/'run_summary.json')
        row['status_sha256']=digest(root/'run_status.json')
        row['input_hashes']=run['input_hashes']
        row['base_commitment_count']=run['base']['commitment_count']
        row['base_commitment_matches_source']=base['commitment']==[g['status_prev'] for _,g in sorted(model['gen'].items(),key=lambda kv:int(kv[0]))]
        assert row['base_commitment_matches_source']
        row['base_method']=run['base']['base_method']
        row['base_nonlinear_refinement']={k:base.get('sparse_ac_economic_refinement',{}).get(k) for k in ('accepted','status','error','total_wall_seconds')}
        row['model_stats']={k:len(model[k]) for k in ('bus','gen','branch','load','shunt')}
        row['model_constants']={k:model[k] for k in ('baseMVA','delta','deltar','deltactg','deltarctg','sm_vio_limit','sm_cost_approx','p_delta_cost_approx','q_delta_cost_approx')}
        row['timings']={k:run.get(k) for k in ('base_algorithm_wall_seconds','base_process_wall_seconds','contingency_parallel_wall_seconds','evaluation_wall_seconds','evaluation_tail_wall_seconds','fast_screen_parallel_wall_seconds')}
        row['routing']={k:run.get(k,status.get(k)) for k in ('fast_workers','workers','fast_screen_fallback_count','fast_screen_feasible_count','ipopt_fallback_count','linearized_contingency_accepted_count','base_pwl_epigraph','base_exact_hessian')}
        row['native_base_objective']=base['base']['objective']
        row['native_base_slacks']={k:sum(base['selected_state'][k]) for k in ('p_delta','q_delta','sm_slack')} if 'selected_state' in base else None
        row['official_objective_components']={k:v for k,v in run['official_evaluation'].items() if any(x in k for x in ('obj','cost','benefit')) and isinstance(v,(int,float))}
        best_by_team={}
        for entry in official[('C2FEN'+family,int(scenario))]:
            if 'benchmark' in entry['team'].lower() or entry['infeasible']!=0 or not isinstance(entry['objective'],(int,float)) or not math.isfinite(entry['objective']):
                continue
            assert entry['contingencies']==case['contingencies'], (family,scenario,entry)
            assert entry['evaluator_infeasible']==0, entry
            prior=best_by_team.get(entry['team'])
            if prior is None or entry['objective']>prior['objective']:
                best_by_team[entry['team']]=entry
        ordered=sorted(best_by_team.values(),key=lambda x:(-x['objective'],x['team']))
        assert len(ordered)>=5
        row['official_top5']=ordered[:5]
        row['eligible_competitors']=len(ordered)
        row['mspp']=ordered[0]['mspp']
        assert all(abs(x['mspp']-row['mspp'])<1e-7 for x in ordered)
        row['credited_gain']=max(0,case['objective']-row['mspp'])
        row['gain_fraction_of_best']=row['credited_gain']/(ordered[0]['objective']-row['mspp'])
        row['shortfall_to_fifth_percent']=100*(ordered[4]['objective']-case['objective'])/abs(ordered[4]['objective'])
        row['meets_top5']=case['objective']>=ordered[4]['objective']
        row['hypothetical_objective_position']=1+sum(x['objective']>case['objective'] for x in ordered)
        cases.append(row)
    families=[]
    for f in summary['families']:
        rows=[c for c in cases if c['buses']==f['buses']]
        families.append(dict(f, median_seconds=statistics.median(c['seconds'] for c in rows),
            min_gap5=min(c['shortfall_to_fifth_percent'] for c in rows),max_gap5=max(c['shortfall_to_fifth_percent'] for c in rows),
            median_gap5=statistics.median(c['shortfall_to_fifth_percent'] for c in rows),
            meets_top5=sum(c['meets_top5'] for c in rows),
            median_gain_capture=statistics.median(c['gain_fraction_of_best'] for c in rows),
            median_base_seconds=statistics.median(c['timings']['base_process_wall_seconds'] for c in rows),
            median_contingency_seconds=statistics.median(c['timings']['contingency_parallel_wall_seconds'] for c in rows),
            median_eval_tail_seconds=statistics.median(c['timings']['evaluation_tail_wall_seconds'] for c in rows)))
    evidence=dict(algorithm_revision=summary['algorithm_revision'],executable_sha256=summary['executable_sha256'],
        config_sha256=summary['config_sha256'],summary_sha256=digest(summary_path),
        workbook_sha256=workbook_hash,workbook_source='https://data.openei.org/files/6197/Challenge_2_Final_Event_Results.xlsx',
        selection='data worksheet; Division 1; feasible; finite raw objective U; ARPA-e benchmark excluded; best feasible row per team; descending unrounded objective; alphabetical exact-tie break',
        cases=cases,families=families)
    out=ROOT/'paper/data'
    out.mkdir(exist_ok=True)
    (out/'evidence.json').write_text(json.dumps(evidence,indent=2)+'\n',encoding='utf-8')
    print('AUDIT: all',len(cases),'saved cases match V34 acceptance and source commitments')
    print('FAMILY_SUMMARY',json.dumps(families,indent=2))
    print('TOP5_TEAMS',sorted({x['team'] for c in cases for x in c['official_top5']}))
    print('OFFICIAL_HARDWARE',Counter((x['code1_resources'],x['code2_resources']) for c in cases for x in c['official_top5']))
    for c in cases:
        print(c['family'],c['scenario'],f"ours={c['objective']:.2f}",f"fifth={c['official_top5'][4]['objective']:.2f}",f"shortfall={c['shortfall_to_fifth_percent']:.2f}%",'base='+c['base_method'])

if __name__ == '__main__':
    main()
