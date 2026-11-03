"""Generate manuscript tables and one vector figure from audited saved results.

No optimization is launched. Requires matplotlib only for the figure.
"""
from pathlib import Path
import json
import statistics
import hashlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
DATA = json.loads((HERE/'data/evidence.json').read_text())
OUT = HERE/'generated'
OUT.mkdir(exist_ok=True)
TEAMS = {
    'Electric Stampede':'ES', 'GMI-GO':'GMI', 'GO-SNIP':'GS',
    'GOT-BSI-OPF':'GOT', 'Gordian Knot':'GK', 'GravityX':'GX',
    'M&H Grid Solutions':'MH', 'Monday Mornings':'MM',
    'NU_Columbia_Artelys':'NCA', 'Pearl Street Technologies':'PS', 'magos':'MA',
}
ROW_END = chr(92)*2

def esc(s):
    return str(s).replace('&',r'\&').replace('_',r'\_').replace('%',r'\%')

def write(name, lines):
    (OUT/name).write_text('\n'.join(lines)+'\n',encoding='utf-8')

def fmt(x, digits=2):
    return f'{x:,.{digits}f}'

rows=[]
timings=[]
allcase=[]
for f in DATA['families']:
    cases=[c for c in DATA['cases'] if c['buses']==f['buses']]
    n=len(cases)
    rows.append(f"{f['buses']:,} & {n}/{n} & {sum(c['contingencies'] for c in cases):,} & "
        f"{min(c['seconds'] for c in cases):.2f} & {statistics.median(c['seconds'] for c in cases):.2f} & "
        f"{max(c['seconds'] for c in cases):.2f} & {f['median_gap5']:.2f} " + ROW_END)
    timings.append(f"{f['buses']:,} & {f['median_base_seconds']:.2f} & {f['median_contingency_seconds']:.2f} & "
        f"{f['median_eval_tail_seconds']:.2f} & {sum(c['routing']['fast_screen_fallback_count'] for c in cases)} & "
        f"{100*f['median_gain_capture']:.2f} " + ROW_END)
    for c in cases:
        r=f"{c['max_residual']:.2e}".split('e')
        residual=rf"${r[0]}\!\times\!10^{{{int(r[1])}}}$"
        allcase.append(f"{c['buses']:,} & {c['scenario']} & {c['normalized_model_buses']:,} & {c['contingencies']:,} & "
            f"{c['seconds']:.2f} & {fmt(c['objective'])} & {residual} & PASS " + ROW_END)
    allcase.append(r'\addlinespace[3pt]')
write('family_rows.tex',rows)
write('timing_rows.tex',timings)
write('case_rows.tex',allcase)

# Each family gets a spacious page of scenario-by-scenario comparisons. Both
# the economic objective and the distinct competition workbook score are shown.
top=[]
for i,f in enumerate(DATA['families']):
    cases=[c for c in DATA['cases'] if c['buses']==f['buses']]
    if i: top.append(r'\clearpage')
    top += [rf"\subsection{{{f['buses']:,}-bus family}}",
        r'\noindent Higher objectives are better. Entries show team code, raw objective $S$, and the published scenario score $Q$ in parentheses. $Q$ is a separate workbook field, not a rounded value recomputed here. The laptop column is the verified V34 objective; no official competition score was awarded to this later experiment.',
        r'\begin{center}\small',r'\setlength{\tabcolsep}{4pt}\renewcommand{\arraystretch}{0.87}',
        r'\begin{tabular}{@{}lrlll@{}}\toprule',
        r'Scen. & V34 objective & Rank & Official team & Official $S$ (published $Q$)\\\midrule']
    for c in cases:
        for k,t in enumerate(c['official_top5'],1):
            first=c['scenario'] if k==1 else ''
            own=fmt(c['objective']) if k==1 else ''
            q=fmt(t['workbook_score'],0)
            top.append(f"{first} & {own} & {k} & {TEAMS[t['team']]} & {fmt(t['objective'])} ({q}) " + ROW_END)
        top.append(r'\addlinespace[5pt]')
    top += [r'\bottomrule\end{tabular}\end{center}',
        r'\noindent Source: Final Event results workbook, Division 1; feasible competitors only, excluding the organizer benchmark. Selection uses unrounded $S$, retaining each team\textquotesingle s best feasible row. Exact values, workbook row numbers, timing fields, and flags are retained in \texttt{data/evidence.json}.']
write('top5_tables.tex',top)

teamrows=[]
for team,abbr in sorted(TEAMS.items(),key=lambda kv:kv[1]):
    entries=[t for c in DATA['cases'] for t in c['official_top5'] if t['team']==team]
    pairs={(t['code1_resources'],t['code2_resources']) for t in entries}
    assert len(pairs)==1, (team,pairs)
    a,b=next(iter(pairs))
    teamrows.append(f"{abbr} & {esc(team)} & \\texttt{{{esc(a)}}} & \\texttt{{{esc(b)}}} " + ROW_END)
write('team_rows.tex',teamrows)

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42})
fig,ax=plt.subplots(figsize=(6.45,2.9))
colors=['#23679c','#147d64','#b85c27','#923e68','#6d65a6','#202d3a']
markers=['o','s','^','D','v','P']
for f,color,marker in zip(DATA['families'],colors,markers):
    cases=[c for c in DATA['cases'] if c['buses']==f['buses']]
    ax.scatter([c['seconds'] for c in cases],[c['shortfall_to_fifth_percent'] for c in cases],
        label=f"{f['buses']:,}",c=color,marker=marker,s=36,edgecolors='white',linewidths=.45)
ax.axvline(300,color='#a52830',linestyle='--',linewidth=1)
ax.text(282,31,'300 s',ha='right',va='top',fontsize=8,color='#a52830')
ax.set_xscale('log')
ax.set_xlim(2.2,370)
ax.set_ylim(0,32)
ax.set_xticks([3,10,30,100,300],labels=['3','10','30','100','300'])
ax.set_xlabel('Laptop end-to-end time (s, logarithmic scale)')
ax.set_ylabel('Shortfall to fifth-best\nreference objective (%)')
ax.grid(axis='y',alpha=.2)
ax.set_axisbelow(True)
ax.spines[['top','right']].set_visible(False)
ax.legend(title='Network family (buses)',ncol=3,loc='upper left',bbox_to_anchor=(0,1.04),
    fontsize=8,title_fontsize=8,frameon=False,columnspacing=1.2,handletextpad=.5)
fig.tight_layout(pad=.6)
fig.savefig(OUT/'runtime_quality.pdf')
plt.close(fig)

# Machine-check coverage rather than trusting tables assembled by hand.
assert len(DATA['cases'])==37
assert sum(c['contingencies'] for c in DATA['cases'])==44517
assert sum(len(c['official_top5']) for c in DATA['cases'])==185
assert all(c['passed'] and not c['meets_top5'] for c in DATA['cases'])
assert all(c['base_commitment_matches_source'] for c in DATA['cases'])
assert sum(c['base_nonlinear_refinement']['status']=='exception_preserved_verified_incumbent' for c in DATA['cases'])==5
print('Generated 37 result rows, 185 official comparison entries, 11 team allocations, and one vector figure.')
