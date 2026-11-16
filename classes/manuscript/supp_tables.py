"""Tables for supplement.tex from results/{mech3,geo,snap}.json -> supp_tables.tex"""
import json, numpy as np
R = '../results/'
m3 = json.load(open(R + 'mech3.json')); geo = json.load(open(R + 'geo.json')); snap = json.load(open(R + 'snap.json'))
NM = {'golden': 'golden, aligned', 'golden-avoid': 'golden, anti', '(5,1)-avoid': '(5,1), anti', '(3,2)-avoid': '(3,2), anti', 'uniform': 'uniform split'}
out = []
# Table S1: T = 1e4 (uniform: mean), EB
out.append(r'\begin{table}[h]\caption{Rigid-joint beam model at $T=10^{4}$ (uniform split: one realisation; over three realisations at $T=3000$ the relative spread is below $2\%$ in $\langle n\rangle$ and $3\%$ in $s$). $\langle n\rangle$: direction-averaged index; $s=\mathrm{tr}\,C^{\rm pin}/\mathrm{tr}\,C$; $s_E$: axial part of the rigid-joint strain energy.}\label{tab:T4}')
out.append(r'\begin{ruledtabular}\begin{tabular}{lccccccccc}')
out.append(r' & \multicolumn{3}{c}{$\bar\rho=0.05$} & \multicolumn{3}{c}{$\bar\rho=0.1$} & \multicolumn{3}{c}{$\bar\rho=0.2$} \\')
out.append(r'rule & $\langle n\rangle$ & $s$ & $s_E$ & $\langle n\rangle$ & $s$ & $s_E$ & $\langle n\rangle$ & $s$ & $s_E$ \\ \hline')
for nm, lab in NM.items():
    row = [lab]
    for rho in (0.05, 0.1, 0.2):
        o = [o for o in m3 if o['rule'] == nm and o['T'] == 10000 and o['rho_t'] == rho and o['model'] == 'EB'][0]
        row += ['%.2f' % o['n_mean'], '%.3f' % o['pin_ratio'], '%.2f' % o['s']]
    out.append(' & '.join(row) + r' \\')
out.append(r'\end{tabular}\end{ruledtabular}\end{table}')
# Table S2: size dependence at rho = 0.1
out.append(r'\begin{table}[h]\caption{Size dependence at $\bar\rho=0.1$: rigid-joint $\langle n\rangle$, pinned ceiling $n_{\max,x}$ and stretching share $s$ (uniform split: mean of three realisations for $T<10^4$).}\label{tab:size}')
out.append(r'\begin{ruledtabular}\begin{tabular}{lcccccccc}')
Ts = (300, 1000, 3000, 10000)
out.append(r'rule & ' + ' & '.join(r'\multicolumn{2}{c}{$T=%d$}' % T for T in Ts) + r' \\')
out.append(' & ' + ' & '.join([r'$\langle n\rangle$, $n_{\max,x}$ & $s$'] * 4) + r' \\ \hline')
for nm, lab in NM.items():
    row = [lab]
    for T in Ts:
        o = [o for o in m3 if o['rule'] == nm and o['T'] == T and o['rho_t'] == 0.1 and o['model'] == 'EB']
        row += ['%.2f, %.2f' % (np.mean([q['n_mean'] for q in o]), np.mean([q['npin_x'] for q in o])), '%.3f' % np.mean([q['pin_ratio'] for q in o])]
    out.append(' & '.join(row) + r' \\')
out.append(r'\end{tabular}\end{ruledtabular}\end{table}')
# Table S3: Timoshenko
out.append(r'\begin{table}[h]\caption{Relative change of $\langle n\rangle$ from Euler--Bernoulli to Timoshenko beams ($\nu=0.3$, $\kappa=5/6$), maximum over $T=300$--$10^{4}$.}\label{tab:timo}')
out.append(r'\begin{ruledtabular}\begin{tabular}{lcc}rule & $\bar\rho=0.2$ & $\bar\rho=0.4$ \\ \hline')
for nm, lab in NM.items():
    row = [lab]
    for rho in (0.2, 0.4):
        ch = []
        for o in m3:
            if o['rule'] == nm and o['rho_t'] == rho and o['model'] == 'timo':
                e = [q for q in m3 if q['rule'] == nm and q['T'] == o['T'] and q['seed'] == o['seed'] and q['rho_t'] == rho and q['model'] == 'EB'][0]
                ch.append(o['n_mean'] / e['n_mean'] - 1)
        row.append('%.1f\\%%' % (100 * max(ch)))
    out.append(' & '.join(row) + r' \\')
out.append(r'\end{tabular}\end{ruledtabular}\end{table}')
# Table S4: geometry
out.append(r'\begin{table}[h]\caption{Geometry in units of $1/L_{\rm tot}$, so that a length $\ell$ equals $\ell L_{\rm tot}/\bar\rho$ thicknesses. Segments between consecutive joints: minimum and median; fraction shorter than $3\lambda$ at $\bar\rho=0.1$ and $0.2$. Offsets $\delta$ between opposite stems on a host: minimum and first percentile.}\label{tab:geo}')
out.append(r'\begin{ruledtabular}\begin{tabular}{lcccccccc}')
out.append(r'rule & $T$ & $\ell_{\min}$ & $\ell_{\rm med}$ & $f_{<3\lambda}(0.1)$ & $f_{<3\lambda}(0.2)$ & $\delta_{\min}$ & $\delta_{1\%}$ & longest$_2$ \\ \hline')
for o in geo:
    lab = NM[o['rule']]
    dq = o['delta_q']
    out.append(' & '.join([lab, r'$10^{%d}$' % round(np.log10(o['T'])), '%.2f' % o['seg_q'][0], '%.2f' % o['seg_q'][3],
                           '%.3f' % o['frac_seg_lt3']['0.1'], '%.3f' % o['frac_seg_lt3']['0.2'],
                           '%.2f' % dq[0], '%.2f' % dq[1], '%.3f' % o['longest2']]) + r' \\')
out.append(r'\end{tabular}\end{ruledtabular}\end{table}')
# Table S5: snap
out.append(r'\begin{table}[h]\caption{Finite-thickness test: largest $\bar\rho$ (step $0.05$) at which merging parallel ligaments closer than $\lambda$ leaves $m=1$.}\label{tab:snap}')
out.append(r'\begin{ruledtabular}\begin{tabular}{lcccc}rule & $T=10^{3}$ & $10^{4}$ & $10^{5}$ & $10^{6}$ \\ \hline')
for nm in ('golden-avoid', '(5,1)-avoid', '(3,2)-avoid'):
    row = [NM[nm]]
    for T in (1000, 10000, 100000, 1000000):
        o = [o for o in snap if o['rule'] == nm and o['T'] == T][0]
        k = 0
        while k < len(o['m_snap']) and o['m_snap'][k] == 1: k += 1
        row.append('%.2f' % o['rho'][k - 1])
    out.append(' & '.join(row) + r' \\')
out.append(r'\end{tabular}\end{ruledtabular}\end{table}')
open('supp_tables.tex', 'w').write('\n'.join(out) + '\n')
print('ok')
