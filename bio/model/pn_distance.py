"""Total-variation distance between model P(n) (results/plant_predictions_N2000.json) and the
shoot-meristem P(n) (results/meristem_topo.json via figures.log values). Run from bio/."""
import json
d = dict(zip(range(3, 11), [0.003, 0.041, 0.267, 0.408, 0.23, 0.046, 0.003, 0.0]))  # fig_plant, 12182 obs
pr = json.load(open('results/plant_predictions_N2000.json'))
for k, v in pr.items():
    m = {int(a): b for a, b in v['Pn'].items()}
    print(f'{k:35s} TV = {0.5 * sum(abs(m.get(n, 0) - d.get(n, 0)) for n in range(2, 16)):.3f}')
