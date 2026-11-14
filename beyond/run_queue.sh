#!/bin/zsh
# sequential queue: waits for the Bloch sweep, then the 3D checks (3 BLAS threads each)
export OMP_NUM_THREADS=3 OPENBLAS_NUM_THREADS=3 MKL_NUM_THREADS=3
cd ~/av2atg/scalerich/beyond
while pgrep -f "python blochsweep.py" >/dev/null; do sleep 20; done
cd d3
python pcheck3d.py  > ../results/pcheck3d.log 2>&1 || echo "pcheck3d failed" >> ../results/queue.log
python skelsweep.py > ../results/skel3d.log 2>&1 || echo "skelsweep failed" >> ../results/queue.log
python tcheck3d.py  > ../results/tcheck3d.log 2>&1 || echo "tcheck3d failed" >> ../results/queue.log
echo "queue done $(date)" >> ../results/queue.log
