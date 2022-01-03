import sys
sys.path.append('/Users/samuel/Documents/mormon/build/src')
import argparse
import mormon as mr
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('--checkpoint')
args = parser.parse_args()

def load_cp(fileName):
    with open(fileName) as file:
        data = file.readlines()

    data = map(lambda x: x.strip().split(), data)
    data = np.array(list(data))

    com = np.array(data[:,0:3], dtype=np.float64)
    pos = np.array(data[:,3:6], dtype=np.float64)
    charge = np.array(data[:,6], dtype=np.float64)
    r = np.array(data[:,7], dtype=np.float64)
    rf = np.array(data[:,8], dtype=np.float64)
    b = np.array(data[:,9], dtype=np.float64)
    name = list(data[:,10])
    return com, pos, charge, r, rf, b, name


infile = args.checkpoint
outfile = "hw_energy_test_temp"
T = 298.0
eps = 78.3
box = [200.0, 200.0, 145.0]
cutoff = 72.5
kVec = [5, 5, 5]
alpha = np.pi / cutoff


sim = mr.Simulator(eps,   T, outfile)

sim.state.set_geometry(2, box)
sim.state.set_energy(2, [cutoff, kVec[0], kVec[1], kVec[2], alpha])

com, pos, charges, r, rf, b, names = load_cp(infile)
b_min = np.zeros(len(b))
b_max = np.zeros(len(b))
sim.state.load_cp(com, pos, charges, r, rf, b, b_min, b_max, names)

sim.add_move(0, [1.0, 0.0]);

sim.finalize()
sim.run(1,    0,     0)

assert np.abs(sim.state.energy + 352.0973589574392) < 1e-12, "Halfwald energy is not correct!"
