import sys
sys.path.append('../build/src')
import argparse
import particlesimulator as ps
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


infile = args.checkpoint#"ewald_test"
outfile = "ewald_energy_test_temp"
T = 298
eps = 78.0
box = [10.0, 10.0, 10.0]
kx = 11
kz = 11
cutoff = 5.0
alpha = 8.0 / 10.0
sigma = 1.0 / (np.sqrt(2.0) * alpha)

ky = kx
kVec = [kx, ky, kz]

sim = ps.Simulator(eps,   T, outfile)

sim.state.set_geometry(0, box) #Have to add geometry before energy!!
sim.state.set_energy(1, [cutoff, kVec[0], kVec[1], kVec[2], alpha, False])

com, pos, charges, r, rf, b, names = load_cp(infile)
b_min = np.zeros(len(b))
b_max = np.zeros(len(b))
sim.state.load_cp(com, pos, charges, r, rf, b, b_min, b_max, names)

sim.add_move(ps.MoveTypes.TRANSLATE, [1.0, 0.0])

#sim.add_sampler(8, 1)
#sim.add_sampler(9, 1)

sim.finalize()
sim.run(1,    0,     0)

assert np.abs(sim.state.energy / 7.188995062152295 + 1.0021255) < 1e-6, "Ewald energy is not correct!"
