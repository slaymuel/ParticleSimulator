# Particle Simulator
ParticleSimulator simulates simple particles in bulk and close to flat surfaces.

## Features

- Simulate simple particles using Markov Chain Monte Carlo (MCMC)
- Support for several thermodynamic ensembles (NVT, NPT and μVT)
- Simulate particles confined by infinitely conducting surfaces
- Several different analysis types
- Parallelization with [OpenMP]
- Python bindings with [pybind11]

## Example usage
Using pybind11, running a simulation is very easy.


```sh
import particlesimulator as ps #Import the package

#Create a simulator object
eps = 78.0 #Relative permittivity
temperature = 298
sim = ps.Simulator(eps,   temperature, outfile)

#Create a geometry (has to be done before setting the energy)
geometryType = 0 #Cubic box
box = [10.0, 10.0, 10.0] #Box dimensions 
sim.state.set_geometry(geometryType, box)

#Set the energy method
energyType = 1 #Ewald summation
cutoff = 5.0 #Distance cutoff when calculating the energy
alpha = 8.0 / 10.0 #Convergence parameter for the Ewald method
kVec = [5, 5, 5] #Reciprocal wave vectors for the Ewald method
sim.state.set_energy(energyType, [cutoff, kVec[0], kVec[1], kVec[2], alpha, False])

#Use the helper method 'load_cp' to load a checkpoint file (see below)
com, pos, charges, r, rf, b, b_min, b_max, names = load_cp(infile)
#Load the particle coordinates etc
sim.state.load_cp(com, pos, charges, r, rf, b, b_min, b_max, names)

#Add Monte Carlo moves
sim.add_move(mr.MoveTypes.TRANSLATE, [1.0, 0.0]);
#Add samplers for analysis
samplerType = 1 #Particle density in the z-dimension
sim.add_sampler(samplerType, 100)

#Have to call finalize() to tell the program that the settings are finished
sim.finalize()
#Run the simulation
microSteps = 10000
macroSteps = 1000
equilibrationSteps = 10
sim.run(macroSteps,    microSteps,     equilibrationSteps)
```
The above script is using the load_cp helper function:
```sh
import numpy as np

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
    b_min = np.array(data[:,10], dtype=np.float64)
    b_max = np.array(data[:,11], dtype=np.float64)
    name = list(data[:,12])
    
    return com, pos, charge, r, rf, b, b_min, b_max, name
```

## License

MIT

[//]: # (These are reference links used in the body of this note and get stripped out when the markdown processor does its job. There is no need to format nicely because it shouldn't be seen. Thanks SO - http://stackoverflow.com/questions/4823468/store-comments-in-markdown-syntax)
 [pybind11]: <https://github.com/pybind/pybind11>
 [OpenMP]: <https://www.openmp.org/>
