import subprocess
import os

DOCKING_PATH = '.'

EXTRA_ATOMS = {'Si': "atom_par Si     4.10  0.200  35.8235  -0.00143  0.0  0.0  0  -1  -1  6",
               'B': "atom_par B      3.84  0.155  29.6478  -0.00152  0.0  0.0  0  -1  -1  0",
               'Pt': "atom_par Pt     2.75  0.080  12.000   -0.00110  0.0  0.0  0  -1  -1  4",
               'Ni': "atom_par Ni     2.83  0.015  12.000   -0.00110  0.0  0.0  0  -1  -1  4"
               }

RDPF = """
parameter_file extra_atoms.dat
autodock_parameter_version 4.2       # used by autodock to validate parameter set
outlev 1                             # diagnostic output level
intelec                              # calculate internal electrostatics
seed pid time                        # seeds for random generator
ligand_types HD C A N NA OA F P SA S Cl Br I
fld r.maps.fld                       # grid_data_file
map r.HD.map
map r.C.map
map r.A.map
map r.N.map
map r.NA.map
map r.OA.map
map r.F.map
map r.P.map
map r.SA.map
map r.S.map
map r.Cl.map
map r.Br.map
map r.I.map
elecmap r.e.map                      # electrostatics map
desolvmap r.d.map                    # desolvation map
move l.pdbqt                         # small molecule
about x y z                          # small molecule center
tran0 random                         # initial coordinates/A or random
quaternion0 random                   # initial orientation
dihe0 random                         # initial dihedrals (relative) or random
torsdof 7                            # torsional degrees of freedom
rmstol 2.0                           # cluster_tolerance/A
extnrg 1000.0                        # external grid energy
e0max 0.0 10000                      # max initial energy; max number of retries
ga_pop_size 150                      # number of individuals in population
ga_num_evals 2500000                 # maximum number of energy evaluations
ga_num_generations 27000             # maximum number of generations
ga_elitism 1                         # number of top individuals to survive to next generation
ga_mutation_rate 0.02                # rate of gene mutation
ga_crossover_rate 0.8                # rate of crossover
ga_window_size 10                    # 
ga_cauchy_alpha 0.0                  # Alpha parameter of Cauchy distribution
ga_cauchy_beta 1.0                   # Beta parameter Cauchy distribution
set_ga                               # set the above parameters for GA or LGA
sw_max_its 300                       # iterations of Solis & Wets local search
sw_max_succ 4                        # consecutive successes before changing rho
sw_max_fail 4                        # consecutive failures before changing rho
sw_rho 1.0                           # size of local search space to sample
sw_lb_rho 0.01                       # lower bound on rho
ls_search_freq 0.06                  # probability of performing local search on individual
set_psw1                             # set the above pseudo-Solis & Wets parameters
unbound_model bound                  # state of unbound ligand
ga_run n                             # do this many hybrid GA-LS runs
analysis                             # perform a ranked cluster analysis
"""

def make_pdbqt_gpf(rc, lg, extra_atoms=[]):
    command = f"mk_prepare_receptor.py --read_pdb {rc} -o r -p -v -g --box_enveloping {lg} --padding 5 -a --default_altloc A"
    subprocess.run(command)

def edit_rgpf(rgpf):
    pass
