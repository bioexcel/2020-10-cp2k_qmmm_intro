## Section 3: Diels-Alder reaction in solution

In this session, we will simulate the Diels Alder reaction in solution. To do so 
we use a QM/MM approach where the reacting part of the molecule is defined 
quantum mechanically, and the rest of the molecule, as well as the solvent, are 
defined classically. This setup requires us to define link atoms. We will 
explain how to do so in this section.

This section is divided into 4 parts:
- System setup
- Classical minimisation and equilibration of the system
- Monitorisation of the QM/MM set up
- Biased simulations runs

For the system setup, minimisation & equilibration, and parts of the 
monitorisation of the QM/MM setup, we will be using the AmberTools package. We 
have recently given a course on using AmberTools to prepare a biological system 
for QM/MM treatment with CP2K. While we will be going over these technicques 
during this session, we will not be doing so in as much detail as in that 
course. If you are interested in using AmberTools, or if you want a more full 
understanding of these methods, please have a look at the
[course video](https://www.youtube.com/watch?v=zilybdb8x-A) and 
[course material](https://www.archer2.ac.uk/training/courses/200609-amber/).

<br/><br/>

### Section 3.1: Setup of the system.

**Ligand parameterisation**

We are using a simple approach to parameterise the reaction product and model 
the reaction backwards ( **P --> R** ). To do so, we use the ambertools 
protocol (using antechamber and prmchk2). For more details on the process, see 
[this link](http://ambermd.org/tutorials/basic/tutorial4b/). 

We are simulating the full molecule this time:

![Full product molecule](https://github.com/salomellabres/CP2K_tutorials_for_biological_simulations/blob/master/GMX_DAA/images/full.product.png) [JS NOTE -- change link]

```
$ antechamber -i product.pdb -fi pdb -o product.mol2 -fo mol2 -nc -1 -c bcc -at gaff2 -rn DAA
$ parmchk2 -i product.mol2 -f mol2 -o product.frcmod
```

For this protocol, we use GAFF2 forcefield parameters to define atom types and 
use AM1-BCC QM level to calculate the charges. 

<br/><br/>

**System preparation**

Using tleap from the ambertools suite, we are going to solvate the system in a 
waterbox, including a Na+ counterion to neutralise the system.

The tleap input file (**leap.in**) looks like this:
```
source leaprc.protein.ff14SB
source leaprc.water.tip3p
source leaprc.gaff2

loadamberparams product.frcmod

LIG = loadmol2 product.mol2

solvatebox LIG TIP3PBOX 15
addions LIG Na+ 1
savepdb LIG system.pdb
saveamberparm LIG system.prmtop system.inpcrd

quit
```

Finally we execute tleap like this:

```
$ tleap -f leap.in
```

> TIP: Always visualise the system to check everything is OK.

<br/><br/>

### Section 3.2: Classical minimisation and equilibration of the system.

---
**NOTE**

For this part of the exercise, we will need the `system.prmtop` and 
`system.inpcrd` files generated during **Section 3.1** of this exercise. If you 
have not managed to generate these files, a copy has been provided 
[here](./exercises/2_minimisation/bckp_missed_steps)

---

Here we are going to leverage the AMBER tools and AMBER tutorial 
[B0](http://ambermd.org/tutorials/basic/tutorial0/). We will use a similar 
protocol to the one in AMBER forcefield to minimise & equilibrate the system 
using classical molecular mechanics. We are not going to explain in detail the 
AMBER input file, you will find a detailed description in the aforementioned 
tutorial. 

For this to work, we will need to copy over the `system.prmtop` and 
`system.inpcrd` generated in the first part of this exercise. We will be using 
the `sander` tool. You will find the input file here: 

- [in.classical_minimisation](./exercises/2_minimisation/in.classical_minimisation)
- [in.classical_heating](./exercises/2_minimisation/in.classical_heating)

**in.classical_minimisation**

```
Initial minimisation of our structure
 &cntrl
  imin=1, maxcyc=4500, ncyc=2000,
  cut=8.0, ntb=1, ntc=2, ntf=2
 /
```

**in.classical_heating**

```
Heat
 &cntrl
  imin=0, ntx=1, irest=0,
  nstlim=10000, dt=0.002,
  ntf=2, ntc=2,
  tempi=0.0, temp0=300.0,
  ntpr=100, ntwx=100,
  cut=8.0,
  ntb=1, ntp=0,
  ntt=3, gamma_ln=2.0,
  nmropt=1, ig=-1,
 /
&wt type='TEMP0', istep1=0, istep2=9000, value1=0.0, value2=300.0 /
&wt type='TEMP0', istep1=9001, istep2=10000, value1=300.0, value2=300.0 /
&wt type='END' /
```

To run the minimisation, we will be sumbmitting a job on ARCHER. This can be 
done by running:

```bash
qsub sub_sander.pbs
```

The lines of interest are equivalent to running the following on the login nodes
(but running on the compute nodes is faster, and will not use up the shared 
login-node resources):

```
$ sander -O -i in.classical_minimisation -o out.classical_minimisation \
            -p system.prmtop -c system.inpcrd -r system.min.r
$ sander -O -i in.classical_heating -o out.classical_heating \
            -p system.prmtop -c system.min.r -r system.md.r \
            -x system.md.nc
```

<br/><br/>

### Section 3.3: Monitorisation of the QM/MM set up.

---

**NOTE**

For this part of the exercise, we will need the `system.prmtop` file generated 
in **Section 3.1**, and the `system.md.r` file generated in **Section 3.2** of 
this exercise. If you have not managed to generate these files, a copy has been 
provided [here](./exercises/3_monitorisation_qmmm/bckp_missed_steps)

---

There are several steps in order to set up a QM/MM system:
- Modifications of the Lennard-Jones parameters, 
- Fix the atomic coordinates format obtained form the equilibration step.
- Decide the QM/MM boundaries and set up the link atoms accordingly.

<br/><br/>

**Modification of Lennard-Jones parameters**

We need to be aware of the limitations of the classical models and how to 
combine them with the QM methods. We have set up the solvent using TIP3P water 
model, which has no Lennard-Jones parameters defined for the hydrogen atoms. We 
have to set them in order to avoid unrealistic interactions with the QM region. 
Here in this example, we are going to use the LJ parameters in the GAFF2 
forcefield.

To modify them, we are going to modify the prmtop file using parmed. Here are 
the PARMED commands to run:

```
$ parmed system.prmtop -i inp.parmed_tip3p
```

The contents of `inp.parmed_tip3p` are:

```
changeLJSingleType :WAT@H1 0.3019 0.047
changeLJSingleType :WAT@H2 0.3019 0.047
outparm system.LJ_mod.prmtop
quit
```

<br/><br/>

**Fix the atomic coordinates format obtained form the equilibration step**

In order to use the restart file from the heat simulations, we need to change 
the format from binary netcdf to a six columns formatted restart file. CPPTRAJ 
(from AMBERtools) can do that for us and simultaneosuly remove periodic boundary 
conditions artifacts.

```
$ cpptraj system.prmtop -i inp.cpptraj_reformat
```

where the contents of `inp.cpptraj_reformat` are:

```
trajin system.md.r
autoimage
trajout system.md.crd restrt
go
quit
```

> **TIP**: CRD AMBER format (10 columns) and RST7 AMBER (6 columns) format 
differ on the number of columns. We are converting the file to restrt (RST7) 
format but naming the file as .crd as CP2K requires. 

For a more detailed explanation of the cpptraj commands, you can look at the 
[CPPTRAJ manual](https://amber-md.github.io/cpptraj/CPPTRAJ.xhtml). 

<br/><br/>

**Definition of the QM/MM boundaries**

After that we set up the QM/MM partition as follows: 

![QM/MM partition](https://github.com/salomellabres/CP2K_tutorials_for_biological_simulations/blob/master/GMX_DAA/images/QMMM.section3.png)

QM region is highlighted in blue and link atoms are shown in purple. The rest 
of the system (water molecules and counterions) is shown in white. 

As shown in the scheme above, the QM/MM boundary divides the ligand and hence, 
covalent bonds. To deal with broken covalent bonds, we are going to use the 
[**IMOMM** scheme](https://www.semanticscholar.org/paper/IMOMM%3A-A-New-Integrated-Ab-Initio-%2B-Molecular-of-Maseras-Morokuma/5e1bb154312e6751c044c8226a9ee73b69db89f6). This scheme adds an additional atomic centre (usually a hydrogen atom) covalently bonded to the last QM atom. This additional atom is not real and its only purpose is to saturate its valency replacing the bond that has been broken. 

To set the link atoms in CP2K, we have to do the following steps.

- Delete the atoms that you want to exclude from the QM part.
- Identify the link atoms.
- Add the link atoms to the system under the QM kind part: 

```
   &LINK
      MM_INDEX  5
      QM_INDEX  38
      LINK_TYPE IMOMM
    &END LINK
```

Here it is an image showing all the indexes of the atoms for this solvated 
system:

![QMMM indexes](https://github.com/salomellabres/CP2K_tutorials_for_biological_simulations/blob/master/GMX_DAA/images/QMMM.section3_indexes.png)

<br/><br/>

You will find the input file here: **GMX_DAA/section3/CP2K/monitorisation_qmmm/qmmm_md.link_atoms.inp**

We have to set up calculation type as `MD` in the `&GLOBAL` section:

```
&GLOBAL 
  PROJECT QMMM
  PRINT_LEVEL LOW
  RUN_TYPE MD
&END GLOBAL
```

The ```&FORCE_EVAL``` section retains part of the definitions we used in the 
vacuum QM calculations. We have to a lot of new parameters:
- ```METHOD QMMM``` : We will use method QMMM
- Similar ```&DFT``` section
- ```&MM``` section: Here we specify all the necessary options to run the MM 
- region of the system. 
  - ```&FORCEFIELD``` section
  - ```&POISSON``` section
- ```&SUBSYS``` section: List all the details of the system:
  - ```&CELL``` section: Size of the system. 
  - ```&TOPOLOGY``` section: Here we specify the AMBER input files. 
  - ```&KIND NA+``` section to specify all the elements that are defined in the 
  - system with a different element name (Na -- NA+)
- ```&QMMM``` section: 
  - ```ECOUPL COULOMB``` electrostatic coupling scheme
  - ```&CELL``` section: size of the QM subsystem. 
  - ```&QM_KIND``` section: We list all the indexes of each element.  
  - ``` &LINK``` section to list the link atom indexes.

```
&FORCE_EVAL
  METHOD QMMM
  STRESS_TENSOR ANALYTICAL
  &DFT
    CHARGE 0
    &QS
      METHOD PM3
    &END QS
    &SCF
      MAX_SCF 2000
      EPS_SCF 1.0E-6
      SCF_GUESS ATOMIC
      &OT
        MINIMIZER DIIS
        PRECONDITIONER FULL_SINGLE_INVERSE
      &END
      &OUTER_SCF
        EPS_SCF 1.0E-6
        MAX_SCF 100
      &END
    &END SCF
  &END DFT
    &MM
    &FORCEFIELD
      PARMTYPE AMBER
      PARM_FILE_NAME post.LJ_mod.prmtop
      &SPLINE
        EMAX_SPLINE 1.0E8
        RCUT_NB [angstrom] 10
      &END SPLINE
    &END FORCEFIELD
    &POISSON
      &EWALD
        EWALD_TYPE SPME
        ALPHA .40
        GMAX 80
      &END EWALD
    &END POISSON
  &END MM
  &SUBSYS
    &CELL
    !Set box dimensions here
      ABC [angstrom] 49.9302020  41.8911420  42.8511700
      ALPHA_BETA_GAMMA 90 90 90
    &END CELL
    &TOPOLOGY
      CONN_FILE_FORMAT AMBER
      CONN_FILE_NAME post.LJ_mod.prmtop
      COORD_FILE_FORMAT CRD
      COORD_FILE_NAME post.md.crd
    &END TOPOLOGY
    &KIND NA+
     ELEMENT Na
    &END KIND
  &END SUBSYS
    &QMMM
    ECOUPL COULOMB
    &CELL
      ABC 40 40 40
      ALPHA_BETA_GAMMA 90 90 90
    &END CELL
    &QM_KIND N
      MM_INDEX 33
    &END QM_KIND
    &QM_KIND O
      MM_INDEX 16 17 35 37
    &END QM_KIND
    &QM_KIND C
      MM_INDEX 27 28 29 30 31 32 34 36 38
    &END QM_KIND
    &QM_KIND H
      MM_INDEX 19 20 21 22
    &END QM_KIND
    &QM_KIND CL
      MM_INDEX 23 24 25 26
    &END QM_KIND
    &QM_KIND S
      MM_INDEX 18
    &END QM_KIND
    &LINK
      MM_INDEX  5
      QM_INDEX  38
      LINK_TYPE IMOMM
    &END LINK
  &END QMMM
&END FORCE_EVAL
```

The main changes are located in the ```&MOTION``` section:

- ```&MD subsection``` : Define all the parameter to run a MD. Including:
  - ```ENSEMBLE``` : NPT_I, NVT...
  - ```&BAROSTAT``` subsection
  - ```&THERMOSTAT``` subsection
  - ```&PRINT``` subsection

```
&MOTION
  &MD
  ENSEMBLE NPT_I
  TIMESTEP [fs] 0.5
  STEPS    20000
  TEMPERATURE 298
  &BAROSTAT
    TIMECON [fs] 100
    PRESSURE [bar] 1.0
  &END BAROSTAT
  &THERMOSTAT
    REGION GLOBAL
    TYPE CSVR
    &CSVR
      TIMECON [fs] 10.
    &END CSVR
  &END THERMOSTAT
  &END MD
  &PRINT
    &RESTART                                    ! This section controls the printing of restart files
      &EACH                                     ! A restart file will be printed every 10000 md steps
        MD 25000
      &END
    &END
    &TRAJECTORY                                 ! Thes section Controls the output of the trajectory
      FORMAT DCD                                ! Format of the output trajectory is DCD
      &EACH                                     ! New trajectory frame will be printed each 100 md steps
        MD 4
      &END
    &END
    &RESTART_HISTORY                            ! This section controls dumping of unique restart files during the run keeping all of them.Most useful if recovery is needed at a later point.
      &EACH                                     ! A new restart file will be printed every 10000 md steps
        MD 5000
      &END
    &END
    &CELL
      &EACH
        MD 100
      &END
    &END
  &END PRINT
&END MOTION
```

<br/><br/>

### Section 3.4: QM/MM enhanced sampling (Metadynamics)

---
**NOTE**

For this part of the exercise, we will need the `system.LJ_mod.prmtop` and the 
`system.md.crd` files generated in **Section 3.3** of this exercise. If you have 
not managed to generate these files, a copy has been provided 
[here](./exercises/4_qmmm_metadynamics/bckp_missed_steps).

---

As we saw in the previous QM/MM simulations, the system remains estable and 
runs smoothly. However, for the Diels Alder reaction to be reversed we should 
extend our equilibrium simulations for a very long time and probably nothing 
would happen given the energy difference between reactants and products. 
Therefore, we need to biase the system in order to sample the Diels Alder 
reaction. To biase the system we need to define collective variables that 
describe the reaction or change we want to observe and use enhanced sampling 
methods to push the system onto that coordinate. 

CP2K has different [collective variables (CV)](https://manual.cp2k.org/trunk/CP2K_INPUT/FORCE_EVAL/SUBSYS/COLVAR.html) 
and several [enhanced sampling methods](https://manual.cp2k.org/trunk/CP2K_INPUT/MOTION/FREE_ENERGY.html) 
implemented. However, for this particular problem none of the implemented CV 
fulfil our needs. Luckily, we can used the 
[PLUMED2 plugin](https://www.plumed.org) to define a custom CV. To do so, you 
will need an additional file (plumed.dat) describing your CV to CP2K.

**plumed.dat**:

``` 
COM ATOMS=30,31 LABEL=com1
COM ATOMS=29,32 LABEL=com2

DISTANCE ATOMS=com1,com2 LABEL=d1
UPPER_WALLS ARG=d1 AT=6.0 KAPPA=150.0 EXP=2 EPS=1 OFFSET=0 LABEL=uwall

#
# Activate metadynamics in d1
# depositing a Gaussian every 500 time steps,
# with height equal to 1.2 kJ/mol,
# and width 0.35 rad for both CVs.
#
METAD ARG=d1 PACE=500 HEIGHT=1.2 SIGMA=0.35,0.35 FILE=HILLS LABEL=metad

# monitor the d1 distance, the upper wall and the metadynamics bias potential
PRINT STRIDE=10 ARG=d1,uwall,metad.bias FILE=COLVAR
```

**qmmm_md.link_atoms.metadynamics.inp**

In the ```&GLOBAL``` section goes as: 

```
&GLOBAL
  PROJECT METAD
  PRINT_LEVEL LOW
  RUN_TYPE MD
&END GLOBAL
```

In the ```&MOTION``` section add: 

```
&FREE_ENERGY
  &METADYN
    USE_PLUMED .TRUE.
    PLUMED_INPUT_FILE  ./plumed.dat
  &END METADYN
&END FREE_ENERGY
  ```

You can also run metadynamics using the collective variables defined in 
[CP2K](https://manual.cp2k.org/trunk/CP2K_INPUT/FORCE_EVAL/SUBSYS/COLVAR.html). 
Here I add a sample version of the ```&FREE_ENERGY``` section:

```
  &FREE_ENERGY
    METHOD METADYN
    &METADYN
      DO_HILLS  .TRUE.
      NT_HILLS 100
      WW [kcalmol] 1.5
      &METAVAR
        WIDTH 0.5 !Also known as scale
        COLVAR 1
      &END METAVAR
    &PRINT
        &COLVAR
           COMMON_ITERATION_LEVELS 3
           &EACH
             MD 10
           &END
        &END
      &END
    &END METADYN
  &END FREE_ENERGY
```

And the COLVAR definition is added to the ```&MM``` section of the 
```&FORCE_EVAL```:

```
    &COLVAR
       &DISTANCE_FUNCTION
          COEFFICIENT 1
          ATOMS 29 30 32 31
       &END DISTANCE_FUNCTION
    &END COLVAR
```

<br/><br/>

---