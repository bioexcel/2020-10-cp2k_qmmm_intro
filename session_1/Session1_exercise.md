## Section 1: Single point calculations of Diels-Alder reaction in vacuum

First, we are going to perform geometry optimisation calculations on 3 
structures: 

- Reactants (**R**), 
- Transition state (**TS**) and
- Products (**P**).

We have provided these 3 structures in XYZ format. These has been inferred 
from the X-ray structure (PDB id 1C1E) and previously optimised using 
Gaussian09. 

<br/>

### 1.0 Basics of CP2K

Below is a description of the main parts of our CP2K input file, with certain 
parts of interest highlighted. In general, CP2K input files are split into three 
main sections. The 
[GLOBAL](https://manual.cp2k.org/cp2k-7_1-branch/CP2K_INPUT/GLOBAL.html) section 
allows you to define general information about the simulation to be run (*e.g.* 
the project name, the type of simulation to be run, the amount of output 
required, the simulation random seed, etc.). The 
[FORCE_EVAL](https://manual.cp2k.org/cp2k-7_1-branch/CP2K_INPUT/FORCE_EVAL.html) 
section is where the simulation parameters (forces, energies, atom types, *etc*.) 
are defined. The final section will depend on the type of simulation you are 
trying to run -- simulations where you expect movement (*e.g.* geometry 
optimisations) will need a 
[MOTION](https://manual.cp2k.org/cp2k-7_1-branch/CP2K_INPUT/MOTION.html) section 
to define the way that particles can move whereas simulations focussing more on 
vibrational modes of freedom or phononics will need a
[VIBRATION_ANALYSIS](https://manual.cp2k.org/cp2k-7_1-branch/CP2K_INPUT/VIBRATIONAL_ANALYSIS.html) 
section. The list of available sections can be found 
[here](https://manual.cp2k.org/cp2k-7_1-branch/CP2K_INPUT.html).

Each of these sections then has a number of subsections that can be used to 
fully define the system being run. In CP2K, a section called `SECTION` with a 
subsection called `SUBSECTION` where a variable `VARIABLE` needs to be defined, 
would look like this:

```bash
&SECTION
  &SUBSECTION
    VARIABLE set_variable
    ...
  &END SUBSECTION
&END SECTION
```

This is recursive, and in a few cases, you will need to define sub-subsections 
*etc.*. Comments can be added using an exclamation mark `!`, and these can be 
set as their own line or following some CP2K instructions (though be aware that 
CP2K instructions that require string inputs will not read `!` as a comment, 
which will cause an irregular termination of the simulation).

<br/><br/>

### 1.1 Setting up a geometry optimisation of the chemical structures of the reactants state

First, we are going to perform a quantum mechanical (QM) geometry optimisation 
for the reactant and product states. We are going to use the [semi-empirical PM3 
method](https://en.wikipedia.org/wiki/PM3_(chemistry))
.

You will find the input files here: 
- [Reactant geometry-optimisation input file](./exercise/1_reactant/inp.reactant_pm3_geoopt)
- [Transition state geometry optimisation input file](./exercise/2_transition_state/inp.transition_state_pm3_geoopt)
- [Transition state frequency calculations input file](./exercise/2_transition_state/inp.transition_state_pm3_freq)
- [Product geomerty optimisation input file](./exercise/3_product/inp.product_pm3_geoopt)

#### 1.1.1 Highlighted regions of inp.reactant_pm3_geoopt

We have to set up calculation type as `GEO_OPT` in the `&GLOBAL` section:

```
&GLOBAL
  PROJECT DA.R        ! Name of the calculation
  PRINT_LEVEL LOW     ! Verbosity of the output
  RUN_TYPE GEO_OPT    ! Calculation type: Geometry optimisation
&END GLOBAL
```

In the `&FORCE_EVAL` section, we define the basic system definitions (such as 
topology and coordinates). Here we are only going to explain the most important 
ones:
- `METHOD QS` : QUICKSTEP is the QM method in CP2K ([link](https://www.cp2k.org/quickstep)).
- In the `&SUBSYS` subsection, we define several parameters of the system. 
  - `&CELL` defines the simulation box size that will contain all QM atoms.
  - `&COORD` defines the starting coordinates of the QM atoms. 
- In the `&DFT` subsection, we define several parameters of the density 
functional theory (DFT) basis set
  - `&CHARGE 0` sets the overall charge of the system.
  - `&QS` subsection where QUICKSTEP parameters are set. Amongst other things, 
we need to specify which QS method we are using (in this case `METHOD PM3`).
  - `&MGRID` subsection where all the parameters for calculating the Gaussian 
plane waves are defined.
  - `&SCF` subsection where parameters for finding a self-consistent solution 
(SCF) of the 
[Kohn-Sham](https://en.wikipedia.org/wiki/Kohn%E2%80%93Sham_equations) DFT 
formalism are defined.

```
&FORCE_EVAL                              ! parameters needed to calculate energy and forces
  METHOD QS                              ! QUICKSTEP method
  &SUBSYS                                ! a subsystem: coordinates, topology, molecules and cell
    &CELL                                ! Input parameters needed to set up the simulation cell
      ABC 12.4138 12.4138 12.4138
      PERIODIC NONE                      ! Define non-periodic boundary conditions
    &END CELL
    &COORD                               ! Coordinates for simple systems specified using explicit XYZ coordinates
  ... R coordinates ...
    &END COORD
  &END SUBSYS  
  &DFT                                   ! Parameter needed by linear combination of atomic orbitals (LCAO) DFT programs
    CHARGE 0                             ! The total overall charge of the system
    &QS                                  ! parameters needed to set up the QUICKSTEP framework
      METHOD PM3                         ! Electronic structure method 
      &SE                                ! Parameters needed to set up the semi-empirical (SE) methods
         &COULOMB                        ! Evaluation of the COULOMB term in SE calculations
           CUTOFF [angstrom] 10.0
         &END
         &EXCHANGE                       ! Evaluation of the EXCHANGE and core Hamiltonian terms in SE calculations
           CUTOFF [angstrom] 10.0
         &END
      &END
    &END QS
    &MGRID                               ! Multigrid information
      CUTOFF 200
      NGRIDS 4
      REL_CUTOFF 30
    &END MGRID
    &SCF                                 ! Parameters needed to perform an SCF run.
      SCF_GUESS ATOMIC
      EPS_SCF 1.0E-05
      MAX_SCF 2000
      &OT                                ! Orbital transformation (OT) method
        MINIMIZER DIIS
        PRECONDITIONER FULL_SINGLE_INVERSEinp.product_geoopt
      &END
      &OUTER_SCF                         ! Parameters controlling the outer SCF loop
        EPS_SCF 1.0E-6
        MAX_SCF 100
      &END
      &PRINT
        &RESTART OFF
        &END RESTART
      &END PRINT
    &END SCF
  &END DFT
&END FORCE_EVAL
```

In the `&MOTION` section, we define the parameters for the geometry 
optimisation. 

```
&MOTION                        
  &GEO_OPT                    ! Environment of the geometry optimizer
    TYPE MINIMIZATION         ! Kind of geometry optimization -- in this case, MINIMIZATION
      MAX_ITER 4000
    OPTIMIZER CG
    &CG                       ! Conjugate gradient (CG) optimization
      MAX_STEEP_STEPS  0
      RESTART_LIMIT 9.0E-01
    &END CG
  &END GEO_OPT
&END MOTION
```
#### 1.1.2 Running CP2K on ARCHER

ARCHER uses modules to store centrally-installed software, and CP2K is one of 
these centrally-installed softwares. Any ARCHER user can load CP2K by running 
`module load cp2k/7.1`. This gives users access to the CP2K executables (as well 
as making sure that all of the correct libraries and dependencies are loaded). 
The ARCHER CP2K module has two executable: `cp2k.sopt` for single-processor jobs 
and `cp2k.popt` for multi-processor jobs. Through this course, we will be using 
the multi-processor `cp2k.popt` executable. Note that, as the ARCHER login nodes 
are shared with other users, running multi-processor jobs is not enabled on them, 
and we will need to submit these simulations to the compute nodes.

We have made sure that all exercises have a script to allow you to submit CP2K 
simulations to the compute nodes with ease -- these sripts can be found in the 
appropriate `exercise` folders, and are called `sub_NAME.pbs`, where the name 
may have been altered to provide a better description. You can submit your 
simulation with the command:

```bash
$ qsub sub_NAME.pbs
```

and you can check on whether your job is queueing with the command:

```bash
$ qstat -u $USER
```

#### 1.1.3 Analysing the results

Once the job has finished, we need to check that the minimisation has converged.
For each optimisation step, we a section that looks similar to this: 

```
  Core-core repulsion energy [eV]:                         111210.31809538803645
  Core Hamiltonian energy [eV]:                             -6070.07521304926922
  Two-electron integral energy [eV]:                      -218128.32380280233338
  Electronic energy [eV]:                                 -115134.23711445042863

  Total energy [eV]:                                        -3923.91901906240992

  Atomic reference energy [eV]:                              3926.78210430983836
  Heat of formation [kcal/mol]:                                66.02431952621455

  outer SCF iter =    6 RMS gradient =   0.72E-06 energy =       -144.2013768850
  outer SCF loop converged in   6 iterations or    8 steps


 ENERGY| Total FORCE_EVAL ( QS ) energy (a.u.):             -144.201376885183862


 *******************************************************************************
 ***                 BRENT   - NUMBER OF ENERGY EVALUATIONS :       1        ***
 *******************************************************************************

 --------  Informations at step =   197 ------------
  Optimization Method        =                   SD
  Total Energy               =      -144.2013768852
  Real energy change         =        -0.0000000005
  Decrease in energy         =                   NO
  Used time                  =                0.857

  Convergence check :
  Max. step size             =         0.0000000000
  Conv. limit for step size  =         0.0030000000
  Convergence in step size   =                  YES
  RMS step size              =         0.0000000000
  Conv. limit for RMS step   =         0.0015000000
  Convergence in RMS step    =                  YES
  Max. gradient              =         0.0079886843
  Conv. limit for gradients  =         0.0004500000
  Conv. for gradients        =                   NO
  RMS gradient               =         0.0015348400
  Conv. limit for RMS grad.  =         0.0003000000
  Conv. for gradients        =                   NO
 ---------------------------------------------------
```

The number of steps it takes for the system to converge can vary. Here we can 
see that it took 197 steps to converge this system. CP2K evaluates the 
convergence using 4 different criteria:

- The maximum calculated displacement for the next step must be essentially 0:

    **Conv. limit for step size**  =        4.50000000E-004
    
    - This default value can be changed by setting `&MOTION` `&GEO_OPT` 
`MAX_DR` to a value of your choosing.

- The root-mean-square (RMS) of the calculated displacement for the next step 
must be also essentially 0:

    **Conv. limit for RMS step**   =        1.50000000E-003
    
    - This default value can be changed by setting `&MOTION` `&GEO_OPT` `RMS_DR`
to a value of your choosing.

- The maximum forces/gradient must be essentially 0: 

    **Conv. limit for gradients**  =        3.00000000E-003
    
    - This default value can be changed by setting `&MOTION` `&GEO_OPT` 
`MAX_FORCE` to a value of your choosing.

- The RMS of the forces/gradient must be also essentially 0:

    **Conv. limit for RMS grad.**  =        3.00000000E-004
    
    - This default value can be changed by setting `&MOTION` `&GEO_OPT` 
`RMS_FORCE` to a value of your choosing.


> **Tip:** Use this command to monitor the runs: grep -a12 "Convergence check" 
out.reactants_pm3_geoopt

#### 1.1.4 Exercise

Now that we have the energy of the geometrically-optimised reactants, we need 
to find the same for the product. You will find all of the files needed in the 
2nd part of the Session 1 exercises. Have a look at 
[inp.product_geoopt](./exercise/2_product/inp.product_geoopt) to ensure that you 
understand what the script is doing and, when you are ready, launch the 
simulation using the submission script provided.

What do you find the total energy of the system to be?

<br/><br/>

### 1.2 Geometry optimisation and Frequency calulations of the transition states

#### Geometry optimisation of the Transition State

Next, we are going to perform a QM geometry optimisation for the transition 
state. We are going to use semi-empirical PM3 again to obtain consistent 
results. The CP2K input file is similar to the geometry optimisation for the 
reactant and product states, but there are several differences:

You will find the input file 
[here](./exercise/3_transition_state/inp.transition_state_pm3_geoopt).

**Highlighted regions of inp.transition_state_pm3_geoopt**

To calculate the geometry for the TS, we need to slightly modify the `&GEO_OPT` 
in the `&MOTION` section. We will use the [DIMER method](https://aip.scitation.org/doi/10.1063/1.480097) 
to calculate the transition states in CP2K. 

```
&MOTION
  &GEO_OPT                          ! Environment of the geometry optimizer
    TYPE TRANSITION_STATE           ! Kind of geometry optimization
    MAX_ITER 2000
    OPTIMIZER CG
    &CG                             ! Conjugate gradient optimization
      MAX_STEEP_STEPS 1000
    &END
    &TRANSITION_STATE               ! Transition state search
      METHOD DIMER                  ! Method to use for locating transition states
      &DIMER                        ! parameters for Dimer Method
        DR 0.0001
        ANGLE_TOLERANCE [deg] 0.5
        INTERPOLATE_GRADIENT
        &ROT_OPT
          OPTIMIZER CG
          MAX_ITER 1000
          &CG
            MAX_STEEP_STEPS 1000
          &END CG
        &END
      &END
    &END
  &END GEO_OPT
&END MOTION
```

<br/><br/>

#### Frequency calculation

Additionally, we need to do an extra calculation to validate the TS structure. 
We need to make sure that there are negative frequencies on the atomic 
vibrations and that they correspond to the formation of the bonds. For the 
vibrational analysis to be correct, we need to use the same QM method used in 
the geometry optimisation.

You will find the input file 
[here](./exercise/3_transition_state/inp.transition_state_freq).

**Highlighted regions of inp._transition_state_freq**

This input script is slightly different from the one for geometric optimisation 
(see above). Here, we highlight the differences between the two.

In `&Global`, we have to set up calculation type as `NORMAL_MODES` instead of 
`GEO_OPT`:

```
&GLOBAL
  PROJECT DA.TS.freq       ! Name of the calculation
  PRINT_LEVEL MEDIUM       ! Verbosity of the output
  RUN_TYPE NORMAL_MODES    ! Calculation type: Normal modes
&END GLOBAL
```

We maintain the same `&FORCE_EVAL` section seen in the geometry optimisations 
runs. We omit the `&MOTION` section and replace it with the following section 
instead:

```
&VIBRATIONAL_ANALYSIS     ! Needed for normal modes, vibrational, or phonon analysis runs.
  NPROC_REP 1             ! Number of processors to be used per replica environment.
  DX 0.01                 ! Increment to construct the HESSIAN using the finite difference method.
  &PRINT
    &MOLDEN_VIB
    &END
    &CARTESIAN_EIGS
    &END
    &PROGRAM_RUN_INFO
      &EACH
        REPLICA_EVAL 1
      &END
    &END
  &END
&END
```

This calculation creates a series of output files. Of note are:

- **TS-VIBRATIONS-1.mol** : Resulting frequencies and structures. 
- **TS-VIBRATIONS-1.eig** : Binary file containing the eigenvectors and 
eigenstates of the frequencies.
- **TS-Hessian-1.hess**   : Binary file containing the Hessian matrix.

In the **DA.TS.freq-VIBRATIONS-1.mol** file, we can see that the frequencies 
listed at the top are negative. The only negative frequency corresponds to the 
formation of the two C-C bonds in a synchonous way. 

```
 [Molden Format]
 [FREQ]
     -605.332321
      -73.720060
      -11.031533
       34.453243
       42.252911

       ...
```

The easiest way to validate the atomic vibrations is to visualise them. The 
default output format in CP2K is MOLDEN format, you might find this software 
difficult to install and to use -- we have therefore provided a python script 
that converts this format to multiple XYZ files, one for each vibration. You can 
visualise these multiple XYZ files using Pymol of VMD.

This image shows the first vibration found using CP2K:

[DA TS vibrations](../Images/session1/vibration.png)

<br/><br/>

### 1.3 Calculation of the energy barrier.

Here we compare the single-point energy calculations between Gaussian and CP2K. 
To do this, we take the final energy output from the simulation and compare it 
to the value obtained for the reactant. Note that CP2K outputs energies in 
arbitraty units (a.u.) -- these are equivalent to Hartree units, and 1 Hartree 
= 627.509 kcal/mol. Here are the results we obtain: 

Energy / kcal/mol | Reactants | Transition States | Product
------------ | ------------- | ------------- | -------------
Structure | [DA R](../Images/session1/DA_reactants.png) | [DA TS](../Images/session1/DA_transition.png) | [DA P](../Images/session1/DA_product.png)
CP2K rel | 0.00 | 47.81 | -8.70
Gaussian09 | 0.00 | 47.81 | -8.69

These results are the same as, or very similar to, those found using Gaussian09.

<br/><br/>

---
