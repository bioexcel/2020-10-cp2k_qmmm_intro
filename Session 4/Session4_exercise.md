## Section 4: Diels-Alder reaction in protein

The protocol used here is very similar to the protocol described for the QM/MM system in solution. The steps are the in order to set up a QM/MM system:

- Modifications of the Lennard-Jones parameters,
- Fix the atomic coordinates format obtained form the equilibration step.
- Decide the QM/MM boundaries and set up the link atoms accordingly.

<br/><br/>

### Section 4.1: Setup of the system. 

**Cleaning the initial structure (from PDB database)**

We start with the catalytic antibody structure (PDB id 1C1E) and we are goingto build the system shown below:

![Protein-ligand system](https://github.com/salomellabres/CP2K_tutorials_for_biological_simulations/blob/master/GMX_DAA/images/QMMM.section4.png)

First, you split the system into 3 parts (protein, ligand, waters). See **fix_1c1e.sh**:

```
# Split water molecules, protein and ligand ENH
grep "HOH" 1c1e.pdb > 1c1e.xray_waters.pdb
grep -v "HOH" 1c1e.pdb | grep -v "ENH" > 1c1e.protein.pdb
grep "ENH" 1c1e.pdb > 1c1e.ligand.pdb
```

Secondly, you use pymol to align the Diels Alder product to the ligand ENH and save the new coordinates (Result: **DAA.aligned.pdb**). 

<br/><br/>

**Ligand parameterisation**

Again, we will use the ambertools protocol (using antechamber and prmchk2). For more details on the process: 
http://ambermd.org/tutorials/basic/tutorial4b/

```
$ antechamber -i DAA.aligned.pdb -fi pdb -o product.mol2 -fo mol2 -nc -1 -c bcc -at gaff2 -rn DAA
$ parmchk2 -i product.mol2 -f mol2 -o product.frcmod
```

Using this protocol, we are going to use GAFF2 forcefield parameters to define atomtupes and use AM1-BCC QM level to calculate the charges. 

<br/><br/>

**System preparation**

Using tleap for the ambertools suite, we are going to prepare the system (which inclused the product solvated in a waterbox with a counterion to neutralise the system. 

The tleap input file (**leap.in**) looks like this:

```
source leaprc.protein.ff14SB
source leaprc.water.tip3p
source leaprc.gaff2

loadamberparams product.frcmod

PROT = loadpdb 1c1e.protein.pdb
LIG = loadmol2 product.mol2
AGUAS = loadpdb 1c1e.xray_waters.pdb

SYS = combine {PROT LIG AGUAS}

solvatebox SYS TIP3PBOX 15
addions SYS Cl- 5
savepdb SYS system.pdb
saveamberparm SYS system.prmtop system.inpcrd

quit
```

Finally we execute tleap like this:

```
$ tleap -f leap.in
```

> TIP: Always visualise the system to check everything is OK.

<br/><br/>


### Section 4.2: Classical minimisation and equilibration of the system.

Here we are going to lto use the same protocol described in the **Section 3.2**. Here are the input files:

- GMX_DAA/section4/CP2K/classical_equilibration/min_classical.in
- GMX_DAA/section4/CP2K/classical_equilibration/heat_classical.in

**min_classical.in**

```
Initial minimisation of our structure
 &cntrl
  imin=1, maxcyc=4500, ncyc=2000,
  cut=8.0, ntb=1, ntc=2, ntf=2
 /
 ```
 
**heat_classical.in**

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

To run the minimisation, you have to use this command:
```
$ sander -O -i min_classical.in -o min_classical.out -p system.prmtop -c system.inpcrd -r system.min.r
$ sander -O -i heat_classical.in -o heat_classical.out -p system.prmtop -c system.min.r -r system.md.r -x system.md.nc
```

<br/><br/>

### Section 4.3: Monitorisation of the QM/MM set up.

There are several steps in order to set up a QM/MM system:
- Modifications of the Lennard-Jones parameters, 
- Fix the atomic coordinates format obtained form the equilibration step.
- Decide the QM/MM boundaries and set up the link atoms accordingly.

<br/><br/>

**Modification of Lennard-Jones parameters**

Again we need to modify the ennard-Jones parameters of TIP3P water model and the hydrogen atom of the hydroxyl groups of the Ser and Tyr residues of the protein. We are going to use the LJ parameters in the GAFF2 forcefield again.

To modify them, we are going to modify the prmtop file using parmed. Here are the PARMED commands to run:

```
$ parmed system.prmtop
```

```
changeLJSingleType :TYR@HH 0.3019 0.047
changeLJSingleType :SER@HG 0.3019 0.047
changeLJSingleType :WAT@H1 0.3019 0.047
changeLJSingleType :WAT@H2 0.3019 0.047
outparm system.LJ_mod.prmtop
quit
```

<br/><br/>

**Fix the atomic coordinates format obtained form the equilibration step**

In order to use the output of the heat simulations, we need to convert the structure from netcdf (binary) to CRD (AMBER) format. We are going to use the CPPTRAJ tool ( provided by the AMBERtools free suite ) to conver the format and recenter the water box. 

```
$ cpptraj system.prmtop
```

```
trajin system.md.r
autoimage
trajout system.md.crd restrt
go
quit
```

For a more detailed explanation of the cpptraj commands, you can look at the [CPPTRAJ manual](https://amber-md.github.io/cpptraj/CPPTRAJ.xhtml). 

<br/><br/>

**Definition of the QM/MM boundaries**

We are going to use the same QM MM partition used in **section 3** and hence, we have to set up the appropiate link atoms. You need to check the atom indexes for each system. 

```
   &LINK
      MM_INDEX  6616
      QM_INDEX  6649
      LINK_TYPE IMOMM
    &END LINK
```

<br/><br/>

The CP2K input resembles the one used in the **Section 3.3**. The index numbers have changed and will need to be updated: **GMX_DAA/section4/CP2K/monitorisation_qmmm/qmmm_md.link_atoms.inp**

```
    &QMMM
    ECOUPL COULOMB
    &CELL
      ABC 40 40 40
      ALPHA_BETA_GAMMA 90 90 90
    &END CELL
    &QM_KIND N
      MM_INDEX 6644
    &END QM_KIND
    &QM_KIND O
      MM_INDEX 6648 6646 6628 6627
    &END QM_KIND
    &QM_KIND C
      MM_INDEX 6638 6639 6640 6641 6642 6643 6649 
    &END QM_KIND
    &QM_KIND H
      MM_INDEX 6630 6631 6632 6633
    &END QM_KIND
    &QM_KIND CL
      MM_INDEX 6634 6635 6636 6637 
    &END QM_KIND
    &QM_KIND S
      MM_INDEX 6629
    &END QM_KIND
    &LINK
      MM_INDEX 6616
      QM_INDEX  6649
      LINK_TYPE IMOMM
    &END LINK
  &END QMMM
```

<br/><br/>


### Section 4.4: QM/MM enhanced sampling (Metadynamics)

As we saw in the previous QM/MM simulations, the system remains estable and runs smoothly. However, we need to set up biased simulations for the Diels Alder reaction to happen. We are going to use the same protocol described before in **section 3.4**

**plumed.dat**:

```
COM ATOMS=6640,6643 LABEL=com1
COM ATOMS=6641,6642 LABEL=com2

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

<br/><br/>
