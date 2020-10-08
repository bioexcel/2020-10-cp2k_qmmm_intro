## Section 2: Diels-Alder reaction in vacuum

To obtain the reaction profile of the reaction in vaccum, we are going to use the [Nundged Elastic Band (NEB)](https://theory.cm.utexas.edu/henkelman/pubs/jonsson98_385.pdf) to estimate the energy profile of the Diels Alder reaction. 

We are going to use the resulting structures from the previous section as guesses for the reaction path. Since the structures were optimised using the semi-empirical PM3 method we will also optimise the BAND using PM3. 

You will find the input file here: **GMX_DAA/section2/CP2K/inp.NEB_pm3**

**Highlighted regions of inp.NED_PM3**

We have to set up calculation type as `BAND` in the `&GLOBAL` section:

```
&GLOBAL
  PROJECT DAA_NEB
  RUN_TYPE BAND          ! Band method
  PRINT_LEVEL MEDIUM
&END GLOBAL
```

We use a similar `&FORCE_EVAL` section of the input file because we want to maintain the same QM theory level.
However, we define the section `&BAND` within the `&MOTION` section which will contain parameters relavent for the NEB. 

- `  BAND_TYPE CI-NEB`
- `  NUMBER_OF_REPLICA` 15
- `  K_SPRING` 0.05
- `  &REPLICA` subsection. Here we list the previous optimised structures as initial guesses for NEB.


```
&MOTION
  &BAND                    ! BAND run
    BAND_TYPE CI-NEB       ! Type of BAND calculation - Climbing Image NEB
    NUMBER_OF_REPLICA 15   ! Number of Replica to use in the BAND
    K_SPRING 0.05          ! Spring constant for the band
    
    &CONVERGENCE_CONTROL.  ! control the convergence criteria for BAND
      MAX_FORCE 0.0010
      RMS_FORCE 0.0050
    &END
    
    ROTATE_FRAMES TRUE.    ! Compute at each BAND step the RMSD and rotate the frames in order to minimize it
    ALIGN_FRAMES TRUE      ! Alignment of the frames at the beginning of a BAND calculation
    
    &CI_NEB                ! CI-NEB type calculation only.
      NSTEPS_IT 5          ! Number of Improved Tangent steps before switching on CI algorithm
    &END
    
    &OPTIMIZE_BAND         ! optimization method for the band 
      OPT_TYPE DIIS        ! DIIS based optimization procedure for BAND
      OPTIMIZE_END_POINTS FALSE
      &DIIS
        MAX_STEPS 500      ! The number of steps to run the NEB
      &END
      
    &END
    &PROGRAM_RUN_INFO
    &END
    &CONVERGENCE_INFO
    &END

    &REPLICA               ! Coordinates of the replica
      COORD_FILE_NAME ./DA.R-pos-1.xyz
    &END
    &REPLICA
      COORD_FILE_NAME ./DA.TS-pos-1.xyz
    &END
    &REPLICA
      COORD_FILE_NAME ./DA.P-pos-1.xyz
    &END
  &END BAND

&END MOTION
```

> **TIP**: The order of the atoms in the xyz files MUST be the same in all the structures. 

> **TIP**: In a production run you would not supply MAX_STEPS and instead run until fully converged.

**Running the NEB calculation**

Submit the following jobscipt to run the calculation.

```
qsub sub-neb.pbs
```

> **TIP**: The NUMBER_OF_REPLICA and NPROC_REP in the input should ideally multiply to give the number of processes used in the job script.


**Understanding the output**


As well as the standard CP2K output in output.neb a series of outputs are generated for each replica in the band. These include:

- ``DAA_NEB-BANDXX.out`` - The geometry optimisation output for each replica.
- ``DAA_NEB-pos-Replica_nr_XX-1.xyz`` - The optimisation trajectory for each replica.
- ``DAA_NEB-r-XX.out`` - The starting point for the geometry optimisation.

By taking the final trajectory printed in each trajectory 
file  you can create an xyz file which will allow you to visualise the optimised transition trajectory.
This can be done by running the following:

```
$ for x in 01 02 03 04 05 06 07 08 09 10 12 13 14 15; do tail -n 26 DAA_NEB-pos-Replica_nr_${x}-1.xyz >> DAA-movie.xyz ; done
```
You may then view this with you chosen xyz viewer e.g. vmd:

```
$ module load vmd
$ vmd DAA-movie.xyz
```
The energy profile can be obtained from the final energies of each of the replicas.
This is printed at the end of output.neb for all replicas, and also in ``DAA_NEB-BANDXX.out`` for the individual replicas.

```
$ grep 'Total Energy' output.neb | tail -n 15
```

As before the energy is given in Hartrees. 

If we write the energy profile with respect to the energy of the reactant, the following profile is obtained:



Replica  |  Energy / kcalmol-1
------------ | ------------- 
1 (R)  | 0.0000
 2 | 0.0584
 3 | 0.1660
 4 | 0.3060
 5 | 0.4751
 6 | 0.6672
 7 | 0.8759
 8 | 1.0982
 9 | 1.3656
10 | 1.9431
11 | 4.8952
12 | 20.587
13 | -6.5748
14 | -7.7414
15 (P)  | -8.6958



![NEB_15replicas](https://github.com/salomellabres/CP2K_tutorials_for_biological_simulations/blob/master/GMX_DAA/section2/CP2K/NEB_15rep.png)

As expected, the energy of the first few images is very close that of the reactant state,
and the energy of the final images is very close to that of the product state.
The transition energy is not as high as that for our predicted transition state
-- this could result from a number of things. Firstly, the transition state used
in Section 1 was a "best first approximation" and might have been slightly off.
Secondly, this simulation shows only 13 steps in the transition (13 snapshots)
-- it's possible that the peak from out approximation can is reached somewhere 
in the transition region. To test this, you could run another nudged elastic band
simulation with *e.g.* images 10 and 12 as your starting and finishing positions.
This will give you a higher resolution for the transition region.


<br/><br/>

---