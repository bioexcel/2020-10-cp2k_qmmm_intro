## Section 2: Diels-Alder reaction in vacuum

To obtain the reaction profile of the reaction in vaccum, we are going to use the [Nundged Elastic Band (NEB)](https://theory.cm.utexas.edu/henkelman/pubs/jonsson98_385.pdf) to estimate the energy profile of the Diels Alder reaction. 

We are going to use the resulting structures from the previous section as guesses for the reaction path. Since the structures were optimised using the semi-empirical PM3. 

You will find the input file here: **GMX_DAA/section2/CP2K/inp.NEB_pm3**

**Highlighted regions of inp.NED_PM3**

We have to set up calculation type as `BAND` in the `&GLOBAL` section:

```
&GLOBAL
  PROJECT DAA_NEB
  RUN_TYPE BAND          ! Band methods
  PRINT_LEVEL MEDIUM
&END GLOBAL
```

We use a similar `&FORCE_EVAL` section of the input file because we want to maintain the same QM level. However, we define a new `&MOTION` section. 

- `BAND_TYPE CI-NEB`
- `NUMBER_OF_REPLICA` 15
- `K_SPRING` 0.05
- `&REPLICA` subsection. Here we list the structures as initial guesses for NEB.

```
&MOTION
  &BAND                    ! BAND run
    BAND_TYPE CI-NEB       ! Type of BAND calculation
    NUMBER_OF_REPLICA 15   ! Number of Replica to use in the BAND
    K_SPRING 0.05          ! Spring constant
    
    &CONVERGENCE_CONTROL.  ! control the convergence criteria for BAND
      MAX_FORCE 0.0010
      RMS_FORCE 0.0050
    &END
    
    ROTATE_FRAMES TRUE.    ! Compute at each BAND step the RMSD and rotate the frames in order to minimize it
    ALIGN_FRAMES TRUE      ! Alignment of the frames at the beginning of a BAND calculation
    
    &CI_NEB                ! CI-NEB type calculation only.
      NSTEPS_IT 5
    &END
    
    &OPTIMIZE_BAND         ! optimization method for the band 
      OPT_TYPE DIIS        ! DIIS based optimization procedure for BAND
      OPTIMIZE_END_POINTS TRUE
      &DIIS
        MAX_STEPS 500
      &END
      
    &END
    &PROGRAM_RUN_INFO
    &END
    &CONVERGENCE_INFO
    &END

    &REPLICA               ! Coordinates and velocities (possibly) of the replica
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

> **TIP**: The order of the atoms MUST be the same on all the structures. 

The energy profile obtained is the following: 

Nodes  |  Energy / kcalmol-1
------------ | ------------- 
R  | 0.0000
 1 | 0.0584
 2 | 0.1660
 3 | 0.3060
 4 | 0.4751
 5 | 0.6672
 6 | 0.8759
 7 | 1.0982
 8 | 1.3656
 9 | 1.9431
10 | 4.8952
11 | 20.587
12 | -6.5748
13 | -7.7414
P  | -8.6958



![NEB_15replicas](https://github.com/salomellabres/CP2K_tutorials_for_biological_simulations/blob/master/GMX_DAA/section2/CP2K/NEB_15rep.png)

As expected, the energy of the first few images is very close that of the reactant state, and the energy of the final images is very close to that of the product state. The transition energy is not as high as that for our predicted transition state -- this could result from a number of things. Firstly, the transition state used in Section 1 was a "best first approximation" and might have been slightly off. Secondly, this simulation shows only 13 steps in the transition (13 snapshots) -- it's possible that the peak from out approximation can is reached somewhere in the transition region. To test this, you could run another nudged elastic banc simulation with *e.g.* images 10 and 12 as your starting and finishing positions. This will give you a higher resolution for the transition region.


<br/><br/>

---