# simulations
Code for generating final state for given decays. 

# running the code 

Use ROOT to run the example codes:

```
root -l generate.C
```

To save the output into lund file use command"

```
root -l -q -e ".x phiproduce.C " | grep -v '^$' > data.lund
```

#Fast Monte-Carlo simulations for CLAS12

The new fastMC package offers event analysis for CLAS12 detector
responses. It can be run on LUND files to mark particles that will
be registered by different detectors. To run fastMC the clas12 
environment must be activated on farms and then run the code
to analyze the LUND file.

```
module load clas12
java -jar /home/gavalian/Software/fastMC/fastMC.jar input.lund output.lund
```

If no output lund file is given, and automatically generated output name
will be used, by appending "_fastmc.lund" to the end of input file name.

In the output file, the 6-th column of the particle list will be changed
to reflect detector responses for the particle, the number 0 will be 
changed to numbers 5,6,7 or 8, with following convention

```
0 - no detector was hit
5 - only DC (all 6 superlayers were hit)
6 - the DC (all 6 superlayers), and FTOF was hit (no ECAL)
7 - the DC and ECAL were hit (no FTOF)
8 - the DC, FTOF and ECAL all were hit by the particle
```


# code output

The code outputs the simlated particles and their decays.
The header prints out the generated angles for production 
process and decay process.

first two numbers are the theta and phi angles in the rest 
frame of the reaction, for proton. where the ep decays to a 
proton and a vector meson.

the second two numbers are the theta and phi angles in of the 
decay products in the rest frame of the vector meson.

The header is followed by the printoutof all the vectors in the 
reaction:

```
header: production [  2.96388,   2.01860], decay [  1.23560   1.61794]
e  in:   0.00000   0.00000  10.60000 [  0.00050]
e out:   1.18395   0.00000   6.40598 [  0.00050]
p out:  -0.10248  -0.12368   0.58169 [  0.93827]
decay:  -1.08148   0.12368   3.61233 [  1.02000]
dgt 1:  -0.58369   0.18683   1.97271 [  0.49368]
dgt 2:  -0.49778  -0.06315   1.63962 [  0.49368]
```

# running for LUND output

run command:

```
root -l -q -b rhoproduce.C
```

or 

```
root -l -q -b phiproduce.C
```

the script prints out events in LUND format on the screen. the header contains
q2,xb, production theta and phi, and decay theta and phi, starting from column 11;
