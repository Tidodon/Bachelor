# Bachelor project code
My transport code for time-dependent current through model systems.
Run program by calling:
$ python transport_code.py 
Specify system description. Be aware that gamma parameter has to be finetuned.

When program is executed it will ask user to specify a list of parameters. I will give here an explanation for what each of them means.
I wrote the code with my own usage in my mind and therefore some of the names used might be hard to understand or only understandable for me.


Size of molecule matrix: 7 
How many atom sites molecule between electrodes is built up of
Size of molecule-lead matrices: 50
Number of sites in each electrode that are taken be a part of 'extended molecule'. 
Size of electrode matrices: 400
Number of bath states
Is molecule benzene? yes or no: no
If molecule is not benzene then molecule will be a polyene/molecular wire. If benzene option is chosen, user will be asked whether it is ortho, meta or para benzene. If Number of of molecular atom sites written to 'Size of molecule matrix: ' is more than 6, then output will be a benzene that is substitutued with a polyene chain in one side.

Hückel theory hopping matrix elements for different system subsections:
Size of direct molecule-lead coupling: -0.3
Insert Molecule-Lead coupling V: -0.3
Above refers to coupling between Molecule-Lead, which is a part of electrode and a given electrode. This value is most probably the same as hopping matrix elements in leads.
Insert hopping matrix element in leads: -0.3
Insert hopping matrix element in molecule: -0.3
Insert bias voltage in V: 0.3
Insert fermi energy in eV: 0
Two options above, fermi energy and bias voltage, come together to constitute electrochemical potential in leads and molecule. In current implementation fermi energy is taken to be the same (equal to zero) in all system subsections.
Insert electronic temperature: 0
Temperature of leads, it steers thermal equilibrium, and intial guess density matrices.
flat, ramp or poisson potential: poisson
Type of potential. Naming can be confusing here. By 'flat' I mean that extended molecule or molecule(depending on system partitioning) is filled up to Fermi level, and electrostatic potential in this system subsection is taken to be equal to zero. By 'ramp' I mean that there is a linear drop in voltage from one lead to another. It results in different local voltages felt by each atomic sites in the inner system section. By 'poisson', I do not mean a potential, but a procedure, where Poisson equation is solved self-consitentely each timestep.
System splitting in 3 or 5 parts. Write 3split or 5split: 5split
Do I split system in Lead - extended molecule - Lead or Lead - molecule lead-molecule - molecule lead - Lead. The difference lays in whether molecule or extended molecule is filled with potential type specified above.
Molecular permativity: 3
Macroscopic quantity, dielectric constant
Interatomic distance in atomic units: 3
In current implementation it is taken to be equal for all atoms in the system. Of course, in reality metal and molecule bond lengths cannot be taken to be identical.
Start Potential type flat or ramp: flat
This option is relevant if self-consistent Poisson equation solving procedure is chosen. In this case a starting guess is specified.
Number of steps: 2000
How many timesteps are to be made.
Size of driving factor: 0.03
driving term, broadening term, damping term. This factor goes by many names. It has to be finetuned. for each specific system. Its value dependes on electrodes' level spacing and effective coupling between molecule and  electrodes'.
Insert non-zero frequency in fs^-1 unit if AC current is desired, else insert 0.0:
Alternating current, it has to be in fs^-1 unit, because calculation are performed in fs units.
Calculation type: Poissonstflat
Above is a user specified name, which will be a part of output files.
plot or test: plot
If 'plot' option is chosen plots, and text files with data will output. If 'test' option is chosen a plot of current vs time will be shown. Not output files.
Produce heatmap across extended molecule. yes or no: no
If 'yes' is written a calculation of current between each of extended molecule atomic sites will calculated. 

When program is done several outputs are produced:
-Currents*.png, plot of bond currents vs time, measured at the center of molecular bridge. If number of molecular sites is uneven(fx. 7) then measurment is taken between site 3 and 4.
-Coherence*.png, plot of relevant coherences from density matrix, used to calculate currents in plot above.
-Sum_e_mol*.png, plot of number of electrons on molecule
-Sum_e_em*.png, plot of number of electrons on extended molecule, both this and above were used to check for mistakes
-ortho*.png, heatmap of currents vs site vs time
-Currents*.txt, a text file where 1 row is naming: Currents, pnn(relevant coherence), p_diag(population of site from which current is measured), sum_e_mol(number of electrons on molecule), sum_e_em(number of electrons on extended molecule). Each subsequent row is data on each of those things in columns. Each row's values are taken after 1 fs timestep.
-Currents_vs_em*.txt, Each column is current between extended molecule sites, first column is current between site 0 and 1, second is current betweeen 1 and 2 and so on. Each row's values are taken each 1 fs timestep.
-Energies_state*.txt, System energies in 'state' representation(in terms of a system subsections eigenvectors)
-Populations_site*.txt, Each row: Number of electrons on each site(columns), each row's values are taken each 100 timesteps.
-Populations_state*.txt, Each row: Number of electrons pr. quantum level of the system, each row's value rare taken each 100 timesteps.


What program can:
-1D system description.
-different lengths of molecular wires and ortho, meta and para benzenes.
-Specification of electronic lead temperature, Fermi level, hopping matrix elements(the same values inside each of 5 system subsection), 
-alternating current, electrostatic potential in form of a linear drop from one electrode to another and selfconsistent Poisson equation solving. Division in extended molecule + leads or molecule, molecule-lead sections and 2 leads.
-System is 1D in its current form and calculations can be performed on Hückel theory level of theory. 
