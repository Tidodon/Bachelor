# Bachelor
My transport code for time-dependent current through model systems.
Run program by calling:
$ python transport_code.py 
Specify system description. Be aware that gamma parameter has to be finetuned.

What program can:
1D system description.
different lengths of molecular wires and ortho, meta and para benzenes.
Specification of electronic lead temperature, Fermi level, hopping matrix elements(the same values inside each of 5 system subsection), 
AC current, ESP potential in form of ramp and selfconsistent Poisson. Division in extended molecule + leads or molecule , molecule-lead sections and 2 leads.
System is 1D in its current form and calculations are intended to be on Hückel theory level of theory. 
