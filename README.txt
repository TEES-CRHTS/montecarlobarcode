Jesse Fine
jfine@tamu.edu
September 12, 2021

This codebase runs MCMatlab, which was created and is distributed by Marti et al. 


SpectraWorkbooks- contains excel sheets of spectra used in Monte Carlo simulations



Phosphorescence Lifetime Decay- These evaluate the photoluminescent output of the near infrared and visible spectrum sensing domains functionalized by phosphorescent dyes.
1) NIRPHOS_20210912
2) VISPHOS_20210912

FRET- These evaluate the photoluminescent output of the near infrared and visible spectrum sensing domains functionalized by the FRET assay described in this work. 
1)NIRFRET_20210912
2)VISFRET_20210912

Autofluorescence

The autofluorescence folder has 2 main files- "Autofluorescence_450nm_20210912.m" and "Autofluorescence_680nm_20210912.m" that determine autofluorescence contribution from each included fluorophore.
The individual contribution is simulated using the appropriate file within the autofluorescence folder. i.e. "Autofluorescence_450nm_FAD_May4.m" was used to determine FAD's autofluorescence contribution when excited at 450nm. 