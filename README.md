# Phanerozoic_decarbonation

author: Evan J. Ramos

#
publication: Ramos, E.J., Lackey, J.S., Barnes, J.D., and Fulton, A.A., 2020. Remnants and rates of metamorphic decarbonation in continental arcs. GSA Today.

#
In this repository, you will find all MATLAB scripts and functions that you can use to replicate analgoue model results found in the above publication. The main scripts are found in the files decarb_Monte_Carlo_script.m and decarb_Monte_Carlo_SNB.m, respectively. The former generates global and North American metamorphic decarbonation rates through the Phanerozoic, using rock information from MACROSTRAT, and the latter computes metamorphic decarbonation rates for the Cretaceous Sierra Nevada batholith using rock type information from local stratigraphy. These scripts call the function decarb_Monte_Cristo.m to compute the metamorphic decarbonation fluxes.

# NOTE 1
Our analogue decarbonation model relies upon numerical model predictions for aureole volumes as a function of intrusion volumes. The numerical model predictions were generated using the script GSAToday_aureole_prediction.m and are post processed in GSAToday_aureole_postprocess, which we have included in this repository. ALL OTHER CODE that these data rely on can be found in another repository at the URL github.com/ejramos/skarn_model.

# NOTE 2
All other functions or text files are used in the scripts. Feel free to tinker with the codes and metadata, but
some unbeknownst errors may arise in part because of this.

# NOTE 3
All code (in scripts and functions) include line-by-line comments on what the line (or section of code) accomplishes.

# NOTE 4
If you have any comments, queries, or suggestions, either contact me directly through GitHub or send me an email at ejramos@utexas.edu with the subject line "Ramos et al 2020 GSA Today model question".
