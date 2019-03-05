# trajectory-analysis-MATLAB
Collection of MATLAB codes to analyse unwinding trajectories

With x,y,z tracked data tracked from output experimental images this readme can be followed to analyse this raw data to determine experimental first pasage time distributions and to find the best model parameters using maximum likelihood estimation.

i)	The x and y positions of each microsphere is selected from a single frame of the LUT data.

ii)	Using these positions the MATLAB software tracks each microsphere to provide the x, y, and z positions for every microscope image frame exported during the experiment.
a)	Initially, only the data representing the 7pN to 0pN to 7pN force alteration is tracked, skipping every 10 frames, to establish which microspheres have a single DNA molecule attached of correct length.
b)	The microspheres that are not correct length or exhibit unexpected behaviour are removed from the list of microsphere positions to be tracked.

iii)	The x, y, and z position is tracked for each correct microsphere on every image frame exported during the experiment.

iv)	Find the starting point before unwinding has taken place by analysing the constant force pre-unwinding data.
a)	Code zl_make_z_only_yyyymmdd.m reads in raw x,y,z data and outputs z only
b)	Code zl_ref_subtract_yyyymmdd.m reads in z data, subtracts reference microsphere position, and outputs reference subtracted z position for all microspheres.
c)	Code zl_filter_yyyymmdd.m reads in reference subtracted microsphere positions and calculates the mean to find an initial microsphere height for each DNA molecule.

v)	Code make_z_only_yyyymmdd.m reads in raw unwinding x,y,z data and outputs z only trajectories.

vi)	Code ref_subtract_yyyymmdd.m reads in z only trajectories, subtracts reference microsphere position, and outputs reference subtracted z positions for all microspheres.

vii)	Code filter_yyyymmdd.m reads in reference subtracted microsphere positions, and takes a sliding window mean to obtain filtered z trajectories.

viii)	Code zero_level_yyyymmdd.m reads in the reference subtracted and filtered data and outputs the same data corrected to start at zero.

ix)	Code to_bp_unwound_yyyymmdd.m reads the initial start position corrected data and converts to base pairs unwound.

x)	Code bp_all_to_individual_cropped_yyyymmdd.m reads in the unwinding trajectories in base pairs unwound and truncates them to ignore the microsphere after it has disappeared. A separate file is output for each microsphere containing the unwinding trajectory.

xi)	Code linear_analysis_yyyymmdd.m performs a linear fit and calculates mean velocity.

xii)	Code fptd_yyyymmdd.m reads in the individual trajectories and calculates the first-passage times.

xiii)	Code mle_fptd_yyymmdd uses maximum likelihood estimation to find the model parameters to describe the data and performs bootstrapping to estimate errors.
