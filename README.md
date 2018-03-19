# trajectory-analysis-MATLAB
Collection of MATLAB codes to analyse unwinding trajectories

With x,y,z tracked data tracked from output experimental images this readme can be followed to analyse this raw data to determine experimental first pasage time distributions and to find the best model parameters using maximum likelihood estimation.

i) Find the starting point before unwinding has taken place by analysing the constant force pre-unwinding data. a) Code ‘zl_make_z_only_yyyymmdd.m’ reads in raw x,y,z data and outputs z only a) Code ‘zl_ref_subtract_yyyymmdd.m’ reads in z data, subtracts reference microsphere position, and outputs reference subtracted z position for all microspheres. a) Code ‘zl_filter_yyyymmdd.m’ reads in reference subtracted microsphere positions and calculates the mean to find an initial microsphere height for each DNA molecule.

ii) Code ‘make_z_only_yyyymmdd.m’ reads in raw unwinding x,y,z data and outputs z only trajectories.

iii) Code ‘ref_subtract_yyyymmdd.m’ reads in z only trajectories, subtracts reference microsphere position, and outputs reference subtracted z positions for all microspheres.

iv) Code ‘filter_yyyymmdd.m’ reads in reference subtracted microsphere positions, and takes a sliding window mean to obtain filtered z trajectories.

v) Code ‘zero_level_yyyymmdd.m’ reads in the reference subtracted and filtered data and outputs the same data corrected to start at zero.

vi) Code ‘to_bp_unwound_yyyymmdd.m’ reads the initial start position corrected data and converts to base pairs unwound.

vii) Code ‘bp_all_to_individual_cropped_yyyymmdd.m’ reads in the unwinding trajectories in base pairs unwound and truncates them to ignore the microsphere after it has disappeared. A separate file is ouput for each microsphere containing the unwinding traejctory.

viii) Code ‘linear_analysis_yyyymmdd.m’ performs a linear fit and calculates mean velocity.

ix) Code ‘fptd_yyymmdd.m’ reads in the individual trajectories and calculates the first-passage times.

x) Code ‘fit_fptd_3_ps_yyymmdd.m’ uses maximum likelihood estimation to find the model parameteres to describe the data and performs bootstrapping to estimate errors.
