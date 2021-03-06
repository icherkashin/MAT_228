#!/usr/bin/env python3

import lib_output_processing as output
import subprocess

# Create HDF5 databse if it does not exist
output.create_output_database()

# Empty task list
task = ["", ""]

# Run computations for the selected initial conditions and grid resolutions
for initial_condition in output.initial_conditions:
 for flux in output.fluxes:
  task[0] = initial_condition
  task[1] = flux
  # Practical range is from k = 6 to k = 17
  for k in range(6,10):
   subprocess.call(["./compute_task", *task, str(k)])

# Process the results for the selected initial conditions and grid resolutions
for initial_condition in output.initial_conditions:
 for flux in output.fluxes:
  task[0] = initial_condition
  task[1] = flux
  
  output.check_conservation(*task)
  output.create_plots(*task)
  output.create_convergence_table(*task)
  output.output_html_content(*task)
  output.include_content_html()
