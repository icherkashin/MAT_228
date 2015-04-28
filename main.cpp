/** 
 * \mainpage MAT 228A Computing Homework
 * 
 * \par Description
 * 
 * Implementation of the conservative finite difference methods for MAT 228A computing homework. 
 * 
 * \par Computation
 * C++
 * 
 * \par Data storage
 * hdf5
 * 
 * \par Data processing
 * Python + h5py + matplotlib
 * 
 * \par Installation
 * 
 * <ol>
 *  <li>
 *   Install libhdf5-dev, python-h5py and python-matplotlib (e.g. using "apt-get install")
 *  </li>
 *  <li>
 *   Adjust the parameters in the CMakeLists.txt file according to your installation directory and location of hdf5 library
 *  </li>
 *  <li>
 *   Run "cmake ." command, which should prepare the makefile
 *  </li>
 *  <li>
 *   Run "make" command, which should compile the "compute_task" executable
 *  </li>
 * </ol>
 * 
 * \par Usage
 *  
 * The syntax for executing a particular computational task is:
 * 
 * <ul><li>./compute_task "Initial Condition" "Flux" "Refinement Exponent"</li></ul>
 * 
 * 
 * The sequence of tasks is programed in the file "execute_tasks.py" using python syntax and functions 
 * defined in the file "lib_output_processing.py."
 * 
 * An example script can be found in the file "execute_task.py" which performs computations for all initial
 * conditions and fluxes and processes all data. 
 * 
 * Based on this example, the user can define behavior that suits his particular needs.
 *
 * \param Initial_Condition Valid initial condition input string @see initial_conditions
 * 
 * \param Flux Valid flux input string @see fluxes
 * 
 * \param Refinement_Exponent Valid grid stepsize; refinement exponent is related to grid stepsize as h = 2^(-Refinement_Exponent)
 *  
 * \par Analysis of the computational results:
 * 
 * <a href="../../../index.html">Theory and Computing Homework</a>
 *
 * \par Author:  
 * @author
 * */

#include <map>
#include "Compute_task.hpp"

int main(int argc, char **argv) {

// Process and validate arguments  
  std::map<std::string, std::string> arguments = process_arguments(argc, argv);

// Select the flux based on the input task and execute the computation  
  if (arguments["flux"] == "Upwind") 
   main_loop<Upwind>(arguments); 
  else if (arguments["flux"] == "Lax_Friedrichs") 
   main_loop<Lax_Friedrichs>(arguments); 
  else if (arguments["flux"] == "Lax_Wendroff") 
   main_loop<Lax_Wendroff>(arguments); 
  else if (arguments["flux"] == "Fromm" or arguments["flux"] == "Fromm_CFL_half") 
   main_loop<Fromm>(arguments); 
  else if (arguments["flux"] == "Fromm_van_Leer" or arguments["flux"] == "Fromm_van_Leer_CFL_half") 
   main_loop<Fromm_van_Leer>(arguments); 
  else if (arguments["flux"] == "Flux_Corrected_Transport") 
   main_loop<Flux_Corrected_Transport>(arguments); 
  else if (arguments["flux"] == "Lax_Wendroff_Fourth_Order") 
   main_loop<Lax_Wendroff_Fourth_Order>(arguments);     
} 