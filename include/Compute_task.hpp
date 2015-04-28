#include <iostream>
#include <cmath>
#include <set>
#include <map>
#include <valarray>
#include <algorithm>
#include <chrono>

#include "boost/math/special_functions/sign.hpp"

#include "H5Cpp.h"

#include "Fluxes.hpp"

/** @brief Valid initial condition input strings.  */
const std::set<std::string> initial_conditions
  {
    "Square_Wave", 
    "Semicircle", 
    "Gaussian_Pulse"    
  };

/** @brief Valid fluxes input strings. */
const std::set<std::string> fluxes 
  {
    "Upwind", 
    "Lax_Friedrichs", 
    "Lax_Wendroff", 
    "Fromm", 
    "Fromm_CFL_half", 
    "Fromm_van_Leer", 
    "Fromm_van_Leer_CFL_half",
    "Flux_Corrected_Transport", 
    "Lax_Wendroff_Fourth_Order"    
  };

  /** 
   * \brief Process and validate input arguments.
   * 
   * \param argv[1] A string containing the name of the initial condition
   * \param argv[2] A string containing the name of the flux
   * \param argv[3] Grid refinement exponent, which is related to grid stepsize h as h = 2^(-Refinement_Exponent)
   * 
   * @return A map (set of key-value pairs) containing validated initial condition name, flux name, and grid refinement exponent.
   */
  
std::map<std::string, std::string> process_arguments(int& argc, char ** & argv) {

// A stack to store one or several error messages that may occur    
  std::string error_messages_stack;

// Error message in case of invalid initial condition input string
    auto valid_initial_conditions = [&] () -> void 
      { 
	std::cout << "###\tERROR:\t" << "\"" << argv[1] << "\" isn't valid initial condition input!" << std::endl;
	
	std::cout << "Valid initial conditions are:" << std::endl;   
	std::for_each(initial_conditions.begin(), initial_conditions.end(), [](std::string initial_condition){std::cout << "\t \""+ initial_condition + "\"" << std::endl;});	
	std::cout << std::endl;
	
	error_messages_stack.append("\n\tInvalid initial condition input");
      };
      
// Error message in case of invalid flux input string     
    auto valid_fluxes = [&] () -> void 
      { 
	std::cout << "###\tERROR:\t" << "\"" << argv[2] << "\" isn't valid flux input!" << std::endl;
	
	std::cout << "Valid fluxes are:" << std::endl;
	std::for_each(fluxes.begin(), fluxes.end(), [](std::string flux){std::cout << "\t \""+ flux + "\"" << std::endl;});	
	std::cout << std::endl;
	
	error_messages_stack.append("\n\tInvalid flux input");
      };

// Error message in case of invalid refinement exponent input 
    auto valid_refinement_exponent_range = [&] () -> void 
      { 
	std::cout << "###\tERROR:\t" << "\"" << argv[3] << "\" isn't valid grid refinement exponent!" << std::endl;
      	
	std::cout << "Valid refinement exponent values are:" << std::endl;
	std::cout << "\tIntegers within the interval [6, 16]" << std::endl;	
	std::cout << std::endl;
	
	error_messages_stack.append("\n\tRefinement exponent out of range");
      };

// Error message in case of missing or extra arguments       
    auto valid_usage = [&] () -> void 
      { 
      std::cout << "###\tERROR:\t" << "Correct usage is:" << std::endl;
      std::cout << "\t./run_computation \"Initial Condition\" \"Method\" \"Refinement Exponent\"" << std::endl;   
      
      std::cout << std::endl;
      
          throw std::out_of_range("\n\tIncorrect input format");
      };

// Map containing validated initial condition name, flux name, and grid refinement exponent      
  std::map<std::string, std::string> arguments;

// Case of missing or extra arguments  
  if (argc != 4)
    valid_usage();

// Case of invalid initial condition input string 
  if (initial_conditions.find(argv[1]) == initial_conditions.end())  
    valid_initial_conditions();
  else  
    arguments.emplace("initial_condition", argv[1]); 
  
// Case of invalid flux input string  
  if (fluxes.find(argv[2]) == fluxes.end())
    valid_fluxes();
  else  
    arguments.emplace("flux", argv[2]); 

// Case of invalid refinement exponent input  
  if (((std::stoi(argv[3]) < 6) || (std::stoi(argv[3]) > 16))) 
     valid_refinement_exponent_range();
  else
    arguments.emplace("refinement_exponent", argv[3]); 

// If any errors occured, throw an exception and print the list of occured errors  
  if (!error_messages_stack.empty())
    throw std::out_of_range(error_messages_stack);
    
  return arguments;
  
};

  /** 
   * \brief Acquire the initial data and computational attributes from the database, 
   * execute the requested computation and output the results to the database.
   * 
   * The choice of fluxes based on a particular task is made through template instantiation.
   * The flux must be defined as a function object, i.e. a class containing data relevant to the 
   * computation of that particular class, and also a method operator() (often called "function call," 
   * "call," or "application" operator) which actually computes the fluxes and stores them in the 
   * designated array.
   * 
   * The other solution to the problem of computing the correct flux according to the task input
   * would be to use pointers to functions. This approach, however, has several disadvanatges 
   * compared to using function objects. For example, it is simple to define private constants
   * in a class that have meaning only for computing the fluxes but not in any other contexts,
   * and such constants can be initialized by a constructor defined in the class. An example of 
   * this is using a constant (1-CFL)/4, which can be given a name and initialized in the class
   * of Fromm's method. This approach allows to define custom private data based on the public
   * data (e.g. private (1-CFL)/4 based on the public CFL number) and it is required to initialize
   * it only once, when the class object is initialized.
   * 
   * In case of a pointer to function, however, one would have to initialize such constant variables
   * every time the function is called, and there is no language capability to separate the private
   * variables of the function apart from the other operations taking place in the function, which
   * is more error prone.
   * 
   * \param arguments Map containing valid initial condition name, flux name, and grid refinement exponent
   */

template<typename Flux>
void main_loop(std::map<std::string, std::string> & arguments) {
 
// Path to the initial data group where the initial data dataset is stored  
  std::string input_group_path = "/" + arguments["initial_condition"];
// Path to the data group where the dataset containing results of the computation will be stored  
  std::string group_path = input_group_path + "/" + arguments["flux"]; 
// Path to the dataset where the results of the computation will be stored    
  std::string dataset_name = "k = " + arguments["refinement_exponent"];
// Path to the initial data dataset containing the initial data 
  std::string dataset_initial_data = dataset_name + " initial_data";

// Open the computational database  
  H5::H5File computations_output_file("output_database/computations_output.hdf5", H5F_ACC_RDWR);

// Data group containing the initial data  
  H5::Group input_group = computations_output_file.openGroup(input_group_path);
// Data group where the dataset containing results of the computation will be stored  
  H5::Group output_group = computations_output_file.openGroup(group_path);

// Grid refinement exponent  
  const unsigned int refinement_exponent = std::stoi(arguments["refinement_exponent"]);

// Retreive the attributes pertaining to the requested computation, e.g. CFL number
  
// Open attributes as objects  
  H5::Attribute S_attribute = output_group.openAttribute("CFL");
  H5::Attribute a_attribute = output_group.openAttribute("a");
  H5::Attribute T_attribute = output_group.openAttribute("T");

// Variables to store the attributes temporarily (buffer)  
  double _S;
  double _a;
  double _T;

// Write attributes to the temporary variables  
  S_attribute.read(H5::PredType::NATIVE_DOUBLE, &_S);
  a_attribute.read(H5::PredType::NATIVE_DOUBLE, &_a);
  T_attribute.read(H5::PredType::NATIVE_DOUBLE, &_T);

// Store attributes as constants  
  const double S = _S;
  const double a = _a;
  const double T = _T;

// Number of cells in the discretization of the grid  
  const unsigned int M = pow(2, refinement_exponent);
// t\h, a constant used in the conservative finite-difference update of the scalar field  
  const double t_over_h = S/a;
// Number of timesteps required to compute the solution with the given parameters: N = T/t = T*M/(t/h)  
  const unsigned int N = std::floor(T*M/t_over_h+0.5);

// Allocate memory for the scalar field plus four ghost cells  
std::valarray<double> _field(M+4);
// Allocate memory for fluxes at the edges of the cells: there are M+1 edges for M cells
std::valarray<double> _fluxes(M+1);
// Shift the pointer to the beginning of the array to allow for indeces in the range [-2, M+1]
  double * field = &_field[0]+2;
// Pointer to the beginning of the array of fluxes: needed in order to use the flux data globally  
  double * fluxes = &_fluxes[0];

// Initialize the function object used to compute the fluxes with the public data   
  Flux update_flux {M, S, a, fluxes, field};
/*
 * Array of dimensions of the datasets stored in the database: in this case, it is one dimensional array of length M. 
 * In general, the i-th component is the length of the i-th dimension of the array that one wishes to store in a HDF5 file.  
*/
  hsize_t dimensions[1] = {M};
/* 
 * Allocate space in the database. 
 * First parameter - rank of the dataset (i.e. number of dimensions)
 * Second: array of dimensions explained above
*/
  H5::DataSpace dataspace(1, dimensions); 

// Dataset object corresponding to the dataset in which the scalar field will be stored  
  H5::DataSet field_dataset = output_group.createDataSet(dataset_name, H5::PredType::NATIVE_DOUBLE, dataspace);
// Dataset object corresponding to the initial data dataset from which initial data will be retrieved   
  H5::DataSet field_initial_dataset = input_group.openDataSet(dataset_initial_data); 

// Retrieval of the initial data and storing in in the scalar field array (i.e. initializing the field array with the initial data)   
  field_initial_dataset.read(field, H5::PredType::NATIVE_DOUBLE);

// Indicate that computation has started to the user  
std::cout << group_path << "/" << dataset_name << ": computation in progress" << std::endl;  
// Mark the time of the beginning of the computation
auto t_0 = std::chrono::system_clock::now();  

// Main computational loop: iterate over all timestep
for (unsigned int t = 0; t < N; ++t) {
  
// Periodic boundary conditions for the cells: update the ghost cells     
  field[-1] = field[M-1];
  field[-2] = field[M-2];
  field[M] = field[0];
  field[M+1] = field[1];

// Compute fluxes across the cells 
  update_flux();

// Conservative finite-difference update of the scalar field (e.g. temperature field field)  
  for(unsigned int i = 0; i < M; ++i) 
    {    
/* 
 * Flux indexing convention: 
 *   flux[i] is flux INTO the i-th cell
 *   -flux[i+1] is flux OUT OF the i-th cell <--> flux[i+1] is flux INTO (i+1)-th cell
 * The flux balance for i-th cell is thus
 *   flux[i] - flux[i+1]
*/
      field[i] += t_over_h*(fluxes[i] - fluxes[i+1]);
    }
}

// Mark the time when computation has been completed
auto t_1 = std::chrono::system_clock::now();

// Compute to time of the computation in seconds
auto execution_time_seconds = std::chrono::duration_cast<std::chrono::seconds>(t_1-t_0).count();
// Compute to time of the computation in minutes
auto execution_time_minutes = std::chrono::duration_cast<std::chrono::minutes>(t_1-t_0).count();

// Write the results of the computations (final updated state of the scalar field) into the database
  field_dataset.write(field, H5::PredType::NATIVE_DOUBLE);

// Indicate that the computation has completed and the time it took to the user   
std::cout << group_path << "/" << dataset_name << ": computation completed in " << execution_time_seconds << " seconds (" << execution_time_minutes << " minutes)" << std::endl;

};