#include <cmath>
#include <set>
#include <algorithm>
#include <valarray>

#include "boost/math/special_functions/sign.hpp"
  
/* NOTE:
 * Flux indexing convention: 
 *   flux[i] is flux INTO the i-th cell
 *   -flux[i+1] is flux OUT OF the i-th cell <--> flux[i+1] is flux INTO (i+1)-th cell
 * The flux balance for i-th cell is thus
 *   flux[i] - flux[i+1]
 * Thus, the boundary condition for fluxes is:
 *   flux[0] = flux[M]
 * since 0-th and M-th cells are identical for periodic boundary conditions
*/

/*
 * NOTE:    
 * Statements involving negative indeces should be computed outside the loop 
 * to avoid segmentation fault, which probably occurs due to the exclusion
 * of memory corresponding to negative indeces from the memory page in 
 * which the operating system allows this process to operate. 
 * That probably happens due to the way the compiler arranges the loop.
*/
 
// Abstract class which serves as a blueprint for other classes of fluxes 
class Flux_base {
public:
// Default constructor which initializes the private data with the public data
  Flux_base(unsigned int M = 0, 
	    double CFL = 0.9, 
	    double a = 3.0,  
	    double * _fluxes = nullptr, 
	    double * _field = nullptr
	   ) : 
	   M(M), 
	   CFL(CFL), 
	   a(a), 
	   _fluxes(_fluxes), 
	   _field(_field) 
	   {};
/*	   
 * operator() is a virtual method, since the classes inheriting from this base class
 * will redefine it in order to compute their particular fluxes
*/
  virtual void operator()() const = 0;
// Empty destructor is sufficient, since no memory will be explicitly allocated on free store by these classes  
  virtual ~Flux_base() {};
// Protected means "will be inherited by children classes"
protected:  
 const unsigned int M; 
 const double CFL;
 const double a;
 double *const _fluxes;
 double *const _field;  
};

class Upwind : public Flux_base {
public:
  Upwind(unsigned int M = 0, 
	 double CFL = 0.9, 
	 double a = 3.0, 
	 double * _fluxes = nullptr, 
	 double * _field = nullptr
	) : 
	Flux_base(M, CFL, a, _fluxes, _field) 
	{};  
	
  ~Upwind() {};  
  
  void operator()() const 
  {
      _fluxes[0] = a*_field[-1];
      for(unsigned int i = 1; i < M; ++i) 
	{
	  _fluxes[i] = a*_field[i-1];
	}
// Periodic boundary condition for fluxes (cells 0 and M are identical) 
      _fluxes[M] = _fluxes[0];
  };
};

class Lax_Friedrichs : public Flux_base {
public:
  Lax_Friedrichs(unsigned int M = 0, 
	 double CFL = 0.9, 
	 double a = 3.0, 
	 double * _fluxes = nullptr, 
	 double * _field = nullptr
	) : 
	Flux_base(M, CFL, a, _fluxes, _field), 
	inv_CFL(1/CFL), 
	half_a(a/2) 
	{};  
	
  ~Lax_Friedrichs() {};  
  
  void operator()() const 
  {
      _fluxes[0] = half_a*( (_field[-1] + _field[0]) + inv_CFL*(_field[-1] - _field[0]) );
      for(unsigned int i = 1; i < M; ++i) 
      {    
       _fluxes[i] = half_a*( (_field[i-1] + _field[i]) + inv_CFL*(_field[i-1] - _field[i]) );
      }
// Periodic boundary condition for fluxes (cells 0 and M are identical) 
      _fluxes[M] = _fluxes[0];
  };
    
private:
  const double inv_CFL;
  const double half_a;
};

class Lax_Wendroff : public Flux_base {
public:
  Lax_Wendroff(unsigned int M = 0, 
	 double CFL = 0.9, 
	 double a = 3.0, 
	 double * _fluxes = nullptr, 
	 double * _field = nullptr
	) : 
	Flux_base(M, CFL, a, _fluxes, _field), 
	half_a(a/2) 
	{};
	
  ~Lax_Wendroff() {};  
  
  void operator()() const 
  {  
       _fluxes[0] = half_a*( (_field[0] + _field[-1]) + CFL*(_field[-1] - _field[0]) );
      for(unsigned int i = 1; i < M; ++i) 
      {    
       _fluxes[i] = half_a*( (_field[i] + _field[i-1]) + CFL*(_field[i-1] - _field[i]) );
      }
// Periodic boundary condition for fluxes (cells 0 and M are identical) 
       _fluxes[M] = _fluxes[0];
  };
    
private:
  const double half_a;
};

class Fromm : public Flux_base{
public:
  Fromm(unsigned int M = 0, 
	 double CFL = 0.9, 
	 double a = 3.0, 
	 double * _fluxes = nullptr, 
	 double * _field = nullptr
	) : 
	Flux_base(M, CFL, a, _fluxes, _field), 
	CFL_expr((1-CFL)/4) 
	{};
	
  ~Fromm() {};  
  
  void operator()() const 
  {   
       _fluxes[0] = a*( _field[-1] + CFL_expr*(_field[0] - _field[-2]) );
       _fluxes[1] = a*( _field[0] + CFL_expr*(_field[1] - _field[-1]) );
      for(unsigned int i = 2; i < M; ++i) 
      {    
       _fluxes[i] = a*( _field[i-1] + CFL_expr*(_field[i] - _field[i-2]) );
      }      
// Periodic boundary condition for fluxes (cells 0 and M are identical) 
       _fluxes[M] = _fluxes[0];
  };
    
private:
  const double CFL_expr;
};

class Fromm_van_Leer : public Flux_base {
public:
  Fromm_van_Leer(unsigned int M = 0, 
	 double CFL = 0.9, 
	 double a = 3.0, 
	 double * _fluxes = nullptr, 
	 double * _field = nullptr
	) : 
	Flux_base(M, CFL, a, _fluxes, _field), 
	CFL_expr(1-CFL), 
	delta(M+1),
	u1(0),
	u2(0),
	u3(0),
	sign_u3(0)
	{};
	
  ~Fromm_van_Leer() {};  

  void update_limiters() const 
  {
     u1 = _field[-1] - _field[-2];
     u2 = _field[0] - _field[-1];
     u3 = 0.25*(u1 + u2);
     sign_u3 = boost::math::sign<double>(u3);
          
     delta[0] = (u1*u2 > 0) ? 
      sign_u3 * CFL_expr * std::min<double>({std::abs(u1), std::abs(u2), std::abs(u3)}) :
      0;
           
     u1 = _field[0] - _field[-1];
     u2 = _field[1] - _field[0];
     u3 = 0.25*(u1 + u2);
     sign_u3 = boost::math::sign<double>(u3);
          
     delta[1] = (u1*u2 > 0) ? 
      sign_u3 * CFL_expr * std::min<double>({std::abs(u1), std::abs(u2), std::abs(u3)}) :
      0;
  
    for(unsigned int i = 2; i < M; ++i) 
    {  
     u1 = _field[i-1] - _field[i-2];
     u2 = _field[i] - _field[i-1];
     u3 = 0.25*(u1 + u2);
     sign_u3 = boost::math::sign<double>(u3);
          
     delta[i] = (u1*u2 > 0) ? 
      sign_u3 * CFL_expr * std::min<double>({std::abs(u1), std::abs(u2), std::abs(u3)}) :
      0;
    }
    
// Periodic boundary condition for fluxes (cells 0 and M are identical) 
     delta[M] = delta[0];
  }
  
  void operator()() const 
  {        
    update_limiters();
          
    _fluxes[0] = a*( _field[-1] + delta[0] );      
     
    for(unsigned int i = 1; i < M+1; ++i) 
     {
      _fluxes[i] = a*( _field[i-1] + delta[i] );
     }     
      
  };
    
private:
  const double CFL_expr;     
  mutable std::valarray<double> delta;
  mutable double u1;
  mutable double u2;
  mutable double u3;
  mutable int sign_u3;
};

class Flux_Corrected_Transport : public Flux_base {
public:
  Flux_Corrected_Transport(unsigned int M = 0, 
	 double CFL = 0.9, 
	 double a = 3.0, 
	 double * _fluxes = nullptr, 
	 double * _field = nullptr
	) : 
	Flux_base(M, CFL, a, _fluxes, _field), 
	CFL_expr(0.5*a*(1-CFL)), 
	_anti_diffusive_fluxes(M+1),
	S(0),
	A1(0),
	A2(0),
	A3(0),
	theta(0),
	t_over_h(CFL/a),
	h_over_t(1/t_over_h), 
	low_order_fluxes{M, CFL, a, _fluxes, _field} 
	{};  
  
  ~Flux_Corrected_Transport() {};  
    
  void anti_diffusive_fluxes() const 
  {
    _anti_diffusive_fluxes[0] = CFL_expr*( _field[0] - _field[-1] );
   for(unsigned int i = 1; i < M; ++i) 
    {
     _anti_diffusive_fluxes[i] = CFL_expr*( _field[i] - _field[i-1] );
    }    
// Periodic boundary condition for fluxes (cells 0 and M are identical)     
   _anti_diffusive_fluxes[M] = _anti_diffusive_fluxes[0];
  }  
 
  void low_order_estimate() const 
  {
  low_order_fluxes();

  for(unsigned int i = 0; i < M; ++i) 
    {
      _field[i] += t_over_h*(_fluxes[i] - _fluxes[i+1]);
    }       
   
  _field[-1] = _field[M-1];
  _field[-2] = _field[M-2];
  _field[M] = _field[0];
  _field[M+1] = _field[1];
  }
  
  void operator()() const 
  {

    anti_diffusive_fluxes();  
    low_order_estimate();
    
      S = boost::math::sign<double>(_anti_diffusive_fluxes[0]);
      A1 = S * h_over_t * (_field[1] - _field[0]);
      A2 = S * h_over_t * (_field[-1] - _field[-2]); 
      A3 = S * _anti_diffusive_fluxes[0];
      theta = std::min<double>({ A1, A2, A3 });
     _fluxes[0] = S * std::max<double>({0, theta});
    
      S = boost::math::sign<double>(_anti_diffusive_fluxes[1]);
      A1 = S * h_over_t * (_field[2] - _field[1]);
      A2 = S * h_over_t * (_field[0] - _field[-1]); 
      A3 = S * _anti_diffusive_fluxes[1];
      theta = std::min<double>({ A1, A2, A3 });
     _fluxes[0] = S * std::max<double>({0, theta});
     
    for(unsigned int i = 2; i < M-1; ++i) 
     {
      S = boost::math::sign<double>(_anti_diffusive_fluxes[i]);
      A1 = S * h_over_t * (_field[i+1] - _field[i]);
      A2 = S * h_over_t * (_field[i-1] - _field[i-2]); 
      A3 = S * _anti_diffusive_fluxes[i]; // abs(A) = sign(A)*A
      theta = std::min<double>({ A1, A2, A3 });
     _fluxes[i] = S * std::max<double>({0, theta});
     }
     
      S = boost::math::sign<double>(_anti_diffusive_fluxes[M-1]);
      A1 = S * h_over_t * (_field[M] - _field[M-1]);
      A2 = S * h_over_t * (_field[M-2] - _field[M-3]); 
      A3 = S * _anti_diffusive_fluxes[M-1];
      theta = std::min<double>({ A1, A2, A3 });
     _fluxes[M-1] = S * std::max<double>({0, theta});
     
// Periodic boundary condition for fluxes (cells 0 and M are identical)       
     _fluxes[M] = _fluxes[0]; 
  };
    
private:
  const double CFL_expr;
  mutable std::valarray<double> _anti_diffusive_fluxes;
  mutable double S;
  mutable double A1;
  mutable double A2;
  mutable double A3;
  mutable double theta;
  const double t_over_h;
  const double h_over_t;   
  Upwind low_order_fluxes;
};

class Lax_Wendroff_Fourth_Order : public Flux_base {
public:
  Lax_Wendroff_Fourth_Order(unsigned int M = 0, 
	 double CFL = 0.9, 
	 double a = 3.0, 
	 double * _fluxes = nullptr, 
	 double * _field = nullptr
	) : 
	Flux_base(M, CFL, a, _fluxes, _field), 
	half_CFL(CFL/2), 
	CFL_expr(a*(4*pow(CFL, 3)+1)/16),
	alpha(a*7/12),
	beta(a*1/12),
	gamma(a*5/4),
	u_n(0),
	F(0),
	D(0) 
	{};  
	
  ~Lax_Wendroff_Fourth_Order() {};  
  
  void operator()() const 
  {
// NOTE: The convergence will be 3-rd order if _fluxes[i] = D;

  u_n = alpha*(_field[0] + _field[-1]) - beta*(_field[1] + _field[-2]);
  F = u_n - half_CFL * (gamma * (_field[0] - _field[-1]) - beta * (_field[1] - _field[-2]) );
  D = CFL_expr * ( ( _field[1] - _field[-2] ) - 3*( _field[0] - _field[-1]) );
  
  _fluxes[0] = F + D;

  u_n = alpha*(_field[1] + _field[0]) - beta*(_field[2] + _field[-1]);
  F = u_n - half_CFL * (gamma * (_field[1] - _field[0]) - beta * (_field[2] - _field[-1]) );
  D = CFL_expr * ( ( _field[2] - _field[-1] ) - 3*( _field[1] - _field[0]) );
  
  _fluxes[1] = F + D;
  
  for(unsigned int i = 2; i < M; ++i) 
    {    
    u_n = alpha*(_field[i] + _field[i-1]) - beta*(_field[i+1] + _field[i-2]);
    F = u_n - half_CFL * (gamma * (_field[i] - _field[i-1]) - beta * (_field[i+1] - _field[i-2]) );
    D = CFL_expr * ( ( _field[i+1] - _field[i-2] ) - 3*( _field[i] - _field[i-1]) );
    
    _fluxes[i] = F + D;
    }
  
// Periodic boundary condition for fluxes (cells 0 and M are identical) 
  _fluxes[M] = _fluxes[0];
  
  }; 
    
private:
  const double half_CFL;
  const double CFL_expr;
  const double alpha;
  const double beta;
  const double gamma;
  mutable double u_n;
  mutable double F;
  mutable double D;
};