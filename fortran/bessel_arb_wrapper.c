// use arb C library to provide accurate Bessel
// functions of complex argument and real order
// http://fredrikj.net/arb/acb_hypgeom.html
#include "acb_hypgeom.h"

double _Complex arb_K(double gcc_nu, double _Complex gcc_z, int kode)
{

  // arb complex types
  acb_t acb_z, acb_res, acb_nu;
  
  // built-in double precision complex
  double _Complex gcc_K;
  
  acb_init(acb_res);
  acb_init(acb_nu);
  acb_zero(acb_nu);
  acb_init(acb_z);

  // gcc_z -> acb_z (argument)
  acb_set_d_d(acb_z, __real__(gcc_z), __imag__(gcc_z));

  // gcc_nu -> acb_nu (order)
  acb_set_d(acb_nu, gcc_nu);

  int extra = 30;
  // compute Bessel function 
  if (kode == 1){ // unscaled
    do {
      acb_hypgeom_bessel_k(acb_res, acb_nu, acb_z, 53 + extra);
      extra *= 2;
    } while (acb_rel_accuracy_bits(acb_res) < 53);
  }
  //else // scaled
    //do {
    //  acb_hypgeom_bessel_k_scaled(acb_res, acb_nu, acb_z, 53 + extra);
    //  extra *= 2;
    //} while (acb_rel_accuracy_bits(acb_res) < 53);    
      
  
  //  arb_res -> gcc_K
  __real__(gcc_K) = arf_get_d(arb_midref(acb_realref(acb_res)), ARF_RND_DOWN);
  __imag__(gcc_K) = arf_get_d(arb_midref(acb_imagref(acb_res)), ARF_RND_DOWN);
  
  acb_clear(acb_res);
  acb_clear(acb_nu);
  acb_clear(acb_z);

  flint_cleanup();
  return gcc_K;
}

double _Complex arb_I(double gcc_nu, double _Complex gcc_z, int kode)
{

  // arb complex types
  acb_t acb_z, acb_res, acb_nu;

  // built-in double precision complex
  double _Complex gcc_I;

  acb_init(acb_res);
  acb_init(acb_nu);
  acb_zero(acb_nu);
  acb_init(acb_z);

  // gcc_z -> acb_z (argument)
  acb_set_d_d(acb_z, __real__(gcc_z), __imag__(gcc_z));

  // gcc_nu -> acb_nu (order)
  acb_set_d(acb_nu, gcc_nu);

  int extra = 30;
  // compute bessel function
  if (kode == 1){
    do {
      acb_hypgeom_bessel_i(acb_res, acb_nu, acb_z, 53 + extra);
      extra *= 2;
    } while (acb_rel_accuracy_bits(acb_res) < 53);
  }
  //else
    //do {
    //  acb_hypgeom_bessel_i_scaled(acb_res, acb_nu, acb_z, 53 + extra);
    //  extra *= 2;
    //} while (acb_rel_accuracy_bits(acb_res) < 53);
    
  
  //  arb_res -> gcc_I
  __real__(gcc_I) = arf_get_d(arb_midref(acb_realref(acb_res)), ARF_RND_DOWN);
  __imag__(gcc_I) = arf_get_d(arb_midref(acb_imagref(acb_res)), ARF_RND_DOWN);
  
  acb_clear(acb_res);
  acb_clear(acb_nu);
  acb_clear(acb_z);

  flint_cleanup();
  return gcc_I;
}


