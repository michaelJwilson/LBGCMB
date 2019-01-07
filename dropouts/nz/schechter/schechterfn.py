import  numpy  as  np


## See https://pythonhosted.org/Astropysics/coremods/models.html

def SchechterLfn(L, phi_star, L_star, alpha):
  '''                                                                                                                                                       
  Schechter function in luminosity.                                                                                                                         
  '''
  LL = L/L_star

  return (phi_star / L_star) * (LL**alpha) * np.exp(-LL)              ## Phi(L) wrt dL;  [h_70/Mpc]^3 per mag for HSC.                                         

def SchechterMfn(M, phi_star, M_star, alpha):
  '''                                                                                                                                                       
  Schechter funtion in absolute magnitude.                                                                                                                  
  '''
  LL  = 10.**(-0.4*(M - M_star))                                      ## LL = L/L_star                                                                        

  return  np.log(10.)*phi_star*np.exp(-LL)*LL**(1. + alpha)/2.5       ## Phi(M_UV) dM_UV, [Phi_star] per mag;


if __name__ == "__main__":
  print("\n\nDone.\n\n")
  
