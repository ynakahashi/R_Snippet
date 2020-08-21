subroutine get_int_parms(sml,eps,big,mnlam,rsqmax,pmin,exmx,itrace)
implicit double precision(a-h,o-z)                                
data sml0,eps0,big0,mnlam0,rsqmax0,pmin0,exmx0,itrace0  /1.0d-5,1.0d-6,9.9d35,5,0.999,1.0d-9,250.0,0/  ! big0 に9.9d35 は 9.9e+35 と同じ？　多分、倍精度実数
sml=sml0                                                          
eps=eps0                                                          
big=big0                                                          
mnlam=mnlam0                                                      
rsqmax=rsqmax0                                                    
pmin=pmin0                                                        
exmx=exmx0                                                        
itrace=itrace0                                                    
return                                                            
entry chg_fract_dev(arg)                                          
sml0=arg                                                          
return                                                            
entry chg_dev_max(arg)                                            
rsqmax0=arg                                                       
return                                                            
entry chg_min_flmin(arg)                                          
eps0=arg                                                          
return                                                            
entry chg_big(arg)                                                
big0=arg                                                          
return                                                            
entry chg_min_lambdas(irg)                                        
mnlam0=irg                                                        
return                                                            
entry chg_min_null_prob(arg)                                      
pmin0=arg                                                         
return                                                            
entry chg_max_exp(arg)                                            
exmx0=arg                                                         
return                                                            
entry chg_itrace(irg)                                             
itrace0=irg                                                       
return                                                            
end