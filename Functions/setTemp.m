function setTemp(T,ii)
  global SPECIFICC_LIQ;
  global SPECIFICC_GAS;
  global mass;
  global vapor_mass;
  global Q;
  global Temp;

  Temp = T;
  Q = Temp * (vapor_mass(ii) * SPECIFICC_GAS(ii) + (mass(ii) - vapor_mass(ii)) * SPECIFICC_LIQ(ii));
endfunction