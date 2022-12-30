function setTemp(T,ii)
  global SPECIFICC_LIQ;
  global SPECIFICC_GAS;
  global mass;
  global vapor_mass;
  global Q;
  global Temp;
  global MAX_SUB;
  global total_mass;

  Temp = T;
  AvgC = 0;
  
  for ii = 1:MAX_SUB
    AvgC += ((vapor_mass(ii) * SPECIFICC_GAS(ii)) + ((mass(ii) - vapor_mass(ii)) * SPECIFICC_LIQ(ii)));
  endfor
  
  if GetMass() > 0;
    AvgC = AvgC / total_mass;
  else
    Temp = 0;
  endif
  
  Q = Temp * AvgC * total_mass;
endfunction