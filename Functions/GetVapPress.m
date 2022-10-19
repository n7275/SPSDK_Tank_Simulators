function vap_press = GetVapPress(ii)
  global mass;
  global ANTIONE_A;
  global ANTIONE_B;
  global Temp;
  
  vap_press = mass(ii)/sum(mass)*exp(ANTIONE_A(ii) - (ANTIONE_B(ii) / Temp))*1E5;
endfunction