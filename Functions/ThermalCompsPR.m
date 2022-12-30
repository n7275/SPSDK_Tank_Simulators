function ThermalCompsPR(dt)

  global mass;
  global MAX_SUB;
  global SPECIFICC_GAS;
  global SPECIFICC_LIQ;
  global MMASS;
  global L_DENSITY;
  global BULK_MOD;
  global VAPPRESS;
  global VAPGRAD;
  global R_CONST;
  global total_mass;
  global vapor_mass;
  global Q;
  global Temp;
  global Volume;
  global Press;
  global ANTIONE_A;
  global ANTIONE_B;
  global PR_Ar;
  global PR_Br;
  global CRITICAL_T;

  ##1. compute average temp

  AvgC = 0;
  Press = 0;
  TotalVapPress = 0;
  n = 0;
  Tcr_equiv = 0;
	

  for ii = 1:MAX_SUB
    AvgC += ((vapor_mass(ii) * SPECIFICC_GAS(ii)) + ((mass(ii) - vapor_mass(ii)) * SPECIFICC_LIQ(ii)));
    n += mass(ii)/MMASS(ii);
    Tcr_equiv += CRITICAL_T(ii) * mass(ii)/sum(mass);
  endfor
    Vm = Volume*1000/n;
    
  if GetMass() > 0;
    AvgC = AvgC / total_mass;
    Temp = Q / AvgC / total_mass;
  else
    Temp = 0;
  endif
  
  for ii = 1:MAX_SUB
    if mass(ii)> 0
      TotalVapPress += GetVapPress(ii) * (mass(ii)/MMASS(ii))/n;
    endif
  endfor

  [A, B, aii, bii] = GetPRCoeffAll(TotalVapPress);
  [x1, x2, x3, ThreeRoots] = SolveCubic(-(1-B),(A-3*B^2-2*B),(-A*B+B^2+B^3))
  
  aii *= 1E6;
  bii *= 1E6;
  
  if(ThreeRoots && Temp < Tcr_equiv)
    #Mixture is either Liquid, VLE, or Vapor
    Zl = x3+(1-B)/3;
    Zv = x1+(1-B)/3;
    
    Vv = Zv*R_CONST*Temp/TotalVapPress*1000;
    Vl = Zl*R_CONST*Temp/TotalVapPress*1000;

    if(Vm > Vv || Vm <Vl)
      Press = (R_CONST/1000*Temp/(Vm-bii) - aii/((Vm*(Vm+bii)+bii*(Vm-bii))))*1E6;
    else
      Press = TotalVapPress;
    endif
  else
    #Mixture is supercritical
    Press = (R_CONST/1000*Temp/(Vm-bii) - aii/((Vm*(Vm+bii)+bii*(Vm-bii))))*1E6;
  endif
  TotalVapPress
  Tcr_equiv
  n
  
##  for ii = 1:1
##    [A, B] = GetPRCoeff(ii);
##    [x1, x2, x3] = SolveCubic(-(1-B),(A-3*B^2-2*B),(-A*B+B^2+B^3));
##    (x1*R_CONST*Temp)/GetVapPress(1)
##  endfor
    
endfunction
