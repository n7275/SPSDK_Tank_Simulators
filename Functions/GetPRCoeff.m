function [A, B] = GetPRCoeff(ii)
	global ACENTRIC;
	global CRITICAL_P;
	global CRITICAL_T;
	global MAX_SUB;
	global R_CONST
  global Temp;
  global mass;
  global Press

  R = R_CONST/1000;
  Pr = Press/CRITICAL_P(ii);
  Tr = Temp/CRITICAL_T(ii);

	kappa = 0.37464 + 1.54226 * ACENTRIC(ii) - 0.26992*ACENTRIC(ii)^2;
  alpha = (1+kappa*(1-sqrt(Tr)))^2;
  
  a = 0.4572355289*(CRITICAL_T(ii)*R)^2*alpha/CRITICAL_P(ii);
  b = 0.0777960739*R*CRITICAL_T(ii)/CRITICAL_P(ii);
  A = a*Press/(R*Temp)^2
  B = b*Press/R/Temp
  

endfunction
