function hV = getVapenth(ii)
  global Temp;
  global CRITICAL_T;
  global ACENTRIC;
  global MMASS;
  global R_CONST;

  if(Temp > CRITICAL_T(ii) || Temp <= 0.0)
    hV = 0.0;
    return;
  else
    hV = (R_CONST/1000*CRITICAL_T(ii)*...
		(7.08*((1-Temp/CRITICAL_T(ii))^0.354) +...
			10.95 * ACENTRIC(ii) * ((1 - Temp / CRITICAL_T(ii))^ 0.456)))/ MMASS(ii);
  endif
endfunction