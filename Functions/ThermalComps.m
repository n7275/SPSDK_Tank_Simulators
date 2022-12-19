function ThermalComps(dt)

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

  ##1. compute average temp

  AvgC = 0;
	vap_press = 0;

  for ii = 1:MAX_SUB
    AvgC += ((vapor_mass(ii) * SPECIFICC_GAS(ii)) + ((mass(ii) - vapor_mass(ii)) * SPECIFICC_GAS(ii)));
  endfor

  if GetMass() > 0;
    AvgC = AvgC / total_mass;
    Temp = Q / AvgC / total_mass;
  else
    Temp = 0;
  endif

  ##2. Compute average Press
  m_i = 0;		##mols
	NV = 0;		##litres
	PNV = 0;
	tNV = 0;

  for ii = 1:MAX_SUB
    m_i += vapor_mass(ii) / MMASS(ii);
    density = L_DENSITY(ii);
##    if(ii == 1)
##      density += 0.56 * Temp * Temp - 134.0 * Temp + 6900.0;
##    elseif(ii == 2)
##      density += 0.03333 * Temp * Temp - 4.3333 * Temp + 73.3333;
##    endif

    tNV = (mass(ii) - vapor_mass(ii)) / density;	##Units of L
		NV += tNV;	##Units of L

		PNV += tNV / BULK_MOD(ii);	##Units of L/Pa

  endfor
  m_i = -m_i * R_CONST * Temp; ##Units of L*Pa
  NV = Volume - NV;	##Units of L

  delta = delta = (NV * NV) - (4.0 * m_i * PNV); ##delta of quadric eq. P^2*PNV+ P*NV + m_i = 0

  if PNV > 0
    Press = (-NV + sqrt(delta)) / (2.0 * PNV);
  else
    Press = 0;
  endif

  for ii = 1:MAX_SUB
    if(Temp <= 0.0)
      vap_press = 0.0;
    else
     vap_press = exp(ANTIONE_A(ii) - (ANTIONE_B(ii) / Temp))*1E5; #this is vapor pressure of current substance
    endif

    #need to boil material if vapor pressure > pressure, otherwise condense
    if vap_press > Press
      Q += Boil(dt,ii);
    else
      Q += Condense(dt,ii);
    endif

  endfor

endfunction
