function [A, B] = GetPRCoeffAll()
	global ACENTRIC;
	global CRITICAL_P;
	global CRITICAL_T;
	global MAX_SUB;
	global R_CONST
  global Temp;
  global mass;

  R = R_CONST/1000;


	A = 0;
	B = 0;

	for ii=1:MAX_SUB;
		Pi = GetVapPress(ii);
    Tri = Temp/CRITICAL_T(ii);
		xi = mass(ii)/sum(mass);
		mi = 0.37464 + 1.54226 * ACENTRIC(ii) - 0.26992*ACENTRIC(ii)^2;

		ai = 0.4572355289*(R*CRITICAL_T(ii))^2/CRITICAL_P(ii)*(1+mi*(1-sqrt(Tri)))^2;
		Ai = ai*Pi/(R*Temp)^2;

    B += (0.0777960739*R*CRITICAL_T(ii)/CRITICAL_P(ii))*(GetVapPress(ii)/R/Temp)*mass(ii)/sum(mass);
		for jj = 1:MAX_SUB;
			Pj = GetVapPress(ii);
      mj = 0.37464 + 1.54226 * ACENTRIC(jj) - 0.26992*ACENTRIC(jj)^2;
      Trj = Temp/CRITICAL_T(jj);
			aj = 0.4572355289*(R*CRITICAL_T(jj))^2/CRITICAL_P(jj)*(1+mj*(1-sqrt(Trj)))^2;
			Aj = aj*Pj/(R*Temp)^2;

			Aij = sqrt(Ai*Aj);
			xj = mass(jj)/sum(mass);

			A += xi*xj*Aij;
		endfor
	endfor
  return;
endfunction