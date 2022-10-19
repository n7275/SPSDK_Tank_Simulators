function [A, B] = GetPRCoeff()
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
		Pri = GetVapPress(ii)/CRITICAL_P(ii);
		Tri = Temp/CRITICAL_T(ii);
		
		B += (0.07780*Pri/Tri)*mass(ii)/sum(mass);
		xi = mass(ii)/sum(mass);
		m = 0.37464 + 1.54226 * ACENTRIC(ii) - 0.26992*ACENTRIC(ii)^2;
		ai = 0.45724*(R*CRITICAL_T(ii))^2/CRITICAL_P(ii)*(1+m*(1-sqrt(Tri)))^2;
		Ai = 0.45724*ai*Pri/Tri^2;
		for jj = 1:MAX_SUB;
			Prj = GetVapPress(jj)/CRITICAL_P(jj);
			Trj = Temp/CRITICAL_T(jj);
			aj = 0.45724*(R*CRITICAL_T(jj))^2/CRITICAL_P(jj)*(1+m*(1-sqrt(Trj)))^2;
			Aj = 0.45724*aj*Prj/Trj^2;
			Aij = sqrt(Ai*Aj);
			xj = mass(jj)/sum(mass);
			
			A += xi*xj*Aij
		endfor
	endfor
  return;
endfunction
