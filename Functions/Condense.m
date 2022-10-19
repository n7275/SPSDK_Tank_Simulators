function q = Condense(dt,ii)
  global vapor_mass;
  global mass;
  global Q;

  vapenth_temp = getVapenth(ii);

  if vapor_mass(ii) < dt
		dt = vapor_mass(ii);
  endif

	vapor_mass(ii) -= dt;
	Q += vapenth_temp * dt;

	q = vapenth_temp * dt;
endfunction