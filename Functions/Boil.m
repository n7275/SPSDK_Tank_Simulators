function q = Boil(dt,ii)
  global vapor_mass;
  global mass;
  global Q;

  vapenth_temp = getVapenth(ii);

  if vapor_mass(ii) + dt > mass(ii) - 1.0
		dt = mass(ii) - 1.0 - vapor_mass(ii);
  endif

	if dt < 0
    q = 0;
		return;
  endif


	if Q < vapenth_temp * dt
		dt = Q / vapenth_temp;
  endif

	vapor_mass(ii) += dt;
	Q -= vapenth_temp * dt;
	q = -vapenth_temp * dt;
endfunction