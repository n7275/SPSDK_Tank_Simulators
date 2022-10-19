function m = GetMass()

  global MAX_SUB;
  global mass = [9];
  global total_mass;

  mass_temporary = 0;
  for ii=1:MAX_SUB
    mass_temporary += mass(ii);
  endfor
  total_mass = mass_temporary;
  m = mass_temporary;
endfunction