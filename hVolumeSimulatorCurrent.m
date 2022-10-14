#NASSP SPSDK Tank simulator;
tic;

clear all;
clc;
close all;

n_itterations =5000;
itteration = 0;
global delta_time = 10;
global mass;
mass = [0,0,0,0,0,0,0,0,0];#g;
global Q;
Q = 190000.0;


global Volume = 800;
global Press = 0;
global Temp = 90;
global total_mass = 0;
global vapor_mass;
vapor_mass = [0,0,0,0,0,0,0,0,0];#g;

global MAX_SUB = 9;

global MMASS;
MMASS	= [31.998,	2.01588,	18.01528,	28.0134,	44.01,		33.434432,		92.146,			92.01,			4.00260		];		#g/mol
global SPECIFICC_GAS;
SPECIFICC_GAS	= [0.658,		10.183,		1.4108,		0.743,		0.6553,		3.625769,		0.48102,		4.6,			3.12		];		#J/g-K .. assume constant
global SPECIFICC_LIQ;
SPECIFICC_LIQ	=	[1.1519,	9.668,		4.184,		1.7848,		0.858,		3.691041,		1.0724,			1.5525,			5.193		];		#J/g-K .. assume constant
global L_DENSITY;
L_DENSITY	= [1141.0,	70.0,		1000.0,		807.0,		1014.0,		1038.5,			899.0,			1450.0,			0.164		];		#g/L @ 103kPa ..assume constant wrt. temp
global BULK_MOD;
BULK_MOD	= [32e6,		24e6,		2.18e6,		32e6,		32e6,		2.55e6,			1.47397e6,		1.362e6,		10e6		];		#Pa .. assume constant and converted from m^3 to L
global CRITICAL_P;
CRITICAL_P = [350115.0,	89631.0,	1523741.0,	234421.0,	508833.0,	3097574.75,		11692906.154,	10132500.0,		226968.0224 ];		#Pa.. critical pressure
global CRITICAL_T;
CRITICAL_T = [154.7,		33.2,		647.3,		126.2,		304.4,		256.9525,		607.15,			431.15,			5.19		];		#K.. critical temperature
global ANTIONE_A;
ANTIONE_A = [9.3199656,	6.59723,	12.490095,	9.0020008,	12.0892,	8.32957,		13.7222,		14.47645,		4.41952 ];			#Antione Equation A constant gives results in bar, must be converter to Pa	[1]
global ANTIONE_B;
ANTIONE_B = [838.91339,	133.793,	4658.1375,	694.78356,	2353.762,	3158.1575,		5309.7973,		4256.07694,		18.65037 ];			#Antione Equation B constant gives results in bar, must be converter to Pa	[2]
global ACENTRIC;
ACENTRIC = [0.022,		-0.216,		0.345,		0.040,		0.288,		0.416,			0.316,			0.0141345,			-0.390];		#[3] Acentric factor

global R_CONST = 8314.4621	##(L*Pa)/(mol*K)

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

function setTemp(T,ii)
  global SPECIFICC_LIQ;
  global SPECIFICC_GAS;
  global mass;
  global vapor_mass;
  global Q;
  global Temp;
  
  Temp = T;
  Q = Temp * (vapor_mass(ii) * SPECIFICC_GAS(ii) + (mass(ii) - vapor_mass(ii)) * SPECIFICC_LIQ(ii));
endfunction

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
    if(ii == 1)
      density += 0.56 * Temp * Temp - 134.0 * Temp + 6900.0;
    elseif(ii == 2)
      density += 0.03333 * Temp * Temp - 4.3333 * Temp + 73.3333;
    endif

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

Press_array = zeros(1,n_itterations);
Temp_array = zeros(1,n_itterations);
time_array = zeros(1,n_itterations);
vapor_fraction_array = zeros(1,n_itterations);
Q_array = zeros(1,n_itterations);

time = 0;
mass = [0,0,0,5000,0,0,0,0,0];

Temps = [80,90,100,110,120,130,140,157,200,300,500];
incrMass = [0,0,0,5000,0,0,0,0,0];
nMassIncr = 150;

P = zeros(nMassIncr,length(Temps));
V = zeros(nMassIncr,length(Temps));
T = zeros(nMassIncr,length(Temps));


for nn = 1:nMassIncr
  for mm = 1:length(Temps)

    setTemp(Temps(mm),4)
    itteration = 1;
    while(itteration < n_itterations)
      itteration = itteration+1;
      ThermalComps(delta_time)

      Press_array(itteration) = Press*0.000145038;
      Temp_array(itteration) = Temp;
      vapor_fraction_array(itteration) = (sum(vapor_mass))/sum(mass);
      Q_array(itteration) = Q/sum(mass);
      time_array(itteration) = time;
      time += delta_time;
      setTemp(Temps(mm),4);
    endwhile

    ##AX = plotyy(time_array,Press_array,time_array,Temp_array);
    ##grid on;
    ##grid minor;
    ##xlabel("Time [sec]");
    ##ylabel(AX(1),"Pressure [PSI]");
    ##ylabel(AX(2),"Temperature [K]");
    ##figure(2)
    ##waitfor(AX)
    ##AX2 = plotyy(time_array,vapor_fraction_array,time_array,Q_array);
    ##grid on;
    ##grid minor;
    ##xlabel("Time [sec]");
    ##ylabel(AX2(1),"Vapor Fraction");
    ##ylabel(AX2(2),"Specific Energy [J/g]");
    ##waitfor(AX2)

    P(nn,mm) = Press*0.000145038;
##    Press*0.000145038;
    V(nn,mm) = Volume/mass(1);
    T(nn,mm) = Temps(mm);

  endfor
  mass = mass + incrMass;
  nn/nMassIncr*100
  toc
endfor

loglog(V,P);
legend('80K','90K','100K','110K','120K','130K','140K','157K','200K','300K','500K');
grid on;
grid minor;
axis([0.001 0.2 1 10000]);
xlabel("Specific Volume [l/g]");
ylabel("Pressure [PSI]");
toc
