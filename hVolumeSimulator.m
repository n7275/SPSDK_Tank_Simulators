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
global SPECIFICC;
SPECIFICC	= [1.669,		9.668,		4.184,		1.040,		0.858,		  3.691041,		2.9056392,		1.270,			5.193		];		#J/g-K .. assume constant
global VAPENTH;
VAPENTH	= [213.13,	445.46,		2260.0,		198.83,		347.0,		1769.195,		991.01556,		414.3,			0.0829		];		#J/g
global VAPPRESS;
VAPPRESS	= [1314841.0,	4925221.0,	39441.0,	1528361.0,	493284.0,	25639.45,		21722.212986,	206782.99342,	14778377.09 ];		#Pa @ 273.00K
global VAPGRAD;
VAPGRAD	= [6556.0,	19045.0,	680.0,		7228.0,		4800.0,		52.87,			111.1,			1754.255683,	874.9447005 ];		#Pa/K.. assume linear dependence of PV / K
global L_DENSITY;
L_DENSITY	= [1141.0,	70.0,		1000.0,		807.0,		1014.0,		1038.5,			899.0,			1450.0,			0.164		];		#g/L @ 103kPa ..assume constant wrt. temp
global BULK_MOD;
BULK_MOD	= [32e6,		24e6,		2.18e6,		32e6,		32e6,		2.55e6,			1.47397e6,		1.362e6,		10e6		];		#Pa .. assume constant and converted from m^3 to L
global CRITICAL_P;
CRITICAL_P = [350115.0,	89631.0,	1523741.0,	234421.0,	508833.0,	3097574.75,		11692906.154,	10132500.0,		226968.0224 ];		#Pa.. critical pressure
global CRITICAL_T;
CRITICAL_T = [154.7,		33.2,		647.3,		126.2,		304.4,		256.9525,		607.15,			431.15,			5.19		];		#K.. critical temperature
##VdW_A		= [1.382E-5,	0.2453E-5,	5.537E-5,	1.37E-5,	3.658E-5,	17.2785E-5,		11.575E-5, 		0,				0.0346E-5	];		//Van der Waals Coefficient Pa*L^2/Mol^2
##VdW_B		= [0.03186,	0.02651,	0.03049,	0.0387,		0.04286,	0.091195,		0.07315,		0,				0.0238		];		//Van der Waals Coefficient B, L/mol
##HvapA		= [0.03186,	0.02651,	0.03049,	0.0387,		0.04286,	0,				0,				0,				0.0238		];

global R_CONST = 8314.4621	##(L*Pa)/(mol*K)

function q = Boil(dt,ii)
  global vapor_mass;
  global mass;
  global VAPENTH;
  global Q;

  if vapor_mass(ii) + dt > mass(ii) - 1.0
		dt = mass(ii) - 1.0 - vapor_mass(ii);
  endif

	if dt < 0
    q = 0;
		return;
  endif


	if Q < VAPENTH(ii) * dt
		dt = Q / VAPENTH(ii);
  endif

	vapor_mass(ii) += dt;
	Q -= VAPENTH(ii) * dt;
	q = -VAPENTH(ii) * dt;
endfunction

function q = Condense(dt,ii)
  global vapor_mass;
  global mass;
  global VAPENTH;
  global Q;

  if vapor_mass(ii) < dt
		dt = vapor_mass(ii);
  endif

	vapor_mass(ii) -= dt;
	Q += VAPENTH(ii) * dt;

	q = VAPENTH(ii) * dt;
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
  global SPECIFICC;
  global mass;
  global Q;

  Q = mass(ii) * T * SPECIFICC(ii);
endfunction

function ThermalComps(dt)

  global mass;
  global MAX_SUB;
  global SPECIFICC;
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

  ##1. compute average temp

  AvgC = 0;
	vap_press = 0;

  for ii = 1:MAX_SUB
    AvgC += mass(ii) * SPECIFICC(ii);
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
    vap_press = VAPPRESS(ii) - (273.0 - Temp) * VAPGRAD(ii);  ##this is vapor pressure of current substance

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
