#NASSP SPSDK Tank simulator;
tic;

clear all;
clc;
close all;

substance2plot = 1;

n_itterations =10000;
itteration = 0;
global delta_time = 0.1;
global mass;
mass = [30,0,0,0,0,0,0,0,0];#g;
global Q;
Q = 2700.0;


global Volume = 10;
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
CRITICAL_P = [4926063.722233855,	1303085.962148109,	19910602.43884709,	3300253.307024634,	7800213.020695794,	3097566.149152389,		14500299.44342058,	9998268.775903003,		228416.3259285132 ];		#Pa.. critical pressure
global CRITICAL_T;
CRITICAL_T = [154.7,		33.2,		647.3,		126.2,		304.4,		645,		607.15,			431.15,			5.19		];		#K.. critical temperature
global ANTIONE_A;
ANTIONE_A = [9.3199656,	6.59723,	12.490095,	9.0020008,	12.0892,	8.32957,		13.7222,		14.47645,		4.41952 ];			#Antione Equation A constant gives results in bar, must be converter to Pa	[1]
global ANTIONE_B;
ANTIONE_B = [838.91339,	133.793,	4658.1375,	694.78356,	2353.762,	3158.1575,		5309.7973,		4256.07694,		18.65037 ];			#Antione Equation B constant gives results in bar, must be converter to Pa	[2]
global ACENTRIC;
ACENTRIC = [0.022,		-0.216,		0.345,		0.040,		0.288,		0.416,			0.316,			0.0141345,			-0.390];		#[3] Acentric factor

global R_CONST = 8314.4621	##(L*Pa)/(mol*K)


PressArray = [];
VolumeArray = [];
for ii = 1:10000
ThermalCompsPR(0.1)
  PressArray(ii) = Press;
  VolumeArray(ii) = Volume;
  Volume -= .000998;
endfor
semilogx(VolumeArray,PressArray,'+')



toc
##Press_array = zeros(1,n_itterations);
##Temp_array = zeros(1,n_itterations);
##time_array = zeros(1,n_itterations);
##vapor_fraction_array = zeros(1,n_itterations);
##Q_array = zeros(1,n_itterations);
##
##
##mass = [50,0,0,0,0,0,0,0,0];
##
##Temps = [50,60,70,80,90,100];
##incrMass = [100,0,0,0,0,0,0,0,0];
##nMassIncr = 1000;
##
##P = zeros(nMassIncr,length(Temps));
##V = zeros(nMassIncr,length(Temps));
##T = zeros(nMassIncr,length(Temps));
##Vf = zeros(nMassIncr,length(Temps));
##
##
##
##for nn = 1:nMassIncr
##  for mm = 1:length(Temps)
##    plast = 0;
##    time = 0;
##    vapor_mass = [0,0,0,0,0,0,0,0,0];
##    setTemp(Temps(mm),substance2plot)
##    itteration = 1;
##    while(true)
##      plast = Press;
##      itteration = itteration+1;
##      ThermalComps(delta_time)
##      if(abs(plast-Press)<100 && itteration> 3)
##        itteration
##        break;
##      endif
####      Press_array(itteration) = Press*0.000145038;
####      Temp_array(itteration) = Temp;
####      vapor_fraction_array(itteration) = (sum(vapor_mass))/sum(mass);
####      Q_array(itteration) = Q/sum(mass);
####      time_array(itteration) = time;
##      time += delta_time;
##      setTemp(Temps(mm),substance2plot);
##      plast = Press;
##    endwhile
##
####    AX = plotyy(time_array,Press_array,time_array,Temp_array);
####    grid on;
####    grid minor;
####    xlabel("Time [sec]");
####    ylabel(AX(1),"Pressure [PSI]");
####    ylabel(AX(2),"Temperature [K]");
####    figure(2)
####    waitfor(AX)
####    AX2 = plotyy(time_array,vapor_fraction_array,time_array,Q_array);
####    grid on;
####    grid minor;
####    xlabel("Time [sec]");
####    ylabel(AX2(1),"Vapor Fraction");
####    ylabel(AX2(2),"Specific Energy [J/g]");
####    waitfor(AX2)
##
##    P(nn,mm) = Press*0.000145038;
####    Press*0.000145038;
##    V(nn,mm) = Volume/mass(substance2plot);
##    T(nn,mm) = Temps(mm);
##    Vf(nn,mm) = (sum(vapor_mass))/sum(mass);
##
##  endfor
##  mass = mass + incrMass;
##  nn/nMassIncr*100
##  toc
##endfor
##figure(1)
##loglog(V,P);
##legend(num2str(Temps'))
##grid on;
##grid minor;
##axis([0.001 1 .01 100]);
##xlabel("Specific Volume [l/g]");
##ylabel("Pressure [PSI]");
##figure(2)
##semilogx(V,Vf);
##legend(num2str(Temps'))
##grid on;
##grid minor;
##axis([0.001 0.1 0 1]);
##xlabel("Specific Volume [l/g]");
##ylabel("Vapor Fraction]");
##toc