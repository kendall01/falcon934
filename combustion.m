function T0 = combustion(phi)
%%
tic
% Calculate the adiabatic flame temperature of stoichiometric 
% methane combustion with engineering air
% CH4 + 2(O2 + 3.76N2) -> 2H2O + CO2 + 2*3.76N2
            
clear all; close all; clc;

% Initialize the gas object -- Only do this once per script
gas = IdealGasMix('gri30.xml');

% Constants
Po = 101300;                % Pa
To = 25 + 273;              % K

% Create the composition vectors
% There are 53 species includd in GRI30, so this will be the size of our array
nsp = nSpecies(gas); 

% we want to know the index of relevant species
iCH4  = speciesIndex(gas,'CH4'); 
iCO2 = speciesIndex(gas,'CO2');
iO2  = speciesIndex(gas,'O2'); 
iN2  = speciesIndex(gas,'N2');
iH2O = speciesIndex(gas,'H2O');

% Note that cantera uses x for mole fractions (and y for mass fractions)
x_r       = zeros(nsp,1);
x_r(iCH4) = 1;
x_r(iO2)  = 2;
x_r(iN2)  = 2*3.76;
x_r = x_r./sum(x_r); %good to normalize, although cantera should do it automatically

x_p       = zeros(nsp,1);
x_p(iH2O) = 2;
x_p(iCO2) = 1;
x_p(iN2)  = 2*3.76;
x_p = x_p./sum(x_p);

% Get the enthalpy of the reactants at To (need 2 independent variables to
% define state)
set(gas,'T',To,'P',Po,'X',x_r);

h_r = enthalpy_mole(gas); % J/kmol

% Test temperatures until prodcut enthalpy matches
h_p = -1e10; % start at a very low enthalpy
T = To;  % start at To and march up in T
dT = 1; % define temp step

while h_p < h_r
    
    T = T + dT;
    set(gas,'T',T,'P',Po,'X',x_p);
    h_p = enthalpy_mole(gas);    % as compared to computing the enthalpy via integrals for each species

end

Adiabatic_Flame_Temp = T

% Now lets try using the equilibrate function
% composition doesnt change unless you equilibrate or define it with x
set(gas,'T',To,'P',Po,'X',x_r);

% Allow the reactants to reach equilibrium (minimize g) at constant
% enthalpy and pressure
equilibrate(gas,'HP');

% The temperature of this gas is the adiabatic flame temperature
% more accurate, because it takes into account minor species
T = temperature(gas)

toc

