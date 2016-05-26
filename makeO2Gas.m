function gas = makeO2Gas()

gas = IdealGasMix('me140_species.xml');

% Create the composition vectors
nsp = nSpecies(gas);

iO2  = speciesIndex(gas,'O2');

% Constants
Po = 6.8e6; % Pa (68 bar)
To = 25 + 273; % K

y_r       = zeros(nsp,1);
y_r(iO2)  = 1;

% Get the enthalpy of the reactants at To (need 2 independent variables to
% define state)
set(gas,'T',To,'P',Po,'Y',y_r);

end