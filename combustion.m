function [T0, gas]= combustion(phi)
    %phi is the mixture ratio = mass of oxygen/mass of fuel
    tic
    
    
    % Initialize the gas object -- Only do this once per script
    gas = IdealGasMix('me140_species.xml');
    
    % Constants
    Po = 1.172E6;  % Pa
    To = 25 + 273; % K
    
    % Create the composition vectors
    nsp = nSpecies(gas);
    
    % we want to know the index of relevant species
    iCO2 = speciesIndex(gas,'CO2');
    iO2  = speciesIndex(gas,'O2');
    iCO  = speciesIndex(gas,'CO');
    iH2O = speciesIndex(gas,'H2O');
    iH2 = speciesIndex(gas,'H2');
    iH = speciesIndex(gas,'H');
    iO = speciesIndex(gas,'O');
    iOH = speciesIndex(gas,'OH');
    iC = speciesIndex(gas,'C');
    iC2H4 = speciesIndex(gas,'C2H4');
    
    %Molar masses
molmass = molarMasses(gas);
    
    % Note that cantera uses x for mole fractions (and y for mass fractions)
    x_r       = zeros(nsp,1);
    x_r(iC2H4) = 1 / molmass(iC2H4);
    x_r(iO2)  = phi / molmass(iO2);
    x_r = x_r./sum(x_r); %good to normalize, although cantera should do it automatically
    
    x_p = zeros(nsp,1);

    % Get the enthalpy of the reactants at To (need 2 independent variables to
    % define state)
    set(gas,'T',To,'P',Po,'X',x_r);
    
    %stoichiometric coefficeints
    z_C2H4 = 1;
    z_O2 = 3;
    z_CO2 = 2;
    z_H2O = 2;
    
    %define LHV and known enthalpies of formation
    LHV_HDPE = 44E6; %J/kg;
    hf_C2H4 = 52.3; % kj/mol @ 298  from http://www.kentchemistry.com/links/Kinetics/EnthalpyFormation.htm
    hf_C2H4 = hf_C2H4 / molmass(iC2H4) * 10^6; %J/kg
    hf_H2O = -241.8; %kj/mol  from kentchem. water vapor
    hf_H2O = hf_H2O / molmass(iH2O) * 10^6; %J/kg
    hf_CO2 = -393.5; %kj/mol  from kentchem. water vapor
    hf_CO2 = hf_CO2 / molmass(iCO2) * 10^6; %J/kg
    hf_O2 = 0;
    
    %Calculate hf_HDPE using chemistry
    hf_HDPE  = (- z_O2 * molmass(iO2) * hf_O2 + z_CO2 * molmass(iCO2) * hf_CO2 + z_H2O * molmass(iH2O) * hf_H2O + molmass(iC2H4) * LHV_HDPE) / (z_C2H4 * molmass(iC2H4));
    %hf_HDPE should be about 1.94e6. 
    h_r = enthalpy_mass(gas);
    hDiff = hf_HDPE - hf_C2H4; %J/kg of fuel
    massfrac = massFractions(gas);
    hDiff = hDiff * massfrac(iC2H4); %J/kg of total gas
    set(gas,'P',Po,'H', h_r + hDiff);

    
%     find enthalpy of formation of hdpe using LHV
%    diff = enthalpy of formation of c2h4 and hdpe
%    add in diff
%     

    
    
    % Now lets try using the equilibrate function
    % composition doesnt change unless you equilibrate or define it with x
%     set(gas,'T',To,'P',Po,'X',x_r);
    
    % Allow the reactants to reach equilibrium (minimize g) at constant
    % enthalpy and pressure
    equilibrate(gas,'HP');
    
    % The temperature of this gas is the adiabatic flame temperature
    % more accurate, because it takes into account minor species
    T0 = temperature(gas);
    33;
    toc
end