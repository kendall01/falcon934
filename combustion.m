function T0 = combustion(phi)
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
    
    x_p       = zeros(nsp,1);

    % Get the enthalpy of the reactants at To (need 2 independent variables to
    % define state)
    set(gas,'T',To,'P',Po,'X',x_r);
    
    h_r = enthalpy_mole(gas); % J/kmol
    LHV = 44E6; %J/kg
    
    % Test temperatures until product enthalpy matches
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
    T0 = temperature(gas);
    
    toc
end