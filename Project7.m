clear all;
close all;

phi = linspace(10,1,100);

T0 = zeros(length(phi),1);
Tt = zeros(length(phi),1);
c = zeros(length(phi),1);

for i = 1:length(phi)
    
    [T0(i),gas] = combustion(phi(i)); %T0, gas must be returned in combustion ([T0(i),gas]= combustion(phi))
    state_ref = gas; %set a reference state 1 returned from combustion for a single phi value
    h_0(i) = enthalpy_mole(state_ref); %define original enthalpy to compare
    P_new = pressure(state_ref); %P_new defined here but will change immediately in loop
    S1 = entropy_mass(state_ref); %define original entropy (mass specific since setState_SP uses that) to reference bc nozzle is isentropic
    RHS = h_0(i) + 100; %start RHS at arbitrary value, anything greater than h_0 which is NOT ALWAYS negative
    x_r = moleFractions(state_ref); %define original mole fractions (may be uneccessary)
    state_new = state_ref; %initialized but should change once the loop is entered
    
    while h_0 < RHS
        % set(gas,'P',P_new); % set the current gas to a new P 
        SP= [S1 P_new]; %vector of constant entropy and new pressure to set the gas to a new state
        state_new = setState_SP(state_new, SP); %should keep all same, but set gas to a new pressure
        % state_new = equilibrate(state_new, 'SP'); %should equilibrate gas
        % to a all new thermo props given constant S and same new P
        % --- (NEWS: don't think we want to equilibrate in the nozzle, don't
        % want a reacting converging nozzle  
        h_t = enthalpy_mole(state_new); %find enthalpy of new gas
        speed = (soundspeed(state_new))^2/2; %find Kinetic E term of new gas
        RHS = h_t + speed %calculate RHS for comparison to h_0
        P_new = P_new - 1000; %change P for next iteration
    end
    
    %why is state_ref changing with state_new everytime?
    Tt(i) = temperature(state_new);  % nozzle throat temp found by equilibrating h0 = ht + Vt^2/2
    
    c(i)= pressure(gas)/(density(gas)*soundspeed(gas));  %c*= P0/(rho*Ut) dependent on new gas mixture (beginning of nozzle not throat) at each mix ratio
    
end

figure(1);
plot(phi, T0, phi, Tt)
xlabel('Mixture Ratio')
ylabel('Temperature [K]') %This is in kelvin, isn't it?
legend('Stagnation Temperature', 'Nozzle Throat Temperature')
figure(2)
plot(phi, c)
xlabel('Mixture Ratio')
ylabel('C*')

%3740 --> temperature that T0 graph should peak at 

