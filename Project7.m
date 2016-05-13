clear all;
close all;

phi = linspace(10,1, 100);

T0 = zeros(length(phi),1);
for i = 1:length(phi)
    
    [T0(i),gas] = combustion(phi(i)); %T0, gas must be returned in combustion ([T0(i),gas]= combustion(phi))
    
    h_0 = enthalpy_mole(gas);
    P_new = pressure(gas);
    S1 = entropy_mole(gas);
    RHS = 0;
    x_r = moleFractions(gas);
    while h_0 < RHS
        set(gas,'S',S1,'P',P_new,'X',x_r);  % As of now, we aren't changing gas based on P-- need to fix this
        P_new = P_new - 15;
        h_t = enthalpy_mole(gas);
        speed = (soundspeed(gas))^2/2;
        RHS = h_t + speed
    end
    Tt = temperature(gas);    % nozzle throat temp found by equilibrating h0= ht + Vt^2/2
    

    c= pressure(gas)/(density(gas)*soundspeed(gas));  %c*= P0/(rho*Ut) dependent on new gas mixture at each mix ratio
    
end

figure(1)
plot(phi, T0)
figure(2)
plot(phi, Tt)
figure(3)
plot(phi, c)
xlabel('Mixture Ratio')
ylabel('Stagnation Temperature')
%3740 %temperature that T0 graph should peak at 

