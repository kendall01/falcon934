clear all;
close all;
tic
phi = linspace(10,1,100);

T0 = zeros(length(phi),1);

rho1 = zeros(length(phi),1);
rho2 = zeros(length(phi),1);
rho2_2 = zeros(length(phi),1);

% %frozen vectors
h_0 = zeros(length(phi),1);
h_e = zeros(length(phi),1);
Tt_frozen = zeros(length(phi),1);
Te_frozen = zeros(length(phi),1);
c_frozen = zeros(length(phi),1);
Cf_frozen = zeros(length(phi),1);
A_ratio_frozen = zeros(length(phi),1);
Ue = zeros(length(phi),1);
Ut = zeros(length(phi),1);
frozen_mole_fracs = zeros(12,length(phi));


%reacted vectors
h_0 = zeros(length(phi),1);
h_e_reacted = zeros(length(phi),1);
Tt_reacted = zeros(length(phi),1);
Te_reacted = zeros(length(phi),1);
c_reacted = zeros(length(phi),1);
Cf_reacted = zeros(length(phi),1);
A_ratio_reacted = zeros(length(phi),1);
Ut_reacted = zeros(length(phi),1);
Ue_reacted = zeros(length(phi),1);
reacted_mole_fracs = zeros(12,length(phi));

% Po = 1.172E6;  % Pa
Po = 6.8e6; % Pa (68 bar)

%parfor runs this loop in parallel. It only works on 2016a with the
%parallels toolbar, otherwise it will still run, but will just run as a
%normal for loop. On the first run, it will take longer as it has to launch
%the parallel pool which takes like a minute. But on subsequent runs it is
%much faster.
parfor i = 1:length(phi)
    i
    % Begin with frozen throat
    [T0(i), gas, y_r] = combustion(phi(i)); %T0, gas must be returned in combustion ([T0(i),gas]= combustion(phi))
    h_0(i) = enthalpy_mass(gas); %define original enthalpy to compare
    P_new = pressure(gas); %P_new defined here but will change immediately in loop
    S1 = entropy_mass(gas); %define original entropy (mass specific since setState_SP uses that) to reference bc nozzle is isentropic
    RHS = h_0(i) + 100; %start RHS at arbitrary value, anything greater than h_0 which is NOT ALWAYS negative
    
    while h_0(i) < RHS
        set(gas,'P',P_new,'S',S1); %should keep all same, but set gas to a new pressure
        h_t = enthalpy_mass(gas); %find enthalpy of new gas
        speed = (soundspeed(gas))^2/2; %find Kinetic E term of new gas
        RHS = h_t + speed; %calculate RHS for comparison to h_0
        P_new = P_new - 1000; %change P for next iteration
    end
    
    Tt_frozen(i) = temperature(gas); %store temps while gas is at throat P and T
    c_frozen(i)= Po/(density(gas)*soundspeed(gas));  %c*= P0/(rho*Ut) dependent on new gas mixture (beginning of nozzle not throat) at each mix ratio
    
    At= 0.0001824; %m^2
   
    rho1(i) = density(gas);
    Ut(i) = soundspeed(gas);
    
    % find frozen exit gas conditions
    P_e = 101325; %Pa
    set(gas,'P',P_e,'S',S1);
    Te_frozen(i) = temperature(gas);
    rho2(i) = density(gas);
    h_e(i) = enthalpy_mass(gas);
    Ue(i) = sqrt(2*(h_0(i)- h_e(i)));
    A_ratio_frozen(i) = rho1(i)*Ut(i)/(rho2(i)*Ue(i));
    Cf_frozen(i) = Ue(i)/c_frozen(i);
    frozen_mole_fracs(:,i) = moleFractions(gas);
    
    % Now find Reacted Gas Throat conditions
    set(gas,'P',Po,'S',S1,'Y',y_r); %must reset to original gas compostion
    equilibrate(gas,'SP'); %making sure the gas is reset (may be uneccesary, idk)
    h_0(i) = enthalpy_mass(gas); %define original enthalpy to compare
    RHS = h_0(i) + 100; %start RHS at arbitrary value, anything greater than h_0 which is NOT ALWAYS negative
    P_new = pressure(gas); %P_new defined here but will change immediately in loop
    
    while h_0(i) < RHS
        set(gas,'P',P_new,'S',S1); %should keep all same, but set gas to a new pressure
        equilibrate(gas,'SP'); %react gas to a new composition with each change in pressure
        h_t = enthalpy_mass(gas); %find enthalpy of new gas
        speed = (soundspeed(gas))^2/2; %find Kinetic E term of new gas
        RHS = h_t + speed; %calculate RHS for comparison to h_0
        %P_new = P_new - 1000; %change P for next iteration
        P_new = P_new - ((RHS - h_0(i)) + 100); % does the same thing but almost three times faster overall
    end
    
    Tt_reacted(i) = temperature(gas); %store temps while gas is at throat P and T
    Ut_reacted(i) = soundspeed(gas); %velocity at throat
    c_reacted(i)= Po/(density(gas)*Ut_reacted(i));  %c*= P0/(rho*Ut) dependent on new gas mixture (beginning of nozzle not throat) at each mix ratio
    
    % Reacted exit
    set(gas,'P',P_e,'S',S1);
    equilibrate(gas,'SP'); %don't know if this is necessary to re-do
    Te_reacted(i) = temperature(gas);
    rho2_2(i) = density(gas);
    h_e_reacted(i) = enthalpy_mass(gas);
    Ue_reacted(i) = sqrt(2*(h_0(i)- h_e_reacted(i)));
    A_ratio_reacted(i) = rho1(i)*Ut(i)/(rho2_2(i)*Ue_reacted(i));
    Cf_reacted(i) = Ue_reacted(i)/c_reacted(i);
    reacted_mole_fracs(:,i) = moleFractions(gas);
end

figure(1);
plot(phi, T0, 'g', phi, Tt_frozen, 'b-.', phi, Tt_reacted, 'b', phi, Te_frozen, 'r-.', phi, Te_reacted, 'r', 'Linewidth', 1.2);
title('Frozen Case');
xlabel('Mixture Ratio');
ylabel('Temperature [K]');
legend('Stagnation Temperature', 'Throat Frozen', 'Throat Reacted', 'Exit Frozen', 'Exit Reacted');
% figure(10);
% plot(phi, T0, phi, Tt_reacted, phi, Te_reacted, 'Linewidth', 1.2);
% title('Reacted Case');
% xlabel('Mixture Ratio');
% ylabel('Temperature [K]');
% legend('Stagnation Temperature', 'Nozzle Throat Temperature', 'Exit Temperature');
figure(2)
plot(phi, c_frozen, phi, c_reacted, 'Linewidth', 1.2);
xlabel('Mixture Ratio');
ylabel('C* (m/s)');
title('C* vs Mixture Ratio');
legend('Frozen', 'Reacted');
% figure(3)
% plot(phi, Ut, phi, Ut_reacted, 'Linewidth', 1.2);
% xlabel('Mixture Ratio');
% ylabel('Throat Velocity (m/s)');
% legend('Frozen', 'Reacted');
% title('Throat Velocity');
figure(6)
plot(phi, Ue, phi, Ue_reacted, 'Linewidth', 1.2);
xlabel('Mixture Ratio');
ylabel('Exit Velocity (m/s)');
legend('Frozen', 'Reacted');
title('Exit Velocity (Ve)');
figure(4)
plot(phi, Cf_frozen, phi, Cf_reacted, 'Linewidth', 1.2);
title('Cf Frozen');
xlabel('Mixture Ratio');
ylabel('Coefficient of Thrust (Cf)');
legend('Frozen', 'Reacted');
figure(5)
plot(phi, A_ratio_frozen, phi, A_ratio_reacted, 'Linewidth', 1.2)
title('Area ratio');
xlabel('Mixture Ratio')
ylabel('Area Ratio (Ae/At)')
legend('Frozen', 'Reacted');
figure(7)
plot(phi, frozen_mole_fracs, 'Linewidth', 1.2);
xlabel('Mixture Ratio');
ylabel('Mole Fractions');
legend('H', 'H2', 'O', 'O2', 'OH', 'C', 'CO', 'CO2', 'H2O', 'C2H4');
title('Mole Fractions of Frozen Nozzle');
figure(8)
plot(phi, reacted_mole_fracs, 'Linewidth', 1.2);
xlabel('Mixture Ratio');
ylabel('Mole Fractions');
legend('H', 'H2', 'O', 'O2', 'OH', 'C', 'CO', 'CO2', 'H2O', 'C2H4');
title('Mole Fractions of Reacted Nozzle');
plotfixer
%3740 --> temperature that T0 graph should peak at 
toc
