clear all;
close all;

N=100;
M=100;
[~, ~, m_dot_oxidizer, m_dot_fuel,iRingD, oRingD, centerD, ~, t_f] = doubleCircleAreaFun(N, M);
t_step = t_f/N;
m_dot_total = m_dot_fuel+m_dot_oxidizer; % kg/s

phi = m_dot_oxidizer./m_dot_fuel;

T0 = zeros(length(phi),1);

rho1 = zeros(length(phi),1);
rho2 = zeros(length(phi),1);

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

Po = 1172110; % Pa (68 bar)
P_e = 101325; %Pa

%parfor runs this loop in parallel. It only works on 2016a with the
%parallels toolbar, otherwise it will still run, but will just run as a
%normal for loop. On the first run, it will take longer as it has to launch
%the parallel pool which takes like a minute. But on subsequent runs it is
%much faster.
parfor i = 1:length(phi)
    i
    [T0(i), gas, y_r] = combustion(phi(i)); %T0, gas must be returned in combustion ([T0(i),gas]= combustion(phi))
    S1 = entropy_mass(gas); %define original entropy (mass specific since setState_SP uses that) to reference bc nozzle is isentropic
    
    % Now find Reacted Gas Throat conditions
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
    rho1(i) = density(gas);
    
    A_t(i) = m_dot_total(i) / rho1(i) / Ut_reacted(i); %Ideal nozzle
    dia_t(i) = sqrt(A_t(i) / pi)*2;
    
%     %calcs for TA nozzle
%     A_t(i) = pi * .684^2/4; %TA nozzle
%     A_e(i) = pi * 1.73^2 /4; 
%     k = 1.4;
% %     A_ratio = A_e(i)/A_t(i);
% %     syms Ma
%     Ma = 1;
%     P_e = Po *  (1 + (k - 1)/2 * Ma ^2 ) ^ (-k/(k - 1));
P_e = P_new/1.5;
    % Reacted exit
    set(gas,'P',P_e,'S',S1);
    Te_reacted(i) = temperature(gas);
    rho2(i) = density(gas);
    h_e_reacted(i) = enthalpy_mass(gas);
    Ue_reacted(i) = sqrt(2*(h_0(i)- h_e_reacted(i)));
    A_ratio_reacted(i) = rho1(i)*Ut_reacted(i)/(rho2(i)*Ue_reacted(i));
    %Cf_reacted(i) = Ue_reacted(i)/c_reacted(i);
    Cf_reacted(i) = Ue_reacted(i)/c_reacted(i); % for TA nozzle
    reacted_mole_fracs(:,i) = moleFractions(gas);
    A_e(i) = A_t(i) * A_ratio_reacted(i); 
    dia_e(i) = sqrt(A_e(i) / pi)*2;
end
g = 9.81; % m/s^2
I_sp = c_reacted.*Cf_reacted/g;
times = t_step*(1:length(I_sp));

figure(1)
plot(times, phi)
xlabel('Time')
ylabel('Phi')
title('Mixture Ratio')

figure(2)
plot(times,I_sp)
xlabel('time (s)')
ylabel('I (specific impulse)');
