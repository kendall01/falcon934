function [CA, SA, mdot_O2, mdot_f, iRingD, oRingD, centerD, ringHoles, t_f] = doubleCircleAreaFun(N, M)

IN_TO_M = .0254;
G_TO_KG = 1000;

%Constants to play with:
O_RING_HOLES = 13; %num holes making up the outer ring of holes
I_RING_HOLES = 7; %num holes making up the inner ring of holes
CENTER_DIA = 15/32 * IN_TO_M; % Diameter of the hole in the center
I_RING_DIA = 7/32 * IN_TO_M; %initial diameter of  holes in inner ring, in meters
O_RING_DIA = 5/32 * IN_TO_M; %initial diameter of  holes in outer ring, in meters

LENGTH = 7.375 * IN_TO_M;
RING_HOLES = O_RING_HOLES + I_RING_HOLES;
n = 0.5; % TODO: Figure out based on data. Has significant effect on mixture ratio
a = .132; %regression coefficient for HDPE from http://www.planete-sciences.org/espace/basedoc/images/d/d7/Aiaa2006HybridFuelRegression.pdf
RHO_HDPE = 970; %kg/m^3
UNITS = 10; % converts from kg/m^2-s to g/cm^2-s
MM_TO_M = 1000;
m_O2 = 0.2491; %mass of O2 from lab burn. Assuming this is all the oxygen in the tank and that it will always be the same and that it is the limiting factor on length of burn.
mdot_O2_init = 0.088; %from lab 1 fire
mdot_O2_fin = 0.06; %from lab 1 fire. Assuming linear.

t_f = m_O2 / mean([mdot_O2_init, mdot_O2_fin]); %number of seconds to run burn.
tstep = t_f/N;
lstep = LENGTH / M;



MM_c2h4 = (2 * 12 + 4 * 1) / G_TO_KG; %[kg/mol]
MM_o2 = (32) / G_TO_KG; %[kg/mol]
stoich_O2 = 3; %C2H4 + 3O2 -> 2 CO2 + 2 H2O

%Runs for time = 3.4 seconds on every run.

%initialize vectors
moldot_f = zeros(N,M);mdot_O2 = zeros(N,M);iRingD = zeros(N,M);oRingD = zeros(N,M);centerD = zeros(N,M);CA = zeros(N,M);SA = zeros(N,M);r = zeros(N,M);mdot_f = zeros(N,M);

%set initial conditions
iRingD(1,:) = I_RING_DIA;
oRingD(1,:) = O_RING_DIA;
centerD(1,:) = CENTER_DIA;
mdot_O2(:,1) = linspace(mdot_O2_init,mdot_O2_fin, N)'; %assuming linear
mdot_gas = mdot_O2;

for i = 1:N-1
    i
    gas = makeO2Gas();
    
    for j = 1:M-1
        CA(i, j) = (I_RING_HOLES * iRingD(i,j)^2 + O_RING_HOLES * oRingD(i,j)^2) * pi / 4 + pi * centerD(i,j)^2 / 4; %Calculate cross sectional area at current step
        SA(i,j) = ((I_RING_HOLES * iRingD(i,j) + O_RING_HOLES * oRingD(i,j)) * pi  + pi * centerD(i,j)) * lstep; %Calculate surface area at current step
        r(i,j) = (a * (mdot_O2(i,j) / CA(i,j) / UNITS) ^ n) / MM_TO_M; % m/s  Burn rate
        mdot_f(i,j) = SA(i,j) * r(i,j) * RHO_HDPE; %[kg/s] per length_step
        moldot_f(i,j) = mdot_f(i,j)/MM_c2h4;
        
        %increment
        iRingD(i+1,j) = iRingD(i,j) + 2 * r(i,j) * tstep;
        oRingD(i+1,j) = oRingD(i,j) + 2 * r(i,j) * tstep;
        centerD(i+1,j) = centerD(i,j) + 2 * r(i,j) * tstep;
        %mdot_O2(i,j+1) = mdot_O2(i,j) - (stoich_O2*moldot_f(i,j)) * MM_o2;
        %mdot_O2(i,j+1) = getO2Flow(mdot_f(i,j),mdot_O2(i,j));
        [mdot_O2(i,j+1),mdot_gas(i,j+1), gas] = getO2Flow(mdot_f(i,j), mdot_gas(i,j), gas);
    end
    j = M;
    CA(i, j) = (I_RING_HOLES * iRingD(i,j)^2 + O_RING_HOLES * oRingD(i,j)^2) * pi / 4 + pi * centerD(i,j)^2 / 4; %Calculate cross sectional area at current step
    SA(i,j) = ((I_RING_HOLES * iRingD(i,j) + O_RING_HOLES * oRingD(i,j)) * pi  + pi * centerD(i,j)) * lstep; %Calculate surface area at current step
    r(i,j) = (a * (mdot_O2(i,j) / CA(i,j) / UNITS) ^ n) / MM_TO_M; % m/s  Burn rate
    mdot_f(i,j) = SA(i,j) * r(i,j) * RHO_HDPE; %[kg/s] per length_step
    moldot_f(i,j) = mdot_f(i,j)/MM_c2h4;
    
    %increment
    iRingD(i+1,j) = iRingD(i,j) + 2 * r(i,j) * tstep;
    oRingD(i+1,j) = oRingD(i,j) + 2 * r(i,j) * tstep;
    centerD(i+1,j) = centerD(i,j) + 2 * r(i,j) * tstep;
end
% start fence post (loop only runs N-1 iterations, this runs the Nth iteration)
i = N;
gas = makeO2Gas();

for j = 1:M-1
    CA(i, j) = (I_RING_HOLES * iRingD(i,j)^2 + O_RING_HOLES * oRingD(i,j)^2) * pi / 4 + pi * centerD(i,j)^2 / 4; %Calculate cross sectional area at current step
    SA(i,j) = ((I_RING_HOLES * iRingD(i,j) + O_RING_HOLES * oRingD(i,j)) * pi  + pi * centerD(i,j)) * lstep; %Calculate surface area at current step
    r(i,j) = (a * (mdot_O2(i,j) / CA(i,j) / UNITS) ^ n) / MM_TO_M; % m/s  Burn rate
    mdot_f(i,j) = SA(i,j) * r(i,j) * RHO_HDPE; %[kg/s] per length_step
    moldot_f(i,j) = mdot_f(i,j)/MM_c2h4;
    
    %increment
    iRingD(i+1,j) = iRingD(i,j) + 2 * r(i,j) * tstep;
    oRingD(i+1,j) = oRingD(i,j) + 2 * r(i,j) * tstep;
    centerD(i+1,j) = centerD(i,j) + 2 * r(i,j) * tstep;
    %mdot_O2(i,j+1) = mdot_O2(i,j) - (stoich_O2*moldot_f(i,j)) * MM_o2;
    %mdot_O2(i,j+1) = getO2Flow(mdot_f(i,j),mdot_O2(i,j));
    [mdot_O2(i,j+1),mdot_gas(i,j+1), gas] = getO2Flow(mdot_f(i,j), mdot_gas(i,j), gas);
end
j = M;
CA(i, j) = (I_RING_HOLES * iRingD(i,j)^2 + O_RING_HOLES * oRingD(i,j)^2) * pi / 4 + pi * centerD(i,j)^2 / 4; %Calculate cross sectional area at current step
SA(i,j) = ((I_RING_HOLES * iRingD(i,j) + O_RING_HOLES * oRingD(i,j)) * pi  + pi * centerD(i,j)) * lstep; %Calculate surface area at current step
r(i,j) = (a * (mdot_O2(i,j) / CA(i,j) / UNITS) ^ n) / MM_TO_M; % m/s  Burn rate
mdot_f(i,j) = SA(i,j) * r(i,j) * RHO_HDPE; %[kg/s] per length_step

%increment
iRingD(i+1,j) = iRingD(i,j) + 2 * r(i,j) * tstep;
oRingD(i+1,j) = oRingD(i,j) + 2 * r(i,j) * tstep;
centerD(i+1,j) = centerD(i,j) + 2 * r(i,j) * tstep;
% end fence post

%fill output variables
mdot_f = sum(mdot_f, 2);
mdot_O2 = mdot_O2(:,1);
ringHoles = [I_RING_HOLES, O_RING_HOLES];