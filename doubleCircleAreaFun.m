function [CA, SA, mdot_O2, mdot_f, iRingD, oRingD, centerD, ringHoles, t_f] = doubleCircleAreaFun(N)

IN_TO_M = .0254;

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
mdot_O2 = linspace(mdot_O2_init,mdot_O2_fin, N)'; %assuming linear

%Runs for time = 3.4 seconds on every run.

%initialize vectors
iRingD = zeros(N,1);oRingD = zeros(N,1);centerD = zeros(N,1);CA = zeros(N,1);SA = zeros(N,1);r = zeros(N,1);mdot_f = zeros(N,1);

%set initial conditions
iRingD(1) = I_RING_DIA;
oRingD(1) = O_RING_DIA;
centerD(1) = CENTER_DIA;

for i = 1:N-1;
    CA(i) = (I_RING_HOLES * iRingD(i)^2 + O_RING_HOLES * oRingD(i)^2) * pi / 4 + pi * centerD(i)^2 / 4; %Calculate cross sectional area at current step
    SA(i) = ((I_RING_HOLES * iRingD(i) + O_RING_HOLES * oRingD(i)) * pi  + pi * centerD(i)) * LENGTH; %Calculate surface area at current step
    r(i) = (a * (mdot_O2(i) / CA(i) / UNITS) ^ n) / MM_TO_M; % m/s  Burn rate
    mdot_f(i) = SA(i) * r(i) * RHO_HDPE; %[kg/s]
    
    %increment
    iRingD(i+1) = iRingD(i) + 2 * r(i) * tstep;
    oRingD(i+1) = oRingD(i) + 2 * r(i) * tstep;
    centerD(i+1) = centerD(i) + 2 * r(i) * tstep;
end
% start fence post (loop only runs N-1 iterations, this runs the Nth iteration)
i = N;
CA(i) = (I_RING_HOLES * iRingD(i)^2 + O_RING_HOLES * oRingD(i)^2) * pi / 4 + pi * centerD(i)^2 / 4; %Calculate cross sectional area at current step
SA(i) = ((I_RING_HOLES * iRingD(i) + O_RING_HOLES * oRingD(i)) * pi  + pi * centerD(i)) * LENGTH; %Calculate surface area at current step
r(i) = (a * (mdot_O2(i) / CA(i) / UNITS) ^ n) / MM_TO_M; % m/s  Burn rate
mdot_f(i) = SA(i) * r(i) * RHO_HDPE; %[kg/s]
% end fence post

%fill output variables
ringHoles = [I_RING_HOLES, O_RING_HOLES];