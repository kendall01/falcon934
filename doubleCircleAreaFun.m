function [CA, SA, mdot_O2, mdot_f, ringD, centerD, ringHoles, t_f] = doubleCircleAreaFun(N)

IN_TO_M = .0254;

%Constants to play with:
OUTER_RING_HOLES = 16; %num holes making up the outer ring of holes
INNER_RING_HOLES = 9; %num holes making up the inner ring of holes
RING_DIA = 3/32 * IN_TO_M; %initial diameter of all holes in ring, in meters
CENTER_DIA = 3/32 * IN_TO_M; % Diameter of the hole in the center

LENGTH = 5.375 * IN_TO_M;
RING_HOLES = OUTER_RING_HOLES + INNER_RING_HOLES;
n = 0.55; % TODO: Figure out based on data. Has significant effect on mixture ratio
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

%initialize vectors
ringD = zeros(N,1);centerD = zeros(N,1);CA = zeros(N,1);SA = zeros(N,1);r = zeros(N,1);mdot_f = zeros(N,1);

%set initial conditions
ringD(1) = RING_DIA;
centerD(1) = CENTER_DIA;

for i = 1:N-1;
    CA(i) = RING_HOLES * pi * ringD(i)^2 / 4 + pi * centerD(i)^2 / 4; %Calculate cross sectional area at current step
    SA(i) = (RING_HOLES * pi * ringD(i) + pi * centerD(i)) * LENGTH; %Calculate surface area at current step
    r(i) = (a * (mdot_O2(i) / CA(i) / UNITS) ^ n) / MM_TO_M; % m/s  Burn rate
    mdot_f(i) = SA(i) * r(i) * RHO_HDPE; %[kg/s]
    
    %increment
    ringD(i+1) = ringD(i) + 2 * r(i) * tstep;
    centerD(i+1) = centerD(i) + 2 * r(i) * tstep;
end
% start fence post (loop only runs N-1 iterations, this runs the Nth iteration)
i = N;
CA(i) = RING_HOLES * pi * ringD(i)^2 / 4 + pi * centerD(i)^2 / 4;
SA(i) = (RING_HOLES * pi * ringD(i) + pi * centerD(i)) * LENGTH;
r(i) = (a * (mdot_O2(i) / CA(i) / UNITS) ^ n) / MM_TO_M; % m/s
mdot_f(i) = SA(i) * r(i) * RHO_HDPE;
% end fence post

%fill output variables
ringHoles = [INNER_RING_HOLES, OUTER_RING_HOLES];