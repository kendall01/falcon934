function [CA, SA, mdot_f, ringD, centerD] = circleAreaFun()

IN_TO_M = .0254;
RING_HOLES = 8; %num holes making up the ring of holes
RING_DIA = 5/16 * IN_TO_M; %initial diameter of holes in ring, in meters
CENTER_DIA = 5/8 * IN_TO_M;
LENGTH = 5.375 * IN_TO_M;
N = 200; %number of time steps
MDOT_O = .0755; %kg/s TODO: Change to vector based off data?
n = 0.55; % TODO: Figure out based on data
a = .132; %regression coefficient for HDPE from http://www.planete-sciences.org/espace/basedoc/images/d/d7/Aiaa2006HybridFuelRegression.pdf
RHO_HDPE = NaN;
UNITS = 10; % converts from kg/m^2-s to g/cm^2-s

t_f = .01; %number of seconds to run burn. 
tstep = t_f/N;
mdot_O2 = linspace(MDOT_O,MDOT_O, N);

ringD = zeros(N,1);centerD = zeros(N,1);CA = zeros(N,1);SA = zeros(N,1);r = zeros(N,1);mdot_f = zeros(N,1);

ringD(1) = RING_DIA;
centerD(1) = CENTER_DIA;

for i = 1:N-1;
    CA(i) = RING_HOLES * pi * ringD(i)^2 / 4 + pi * centerD(i)^2 / 4;
    SA(i) = (RING_HOLES * pi * ringD(i) + pi * centerD(i)) * LENGTH;
    r(i) = a * (mdot_O2(i) / CA(i) / UNITS) ^ n;
    mdot_f(i) = SA(i) * r(i) * tstep * RHO_HDPE;
    
    %increment
    ringD(i+1) = ringD(i) + 2 * r(i) * tstep;
    centerD(i+1) = centerD(i) + 2 * r(i) * tstep;
end
%fence post
CA(N) = RING_HOLES * pi * ringD(N)^2 / 4 + pi * centerD(N)^2 / 4;
SA(N) = RING_HOLES * pi * ringD(N) + pi * centerD(N);
r(N) = (mdot_O2(N) / CA(N)) ^ n;
mdot_f(N) = SA(N) * r(N) * tstep * RHO_HDPE;
    