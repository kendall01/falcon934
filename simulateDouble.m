function [] = simulateDouble(ringD, centerD, numHoles, N)

OUT = 2; % pound define for outer ring
IN = 1;  % pound define for inner ring
IN_TO_M = .0254;
RING_FACTOR = numHoles(OUT)/numHoles(IN); % ratio of radius of outer ring to inner ring. Also, of number of holes since circumference is proportional to radius.
IR = .0115; % [m] radius of inner ring. .018 works well. how far from the center of the grain does the circle of ring holes lie
R = RING_FACTOR * IR; % [m] radius of outer ring
OR = 1.995/2 * IN_TO_M; % [m] outer radius of fuel grain

figure(1)
for i = 1:length(centerD)
    ezpolar(@(x)OR) %Draws exterior of fuel grain
    hold on
    [center(:,1), center(:,2)] = pol2cart(linspace(-2*pi, 2*pi,N)', linspace(centerD(i)/2,centerD(i)/2,N)'); %Calculates coordinates of center hole by expressing it in polar and then converting to cartesian
    plot(center(:,1), center(:,2)) % plots center hole
    hold on
    % Outer Ring
    [ring(:,1), ring(:,2)] = pol2cart(linspace(-2*pi, 2*pi,N)', linspace(ringD(i)/2,ringD(i)/2,N)'); %Calculates polar coordinates of outer ring hole if it were at the origin and converts to cartesian coordinates
    for j = 1:numHoles(OUT)
        %iterates through each hole in the outer ring, moves it from the
        %origin to its location on the ring, then plots it
        hold on
        x_off = R * cos((2*pi/numHoles(OUT)) * j);
        y_off = R * sin((2*pi/numHoles(OUT)) * j);
        plot(ring(:,1) + x_off, ring(:,2) + y_off)
    end
    hold on
    % Inner Ring
    [inner_ring(:,1), inner_ring(:,2)] = pol2cart(linspace(-2*pi, 2*pi,N)', linspace(ringD(i)/2,ringD(i)/2,N)'); %Calculates polar coordinates of inner ring hole if it were at the origin and converts to cartesian coordinates
    for j = 1:numHoles(IN)
        %iterates through each hole in the inner ring, moves it from the
        %origin to its location on the ring, then plots it
        hold on
        x_off = IR * cos((2*pi/numHoles(IN)) * j);
        y_off = IR * sin((2*pi/numHoles(IN)) * j);
        plot(inner_ring(:,1) + x_off, inner_ring(:,2) + y_off)
    end
    axis tight
    hold off
    pause(.0001); %Ideally it simulates it at real world time. The time step of the simulation when N=50 happens to be about the length of time it takes to run the above for loop of plotting. Thus, a tiny pause is put in here so that the plot actually appears because without any pause it doesn't appear.
end


end