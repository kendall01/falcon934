clear all;
close all;

phi = linspace(10,1, 100);

T0 = zeros(length(phi),1);
for i = 1:length(phi)
    T0(i) = combustion(phi(i));
end

figure(1)
plot(phi, T0)
xlabel('Mixture Ratio')
ylabel('Stagnation Temperature')
%3740 %temperature that T0 graph should peak at 


%% Take T0 and do things