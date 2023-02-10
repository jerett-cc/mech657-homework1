close all;
A = csvread("problem3_mach.csv");

figure()
plot(A(:,1), A(:,2))
ylabel("Mach number")

B = csvread("problem3_density.csv");

figure()
plot(B(:,1), B(:,2))
ylabel("density")
ylim([0 1.1])
yticks(0:0.1:1.1)

C = csvread("problem3_pressure.csv");

figure()
plot(C(:,1), C(:,2))
ylabel("pressure in Pa")
ylim([1 11e4])