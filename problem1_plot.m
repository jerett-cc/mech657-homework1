C = csvread("problem1_mach.csv");

figure()
plot(C(:,1), C(:,2));
ylabel("mach number")

B = csvread("problem1_pressure.csv");

figure()
plot(B(:,1), B(:,2));
ylabel("pressure in kPa")

C = csvread("problem1_temp.csv");

figure()
plot(C(:,1), C(:,2));
ylabel("temperature in K")

D = csvread("problem1_density.csv");

figure()
plot(D(:,1), D(:,2));
ylabel("density")
