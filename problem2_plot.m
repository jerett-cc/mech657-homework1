C1 = csvread("problem2_before_shock_mach.csv");
C2 = csvread("problem2_after_shock_mach.csv");
C = [C1;C2];

figure()
plot(C(:,1), C(:,2));
ylabel("mach number")



B1 = csvread("problem2_before_shock_pressure.csv");
B2 = csvread("problem2_after_shock_pressure.csv");
B = [B1;B2];

figure()
plot(B(:,1), B(:,2));
ylabel("pressure in Pa")