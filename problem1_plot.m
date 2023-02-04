A = csvread("problem1_mach.txt");

figure()
plot(A(:,1), A(:,2));
ylabel("mach number")