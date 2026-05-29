clear

B0 = 31200/5.5^3*1e-9;     % background magnetic field, T
ne = 50e6;      % cold electron density, m^-3, 1 cm^-3 = 1e6 m^-3

dd_cold_ion=[
    0.7, 0.1, 0.1;
      1,   1,   1;
      1,   4,  16;
    ]; % [ns/ne;qs/qe;ms/mp];
dd_warm_ion=[
    0.05, 0.05;
       1,    1;
       1,    4;
     1.5,    0;
      10,   10;
    ];% [ns/ne;qs/qe;ms/mp;Tp/Tz-1;Tz(keV)]
%% Scan normalized wave number
%  K = k c / (omega_pe * sqrt(eps))
K = logspace(-2, 1, 2000);

[K,xH]=PHEMICS(K,0.5+1e-3i,B0,ne,dd_cold_ion,dd_warm_ion);
[K,xHe]=PHEMICS(K,0.1+1e-3i,B0,ne,dd_cold_ion,dd_warm_ion);
[K,xO]=PHEMICS(K,0.01+1e-3i,B0,ne,dd_cold_ion,dd_warm_ion);

%%
figure('Color','w')
subplot(3,1,1)
plot(K, real(xH), '--r', 'LineWidth', 1.5); 
hold on;
plot(K, real(xHe), '-r', 'LineWidth', 1.5);
plot(K, real(xO), '-.r', 'LineWidth', 1.5);
xlabel('K = k c / (\omega_{pe}\epsilon^{1/2})');
ylabel('\omega_r / \Omega_p');
legend('H band', 'He band');
grid on;
axis([0.001 100 0 1])
set(gca,'XScale','log','FontName','times new roman','FontSize',14)

%% Plot growth / damping rate
subplot(3,1,2)
plot(K, imag(xH), '--r', 'LineWidth', 1.5); 
hold on;
plot(K, imag(xHe), '-r', 'LineWidth', 1.5);
plot(K, imag(xO), '-.r', 'LineWidth', 1.5);
yline(0, 'k-');
xlabel('K = k c / (\omega_{pe}\epsilon^{1/2})');
ylabel('\gamma / \Omega_p');
legend('H band', 'He band');
grid on;
axis([0.001 100 -0.0025 inf])
set(gca,'XScale','log','FontName','times new roman','FontSize',14)

subplot(3,1,3)
plot(real(xH), imag(xH), '--r', 'LineWidth', 1.5); 
hold on;
plot(real(xHe), imag(xHe), '-r', 'LineWidth', 1.5);
plot(real(xO), imag(xO), '-.r', 'LineWidth', 1.5);
xline(1/4, 'k--','\Omega_{He+}');
xline(1/16, 'k--','\Omega_{O+}');
yline(0, 'k--');
xlabel('\omega / \Omega_p');
ylabel('\gamma / \Omega_p');
legend('H band', 'He band');
grid on;
axis([0 1 -0.0025 1e-2])
set(gca,'XScale','linear','FontName','times new roman','FontSize',14)