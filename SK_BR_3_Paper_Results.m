function ode_plot()
%SK-BR-3
% Define parameters
kon = 0.3684e9;
koff = 0.014;
Vs = 35.855e4;
ke = 0.09;
mu = 0.011;
kdeg = 0.038;
kout = 0.015;
NAv = 6.02*10^23;

% Define initial conditions
R0 = 3.55e6;
C0 = 1;
I0 = 0;
D0 = 0;
Ab0 = 1e-18;
N0 = 1e8;

% Time span
tspan1 = [0 130]; % from 0 to 130 hours
tspan2 = [0 170]; % from 0 to 170 hours

% Solve ODEs
[t1, y1] = ode45(@odes1, tspan1, [R0 C0 I0 D0 Ab0 N0]);
[t2, y2] = ode45(@odes2, tspan2, [R0 C0 I0 D0 Ab0 N0]);

% Plot results
figure;
plot(t1, y1(:, 2) + y1(:, 3));
hold on
plot(t1, y1(:, 2));
hold on;
plot(t1, y1(:, 3));
hold off;
xlim([0 140])
ylim([0 1])
set(gca, 'ytick', 0:0.2:1);
xlabel('Time (h)');
ylabel('Fraction of Full Antibody After wash');
legend('Total Full ADC (C + I)', 'Complex (C)', 'Internalized, Intact ADC (I)')
title('SK BR 3');

figure;
plot(t2, y2(:, 2) + y2(:, 3) + y2(:, 4));
hold on
plot(t2, y2(:, 2));
hold on;
plot(t2, y2(:, 3));
hold on;
plot(t2, y2(:, 4));
hold off;
xlim([0 170])
ylim([0 1])
set(gca, 'ytick', 0:0.2:1);
xlabel('Time (h)');
ylabel('Fraction of 647 Signal');
legend('Total Alexa 647 Signal (C + I + D)', 'Complex (C)', 'Internalized, Intact ADC (I)', 'Degraded (D)')
title('SK BR 3');

% Define ODEs
function dydt = odes1(t1, y1)
    dydt = zeros(6, 1);
    dydt(1) = -kon*y1(5)*y1(1) + koff*y1(2) + Vs - ke*y1(1) - mu*y1(1);
    dydt(2) = kon*y1(5)*y1(1) - koff*y1(2) - ke*y1(2) - mu*y1(2);
    dydt(3) = ke*y1(2) - kdeg*y1(3) - mu*y1(3);
    dydt(4) = kdeg*y1(3) - kout*y1(4) - mu*y1(4);
    dydt(5) = (koff*y1(2) - kon*y1(5)*y1(1))*(y1(6)/NAv);
    dydt(6) = mu*y1(6);
end

function dydt = odes2(t2, y2)
    dydt = zeros(6, 1);
    dydt(1) = -kon*y2(5)*y2(1) + koff*y2(2) + Vs - ke*y2(1) - mu*y2(1);
    dydt(2) = kon*y2(5)*y2(1) - koff*y2(2) - ke*y2(2) - mu*y2(2);
    dydt(3) = ke*y2(2) - kdeg*y2(3) - mu*y2(3);
    dydt(4) = kdeg*y2(3) - kout*y2(4) - mu*y2(4);
    dydt(5) = (koff*y2(2) - kon*y2(5)*y2(1))*(y2(6)/NAv);
    dydt(6) = mu*y2(6);
end

end