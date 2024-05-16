function ode_plot()
%NCI-N87
% Define parameters
kon = 0.3684e9;
koff = 0.014;
Vs = 22.75e4;
ke = 0.051;
mu = 0.019;
kdeg = 0.027;
kout = 0.022;
NAv = 6.02*10^23;

% Define initial conditions
R0 = 3.25e6;
C0 = 0;
I0 = 0;
D0 = 0;
Ab0 = 1e-8;
N0 = 1e8;

% Time span
tspan = [0 250]; % from 0 to 250 hours

% Solve ODEs
[t, y] = ode45(@odes, tspan, [R0 C0 I0 D0 Ab0 N0]);

% Plot results
figure;
subplot(3, 2, 1);
plot(t, y(:, 1));
xlabel('Time (hours)');
ylabel('R');
title('R vs Time');

subplot(3, 2, 2);
plot(t, y(:, 2));
xlabel('Time (hours)');
ylabel('C');
title('C vs Time');

subplot(3, 2, 3);
plot(t, y(:, 3));
xlabel('Time (hours)');
ylabel('I');
title('I vs Time');

subplot(3, 2, 4);
plot(t, y(:, 4));
xlabel('Time (hours)');
ylabel('D');
title('D vs Time');

subplot(3, 2, 5);
plot(t, y(:, 5));
xlabel('Time (hours)');
ylabel('[Ab]');
title('[Ab] vs Time');

subplot(3, 2, 6);
plot(t, y(:, 6));
xlabel('Time (hours)');
ylabel('N');
title('N vs Time');

sgtitle("NCI-N87")

% Define ODEs
function dydt = odes(t, y)
    dydt = zeros(6, 1);
    dydt(1) = -kon*y(5)*y(1) + koff*y(2) + Vs - ke*y(1) - mu*y(1);
    dydt(2) = kon*y(5)*y(1) - koff*y(2) - ke*y(2) - mu*y(2);
    dydt(3) = ke*y(2) - kdeg*y(3) - mu*y(3);
    dydt(4) = kdeg*y(3) - kout*y(4) - mu*y(4);
    dydt(5) = (koff*y(2) - kon*y(5)*y(1))*(y(6)/NAv);
    dydt(6) = mu*y(6);
end

end