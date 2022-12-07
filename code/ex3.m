
close all;
clear;

ts = 0.02:0.02:4;

% ys = zeros(200, 1);

% for i=1:200
%     ys(i) = model(ts(i), 3, 150, 2);
% end

ys = model(ts, 3, 150, 2);

stdv = 0.05 * max(abs(ys));

noises = normrnd(0, stdv, 1, 200);

yst = ys + noises;

% yst = ys;


f = @(x) computeF(ts, yst, x(1), x(2), x(3));
df =  @(x) computeDf(ts,yst,x(1), x(2), x(3));
d2f = @(x) computeD2f(ts,yst,x(1), x(2), x(3));

J = @(x) computeJacobian(ts, x(1), x(2), x(3));
r = @(x) (model(ts, x(1), x(2), x(3)) - yst)';

F = struct('f', f, 'df', df, 'd2f', d2f, 'J', J, 'r', r);
descent = 'gauss';

ls = @(x, p, z) z;

alpha0 = 0.05;

x0 = [1,1,1]';

tol = 0.00001;

maxIter = 10000;

[GN_x, GN_f, ~, ~] = descentLineSearch(F, descent, ls, alpha0, x0, tol, maxIter)


x0 = [1,1,1]';
Delta = 1;
eta = 0.001;
debug = 0;
F2 = 0;
solver = @(F, x_k, Delta) solverCMlevenberg(F, x_k, Delta, 100);

[LM_x, LM_f, ~, ~] = trustRegion(F, x0, solver, Delta, eta, tol, maxIter, debug, F2)


GN_estimated = model(ts, GN_x(1), GN_x(2), GN_x(3));
LM_estimated = model(ts, LM_x(1), LM_x(2), LM_x(3));



f = figure();
plot(ts, yst);
hold on;
plot(ts, GN_estimated);
legend('Measurements', 'Gauss-Newton Estimate');
xlabel('t');
ylabel('signal');
title('Measurements Vs Estimated Signal by Gauss-Newton');
saveas(f, 'GN.png');

f = figure();
plot(ts, yst);
hold on;
plot(ts, LM_estimated);
legend('Measurements','Levenberg-Marquardt Estimate');
xlabel('t');
ylabel('signal');
title('Measurements Vs Estimated Signal by Levenberg-Marquardt');
saveas(f, 'LM.png');


















