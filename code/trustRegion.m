
function [x_k, f_k, k, info] = trustRegion(F, x0, solverCM, Delta, eta, tol, maxIter, debug, F2)
% TRUSTREGION Trust region iteration
% [x_k, f_k, k, info] = trustRegion(F, x0, solverCM, Delta, eta, tol, maxIter, debug, F2)
% INPUTS
% F: structure with fields
%   - f: function handler
%   - df: gradient handler
%   - d2f: Hessian handler
% x_k: current iterate
% solverCM: handle to solver to quadratic constraint trust regio problem
% Delta: upper limit on trust region radius
% eta: step acceptance relative progress threshold
% tol: stopping condition on minimal allowed step
%      norm(x_k - x_k_1)/norm(x_k) < tol;
% maxIter: maximum number of iterations
% debug: debugging parameter switches on visualization of quadratic model
%        and various step options. Only works for functions in R^2
% F2: needed if debug == 1. F2 is equivalent of F but formulated as function of (x,y)
%     to enable meshgrid evaluation
% OUTPUT
% x_k: minimum
% f_k: objective function value at minimum
% k: number of iterations
% info: structure containing iteration history
%   - xs: taken steps
%   - xind: iterations at which steps were taken
%   - stopCond: shows if stopping criterium was satisfied, otherwsise
% k = maxIter 
% Reference: Algorithm 4.1 in Nocedal Wright
%
% Copyright (C) 2017 Marta M. Betcke, Kiko Rullan
% Parameters
% Choose stopping condition {'step', 'grad'}

stopType = 'grad';
% Initialisation
Delta_k = 0.5*Delta;
stopCond = false;
k = 0;
x_k = x0;
nTaken = 0;
info.xs = zeros(length(x0), maxIter);
info.xs(:,1) = x0;
info.xind = zeros(1,maxIter);
info.xind(1) = 1;
while ~stopCond && (k < maxIter)
  k = k+1;
  % Construct and solve quadratic model
  Mk.m = @(p) F.f(x_k) + F.df(x_k)'*p + 0.5*p'*F.d2f(x_k)*p;
  Mk.dm = @(p) F.df(x_k) + F.d2f(x_k)*p;
  Mk.d2m = @(p) F.d2f(x_k);
  p = solverCM(F, x_k, Delta_k);
  if debug
    % Visualise quadratic model and various steps
    figure(1); clf;
    plotTaylor(F2, x_k, [x_k - 4*Delta_k, x_k + 4*Delta_k], Delta_k, p);
    hold on,
    g = -F.df(x_k);
    gu = -F.d2f(x_k)\g;

    plot(x_k(1) + g(1)*Delta_k/norm(g), x_k(2) + g(2)*Delta_k/norm(g), 'rs')
    plot(x_k(1) + gu(1)*Delta_k/norm(gu), x_k(2) + gu(2)*Delta_k/norm(gu), 'bo')
    pause
  end
  % Evaluate actual to predicted reduction ratio
  rho_k = (F.f(x_k) - F.f(x_k + p)) / (Mk.m(0*p) - Mk.m(p)) ;
  if (Mk.m(0*p) < Mk.m(p))
    disp(strcat('Ascent - iter', num2str(k)))
  end
  % Record iteration information
  info.rhos(k) = rho_k;
  info.Deltas(k) = Delta_k;
  if rho_k < 0.25
    % Shink trust region
    Delta_k = 0.25*Delta_k;
  else
    if rho_k > 0.75 && abs(p'*p - Delta_k^2) < 1e-12
      % Expand trust region
      Delta_k = max(2*Delta_k, Delta);
    end
  end
  % Accept step if rho_k > eta
  if rho_k > eta
    x_k_1 = x_k;
    x_k = x_k + p;
    % Record all taken steps including iteration index
    nTaken = nTaken + 1;
    info.xs(:,nTaken+1) = x_k;
    info.xind(nTaken+1) = k;
    % Evaluate stopping condition:
    switch stopType
      case 'step'
        % relative step length
        stopCond = (norm(x_k - x_k_1)/norm(x_k_1) < tol);
      case 'grad'
        % gradient norm
        %stopCond = (norm(F.df(x_k)) < tol);
        stopCond = (norm(F.df(x_k), 'inf') < tol*(1 + abs(F.f(x_k))));
    end
  elseif Delta_k < 1e-6*Delta
    % Stop iteration if Delta_k shrank below 1e-6*Delta. Otherwise, if the model
    % does not improve inspite of shinking, the algorithm would shrink Delta_k indefinitely.
    warning('Region of interest is to small. Terminating iteration.')
    break;
  
  end
end

f_k = F.f(x_k);
info.stopCond = stopCond;
info.xs(:,nTaken+2:end) = [];
info.xind(nTaken+2:end) = [];
info.rhos(k+1:end) = [];
info.Deltas(k+1:end) = [];

end





