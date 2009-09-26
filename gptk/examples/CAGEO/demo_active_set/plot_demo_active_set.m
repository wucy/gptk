% CAGEO - Demonstration of projection onto active set
% as more and more observation are seen. This simply plots 
% the data generated from the C++ example.

addpath ..

% Load data
active_set = csvread('demo_set_active_set.csv');
gp = csvread('demo_set_gp.csv');
psgpmean = csvread('demo_set_psgp_mean.csv');
psgpvar = csvread('demo_set_psgp_var.csv');
obs = csvread('demo_set_obs.csv');
test = csvread('demo_set_test.csv');


Xtest = test(1,:);
Ytest = test(2,:);

Xobs = obs(1, :);
Yobs = obs(2, :);

N = 4;

for i=1:N
  k=2^(i-1)*8;   % 8, 16, 32, 64
  
  m = psgpmean(i,:);
  v = psgpvar(i,:);
  
  Xactive = active_set(i,1:k);
  Yactive = active_set(N+i,1:k);
  
  subplot(2,ceil(N/2)+1,i); hold on;
  error_fill(Xtest', m'-2*sqrt(v'), m'+2*sqrt(v'));
  
  plot(Xtest, Ytest, '-', 'Color', [0.5,0.5,0.5]);
  
  plot(Xobs, Yobs, '+', 'Color', [0.7,0.7,0.7]);
  plot(Xtest, m, '-', 'Color', [0.5,0.5,0.5], 'LineWidth', 2);
  plot(Xactive, Yactive, '.k');
  
  axis([-20.5,20.5,-3,3]);
  
  axis square;
  %set(gca,'ytick',[]);
  %set(gca,'xtick',[]);
  box on;
  
  title([num2str(k) ' active points']);
  
end

legend('PSGP/GP variance', 'true function', 'observations', 'PSGP/GP mean', 'active points');

m = gp(1,:);
v = gp(2,:);

subplot(2,ceil(N/2)+1,N+1); hold on;
error_fill(Xtest', m'-2*sqrt(v'), m'+2*sqrt(v'));
plot(Xtest, Ytest, '-', 'Color', [0.5,0.5,0.5]);
plot(Xobs, Yobs, '+', 'Color', [0.7,0.7,0.7]);
plot(Xtest, m, '-', 'Color', [0.5,0.5,0.5], 'LineWidth', 2);
axis([-20.5,20.5,-3,3]);
axis square;
box on;
title('full GP');

