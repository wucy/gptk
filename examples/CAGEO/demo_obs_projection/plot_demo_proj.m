% CAGEO - Demonstration of projection onto active set
% as more and more observation are seen. This simply plots 
% the data generated from the C++ example.

% Load data
active_set = csvread('demo_proj_active_set.csv');
psgpmean = csvread('demo_proj_psgp_mean.csv');
psgpvar = csvread('demo_proj_psgp_var.csv');
obs = csvread('demo_proj_obs.csv');
test = csvread('demo_proj_test.csv');

Xactive = active_set(1,:);
Yactive = active_set(2,:);

Xtest = test(1,:);
Ytest = test(2,:);

Xobs = obs(1:7, :);
Yobs = obs(8:end, :);

N = 7;

for i=1:N
  k=2^(i-1);
  
  m = psgpmean(i,:);
  v = psgpvar(i,:);
  
  subplot(2,ceil(N/2),i); hold on;
  error_fill(Xtest', m'-2*sqrt(v'), m'+2*sqrt(v'));
  
  plot(Xtest, Ytest, '-', 'Color', [0.5,0.5,0.5]);
  
  plot(Xobs(i,1:k), Yobs(i,1:k), '+', 'Color', [0.7,0.7,0.7]);
  plot(Xtest, m, '-', 'Color', [0.5,0.5,0.5], 'LineWidth', 2);
  plot(Xactive, Yactive, '.k');
  
  axis([-20.5,20.5,-5,5]);
  
  axis square;
  %set(gca,'ytick',[]);
  %set(gca,'xtick',[]);
  box on;
  
  title([num2str(k) ' observations']);
  
end

legend('PSGP variance', 'true function', 'observations', 'PSGP mean', ...
'active points');
