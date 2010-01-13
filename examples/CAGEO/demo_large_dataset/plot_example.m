% Load observations
obs = load('INTAMP_scenario6_1080.csv');

subplot(1,3,1);
scatter(obs(:,1), obs(:,2), 2, obs(:,3));
colorbar
title('Observations');
axis image;

% Load gridded data
pred = csvread('demo_scenario6_psgp_grid.csv');  

% Grid size
ngridx = size(pred,2)/100;     
ngridy = 100;

x = reshape( pred(1, :), ngridx, ngridy );
y = reshape( pred(2, :), ngridx, ngridy );
m = reshape( pred(3, :), ngridx, ngridy );
v = reshape( pred(4, :), ngridx, ngridy );

xrange = x(1,:);
yrange = y(:,1);

subplot(1,3,2);
imagesc(xrange, yrange, m);
colorbar
title('PSGP mean');
axis image;

subplot(1,3,3);
imagesc(xrange, yrange, v);
colorbar
title('PSGP var');
axis image;
