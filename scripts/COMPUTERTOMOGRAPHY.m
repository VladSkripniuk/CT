%% Filtered Backprojection (choose filter and approach)

filter = 'Ram-Lak';
% filter = 'Cosine';
% filter = 'Shepp-Logan';
% filter = 'no filter';

approach = 'analytic';
% approach = 'integral';

method = 'fb';
%method = 'art';

resultsfigure = figure;
set(resultsfigure,'PaperPositionMode','auto')

objects = InitObjects();
res = 255;
%res = 128;

p = 300;
%p = 150;
q = 100;
%q = 50;

numberofiterations = 10;
lambda = 1;

tstart = tic;
phantom = PhantomCircle(objects, res);
fprintf("phantom generation: %0.2e s\n", toc(tstart));

% img = imread('../pics/lung.jpg');
% img = rgb2gray(img);
% phantom = imresize(img, [255*2+1, 255*2+1]);

tstart = tic;
if strcmp(approach, 'analytic')
    sinogram = GenerateMeasuredData(p, q, objects);
elseif strcmp(approach, 'integral')
    sinogram = GenerateMeasuredDataIntegration(p, q, phantom);
else
    print('invalid data generation approach');
end
fprintf("sinogram generation: %0.2e s\n", toc(tstart));



if strcmp(method, 'art')
    % reshape sinogram into a vector
    g = reshape(sinogram, [size(sinogram, 1)*size(sinogram, 2), 1]);
    
    tstart = tic;
    A = CalculateA(p, q, (2*res+1)*(2*res+1));
    fprintf("calculation of matrix A: %0.2e s\n", toc(tstart));
    
    tstart = tic;
    F = Kaczmarz(A, g, lambda, numberofiterations);
    fprintf("Kaczmarz Method: %0.2e s\n", toc(tstart));

    % reshape result vector into a square image
    fFBI = reshape(F, [2*res+1, 2*res+1]);
else
    s = linspace(-1, 1, 2*q+1);
    b = 100*pi;
    wb = CalculateFilter(s, b, filter);
    %wb = CalculateFilter(s, b, 'no filter');

    tstart = tic;
    convolution = CalculateConvolution(sinogram, p, q, wb);
    fprintf("precompute convolutions: %0.2e s\n", toc(tstart));

    tstart = tic;
    fFBI = CalculateBackprojection(convolution, p, q, res);
    fprintf("FBI: %0.2e s\n", toc(tstart));
end

% fFBI = fFBI * 100;

error = abs(phantom - fFBI);
relative_error = error ./ (phantom + 1e-4);

% Plot the results:
if strcmp(method, 'art')
    subplot(1,3,1);
    imshow(phantom, [min(min(phantom)) max(max(phantom))]), colorbar;
    title('Phantom');
    
    subplot(1,3,2);
    imshow(sinogram,[min(min(sinogram)),max(max(sinogram))]), colorbar;
    title('Sinogram');
    
    subplot(1,3,3);
    imshow(fFBI,[min(min(fFBI)) max(max(fFBI))]), colorbar;
    title('fFBI');
    
    set(gcf,'position',[100 100 1920 1080])
    print(resultsfigure, '-dpng', sprintf('../pics/Kaczmarz_p%d_q%d_b100pi_w_%s_lam%d_it%d.png', p, q, approach, lambda, numberofiterations), '-r300');
else
    h1 = subplot(2,3,1);
    disp(get(h1, 'Position'));
    imshow(phantom, [min(min(phantom)) max(max(phantom))]), colorbar;
    disp(get(h1, 'Position'));
    title('Phantom', 'FontSize',20);

    subplot(2,3,2);
    imshow(sinogram,[min(min(sinogram)),max(max(sinogram))]), colorbar;
    title('Sinogram', 'FontSize',20);

    subplot(2,3,3);
    imshow(convolution,[min(min(convolution)), max(max(convolution))]), colorbar;
    title('Convolution', 'FontSize',20);

    subplot(2,3,4);
    imshow(fFBI,[min(min(fFBI)) max(max(fFBI))]), colorbar;
    title('fFBI reconstruction', 'FontSize',20);

    subplot(2,3,5);
    imshow(error,[min(min(error)) max(max(error))]), colorbar;
    title('error', 'FontSize',20);
    
    subplot(2,3,6);
    imshow(relative_error,[min(min(relative_error)) max(max(relative_error))]), colorbar;
    title('relative error', 'FontSize',20);
    
    set(gcf,'position',[100 100 1920 1080])
    print(resultsfigure, '-dpng', sprintf('../pics/reconstruction_p%d_q%d_b100pi_w_%s_%s_TEST.png', p, q, approach, filter), '-r300');
end
    
%% Noise Level different filters
filter = 'Ram-Lak';

approach = 'analytic';

resultsfigure = figure;
set(resultsfigure,'PaperPositionMode','auto')

objects = InitObjects();
res = 64;

p = 150;
q = 50;

tstart = tic;
phantom = PhantomCircle(objects, res);
fprintf("phantom generation: %0.2e s\n", toc(tstart));

tstart = tic;
if strcmp(approach, 'analytic')
    sinogram = GenerateMeasuredData(p, q, objects);
elseif strcmp(approach, 'integral')
    sinogram = GenerateMeasuredDataIntegration(p, q, phantom);
else
    print('invalid data generation approach');
end
fprintf("sinogram generation: %0.2e s\n", toc(tstart));

s = linspace(-1, 1, 2*q+1);
b = q*pi;
wb = CalculateFilter(s, b, "Ram-Lak");

noise_levels = [0, 0.01, 0.05, 0.2];
n_levels = size(noise_levels, 2);

sinogram_clean = sinogram;

for noise_level_index=1:n_levels

    sinogram = addNoise(sinogram_clean, noise_levels(noise_level_index));
    convolution = CalculateConvolution(sinogram, p, q, wb);

    subplot(3, n_levels, noise_level_index);
    fFBI = CalculateBackprojection(convolution, p, q, res);
    imshow(fFBI,[min(min(fFBI)),max(max(fFBI))]), colorbar;
    title(sprintf('noise: %.0f %%', 100*noise_levels(noise_level_index)), 'FontSize',20);

    
    
end

wb = CalculateFilter(s, b, "Cosine");

for noise_level_index=1:n_levels
    
    sinogram = addNoise(sinogram_clean, noise_levels(noise_level_index));
    convolution = CalculateConvolution(sinogram, p, q, wb);
    
    subplot(3, n_levels, n_levels + noise_level_index);
    fFBI = CalculateBackprojection(convolution, p, q, res);
    imshow(fFBI,[min(min(fFBI)),max(max(fFBI))]), colorbar;
%     title(sprintf('noise: %.0f %%', 100*noise_levels(noise_level_index)), 'FontSize',20);
end


wb = CalculateFilter(s, b, "Shepp-Logan");

for noise_level_index=1:n_levels

    sinogram = addNoise(sinogram_clean, noise_levels(noise_level_index));
    convolution = CalculateConvolution(sinogram, p, q, wb);
    
    subplot(3, n_levels, 2*n_levels + noise_level_index);
    fFBI = CalculateBackprojection(convolution, p, q, res);
    imshow(fFBI,[min(min(fFBI)),max(max(fFBI))]), colorbar;
%     title(sprintf('noise: %.0f %%', 100*noise_levels(noise_level_index)), 'FontSize',20);
 
    
end

set(gcf,'position',[100 100 1920 1080])
print(resultsfigure, '-dpng', sprintf('../pics/filters_noise.png'), '-r300');


%% Noise Level FBP vs ART
approach = 'analytic';

resultsfigure = figure;
set(resultsfigure,'PaperPositionMode','auto')

objects = InitObjects();
res = 64;

p = 150;
q = 50;

tstart = tic;
phantom = PhantomCircle(objects, res);
fprintf("phantom generation: %0.2e s\n", toc(tstart));

tstart = tic;
if strcmp(approach, 'analytic')
    sinogram = GenerateMeasuredData(p, q, objects);
elseif strcmp(approach, 'integral')
    sinogram = GenerateMeasuredDataIntegration(p, q, phantom);
else
    print('invalid data generation approach');
end
fprintf("sinogram generation: %0.2e s\n", toc(tstart));

s = linspace(-1, 1, 2*q+1);
b = q*pi;
wb = CalculateFilter(s, b, "Cosine");

noise_levels = [0, 0.01, 0.05, 0.2];
n_levels = size(noise_levels, 2);

sinogram_clean = sinogram;

for noise_level_index=1:n_levels

    sinogram = addNoise(sinogram_clean, noise_levels(noise_level_index));
    convolution = CalculateConvolution(sinogram, p, q, wb);

    subplot(2, n_levels, noise_level_index);
    fFBI = CalculateBackprojection(convolution, p, q, res);
    imshow(fFBI,[min(min(fFBI)),max(max(fFBI))]), colorbar;
    title(sprintf('noise: %.0f %%', 100*noise_levels(noise_level_index)), 'FontSize',20);

    
    
end


A = CalculateA(p, q, M);

for noise_level_index=1:n_levels

    sinogram = addNoise(sinogram_clean, noise_levels(noise_level_index));
    g = reshape(sinogram', p * (2*q+1), 1);

    num_iter = 5;
    lambda_schedule = 2 * ones(1, num_iter) ./ sqrt(1:num_iter);
    
    tstart = tic;
    F = Kaczmarz(A, g, lambda_schedule, num_iter, "dont-saveplots", "randomized", "dont-save");
    fprintf("Kaczmarz wall time: %0.2e s\n", toc(tstart));

    subplot(2, n_levels, n_levels + noise_level_index);
    imshow(reshape(F, sqrt(M), sqrt(M)),[min(min(F)),max(max(F))]), colorbar;        
    
end


set(gcf,'position',[100 100 1920 1080])
print(resultsfigure, '-dpng', sprintf('../pics/ARTvsFBP_noise.png'), '-r300');



%% Noise Level Experiment (Ex 2c))
filter = 'Ram-Lak';
% filter = 'Cosine';

approach = 'analytic';

resultsfigure = figure;
set(resultsfigure,'PaperPositionMode','auto')

objects = InitObjects();
res = 64;

p = 150;
q = 50;

tstart = tic;
phantom = PhantomCircle(objects, res);
fprintf("phantom generation: %0.2e s\n", toc(tstart));

tstart = tic;
if strcmp(approach, 'analytic')
    sinogram = GenerateMeasuredData(p, q, objects);
elseif strcmp(approach, 'integral')
    sinogram = GenerateMeasuredDataIntegration(p, q, phantom);
else
    print('invalid data generation approach');
end
fprintf("sinogram generation: %0.2e s\n", toc(tstart));

s = linspace(-1, 1, 2*q+1);
b = q*pi;
wb = CalculateFilter(s, b, filter);
%wb = CalculateFilter(s, b, 'no filter');

noise_levels = [0, 0.01, 0.05, 0.2, 0.5, 1.0];
n_levels = size(noise_levels, 2);

for noise_level_index=1:n_levels

    fprintf('noise level: %.2f \n', noise_levels(noise_level_index));
    sinogram = addNoise(sinogram, noise_levels(noise_level_index));

    subplot(3, n_levels, noise_level_index);
    imshow(sinogram,[min(min(sinogram)),max(max(sinogram))]), colorbar;
    title(sprintf('noise: %.0f %%', 100*noise_levels(noise_level_index)), 'FontSize',20);
    
    subplot(3, n_levels, n_levels + noise_level_index);
    convolution = CalculateConvolution(sinogram, p, q, wb);
    imshow(convolution,[min(min(convolution)),max(max(convolution))]), colorbar;

    subplot(3, n_levels, 2*n_levels + noise_level_index);
    fFBI = CalculateBackprojection(convolution, p, q, res);
    imshow(fFBI,[min(min(fFBI)),max(max(fFBI))]), colorbar;
    
end

set(gcf,'position',[100 100 1920 1080])
print(resultsfigure, '-dpng', sprintf('../pics/reconstruction_p%d_q%d_b100pi_w_%s_%s_noise.png', p, q, approach, filter), '-r300');

%% Cut-Off Frequency Experiment (Ex 2d))
filter = 'Shepp-Logan';
approach = 'analytic';
noise_level = 0.2;

resultsfigure = figure;
set(resultsfigure,'PaperPositionMode','auto')

objects = InitObjects();
res = 255;

p = 300;
q = 100;

tstart = tic;
phantom = PhantomCircle(objects, res);
fprintf("phantom generation: %0.2e s\n", toc(tstart));

tstart = tic;
if strcmp(approach, 'analytic')
    sinogram = GenerateMeasuredData(p, q, objects);
elseif strcmp(approach, 'integral')
    sinogram = GenerateMeasuredDataIntegration(p, q, phantom);
else
    print('invalid data generation approach');
end
fprintf("sinogram generation: %0.2e s\n", toc(tstart));

fprintf('noise level: %.2f \n', noise_level);
sinogram = addNoise(sinogram, noise_level);

s = linspace(-1, 1, 2*q+1);

%wb = CalculateFilter(s, b, 'no filter');

b_values = [q, q/2, q/4, q/8, q/16, q/32];
n_b = size(b_values, 2);

for b_value_index=1:n_b

    b = b_values(b_value_index)*pi;
    wb = CalculateFilter(s, b, filter);
    
    fprintf('b value: %.2f*pi \n', b_values(b_value_index));
    subplot(3, n_b, b_value_index);
    imshow(sinogram,[min(min(sinogram)),max(max(sinogram))]), colorbar;
    title(sprintf('b value: %.2f*pi', b_values(b_value_index)));
    
    subplot(3, n_b, n_b + b_value_index);
    convolution = CalculateConvolution(sinogram, p, q, wb);
    imshow(convolution,[min(min(convolution)),max(max(convolution))]), colorbar;
    
    subplot(3, n_b, 2*n_b + b_value_index);
    fFBI = CalculateBackprojection(convolution, p, q, res);
    imshow(fFBI,[min(min(fFBI)),max(max(fFBI))]), colorbar;
    
end

set(gcf,'position',[100 100 1920 1080])
print(resultsfigure, '-dpng', sprintf('../pics/reconstruction_p%d_q%d_w_%s_%s_blevels.png', p, q, approach, filter), '-r300');


%% FBPTest (Ex 2e))
% filter = 'Ram-Lak';
% filter = 'Cosine';
filter = 'Shepp-Logan';
% filter = 'no filter';

approach = 'integral';

resultsfigure = figure;
set(resultsfigure,'PaperPositionMode','auto')

objects = InitObjects();
res = 128;

p = 180;
q = 100;

tstart = tic;
phantom_struct = load('../pics/FBPTest.mat');
phantom = phantom_struct.phantom;
fprintf("phantom generation: %0.2e s\n", toc(tstart));

h1 = subplot(2,3,1);
disp(get(h1, 'Position'));
imshow(phantom, [min(min(phantom)) max(max(phantom))]), colorbar;
disp(get(h1, 'Position'));
title('Phantom');

tstart = tic;
if strcmp(approach, 'analytic')
    sinogram = GenerateMeasuredData(p, q, objects);
elseif strcmp(approach, 'integral')
    sinogram = GenerateMeasuredDataIntegration(p, q, phantom);
else
    print('invalid data generation approach');
end
fprintf("sinogram generation: %0.2e s\n", toc(tstart));

subplot(2,3,2);
imshow(sinogram,[min(min(sinogram)),max(max(sinogram))]), colorbar;
title('Sinogram');

s = linspace(-1, 1, 2*q+1);
b = 100*pi;
wb = CalculateFilter(s, b, filter);
%wb = CalculateFilter(s, b, 'no filter');

tstart = tic;
convolution = CalculateConvolution(sinogram, p, q, wb);
fprintf("precompute convolutions: %0.2e s\n", toc(tstart));

subplot(2,3,3);
imshow(convolution,[min(min(convolution)), max(max(convolution))]), colorbar;
title('Convolution');

tstart = tic;
fFBI = CalculateBackprojection(convolution, p, q, res);
fprintf("FBI: %0.2e s\n", toc(tstart));


subplot(2,3,4);
imshow(fFBI,[min(min(fFBI)) max(max(fFBI))]), colorbar;
title('fFBI');

error = abs(phantom - fFBI);
subplot(2,3,5);
imshow(error,[min(min(error)) max(max(error))]), colorbar;
title('error');

set(gcf,'position',[100 100 1920 1080])
print(resultsfigure, '-dpng', sprintf('../pics/reconstruction_FBPTest_p%d_q%d_b100pi_w_%s_%s.png', p, q, approach, filter), '-r300');


%% Sinogram Comparison (Ex 2f))

resultsfigure = figure;
set(resultsfigure,'PaperPositionMode','auto')

objects = InitObjects();
%res = 255;
res = 255;

p = 600;
%p = 150;
q = 200;
%q = 50;

tstart = tic;
phantom = PhantomCircle(objects, res);
fprintf("phantom generation: %0.2e s\n", toc(tstart));

tstart = tic;
sinogram_analytic = GenerateMeasuredData(p, q, objects);
analytic_time = toc(tstart);
fprintf("sinogram generation (analytic): %0.2e s\n", analytic_time);
tstart = tic;
sinogram_integration = GenerateMeasuredDataIntegration(p, q, phantom);
integration_time = toc(tstart);
fprintf("sinogram generation (integration): %0.2e s\n", integration_time);

%Plot Sinograms only
subplot(1, 2, 1);
imshow(sinogram_analytic,[min(min(sinogram_analytic)),max(max(sinogram_analytic))]), colorbar;
title(sprintf('Analytic Calculation, time: %f', analytic_time));
subplot(1, 2, 2);
imshow(sinogram_integration,[min(min(sinogram_integration)),max(max(sinogram_integration))]), colorbar;
title(sprintf('Calculation via Integration, time: %f', integration_time));
set(gcf,'position',[100 100 1920 1080])
print(resultsfigure, '-dpng', sprintf('../pics/sinogram_comparison_p%d_q%d_res%d.png', p, q, res), '-r300');


%Include Plot of Phantom
% subplot(1, 3, 1);
% imshow(phantom,[min(min(phantom)),max(max(phantom))]), colorbar;
% title(sprintf('Shepp-Logan Phantom, Resolution: %d', res));
% subplot(1, 3, 2);
% imshow(sinogram_analytic,[min(min(sinogram_analytic)),max(max(sinogram_analytic))]), colorbar;
% title(sprintf('Analytic Calculation, time: %f', analytic_time));
% subplot(1, 3, 3);
% imshow(sinogram_integration,[min(min(sinogram_integration)),max(max(sinogram_integration))]), colorbar;
% title(sprintf('Calculation via Integration, time: %f', integration_time));
% set(gcf,'position',[100 100 1920 1080])
% print(resultsfigure, '-dpng', sprintf('../pics/sinogram_comparison_p%d_q%d_res%d_phantom.png', p, q, res), '-r300');

%% Limited Angles (Ex 2g))

%filter = 'Ram-Lak';
% filter = 'Cosine';
filter = 'Shepp-Logan';

approach = 'analytic';

resultsfigure = figure;
set(resultsfigure,'PaperPositionMode','auto')

objects = InitObjects();
res = 64;

p = 150;
q = 50;

% tstart = tic;
% phantom = PhantomCircle(objects, res);
% fprintf("phantom generation: %0.2e s\n", toc(tstart));

s = linspace(-1, 1, 2*q+1);
b = q*pi;
wb = CalculateFilter(s, b, filter);
%wb = CalculateFilter(s, b, 'no filter');

angle_restr = [0, pi/9, 2*pi/9, 3*pi/9, 4*pi/9];

prettystrings = ["$$\phi \in [0, \pi]$$";
    "$$\phi \in \left[\frac{\pi}{9}, \frac{8\pi}{9}\right]$$";
    "$$\phi \in \left[\frac{2\pi}{9}, \frac{7\pi}{9}\right]$$"; 
    "$$\phi \in \left[\frac{3\pi}{9}, \frac{6\pi}{9}\right]$$"; 
    "$$\phi \in \left[\frac{4\pi}{9}, \frac{5\pi}{9}\right]$$"];
n_levels = size(angle_restr, 2);

for angle_index=1:n_levels
    
    tstart = tic;

    sinogram = GenerateMeasuredData(p, q, objects, angle_restr(angle_index));

    fprintf("sinogram generation: %0.2e s\n", toc(tstart));

    subplot(3, n_levels, angle_index);
    imshow(sinogram,[min(min(sinogram)),max(max(sinogram))]), colorbar;
    title(prettystrings(angle_index), 'interpreter', 'latex', 'FontSize',20);
    
    subplot(3, n_levels, n_levels + angle_index);
    convolution = CalculateConvolution(sinogram, p, q, wb);
    imshow(convolution,[min(min(convolution)),max(max(convolution))]), colorbar;

    subplot(3, n_levels, 2*n_levels + angle_index);
    fFBI = CalculateBackprojection(convolution, p, q, res, angle_restr(angle_index));
    imshow(fFBI,[min(min(fFBI)),max(max(fFBI))]), colorbar;
    
end

set(gcf,'position',[100 100 1920 1080])
print(resultsfigure, '-dpng', sprintf('../pics/reconstruction_p%d_q%d_b100pi_w_%s_%s_limangles.png', p, q, approach, filter), '-r300');


%% Limited Angles ART

resultsfigure = figure;
set(resultsfigure,'PaperPositionMode','auto')

objects = InitObjects();
res = 64;
M = (2*res+1)^2;

p = 150;
q = 50;

prettystrings = ["$$\phi \in [0, \pi]$$";
    "$$\phi \in \left[\frac{\pi}{9}, \frac{8\pi}{9}\right]$$";
    "$$\phi \in \left[\frac{2\pi}{9}, \frac{7\pi}{9}\right]$$"; 
    "$$\phi \in \left[\frac{3\pi}{9}, \frac{6\pi}{9}\right]$$"; 
    "$$\phi \in \left[\frac{4\pi}{9}, \frac{5\pi}{9}\right]$$"];
angle_restr = [0, pi/9, 2*pi/9, 3*pi/9, 4*pi/9];
n_levels = size(angle_restr, 2);

for angle_index=1:n_levels
    
    sinogram = GenerateMeasuredData(p, q, objects, angle_restr(angle_index));
    g = reshape(sinogram', p * (2*q+1), 1);
    
    subplot(1, n_levels, angle_index);
    
    A = CalculateA(p, q, M, angle_restr(angle_index));

    num_iter = 5;
    lambda_schedule = 2 * ones(1, num_iter) ./ sqrt(1:num_iter);
    
    tstart = tic;
    F = Kaczmarz(A, g, lambda_schedule, num_iter, "dont-saveplots", "randomized", "dont-save");
    fprintf("Kaczmarz wall time: %0.2e s\n", toc(tstart));

    subplot(1, n_levels, angle_index);
    imshow(reshape(F, sqrt(M), sqrt(M)),[min(min(F)),max(max(F))]), colorbar;
    title(prettystrings(angle_index), 'interpreter', 'latex', 'FontSize',20);
    
end

set(gcf,'position',[100 100 1920 1080])
print(resultsfigure, '-dpng', sprintf('../pics/ART_limangles.png'), '-r300');


%% Sheet 3, Ex 4 a)

objects = InitObjects();

res = 127;
M = (2*res+1)^2;

p = 150;
q = 50;

sinogram = GenerateMeasuredData(p, q,objects);
g = reshape(sinogram', p * (2*q+1), 1);

tstart = tic;
A = CalculateA(p, q, M);
fprintf("Matrix A calculation time: %0.2e s\n", toc(tstart));

tstart = tic;
F = Kaczmarz(A, g, 1, 1, "saveplots", "not-randomized", "s3ex4a");
% F = Kaczmarz_highres(M, p, q, g, 1, 1, "saveplots", "randomized");
fprintf("Kaczmarz 1 iteration time: %0.2e s\n", toc(tstart));

%% Sheet 3, Ex 4 b) const lambda
objects = InitObjects();

res = 64;
M = (2*res+1)^2;

p = 150;
q = 50;

sinogram = GenerateMeasuredData(p, q, objects);
g = reshape(sinogram', p * (2*q+1), 1);

A = CalculateA(p, q, M);

num_iter = 15;

lambda_schedule = ones(num_iter, 1);
tstart = tic;
F = Kaczmarz(A, g, lambda_schedule, num_iter, "saveplots", "not-randomized", "s3ex4b");
fprintf("Kaczmarz wall time: %0.2e s\n", toc(tstart));

%% Sheet 3, Ex 4 b) decay lambda
objects = InitObjects();

res = 64;
M = (2*res+1)^2;

p = 150;
q = 50;

sinogram = GenerateMeasuredData(p, q, objects);
g = reshape(sinogram', p * (2*q+1), 1);

A = CalculateA(p, q, M);

num_iter = 15;

lambda_schedule = 2 * ones(1, num_iter) ./ sqrt(1:num_iter);
tstart = tic;
F = Kaczmarz(A, g, lambda_schedule, num_iter, "saveplots", "not-randomized", "s3ex4bdecay");
fprintf("Kaczmarz wall time: %0.2e s\n", toc(tstart));


%% Sheet 3, Ex 4 d)
objects = InitObjects();

res = 64;
M = (2*res+1)^2;

% p = 90;
% q = 30;
p = 300;
q = 100;


sinogram = GenerateMeasuredData(p, q, objects);
g = reshape(sinogram', p * (2*q+1), 1);

A = CalculateA(p, q, M);

num_iter = 15;

lambda_schedule = 2 * ones(1, num_iter) ./ sqrt(1:num_iter);
tstart = tic;
F = Kaczmarz(A, g, lambda_schedule, num_iter, "saveplots", "not-randomized", "s3ex4dp300q100");
fprintf("Kaczmarz wall time: %0.2e s\n", toc(tstart));

%% Sheet 3, Ex 4 randomized
objects = InitObjects();

res = 64;
M = (2*res+1)^2;

p = 150;
q = 50;


sinogram = GenerateMeasuredData(p, q, objects);
g = reshape(sinogram', p * (2*q+1), 1);

A = CalculateA(p, q, M);

num_iter = 15;

lambda_schedule = 2 * ones(1, num_iter) ./ sqrt(1:num_iter);
tstart = tic;
F = Kaczmarz(A, g, lambda_schedule, num_iter, "saveplots", "randomized", "s3ex4rand");
fprintf("Kaczmarz wall time: %0.2e s\n", toc(tstart));

%% Sheet 3, Ex 4 randomized over- and underdetermined systems
objects = InitObjects();

res = 64;
M = (2*res+1)^2;

p = 300;
q = 100;


sinogram = GenerateMeasuredData(p, q, objects);
g = reshape(sinogram', p * (2*q+1), 1);

A = CalculateA(p, q, M);

num_iter = 15;

lambda_schedule = 2 * ones(1, num_iter) ./ sqrt(1:num_iter);
tstart = tic;
F = Kaczmarz(A, g, lambda_schedule, num_iter, "saveplots", "randomized", "s3ex4randp300q100");
fprintf("Kaczmarz wall time: %0.2e s\n", toc(tstart));


%% Sheet 3, Ex 4 FBP vs ART with smaller p q
resultsfigure = figure;
set(resultsfigure,'PaperPositionMode','auto')

objects = InitObjects();

res = 64;
M = (2*res+1)^2;


q_list = [10 20 30 50];

for i=1:4
    q = q_list(i);
    p = 3*q;

    s = linspace(-1, 1, 2*q+1);
    b = 50*pi;
    %b = q*pi;
    wb = CalculateFilter(s, b, "Cosine");

    
    sinogram = GenerateMeasuredData(p, q, objects);
    
    convolution = CalculateConvolution(sinogram, p, q, wb);
    fFBI = CalculateBackprojection(convolution, p, q, res);
    
    subplot(2, 4, i);
    imshow(fFBI,[min(min(fFBI)),max(max(fFBI))]), colorbar;
    title(sprintf("p=%d, q=%d", p, q), 'FontSize',20);
    
%     g = reshape(sinogram', p * (2*q+1), 1);    
%     A = CalculateA(p, q, M);
% 
%     num_iter = 15;
% 
%     lambda_schedule = 2 * ones(1, num_iter) ./ sqrt(1:num_iter);
%     tstart = tic;
%     F = Kaczmarz(A, g, lambda_schedule, num_iter, "dont-save", "randomized", "dont-save");
%     fprintf("Kaczmarz wall time: %0.2e s\n", toc(tstart));
%     
%     subplot(2, 4, i+4);
%     imshow(reshape(F, sqrt(M), sqrt(M)),[min(min(F)),max(max(F))]), colorbar;
%     
    
end


%% Sheet 3, Ex 4 limit angle
objects = InitObjects();

res = 64;
M = (2*res+1)^2;

p = 150;
q = 50;

angle_restr = 4*pi/9;

sinogram = GenerateMeasuredData(p, q, objects,angle_restr);
g = reshape(sinogram', p * (2*q+1), 1);

A = CalculateA(p, q, M, angle_restr);

num_iter = 15;

lambda_schedule = 2 * ones(1, num_iter) ./ sqrt(1:num_iter);
tstart = tic;
F = Kaczmarz(A, g, lambda_schedule, num_iter, "saveplots", "randomized", "s3ex44pi9");
fprintf("Kaczmarz wall time: %0.2e s\n", toc(tstart));


%% Sheet 3, Ex 4 high resolution
objects = InitObjects();

res = 128;
M = (2*res+1)^2;
p = 180;
q = 100;

angle_restr = 0;

sinogram = GenerateMeasuredDataIntegration(p, q, phantom);
g = reshape(sinogram', p * (2*q+1), 1);

% A = CalculateA(p, q, M, angle_restr);

num_iter = 1;

lambda_schedule = 2 * ones(1, num_iter) ./ sqrt(1:num_iter);
tstart = tic;
F = Kaczmarz_highres(M, p, q, g, lambda_schedule, num_iter, "saveplots", "randomized", "s3ex4res255");
fprintf("Kaczmarz wall time: %0.2e s\n", toc(tstart));


%% Sheet 2, Ex 3d FBP cut off freq
resultsfigure = figure;
set(resultsfigure,'PaperPositionMode','auto')

objects = InitObjects();

res = 64;
M = (2*res+1)^2;

prettystrings = ["$$b = 10 \pi$$";
    "$$b = 20 \pi$$";
    "$$b = 50 \pi$$";
    "$$b = 100 \pi$$"];

b_list = [10*pi, 20*pi, 50*pi, 100*pi];

sinogram = GenerateMeasuredData(p, q, objects);
sinogram = addNoise(sinogram, 0.2);

for i=1:4
    q = 50;
    p = 150;

    s = linspace(-1, 1, 2*q+1);
    b = b_list(i);
    wb = CalculateFilter(s, b, "Cosine");
    
    convolution = CalculateConvolution(sinogram, p, q, wb);
    fFBI = CalculateBackprojection(convolution, p, q, res);
    
    subplot(1, 4, i);
    imshow(fFBI,[min(min(fFBI)),max(max(fFBI))]), colorbar;
    title(prettystrings(i), 'interpreter', 'latex', 'FontSize',20);
    
end

set(gcf,'position',[100 100 1920 1080])
print(resultsfigure, '-dpng', sprintf('../pics/cutoffs.png'), '-r300');

%% Sheet 2, Ex 2 e
resultsfigure = figure;
set(resultsfigure,'PaperPositionMode','auto')

objects = InitObjects();

res = 128;
M = (2*res+1)^2;

p = 180;
q = 100;


s = linspace(-1, 1, 2*q+1);
b = q*pi;
wb = CalculateFilter(s, b, "Cosine");

phantom = load('FBPTest.mat');
phantom = phantom.phantom;

% phantom = PhantomCircle(objects, res);


sinogram = GenerateMeasuredDataIntegration(p, q, phantom);

tstart = tic;
convolution = CalculateConvolution(sinogram, p, q, wb);
fFBI = CalculateBackprojection(convolution, p, q, res);
disp(toc(tstart))

subplot(1, 2, 1);
imshow(fFBI,[min(min(fFBI)),max(max(fFBI))]), colorbar;
% title(sprintf("p=%d, q=%d", p, q), 'FontSize',20);

tstart = tic;
g = reshape(sinogram', p * (2*q+1), 1);    
% A = CalculateA(p, q, M);

num_iter = 15;

lambda_schedule = 2 * ones(1, num_iter) ./ sqrt(1:num_iter);

F = Kaczmarz_highres(M, p, q, g, lambda_schedule, num_iter, "dontsave", "randomized", "dontsave");
fprintf("Kaczmarz wall time: %0.2e s\n", toc(tstart));

subplot(1, 2, 2);
imshow(reshape(F, sqrt(M), sqrt(M)),[min(min(F)) max(max(F))]), colorbar;

set(gcf,'position',[100 100 1920 1080])
print(resultsfigure, '-dpng', '../pics/bird.png', '-r300');
