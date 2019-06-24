%% Filtered Backprojection (choose filter and approach)

% filter = 'Ram-Lak';
% filter = 'Cosine';
filter = 'Shepp-Logan';
% filter = 'no filter';

% approach = 'analytic';
approach = 'integral';

resultsfigure = figure;
set(resultsfigure,'PaperPositionMode','auto')

objects = InitObjects();
res = 255;

p = 300;
q = 100;

tstart = tic;
phantom = PhantomCircle(objects, res);
fprintf("phantom generation: %0.2e s\n", toc(tstart));

% img = imread('../pics/lung.jpg');
% img = rgb2gray(img);
% phantom = imresize(img, [255*2+1, 255*2+1]);

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
print(resultsfigure, '-dpng', sprintf('../pics/reconstruction_p%d_q%d_b100pi_w_%s_%s.png', p, q, approach, filter), '-r300');


%% Noise Level Experiment
filter = 'Ram-Lak';
% filter = 'Cosine';

approach = 'analytic';

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

s = linspace(-1, 1, 2*q+1);
b = 100*pi;
wb = CalculateFilter(s, b, filter);
%wb = CalculateFilter(s, b, 'no filter');

noise_levels = [0, 0.01, 0.05, 0.2, 0.5, 1.0];
n_levels = size(noise_levels, 2);

for noise_level_index=1:n_levels

    fprintf('noise level: %.2f \n', noise_levels(noise_level_index));
    sinogram = addNoise(sinogram, noise_levels(noise_level_index));

    subplot(3, n_levels, noise_level_index);
    imshow(sinogram,[min(min(sinogram)),max(max(sinogram))]), colorbar;
    title(sprintf('noise level: %.2f', noise_levels(noise_level_index)));
    
    subplot(3, n_levels, n_levels + noise_level_index);
    convolution = CalculateConvolution(sinogram, p, q, wb);
    imshow(convolution,[min(min(convolution)),max(max(convolution))]), colorbar;

    subplot(3, n_levels, 2*n_levels + noise_level_index);
    fFBI = CalculateBackprojection(convolution, p, q, res);
    imshow(fFBI,[min(min(fFBI)),max(max(fFBI))]), colorbar;
    
end

set(gcf,'position',[100 100 1920 1080])
print(resultsfigure, '-dpng', sprintf('../pics/reconstruction_p%d_q%d_b100pi_w_%s_%s_noise.png', p, q, approach, filter), '-r300');

%% Cut-Off Frequency Experiment
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


%% FBPTest
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
