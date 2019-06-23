% Exercise 1

resultsfigure = figure;
set(resultsfigure,'PaperPositionMode','auto')

objects = InitObjects();
res = 255;

tstart = tic;
phantom = PhantomCircle(objects, res);
fprintf("phantom generation: %0.2e s\n", toc(tstart));

h1 = subplot(2,3,1);
disp(get(h1, 'Position'));
imshow(phantom, [min(min(phantom)) max(max(phantom))]), colorbar;
disp(get(h1, 'Position'));
title('Phantom');

p = 300;
q = 100;

tstart = tic;
sinogram = GenerateMeasuredData(p, q, objects);
fprintf("sinogram generation: %0.2e s\n", toc(tstart));

subplot(2,3,2);
imshow(sinogram,[min(min(sinogram)),max(max(sinogram))]), colorbar;
title('Sinogram');

s = linspace(-1, 1, 2*q+1);
b = 100*pi;
wb = CalculateFilter(s, b, 'Cosine');
% wb = CalculateFilter(s, b, 'no filter');


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
print(resultsfigure, '-dpng', sprintf('../pics/reconstruction_p%d_q%d_b100pi_w_cosine.png', p, q), '-r300');