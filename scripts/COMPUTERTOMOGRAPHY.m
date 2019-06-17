% Exercise 1

objects = InitObjects();

imshow(PhantomCircle(objects, 255),[0 255]);

p = 300;
q = 100;
sinogram = GenerateMeasuredData(p, q, objects);

s = linspace(-1, 1, 2*q+1);
b = pi;
wb = CalculateFilter(s, b, 'Shepp-Logan');

convolution = CalculateConvolution(sinogram, p, q, wb);
fFBI = CalculateBackprojection(convolution, p, q, 255);
disp(1);
