Nx = 512; Ny = 512;
dx = 0.01e-3; dy = dx;
X = 1e-3; Y = X;
NX = 200; NY = 200;
lambda = 632.8e-9;
k = 2*pi/lambda;
z0 = 1e-2;
alpha = 0;
beta = 0.001;
oX = lambda*z0/dx;
obj = zeros(NX, NY);

i = NX/2-10:NX/2+10;
obj(i,NY/2-10) = 1;
obj(i,NY/2+10) = 1;

j = NY/2-10:NY/2+10;
obj(NX/2-10, j) = 1;
obj(NX/2+10, j) = 1;

subplot(1,3,1);
imshow(obj);

EXP = zeros(NX, NY);
R = zeros(NX,NY);

for i = 1:NX
    for j = 1:NY
        EXP(i,j) = exp(sqrt(-1)*2*pi/lambda/2/z0*((i-Nx/2)^2*dx^2+(j-Nx/2)^2*dy^2));
        R(i,j) = exp(-sqrt(-1)*k*(i*sin(alpha) + j*sin(beta)));
    end
end


Obj = EXP .* fft2(obj.*EXP);
subplot(1,3,2);

H = abs(Obj+R).^2;
H = H/max(max(H));
imshow(H);

subplot(1,3,3);
fftQxI = fft2(H);
tmp = log(abs(fftQxI));
imshow(tmp/max(max(tmp)));

CH = conj(R) .* H;
cstrObj = EXP .* fft2(CH .* EXP);
cstrObj = abs(cstrObj).^2;
figure;
imshow(cstrObj/max(max(cstrObj)));