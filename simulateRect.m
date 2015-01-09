%% ����ȫϢ�������ֵ����
% author: ��Ԫ
% email: zuoyuan@mail.ustc.edu.cn
% ƽ�й��������ؽ�

%% һЩ�����Ķ���
Nx = 512; Ny = 512; 
dx = 0.01e-3; dy = dx;
Lx = Nx*dx; Ly = Ny*dy;
X = 1e-3; Y = X;

lambda = 632.8e-9;
k = 2*pi/lambda;
z0 = max((4*X+2*Lx)/lambda*dx, (4*Y+2*Ly)/lambda*dy);
dx0 = lambda*z0/Nx/dx;
dy0 = lambda*z0/Ny/dy;

alpha = 0.031416;asin((3*X+Lx)/2/z0); % �ο���Ƕ�ƫ�ã���Ӱ�����ͼ����Ƶ�׷ֲ�����������������Թ۲쵽Ƶ�׵ı仯
beta = asin((3*Y+Ly)/2/z0);


%% Ŀ����
obj = im2double(rgb2gray(imread('kuang.jpg')));
subplot(2,3,1);
imshow(obj);

%% �ο���
R = zeros(Nx,Ny);

for i = 1:Nx
    for j = 1:Ny
        R(i,j) = exp(-sqrt(-1)*k*((i-Nx/2)*dx*sin(alpha) + (j-Nx/2)*dy*sin(beta)));
    end
end


%% �õ�CCD���ܵ��ĸ���ͼ��
Obj = fresnelDiffraction(obj, z0, dx0, lambda);
% Obj = angularSpectrumDiffraction(obj, z0, dx0, lambda);

subplot(2,3,2);

Obj = Obj / max(max(abs(Obj))) * 5; % ��ʵ��ֵ��һ������������ʾ��Բο���R�����ǿ��
H = abs(Obj+R).^2;
H = H/max(max(H));
imshow(H);

%% ����ͼ��Ƶ��
subplot(2,3,3);
fftQxI = fftshift(fft2(H));
tmp = abs(fftQxI);
imshow(tmp/max(max(tmp))*20);

%% 
CH = H;
subplot(2,3,4);

animate = 1;
if animate
    for zi = linspace(1e-5, z0, 20)
        constructedObj = fresnelDiffraction(CH, zi, dx, lambda);
        coI = abs((constructedObj));
        imshow(coI/max(max(coI))*5);
        pause(.1);
    end
else
    constructedObj = angularSpectrumDiffraction(CH, z0, dx, lambda);
    coI = abs((constructedObj));
    imshow(coI/max(max(coI)));
end
subplot(2,3,5);
plot(abs(coI(Nx/2,:))/max(max(coI)));
axis([0 512 0 1]);
% figure;
% imshow(cstrObj/max(max(cstrObj)));