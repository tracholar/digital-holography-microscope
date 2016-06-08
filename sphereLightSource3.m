%% ����ȫϢ�������ֵ����
% author: ��Ԫ
% email: zuoyuan@mail.ustc.edu.cn
% ���沨�������ؽ���������͸������Ҷ�任ȫϢ����һ��ͼƬ

%% һЩ�����Ķ���
Nx = 512; Ny = 512; 
dx = 0.01e-3; dy = dx;
Lx = Nx*dx; Ly = Ny*dy;
X = dx*400; Y = X;

lambda = 632.8e-9;
k = 2*pi/lambda;
z0 = max(4*X/lambda*dx, 4*Y/lambda*dy); % ������͸������Ҷ�任ȫϵ����С��¼����
zr = z0;    
dx0 = lambda*z0/Nx/dx;
dy0 = lambda*z0/Ny/dy;

xr = 3/2*X; % �ο���Ƕ�ƫ�ã���Ӱ�����ͼ����Ƶ�׷ֲ�����������������Թ۲쵽Ƶ�׵ı仯
yr = 3/2*Y;



%% Ŀ����
obj = im2double(rgb2gray(imread('zkd.jpg')));
subplot(2,3,1);
imshow(obj);

%% �ο��⣬���沨������������Ϊ����������
R = zeros(Nx,Ny);
for i = 1:Nx
    for j = 1:Ny
        x = (i-Nx/2)*dx;
        y = (j-Ny/2)*dy;
        R(i,j) = exp(sqrt(-1)*k/2/zr*((x-xr)^2 + (y-yr)^2));
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


%% ���ֹ⣬ͬ�����沨
R2 = zeros(Nx,Ny);
for i = 1:Nx
    for j = 1:Ny
        x = (i-Nx/2)*dx;
        y = (j-Ny/2)*dy;
        R2(i,j) = exp(-sqrt(-1)*k/2/zr*((x)^2 + (y)^2));
    end
end

%% ����
CH = R2 .* H;
subplot(2,3,4);

animate = 0;
if animate
    for zi = linspace(1e-5, z0, 20)
        constructedObj = fresnelDiffraction(CH, zi, dx, lambda);
        coI = abs((constructedObj));
        imshow(coI/max(max(coI))*20);
        pause(.1);
    end
else
    constructedObj = fresnelDiffraction(CH, z0, dx, lambda);
    coI = abs((constructedObj));
    imshow(coI/max(max(coI))*20);
end
subplot(2,3,5);
plot(abs(coI(Nx/2,:))/max(coI(Nx/2,:)));
axis([0 512 0 1]);
% figure;
% imshow(cstrObj/max(max(cstrObj)));