%% 数字全息成像的数值仿真
% author: 左元
% email: zuoyuan@mail.ustc.edu.cn
% 球面波照明和重建，同轴无透镜傅里叶变换全息

%% 一些常数的定义
Nx = 512; Ny = 512; 
dx = 0.01e-3; dy = dx;
Lx = Nx*dx; Ly = Ny*dy;
X = 1e-3; Y = X;

lambda = 632.8e-9;
k = 2*pi/lambda; m = 1;
z0 = max((X+(1-1/m)*Lx)/lambda*dx, (Y+(1-1/m)*Ly)/lambda*dy); % 同轴无透镜傅里叶变换全系，最小记录距离
zr = m*z0;    
dx0 = lambda*z0/Nx/dx;
dy0 = lambda*z0/Ny/dy;

xr = 3/2*X; % 参考光角度偏置，将影响干涉图样的频谱分布，调整这个参数可以观察到频谱的变化
yr = 3/2*Y;



%% 目标物
obj = im2double(rgb2gray(imread('kuang.jpg')));
subplot(2,3,1);
imshow(obj);

%% 参考光，球面波，菲涅尔近似为二次相因子
R = zeros(Nx,Ny);
for i = 1:Nx
    for j = 1:Ny
        x = (i-Nx/2)*dx;
        y = (j-Ny/2)*dy;
        R(i,j) = exp(sqrt(-1)*k/2/zr*((x)^2 + (y)^2));
    end
end


%% 得到CCD接受到的干涉图样
Obj = fresnelDiffraction(obj, z0, dx0, lambda);
subplot(2,3,2);

Obj = Obj / max(max(abs(Obj))) * 5; % 对实际值归一化，倍乘数表示相对参考光R的相对强度
H = abs(Obj+R).^2;
H = H/max(max(H));
imshow(H);

%% 干涉图样频谱
subplot(2,3,3);
fftQxI = fftshift(fft2(H));
tmp = abs(fftQxI);
imshow(tmp/max(max(tmp))*20);


%% 再现光，离轴球面波，增加相移
R2 = zeros(Nx,Ny);
for i = 1:Nx
    for j = 1:Ny
        x = (i-Nx/2)*dx;
        y = (j-Ny/2)*dy;
        R2(i,j) = exp(-sqrt(-1)*k/2/z0*((x-xr)^2 + (y-yr)^2));
    end
end

%% 再现
CH = R2 .* H;
subplot(2,3,4);

animate = 1;
if animate
    for zi = linspace(1e-5, z0, 40)
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
plot(abs(coI(Nx/2,:))/max(max(coI)));
axis([0 512 0 1]);
% figure;
% imshow(cstrObj/max(max(cstrObj)));