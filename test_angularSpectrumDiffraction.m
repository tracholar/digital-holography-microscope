N = 512;
obj = im2double(rgb2gray(imread('kuang.jpg')));
z0 = 1e1;
dx = 1.13e-4;
lambda = 650e-9;

Obj = fresnelDiffraction(obj, z0, dx, lambda);
Obj = Obj/max(max(abs(Obj)));
subplot(2,3,1);
imshow(obj);
subplot(2,3,2);
imshow(abs(Obj));
subplot(2,3,3);
imshow(angle(Obj));

Obj2 = angularSpectrumDiffraction(obj, z0, dx, lambda);
Obj2 = Obj2/max(max(abs(Obj2)));
subplot(2,3,4);
imshow(obj);
subplot(2,3,5);
imshow(abs(Obj2));
subplot(2,3,6);
imshow(angle(Obj2));
