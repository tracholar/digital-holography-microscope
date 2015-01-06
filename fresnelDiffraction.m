%% Simulate fresnell diffraction

function Y = fresnelDiffraction(X,z0,d0,lambda)
    N = length(X);
    if size(X,2)~=N
        error('X is not square');
    end
    
    k = 2*pi/lambda;
    
    dx0 = d0; dy0 = d0;
    dx = lambda*z0/N/d0;    %源目标的采样到接受平面采样间隔的换算
    dy = dx;
    
    Y = zeros(size(X));
    EXP0 = zeros(size(X));
    EXP = zeros(size(X));
    for i = 1:N
        for j = 1:N
            EXP0(i,j) = exp(sqrt(-1)*k/2/z0*((i-N/2)^2*dx0^2+(j-N/2)^2*dy0^2));
            EXP(i,j) = exp(sqrt(-1)*k/2/z0*((i-N/2)^2*dx^2+(j-N/2)^2*dy^2));
        end
    end
    Y(:) = EXP .* fftshift(fft2(X .* EXP0));
end