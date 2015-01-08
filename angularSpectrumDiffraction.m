%% 利用角谱理论，计算衍射以及全息重建
function Y = angularSpectrumDiffraction(X, z0, d0, lambda)
    fftX = fft2(X);
    Gb = zeros(size(X));
    N = length(X);
    ldf = lambda/N/d0;  % lambda * df
    for i = 1:N
        for j = 1:N
            Gb(i,j) = exp(sqrt(-1)*2*pi/lambda*z0*sqrt(1 - (i-N/2)^2*ldf^2 - (j-N/2)^2*ldf^2));
            
            % 菲涅尔近似
%             Gb(i, j) = exp(sqrt(-1)*pi/lambda*z0 *(2 - (i-N/2)^2*ldf^2 - (j-N/2)^2*ldf^2));   
            
        end
    end
    fftY = fftX .* fftshift(Gb);
    Y = ifft2(fftY);
    
end