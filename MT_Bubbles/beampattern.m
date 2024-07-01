array = ones([8,8]);

arrayFFT = fftshift(abs(fft2(array, 512,512)));

figure;

surf(arrayFFT, 'LineStyle','none');