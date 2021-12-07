function s = Perlin_Noise_Sparse(x,y,iterations,saturation)
    s = full(sprandn(x,y,0.4));     % Prepare output image (size: m x m)
    s(s<0) = 0; s(s>1) = 1;
    for i=iterations     
        s_new = full(sprandn(ceil(x/i)+i,ceil(y/i)+i,0.4)); s_new(s_new<0) = 0; s_new(s_new>1) = 1;
        d = interp2(s_new, i-1, 'spline');
        s = s + i * d(1:x, 1:y);
    end
    s = (s - min(min(s(:,:)))) ./ (max(max(s(:,:))) - min(min(s(:,:))));
    s = s.*(1+saturation*2);s = s-saturation; s(s<0) = 0; s(s>1) = 1;
end