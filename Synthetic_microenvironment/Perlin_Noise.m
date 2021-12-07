function s = Perlin_Noise(x,y,iterations,saturation)
    s = randn([x,y]);     % Prepare output image (size: m x m)
    for i=iterations   
        if i<5
            d = interp2(randn([ceil(x/i)+i,ceil(y/i)+i]), i-1, 'spline');
        elseif i<7
            d = interp2(randn([ceil(x/(i+3)),ceil(y/(i+3))]), i-1, 'spline');
        elseif i<9
            d = interp2(randn([ceil(x/(i+10)),ceil(y/(i+10))]), i-1, 'spline');
        elseif i<11
            d = interp2(randn([ceil(x/(i+100)),ceil(y/(i+100))]), i-1, 'spline');
        elseif i<20
            d = interp2(randn([ceil(x/(i+400)),ceil(y/(i+400))]), i-1, 'spline');
        end
        s = s + i * d(1:x, 1:y);
    end
    s = (s - min(min(s(:,:)))) ./ (max(max(s(:,:))) - min(min(s(:,:))));
    s = s.*(1+saturation*2);s = s-saturation; s(s<0) = 0; s(s>1) = 1;
end