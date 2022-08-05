function imRGB = to_rgb(im)

im_s = size(im);
if im_s(3)==6;
    imRGB = zeros([im_s(1),im_s(2),3]);
    imRGB(:,:,3) = imRGB(:,:,3) + im(:,:,1); % Blue
    imRGB(:,:,2) = imRGB(:,:,2) + im(:,:,2); % Green
    imRGB(:,:,1) = imRGB(:,:,1) + im(:,:,3); % Red
    imRGB(:,:,1) = imRGB(:,:,1) + im(:,:,4); imRGB(:,:,2) = imRGB(:,:,2) + im(:,:,4); % Yellow    
    imRGB(:,:,1) = imRGB(:,:,1) + im(:,:,5); imRGB(:,:,3) = imRGB(:,:,3) + im(:,:,5); % Magenta    
    imRGB(:,:,2) = imRGB(:,:,2) + im(:,:,6); imRGB(:,:,3) = imRGB(:,:,3) + im(:,:,6); % Cyan   
    figure; imshow(imRGB,[])

    for ii=1:size(im,3); im(:,:,ii)=im(:,:,ii)./max(max(im(:,:,ii))); end;
    im_s = size(im);
    imRGB = zeros([im_s(1),im_s(2),3]);
    imRGB(:,:,3) = imRGB(:,:,3) + im(:,:,1); % Blue
    imRGB(:,:,2) = imRGB(:,:,2) + im(:,:,2); % Green
    imRGB(:,:,1) = imRGB(:,:,1) + im(:,:,3); % Red
    imRGB(:,:,1) = imRGB(:,:,1) + im(:,:,4); imRGB(:,:,2) = imRGB(:,:,2) + im(:,:,4); % Yellow    
    imRGB(:,:,1) = imRGB(:,:,1) + im(:,:,5); imRGB(:,:,3) = imRGB(:,:,3) + im(:,:,5); % Magenta    
    imRGB(:,:,2) = imRGB(:,:,2) + im(:,:,6); imRGB(:,:,3) = imRGB(:,:,3) + im(:,:,6); % Cyan   
    figure; imshow(imRGB,[])
    
else if im_s(3)==5;       
    imRGB = zeros([im_s(1),im_s(2),3]);
    imRGB(:,:,3) = imRGB(:,:,3) + im(:,:,1); % Blue
    imRGB(:,:,2) = imRGB(:,:,2) + im(:,:,2); % Green
    imRGB(:,:,1) = imRGB(:,:,1) + im(:,:,3); % Red
    imRGB(:,:,1) = imRGB(:,:,1) + im(:,:,4); imRGB(:,:,2) = imRGB(:,:,2) + im(:,:,4); % Yellow    
    imRGB(:,:,1) = imRGB(:,:,1) + im(:,:,5); imRGB(:,:,3) = imRGB(:,:,3) + im(:,:,5); % Magenta    
    figure; imshow(imRGB,[])

    for ii=1:size(im,3); im(:,:,ii)=im(:,:,ii)./max(max(im(:,:,ii))); end;
    im_s = size(im);
    imRGB = zeros([im_s(1),im_s(2),3]);
    imRGB(:,:,3) = imRGB(:,:,3) + im(:,:,1); % Blue
    imRGB(:,:,2) = imRGB(:,:,2) + im(:,:,2); % Green
    imRGB(:,:,1) = imRGB(:,:,1) + im(:,:,3); % Red
    imRGB(:,:,1) = imRGB(:,:,1) + im(:,:,4); imRGB(:,:,2) = imRGB(:,:,2) + im(:,:,4); % Yellow    
    imRGB(:,:,1) = imRGB(:,:,1) + im(:,:,5); imRGB(:,:,3) = imRGB(:,:,3) + im(:,:,5); % Magenta    
    figure; imshow(imRGB,[])
    
else if im_s(3)==4;       
    imRGB = zeros([im_s(1),im_s(2),3]);
    imRGB(:,:,3) = imRGB(:,:,3) + im(:,:,1); % Blue
    imRGB(:,:,2) = imRGB(:,:,2) + im(:,:,2); % Green
    imRGB(:,:,1) = imRGB(:,:,1) + im(:,:,3); % Red
    imRGB(:,:,1) = imRGB(:,:,1) + im(:,:,4); imRGB(:,:,2) = imRGB(:,:,2) + im(:,:,4); % Yellow    
    figure; imshow(imRGB,[])

    for ii=1:size(im,3); im(:,:,ii)=im(:,:,ii)./max(max(im(:,:,ii))); end;
    im_s = size(im);
    imRGB = zeros([im_s(1),im_s(2),3]);
    imRGB(:,:,3) = imRGB(:,:,3) + im(:,:,1); % Blue
    imRGB(:,:,2) = imRGB(:,:,2) + im(:,:,2); % Green
    imRGB(:,:,1) = imRGB(:,:,1) + im(:,:,3); % Red
    imRGB(:,:,1) = imRGB(:,:,1) + im(:,:,4); imRGB(:,:,2) = imRGB(:,:,2) + im(:,:,4); % Yellow    
    figure; imshow(imRGB,[])
    
else im_s(3)==3;
    imRGB = zeros([im_s(1),im_s(2),3]);
    imRGB(:,:,3) = imRGB(:,:,3) + im(:,:,1); % Blue
    imRGB(:,:,2) = imRGB(:,:,2) + im(:,:,2); % Green
    imRGB(:,:,1) = imRGB(:,:,1) + im(:,:,3); % Red
    figure; imshow(imRGB,[])
    
    for ii=1:size(im,3); im(:,:,ii)=im(:,:,ii)./max(max(im(:,:,ii))); end;
    im_s = size(im);
    imRGB = zeros([im_s(1),im_s(2),3]);
    imRGB(:,:,3) = imRGB(:,:,3) + im(:,:,1); % Blue
    imRGB(:,:,2) = imRGB(:,:,2) + im(:,:,2); % Green
    imRGB(:,:,1) = imRGB(:,:,1) + im(:,:,3); % Red
    figure; imshow(imRGB,[])
end
    end
end
end

