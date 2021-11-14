function [Image_Pre, Image_Seg, Image_lbp, DataLBP, Class] = TaskHMI(I)

I_G = rgb2gray(I);
I_G = im2double(I_G);

Fft_I = fft2(log(I_G + 0.01));
Fil_I = Butterworth_HighPassFilter(I_G, 500, 3);
temp = Fft_I .* Fil_I;
h = real(ifft2(temp));
h = exp(h);
Img = ifftshow(h);
Ventana = [80 16 630 550];
Img1 = imcrop(I_G,Ventana);
Image_Pre = Img1;

[Img] = segment(Img);
I_R = I(:,:,1);
I_G = I(:,:,2);
I_B = I(:,:,3);
Img_seg = uint8(Img);
I_R = I_R .* Img_seg;
I_G = I_G .* Img_seg;
I_B = I_B .* Img_seg;
Img_seg = cat(3, I_R, I_G, I_B);
Image_Seg = Img_seg;


%++++++Algortimo LBP
grayImage=Image_Pre;
% Get the dimensions of the image.  numberOfColorBands should be = 1.
[rows columns numberOfColorBands] = size(grayImage);
% Preallocate/instantiate array for the local binary pattern.
localBinaryPatternImage1 = zeros(size(grayImage), 'uint8');
localBinaryPatternImage2 = zeros(size(grayImage), 'uint8');
localBinaryPatternImage3 = zeros(size(grayImage), 'uint8');
localBinaryPatternImage4 = zeros(size(grayImage), 'uint8');
localBinaryPatternImage5 = zeros(size(grayImage), 'uint8');
localBinaryPatternImage6 = zeros(size(grayImage), 'uint8');
localBinaryPatternImage7 = zeros(size(grayImage), 'uint8');
localBinaryPatternImage8 = zeros(size(grayImage), 'uint8');
tic;
for row = 2 : rows - 1
	for col = 2 : columns - 1
		centerPixel = grayImage(row, col);
		pixel7=grayImage(row-1, col-1) > centerPixel;
		pixel6=grayImage(row-1, col) > centerPixel;
		pixel5=grayImage(row-1, col+1) > centerPixel;
		pixel4=grayImage(row, col+1) > centerPixel;
		pixel3=grayImage(row+1, col+1) > centerPixel;
		pixel2=grayImage(row+1, col) > centerPixel;
		pixel1=grayImage(row+1, col-1) > centerPixel;
		pixel0=grayImage(row, col-1) > centerPixel;
		
		% Create LBP image with the starting, LSB pixel in the upper left.
		eightBitNumber = uint8(...
			pixel7 * 2^7 + pixel6 * 2^6 + ...
			pixel5 * 2^5 + pixel4 * 2^4 + ...
			pixel3 * 2^3 + pixel2 * 2^2 + ...
			pixel1 * 2 + pixel0);
		% Or you can use the built-in function bwpack(), which is somewhat simpler but a lot slower.
		% 		eightBitNumber = uint8(bwpack([pixel0; pixel1; pixel2; pixel3; pixel4; pixel5; pixel6; pixel7]));
		localBinaryPatternImage1(row, col) = eightBitNumber;
		
		% Create LBP image with the starting, LSB pixel in the upper middle.
		eightBitNumber = uint8(...
			pixel6 * 2^7 + pixel5 * 2^6 + ...
			pixel5 * 2^4 + pixel3 * 2^4 + ...
			pixel3 * 2^2 + pixel1 * 2^2 + ...
			pixel0 * 2 + pixel7);
		% Or you can use the built-in function bwpack(), which is somewhat simpler but a lot slower.
		% 		eightBitNumber = uint8(bwpack([pixel0; pixel1; pixel2; pixel3; pixel4; pixel5; pixel6; pixel7]));
		localBinaryPatternImage2(row, col) = eightBitNumber;
		
		% Create LBP image with the starting, LSB pixel in the upper right.
		eightBitNumber = uint8(...
			pixel5 * 2^7 + pixel4 * 2^6 + ...
			pixel3 * 2^5 + pixel2 * 2^4 + ...
			pixel1 * 2^3 + pixel0 * 2^2 + ...
			pixel7 * 2 + pixel6);
		% Or you can use the built-in function bwpack(), which is somewhat simpler but a lot slower.
		% 		eightBitNumber = uint8(bwpack([pixel0; pixel1; pixel2; pixel3; pixel4; pixel5; pixel6; pixel7]));
		localBinaryPatternImage3(row, col) = eightBitNumber;
		
		% Create LBP image with the starting, LSB pixel in the center right.
		eightBitNumber = uint8(...
			pixel4 * 2^7 + pixel3 * 2^6 + ...
			pixel2 * 2^5 + pixel1 * 2^4 + ...
			pixel0 * 2^3 + pixel7 * 2^2 + ...
			pixel6 * 2 + pixel5);
		% Or you can use the built-in function bwpack(), which is somewhat simpler but a lot slower.
		% 		eightBitNumber = uint8(bwpack([pixel0; pixel1; pixel2; pixel3; pixel4; pixel5; pixel6; pixel7]));
		localBinaryPatternImage4(row, col) = eightBitNumber;
		
		% Create LBP image with the starting, LSB pixel in the lower right.
		eightBitNumber = uint8(...
			pixel3 * 2^7 + pixel2 * 2^6 + ...
			pixel1 * 2^5 + pixel0 * 2^4 + ...
			pixel7 * 2^3 + pixel6 * 2^2 + ...
			pixel5 * 2 + pixel0);
		% Or you can use the built-in function bwpack(), which is somewhat simpler but a lot slower.
		% 		eightBitNumber = uint8(bwpack([pixel0; pixel1; pixel2; pixel3; pixel4; pixel5; pixel6; pixel7]));
		localBinaryPatternImage5(row, col) = eightBitNumber;
		
		% Create LBP image with the starting, LSB pixel in the lower center.
		eightBitNumber = uint8(...
			pixel2 * 2^7 + pixel1 * 2^6 + ...
			pixel0 * 2^5 + pixel7 * 2^4 + ...
			pixel6 * 2^3 + pixel5 * 2^2 + ...
			pixel4 * 2 + pixel3);
		% Or you can use the built-in function bwpack(), which is somewhat simpler but a lot slower.
		% 		eightBitNumber = uint8(bwpack([pixel0; pixel1; pixel2; pixel3; pixel4; pixel5; pixel6; pixel7]));
		localBinaryPatternImage6(row, col) = eightBitNumber;
		
		% Create LBP image with the starting, LSB pixel in the lower left.
		eightBitNumber = uint8(...
			pixel1 * 2^7 + pixel0 * 2^6 + ...
			pixel7 * 2^5 + pixel6 * 2^4 + ...
			pixel5 * 2^3 + pixel4 * 2^2 + ...
			pixel3 * 2 + pixel2);
		% Or you can use the built-in function bwpack(), which is somewhat simpler but a lot slower.
		% 		eightBitNumber = uint8(bwpack([pixel0; pixel1; pixel2; pixel3; pixel4; pixel5; pixel6; pixel7]));
		localBinaryPatternImage7(row, col) = eightBitNumber;
		
		% Create LBP image with the starting, LSB pixel in the center left.
		eightBitNumber = uint8(...
			pixel0 * 2^7 + pixel7 * 2^6 + ...
			pixel6 * 2^5 + pixel5 * 2^4 + ...
			pixel4 * 2^3 + pixel3 * 2^2 + ...
			pixel2 * 2 + pixel1);
		% Or you can use the built-in function bwpack(), which is somewhat simpler but a lot slower.
		% 		eightBitNumber = uint8(bwpack([pixel0; pixel1; pixel2; pixel3; pixel4; pixel5; pixel6; pixel7]));
		localBinaryPatternImage8(row, col) = eightBitNumber;
		
	end
end
toc;
% Outer layer of pixels will be zero because they didn't have 8 neighbors.
% So, to avoid a huge spike in the histogram at zero, replace the outer layer of pixels with the next closest layer.
localBinaryPatternImage1(1, :) = localBinaryPatternImage1(2, :);
localBinaryPatternImage1(end, :) = localBinaryPatternImage1(end-1, :);
localBinaryPatternImage1(:, 1) = localBinaryPatternImage1(:, 2);
localBinaryPatternImage1(:, end) = localBinaryPatternImage1(:, end-1);

%*****Devuelvo Valores LBP*******
Image_lbp = localBinaryPatternImage1;
DataLBP= extractLBPFeatures(Image_Pre,'NumNeighbors',8,'Radius',2);

%Clasificacion
B = readmatrix('D:\VIU\04_Asignaturas\21_Trabajo de Final\Programa\Modelo_MRL.xlsx');;
mtest1 = 1;
xtest1 = [ones(mtest1,1) DataLBP];
ztest1 = xtest1*B(:,1);
ztest2 = xtest1*B(:,2);
htt1 = exp(ztest1)./(1.0+exp(ztest1)+exp(ztest2));
htt2 = exp(ztest2)./(1.0+exp(ztest1)+exp(ztest2));
htt3 = 1.0./(1.0+exp(ztest1)+exp(ztest2));
Aux = 0;
if ((htt1 > htt2) &&(htt1 > htt3));
Aux = 1; 
end
if ((htt2 > htt1) &&(htt2 > htt3));
Aux = 2; 
end
if ((htt3 > htt1) &&(htt3 > htt2));
Aux= 3; 
end
Class = Aux;

end

function [t] = otsu_impl(counts)
num_bins = numel(counts);
counts = double( counts(:) );
p = counts / sum(counts);
omega = cumsum(p);
mu = cumsum(p .* (1:num_bins)');
mu_t = mu(end);
sigma_b_squared = (mu_t * omega - mu).^2 ./ (omega .* (1 - omega));
maxval = max(sigma_b_squared);
isfinite_maxval = isfinite(maxval);
if isfinite_maxval
    idx = mean(find(sigma_b_squared == maxval));
    t = (idx - 1) / (num_bins - 1);
else
    t = 0.0;
end
end

function [level] = thresh_impl(I)
if ~isempty(I)
    I = im2uint8(I(:));
    num_bins = 256;
    counts = imhist(I,num_bins);
    level = otsu_impl(counts);
else
    level = 0.0;
end
end

function [out] = Butterworth_HighPassFilter(I, d, n)
h = size(I, 1);
w = size(I, 2);
[x, y] = meshgrid(-floor(w / 2):floor(w - 1) / 2, -floor(h / 2):floor(h - 1) / 2);
out = 1 ./ (1 + (d ./ (x .^ 2 + y .^ 2) .^ 0.5) .^ (2 * n));
end

function [I] = ifftshow(f)
fabs = abs(f);
fmax = max(fabs(:));
I = (fabs / fmax);
end

function [outImage] = segment(I)

[m, n] = size(I);
level = thresh_impl(I);     
I_Bin = imbinarize(I,level);   
K = medfilt2(I_Bin);
se = strel('disk', 4);
hairs = imbothat(K,se);      
woHair = K;
for i = 1:m
    for j =  1:n
        if hairs(i, j) == 1
            woHair(i,j) = 1;
        end
    end
end
woHair = imerode(woHair, se);
woHair_Edge = edge(woHair,'Canny');
se90 = strel('line', 2, 90);
se0 = strel('line', 2, 0);
woHair_Edge_Dil = imdilate(woHair_Edge, [se90 se0]);
outImage = imfill(woHair_Edge_Dil, 'holes');
outImage = bwareafilt(outImage, 1);
end