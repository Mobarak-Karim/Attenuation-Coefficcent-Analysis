

% Step 1: Remove fixed-pattern noise by subtracting mean of many A-lines
% from each A-line.
figure; imagesc(image);colormap(gray)
for k = 1:size(image,2)
    newImage(:,k) = image(:,k)-mean(image,2);
end

figure; imagesc(newImage); colormap(gray)

% Step 2: Remove noise floor signal N(z)

load('D:\BOL\Project2\2nd code\NoiseFloor.mat','NoiseFloor_N_z')
N = NoiseFloor_N_z;
% N = mean(background_signal,2);

% for k = 1:size(image,2)
%     N(:,k) = image(:,k)-noiseFloorSignal;
% end

% Step 3: Model signal decay S(z)

zRange = 1:size(image,1);
sigma =4.7* 10^-3;
for k = 1:length(zRange)
S(:,k) = exp((-zRange(k)^2)/(sigma^2))
end


% Step 4: Calculate corrected data I(z)

for k = 1:size(image,2)
    I(:,k) = ((image(:,k) - N) / S(:,k));
    
end

figure; imagesc(I); colormap(gray);
