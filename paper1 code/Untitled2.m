clear; close all force; clc;

%for i= 0: 11
    %start_path = ...
        %['D:\BOL\Project2\LA2021.06.24_600x1x15\mirror',num2str(i)];
start_path = ...
    'D:\BOL\Project2\LA2021.06.24_600x1x15\mirror0';
pathname = uigetdir(start_path,'Select directory with the phase files');
data_type = 'uint16';

files = dir(fullfile(pathname,'*.dat'));
spec_filename = (dir(fullfile(start_path,'*.spectrum')));
spectrum = (double(load(fullfile(start_path,spec_filename(1).name))));
spec_len = length(spectrum);
[xData, yData] = prepareCurveData( [], spectrum );
% Set up fittype and options.
ft = fittype( 'poly2' );
opts = fitoptions( 'Method', 'LinearLeastSquares' );
opts.Robust = 'Bisquare';
% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
fit_spec = feval(fitresult,1:spec_len);
k_space = (2*pi./fit_spec);
new_ks = ((k_space(1) -(0:(spec_len-1))*(k_space(1)-k_space(spec_len))/(spec_len-1)));

BScanWidth = 600;
top_depth = 10;
how_many_depths = 600;
hann_repmat = (double(repmat(hann(spec_len),[1,BScanWidth])));
images=zeros(how_many_depths,BScanWidth,length(files));
% mkdir(pathname,'images_2');
for i = 1:length(files)
    filename = files(i).name;
    fid = fopen(fullfile(pathname,filename));
    
    switch data_type
        
        case 'uint8' % 8 bit data from 8 bit camera
            % load the raw fringe data from *.dat file
            raw_fringes = double(fread(fid,[spec_len,BScanWidth],'uint8',0,'b'));
            fclose(fid);
            % ensure that ks are increasing
            if new_ks(end) < new_ks(1)
                % interpolate/resample raw fringe data
                temp = single(interp1(flip(k_space),double(flip(raw_fringes)),...
                    flip(new_ks),'linear')); %#ok<*PFTIN>
            else
                % interpolate/resample raw fringe data
                temp = single(interp1(k_space,double(raw_fringes),...
                    new_ks,'linear'));
            end
            
        case 'int16' % 16 bit data from 16 bit cameras
            % load the raw fringe data from *.dat file
            raw_fringes = double(fread(fid,[spec_len,BScanWidth],'int16',0,'b'));
            fclose(fid);
            % ensure that ks are increasing
            if new_ks(end) < new_ks(1)
                % interpolate/resample raw fringe data
                temp = single(interp1(flip(k_space),double(flip(raw_fringes)),...
                    flip(new_ks),'linear'));
            else
                % interpolate/resample raw fringe data
                temp = single(interp1(k_space,double(raw_fringes),...
                    new_ks,'linear'));
            end
            
        case 'uint16' % 16 bit data from PhS-SSOCT
            % load the raw fringe data from *.dat file
            raw_fringes = double(fread(fid,[spec_len,BScanWidth],'uint16',0,'b'));
            fclose(fid);
            % ensure that ks are increasing
            if new_ks(end) < new_ks(1)
                % interpolate/resample raw fringe data
                temp = single(interp1(flip(k_space),double(flip(raw_fringes)),...
                    flip(new_ks),'linear'));
            else
                % interpolate/resample raw fringe data
                temp = single(interp1(k_space,double(raw_fringes),...
                    new_ks,'linear'));
            end
            
    end  
    bg_fringe = median(temp,2);
    
    fft_1 = fft(bsxfun(@times,hann_repmat,temp-bg_fringe),[],1);
    fft_1 = fft_1(1:1600,:);
    image = 20*log10(abs (fft_1(top_depth:top_depth+how_many_depths-1,:)));
    images(:,:,i)=image; 
%     write_image = mat2gray(image,[30,80]);
%     
%     figure(1); imagesc(image); caxis([30,80]); colormap(gray);
%     set(gcf,'units','normalized','outerposition',[0,0,1,1]);
%     title([num2str(i),' of ',num2str(length(files))]);pause(0.01);
%     imwrite(write_image,fullfile(pathname,'images_2',['frame',num2str(i,'%05d'),'.jpg']),'jpg');
    
end

figure(1); imagesc(image,[40 120]); colormap(gray);
% imagesc(image(:, 1));
% aline= (image(:, 1));
% figure(2); plot(aline(:, 1));
surface= zeros(1, BScanWidth);

% for i = 1: length(aline)
%     
%     surface(i) = find(aline>80, 1);
%     
% end


for i = 1: 600
   
    aline= image(:, i);
    surface(i) = find(aline>74.645, 1);
    
end
plot(aline);
s_surface=smooth(surface,0.1);

hold on; plot(surface,'.r','MarkerSize',5); hold off;
mirror_peak= surface;

[~,b,~]=fileparts(pathname);
save(fullfile(pathname,['surface_',b,'.mat']),'surface');
num2clip
