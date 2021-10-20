
clear; close all force; clc;

start_path = ...
    'D:\BOL\Project2\2021.06.15\LA2021.06.15\.25percent';
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

slope_all = size(zeros(length(files),BScanWidth));
exp_decay_all= size(zeros(length(files),BScanWidth));
sd_aline_all = size(zeros(length(files),BScanWidth));

for i=1:length(files)
%     test_image=images(:,:,1);
    test_image=images(:,:,i);
    for pos=1:BScanWidth
        
        figure(1); imagesc(test_image,[30,120]); colormap(gray);
%         aline=test_image(:,41);
         aline=test_image(:,pos);
        figure(2); plot(aline);
        
        window=100:500;
        
        aline_cropped=aline(window);
        
        x=1:length(aline_cropped);
        
        [xData, yData] = prepareCurveData( x, aline_cropped );
        
        % Set up fittype and options.
        ft = fittype( 'poly1' );
        opts = fitoptions( 'Method', 'LinearLeastSquares' );
        opts.Robust = 'Bisquare';
        
        % Fit model to data.
        [fitresult, gof] = fit( xData, yData, ft, opts );
        
        linear_fit_slope=fitresult.p1;
        slope_all(i,pos)=linear_fit_slope;
        
        fittedLine=feval(fitresult,x);
        %figure(3); plot(x,fittedLine);
        
        aline_subtracted=aline_cropped-fittedLine;
        figure(3); plot(aline_subtracted);
        sd_aline=std(aline_subtracted);
        sd_aline_all(i,pos)=sd_aline;
        
        mean_aline=mean(aline_subtracted);
        
         figure(4); plot(aline_cropped); hold on; plot(fittedLine); hold off;
        
        
        X=aline_subtracted;
        x_corr=x*5.68;                  %per pixel resulation is 5.68
        Fs = 1/5.68;                    % Sampling frequency(distance per sampling)
        T = 5.68;                       % Sampling period 
        L = length(x_corr);             % Length of signal
        t = x_corr;                     % Time vector
        
        
        figure(4);
        plot(x,X)
        title('Signal Corrupted with Zero-Mean Random Noise')
        xlabel('t (microns)')           % i have to covert to db
        ylabel('X(t)')
        
        
        Y = fft(X);
        P2 = abs(Y/L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
        
        f = Fs*(0:(L/2))/L;
        figure(5); plot(f,P1)
        title('Single-Sided Amplitude Spectrum of X(t)')
        xlabel('f (Hz)')
        ylabel('|P1(f)|')
        
        % Fit: 'untitled fit 1'.
        [xData, yData] = prepareCurveData( f, P1 );
        
        % Set up fittype and options.
        ft = fittype( 'exp1' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Display = 'Off';
        opts.StartPoint = [0.868454356691457 -21.1707637609794];
        
        % Fit model to data.
        [fitresult2, ~] = fit( xData, yData, ft, opts );
        y_exp=feval(fitresult2,f);
        figure(6); hold on; plot(f,y_exp); hold off;
        exp_decay=fitresult2.a;
         
         exp_decay_all(i,pos)=exp_decay;
    end
 end

 
        


close all force;
std_all= mean(sd_aline_all, 2);
slope= mean(slope_all, 2);
exponential_decay= mean(exp_decay_all, 2);

boxplot(all_data);
figure(7); boxplot(slope);
figure(8); boxplot(std_all);
figure(9); boxplot(exponential_decay);