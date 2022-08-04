clear all;
close all;

addpath('functions');
addpath('data');

%% Reverse correlation with white noise

%n_movies  = 7;
imgSiz  = 32; 
nFrames = 375;
nPixels = imgSiz^2; 

figure(); 

N_movies=[1,10,100,1000];

for N_movies_=1:4
    n_movies  = N_movies(N_movies_);

    [rfModel,rfModelVec,rfMap2d_white,rfMapVec_white,stimMovie_white,resp_white,VAF,output_actual,stimMovie]=reverseCorr(n_movies,'white');

    subplot(1,5,1);
    imagescZadj(rfModel);  
    title('Actual RF')

    subplot(1,5,N_movies_+1);
    imagescZadj(rfMap2d_white);
    title(sprintf('Ns = %.0f',n_movies*375))

    mse(N_movies_)=sum((rfModelVec'-rfMapVec_white).^2)/1024;
    vaf(N_movies_)=VAF;
end

N_frames=N_movies*375; % each movie has 375 frames

% plot MSE, VAF variation with number of stimuli
figure()
plot(N_frames,mse,'Linewidth',1.5)
set(gca, 'XScale', 'log'); grid on;
xlim([375,375000])
xlabel('Number of stimuli (Ns)'); ylabel('Mean squarred error (MSE)'); title('MSE variation with number of stimuli'); 
figure();
plot(N_frames,vaf,'Linewidth',1.5)
set(gca, 'XScale', 'log'); grid on;
xlim([375,375000])
xlabel('Number of stimuli (Ns)'); ylabel('VAF (%)');  ; title('VAF variation with number of stimuli'); 


% check the pixel wise correlations in white noise
cross_corr_white=corrcoef(stimMovie_white'); 
figure();
imagesc(cross_corr_white);
xlabel('Pixel #'); ylabel('Pixel #'); title('Pixel-wise correlations of white noise');

% visualize few example white noise images used
figure();
for example=1:4
    stimImgVec = stimMovie_white(:,example);
    stimImg = reshape(stimImgVec,imgSiz,imgSiz);
    subplot(2,2,example);
    imagescZadj(stimImg); colormap('gray'); xlabel('Pixel #'); ylabel('Pixel #');
end
sgtitle('Example white noise stimuli')

%% Reverse correlation with natural image movies

n_movies=7;
[rfModel,rfModelVec,rfMap2d_natural,rfMapVec_natural,stimMovie_natural,resp_natural,vaf,output_actual,stimMovie]=reverseCorr(n_movies,'McGill_clips');

figure(); 
subplot(1,2,1);
imagescZadj(rfModel);  title('RF model filter'); 
subplot(1,2,2);
imagescZadj(rfMap2d_natural); title('ReverseCorr-Natural images');  


% check the pixel wise correlations in natural images
cross_corr_natural=corrcoef(stimMovie_natural'); 
figure();
imagesc(cross_corr_natural);
xlabel('Pixel #'); ylabel('Pixel #'); title('Pixel-wise correlations of natural images');

% visualize few example natural images used
figure();
for example=1:4
    stimImgVec = stimMovie_natural(:,example+50); %take images 50,100,150 and 200
    stimImg = reshape(stimImgVec,imgSiz,imgSiz);
    subplot(2,2,example);
    imagescZadj(stimImg); colormap('gray'); xlabel('Pixel #'); ylabel('Pixel #');
end
sgtitle('Example natural image stimuli')

%% Testing the relationship between the corr matrix and receptive field estimate
% Note: This is the euqation5 of B Willmore and D Smyth, 2003 paper
% for non-orthogonal stimuli,the reverse-correlation receptive-field estimate
% will be related to the true receptive field as below.

est_rf=cross_corr_natural*rfModelVec';
rfMap2d_est = reshape(est_rf,imgSiz,imgSiz);

figure(); 
subplot(1,2,1);
imagescZadj(rfMap2d_natural); title('RF-ReverseCorr with NaturalImages'); colorbar; 
subplot(1,2,2);
imagescZadj(rfMap2d_est); title('RF-construct with equation'); colorbar; 

%% Correction of the RF-reverseCorr of natural images using the power spectrum
% this is one way to correct the RF estimates from natural images using the
% power spectrum of natural images. Ref: Section 2.3 of B Willmore and D Smyth, 2003

power=0;
for i=1:n_movies*nFrames
    fourier_trans=fft(stimMovie_natural(:,i)); % fft
    power=power+abs(fourier_trans.^2); % power spectrum 
end

% fourier transform of 1d RF-reverseCorr of natural images
fourier_rf=fft(rfMapVec_natural);

% correct the fourier transform of 1d RF-reverseCorr of natural images
corrected_fourier=(n_movies*nFrames)*fourier_rf./power; 

% inverse fourier tranform of the corrected
corrected_rf = ifft(corrected_fourier);
corrected_rf2d = reshape(corrected_rf,imgSiz,imgSiz);

% visualize
figure()
subplot(1,3,1)
imagescZadj(rfModel); title('RF model filter')
subplot(1,3,2)
imagescZadj(rfMap2d_natural); title('ReverseCorr-NaturalImages')
subplot(1,3,3)
imagescZadj(corrected_rf2d); title('ReverseCorr-NaturalImages-corrected')

%% Alternative spectrum Correction 
% this is another way to correct the RF estimates from natural images using the
% amplitude spectrum of natural images. Ref: Section 2.3.1 of B Willmore and D Smyth, 2003

stim_ampSpectrum=0;
for i=1:n_movies*nFrames
    fourier_trans=fft(stimMovie_natural(:,i)); % fft
    stim_ampSpectrum=stim_ampSpectrum+abs(fourier_trans); % amplitude spectrum 
end

% correct the fourier transform of 1d RF-reverseCorr of natural images
corrected_fourier2=(n_movies*nFrames)*fourier_rf./stim_ampSpectrum; 

% inverse fourier tranform of the corrected
corrected_rf2 = ifft(corrected_fourier2);
corrected_rf2d2 = reshape(corrected_rf2,imgSiz,imgSiz);

% visualize
figure()
subplot(1,3,1)
imagescZadj(rfModel); title('RF model filter')
subplot(1,3,2)
imagescZadj(rfMap2d_natural); title('ReverseCorr-NaturalImages')
subplot(1,3,3)
imagescZadj(corrected_rf2d2); title('ReverseCorr-NaturalImages-corrected')


%% Testing RF correction with increasing number of movies

figure()
for n_movies = 3:1:7
    iteration=n_movies-2;
    [rfModel,rfModelVec,rfMap2d_natural,rfMapVec_natural,stimMovie_natural,resp_natural,vaf,output_actual,stimMovie_test]=reverseCorr(n_movies,'McGill_clips');
    VAF_before(iteration)=vaf;

    mse_before(iteration)=sum((rfModelVec'-rfMapVec_natural).^2)/1024;

    subplot(3,5,iteration)
    imagescZadj(rfMap2d_natural); 
    title(sprintf('Ns = %.0f',n_movies*375))

    %% RF correction
    power=0;
    for i=1:n_movies*nFrames
        fourier_trans=fft(stimMovie_natural(:,i)); % fft
        power=power+abs(fourier_trans.^2); % power spectrum 
    end

    % fourier transform of 1d RF-reverseCorr of natural images
    fourier_rf=fft(rfMapVec_natural);

    % correct the fourier transform of 1d RF-reverseCorr of natural images
    corrected_fourier=(n_movies*nFrames)*fourier_rf./power; 

    % inverse fourier tranform of the corrected
    corrected_rf = ifft(corrected_fourier);
    corrected_rf=corrected_rf/max(max(abs(corrected_rf)));
    corrected_rf2d = reshape(corrected_rf,imgSiz,imgSiz);
    
    subplot(3,5,5+iteration)
    imagescZadj(corrected_rf2d);
    
    % MSE 
    mse_1(iteration)=sum((rfModelVec'-corrected_rf).^2)/1024;
    
    % VAF
    output_pred=corrected_rf'*stimMovie_test;
    vaf = corrcoef(output_pred,output_actual);  % -> 2x2 matrix, ones on diagonal
    vaf = vaf(1,2);
    vaf= vaf^2.;
    vaf=vaf*100;
    VAF_test1(iteration)=vaf;
    
    %% RF correction method2
    
    stim_ampSpectrum=0;
    for i=1:n_movies*nFrames
        fourier_trans=fft(stimMovie_natural(:,i)); % fft
        stim_ampSpectrum=stim_ampSpectrum+abs(fourier_trans); % power spectrum 
    end

    % correct the fourier transform of 1d RF-reverseCorr of natural images
    corrected_fourier2=(n_movies*nFrames)*fourier_rf./stim_ampSpectrum; 

    % inverse fourier tranform of the corrected
    corrected_rf2 = ifft(corrected_fourier2);
    corrected_rf2=corrected_rf2/max(max(abs(corrected_rf2)));
    corrected_rf2d2 = reshape(corrected_rf2,imgSiz,imgSiz);
    
    subplot(3,5,10+iteration)
    imagescZadj(corrected_rf2d2); 
    
    % MSE 
    mse_2(iteration)=sum((rfModelVec'-corrected_rf2).^2)/1024;
    % VAF
    output_pred=corrected_rf2'*stimMovie_test;
    vaf2 = corrcoef(output_pred,output_actual);  % -> 2x2 matrix, ones on diagonal
    vaf2 = vaf2(1,2);
    vaf2= vaf2^2.;
    vaf2=vaf2*100;
    VAF_test2(iteration)=vaf2;
    
end

figure();
plot((3:1:7)*375,mse_1,'Linewidth',1.5)
hold on; grid on;
plot((3:1:7)*375,mse_2,'Linewidth',1.5)
legend('Correction method 1','Correction method 2')
xlabel('Number of stimuli (Ns)'); ylabel('MSE');  ; title('MSE variation with number of stimuli'); 

figure();
plot((3:1:7)*375,VAF_test1,'Linewidth',1.5)
hold on; grid on;
plot((3:1:7)*375,VAF_test2,'Linewidth',1.5)
legend('Correction method 1','Correction method 2')
xlabel('Number of stimuli (Ns)'); ylabel('VAF %');  ; title('VAF variation with number of stimuli'); 

