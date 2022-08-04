%% Reverse correlation

function [rfModel,rfModelVec,rfMap2d,rfMapVec,stimMovie_all,resp_all,vaf,output_actual,stimMovie]=reverseCorr(n_movies,stimulus_type)

imgSiz  = 32;     
nPixels = imgSiz^2; 
nFrames = 375;

% specify model receptive field (Gabor function followed by half-power law)
model.lambda = 8;
model.phase  = 0;
model.ori    = 0;  
model.pwrExp = 2; 

%  create model filter
rfModel = makeModelRF(model,imgSiz);    % creates model filter (Gabor function)
rfModel=rfModel/max(max(abs(rfModel)));
rfModelVec  = reshape(rfModel,1,nPixels);      % make a 1d version, for later use            

stimMovie_all = [];
resp_all  = [];

for iMovie=1:n_movies
    getStimulusMovies;            % -> stimMovie = nPixels x nFrames, range -1 to +1
    output = rfModelVec*stimMovie;    % linear filter response to the stimulus    
    output = hwr(output);  % half-wave rectify (set negative values to zero)
    %output = output.^model.pwrExp;  % power-low for positive values
                   
    stimMovie_all = [stimMovie_all stimMovie]; % nPixels x nFrames*nMovies.train
    resp_all    = [resp_all    output];        %     1   x nFrames*nMovies.train  
                   
end  % end of iMovie-loop

crossCorrAll = zeros(nPixels,n_movies*nFrames);
for iPix=1:nPixels
    crossCorrAll(iPix,:) = xcorr(resp_all,stimMovie_all(iPix,:),0,'unbiased');  % cross-correlation
end
rfMap = crossCorrAll(:,end);  % only take positive lagged values
rfMap2d = reshape(rfMap,imgSiz,imgSiz);  % reshape into 2d
rfMap2d=rfMap2d/max(max(abs(rfMap2d)));
rfMapVec  = reshape(rfMap2d,nPixels,1);    % make a 1-d version

% test movie
iMovie=8;
getStimulusMovies;  
output = rfModelVec*stimMovie;
output_actual = hwr(output); 


output_pred=rfMapVec'*stimMovie;

vaf.R_matrix = corrcoef(output_pred,output_actual);  % -> 2x2 matrix, ones on diagonal
vaf.offDiag = vaf.R_matrix(1,2);
vaf= vaf.offDiag^2.;
vaf=vaf*100;
end