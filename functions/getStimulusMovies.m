%% getStimulusMovies.m
% script to create a stimulus movie, either white noise or series of natural images


if strcmp(stimulus_type,'white') 
    stimMovie  = 2*(rand(nPixels,nFrames) - 0.5);  % white noise, range -1 to +1
    
elseif strcmp(stimulus_type,'McGill_clips')   % these are 32^2 x 375
    % read file
    files.thisFileName = ['McGill_clips_0' int2str(iMovie) '.mat'];
    fprintf(1,'loading stimulus file: %s\n',files.thisFileName);
    load(files.thisFileName); 
    xMovie = thisNewMovie;
    % create stimMovie, range -1 to +1, imgSiz^2 x 375
    stimMovie = double(reshape(xMovie,nPixels,nFrames));
    stimMovie = stimMovie/128;  % range -1 to +1 
end
    





