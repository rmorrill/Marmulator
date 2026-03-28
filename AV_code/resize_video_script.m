% script to resize videos using ffmpeg 
movie_folder = '/home/ryan/Documents/Stimuli/AV_samples/1121Stimuli/resized_fixed_group'; 
new_movie_folder = '/home/ryan/Documents/Stimuli/AV_samples/1121Stimuli/resized_fixed_group/resized_fixed_group_small'; 
fnames = dir([movie_folder '/*.mov']);
disp(fnames)
%%
for i = 1:numel(fnames)
    curr_fname = fnames(i).name; 
    new_fname = strrep(curr_fname, '.mov', '_small.mov')
    cmdResize = sprintf('ffmpeg -i "%s" -vf "scale=500:500" "%s"', ...
        fullfile(movie_folder,curr_fname), fullfile(new_movie_folder, new_fname)); 
    system(cmdResize); 
end

