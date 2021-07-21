function [fname] = h5Animation2D(input_foldername, extra_params)
%H5ANIMATION2D makes an animation of the results from an h5 file
%   INPUTS:
%       input_foldername  (string): the name of the folder where the results are
%       extra_params (struct): all extra parameters for making the animation
%       extra_params.real_time (float): how long the animation should go
%           for in real time in seconds [DEFAULT: 10.0 seconds]
%       extra_params.sim_interval (float, size = (2,)): time intervals
%           between which the result should be animated (inclusive). [DEFAULT:
%           full range of t values]
%       extra_params.fps (float): the frames per second of the output
%           animation [DEFAULT: 20 fps]
%       extra_params.dpi (int): the dots per inch of the output animation
%           [DEFAULT: 200]
%       extra_params.interpolate_resolution (int, size = (2,)): whether to interpolate
%           the result to a finer mesh before plotting. Set to 0, or set one of
%           the values to 0 to prevent interpolation. [DEFAULT: 0]
%       extra_params.range_max (float, size =({# of chemical species},)):
%           sets the scale of the plot. Each value is the maximum of the range
%           for the respective chemical species.
%       extra_params.range_min (float, size =({# of chemical species},)):
%           sets the scale of the plot. Each value is the minimum of the range
%           for the respective chemical species.

   
    fname = strcat(input_foldername, '/results.h5');
    if nargin <2
        extra_params.fps=20;
    end
    
    try
        t = h5read(fname, "/t");
    catch e
        error("Error in reading full range of time values")
    end
    % find the tange of t values
    if ~isfield(extra_params, 'sim_interval')
        time_index_start = 1; % time index start
        
        info = h5info(fname, "/t");
        time_index_end = info.Dataspace.Size;
    else
        % check that time values are in range
        if extra_params.sim_interval(1) < 0
            error("extra_params.sim_interval(1) out of range")
        end
        
        info = h5info(fname, "/t");
        
        
        if extra_params.sim_interval(2) > h5read(fname, "/t", double(info.Dataspace.Size), double(1))
            error("extra_params.sim_interval(2) out of range")
        end
        [~, time_index_start] = max(t >= extra_params.sim_interval(1));
        [~, time_index_end] = max(t >= extra_params.sim_interval(2));
    end
    
    
    % find the number of frames
    if ~isfield(extra_params, "fps")
        extra_params.fps = 20;
    end
    if ~isfield(extra_params, "real_time")
        extra_params.real_time = 10;
    end
    
    number_of_frames = extra_params.fps * extra_params.real_time;
    t_indices = linspace(time_index_start, time_index_end, number_of_frames);
    
    % set up interpolation
    if ~isfield(extra_params, 'interpolate_resolution')
        extra_params.interpolate_resolution = 0;
    end
    
    if any(extra_params.interpolate_resolution == 0)
        interp = false;
    else
        interp = true;
    end
    
    % load x and y values, and calculate meshgrid
    try
        x = h5read(fname, "/x");
        y = h5read(fname, "/y");
        [xx, yy] = meshgrid(x, y);
    catch e
       error("error in reading x and y values") 
    end
    
    % set up interpolation grid
    if interp
        x_big = linspace(x(1), x(end), extra_params.interpolate_resolution(1));
        y_big = linspace(y(1), y(end), extra_params.interpolate_resolution(2));
        [xx_big, yy_big] = meshgrid(x_big, y_big);
    else
        x_big = x; y_big = y;
        xx_big = xx; yy_big = yy;
    end
    
    % check if dpi is found
    if ~isfield(extra_params, "dpi")
        extra_params.dpi = 200;
    end
    
    myVideo = VideoWriter(fname);
    myVideo.FrameRate = extra_params.fps;
    myVideo.Quality = 100;
    open(myVideo);
    
    
    
    %timeskip = ceil(t_max / (real_time * 20 * (t(2) - t(1))));
    %file_len = length(t);
    
    u_info = h5info(fname, "/uu");
    u_size = double(u_info.Dataspace.Size);
    page_size = u_size;
    page_size(3) = double(1);
    
    % find minimum and maximum values
    if ~isfield(extra_params, "range_max") || ~isfield(extra_params, "range_min")
        fprintf("range_max or range_min not found in extra_params. Finding range...\n")
        extra_params.range_max = -Inf * ones( [1 u_size(2)] );
        extra_params.range_min = Inf * ones( [1 u_size(2)] );
        for frame_i = 1:length(t_indices)
            u = h5read(fname, "/uu", double([1 1 t_indices(frame_i)]), page_size);
            extra_params.range_max = max(extra_params.range_max, max(u));
            extra_params.range_min = min(extra_params.range_min, min(u));
        end
        fprintf("Range found:");
        extra_params.range_max
        extra_params.range_min
    end
    
    for frame_i = 1:length(t_indices)
        u = h5read(fname, "/uu", double([1 1 t_indices(frame_i)]), page_size);
        clf;
        for species = 1:length(u(1, :, 1))
            subplot(1, length(u(1, :, 1)), species)
            uu = reshape(u(:, species, 1), size(xx));
            
            if interp
                uu_big = interp2(xx,yy,uu,xx_big,yy_big,'spline');
            else
                uu_big = uu;
            end
            crange = [extra_params.range_min(species) extra_params.range_max(species)];
            hold on;
            ylim([min(y_big) max(y_big)])
            imagesc(x_big, y_big, uu_big, crange );
            title(sprintf("Chemical %i.\nTime = %.2f", species, t(floor(t_indices(frame_i)))));
            axis square;
            axis tight
            colorbar;
            hold off;
        end
        
        %frame = getframe(gcf);
        
        print('-dpng',strcat("-r", string(floor(extra_params.dpi))), strcat(input_foldername, "/temp"))
        
        I = imread(strcat(input_foldername, "/temp.png"));       % read saved image
        frame = im2frame(I);              % convert image to frame
        writeVideo(myVideo, frame);
        
    end
    close(myVideo)
end
