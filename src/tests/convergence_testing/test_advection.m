system("nard2 new $(dirname $(whereis nard2 | awk '{print $2}'))/src/user_functions/adv_test.f90 _convtesting_adv");

Nx_range = 2.^(5:9);
foldername = sprintf("%s/test_adv", pwd());
mkdir(foldername);
if false
parfor i = 1:numel(Nx_range)
    conf = default_config_ADV(Nx_range(i));
    %conf.rparams(1) = dt_range(i);
    conf.iparams(1:8) = [1000 10000000000 1 1 0 1 1 1];
    % savenum, max_save_size, BCx, BCy, (BCz), timestepping_method, NAD?
    
    conf.savefilename = sprintf("%s/results_%i_adv.h5",foldername, i);
    conf.plotfilename = sprintf("%s/temp%i_adv.png", foldername, i);
    conf_filename = sprintf('%s/config_%i_adv.h5', foldername, i);
    write_config(conf, conf_filename, true); 
    nard2('_convtesting_adv', conf_filename);
end
end

errors = zeros(size(dt_range));
for i = 1:length(dt_range)
    
    savefilename = sprintf("%s/test_adv/results_%i_adv.h5",pwd(), i);
    x = h5read(savefilename, "/x");
    y = h5read(savefilename, "/y");
    uu = h5read(savefilename, "/uu");
    t = h5read(savefilename, "/t");
    
    [xx, yy] = meshgrid(x, y);
    integration_weights = xx * 0 + max(x) * max(y) / numel(xx);
    
    current_index = 1;
    for current_index = 1:1001
        current_time = t(current_index);
        if current_time < 0 
            break
        end
        current_u = uu(:, 1, current_index);
        exact_solution = cos(2 * pi * xx - 2 * pi * 0.05*current_time*yy);
        error_sq = ( exact_solution - reshape(current_u, size(xx)) ).^2;
        
        %{%
        figure(1)
        surf(xx, yy, reshape(current_u, size(xx)));
        drawnow;
        
        figure(2)
        surf(xx, yy, exact_solution);
        drawnow;
        %}
        
        errors(i) = max(sqrt( sum(error_sq .* integration_weights, 'all') ), errors(i));
        
    end
    loglog(dt_range, errors)
end

function [conf] = default_config_ADV(N)
    
    dt = 1/(512 * 64 * 0.05);

    conf.x = linspace(0, 1, N+1);
    conf.x = conf.x(2:end);
    conf.y = conf.x;
    conf.diffusion_consts = [0 ];
    conf.user_params = [-0.05];


    conf.rparams = zeros(64, 1);
    conf.rparams(1:3) = [dt 1 60.0];
    % dt, max time, plot_interval
    conf.iparams = zeros(64, 1);
    conf.iparams(1:8) = [1000 10000000000 1 1 0 2 1 1];
    % savenum, max_save_size, BCx, BCy, (BCz), timestepping_method,
    % non-autonomous diffusion? (1 for NA, 0 for A);
    % advection
    conf.iparams = int64(conf.iparams);

    conf.DBCx_plus = [0.0 0.0];
    conf.DBCx_minus = [0.0 0.0];
    conf.DBCy_plus = [0.0 0.0];
    conf.DBCy_minus = [0.0 0.0];

    conf.DBCx_plus_mask = int64([0 0]);
    conf.DBCx_minus_mask = int64([0 0]);
    conf.DBCy_plus_mask = int64([0 0]);
    conf.DBCy_minus_mask = int64([0 0]);

    [xx, yy] = meshgrid(conf.x, conf.y);

    temp = cos(2 * pi * xx);
    conf.IC(:, 1) = temp(:);
    
  
end