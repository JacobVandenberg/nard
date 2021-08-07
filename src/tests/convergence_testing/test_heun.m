system("nard2 new $(dirname $(whereis nard2 | awk '{print $2}'))/src/user_functions/heun_test.f90 _convtesting_heun");

dt_range = 2.^(-4:-1:-9);
foldername = sprintf("%s/test_heun", pwd());
mkdir(foldername);
parfor i = 1:numel(dt_range)
    conf = default_config_CNNAD;
    conf.rparams(1) = dt_range(i);
    conf.iparams(1:7) = [1000 10000000000 1 1 0 2 1];
    % savenum, max_save_size, BCx, BCy, (BCz), timestepping_method, NAD?
    
    conf.savefilename = sprintf("%s/results_%i_heun.h5",foldername, i);
    conf.plotfilename = sprintf("%s/temp%i_heun.png", foldername, i);
    conf_filename = sprintf('%s/config_%i_heun.h5', foldername, i);
    write_config(conf, conf_filename, true);
    nard2('_convtesting_heun', conf_filename);
end

errors = zeros(size(dt_range));
for i = 1:length(dt_range)
    
    savefilename = sprintf("%s/test_heun/results_%i_heun.h5",pwd(), i);
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
        exact_solution = xx*0 + sqrt( 2/3 * current_time^3 + 1 );
        error_sq = ( exact_solution - reshape(current_u, size(xx)) ).^2;
        errors(i) = max(sqrt( sum(error_sq .* integration_weights, 'all') ), errors(i));
        
    end
    loglog(dt_range, errors)
end

function [conf] = default_config_CNNAD()
    
    dt = 0.0001;
    N = 32;

    conf.x = linspace(-1, 1, N+1);
    conf.x = conf.x(2:end);
    conf.y = conf.x;
    conf.diffusion_consts = [1.0 ];
    conf.user_params = [0];


    conf.rparams = zeros(64, 1);
    conf.rparams(1:3) = [dt 1.0 60.0];
    % dt, max time, plot_interval
    conf.iparams = zeros(64, 1);
    conf.iparams(1:7) = [1000 10000000000 1 1 0 1 1];
    % savenum, max_save_size, BCx, BCy, (BCz), timestepping_method,
    % non-autonomous diffusion? (1 for NA, 0 for A)
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

    temp = 1 + xx*0;
    conf.IC(:, 1) = temp(:);
    
  
end