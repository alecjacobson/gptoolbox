function interactive_fit(filename,x_name,y_name,skip_interaction)
  % function interactive_fit(filename,x_name,y_name,skip_interaction)
  % Reads in plot X and Y data from filename stored in variables x_name, and
  % y_name respectively. Then interactively lets the user pick axes, slope and
  % y_intercept to match the data. Those variables are then saved in the file.
  % load the data from the .mat file
  load(filename,x_name,y_name);
  x = eval(x_name);
  y = eval(y_name);

  % none of these should exist or they will be overwritten
  if( ...
    exist('x_max') | ...
    exist('y_max') | ...
    exist('x_min') | ...
    exist('y_min') | ...
    exist('m') | ...
    exist('b') | ...
    exist('x_label') | ...
    exist('y_label') | ...
    exist('print_name') | ...
    exist('log_base'))
    error(['Error: x_max,  y_max,  x_min, y_min,x_label,y_label, m, b, and log_base' ...
      ' must be undefined.']);
  end


  % load in axis min and maxes, and line fit info if they exist
  warning off MATLAB:load:variableNotFound;
  load(filename,'x_max','y_max','x_min','y_min','m','b','log_base','x_label','y_label','print_name');
  warning on MATLAB:load:variableNotFound;
  %fprintf('OVERRIDE: \n');
  %x_min = -2.5;
  %x_max= -0.5;
  %y_min = -7;
  %y_max = -4;

  plot_figure = figure;

  if(~exist('skip_interaction') || ~skip_interaction)
    % prompt for whether plot is log or not
    if(~exist('log_base'))
      %log_base = str2double(input('Log plot? [Enter base or No]: ', 's'));
      log_base = str2double(input('Log plot? [10]: ', 's'));

      % ignore nonsense
      if(isnan(log_base) || log_base<=1)
        %log_base = false;
        log_base = 10;
      end
    end

  end

    if(log_base)
      fprintf('Log-log base is %d\n',log_base);
      % take log of x and y
      x = log(x)/log(log_base);
      y = log(y)/log(log_base);
    end

    
    % auto-set axes if not already existing
  if(~exist('skip_interaction') || ~skip_interaction)
    if( ...
      ~exist('x_max') | ...
      ~exist('y_max') | ...
      ~exist('x_min') | ...
      ~exist('y_min'))
      x_min = floor(min(x));
      y_min = floor(min(y));
      x_max = ceil(max(x));
      y_max = ceil(max(y));
      %x_min = (min(x));
      %y_min = (min(y));
      %x_max = (max(x));
      %y_max = (max(y));
    end

    if(~exist('x_label'))
      %x_label = '';
      x_label = 'h';
    end

    %fprintf('OVERRIDE: x_label = "h"\n');
    %x_label = 'h';

    if(~exist('y_label'))
      %y_label = '';
      y_label = 'error';
    end


    figure(plot_figure);
    interactive_fit_plot(x,y,x_min,x_max,y_min,y_max,x_label,y_label,log_base);

    %fprintf('Current axes are: %f < x < %f, %f < y < %f\n', ...
    %  x_min,x_max,y_min,y_max
    reply = str2double(input(['New x_min (' num2str(x_min) '): '], 's'));
    if(~isnan(reply))
      x_min = reply;
    end

    figure(plot_figure);
    interactive_fit_plot(x,y,x_min,x_max,y_min,y_max,x_label,y_label,log_base);

    reply = str2double(input(['New x_max (' num2str(x_max) '): '], 's'));
    if(~isnan(reply))
      x_max= reply;
    end

    figure(plot_figure);
    interactive_fit_plot(x,y,x_min,x_max,y_min,y_max,x_label,y_label,log_base);

    reply = str2double(input(['New y_min (' num2str(y_min) '): '], 's'));
    if(~isnan(reply))
      y_min = reply;
    end

    figure(plot_figure);
    interactive_fit_plot(x,y,x_min,x_max,y_min,y_max,x_label,y_label,log_base);

    reply = str2double(input(['New y_max (' num2str(y_max) '): '], 's'));
    if(~isnan(reply))
      y_max = reply;
    end

    figure(plot_figure);
    interactive_fit_plot(x,y,x_min,x_max,y_min,y_max,x_label,y_label,log_base);

    reply = input(['New x_label ("' x_label '"): '], 's');
    if(~isempty(reply))
      x_label = reply;
    end
    reply = input(['New y_label ("' y_label '"): '], 's');
    if(~isempty(reply))
      y_label = reply;
    end

    figure(plot_figure);
    interactive_fit_plot(x,y,x_min,x_max,y_min,y_max,x_label,y_label,log_base);

    % if slope and intercept are not present then solve for them
    if(~exist('m') || ~exist('b'))
      linear_fit = polyfit(x,y,1);
      b = linear_fit(2);
      m = linear_fit(1);
    end

    figure(plot_figure);
    interactive_fit_plot(x,y,x_min,x_max,y_min,y_max,x_label,y_label,log_base,m,b);

    % keep prompting for a new slope and y-intercept until user specifies
    % something other than a number or nothing for both
    while(true)
      reply = str2double(input(['New slope (' num2str(m) '): '], 's'));
      if(~isnan(reply))
        m = reply;
        nan_m = false;
      else
        nan_m = true;
      end
      figure(plot_figure);
      interactive_fit_plot(x,y,x_min,x_max,y_min,y_max,x_label,y_label,log_base,m,b);

      reply = str2double(input(['New y-intercept (' num2str(b) '): '], 's'));
      if(~isnan(reply))
        b = reply;
      elseif(nan_m)
        % break loop if user has answer nothing for both m and b
        break;
      end
      figure(plot_figure);
      interactive_fit_plot(x,y,x_min,x_max,y_min,y_max,x_label,y_label,log_base,m,b);
    end


    if(~exist('print_name'))
      print_name = [filename '.eps'];
    end
    reply = input(['New print_name ("' print_name '"): '], 's');
    if(~isempty(reply))
      print_name = reply;
    end

  end

  figure(plot_figure);
  interactive_fit_plot(x,y,x_min,x_max,y_min,y_max,x_label,y_label,log_base,m,b);
  fprintf('Plot printed to %s\n',print_name);
  print('-depsc',print_name);
  close(plot_figure);
  % append new variables to file with data left untouched
  save(filename,'log_base','x_max','x_min','y_max','y_min','x_label','y_label','print_name','m','b','-append');

end
