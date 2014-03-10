function [t, measurement_overhead, measurement_details] = timeit(f, num_outputs)
%TIMEIT Measure time required to run function.
%   T = TIMEIT(F) measures the time (in seconds) required to run F, which is a
%   function handle.  TIMEIT calls F with either no output arguments or one
%   output argument depending on nargout(F).
%
%   T = TIMEIT(F,N) calls F with N output arguments.  N can be 0, 1, 2, 3, 4, or 5.
%
%   TIMEIT handles automatically the usual benchmarking procedures of "warming
%   up" F, figuring out how many times to repeat F in a timing loop, etc.
%   TIMEIT also compensates for the estimated time-measurement overhead
%   associated with tic/toc and with calling function handles.  TIMEIT returns
%   the median of several repeated measurements.
%
%   Examples
%   --------
%   How much time does it take to compute sum(A.' .* B, 1), where A is
%   12000-by-400 and B is 400-by-12000?
%
%       A = rand(12000, 400);
%       B = rand(400, 12000);
%       f = @() sum(A.' .* B, 1);
%       timeit(f)
%
%   How much time does it take to call svd with three output arguments?
%
%       X = [1 2; 3 4; 5 6; 7 8];
%       f = @() svd(X);
%       timeit(f, 3)
%
%   How much time does it take to dilate the text.png image with
%   a 25-by-25 all-ones structuring element? (This example uses Image Processing
%   Toolbox functions.)
%
%       bw = imread('text.png');
%       se = strel(ones(25, 25));
%       g = @() imdilate(bw, se);
%       timeit(g)

%   Steve Eddins
%   Copyright 2008-2010 The MathWorks, Inc.

if nargin < 2
    num_outputs = min(numOutputs(f), 1);
else
    if num_outputs > 5
        warning('MATLAB:timeit:tooManyOutputs', ...
            'Too many function output arguments specified. timeit will call your function with 5 output arguments.');
    end
end

t_rough = roughEstimate(f, num_outputs);

% Calculate the number of inner-loop repetitions so that the inner for-loop
% takes at least about 1ms to execute.
desired_inner_loop_time = 0.001;
num_inner_iterations = max(ceil(desired_inner_loop_time / t_rough), 1);

% Run the outer loop enough times to give a reasonable set of inputs to median.
num_outer_iterations = 11;

% If the estimated running time for the timing loops is too long,
% reduce the number of outer loop iterations.
estimated_running_time = num_outer_iterations * num_inner_iterations * t_rough;
long_time = 15;
min_outer_iterations = 3;
if estimated_running_time > long_time
    num_outer_iterations = ceil(long_time / (num_inner_iterations * t_rough));
    num_outer_iterations = max(num_outer_iterations, min_outer_iterations);
end

times = zeros(num_outer_iterations, 1);

for k = 1:num_outer_iterations
    % Coding note: An earlier version of this code constructed an "outputs" cell
    % array, which was used in comma-separated form for the left-hand side of
    % the call to f().  It turned out, though, that the comma-separated output
    % argument added significant measurement overhead.  Therefore, the cases
    % for different numbers of output arguments are hard-coded into the switch
    % statement below.
    switch num_outputs
        case 0
            tic();
            for p = 1:num_inner_iterations
                f();
            end
            times(k) = toc();
            
        case 1
            tic();
            for p = 1:num_inner_iterations
                output = f();
            end
            times(k) = toc();
            
        case 2
            tic();
            for p = 1:num_inner_iterations
                [output1, output2] = f();
            end
            times(k) = toc();
            
        case 3
            tic();
            for p = 1:num_inner_iterations
                [output1, output2, output3] = f();
            end
            times(k) = toc();
            
        case 4
            tic();
            for p = 1:num_inner_iterations
                [output1, output2, output3, output4] = f();
            end
            times(k) = toc();
            
        otherwise
            tic();
            for p = 1:num_inner_iterations
                [output1, output2, output3, output4, output5] = f();
            end
            times(k) = toc();
    end

end

t = median(times) / num_inner_iterations;

measurement_details.EmptyFunctionCallTime = emptyFunctionCallTime();
measurement_details.SimpleFunctionHandleCallTime = simpleFunctionHandleCallTime();
measurement_details.AnonymousFunctionHandleCallTime = anonFunctionHandleCallTime();
measurement_details.TicTocCallTime = tictocCallTime();
measurement_overhead = (tictocCallTime() / num_inner_iterations) + ...
    functionHandleCallOverhead(f);

t = max(t - measurement_overhead, 0);

if t < (5 * measurement_overhead)
    warning('MATLAB:timeit:HighOverhead', 'The measured time for F may be inaccurate because it is close to the estimated time-measurement overhead (%.1e seconds).  Try measuring something that takes longer.', measurement_overhead);
end

function t = roughEstimate(f, num_f_outputs)
%   Return rough estimate of time required for one execution of
%   f().  Basic warmups are done, but no fancy looping, medians,
%   etc.

% Warm up tic/toc.
tic();
elapsed = toc();
tic();
elapsed = toc();

% Call f() in a loop for at least a millisecond.
times = [];
time_threshold = 3;
iter_count = 0;
while sum(times) < 0.001
    iter_count = iter_count + 1;
    
    switch num_f_outputs
        case 0
            tic();
            f();
            times(end+1) = toc();
            
        case 1
            tic();
            output1 = f();
            times(end+1) = toc();
            
        case 2
            tic();
            [output1, output2] = f();
            times(end+1) = toc();
            
        case 3
            tic();
            [output1, output2, output3] = f();
            times(end+1) = toc();
            
        case 4
            tic();
            [output1, output2, output3, output4] = f();
            times(end+1) = toc();
            
        otherwise
            tic();
            [output1, output2, output3, output4, output5] = f();
            times(end+1) = toc();
    end
    
    if iter_count == 1
        if times > time_threshold
            % If the first call to f() takes more than time_threshold to run,
            % then just use the result from that call.  The assumption is that
            % first-time effects are negligible compared to the running time for
            % f().
            break;
        else
            % Discard first timing.
            times = [];
        end
    end
end

t = median(times);
        
function n = numOutputs(f)
%   Return the number of output arguments to be used when calling the function
%   handle f.  
%   * If nargout(f) > 0, return 1.
%   * If nargout(f) == 0, return 0.
%   * If nargout(f) < 0, use try/catch to determine whether to call f with one
%     or zero output arguments.
%     Note: It is not documented (as of R2008b) that nargout can return -1.
%     However, it appears to do so for functions that use varargout and for
%     anonymous function handles.  

n = nargout(f);
if n < 0
   try
      a = f();
      % If the line above doesn't throw an error, then it's OK to call f() with
      % one output argument.
      n = 1;
      
   catch %#ok<CTCH>
      % If we get here, assume it's because f() has zero output arguments.  In
      % recent versions of MATLAB we could catch the specific exception ID
      % MATLAB:maxlhs, but that would limit the use of timeit to MATLAB versions
      % since the introduction of MExceptions.
      n = 0;
   end
end

function t = tictocCallTime
% Return the estimated time required to call tic/toc.

persistent ttct
if ~isempty(ttct)
    t = ttct;
    return
end

% Warm up tic/toc.
temp = tic(); elapsed = toc();
temp = tic(); elapsed = toc();
temp = tic(); elapsed = toc();

num_repeats = 11;
times = zeros(1, num_repeats);

for k = 1:num_repeats
   times(k) = tictocTimeExperiment();
end

t = min(times);
ttct = t;

function t = tictocTimeExperiment
% Call tic/toc 100 times and return the average time required.

elapsed = 0;
% Call tic/toc 100 times.
tic(); elapsed = elapsed +  toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();
tic(); elapsed = elapsed + toc();

t = elapsed / 100;

function emptyFunction()

function t = functionHandleCallOverhead(f)
% Return the estimated overhead, in seconds, for calling a function handle
% compared to calling a normal function.

fcns = functions(f);
if strcmp(fcns.type, 'anonymous')
    t = anonFunctionHandleCallTime();
else
    t = simpleFunctionHandleCallTime();
end

t = max(t - emptyFunctionCallTime(), 0);

function t = simpleFunctionHandleCallTime
% Return the estimated time required to call a simple function handle to a
% function with an empty body.
%
% A simple function handle fh has the form @foo.

persistent sfhct
if ~isempty(sfhct)
    t = sfhct;
    return
end

num_repeats = 101;
% num_repeats chosen to take about 100 ms, assuming that
% timeFunctionHandleCall() takes about 1 ms.
times = zeros(1, num_repeats);

fh = @emptyFunction;

% Warm up fh().
fh();
fh();
fh();

for k = 1:num_repeats
   times(k) = functionHandleTimeExperiment(fh);
end

t = min(times);
sfhct = t;

function t = anonFunctionHandleCallTime
% Return the estimated time required to call an anonymous function handle that
% calls a function with an empty body.
%
% An anonymous function handle fh has the form @(arg_list) expression. For
% example:
%
%       fh = @(thetad) sin(thetad * pi / 180)

persistent afhct
if ~isempty(afhct)
    t = afhct;
    return
end

num_repeats = 101;
% num_repeats chosen to take about 100 ms, assuming that timeFunctionCall()
% takes about 1 ms.
times = zeros(1, num_repeats);

fh = @() emptyFunction();

% Warm up fh().
fh();
fh();
fh();

for k = 1:num_repeats
   times(k) = functionHandleTimeExperiment(fh);
end

t = min(times);
afhct = t;

function t = functionHandleTimeExperiment(fh)
% Call the function handle fh 2000 times and return the average time required.

% Record starting time.
tic();

fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();

fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();

fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();

fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();

fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();

fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();

fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();

fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();

fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();

fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();

fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();

fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();

fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();

fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();

fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();

fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();

fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();

fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();

fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();

fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();
fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh(); fh();

t = toc() / 2000;

function t = emptyFunctionCallTime()
% Return the estimated time required to call a function with an empty body.

persistent efct
if ~isempty(efct)
    t = efct;
    return
end

% Warm up emptyFunction.
emptyFunction();
emptyFunction();
emptyFunction();

num_repeats = 101;
% num_repeats chosen to take about 100 ms, assuming that timeFunctionCall()
% takes about 1 ms.
times = zeros(1, num_repeats);

for k = 1:num_repeats
   times(k) = emptyFunctionTimeExperiment();
end

t = min(times);
efct = t;

function t = emptyFunctionTimeExperiment()
% Call emptyFunction() 2000 times and return the average time required.

% Record starting time.
tic();

emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();

emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();

emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();

emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();

emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();

emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();

emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();

emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();

emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();

emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();

emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();

emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();

emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();

emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();

emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();

emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();

emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();

emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();

emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();

emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();

emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();
emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction(); emptyFunction();

t = toc() / 2000;

% Copyright (c) 2010, The MathWorks, Inc.
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
%     * Neither the name of the The MathWorks, Inc. nor the names
%       of its contributors may be used to endorse or promote products derived
%       from this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
