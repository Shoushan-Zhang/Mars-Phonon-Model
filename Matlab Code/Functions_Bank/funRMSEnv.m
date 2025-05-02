function rms_env = funRMSEnv(signal, winLen_num)

% FUNRMSENV Compute the RMS envelope of a signal using a moving window.
%
%   rms_env = funRMSEnv(signal, winLen_num)
%
%   Inputs:
%       signal      - Input signal vector (1D array)
%       winLen_num  - Window length in number of samples (positive scalar)
%
%   Output:
%       rms_env     - RMS envelope of the input signal

    % Input validation
    if nargin < 2
        error('Two inputs are required: signal and winLen_num.');
    end

    if ~isvector(signal)
        error('Input signal must be a 1D array.');
    end

    if ~isscalar(winLen_num) || winLen_num <= 0
        error('Window length must be a positive scalar (number of samples).');
    end

    if winLen_num > length(signal)
        error('Window length exceeds signal length.');
    end

    % Ensure signal is a column vector
    signal = signal(:);

    % Square the signal
    squared_signal = signal .^ 2;

    % Apply moving average (smoothing)
    smoothed_signal = movmean(squared_signal, winLen_num);

    % Take the square root
    rms_env = sqrt(smoothed_signal);
end
