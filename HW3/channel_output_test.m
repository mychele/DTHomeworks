function [  ] = channel_output_test( x, T, Tc)

Q0 = T/Tc; % interpolation factor
q = [0,0,0,0,0,0,0,0,0.19*exp(-1i*2.21), 0.09*exp(1i*1.64), 0.7*exp(-1i*2.57), ...
    0.45, 0.6*exp(-1i*2.26), 0.35*exp(1i*3.15), 0.24*exp(1i*1.34), 0.37*exp(1i*2.6), ...
    0.34*exp(-1i*1.17), 0, 0.15*exp(-1i*2.66), 0.15*exp(1i*3.27), 0.17*exp(1i*2.13), ...
    0.4*exp(1i*2.06), 0.58*exp(-1i*1.51), 0.03*exp(1i*2.15), 0.18*exp(1i*3.6), ...
    0.29*exp(1i*3.17), 0.4*exp(-1i*1.63), 0.07*exp(-1i*3.16)];

out = channel_output(x, T, Tc, 100000);
out_control = conv(interp_nofilter(x, Q0), q);

stem(abs(out))
hold on, stem(abs(out_control), 'x')

%errsquared = sum(abs(out-out_control).^2);

end