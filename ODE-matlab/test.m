mu = @(t) (square(t, 50)+1)/2;
mu(10)

%%
t = linspace(-pi,2*pi,300);
x = 1.15*square(2*t);

plot(t/pi,x,'.-',t/pi,1.15*sin(2*t))
xlabel('t / \pi')
grid on

%%
start_time = 0;
end_time = 20;
precision = 0.1; % lower = more precise
total_points = (end_time - start_time) / precision;

time = linspace(start_time, end_time, total_points);
mu = ~mod(floor(time./(pi)),2)