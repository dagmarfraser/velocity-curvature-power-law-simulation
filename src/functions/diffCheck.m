
N = 9
m = 1

samplerate = 10;

x = 1:1000;

x = sin(x);

xdiff = diff(x)

[x0,x1,x2]=HaS_smoothDifferentiate(x,N,m,samplerate);
clf
plot(xdiff)
hold on;
plot(x1 / samplerate)