function out = LTE_setclock(in,clock)

for ii = 1:length(in)
    in(ii).clock = clock;
end
out = in;