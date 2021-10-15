function [idx,idy] = pos2grid(x,y,step)
idx = round(x/step)+1;
idy = round(y/step)+1;
end