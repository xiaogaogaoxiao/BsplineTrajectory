function [x,y] = grid2pos(idx,idy,step)
x = (idx-1)*step;
y = (idy-1)*step;
end