function X = testderX(sole,x,sX)
sole.coor = x;
sole.stiffness();
X = sole.(sX);
end