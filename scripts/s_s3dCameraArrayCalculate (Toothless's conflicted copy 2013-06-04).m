
%calculates the new lookAt positions for a camera array

p1 = [359.4034 -341.6905 28 ];
p2 = [357.9510 -340.3236 28 ];
v = p2 - p1
vp = [v(2) -v(1) 0]

vVertical = [ 0 0 1];

horIncrement = 2
verIncrement = 2

newP1 = p1 + vp * horIncrement + vVertical * verIncrement
newP2 = newP1 + v