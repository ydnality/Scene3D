
%%calculates the new lookAt positions for a camera array

p1 = [-129.027725 -87.079384 80.930077   ];
p2 = [-82.804611 -61.262680 72.396408  ];
v = p2 - p1
vp = [v(2) -v(1) 0]
vp = vp./sqrt(sum(vp .* vp))

vVertical = [ 0 0 1];

horIncrement = 3
verIncrement = 0

newP1 = p1 + vp * horIncrement + vVertical * verIncrement
newP2 = newP1 + v

%% calculate new lookAt positions that are further back (for flash experiment)

p1 = [-56.914787 -105.385544 35.0148 ];
p2 = [-56.487434 -104.481461 34.8    ];
v = p2 - p1
v = v./sqrt(sum(v .* v))

depthIncrement = -50;

newP1 = p1 + depthIncrement .* v
newP2 = p2 + depthIncrement .* v



%% calculate new lookAt positions that are to the right (for flash experiment)

p1 = [-36.914787 -105.385544 13.014802       ];
p2 = [ -36.914787 -104.481461 12.95     ];
v = p2 - p1
vp = [v(2) -v(1) 0]
vp = vp./sqrt(sum(vp .* vp))

vVertical = [ 0 0 1];

horIncrement = 25
verIncrement = 0

newP1 = p1 + vp * horIncrement + vVertical * verIncrement
newP2 = newP1 + v
