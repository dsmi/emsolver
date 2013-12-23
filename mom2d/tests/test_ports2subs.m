function test_ports2subs()

ports = { [ 2 4 ] [ 8 9 10 ] };

[ faceidx portidx ] = ports2subs(ports);

test_faceidx = [ 2 4 8 9 10 ];
test_portidx = [ 1 1 2 2 2 ];

assertEquals(test_faceidx, faceidx);
assertEquals(test_portidx, portidx);
