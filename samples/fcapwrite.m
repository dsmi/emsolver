function fcapwrite(fileName, tri, x, y, z, conductors);
% tswrite(fileName, freq, Y, t, r0)
%
% Save geometry in FastCap general format.
%

f = fopen(fileName, 'wt');

fprintf(f, '0 Exported from emsolver\n');

nc = length(conductors);
for i = 1:nc
    for j = conductors{i}
        x1 = x( tri(j, 1) );
        y1 = y( tri(j, 1) );
        z1 = z( tri(j, 1) );
        x2 = x( tri(j, 2) );
        y2 = y( tri(j, 2) );
        z2 = z( tri(j, 2) );
        x3 = x( tri(j, 3) );
        y3 = y( tri(j, 3) );
        z3 = z( tri(j, 3) );
        fspec = 'T %i %.11e %.11e %.11e %.11e %.11e %.11e %.11e %.11e %.11e\n';
        fprintf(f, fspec, i, x1, y1, z1, x2, y2, z2, x3, y3, z3);
    end
end

fclose(f);
