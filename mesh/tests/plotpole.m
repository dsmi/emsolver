
addpath("..");

[ tri, x, y, z ] = mkpole(2, 1, 4, 5, 2, 0.5);

trisurf(tri, x, y, z);
xlabel('X');
ylabel('Y');
zlabel('Z'); 

hold on;

[ nx, ny, nz, cx, cy, cz ] = calc_triangle_normals(tri, x, y, z);

quiver3(cx,cy,cz,nx,ny,nz);

hold off;
