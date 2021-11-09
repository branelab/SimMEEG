function iV = idx3Dto1D(idx3D, dims, bMatlab)
% Convert a vector of 3 coordinate indecies to a one-dimentional index 
nx = dims(1); ny = dims(2); nz = dims(3);
ix = idx3D(1); iy = idx3D(2); iz = idx3D(3);

if bMatlab == 0
    iV = (ix - 1)*ny*nz + (iy - 1)*nz + iz;
else
    iV = (iz - 1)*ny*nx + (iy - 1)*nx + ix;
end
