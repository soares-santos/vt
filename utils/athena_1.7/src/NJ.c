#include <stdlib.h>
#include <stdio.h>
#include <math.h>

typedef unsigned int uint;

/* Calculates Nx, Ny with Nx * Ny = N, and Ny/Nx = y/x = r */
void Nxy_jackknife(uint N, double r, uint *Nx, uint *Ny)
{
   double dNx, dNy;

   dNx = sqrt((double)N / r);
   dNy = sqrt((double)N * r);

   *Nx  = (uint)round(dNx);
   *Ny  = (uint)round(dNy);

   if (*Nx == 0) (*Nx) ++;
   if (*Ny == 0) (*Ny) ++;

   printf("%g %g %u %u\n", dNx, dNy, *Nx, *Ny);
}

int main(int argc, char *argv[])
{
   uint N, Nx, Ny;
   double r;

   N = atoi(argv[1]);
   r = atof(argv[2]);

   Nxy_jackknife(N, r, &Nx, &Ny);

   return 0;
}
