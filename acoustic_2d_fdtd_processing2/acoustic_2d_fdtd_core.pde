/* 2-D Acoustic FDTD Simulation Solver for Processing 2 (rev. Nov. 27th, 2014)
    by Yoshiki NAGATANI
      http://ultrasonics.jp/nagatani/e/
      https://twitter.com/nagataniyoshiki

  ==================================================================================
  [Variables]

   The followings should be defined as global variables in your main script.

   - dx: Spatical Resolution [m]
   - dt: Temporal Resolution [s]
   - dt_over_dx: dt/dx (for more efficient calculation)

   - NX: Spatial Size for X direction [pixels]
   - NY: Spatial Size for Y direction [pixels]

   - Vx: Particle Velocity for X direction [m/s] ==> 2-D Array: (NX+1 x NY)
   - Vy: Particle Velocity for Y direction [m/s] ==> 2-D Array: (NX x NY+1)
   - P:  Sound Pressure [Pa] ==> 2-D Array: (NX x NY)

   - rho:        Densities [kg/m^3] ==> 1-D Array: [Medium1, Medium2, ...]
   - kappa:      Bulk Moduli [Pa] ==> 1-D Array: [Medium1, Medium2, ...]
  ==================================================================================

  *** Technical Notes ***
   This program has been tested with Processing 2.2.1 (64bit) both on Windows and OSX.
   The speed of calculation depends on the CPU power of your PC. Intel Core i3 series
   or higher is recommended.
   For more detailed information about FDTD method (including 3-D elastic simulation),
   please refer our papers on simulation. ==> http://ultrasonics.jp/nagatani/fdtd/
   Thank you.

*/



/* Initialize Field Values */
void InitField(){
	for(int i=0;i<NX+1;i++){
		for(int j=0;j<NY;j++){
			Vx[i][j] = 0.0;
		}
	}
	for(int i=0;i<NX;i++){
		for(int j=0;j<NY+1;j++){
			Vy[i][j] = 0.0;
		}
	}
	for(int i=0;i<NX;i++){
		for(int j=0;j<NY;j++){
			P[i][j]  = 0.0;
		}
	}
}

/* Update Particle Velocities */
void UpdateV()
{
	for(int i=1;i<NX;i++){
		for(int j=0;j<NY;j++){
			Vx[i][j] += - dt_over_dx / ( (rho[Model[i][j]]+rho[Model[i-1][j]])/2.0 ) * ( P[i][j] - P[i-1][j] );
//			Vx[i][j] += - dt / ( dx * (rho[Model[i][j]]+rho[Model[i-1][j]])/2.0 ) * ( P[i][j] - P[i-1][j] );
		}
	}

	for(int i=0;i<NX;i++){
		for(int j=1;j<NY;j++){
			Vy[i][j] += - dt_over_dx / ( (rho[Model[i][j]]+rho[Model[i][j-1]])/2.0 ) * ( P[i][j] - P[i][j-1] );
		}
	}
}

/* Update Sound Pressure */
void UpdateP()
{
	for(int i=0;i<NX;i++){
		for(int j=0;j<NY;j++){
			P[i][j] += - ( kappa[Model[i][j]] * dt_over_dx )
			            * ( ( Vx[i+1][j] - Vx[i][j] ) + ( Vy[i][j+1] - Vy[i][j] ) );
//			P[i][j] += - ( kappa[Model[i][j]] * dt / dx )
//			            * ( ( Vx[i+1][j] - Vx[i][j] ) + ( Vy[i][j+1] - Vy[i][j] ) );
		}
	}
}
