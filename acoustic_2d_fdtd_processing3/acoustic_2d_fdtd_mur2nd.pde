/* Mur's 2nd order absorbing boundary for 2-D Acoustic FDTD (rev. Nov. 27th, 2014)
    by Yoshiki NAGATANI
      http://ultrasonics.jp/nagatani/e/
      https://twitter.com/nagataniyoshiki

  ==================================================================================
  [Variables]

   The followings should be defined as global variables in your main script.

   - dx: Spatical Resolution [m]
   - dt: Temporal Resolution [s]

   - NX: Spatial Size for X direction [pixels]
   - NY: Spatial Size for Y direction [pixels]

   - P:  Sound Pressure [Pa] ==> 2-D Array: (NX x NY)

   - rho:   Densities [kg/m^3] ==> 1-D Array: [Medium1, Medium2, ...]
   - kappa: Bulk Moduli [Pa] ==> 1-D Array: [Medium1, Medium2, ...]

   - Mur_X1: Mur's 1st order absorbing boundary ==> 2-D Array: (4 x NY)
   - Mur_X2: Mur's 2nd order absorbing boundary ==> 2-D Array: (4 x NY)
   - Mur_Y1: Mur's 1st order absorbing boundary ==> 2-D Array: (NX x 4)
   - Mur_Y2: Mur's 2nd order absorbing boundary ==> 2-D Array: (NX x 4)
  ==================================================================================

  *** Technical Notes ***
   This program has been tested with Processing 2.2.1 and 3.1.2 (64bit) both on Windows and OSX.
   The speed of calculation depends on the CPU power of your PC. Intel Core i3 series
   or higher is recommended.
   For more detailed information about FDTD method (including 3-D elastic simulation),
   please refer our papers on simulation. ==> http://ultrasonics.jp/nagatani/fdtd/
   Thank you.

*/



/* Initialize Mur's 2nd Order Absorption Values */
void InitMur2nd(){
	for(int i=0;i<NX;i++){
		Mur_Y1[i][0] = 0.0;
		Mur_Y1[i][1] = 0.0;
		Mur_Y1[i][2] = 0.0;
		Mur_Y1[i][3] = 0.0;
		Mur_Y2[i][0] = 0.0;
		Mur_Y2[i][1] = 0.0;
		Mur_Y2[i][2] = 0.0;
		Mur_Y2[i][3] = 0.0;
	}
	for(int j=0;j<NY;j++){
		Mur_X1[0][j] = 0.0;
		Mur_X1[1][j] = 0.0;
		Mur_X1[2][j] = 0.0;
		Mur_X1[3][j] = 0.0;
		Mur_X2[0][j] = 0.0;
		Mur_X2[1][j] = 0.0;
		Mur_X2[2][j] = 0.0;
		Mur_X2[3][j] = 0.0;
	}
}

/* Mur's 2nd Order Absorption */
void Mur2nd()
{
	float vel;	// Wave velocity at the boundary
	int i,j;

	/* Mur's 2nd Order Absorption */
	for(i=2;i<NX-2;i++){
		vel = sqrt(kappa[Model[i][0]]/rho[Model[i][0]]);
		P[i][0] = - Mur_Y2[i][1]
		          + (vel*dt-dx)/(vel*dt+dx) * ( P[i][1] + Mur_Y2[i][0] )
		          + (2.0*dx)/(vel*dt+dx) * ( Mur_Y1[i][0] + Mur_Y1[i][1] )
		          + (dx*vel*vel*dt*dt)/(2.0*dx*dx*(vel*dt+dx))
		             * (   Mur_Y1[i+1][0] - 2.0 * Mur_Y1[i][0]
		                 + Mur_Y1[i-1][0] + Mur_Y1[i+1][1]
		                 - 2.0 * Mur_Y1[i][1] + Mur_Y1[i-1][1] );
		vel = sqrt(kappa[Model[i][NY-1]]/rho[Model[i][NY-1]]);
		P[i][NY-1] = - Mur_Y2[i][2]
		          + (vel*dt-dx)/(vel*dt+dx) * ( P[i][NY-2] + Mur_Y2[i][3] )
		          + (2.0*dx)/(vel*dt+dx) * ( Mur_Y1[i][3] + Mur_Y1[i][2] )
		          + (dx*vel*vel*dt*dt)/(2.0*dx*dx*(vel*dt+dx))
		             * (   Mur_Y1[i+1][3] - 2.0 * Mur_Y1[i][3]
		                 + Mur_Y1[i-1][3] + Mur_Y1[i+1][2]
		                 - 2.0 * Mur_Y1[i][2] + Mur_Y1[i-1][2] );
	}
	for(j=2;j<NY-2;j++){
		vel = sqrt(kappa[Model[0][j]]/rho[Model[0][j]]);
		P[0][j] = - Mur_X2[1][j]
		          + (vel*dt-dx)/(vel*dt+dx) * ( P[1][j] + Mur_X2[0][j] )
		          + (2.0*dx)/(vel*dt+dx) * ( Mur_X1[0][j] + Mur_X1[1][j] )
		          + (dx*vel*vel*dt*dt)/(2.0*dx*dx*(vel*dt+dx))
		             * (   Mur_X1[0][j+1] - 2.0 * Mur_X1[0][j]
		                 + Mur_X1[0][j-1] + Mur_X1[1][j+1]
		                 - 2.0 * Mur_X1[1][j] + Mur_X1[1][j-1] );
		vel = sqrt(kappa[Model[NX-1][j]]/rho[Model[NX-1][j]]);
		P[NX-1][j] = - Mur_X2[2][j]
		          + (vel*dt-dx)/(vel*dt+dx) * ( P[NX-2][j] + Mur_X2[3][j] )
		          + (2.0*dx)/(vel*dt+dx) * ( Mur_X1[3][j] + Mur_X1[2][j] )
		          + (dx*vel*vel*dt*dt)/(2.0*dx*dx*(vel*dt+dx))
		             * (   Mur_X1[3][j+1] - 2.0 * Mur_X1[3][j]
		                 + Mur_X1[3][j-1] + Mur_X1[2][j+1]
		                 - 2.0 * Mur_X1[2][j] + Mur_X1[2][j-1] );
	}

	/* Mur's 1st Order Absorption for 4 corners*/
	i = 1;
		vel = sqrt(kappa[Model[i][0]]/rho[Model[i][0]]);
		P[i][0] = Mur_Y1[i][1] + (vel*dt-dx)/(vel*dt+dx) * (P[i][1] - Mur_Y1[i][0]);
		vel = sqrt(kappa[Model[i][NY-1]]/rho[Model[i][NY-1]]);
		P[i][NY-1] = Mur_Y1[i][2] + (vel*dt-dx)/(vel*dt+dx) * (P[i][NY-2] - Mur_Y1[i][3]);
	i = NX-2;
		vel = sqrt(kappa[Model[i][0]]/rho[Model[i][0]]);
		P[i][0] = Mur_Y1[i][1] + (vel*dt-dx)/(vel*dt+dx) * (P[i][1] - Mur_Y1[i][0]);
		vel = sqrt(kappa[Model[i][NY-1]]/rho[Model[i][NY-1]]);
		P[i][NY-1] = Mur_Y1[i][2] + (vel*dt-dx)/(vel*dt+dx) * (P[i][NY-2] - Mur_Y1[i][3]);
	j = 1;
		vel = sqrt(kappa[Model[0][j]]/rho[Model[0][j]]);
		P[0][j] = Mur_X1[1][j] + (vel*dt-dx)/(vel*dt+dx) * (P[1][j] - Mur_X1[0][j]);
		vel = sqrt(kappa[Model[NX-1][j]]/rho[Model[NX-1][j]]);
		P[NX-1][j] = Mur_X1[2][j] + (vel*dt-dx)/(vel*dt+dx) * (P[NX-2][j] - Mur_X1[3][j]);
	j = NY - 2;
		vel = sqrt(kappa[Model[0][j]]/rho[Model[0][j]]);
		P[0][j] = Mur_X1[1][j] + (vel*dt-dx)/(vel*dt+dx) * (P[1][j] - Mur_X1[0][j]);
		vel = sqrt(kappa[Model[NX-1][j]]/rho[Model[NX-1][j]]);
		P[NX-1][j] = Mur_X1[2][j] + (vel*dt-dx)/(vel*dt+dx) * (P[NX-2][j] - Mur_X1[3][j]);

	/* Copy Previous Values */
	Mur2ndCopy();
}

/* Copy the Filed Values for Mur's 2nd Order Absorption */
void Mur2ndCopy()
{
	for(int i=0;i<NX;i++){
		/* Copy 1st Old Values to 2nd Old Values*/
		Mur_Y2[i][0] = Mur_Y1[i][0];
		Mur_Y2[i][1] = Mur_Y1[i][1];
		Mur_Y2[i][2] = Mur_Y1[i][2];
		Mur_Y2[i][3] = Mur_Y1[i][3];

		/* Copy Present Values */
		Mur_Y1[i][0] = P[i][0];
		Mur_Y1[i][1] = P[i][1];
		Mur_Y1[i][2] = P[i][NY-2];
		Mur_Y1[i][3] = P[i][NY-1];
	}
	for(int j=0;j<NY;j++){
		/* Copy 1st Old Values to 2nd Old Values*/
		Mur_X2[0][j] = Mur_X1[0][j];
		Mur_X2[1][j] = Mur_X1[1][j];
		Mur_X2[2][j] = Mur_X1[2][j];
		Mur_X2[3][j] = Mur_X1[3][j];

		/* Copy Present Values */
		Mur_X1[0][j] = P[0][j];
		Mur_X1[1][j] = P[1][j];
		Mur_X1[2][j] = P[NX-2][j];
		Mur_X1[3][j] = P[NX-1][j];
	}
}

/* Initialize Mur's 1st Order Absorption Values */
void InitMur1st(){
	for(int i=0;i<NX;i++){
		Mur_Y1[i][0] = 0.0;
		Mur_Y1[i][1] = 0.0;
		Mur_Y1[i][2] = 0.0;
		Mur_Y1[i][3] = 0.0;
	}
	for(int j=0;j<NY;j++){
		Mur_X1[0][j] = 0.0;
		Mur_X1[1][j] = 0.0;
		Mur_X1[2][j] = 0.0;
		Mur_X1[3][j] = 0.0;
	}
}

/* Mur's 1st Order Absorption */
void Mur1st()
{
	float vel;

	for(int i=1;i<NX-1;i++){
		vel = sqrt(kappa[Model[i][0]]/rho[Model[i][0]]);
		P[i][0] = Mur_Y1[i][1] + (vel*dt-dx)/(vel*dt+dx) * (P[i][1] - Mur_Y1[i][0]);
		vel = sqrt(kappa[Model[i][NY-1]]/rho[Model[i][NY-1]]);
		P[i][NY-1] = Mur_Y1[i][2] + (vel*dt-dx)/(vel*dt+dx) * (P[i][NY-2] - Mur_Y1[i][3]);
	}
	for(int j=1;j<NY-1;j++){
		vel = sqrt(kappa[Model[0][j]]/rho[Model[0][j]]);
		P[0][j] = Mur_X1[1][j] + (vel*dt-dx)/(vel*dt+dx) * (P[1][j] - Mur_X1[0][j]);
		vel = sqrt(kappa[Model[NX-1][j]]/rho[Model[NX-1][j]]);
		P[NX-1][j] = Mur_X1[2][j] + (vel*dt-dx)/(vel*dt+dx) * (P[NX-2][j] - Mur_X1[3][j]);
	}

	/* Copy Previous Values */
	Mur1stCopy();
}

/* Copy the Filed Values for Mur's 1st Order Absorption */
void Mur1stCopy()
{
	/* Copy Previous Values */
	for(int i=0;i<NX;i++){
		Mur_Y1[i][0] = P[i][0];
		Mur_Y1[i][1] = P[i][1];
		Mur_Y1[i][2] = P[i][NY-2];
		Mur_Y1[i][3] = P[i][NY-1];
	}
	for(int j=0;j<NY;j++){
		Mur_X1[0][j] = P[0][j];
		Mur_X1[1][j] = P[1][j];
		Mur_X1[2][j] = P[NX-2][j];
		Mur_X1[3][j] = P[NX-1][j];
	}
}
