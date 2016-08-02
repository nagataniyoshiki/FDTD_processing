/* 2-D Acoustic FDTD Simulation Demo for Processing 3 (rev. Aug. 2nd, 2016)
    by Yoshiki NAGATANI
      http://ultrasonics.jp/nagatani/e/
      https://twitter.com/nagataniyoshiki

  This is a truly physical simulation program of the sound wave propagation
  in a two-dimensional field filled with fluid media and Mur's 2nd order absorbing
  boundary or a total reflecting wall.
  This program solves the "2D Acoustic FDTD (finite-difference time-domain) method".
  The dark part consists of air (332 m/s) and the brighter part consists of
  a mixture of air and helium (468 m/s). The grid resolution is 10 mm/pixel and
  the time step is 15 us/step (us = micro second = 1/1000000 second). A single pulse
  of sinusoidal wave at 1 kHz with Hann window is transmitted.
  Have fun!!

  ==================================================================================
  [Usage]

   - Click: Create a point source
   - Drag:  Create a line source (i.e. an array of point sources)
   - 'r': Reset the sound field
   - 'c': Toggle Color Map (Jet and Gray scale)
   - 'm': Toggle Displaying Help Messages

  [This program requires:]
   - acoustic_2d_fdtd_core.pde  : Core Solver of 2-D Acoustic FDTD method
   - acoustic_2d_fdtd_mur2nd.pde: Mur's 2nd order absorbing boundary

  ==================================================================================

  *** Technical Notes ***
   This program has been tested with Processing 3.1.2 (64bit) both on Windows and OSX.
   The speed of calculation depends on the CPU power of your PC. Intel Core i3 series
   or higher is recommended.
   For more detailed information about FDTD method (including 3-D elastic simulation),
   please refer our papers on simulation. ==> http://ultrasonics.jp/nagatani/fdtd/
   Thank you.

*/


String photo_file = "lena.png";		// Arbitrary photo file (.png/.jpg/.gif etc.)
int Thresh = 192;			// Threshold for Binarizing Photo Image (0 to 255)

float image_intensity = 400;		// Brightness of Acoustic Field (>0)
float model_intensity = 0.4;		// Brightness of Model (>0)
boolean DrawJetColormap = true;		// true: MATLAB-like jet colormap / false: Gray scale
boolean DisplayHelpMessage = true;

int RandomPointSourceInterval = 2000;	// Interval of Rondom Source Transmission [steps]
int MaxLineSourceNumber = 1000;		// Maximum Number of Points of Line Source

float dx = 10.0e-3;			// Spatial Resolution [m/pixel]
float dt = 15.0e-6;			// Temporal Resolution [s/step]
float dt_over_dx = dt/dx;		// for more efficient calculation


float rho[]   = {1.29,	0.73};		// Densities [kg/m^3]
float kappa[] = {142.0e3, 160.0e3};	// Bulk Moduli [Pa]

float freq = 1.0e3;			// Frequency of Initial Waveform [Hz]

int NX;					// Spatial Size for X direction [pixels]
int NY;					// Spatial Size for Y direction [pixels]

float[][] Vx;				// Particle Velocity for X direction [m/s]
float[][] Vy;				// Particle Velocity for Y direction [m/s]
float[][] P;				// Sound Pressure [Pa]
int[][] Model;				// Model

float[][] Mur_X1;			// Mur's 2nd-Order Absorption Layer
float[][] Mur_X2;			// Mur's 2nd-Order Absorption Layer
float[][] Mur_Y1;			// Mur's 2nd-Order Absorption Layer
float[][] Mur_Y2;			// Mur's 2nd-Order Absorption Layer

int PointSourceX;			// X coordinate of Point Sound Source [pixel]
int PointSourceY;			// Y coordinate of Point Sound Source [pixel]
int[] LineSourceX;			// X coordinate of Line Sound Source [pixel]
int[] LineSourceY;			// Y coordinate of Line Sound Source [pixel]

int n_point = 0;
int n_line = 1000;
float sig_point;
float sig_line;
int LineSourceNumber = 0;
boolean LineSourceDragging = false;
boolean LineSourceRunning = false;
float col,col_r,col_g,col_b;
PImage photo;


/* Preparation *********************************************************/

/* setting() function is for Processing 3.x */
void settings() {
  /* Read Photo */
  photo = loadImage(photo_file);
  NX = photo.height;
  NY = photo.width;
  size(NY,NX);
}

void setup()
{
	/* Read Photo (for Processing 2.x) */
	//photo = loadImage(photo_file);
	//NX = photo.height;
	//NY = photo.width;
	//size(NY,NX);

	PointSourceX = int(NX*0.3);
	PointSourceY = int(NY*0.6);

	/* Allay Allocation */
	Vx = new float[NX+1][NY];
	Vy = new float[NX][NY+1];
	P  = new float[NX][NY];
	Model= new int[NX][NY];
	LineSourceX  = new int[MaxLineSourceNumber];
	LineSourceY  = new int[MaxLineSourceNumber];

	Mur_X1 = new float[4][NY];
	Mur_X2 = new float[4][NY];
	Mur_Y1 = new float[NX][4];
	Mur_Y2 = new float[NX][4];

	/* Initialize Field Values */
	InitField();
	InitMur2nd();

	/* Create Model */
	for(int i=0;i<NX;i++){
		for(int j=0;j<NY;j++){
			color c = photo.pixels[i*NY+j];
			if ( brightness(c) > Thresh )
				Model[i][j] = 1;	// air and helium
			else
				Model[i][j] = 0;	// air
		}
	}

}


/* Main Loop *********************************************************/
void draw()
{
	PImage img = createImage( NY, NX, RGB );

	/* Update the Acoustic Field (This is the main part of FDTD !!)  */
	UpdateV();
	UpdateP();

	/* Mur's 2nd Order Absorption */
	Mur2nd();

	/* Initial Waveform from a Point Source (1 pulse of sinusoidal wave with Hann window) */
	if( n_point < (1.0/freq)/dt ){
		sig_point = (1.0-cos((2.0*PI*freq*n_point*dt)))/2.0 * sin((2.0*PI*freq*n_point*dt));
		P[PointSourceX][PointSourceY] += sig_point;
	}

	/* Initial Waveform from a Line Source (1 pulse of sinusoidal wave with Hann window) */
	if( LineSourceRunning == true){
		if( n_line < (1.0/freq)/dt ){
			sig_line = (1.0-cos((2.0*PI*freq*n_line*dt)))/2.0 * sin((2.0*PI*freq*n_line*dt)) / sqrt(LineSourceNumber);
			for(int i=0; i<LineSourceNumber; i++){
				P[LineSourceX[i]][LineSourceY[i]] += sig_line;
			}
		}
		else{
			LineSourceRunning = false;
		}
	}

	/* Output Step Numbers to Console */
	println(n_point, n_line);

	/* Draw the Acoustic Field */
	for(int i=0; i<NX; i++){
		for(int j=0; j<NY; j++){
			col = abs((float)(P[i][j]*image_intensity + Model[i][j]*model_intensity));	// Value in Gray scale
			if( DrawJetColormap ){	// conversion into Jet colormap
				col = min(4.0,col);	// avoid black color
				col_r = 255*min(max(min(col-1.5,-col+4.5),0),1);
				col_g = 255*min(max(min(col-0.5,-col+3.5),0),1);
				col_b = 255*min(max(min(col+0.5,-col+2.5),0),1);
			}
			else{	// conversion into 255 scales
				col_r = 255*min(max(col,0),1);
				col_g = col_r;
				col_b = col_r;
			}
			img.pixels[i*NY+j] = color(col_r,col_g,col_b);
		}
	}
	/* Draw the Line Source while Mouse Dragging */
	if( LineSourceDragging == true ){
		for(int i=0; i<LineSourceNumber; i++){
			img.pixels[LineSourceX[i]*NY+LineSourceY[i]] = color(255,127,127);
		}
	}

	/* Output the Image */
	image(img, 0, 0);

	/* Text Messages */
	fill(192);
	if( DisplayHelpMessage )
		text("@nagataniyoshiki  [ Click or Drag anywhere / 'r' to reset / 'c' to toggle color map ]", 10, NX-5);
	fill(224);

	n_point++;
	n_line++;

	/* Random Source at given time intervals */
	if( n_point > RandomPointSourceInterval && n_line > RandomPointSourceInterval )
		RandomPointSource();

}


/* Sub Routines ************************************************************/


/* Create a Random Point Source */
void RandomPointSource()
{
	PointSourceX = int(random(NX-1));
	PointSourceY = int(random(NY-1));
	n_point = 0;
}

/* Create a Point Source at Mouse Clicking Point */
void mouseClicked()
{
	if( n_point >= (1.0/freq)/dt ){	// Only one point source can exist simultaneously
		if( 0 < mouseY && mouseY < NX && 0 < mouseX && mouseX < NY ){
			PointSourceX = mouseY;
			PointSourceY = mouseX;
			n_point = 0;
		}
	}
}

/* Create a Line Source on Mouse Dragging Trace */
void mouseDragged()
{
	if( n_line >= (1.0/freq)/dt ){	// Only one line source can exist simultaneously
		// Start Dragging
		if( LineSourceDragging == false ){
			LineSourceDragging = true;
			LineSourceNumber = 0;
		}
		// Add a Point to the Line Source
		if( LineSourceNumber < MaxLineSourceNumber) {
			if( 0 < mouseY && mouseY < NX && 0 < mouseX && mouseX < NY ){
				LineSourceX[LineSourceNumber] = mouseY;
				LineSourceY[LineSourceNumber] = mouseX;
				LineSourceNumber++;
			}
		}
	}
}

/* End Dragging */
void mouseReleased()
{
	if( LineSourceDragging == true && LineSourceRunning == false && n_line >= (1.0/freq)/dt ){
		LineSourceDragging = false;
		LineSourceRunning = true;
		n_line = 0;
	}
}

void keyPressed()
{
	/* Reset Sound Field */
	if (key == 'r') {
		InitField();
		InitMur2nd();
		n_line = int((1.0/freq)/dt + 1);	// Disable Point Sound Source
		n_point = int((1.0/freq)/dt + 1);	// Disable Line Sound Source
	}

	/* Toggle Color map (Jet and Gray scale) */
	if (key == 'c' ) {
		DrawJetColormap = ! DrawJetColormap;
	}

	/* Toggle Displaying Help Message */
	if (key == 'm' ) {
		DisplayHelpMessage = ! DisplayHelpMessage;
	}


}