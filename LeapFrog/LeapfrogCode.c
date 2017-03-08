#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAXPNT 20000

//Function Declarations
void PerformIntergration();
void LeapFrogStep(double x[], double y[], double z[], double vx[], double vy[], double vz[], int n, double dt);
void WriteStateToFile();

double AccelerationX(double ax[], double x[], double y[], double z[], int n, int i);
double AccelerationY(double ax[], double x[], double y[], double z[], int n, int i);
double AccelerationZ(double ax[], double x[], double y[], double z[], int n, int i);



int main() {
FILE *fp;

// Constants
const double  _dt; /* timestep for integration */

// Variables
int _n /* number of points */, _mstep /* number of steps to take (max) */,
		_nout /* steps between outputs */, _nstep /* current step */;

double _x[MAXPNT]; //x-Position
double _y[MAXPNT]; //y-Position
double _z[MAXPNT]; //z-Position
double _r[MAXPNT]; //radial position

double _vx[MAXPNT]; //Velocity in x Component
double _vy[MAXPNT]; //Velocity in y Component
double _vz[MAXPNT]; //Velocity in z Component

double _ax[MAXPNT]; //Acceleration in x Component
double _ay[MAXPNT]; //Acceleration in y Component
double _az[MAXPNT]; //Acceleration in z Component

double _tNow; // time at that point


//Main Body Code

fp = fopen("OutputFile.data in shell","w+");

printf("Hello world");
fprintf(fp, "Hello world in file");

fclose(fp);

return 0;
}


double AccelerationX(double ax[], double x[], double y[], double z[], int n, int i) {
	// IMPORTANT: we are going to define a half step as i+=1 & a full step  by i+=2
		double w = sqrt((_k / (_a*_a*_a)));
		_r[i] = sqrt((x[i]*x[i]) + (y[i]*y[i]) + (z[i]*z[i])); // definition of r
		return -(((_a*_a*_a)*(w*w*w)) / (_r[i]*_r[i]*_r[i]))* (x[i] - _a*_e);
}

double AccelerationY(double ay[], double x[], double y[], double z[], int n, int i) {
	// IMPORTANT: we are going to define a half step as i+=1 & a full step  by i+=2
		double w = sqrt((_k / (_a*_a*_a)));
		_r[i] = sqrt((x[i]*x[i]) + (y[i]*y[i]) + (z[i]*z[i])); // definition of r
		return -(((_a*_a*_a)*(w*w*w)) / (_r[i]*_r[i]*_r[i]))* (y[i]);
		//Y component of Second derivative of r, from eliptical orbit equations
}

double AccelerationZ(double az[], double x[], double y[], double z[], int n, int i) {
	// IMPORTANT: we are going to define a half step as i+=1 & a full step  by i+=2
		double w = sqrt((_k / (_a*_a*_a)));
		_r[i] = sqrt((x[i]*x[i]) + (y[i]*y[i]) + (z[i]*z[i])); // definition of r
		return -(((_a*_a*_a)*(w*w*w)) / (_r[i]*_r[i]*_r[i]))* (y[i]);
		//Y component of Second derivative of r, from eliptical orbit equations
}


void PerformIntergration() {

	for (_nstep = 0; _nstep < _mstep; _nstep++) { // Loop mstep times
		if (_nstep % _nout == 0) {
			//WriteStateToFile();
		}
		LeapFrogStep(_x, _y, _z, _vx, _vy, _vz, _n, _dt); //take intergration step
		_tNow += _dt; //update time
	}

	if (_mstep % _nout == 0) { // if last output wanted
		//WriteStateToFile();
	}
}

/* LeapFrogStep: map t to t+dt. Warning: not accurate unless the
timestep dt is fixed from one call to another */
void LeapFrogStep(double x[], double y[], double z[], double vx[], double vy[], double vz[], int n, double dt) {
	int i;
	double ax[MAXPNT];
	double ay[MAXPNT];
	double az[MAXPNT];
	x[0] = _x[0];
	vx[0] = _vx[0];
	y[0] = _y[0];
	vy[0] = _vy[0];
	z[0] = _z[0];
	vz[0] = _vz[0];

	_outputFile << "_nstep: " << _nstep << "\n";

	if (_nstep == 0) {
		// IMPORTANT: we are going to define a half step as i+=1 & a full step by i+=2
		for (i = 0; i < 2 * n; i += 2) {
			ax[2 * i] = AccelerationX(ax, x, y, z, n, 2 * i); //__1.0__
			ay[2 * i] = AccelerationY(ay, x, y, z, n, 2 * i);
			az[2 * i] = AccelerationZ(az, x, y,z, n, 2 * i);
			}

		vx[1] = vx[0] + ax[0] * dt * 0.5; //Taylor Approximate first half step __2__
	    vy[1] = vy[0] + ay[0] * dt * 0.5;
        vz[1] = vz[0] + az[0] * dt * 0.5;
		}

	for (i = (1 + _nstep); i < (n + _nstep) + 1; i += 2) {

		x[2 * i] = x[(2 * i) - 2] + dt * vx[(2 * i) - 1]; //advance position by full-step __3__
		y[2 * i] = y[(2 * i) - 2] + dt * vy[(2 * i) - 1]; //advance position by full-step
		z[2 * i] = z[(2 * i) - 2] + dt * vz[(2 * i) - 1]; //advance position by full-step

        ax[2 * i] = AccelerationX(ax, x, y, z, n, 2 * i); //__1.1__
        ay[2 * i] = AccelerationY(ay, x, y, z, n, 2 * i);
        az[2 * i] = AccelerationZ(az, x, y,z, n, 2 * i);

		vx[(2 * i) + 1] = vx[(2 * i) - 1] + dt * ax[2 * i]; //advance position by full-step  __4__
		vy[(2 * i) + 1] = vy[(2 * i) - 1] + dt * ay[2 * i]; //advance position by full-step
		vz[(2 * i) + 1] = vz[(2 * i) - 1] + dt * az[2 * i]; //advance position by full-step

	}

}

void WriteStateToFile() {
            fprintf(fp, "%f \t", _tNow);
            fprintf(fp, "%f \t",_x[2 * (_nstep)]);
            fprintf(fp, "%f \t", _vx[(2 * (_nstep)) + 1]);
            fprintf(fp, "%f \t",_y[2 * (_nstep)]);
            fprintf(fp, "%f \t", _vy[(2 * (_nstep)) + 1]);
            fprintf(fp, "%f \t",_z[2 * (_nstep)]);
            fprintf(fp, "%f \t", _vz[(2 * (_nstep)) + 1]);
}
