#ifndef RESULT_SET_H
#define RESULT_SET_H
struct result_set {
	QString casename;
	int nsteps;
 	double *tnow;      // simulation time (mins)
    double *pData[maxGraphs];
    double maxValue[maxGraphs];
};
typedef result_set RESULT_SET;
#endif
