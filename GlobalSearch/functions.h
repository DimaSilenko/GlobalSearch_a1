#include "evolvent.h"
#include "problem_manager.h"
#include "extended.h"
#include "mpi.h"

//double Perebor(int maxTrial, int demension, IProblem* problem);
void QuickSort(double* mas, int first, int last);
double AGP(int maxTrial, double accuracy, int demension, IProblem* problem);