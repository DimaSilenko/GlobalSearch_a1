#ifndef __ISOLINES_PLOTTER_H__
#define __ISOLINES_PLOTTER_H__

#include <string>
#include "../GlobalSearch/problem_interface.h"

void plot2dProblemIsolines(IProblem* problem, const std::string& outputPath,
  const std::string& pointsPath, bool isUPPlot = false, double* point = 0);

#endif
