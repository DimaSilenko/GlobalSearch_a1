#include "isolinesPlotter.h"

#include <dislin/include/discpp.h>
#include <algorithm>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>

#define GRID_SIZE 300

typedef std::pair<double, double> point2d;

size_t readPoints(const std::string& filename, std::vector<point2d>& points)
{
  std::ifstream input;
  std::string currentLine(512, ' ');

  input.open(filename, std::ios_base::in);

  if (!input.is_open())
    return 0;

  input.getline(&currentLine[0], currentLine.size());
  size_t numberOfPoints = std::stoi(currentLine, NULL);
  points.reserve(numberOfPoints + 2);

  while (!input.eof()) {
    size_t nextPosition;
    point2d currentPoint;
    input.getline(&currentLine[0], currentLine.size());
    currentPoint.first = std::stod(currentLine, &nextPosition);
    currentPoint.second = std::stod(currentLine.substr(nextPosition));
    points.push_back(currentPoint);
  }

  input.close();

  return numberOfPoints;
}

void plot2dProblemIsolines(IProblem* problem, const std::string& outputPath,
  const std::string& pointsPath, bool isUPPlot, double* point)
{
  /// Количество линий уровней
  const int linesNumber = 50;
  /// Максимальная размерность задачи
  const int maxDim = 200;
  // Расскоментировать если в плоттере необходимо строить более полотные линии уровней
  //isUPPlot = true;

  std::vector<point2d> points;
  size_t totalPoints = readPoints(pointsPath, points);

  int n = GRID_SIZE, i, j, width = 900, height = 900;
  double  x_step, y_step, x_left, y_left, x_right, y_right;
  double arg[maxDim];


  double *xray, *yray, **zmat;

  xray = new double[GRID_SIZE];
  yray = new double[GRID_SIZE];
  double* zmatData = new double[GRID_SIZE * GRID_SIZE];
  zmat = new double*[GRID_SIZE];

  for (int i = 0; i < GRID_SIZE; i++)
    zmat[i] = &zmatData[i * GRID_SIZE];

  double lb[maxDim], ub[maxDim];
  problem->GetBounds(lb, ub);

  double nArg[maxDim];
  int dim = problem->GetDimension();

  if (point == 0)
    problem->GetOptimumPoint(nArg);
    //for (int q = 0; q < dim; q++)
    //  nArg[q] = (ub[q] + lb[q]) / 2.0;
  else
    for (int q = 0; q < dim; q++)
      nArg[q] = point[q];

  //problem->GetOptimumPoint(nArg);

  int startCoordinate = dim - 2;

  x_left = lb[startCoordinate];
  y_left = lb[startCoordinate + 1];
  x_right = ub[startCoordinate];
  y_right = ub[startCoordinate + 1];

  x_step = (x_right - x_left) / (n - 1);
  y_step = (y_right - y_left) / (n - 1);

  for (i = 0; i < n; i++)
  {
    xray[i] = x_left + i * x_step;
    yray[i] = y_left + i * y_step;
  }

  Dislin g;
  g.metafl("png");
  g.winsiz(width, height);
  g.pagmod("LAND");
  g.page(2400, 2400);
  g.setfil(outputPath.c_str());
  g.sclfac(2.0);
  g.filmod("VERSION");
  g.scrmod("revers");
  g.disini();
  g.complx();
  g.name("X-axis", "x");
  g.name("Y-axis", "y");

  g.axspos(240, 2200);
  g.axslen(2100, 2100);

  g.graf(x_left, x_right, x_left, (x_right - x_left) / 4,
    y_left, y_right, y_left, (y_right - y_left) / 4);

  g.height(30);

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
    {
      arg[0] = xray[i]; arg[1] = yray[j];
      nArg[startCoordinate] = arg[0];
      nArg[startCoordinate + 1] = arg[1];

      zmat[i][j] = -HUGE_VAL;
      for (int k = 0; k < problem->GetNumberOfConstraints(); k++)
      {
        zmat[i][j] = problem->CalculateFunctionals(nArg, k);
        if (zmat[i][j] > 0.0)
          break;
      }
    }
  int colors[] = { 150 };
  g.shdmod("LOWER", "CELL");
  g.shdmod("UPPER", "COLOR");
  g.conclr(colors, 1);
  if (problem->GetNumberOfFunctions() != 1)
  {
    double levels[] = { 0.0 };
    g.conshd(xray, n, yray, n, zmatData, levels, 1);
  }

  int targetIndex = problem->GetNumberOfFunctions() - 1;
  double minZ = HUGE_VAL;
  double maxZ = -HUGE_VAL;
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
    {
      arg[0] = xray[i]; arg[1] = yray[j];
      nArg[startCoordinate] = arg[0];
      nArg[startCoordinate + 1] = arg[1];
      zmat[i][j] = problem->CalculateFunctionals(nArg, targetIndex);
      minZ = std::min(zmat[i][j], minZ);
      maxZ = std::max(zmat[i][j], maxZ);
    }

  double step = (maxZ - minZ) / linesNumber;
  double* zlevs = new double[linesNumber];
  zlevs[0] = minZ;
  zlevs[linesNumber - 1] = maxZ;

  double d = maxZ - minZ;
  for (int t = 1; t < linesNumber - 1; t++)
  {
    double xt = ((double)t) / (linesNumber - 1);
    double xxt = xt * xt * xt;
    zlevs[t] = xxt * d + minZ;
  }

  if (isUPPlot)
  {
    for (i = 0; i < linesNumber; i++)
    {
      double zlev = zlevs[i];
      g.contur(xray, n, yray, n, zmatData, zlev);
    }
  }
  else
  {
    for (i = 0; i < linesNumber; i++)
    {
      double zlev = minZ + i*step + 0.05;
      g.contur(xray, n, yray, n, zmatData, zlev);
    }
  }

  for (int k = 0; k < problem->GetNumberOfConstraints(); k++)
  {
    for (i = 0; i < n; i++)
      for (j = 0; j < n; j++)
      {
        arg[0] = xray[i]; arg[1] = yray[j];
        nArg[startCoordinate] = arg[0];
        nArg[startCoordinate + 1] = arg[1];
        zmat[i][j] = problem->CalculateFunctionals(nArg, k);
      }
    g.contur(xray, n, yray, n, zmatData, 0.0);
  }

  if (points.size())
  {
    g.color("white");
    g.hsymbl(15);

    for (size_t k = 0; k < totalPoints; k++)
      g.rlsymb(21, points[k].first, points[k].second);
    g.color("red");
    g.hsymbl(25);
    g.rlsymb(21, points[totalPoints].first, points[totalPoints].second);
    g.hsymbl(20);
    if(totalPoints == points.size() - 2)
    {
      g.color("blue");
      g.rlsymb(21, points[totalPoints + 1].first, points[totalPoints + 1].second);
    }
  }

  g.height(50);
  g.color("fore");
  g.title();
  g.disfin();
}
