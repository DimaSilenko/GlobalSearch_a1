#include "evolvent.h"
#include "problem_manager.h"
#include "extended.h"
#include "mpi.h"
#include "omp.h"
#include <ctime>

//-----------------------------------------------
#include <iostream>
#include "../SRC/sample_src/HansenProblem.hpp"
#include "../SRC/sample_src/HansenProblemFamily.hpp"

#include "../SRC/sample_src/HillProblem.hpp"
#include "../SRC/sample_src/HillProblemFamily.hpp"

#include "../SRC/sample_src/ShekelProblem.hpp"
#include "../SRC/sample_src/ShekelProblemFamily.hpp"

#include "../SRC/sample_src/OptSqConstrProblem.hpp"

#include "../SRC/sample_src/grishagin_function.hpp"
#include "../SRC/sample_src/GrishaginProblemFamily.hpp"
#include "../SRC/sample_src/GrishaginConstrainedProblem.hpp"
#include "../SRC/sample_src/GrishaginConstrainedProblemFamily.hpp"

#include "../SRC/sample_src/GKLSProblem.hpp"
#include "../SRC/sample_src/GKLSProblemFamily.hpp"
#include "../SRC/sample_src/GKLSConstrainedProblem.hpp"
#include "../SRC/sample_src/GKLSConstrainedProblemFamily.hpp"
//---------------------------------------------------------------

class Interval
{
protected:
  std::pair <double, double>* zip_x;
  std::pair <double, double>* z;
  double* R;
  int HeapSize;
public:
  Interval();
  Interval(int maxTrials);

  void out();

  //Добавление нового элемента
  void AddElem(double R, double x_l, double x_r, double z_l, double z_r);
  //Обновляем все остальные R
  void Update(double M, double r, int demension);
  //Извлекаем левую границу интервала с макс R
  double GetMaxPointLeft();
  //Извлекаем правую границу интервала с макс R
  double GetMaxPointRight();
  //Извлекаем значение функции, соответствующее левой границе интервала с макс R
  double GetMaxZLeft();
  //Извлекаем значение функции, соответствующее правой границе интервала с макс R
  double GetMaxZRight();
  //Извлекаем верхушку кучи (макс R)
  double GetMax();
  //Получаем размер кучи
  int GetSize();
  //Восстанавливаем свойство кучи
  void Heapify(int i);
};

double Perebor(int maxTrial, int demension, IProblem* problem, double* *bestX);
void QuickSort(double* mas, int first, int last);
void QuickSort_D(double* mas, int first, int last, double* z);
double AGP(int maxTrial, double accuracy, int demension, IProblem* problem);
double AGP_Space(int maxTrial, double accuracy, int demension, IProblem* problem,
  double* *BestX, int *iteration, double r, double* Truemin, int numThr, int flag, char* argv[]);