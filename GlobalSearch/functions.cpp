#include "functions.h"

//double Perebor(int maxTrial, int demension, IProblem* problem)
//{
//  int size;
//  int rank;
//  int delta;
//  int ost;
//  MPI_Comm_size(MPI_COMM_WORLD, &size);
//  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//
//  double min, sum;
//
//  double* lower_bounds = new double[demension];
//  double* upper_bounds = new double[demension];
//  problem->GetBounds(lower_bounds, upper_bounds);
//  double shag = lower_bounds[0];
//  double* masX = new double[maxTrial + 1];
//  masX[0] = shag;
//
//  //Подсчитаем координаты X у точек, в которых потом будем считать значение функции
//  for (int i = 1; i < maxTrial + 1; i++)
//  {
//    if (shag >= upper_bounds[0] + 0.1)
//      break;
//    shag += (upper_bounds[0] - lower_bounds[0]) / maxTrial;
//    masX[i] = shag;
//  }
//
//  int flag = 0;
//  int tmpsize = 0;
//  ost = (maxTrial + 1) % size;
//  if ((maxTrial + 1) <= size)
//  {
//    tmpsize = size;
//    size = 1;
//    ost = 0;
//    flag = -1;
//  }
//
//  delta = (maxTrial + 1) / size;
//
//  if (ost == 0)
//  {
//    if (rank == 0)
//    {
//      for (int i = delta; i < (maxTrial + 1) - delta; i += delta)
//      {
//        MPI_Send(&masX[i], delta, MPI_DOUBLE, i / delta, 0, MPI_COMM_WORLD);
//      }
//    }
//  }
//  else
//  {
//    if (rank == 0)
//    {
//      for (int i = delta + ost; i <= (maxTrial + 1) - delta; i += delta)
//      {
//        MPI_Send(&masX[i], delta, MPI_DOUBLE, (i - ost) / delta, 0, MPI_COMM_WORLD);
//      }
//    }
//  }
//
//  MPI_Status status;
//  double* b = new double[delta];
//  for (int i = 0; i < delta; i++)
//    b[i] = 0;
//  double prom_min;
//  double* tmp = new double[demension];
//  for (int i = 0; i < demension; i++)
//    tmp[i] = 0;
//  if (rank == 0)
//  {
//    for (int i = 0; i < delta + ost; i++)
//    {
//      for (int j = 0; j < demension; j++)
//        tmp[j] = masX[i];
//      sum = problem->CalculateFunctionals(tmp, 0);
//      //printf("Now Func %lf\n", sum);
//
//      if (i == 0)
//        prom_min = sum;
//      else
//        if (sum < prom_min)
//          prom_min = sum;
//    }
//  }
//  else
//  {
//    if (flag == 0)
//    {
//      MPI_Recv(&b[0], delta, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
//      for (int i = 0; i < delta; i++)
//      {
//        for (int j = 0; j < demension; j++)
//          tmp[j] = b[i];
//        sum = problem->CalculateFunctionals(tmp, 0);
//        //printf("Now Func %lf\n", sum);
//
//        if (i == 0)
//          prom_min = sum;
//        else
//          if (sum < prom_min)
//            prom_min = sum;
//      }
//    }
//  }
//
//  if (flag == 0)
//  {
//    MPI_Reduce(&prom_min, &min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
//  }
//  else
//  {
//    if (rank == 0)
//    {
//      return prom_min;
//    }
//    else
//    {
//      return -1;
//    }
//  }
//  
//  MPI_Barrier(MPI_COMM_WORLD);
//  delete[] masX;
//  delete[] tmp;
//  delete[] lower_bounds;
//  delete[] upper_bounds;
//  delete[] b;
//  return min;
//}

void QuickSort(double* mas, int first, int last)
{
  double count;
  int f = first, l = last;
  double mid = mas[(f + l) / 2]; //вычисление опорного элемента

  do
  {
    while (mas[f] < mid)
      f++;
    while (mas[l] > mid)
      l--;

    if (f <= l) //перестановка элементов
    {
      count = mas[f];
      mas[f] = mas[l];
      mas[l] = count;

      f++;
      l--;
    }
  } while (f < l);
  if (first < l)
    QuickSort(mas, first, l);
  if (f < last)
    QuickSort(mas, f, last);
}

//Вопрос: как работать с точностью? Мы же не знаем, сколько итераций потребуется, соответственно не знаем
//массив какого размера создавать
double AGP(int maxTrial, double accuracy, int demension, IProblem* problem)
{
  double* lower_bounds = new double[demension];
  double* upper_bounds = new double[demension];
  problem->GetBounds(lower_bounds, upper_bounds);

  double* X = new double[maxTrial];
  for (int i = 0; i < maxTrial; i++)
    X[i] = 0;
  X[0] = lower_bounds[0];
  X[1] = upper_bounds[0];

  //Для подсчета минимума
  double* tmpmin = new double[demension];
  tmpmin[0] = X[0];
  double min = problem->CalculateFunctionals(tmpmin, 0);
  double prom_min;
  int iteration = 0;

  int first = 0;
  int last = 1;
  int k = 2;

  double* tmp = new double[demension];
  double z0;
  double z1;
  double M;

  int t = 1;
  double r = 4;
  double delta;
  double second;
  double third;
  double R;

  for (int i = 0; i < maxTrial; i++)
  {
    //упорядочить точки по координате
    QuickSort(X, first, last);

    t = 1;

    //вычислить оценку M для неизвестной константы Липшеца
    for (int j = 1; j < k; j++)
    {
      double m;
      tmp[0] = X[j - 1];
      z0 = problem->CalculateFunctionals(tmp, 0);
      tmp[0] = X[j];
      z1 = problem->CalculateFunctionals(tmp, 0);
      m = abs(z1 - z0) / (X[j] - X[j - 1]);
      if (j == 1)
        M = m;
      if (m > M)
        M = m;
    }

    //Вычислим характеристику R(i);
    for (int j = 1; j < k; j++)
    {
      double Rtmp;
      tmp[0] = X[j - 1];
      z0 = problem->CalculateFunctionals(tmp, 0);
      tmp[0] = X[j];
      z1 = problem->CalculateFunctionals(tmp, 0);
      delta = X[j] - X[j - 1];
      second = (z1 - z0)*(z1 - z0) / (r*r*M*M*delta);
      third = 2 * (z1 + z0) / (r*M);
      Rtmp = delta + second - third;
      if (j == 1)
        R = Rtmp;
      if (Rtmp > R)
      {
        R = Rtmp;
        t = j;
      }
    }

    k++;
    tmp[0] = X[t - 1];
    z0 = problem->CalculateFunctionals(tmp, 0);
    tmp[0] = X[t];
    z1 = problem->CalculateFunctionals(tmp, 0);
    second = (X[t] + X[t - 1]) / 2;
    third = (z1 - z0) / (2 * r*M);
    X[k - 1] = second - third;

    tmpmin[0] = X[k - 1];
    prom_min = problem->CalculateFunctionals(tmpmin, 0);
    if (prom_min < min)
      min = prom_min;
    iteration++;

    if (abs(X[k - 1] - X[t]) < accuracy)
      break;
    else if (abs(X[k - 1] - X[t - 1]) < accuracy)
      break;
    last++;
  }

  printf("Number of iteration = %d\n", iteration);
  return min;
}