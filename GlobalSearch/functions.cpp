#include "functions.h"

using namespace std;

Interval::Interval()
{
  int s = 10000;
  zip_x = new pair<double, double>[s-1];
  z = new pair<double, double>[s - 1];
  R = new double[s - 1];
  HeapSize = 0;
}

Interval::Interval(int maxTrials)
{
  zip_x = new pair<double, double>[maxTrials - 1];
  z = new pair<double, double>[maxTrials - 1];
  R = new double[maxTrials - 1];
  HeapSize = 0;
}

void Interval::out()
{
  int i = 0;
  int k = 1;
  while (i < HeapSize)
  {
    while ((i < k) && (i < HeapSize))
    {
      cout << R[i] << " ";
      i++;
    }
    cout << endl;
    k = k * 2 + 1;
  }
}

void Interval::AddElem(double _R, double x_l, double x_r, double z_l, double z_r)
{
  int i = HeapSize;
  R[i] = _R;
  zip_x[i].first = x_l;
  zip_x[i].second = x_r;
  z[i].first = z_l;
  z[i].second = z_r;
  int parent = (i - 1) / 2;

  while (parent >= 0 && i > 0)
  {
    if (R[i] > R[parent])
    {
      double temp = R[i];
      R[i] = R[parent];
      R[parent] = temp;

      temp = zip_x[i].first;
      zip_x[i].first = zip_x[parent].first;
      zip_x[parent].first = temp;

      temp = zip_x[i].second;
      zip_x[i].second = zip_x[parent].second;
      zip_x[parent].second = temp;

      temp = z[i].first;
      z[i].first = z[parent].first;
      z[parent].first = temp;

      temp = z[i].second;
      z[i].second = z[parent].second;
      z[parent].second = temp;
    }
    i = parent;
    parent = (i - 1) / 2;
  }
  HeapSize++;
}

void Interval::Update(double M, double r, int demension)
{
  double N = 1.0 / demension;
  double delta = 0;
  double second = 0;
  double third = 0;

  for (int i = 0; i < HeapSize; i++)
  {
    if (zip_x[i].first == 0)
    {
      delta = pow((zip_x[i].second - zip_x[i].first), N);
      R[i] = 2 * delta - 4 * z[i].second / (r*M);
    }
    else if (zip_x[i].second == 1)
    {
      delta = pow((zip_x[i].second - zip_x[i].first), N);
      R[i] = 2 * delta - 4 * z[i].first / (r*M);
    }
    else
    {
      delta = pow((zip_x[i].second - zip_x[i].first), N);
      second = (z[i].second - z[i].first)*(z[i].second - z[i].first) / (r*r*M*M*delta);
      third = 2 * (z[i].second + z[i].first) / (r*M);
      R[i] = delta + second - third;
    }
    
    int parent = (i - 1) / 2;
    int j = i;

    while (parent >= 0 && j > 0)
    {
      if (R[j] > R[parent])
      {
        double temp = R[j];
        R[j] = R[parent];
        R[parent] = temp;

        temp = zip_x[j].first;
        zip_x[j].first = zip_x[parent].first;
        zip_x[parent].first = temp;

        temp = zip_x[j].second;
        zip_x[j].second = zip_x[parent].second;
        zip_x[parent].second = temp;

        temp = z[j].first;
        z[j].first = z[parent].first;
        z[parent].first = temp;

        temp = z[j].second;
        z[j].second = z[parent].second;
        z[parent].second = temp;
      }
      j = parent;
      parent = (j - 1) / 2;
    }
  }
}

double Interval::GetMaxPointLeft()
{
  return zip_x[0].first;
}

double Interval::GetMaxPointRight()
{
  return zip_x[0].second;
}

double Interval::GetMaxZLeft()
{
  return z[0].first;
}

double Interval::GetMaxZRight()
{
  return z[0].second;
}

double Interval::GetMax()
{
  double x = R[0];
  R[0] = R[HeapSize - 1];
  zip_x[0].first = zip_x[HeapSize - 1].first;
  zip_x[0].second = zip_x[HeapSize - 1].second;
  z[0].first = z[HeapSize - 1].first;
  z[0].second = z[HeapSize - 1].second;
  HeapSize--;
  Heapify(0);
  return(x);
}

int Interval::GetSize()
{
  return HeapSize;
}

void Interval::Heapify(int i)
{
  int left = 2 * i + 1;
  int right = 2 * i + 2;
  double temp;
  if (left < HeapSize)
  {
    if (R[i] < R[left])
    {
      temp = R[i];
      R[i] = R[left];
      R[left] = temp;

      temp = zip_x[i].first;
      zip_x[i].first = zip_x[left].first;
      zip_x[left].first = temp;

      temp = zip_x[i].second;
      zip_x[i].second = zip_x[left].second;
      zip_x[left].second = temp;

      temp = z[i].first;
      z[i].first = z[left].first;
      z[left].first = temp;

      temp = z[i].second;
      z[i].second = z[left].second;
      z[left].second = temp;

      Heapify(left);
    }
  }

  if (right < HeapSize)
  {
    if (R[i] < R[right])
    {
      temp = R[i];
      R[i] = R[right];
      R[right] = temp;

      temp = zip_x[i].first;
      zip_x[i].first = zip_x[right].first;
      zip_x[right].first = temp;

      temp = zip_x[i].second;
      zip_x[i].second = zip_x[right].second;
      zip_x[right].second = temp;

      temp = z[i].first;
      z[i].first = z[right].first;
      z[right].first = temp;

      temp = z[i].second;
      z[i].second = z[right].second;
      z[right].second = temp;

      Heapify(right);
    }
  }
}
//----------------------Конец описания класса--------------------------------------------------

double Perebor(int maxTrial, int demension, IProblem* problem, double* *bestX)
{
  int size;
  int rank;
  int delta;
  int ost;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  double min, sum;

  double* lower_bounds = new double[demension];
  double* upper_bounds = new double[demension];
  problem->GetBounds(lower_bounds, upper_bounds);
  double* shag = new double[demension];
  for (int i =0; i < demension; i++)
    shag[i] = lower_bounds[i];
  double** masX = new double*[maxTrial + 1];
  for (int i = 0; i < maxTrial + 1; i++)
    masX[i] = new double[demension];
  for (int i = 0; i < demension; i++)
    masX[0][i] = shag[i];

  //Подсчитаем координаты X у точек, в которых потом будем считать значение функции
  for (int i = 1; i < maxTrial + 1; i++)
  {
    for (int j = 0; j < demension; j++) {
      if (shag[j] >= upper_bounds[j] + 0.1)
        break;
      shag[j] += (upper_bounds[j] - lower_bounds[j]) / (maxTrial + 1);
      masX[i][j] = shag[j];
    }
  }

  int flag = 0;
  int tmpsize = 0;
  ost = (maxTrial + 1) % size;
  if ((maxTrial + 1) <= size)
  {
    tmpsize = size;
    size = 1;
    ost = 0;
    flag = -1;
  }

  delta = (maxTrial + 1) / size;

  double prom_min = 0;
  if (rank == 0)
  {
    double* tmp = new double[demension];

    for (int i = 0; i < delta + ost; i++)
    {
      if (demension == 2) {
        tmp[0] = masX[i][0];
        for (int j = 0; j < delta + ost; j++) {
          tmp[1] = masX[j][1];

          sum = problem->CalculateFunctionals(tmp, 0);

          if ((j == 0)&&(i==0)) {
            prom_min = sum;
            (*bestX)[0] = tmp[0];
            (*bestX)[1] = tmp[1];
          }
          else
            if (sum < prom_min) {
              prom_min = sum;
              (*bestX)[0] = tmp[0];
              (*bestX)[1] = tmp[1];
            }
        }
      }
      if (demension == 1) {
        for (int j = 0; j < demension; j++)
          tmp[j] = masX[i][j];
        sum = problem->CalculateFunctionals(tmp, 0);

        if (i == 0) {
          prom_min = sum;
          (*bestX)[0] = tmp[0];
        }
        else
          if (sum < prom_min) {
            prom_min = sum;
            (*bestX)[0] = tmp[0];
          }
      }
    }
    min = prom_min;
    //delete[] tmp;
  }
  //else
  //{
  //  if (flag == 0)
  //  {
  //    double* tmp = new double[demension];
  //    for (int i = delta*rank + ost; i < delta* rank + ost + delta; i++)
  //    {
  //      if (demension == 2) {
  //        tmp[0] = masX[i][0];
  //        for (int j = delta * rank + ost; j < delta* rank + ost + delta; j++) {
  //          tmp[1] = masX[j][1];

  //          sum = problem->CalculateFunctionals(tmp, 0);

  //          if (i == delta * rank + ost) {
  //            prom_min = sum;
  //          }
  //          else
  //            if (sum < prom_min)
  //              prom_min = sum;
  //        }
  //      }
  //      if (demension == 1) {
  //        for (int j = 0; j < demension; j++)
  //          tmp[j] = masX[i][j];
  //        sum = problem->CalculateFunctionals(tmp, 0);

  //        if (i == 0) {
  //          prom_min = sum;
  //        }
  //        else
  //          if (sum < prom_min)
  //            prom_min = sum;
  //      }
  //    }
  //    //delete[] tmp;
  //  }
  //}

  if (flag == 0)
  {
    MPI_Reduce(&prom_min, &min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  }
  else
  {
    if (rank == 0)
    {
      return prom_min;
    }
    else
    {
      return -1;
    }
  }
  
  delete[] masX;
  delete[] lower_bounds;
  delete[] upper_bounds;
  return min;
}

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

void QuickSort_D(double* mas, int first, int last, double* z)
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

      count = z[f];
      z[f] = z[l];
      z[l] = count;

      f++;
      l--;
    }
  } while (f < l);
  if (first < l)
    QuickSort_D(mas, first, l, z);
  if (f < last)
    QuickSort_D(mas, f, last, z);
}

int sng(double x) {
  if (x > 0)
    return 1;
  else
    return -1;
}

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
  //delete[] tmp;
  //delete[] tmpmin;
  delete[] lower_bounds;
  delete[] upper_bounds;

  return min;
}

// Для многомерного случая------------------------------------------------------------
double AGP_Space(int maxTrial, double accuracy, int demension, IProblem* problem, double* *BestX, int *iteration, double r, double* Truemin, int NumThr, int flag, char* argv[])
{
  double N = 1.0 / demension;

  double* lower_bounds = new double[demension];
  double* upper_bounds = new double[demension];
  problem->GetBounds(lower_bounds, upper_bounds);

  //Многомерная точка
  double* X = new double[demension];
  for (int j = 0; j < demension; j++)
    X[j] = 0;

  double z_l, z_r, z;

  Interval inter(maxTrial);

  //Нам надо вычислить значение на левой гарнице и сжать в отрезок [0; 1] и то же самое с правой границей
  double zip_x_l = 0;
  double zip_x_r = 1;
  double zip_x = 0;

  //Для подсчета минимума
  double min;
  GetImage(zip_x_l, X, lower_bounds, upper_bounds, demension, 10);
  min = problem->CalculateFunctionals(X, 0);
  z_l = min;

  GetImage(zip_x_r, X, lower_bounds, upper_bounds, demension, 10);
  z_r = problem->CalculateFunctionals(X, 0);

  for (int i = 0; i < demension; i++)
    (*BestX)[i] = X[i];
  double prom_min;

  double M = 1;
  
  // Флаг на обновление M
  int f = 0;

  double delta;
  double second;
  double third;
  double R = DBL_MIN;
  double p;

  //Открываем файл для записи-------------------------------------------------------------------
 /* std::string pointLogName = argv[3];
  FILE* pointLog = fopen(pointLogName.c_str(), "w");*/
  //--------------------------------------------------------------------------------------------

  for (int i = 0; i < maxTrial; i++)
  {
    inter.out();
    cout << "\n";
    if (inter.GetSize() != 0)
    {
      zip_x_l = inter.GetMaxPointLeft();
      zip_x_r = inter.GetMaxPointRight();
      z_l = inter.GetMaxZLeft();
      z_r = inter.GetMaxZRight();
    }

    if (zip_x_l != 0 && zip_x_r != 1) {
      second = (zip_x_r + zip_x_l) / 2;
      p = pow((abs(z_r - z_l) / M), demension);
      third = p / 2 / r;
      zip_x = second - sng(z_r - z_l) * third;
    }
    else
      zip_x = (zip_x_r + zip_x_l) / 2;


    GetImage(zip_x, X, lower_bounds, upper_bounds, demension, 10);
    prom_min = problem->CalculateFunctionals(X, 0);
    z = prom_min;
    if (prom_min < min) {
      min = prom_min;
      for (int kk = 0; kk < demension; kk++)
        (*BestX)[kk] = X[kk];
    }
    (*iteration)++;

    if ((zip_x >= zip_x_r) || (zip_x <= zip_x_l))
      std::cout << "ERROR!!!!!!!\n";

    if (flag == 0) {
      p = pow(abs(zip_x - zip_x_r), N);
      if (p < accuracy)
        break;
      p = pow(abs(zip_x - zip_x_l), N);
      if (p < accuracy)
        break;
    }
    else {
      double maxabs = abs(Truemin[0] - (*BestX)[0]);
      if (abs(Truemin[1] - (*BestX)[1]) > maxabs)
        maxabs = abs(Truemin[1] - (*BestX)[1]);

      if ((abs(Truemin[0] - (*BestX)[0]) < accuracy) && (abs(Truemin[1] - (*BestX)[1]) < accuracy))
        break;
    }

    //вычислить оценку M для неизвестной константы Липшеца
    double m = 1;
    p = pow(abs(zip_x - zip_x_r), N);
    m = abs(z - z_r) / p;
    if (m > M)
    {
      M = m;
      f = 1;
    }
    p = pow(abs(zip_x - zip_x_l), N);
    m = abs(z - z_l) / p;
    if (m > M)
    {
      M = m;
      f = 1;
    }

    if (i == 0) {
      M = 1;
      f = 1;
    }

    //Вычислим характеристику R(i)
    R = DBL_MIN;

    if (inter.GetSize() != 0)
      inter.GetMax();

    if (zip_x_l == 0)
    {
      delta = pow((zip_x - zip_x_l), N);
      R = 2 * delta - 4 * z / (r*M);
    }
    else
    {
      delta = pow((zip_x - zip_x_l), N);
      second = (z - z_l)*(z - z_l) / (r*r*M*M*delta);
      third = 2 * (z + z_l) / (r*M);
      R = delta + second - third;
    }

    inter.AddElem(R, zip_x_l, zip_x, z_l, z);

    //-----------------------------------------------
    R = DBL_MIN;

    if (zip_x_r == 1)
    {
      delta = pow((zip_x_r - zip_x), N);
      R = 2 * delta - 4 * z / (r*M);
    }
    else
    {
      delta = pow((zip_x_r - zip_x), N);
      second = (z_r - z)*(z_r - z) / (r*r*M*M*delta);
      third = 2 * (z_r + z) / (r*M);
      R = delta + second - third;
    }

    inter.AddElem(R, zip_x, zip_x_r, z, z_r);

    if (f == 1)
    {
      inter.Update(M, r, demension);
      f = 0;
    }

  }

  ////------------------------------------------------------------------------------------------
  //fprintf(pointLog, "%d\n", (*iteration));
  ////Записываем в файл точки----------------------------------------------------------------------------------
  //for (int i = 1; i < (*iteration); i++)
  //{
  //  for (int j = 0; j < demension; j++)
  //    fprintf(pointLog, "%lf ", X[i][j]);
  //  fprintf(pointLog, "\n");
  //}
  ////---------------------------------------------------------------------------------------------------------------

  ///// Печатаем найденную точку
  //for (int j = 0; j < demension; j++)
  //  fprintf(pointLog, "%lf ", (*BestX)[j]);
  //fprintf(pointLog, "\n");

  ///// Известная точка глобального минимума
  //for (int j = 0; j < demension; j++)
  //  fprintf(pointLog, "%lf ", Truemin[j]);
  //fclose(pointLog);
  ////------------------------------------------------------------------------------------------

  delete[] lower_bounds;
  delete[] upper_bounds;

  return min;
}
