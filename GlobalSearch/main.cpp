//#define OMPI_IMPORTS
#include "evolvent.h"
#include "problem_manager.h"
#include "extended.h"
#include "mpi.h"

#define _USE_MATH_DEFINES
#include <math.h>
#include <ctime>
#include <stdio.h>
#include <omp.h>
//#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include "functions.h"

//TEST(AGP, test_AGP_ra) {
//  double ra = 0;
//
//  TProblemManager manager;
//  std::string problemName = "ra";
//  manager.LoadProblemLibrary(problemName);
//  IProblem* problem = manager.GetProblem();
//  int dimension = problem->GetDimension();
//  double res = AGP(1000, 0.00001, dimension, problem);
//
//  EXPECT_EQ(ra, res);
//}


int  main(int argc, char* argv[])
{
  /// мэнаджер подключения задач
  TProblemManager manager;
  /// Число проведенных испытаний
  int trialCount = 0;
  /// Найдена ли точка глобального минимума
  bool isFindOptimalPoint = false;
  /// 
  //int iter = 0;
  /// Координаты точек испытания
  double** points = 0;
  /// Максимальное число испытаний
  int maxTrial = 1000;

  int dimension = 0;
  IProblem* problem;

  //Инициализация тестов и mpi
  //::testing::InitGoogleTest(&argc, argv);
  //MPI_Init(&argc, &argv);
  //int rankM;
  //MPI_Comm_rank(MPI_COMM_WORLD, &rankM);
  //::testing::AddGlobalTestEnvironment(new GTestMPIListener::MPIEnvironment);
  //::testing::TestEventListeners& listeners =
    //::testing::UnitTest::GetInstance()->listeners();

  //listeners.Release(listeners.default_result_printer());
  //listeners.Release(listeners.default_xml_generator());

  //listeners.Append(new GTestMPIListener::MPIMinimalistPrinter);

  if (argc > 1)
  {
    /// Имя задачи
    std::string problemName = argv[1];
    /// Конфигурационный файл
    std::string problepConf = "";
    if (argc > 2)
    {
      problepConf = argv[2];
    }

    int err;
    /// Загружаем задачу
    if (manager.LoadProblemLibrary(problemName) != TProblemManager::OK_)
      return 0;
    /// Получаем задачу
    problem = manager.GetProblem();
    /// Устанавливаем размерность
    err = problem->SetDimension(1);
    if (err != TProblemManager::OK_)
    {
      printf("Error SetDimension!\n");
      return 0;
    }
    /// Задаем файл с конфигурацией задачи
    err = problem->SetConfigPath(problepConf);
    if (err != TProblemManager::OK_)
    {
      printf("Error SetConfigPath!\n");
      return 0;
    }
    /// Инициализация
    err = problem->Initialize();
    if (err != TProblemManager::OK_)
    {
      printf("Error Initialize!\n");
      return 0;
    }
    /// Получаем размерность из задачи
    dimension = problem->GetDimension();
    /// Граници области поиска
    double* lower_bounds = new double[dimension];
    double* upper_bounds = new double[dimension];
    /// Получаем Граници области
    problem->GetBounds(lower_bounds, upper_bounds);

    /// Координаты найденного  минимума
    double* minx = new double[dimension];
    /// Значение минимума
    double minf = -std::numeric_limits<double>::infinity();

    if (argc > 3)
    {
      points = new double*[maxTrial * 2];
    }

    if (/*rankM == 0*/1)
    {

      //if (Extended::GetPrecision() == 0.01)
      //{
      //  if (parameters.m * (pTask->GetFreeN() - pTask->GetNumberOfDiscreteVariable()) <= 50)
      //  {
      //    Extended::SetTypeID(etDouble);
      //  }
      //  Extended::SetPrecision(1 / (::pow(2., parameters.m * (pTask->GetFreeN() -
      //    pTask->GetNumberOfDiscreteVariable()))));
      //}


      /////////////////////
      ///Вычисоения
      ////////////////////

      minx[0] = 0.5;
      minx[1] = 0;

      //    int n = 100;
      //
      //    Extended* ex = new Extended [n];
      //    Extended* ey = new Extended[n];
      //    Extended* ez = new Extended[n];
      //
      //    for (int i = 0; i < n; i++)
      //    {
      //      ex[i] = 1.0 / (i + 1.0);
      //      ey[i] = i * 100.0;
      //      ez[i] = (ex[i] + ey[i]) / 1000000.0;
      //
      //      //printf("ez[%d] = %lf\n", i, ez[i].toDouble());
      //    }
      //
      //    int t = 4;
      //#pragma omp parallel for num_threads(t)
      //    for (int i = 0; i < t; i++)
      //    {
      //      int thread = omp_get_thread_num();
      //      printf("%d\n", thread);
      //    }
      //
      //    Extended* exp = new Extended[n];
      //    Extended* eyp = new Extended[n];
      //    Extended* ezp = new Extended[n];
      //
      //#pragma omp parallel for num_threads(t)
      //    for (int i = 0; i < n; i++)
      //    {
      //      exp[i] = 1.0 / (i + 1.0);
      //      eyp[i] = i * 100.0;
      //      ezp[i] = (exp[i] + eyp[i]) / 1000000.0;
      //
      //      //printf("ez[%d] = %lf\n", i, ez[i].toDouble());
      //    }
      //
      //    for (int i = 0; i < n; i++)
      //    {
      //      if ((ex[i] != exp[i]) || (ey[i] != eyp[i]) || (ez[i] != ezp[i]))
      //        printf("ez[%d] = %lf\tezp[%d] = %lf\n", i, ez[i].toDouble(), i, ezp[i].toDouble());
      //    }
      //
      //
      minf = problem->CalculateFunctionals(minx, 0);

      double x = 0.5;
      GetInverseImage(minx, x, lower_bounds, upper_bounds, dimension, 10);
      GetImage(x, minx, lower_bounds, upper_bounds, dimension, 10);

      trialCount++;

      if (argc > 3)
      {
        points[0] = new double[2];

        points[0][0] = minx[0];
        points[0][1] = minx[1];
      }

      ////////////////////
      ///
      ///////////////////

      for (int i = 0; i < dimension; i++)
        printf("x[%d] = %lf\n", i, minx[i]);
      printf("min = %lf\n", minf);
      printf("Point count = %d\n\n", trialCount);

      /// Известная точка глобального минимума
      double* BestTrialy = new double[dimension];
      if (problem->GetOptimumPoint(BestTrialy) == TProblemManager::OK_)
      {
        for (int j = 0; j < dimension; j++)
          printf("Optimal x[%d] = %lf \n", j, BestTrialy[j]);
      }

      // Проверяем попадание в окрестность глобального оптимума
      double Epsilon = 0.01;
      isFindOptimalPoint = true;
      for (int i = 0; i < dimension; i++)
      {
        double fabsx = fabs(BestTrialy[i] - minx[i]);
        double fm = Epsilon * (upper_bounds[i] - lower_bounds[i]);
        if (fabsx > fm)
        {
          isFindOptimalPoint = false;
          break;
        }
      }

      if (isFindOptimalPoint)
      {
        printf("\nGlobal optimum FOUND!\n\n");
      }

      if (argc > 3)
      {
        /// Имя файла с координатами точек
        std::string pointLogName = argv[3];
        FILE* pointLog = fopen(pointLogName.c_str(), "w");
        /// Печатаем точки испытания
        fprintf(pointLog, "%d\n", trialCount);
        for (int i = 0; i < trialCount; i++)
        {
          for (int j = 0; j < dimension; j++)
            fprintf(pointLog, "%lf ", points[i][j]);
          fprintf(pointLog, "\n");
        }

        /// Печатаем найденную точку
        for (int j = 0; j < dimension; j++)
          fprintf(pointLog, "%lf ", minx[j]);
        fprintf(pointLog, "\n");

        /// Известная точка глобального минимума
        if (problem->GetOptimumPoint(BestTrialy) == TProblemManager::OK_)
        {
          for (int j = 0; j < dimension - 1; j++)
            fprintf(pointLog, "%lf ", BestTrialy[j]);
          fprintf(pointLog, "%lf", BestTrialy[dimension - 1]);
        }

        fclose(pointLog);
        delete[] points;


      }
    }
  }
  
  printf("\n");
  /*MPI_Barrier(MPI_COMM_WORLD);

  unsigned int start_time = clock();
  double per = Perebor(1000, dimension, problem);
  unsigned int end_time = clock();

  if (rankM == 0) 
  {
  double seconds = (double)(end_time - start_time) / CLOCKS_PER_SEC;
  printf("Runtime = %lf\n", seconds);

  printf("Try to do perebor = %lf\n", per);
  }
  MPI_Finalize();*/

  double start =  clock();
  double agp = AGP(10000, 0.000001, dimension, problem);
  double end = clock();
  double seconds = (end - start) / CLOCKS_PER_SEC;
  printf("Runtime AGP = %lf\n", seconds);
  printf("Try to do AGP = %lf\n", agp);

  //return RUN_ALL_TESTS();
  return 0;
}