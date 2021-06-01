//#define OMPI_IMPORTS
//#pragma comment(lib, "D:\\GlobalSearch\\Test\\GlobalSearch_a1\\Drawing\\dislin\\lib\\win\\x86\\discpp.lib")

#include "evolvent.h"
#include "problem_manager.h"
#include "extended.h"
#include "mpi.h"
#include "pugixml.hpp"
#include "isolinesPlotter.h"

#define _USE_MATH_DEFINES
#include <math.h>
#include <ctime>
#include <stdio.h>
#include <omp.h>
#include <gtest/gtest.h>
#include "functions.h"


int  main(int argc, char* argv[])
{
  double eps = 0;
  double r = 0;
  /// �������� ����������� �����
  TProblemManager manager;
  /// ����� ����������� ���������
  int trialCount = 0;
  /// ������� �� ����� ����������� ��������
  bool isFindOptimalPoint = false;
  /// 
  //int iter = 0;
  /// ���������� ����� ���������
  double** points = 0;
  /// ������������ ����� ���������
  int maxTrial = 1000;

  int dimension = 0;
  IProblem* problem;

  //������������� mpi
  MPI_Init(&argc, &argv);
  int rankM;
  MPI_Comm_rank(MPI_COMM_WORLD, &rankM);

  /// ��������� ����� ����������� ��������
  double* BestTrialy;

  if (argc > 1)
  {
    /// ��� ������
    std::string problemName = argv[1];
    /// ���������������� ����
    std::string problepConf = "";
    if (argc > 2)
    {
      problepConf = argv[2];
    }

    int err;
    /// ��������� ������
    if (manager.LoadProblemLibrary(problemName) != TProblemManager::OK_)
      return 0;
    /// �������� ������
    problem = manager.GetProblem();
    /// ������������� �����������
    err = problem->SetDimension(2);
    if (err != TProblemManager::OK_)
    {
      printf("Error SetDimension!\n");
      return 0;
    }
    /// ������ ���� � ������������� ������
    err = problem->SetConfigPath(problepConf);
    if (err != TProblemManager::OK_)
    {
      printf("Error SetConfigPath!\n");
      return 0;
    }
    /// �������������
    err = problem->Initialize();
    if (err != TProblemManager::OK_)
    {
      printf("Error Initialize!\n");
      return 0;
    }
    /// �������� ����������� �� ������
    dimension = problem->GetDimension();
    printf("Dimension = %d\n", dimension);

    pugi::xml_document doc;
    pugi::xml_parse_result result = doc.load_file(problepConf.c_str());
    pugi::xml_node config = doc.child("config");
    
    eps = std::stod(config.child("epsilon").child_value());
    printf("Epsilon = %lf\n", eps);
    r = std::stod(config.child("r").child_value());
    //��������� ����� ������
    int funcNumber = std::stoi(config.child("function_number").child_value());
    printf("Function_number = %d\n", problem->GetNumberOfFunctions());
    //��� GKLS-----------------------------------------------------------------------------------------------
    double global_dist = std::stod(config.child("global_dist").child_value());
    double global_radius = std::stod(config.child("global_radius").child_value());
    int numMinuma = std::stoi(config.child("num_minima").child_value());


    printf("Task_number = %d\n", funcNumber);
    printf("global_dist = %lf\n", global_dist);
    printf("global_radius = %lf\n", global_radius);
    printf("numMinuma = %d\n\n", numMinuma);
    /// ������� ������� ������
    double* lower_bounds = new double[dimension];
    double* upper_bounds = new double[dimension];
    /// �������� ������� �������
    problem->GetBounds(lower_bounds, upper_bounds);
    /// ���������� ����������  ��������
    double* minx = new double[dimension];
    /// �������� ��������
    double minf = -std::numeric_limits<double>::infinity();

    if (argc > 3)
    {
      points = new double*[maxTrial * 2];
    }

    if (rankM == 0)
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
      ///����������
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

      //-----------------------------------------------------------------------------------------------------------------
      /*for (int i = 0; i < dimension; i++)
        printf("x[%d] = %lf\n", i, minx[i]);
      printf("min = %lf\n", minf);
      printf("Point count = %d\n\n", trialCount);*/
      //-----------------------------------------------------------------------------------------------------------------

      BestTrialy = new double[dimension];

      if (problem->GetOptimumPoint(BestTrialy) == TProblemManager::OK_)
      {
        for (int j = 0; j < dimension; j++)
          printf("Optimal x[%d] = %lf \n", j, BestTrialy[j]);
      }

      // ��������� ��������� � ����������� ����������� ��������
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
        /// ��� ����� � ������������ �����
        std::string pointLogName = argv[3];
        FILE* pointLog = fopen(pointLogName.c_str(), "w");
        /// �������� ����� ���������
        fprintf(pointLog, "%d\n", trialCount);
        for (int i = 0; i < trialCount; i++)
        {
          for (int j = 0; j < dimension; j++)
            fprintf(pointLog, "%lf ", points[i][j]);
          fprintf(pointLog, "\n");
        }

        /// �������� ��������� �����
        for (int j = 0; j < dimension; j++)
          fprintf(pointLog, "%lf ", minx[j]);
        fprintf(pointLog, "\n");

        /// ��������� ����� ����������� ��������
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
 
  MPI_Barrier(MPI_COMM_WORLD);
  int iter = 1000000;

  double* bestX = new double[dimension];
  double start_time = MPI_Wtime();
  double per = 0;// Perebor(iter, dimension, problem, &bestX);
  double end_time = MPI_Wtime();

  double By = problem->CalculateFunctionals(BestTrialy, 0);
  printf("min = %lf\n", By);

  if (rankM == 0) 
  {
    //������� ��� ���� �� ����� ������
    /*printf("\nPerebor:\n");

    double seconds = end_time - start_time;
    printf("Iteration Perebor = %d\n", iter);
    printf("Solve time Perebor = %lf\n", seconds);
    printf("Try to do perebor: Function = %lf X = (%lf, %lf)\n", per, bestX[0], bestX[1]);
    printf("Accuracy Function Perebor = %lf, Accuracy X Perebor = (%lf, %lf)\n", (abs(per - By)),
      (abs(bestX[0] - BestTrialy[0])), (abs(bestX[1] - BestTrialy[1])));

    bool answer = ((abs(bestX[0] - BestTrialy[0]) < eps) && (abs(bestX[1] - BestTrialy[1]) < eps));
    if (answer)
      printf("Perebor found Global optimum!\n");
    else
      printf("Something went wrong with Perebor!!\n");*/

    int iteration = 0;
    printf("\n\nAGP:\n");
    printf("r = %lf\n", r);
    int numThr = 2;
    int flag = 1;
    clock_t start = clock();
    double agp = AGP_Space(iter, eps, dimension, problem, &bestX, &iteration, r, BestTrialy, numThr, flag, argv);
    clock_t end = clock();
    
    double seconds = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Iteration AGP = %d\n", iteration);
    printf("Solve time AGP = %lf\n", seconds);
    printf("Try to do AGP: Function = %lf X = (%lf, %lf)\n", agp, bestX[0], bestX[1]);

    printf("Accuracy Function AGP = %lf, Accuracy X AGP = (%lf, %lf)\n", (abs(agp - By)),
      (abs(bestX[0] - BestTrialy[0])), (abs(bestX[1] - BestTrialy[1])));

    bool answer = ((abs(bestX[0] - BestTrialy[0]) < eps) && (abs(bestX[1] - BestTrialy[1]) < eps));
    if (answer)
      printf("AGP found Global optimum!\n");
    else
      printf("Something went wrong with AGP!!\n");

    //������������ �����, ��� ������� ���, ����� ������ ��� �� ���
    /*plot2dProblemIsolines(problem, "D:\\GlobalSearch\\Test\\GlobalSearch_a1\\points.png",
      "D:\\GlobalSearch\\Test\\GlobalSearch_a1\\points.txt", false, 0);*/
  }
  MPI_Finalize();

  return 0;
}