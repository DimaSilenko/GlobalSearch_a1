#ifndef EVOLVENT_H
#define  EVOLVENT_H

#define MaxDim 50

// ������ �������� _y (������������� ��������� �� ������ x (�� ������� 0..1)
// ������ ������� ������� ������ - A
// ������� ������� ������� ������ - B
// ����������� - N
// ��������� ���������� ��������� - m
void GetImage(double x, double* _y, const double* A, const double* B, int N, int m);
// ������ ����� � (�� ������� 0..1) �� ���������  _y (������������� ���������)
// ������ ������� ������� ������ - A
// ������� ������� ������� ������ - B
// ����������� - N
// ��������� ���������� ��������� - m
void GetInverseImage( double* _y, double& x, const double* A, const double* B, int N, int m);

#endif //EVOLVENT_H