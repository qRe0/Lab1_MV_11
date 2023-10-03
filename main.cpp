#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <ctime>

using namespace std;

const int N = 3;

void printMatrix(double** matrix, int N);
void printVector(double* vector, int N);

double* solveWithGauss(double** A, double* b, int N);
double* solveWithCholesky(double** A, double* b, int N);
double computeDeterminant(double** A, int N);

int main() {

    srand(time(NULL));

    double** A = new double*[N];
    for(int i = 0; i < N; i++) {
        A[i] = new double[N];
    }

    double* b = new double[N];

    // Заполнение матрицы
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++) {
            if(i == j) {
                A[i][j] = 5 * sqrt(i+1);
            } else {
                A[i][j] = sqrt((i+1)*(j+1));
            }
        }
    }

    // Заполнение вектора случайными значениями
    for(int i = 0; i < N; i++) {
        b[i] = (double)rand() / RAND_MAX;
    }

    cout << "Matrix:" << endl;
    printMatrix(A, N);

    cout << "Vector b:" << endl;
    printVector(b, N);

    // Решение системы методом Гаусса
    double* xGauss = solveWithGauss(A, b, N);

    // Вычисление определителя
    double determinantA = computeDeterminant(A, N);

    cout << "Solve with Gauss" << endl;
    printVector(xGauss, N);

    cout << "Det A: " << determinantA << endl;

    // Решение системы методом Холецкого
    double* xCholesky = solveWithCholesky(A, b, N);

    cout << "Solve with Cholesky:" << endl;
    printVector(xCholesky, N);

    // Освобождение памяти
    for(int i = 0; i < N; i++) {
        delete[] A[i];
    }
    delete[] A;
    delete[] b;
    delete[] xGauss;
    delete[] xCholesky;

    return 0;
}

// Метод Гаусса
double* solveWithGauss(double** A, double* b, int N) {

    double* x = new double[N];

    // Прямой ход
    for(int i = 0; i < N; i++) {

        // Выбор главного элемента
        int maxRow = i;
        for(int k = i+1; k < N; k++) {
            if(abs(A[k][i]) > abs(A[maxRow][i]))
                maxRow = k;
        }

        // Перестановка строк
        if(maxRow != i) {
            double* temp = A[i];
            A[i] = A[maxRow];
            A[maxRow] = temp;

            double tempB = b[i];
            b[i] = b[maxRow];
            b[maxRow] = tempB;
        }

        // Обнуление элементов ниже главного
        for(int k = i+1; k < N; k++) {
            double factor = A[k][i] / A[i][i];
            b[k] -= factor * b[i];
            for(int j = i; j < N; j++) {
                A[k][j] -= factor * A[i][j];
            }
        }
    }

    // Обратный ход
    for(int i = N-1; i >= 0; i--) {
        x[i] = b[i] / A[i][i];
        for(int j = i+1; j < N; j++) {
            x[i] -= A[i][j] * x[j] / A[i][i];
        }
    }

    return x;
}

// Вычисление определителя
double computeDeterminant(double** A, int N) {

    double determinant = 1.0;

    for(int i = 0; i < N; i++) {
        determinant *= A[i][i];
    }

    return determinant;
}

// Метод Холецкого
double* solveWithCholesky(double** A, double* b, int N) {

    // Разложение на нижнюю L и верхнюю U матрицы
    double** L = new double*[N];
    double** U = new double*[N];
    for(int i = 0; i < N; i++) {
        L[i] = new double[N];
        U[i] = new double[N];
    }

    for(int i = 0; i < N; i++) {
        for(int j = 0; j <= i; j++) {
            double s = 0;
            for(int k = 0; k < j; k++) {
                s += L[i][k] * U[k][j];
            }
            L[i][j] = A[i][j] - s;
        }

        for(int j = i; j < N; j++) {
            double s = 0;
            for(int k = 0; k < i; k++) {
                s += L[i][k] * U[k][j];
            }
            U[i][j] = (A[i][j] - s) / L[i][i];
        }
    }

    // Решение Ly = b
    double* y = new double[N];
    y[0] = b[0] / L[0][0];

    for(int i = 1; i < N; i++) {
        double s = 0;
        for(int j = 0; j < i; j++) {
            s += L[i][j] * y[j];
        }
        y[i] = (b[i] - s) / L[i][i];
    }

    // Решение Ux = y
    double* x = new double[N];
    x[N-1] = y[N-1] / U[N-1][N-1];

    for(int i = N-2; i >= 0; i--) {
        double s = 0;
        for(int j = i+1; j < N; j++) {
            s += U[i][j] * x[j];
        }
        x[i] = (y[i] - s) / U[i][i];
    }

    // Освобождение памяти
    for(int i = 0; i < N; i++) {
        delete[] L[i];
        delete[] U[i];
    }
    delete[] L;
    delete[] U;
    delete[] y;

    return x;
}

// Вывод матрицы
void printMatrix(double** matrix, int N) {

    cout << fixed << setprecision(4);

    for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++) {
            cout << matrix[i][j] << " ";
        }
        cout << endl;
    }
}

// Вывод вектора
void printVector(double* vector, int N) {

    cout << fixed << setprecision(4);

    for(int i = 0; i < N; i++) {
        cout << vector[i] << endl;
    }
}