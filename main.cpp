#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <ctime>

using namespace std;

const int N = 5; // Размерность матрицы
const float detCalc = 15540.; // Определитель матрицы A 5x5, вычисленный в WolframAlpha

// Метод Гаусса
float *solveWithGauss(float **A, float *b, int N) {

    auto *x = new float[N];

    // Прямой ход
    for (int i = 0; i < N; i++) {

        // Выбор главного элемента
        int maxRow = i;
        for (int k = i + 1; k < N; k++) {
            if (abs(A[k][i]) > abs(A[maxRow][i]))
                maxRow = k;
        }

        // Перестановка строк
        if (maxRow != i) {
            float *temp = A[i];
            A[i] = A[maxRow];
            A[maxRow] = temp;

            float tempB = b[i];
            b[i] = b[maxRow];
            b[maxRow] = tempB;
        }

        // Обнуление элементов ниже главного
        for (int k = i + 1; k < N; k++) {
            float factor = A[k][i] / A[i][i];
            b[k] -= factor * b[i];
            for (int j = i; j < N; j++) {
                A[k][j] -= factor * A[i][j];
            }
        }
    }

    // Обратный ход
    for (int i = N - 1; i >= 0; i--) {
        x[i] = b[i] / A[i][i];
        for (int j = i + 1; j < N; j++) {
            x[i] -= A[i][j] * x[j] / A[i][i];
        }
    }

    return x;
}

// Вычисление определителя
float computeDeterminant(float **A, int N) {

    float determinant = 1.0;

    for (int i = 0; i < N; i++) {
        determinant *= A[i][i];
    }

    return determinant;
}

// Метод Холецкого
float *solveWithCholesky(float **A, const float *b, int N) {

    // Разложение на нижнюю L и верхнюю U матрицы
    auto **L = new float *[N];
    auto **U = new float *[N];
    for (int i = 0; i < N; i++) {
        L[i] = new float[N];
        U[i] = new float[N];
    }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j <= i; j++) {
            float s = 0;
            for (int k = 0; k < j; k++) {
                s += L[i][k] * U[k][j];
            }
            L[i][j] = A[i][j] - s;
        }

        for (int j = i; j < N; j++) {
            float s = 0;
            for (int k = 0; k < i; k++) {
                s += L[i][k] * U[k][j];
            }
            U[i][j] = (A[i][j] - s) / L[i][i];
        }
    }

    // Решение Ly = b
    // L - нижняя треугольная матрица
    auto *y = new float[N];
    y[0] = b[0] / L[0][0];

    for (int i = 1; i < N; i++) {
        float s = 0;
        for (int j = 0; j < i; j++) {
            s += L[i][j] * y[j];
        }
        y[i] = (b[i] - s) / L[i][i];
    }

    // Решение Ux = y
    // U - верхняя треугольная матрица
    auto *x = new float[N];
    x[N - 1] = y[N - 1] / U[N - 1][N - 1];

    for (int i = N - 2; i >= 0; i--) {
        float s = 0;
        for (int j = i + 1; j < N; j++) {
            s += U[i][j] * x[j];
        }
        x[i] = (y[i] - s) / U[i][i];
    }

    // Освобождение памяти
    for (int i = 0; i < N; i++) {
        delete[] L[i];
        delete[] U[i];
    }
    delete[] L;
    delete[] U;
    delete[] y;

    return x;
}

// Вывод матрицы
void printMatrix(float **matrix, int N) {

//    cout << fixed << setprecision(4);

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            cout << setw(10) << matrix[i][j];
        }
        cout << endl << endl;
    }
}

void printMatrix(float **matrix, float *b, int N, int max_rows = 10, int max_cols = 10) {
    cout << fixed << setprecision(1);

    for (int i = 0; i < min(N, max_rows); i++) {
        for (int j = 0; j < min(N, max_cols); j++) {
            cout << setw(5) << matrix[i][j] << " ";
        }

        if (N > max_cols) {
            cout << " ...";
        }

        cout << " | " << b[i] << endl;
    }

    if (N > max_rows) {
        for (int i = 0; i < max_cols; i++) {
            cout << setw(5) << "..." << " ";
        }

        if (N > max_cols) {
            cout << " ...";
        }

        cout << " | ..." << endl;
    }

    cout << endl;
    cout << fixed << setprecision(8);
}

// Вывод вектора
void printVector(float *vector, int N) {

//    cout << fixed << setprecision(4);

    for (int i = 0; i < N; i++) {
        cout << "| " << vector[i] << " |" << endl;
    }
}

float computeMaxNorm(float **A, float *x, float *b, int N) {
    float maxNorm = 0.0;

    for (int i = 0; i < N; i++) {
        float sum = 0.0;
        for (int j = 0; j < N; j++) {
            sum += A[i][j] * x[j];
        }
        float residual = abs(sum - b[i]);
        if (residual > maxNorm) {
            maxNorm = residual;
        }
    }

    return maxNorm;
}

int main() {

    srand(time(nullptr));

    auto **A = new float *[N];
    for (int i = 0; i < N; i++) {
        A[i] = new float[N];
    }

    auto *b = new float[N];

    // Заполнение матрицы
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (i == j) {
                A[i][j] = 5 * sqrt(i + 1);
            } else {
                A[i][j] = sqrt((i + 1) * (j + 1));
            }
        }
    }

    // Заполнение вектора случайными значениями
    for (int i = 0; i < N; i++) {
        b[i] = (float) rand() / RAND_MAX;
    }

    cout << endl;
    cout << "Matrix:" << endl;
    printMatrix(A, N);

    cout << endl;
    cout << "Vector b:" << endl;
    printVector(b, N);

    // Решение системы методом Гаусса
    auto start_time = clock();
    float *xGauss = solveWithGauss(A, b, N);
    cout << endl;
    cout << "Solve with Gauss" << endl;
    printVector(xGauss, N);
    cout << "Work time of Gauss: " << (double)(clock() - start_time) / CLOCKS_PER_SEC << endl << endl;

    // Решение системы методом Холецкого
    start_time = clock();
    float *xCholesky = solveWithCholesky(A, b, N);
    cout << endl;
    cout << "Solve with Cholesky:" << endl;
    printVector(xCholesky, N);
    cout << "Work time of Cholesky: " << (double)(clock() - start_time) / CLOCKS_PER_SEC << endl << endl;

    // Вычисление максимум-нормы невязки для метода Гаусса
    float maxNormGauss = computeMaxNorm(A, xGauss, b, N);
    cout << endl;
    cout << "Max Norm Residual (Gauss): " << maxNormGauss;

    // Вычисление максимум-нормы невязки для метода Холецкого
    float maxNormCholesky = computeMaxNorm(A, xCholesky, b, N);
    cout << endl;
    cout << "Max Norm Residual (Cholesky): " << maxNormCholesky << endl;

    // Вычисление определителя
    float determinantA = computeDeterminant(A, N);
    cout << endl;
    cout << "Det A: " << determinantA << endl;
    cout << "DeFacto Det A: " << detCalc << endl;
    cout << "Absolut Error: " << abs(determinantA - detCalc) << endl;
    cout << "Relative Error: " << abs(determinantA - detCalc) / abs(determinantA) * 100 <<"%" << endl;


    // Освобождение памяти
    for (int i = 0; i < N; i++) {
        delete[] A[i];
    }
    delete[] A;
    delete[] b;
    delete[] xGauss;
    delete[] xCholesky;

    return 0;
}