#define PRECISION 10 // number of decimal places
#define A -3
#define B 6
#define FUNCTION pow(cos(x - 5), 2) + x * sqrt(e)
#define METHOD Solve_Gauss // Solve_Midpoint, Solve_Trapeze, Solve_Simpson, Solve_Gauss

#define e 2.71828182845904523536
#define pi 3.14159265358979323846

#define _STRINGIFY(s) #s
#define STRINGIFY(s) _STRINGIFY(s)

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
using namespace std;

double f(double x)
{
    return FUNCTION;
}

double Solve_Midpoint(double a, double b, unsigned int n, double h)
{
    double integral = 0;
    for (int i = 0; i < n; i++)
    {
        double x = a + h * i;
        integral += f(x - 0.5 * h);
    }  
    integral *= h;
    return integral;
}

double Solve_Trapeze(double a, double b, unsigned int n, double h)
{
    double integral = 0;
    for (int i = 1; i <= n; i++)
    {
        double x = a + h * i;
        integral += (f(x - h) + f(x)) / 2;
    }
    integral *= h;
    return integral;
}

double Solve_Simpson(double a, double b, unsigned int n, double h)
{
    double integral = 0;
    for (int i = 1; i <= n; i++)
    {
        double x = a + h * i;
        integral += f(x - h) + 4 * f(x - 0.5 * h) + f(x);
    }
    integral *= h / 6;
    return integral;
}

double Solve_Gauss(double a, double b, unsigned int N, double h)
{
    // Polynomial degree = 3
    const int n = 3;
    const double x[] = { -0.7745966692, 0.0, 0.7745966692 };
    const double q[] = { 0.5555555556, 0.8888888888, 0.5555555556 };

    double sum_j = 0;
    for (int j = 1; j <= n; j++)
    {
        double sum_i = 0;
        for (int i = 1; i <= N; i++)
        {
            double t_i = a + i * h;
            double t_ij = (t_i - h + t_i) / 2 + h / 2 * x[j - 1];
            sum_i += f(t_ij);
        }
        sum_j += q[j - 1] * sum_i;
    }
    sum_j *= h / 2;
    return sum_j;
}

int main()
{
    setlocale(LC_ALL, "");

    cout << "Chosen solution method: " << STRINGIFY(METHOD) << "\n";
    cout << "Precision: " << PRECISION << " decimal places\n";
    cout << "Lower integration limit: " << A << "\n";
    cout << "Upper integration limit: " << B << "\n";
    cout << "Function: f(x) = " << STRINGIFY(FUNCTION) << "\n\n";

    unsigned int iterCount = 1;
    double eps = 1 / pow(10, PRECISION + 1);
    unsigned int n = 2;
    double h = abs(double(B - A)) / n;

    double I_1N = METHOD(A, B, n, h);
    double integral;
    while (true)
    {
        n *= 2;
        h /= 2;
        cout << "Iteration " << setw(2) << iterCount << ": ";

        double I_2N = METHOD(A, B, n, h);
        if (abs(I_1N - I_2N) < eps)
        {
            integral = I_2N;
            cout << "integral has been solved (" << setprecision(17) << I_2N << ")\n";
            break;
        }
        else
        {
            I_1N = I_2N;
            iterCount++;
            cout << "insufficient precision   (" << setprecision(17) << I_2N << ")\n";
        }
    }

    cout << "\nThe integral of the given function from " << A << " to " << B;
    cout << "\nwith a precision of " << PRECISION << " decimal places";
    cout << " is " << fixed << setprecision(PRECISION) << integral;
    cout << "\nIterations count: " << iterCount << "\n";
}