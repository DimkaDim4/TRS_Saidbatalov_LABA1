#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;

double dxdt(double t, double x, double v)
{
    return v;
}

double dvdt(double t, double x, double v)
{
    return -2. / 9.1 * x * (1 - 2 * x * x);
}

double dxdt_SM(double x, double u, double v)
{
    return v;
}

double dvdt_SM(double x, double u, double v)
{
    return -u + 4. * exp(x);
}

void Euler(double* t, double* x, double* v, double x_0, double v_0, double a, double b, int N)
{
    double h = (b - a) / N;

    for (int i = 0; i < N + 1; i++)
    {
        t[i] = i * h;
    }

    x[0] = x_0; v[0] = v_0;

    for (int i = 0; i < N; i++) {
        x[i + 1] = x[i] + h * dxdt(t[i], x[i], v[i]);
        v[i + 1] = v[i] + h * dvdt(t[i], x[i], v[i]);
    }
}

void Adams(double* t, double* x, double* v, double x_0, double v_0, double a, double b, int N)
{
    double h = (b - a) / N;

    for (int i = 0; i < N + 1; i++)
    {
        t[i] = i * h;
    }

    x[0] = x_0; v[0] = v_0;
    x[1] = x[0] + h * dxdt(t[0], x[0], v[0]);
    v[1] = v[0] + h * dvdt(t[0], x[0], v[0]);

    for (int i = 1; i < N; i++)
    {
        x[i + 1] = x[i] + h / 2. * (3. * dxdt(t[i], x[i], v[i]) - dxdt(t[i - 1], x[i - 1], v[i - 1]));
        v[i + 1] = v[i] + h / 2. * (3. * dvdt(t[i], x[i], v[i]) - dvdt(t[i - 1], x[i - 1], v[i - 1]));
    }
}

void RungeKutta4(double* t, double* x, double* v, double x_0, double v_0, double a, double b, int N)
{
    double h = (b - a) / N;

    for (int i = 0; i < N + 1; i++)
    {
        t[i] = i * h;
    }

    x[0] = x_0; v[0] = v_0;

    double k1, k2, k3, k4, n1, n2, n3, n4;

    for (int i = 0; i < N; i++) {
        k1 = dxdt(t[i], x[i], v[i]);
        n1 = dvdt(t[i], x[i], v[i]);

        k2 = dxdt(t[i] + h / 2., x[i] + h * k1 / 2., v[i] + h * n1 / 2.);
        n2 = dvdt(t[i] + h / 2., x[i] + h * k1 / 2., v[i] + h * n1 / 2.);

        k3 = dxdt(t[i] + h / 2., x[i] + h * k2 / 2., v[i] + h * n2 / 2.);
        n3 = dvdt(t[i] + h / 2., x[i] + h * k2 / 2., v[i] + h * n2 / 2.);

        k4 = dxdt(t[i] + h, x[i] + h * k3, v[i] + h * n3);
        n4 = dvdt(t[i] + h, x[i] + h * k3, v[i] + h * n3);

        x[i + 1] = x[i] + (h / 6.) * (k1 + 2. * k2 + 2. * k3 + k4);
        v[i + 1] = v[i] + (h / 6.) * (n1 + 2. * n2 + 2. * n3 + n4);
    }
}

void write(std::string str, double* t, double* x, double* v, int N)
{
    ofstream file;
    file.open(str);
    for (int i = 0; i < N; i++)
    {
        file << t[i] << "\t" << x[i] << "\t" << v[i] << endl;
    };
    file.close();
}

void write(std::string str, double* t, double* x, int N)
{
    ofstream file;
    file.open(str);
    for (int i = 0; i < N; i++)
    {
        file << t[i] << "\t" << x[i] << endl;
    };
    file.close();
}

void read(std::string str, double* x)
{
    ifstream _file_1(str);
    double x_i;
    int i = 0;
    if (_file_1.is_open())
    {
        while (_file_1 >> x_i)
        {
            x[i] = x_i;
            i++;
        }
        _file_1.close();
    }
}

double aitken_process(std::string str, double* x_1, double* x_2, double* x_3, int N, int k)
{
    ofstream file;
    double delta_max_1 = 0;
    double delta_max_2 = 0;

    file.open(str);
    for (int i = 0; i < N + 1; i++)
    {
        double delta_1 = fabs(x_1[i] - x_2[i * k]);
        double delta_2 = fabs(x_2[i * k] - x_3[i * k * k]);
        if (delta_1 > delta_max_1)
            delta_max_1 = delta_1;
        if (delta_2 > delta_max_2)
            delta_max_2 = delta_2;
    }
    file << log(N) << "\t" << log(delta_max_1) << endl;
    file << log(N * k) << "\t" << log(delta_max_2) << endl;
    file.close();
    double p = (log(delta_max_2) - log(delta_max_1)) / (log(N * k) - log(N));
    return p;
}

void Runge(double* x_1, double* x_2, double* res, int N, int k, double p)
{
    double G = pow(k, p) / (pow(k, p) - 1.);
    for (int i = 0; i < N + 1; i++)
    {
        res[i] = G * x_1[i] + (1. - G) * x_2[i * k];
    }
}

double residu(double* x_1, double* x_2, double* res, int N, int k)
{
    double delta = 0.;
    for (int i = 0; i < N + 1; i++)
    {
        res[i] = x_1[i] - x_2[i * k];
        if (fabs(res[i]) > delta)
            delta = fabs(res[i]);
    }
    return delta;
}

void sweep(double* a, double* b, double* c, double* d, double* res, int N) // метод прогонки
{
    double* alpha = new double[N];
    double* beta = new double[N + 1];
    double* y = new double[N + 1];

    y[0] = b[0];
    alpha[0] = -c[0] / y[0];
    beta[0] = d[0] / y[0];

    for (int i = 1; i < N; i++)
    {
        y[i] = b[i] + a[i] * alpha[i - 1];
        alpha[i] = -c[i] / y[i];
        beta[i] = (d[i] - a[i] * beta[i - 1]) / y[i];
    }
    y[N] = b[N] + a[N] * alpha[N - 1];
    beta[N] = (d[N] - a[N] * beta[N]) / y[N];

    res[N] = beta[N];
    for (int i = N - 1; i >= 0; i--)
    {
        res[i] = alpha[i] * res[i + 1] + beta[i];
    }
}

void FDM(double* x, double* u, double A, double B, int N) // конечно-разностный метод
{
    double* a = new double[N + 1];
    double* b = new double[N + 1];
    double* c = new double[N + 1];
    double* d = new double[N + 1];

    double h = (B - A) / N;
    x[0] = A;
    b[0] = -1. / h + (h - 2.) / 2;
    c[0] = 1. / h;
    d[0] = h * 2. * exp(x[0]) - 4.;

    for (int i = 1; i < N + 1; i++)
    {
        x[i] = x[i - 1] + h;
    }

    for (int i = 1; i < N; i++)
    {
        a[i] = 1. / (h * h);
        b[i] = 1. - 2. * a[i];
        c[i] = a[i];
        d[i] = 4. * exp(x[i]);
    }
    a[N] = 0.;
    b[N] = 1.;
    d[N] = -3.;

    sweep(a, b, c, d, u, N);
}

void RungeKutta_SM(double* x, double* u, double u_0, double v_0, double a, double b, int N)
{
    double* v = new double[N + 1];
    double h = (b - a) / N;
    u[0] = u_0; v[0] = v_0;

    double k1, k2, k3, k4, n1, n2, n3, n4;

    for (int i = 0; i < N; i++) {
        k1 = dxdt_SM(x[i], u[i], v[i]);
        n1 = dvdt_SM(x[i], u[i], v[i]);

        k2 = dxdt_SM(x[i] + h / 2., u[i] + h * k1 / 2., v[i] + h * n1 / 2.);
        n2 = dvdt_SM(x[i] + h / 2., u[i] + h * k1 / 2., v[i] + h * n1 / 2.);

        k3 = dxdt_SM(x[i] + h / 2., u[i] + h * k2 / 2., v[i] + h * n2 / 2.);
        n3 = dvdt_SM(x[i] + h / 2., u[i] + h * k2 / 2., v[i] + h * n2 / 2.);

        k4 = dxdt_SM(x[i] + h / 2., u[i] + h * k3, v[i] + h * n3);
        n4 = dvdt_SM(x[i] + h / 2., u[i] + h * k3, v[i] + h * n3);

        u[i + 1] = u[i] + (h / 6.) * (k1 + 2. * k2 + 2. * k3 + k4);
        v[i + 1] = v[i] + (h / 6.) * (n1 + 2. * n2 + 2. * n3 + n4);
    }
    delete[] v;
}

void SM(double* x, double* u_1, double a, double b, int N, double eps) // метод стрельбы
{
    double eta_1 = 0.;
    double eta_2 = 0.;
    double eta_3 = 0.;
    double* u_2 = new double[N + 1];
    double* u_3 = new double[N + 1];
    double h = (b - a) / N;
    x[0] = a;
    for (int i = 1; i < N + 1; i++)
    {
        x[i] = x[i - 1] + h;
    }

    h = 0.1;

    do
    {
        if (eta_1 >=0 )
            eta_1 += h;
        if (eta_1 > 100.)
        {
            eta_1 = -h;
        }

        if (eta_1 < 0)
            eta_1 -= h;
        if (eta_1 < -100.)
        {
            throw exception("Parametr eta_1 not found");
        }
        RungeKutta_SM(x, u_1, eta_1, eta_1 - 4., a, b, N);
    } while ((u_1[N] + 3.) < 0.);

    do
    {
        if (eta_2 >= 0.)
            eta_2 += h;
        if (eta_2 > 100.)
        {
            eta_2 = -h;
        }

        if (eta_2 < 0.)
            eta_2 -= h;
        if (eta_2 < -100.)
        {
            throw exception("Parametr eta_2 not found");
        }
        RungeKutta_SM(x, u_1, eta_2, eta_2 - 4., a, b, N);
    } while ((u_1[N] + 3.) > 0.);

    if (eta_1 > eta_2)
    {
        swap(eta_1, eta_2);
    }

    do
    {
        RungeKutta_SM(x, u_1, eta_1, eta_1 - 4., a, b, N);
        RungeKutta_SM(x, u_2, eta_2, eta_2 - 4., a, b, N);

        eta_3 = (eta_1 + eta_2) / 2.;
        RungeKutta_SM(x, u_3, eta_3, eta_3 - 4., a, b, N);

        double res = u_3[N];

        if ((u_1[N] + 3.) * (u_3[N] + 3.) <= 0.)
            eta_2 = eta_3;
        if ((u_3[N] + 3.) * (u_2[N] + 3.) <= 0.)
            eta_1 = eta_3;
    } while (fabs(u_3[N] + 3.) >= eps);

    delete[] u_2;
    delete[] u_3;
}

int main()
{
    double a = 0., b = 50.;
    double x_0 = 0., v_0 = 0.1;
    double p;
    int N = 500;
    int k = 10;

    //double* res_1 = new double[N + 1];
    //double* res_2 = new double[N * k + 1];
    //double* res_3 = new double[N * k * k + 1];
    //double* residu_1 = new double[N + 1];
    //double* residu_2 = new double[N * k + 1];
    //double* residu_3 = new double[N * k * k + 1];
    //double* runge = new double[N + 1];
    //ofstream file;

    //read("Solution1.txt", res_1); //501
    //read("Solution2.txt", res_2); //5001
    //read("Solution3.txt", res_3); //50001

    //double* t_1E = new double[N + 1], * t_2E = new double[N * k + 1], * t_3E = new double[N * k * k + 1];
    //double* x_1E = new double[N + 1], * x_2E = new double[N * k + 1], * x_3E = new double[N * k * k + 1];
    //double* v_1E = new double[N + 1], * v_2E = new double[N * k + 1], * v_3E = new double[N * k * k + 1];

    //double* t_1 = new double[N + 1], * t_2 = new double[N * k + 1], * t_3 = new double[N * k * k + 1];
    //double* x_1 = new double[N + 1], * x_2 = new double[N * k + 1], * x_3 = new double[N * k * k + 1];
    //double* v_1 = new double[N + 1], * v_2 = new double[N * k + 1], * v_3 = new double[N * k * k + 1];

    /////////////////////////////////////////////////////////
    //// задание 1
    /////////////////////////////////////////////////////////
    //Euler(t_1E, x_1E, v_1E, x_0, v_0, a, b, N);
    //Euler(t_2E, x_2E, v_2E, x_0, v_0, a, b, N * k);
    //Euler(t_3E, x_3E, v_3E, x_0, v_0, a, b, N * k * k);

    //write("Euler1.txt", t_1E, x_1E, N + 1);
    //write("Euler2.txt", t_2E, x_2E, N * k + 1);
    //write("Euler3.txt", t_3E, x_3E, N * k * k + 1);

    //p = aitken_process("DELTA_Euler.txt", x_1E, x_2E, x_3E, N, k);
    //cout << "Euler: p = " << p <<endl;

    //Runge(x_1E, x_2E, runge, N, k, p);
    //write("Runge_Euler.txt", t_1E, runge, N + 1);

    //double delta_1 = residu(x_1E, res_3, residu_1, N, k * k);
    //double delta_2 = residu(x_2E, res_3, residu_2, N * k, k);
    //double delta_3 = residu(x_3E, res_3, residu_3, N * k * k, 1);

    //file.open("Delta_Euler.txt");
    //file << N << "\t" << delta_1 << endl;
    //file << N * k << "\t" << delta_2 << endl;
    //file << N * k * k << "\t" << delta_3 << endl;
    //file.close();

    //write("Residu_Euler_1.txt", t_1E, residu_1, N + 1);
    //write("Residu_Euler_2.txt", t_2E, residu_2, N * k + 1);
    //write("Residu_Euler_3.txt", t_3E, residu_3, N * k * k + 1);

    /////////////////////////////////////////////////////////
    //// задание 2
    /////////////////////////////////////////////////////////
    //Adams(t_1, x_1, v_1, x_0, v_0, a, b, N);
    //Adams(t_2, x_2, v_2, x_0, v_0, a, b, N * k);
    //Adams(t_3, x_3, v_3, x_0, v_0, a, b, N * k * k);

    //write("Adams_1.txt", t_1, x_1, N + 1);
    //write("Adams_2.txt", t_2, x_2, N * k + 1);
    //write("Adams_3.txt", t_3, x_3, N * k * k + 1);

    //p = aitken_process("DELTA_Adams.txt", x_1, x_2, x_3, N, k);
    //cout << "Adams: p = " << p << endl;

    //delta_1 = residu(x_1, x_1E, residu_1, N, 1);
    //delta_2 = residu(x_2, x_2E, residu_2, N * k, 1);
    //delta_3 = residu(x_3, x_3E, residu_3, N * k * k, 1);

    //write("Residu_Adams_Euler_1.txt", t_1, residu_1, N + 1);
    //write("Residu_Adams_Euler_2.txt", t_2, residu_2, N * k + 1);
    //write("Residu_Adams_Euler_3.txt", t_3, residu_3, N * k * k + 1);

    //file.open("Delta_Euler_Adams.txt");
    //file << N << "\t" << delta_1 << endl;
    //file << N * k << "\t" << delta_2 << endl;
    //file << N * k * k << "\t" << delta_3 << endl;
    //file.close();

    //residu(x_1, res_3, residu_1, N, k * k);
    //residu(x_2, res_3, residu_2, N * k, k);
    //residu(x_3, res_3, residu_3, N * k * k, 1);

    //write("Residu_Adams_1.txt", t_1, residu_1, N + 1);
    //write("Residu_Adams_2.txt", t_2, residu_2, N * k + 1);
    //write("Residu_Adams_3.txt", t_3, residu_3, N * k * k + 1);

    /////////////////////////////////////////////////////////
    //// задание 3
    /////////////////////////////////////////////////////////
    //RungeKutta4(t_1, x_1, v_1, x_0, v_0, a, b, N);
    //RungeKutta4(t_2, x_2, v_2, x_0, v_0, a, b, N * k);
    //RungeKutta4(t_3, x_3, v_3, x_0, v_0, a, b, N * k * k);

    //write("RungeKutta_1.txt", t_1, x_1, N + 1);
    //write("RungeKutta_2.txt", t_2, x_2, N * k + 1);
    //write("RungeKutta_3.txt", t_3, x_3, N * k * k + 1);

    //p = aitken_process("DELTA_Runge.txt", x_1, x_2, x_3, N, k);
    //cout << "Runge: p = " << p << endl;

    //delta_1 = residu(x_1, res_3, residu_1, N, k * k);
    //delta_2 = residu(x_2, res_3, residu_2, N * k, k);
    //delta_3 = residu(x_3, res_3, residu_3, N * k * k, 1);

    //file.open("Delta_RungeKutta.txt");
    //file << N << "\t" << delta_1 << endl;
    //file << N * k << "\t" << delta_2 << endl;
    //file << N * k * k << "\t" << delta_3 << endl;
    //file.close();

    //write("Residu_RungeKutta_1.txt", t_1, residu_1, N + 1);
    //write("Residu_RungeKutta_2.txt", t_2, residu_2, N * k + 1);
    //write("Residu_RungeKutta_3.txt", t_3, residu_3, N * k * k + 1);

    //delete[] res_1;
    //delete[] res_2;
    //delete[] res_3;
    //delete[] residu_1;
    //delete[] residu_2;
    //delete[] residu_3;
    //delete[] runge;
    //delete[] t_1;
    //delete[] t_2;
    //delete[] t_3;
    //delete[] x_1;
    //delete[] x_2;
    //delete[] x_3;
    //delete[] v_1;
    //delete[] v_2;
    //delete[] v_3;
    //delete[] t_1E;
    //delete[] t_2E;
    //delete[] t_3E;
    //delete[] x_1E;
    //delete[] x_2E;
    //delete[] x_3E;
    //delete[] v_1E;
    //delete[] v_2E;
    //delete[] v_3E;

    ///////////////////////////////////////////////////////
    // задание 4
    ///////////////////////////////////////////////////////
    N = 500;
    k = 4;
    a = -20.;
    b = -1.;
    double* x_1 = new double[N + 1], * x_2 = new double[N * k + 1], * x_3 = new double[N * k * k + 1];
    double* u_1 = new double[N + 1], * u_2 = new double[N * k + 1], * u_3 = new double[N * k * k + 1];

    double* residu_1 = new double[N + 1];
    double* residu_2 = new double[N * k + 1];
    double* residu_3 = new double[N * k * k + 1];

    double* solution_1 = new double[N + 1];
    double* solution_2 = new double[N * k + 1];
    double* solution_3 = new double[N * k * k + 1];

    read("solu1.txt", solution_1); //501
    read("solu2.txt", solution_2); //2001
    read("solu3.txt", solution_3); //8001

    FDM(x_1, u_1, a, b, N);
    FDM(x_2, u_2, a, b, N * k);
    FDM(x_3, u_3, a, b, N * k * k);

    write("FDM_500.txt", x_1, u_1, N + 1);
    write("FDM_2000.txt", x_2, u_2, N * k + 1);
    write("FDM_8000.txt", x_3, u_3, N * k * k + 1);

    double delta_1 = residu(u_1, solution_1, residu_1, N, 1);
    double delta_2 = residu(u_2, solution_2, residu_2, N * k, 1);
    double delta_3 = residu(u_3, solution_3, residu_3, N * k * k, 1);

    write("Residu_FDM_1.txt", x_1, residu_1, N + 1);
    write("Residu_FDM_2.txt", x_2, residu_2, N * k + 1);
    write("Residu_FDM_3.txt", x_3, residu_3, N * k * k + 1);

    delete[] solution_1;
    delete[] solution_2;
    delete[] solution_3;

    ///////////////////////////////////////////////////////
    // задание 5
    ///////////////////////////////////////////////////////

    double* x_1_SM = new double[N + 1], * x_2_SM = new double[N * k + 1], * x_3_SM = new double[N * k * k + 1];
    double* u_1_SM = new double[N + 1], * u_2_SM = new double[N * k + 1], * u_3_SM = new double[N * k * k + 1];

    SM(x_1_SM, u_1_SM, a, b, N, pow((b - a) / N, 4));
    SM(x_2_SM, u_2_SM, a, b, N * k, pow((b - a) / (N * k), 4));
    SM(x_3_SM, u_3_SM, a, b, N * k * k, pow((b - a) / (N * k * k), 4));

    write("SM_1.txt", x_1_SM, u_1_SM, N + 1);
    write("SM_2.txt", x_2_SM, u_2_SM, N* k + 1);
    write("SM_3.txt", x_3_SM, u_3_SM, N* k* k + 1);

    residu(u_1, u_1_SM, residu_1, N, 1);
    residu(u_2, u_2_SM, residu_2, N * k, 1);
    residu(u_3, u_3_SM, residu_3, N * k * k, 1);

    write("Residu_FDM_SM_1.txt", x_1, residu_1, N + 1);
    write("Residu_FDM_SM_2.txt", x_2, residu_2, N * k + 1);
    write("Residu_FDM_SM_3.txt", x_3, residu_3, N * k * k + 1);

    delete[] residu_1;
    delete[] residu_2;
    delete[] residu_3;
    delete[] x_1_SM;
    delete[] x_2_SM;
    delete[] x_3_SM;
    delete[] u_1_SM;
    delete[] u_2_SM;
    delete[] u_3_SM;
    delete[] x_1;
    delete[] x_2;
    delete[] x_3;
    delete[] u_1;
    delete[] u_2;
    delete[] u_3;
}