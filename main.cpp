#include <iostream>
#include <fstream>
#include <math.h>

const double a = 0.5, b = 2, eps = 1.0E-5,eps2= 1.0E-8;
const double c = 0, d = M_PI,dopEps = 1.0E-15;

//ряд Тейлора
double J(double x)
{
    if (x == 0)
        return 0;
    double sum = (-x*x)/4, an = (-x*x)/4;
    int k = 1;
    while (fabs(an) > eps)
    {
        an = an * ( (-x * x  * k) / ( (2 * k + 2)  * (2 * k + 1) * (k + 1)));
        sum += an;
        k++;
    }

    return sum;
}

std::pair<double,int> J1(double x)
{
    if (x == 0)
        return {0,1};
    double sum = (-x*x)/4, an = (-x*x)/4;
    int k = 1;
    while (fabs(an) > eps)
    {
        an = an * ( (-x * x  * k) / ( (2 * k + 2)  * (2 * k + 1) * (k + 1)));
        sum += an;
        k++;
    }
    return std::pair{sum,k};
}
int*  Tab1(double* x, int n) {
    int *k = new int[n + 1];
    for (int i = 0; i <= n; i++)
    {
        k[i] = J1(x[i]).second;
    }
    return k;
}

//табуляция
double* Tab(double* x, int n)
{
    double* y = new double[n + 1];
    for (int i = 0; i <= n; i++)
        y[i] = J(x[i]);
    return y;
}


// Равномерные узлы
double* Ravn(int n)
{
    double* x = new double[n + 1];
    double h = 0.2;//(b - a) / n;
    for (int i = 0; i <= n; i++)
        x[i] = a + i * h;
    return x;
}

// Корни полинома Чебышева
double* Cheb(int n)
{
    double* x = new double[n + 1];
    for (int i = 0; i <= n; i++)
        x[i] = (b - a) / 2 * cos((2 * i + 1) * M_PI / (2 * n + 2)) + (b + a) / 2;
    return x;
}

//Разделенные разности для формулы Ньютона
double* Dif(double* x, int m)
{
    double* f = new double[m + 1];
    for (int i = 0; i <= m; i++)
        f[i] = J(x[i]);
    for (int k = 1; k <= m; k++)
        for (int i = m; i >= k; i--)
            f[i] = (f[i] - f[i - 1]) / (x[i] - x[i - k]);
    return f;
}

//Вычисление сумм
double Ln(double* f, double* x1, double x, int n)
{
    double q = 1, sum = f[0];
    for (int i = 1; i <= n; i++)
    {
        q *= (x - x1[i - 1]);
        sum += q * f[i];
    }
    return sum;
}

//Полином Ньютона, итоговая подстановка
double* polNewton(double* x, double* x1, int n, int m)
{
    double* L = new double[n + 1];
    double* f = Dif(x1, m);
    for (int i = 0; i <= n; i++)
        L[i] = Ln(f, x1, x[i], m);
    delete[] f;
    return L;
}

double* Eps(double* L, double* J, int n)
{
    double* eps = new double[n + 1];
    for (int i = 0; i <= n; i++)
        eps[i] = fabs(L[i] - J[i]);
    return eps;
}

//Функция поиска максимума
double Max(double* x, int n)
{
    double max = x[0];
    for (int i = 1; i <= n; i++)
        if (max < x[i])
            max = x[i];
    return max;
}

//Погрешность интерполяции
double maxEps(double* J, double* x, int n, int m, char c)
{
    double* x1;
    if (c == 'r')
        x1 = Ravn(m);
    else x1 = Cheb(m);
    double* L = polNewton(x, x1, n, m);
    double* eps = Eps(L, J, n);
    double max = Max(eps, n);
    delete[] x1, L, eps;
    return max;
}

void writeMaxEps(double* J, double* x, int n, int m, char c, std::ofstream& f)
{
    double max1 = maxEps(J, x, n, m, c);
    f << m << '\t' << max1 << '\n';
    m += 5;
    double max2 = maxEps(J, x, n, m, c);
    f << m << '\t' << max2 << '\n';
    while (m <= 70 )
    {
        max1 = max2;
        m += 5;
        max2 = maxEps(J, x, n, m, c);
        f << m << '\t' << max2 << '\n';
    }
    f << '\n';
}

void writeTab(double* y, double* x, int n, std::ofstream& f,int* kol)
{
    f <<  "X[i]" << "\t" << "J" <<  "\t" << "K"<<'\n';
    for (int i = 0; i <= n; i++)
        f << x[i] << "\t" << y[i]<< "\t" << kol[i] << '\n';
    f << '\n';
}

void write(double* y, double* x, int n, std::ofstream& f)
{
    f << "X[i]" << "\t" << "J" << '\n';
    for (int i = 0; i <= n; i++)
        f << x[i] << "\t" << y[i] << '\n';
    f << '\n';
}

//3 задание
//квадратурная формула центральных прямоугольников
double SumCenterRect(double x, int N)
{
    double s = 0, h = (x - c) / N;
    for (int i = 0; i <= N; i++)
    {
        double a = (cos(c + i * h - h/2)-1)/(c + i * h - h/2);
        s += a * h;
    }
    return s;
}

//квадратурная формула трапеции
double SumTrap(double x, int N)
{
    double s = 0, h = (x - c) / N;
    for (int i = 2 ;i <= N ; i++)
    {
        double a = (cos(c + (i-1) * h) - 1)/(c + (i-1) * h) +  (cos(c + i * h)-1)/(c + i * h);
        s += a * h/2;
    }
    return s ;
}

//квадратурная формула Симпсона
double SumSim(double x, int N)
{
    double s = 0, h = (x - c) / N;
    for(int i = 2; i <= N; ++i)
    {
        double a = (cos(c + (i-1) * h) - 1)/(c + (i-1) * h) + 4*(cos(c + i * h - h/2)-1)/(c + i * h - h/2) +(cos(c + i * h) - 1)/(c + i * h);
        s += a * h/6;
    }
    return s ;
}
//квадратурная формула Гаусса
double SumGauss(double x, int N)
{
    double s = 0, h = (x - c) / N;
    for (int i = 0; i <= N - 1; i++)
    {
        double a =  (cos(c + i * h + (h / 2) * (1 - 1 / sqrt(3.0)))-1)/(c + i * h + (h / 2) * (1 - 1 / sqrt(3.0))) +  (cos(c + i * h + (h / 2) *
                                                                                                                                       (1 + 1 / sqrt(3.0)))-1)/(c + i * h + (h / 2) *
                                                                                                                                                                            (1 + 1 / sqrt(3.0)));
        s += a * h / 2;
    }
    return s ;
}

//вывод
void kvad(double* J, double* x, int n, char c, std::ofstream& f)
{
    for (int i = 1; i <= n; i++)
    {
        //кол-во разбиения отрезка,начинаем с двух
        int N = 2;
        double Sn, S2n;

        //t-трапеция, c-центральная,g-гаусса
        if (c == 't') {
            Sn = SumTrap(x[i], N);
        }
        else if (c == 'c') {
            Sn = SumCenterRect(x[i], N);
        }
        else if(c == 's'){
            Sn = SumSim(x[i], N);
        }
        else Sn = SumGauss(x[i], N);

        N *= 2;

        if (c == 't') {
            S2n = SumTrap(x[i], N);
        }
        else if (c == 'c') {
            S2n = SumCenterRect(x[i], N);
        }
        else if(c == 's'){
            S2n = SumSim(x[i], N);
        }
        else S2n = SumGauss(x[i], N);

        //условия достижения точности eps
        while (fabs(Sn - S2n) > eps2)
        {
            Sn = S2n; N *= 2;

            if (c == 't') {
                S2n = SumTrap(x[i], N);
            }
            else if (c == 'c') {
                S2n = SumCenterRect(x[i], N);
            }
            else if(c == 's'){
                S2n = SumSim(x[i], N);
            }
            else S2n = SumGauss(x[i], N);
        }
        //пишем значение в файл
        f << x[i] << '\t' << S2n << '\t' << N << '\t' << fabs(J[i] - S2n) << '\n';
    }
    f << '\n';
}

//4 задание
double dJ(double x)
{
    if (x == 0)
        return 0;
    double an = -x / 2, sum = an;
    int k = 1;
    while (fabs(an) > eps)
    {
        an *= (- x * x) / ((2 * k + 2 ) * (2 * k + 1));
        sum += an;
        k++;
    }
    return sum;
}

double* F(int n)
{
    double* x = new double[n + 1];
    x[0] = J(a);
    x[n] = J(b);
    for (int i = 1; i < n; i++)
        x[i] = x[0] + i * (x[n] - x[0]) / n;
    return x;
}

// Метод хорд(обратная функция)
void MethodChord(double* x, int n, double epsilon, std::ofstream& file) {
    double* f = F(n);
    file << "Метод хорд\n";
    file << f[0] << '\t' << 0 << '\t' << 1 << '\n';

    for (int i = 1; i <=n; i++) {

        double z0 = x[i];
        double z1 = z0 - (J(z0) - f[i]) / dJ(z0);
        int k = 1;

        while (fabs(z1 - z0) > epsilon) {

            double z2 = (z0 * (J(z1) - f[i]) - z1 * J(z0)) / ((J(z1) - f[i]) - J(z0));
            z0 = z1;
            z1 = z2;
            k++;
        }
        file <<x[i]<< '\t'<< f[i] << '\t' << z1 << '\t' << k << '\n';
    }
}

// Метод касательных(обратная функция)
void MethodTanget(double* x, int n, double epsilon, std::ofstream& file) {
    double* f = F(n);
    file << "Метод касательных\n";
    file << f[0] << '\t' << 0 << '\t' << 1 << '\n';

    for (int i = 1; i <= n; i++) {

        double z0 = x[i];
        double z1 = z0 - (J(z0) - f[i]) / dJ(z0);
        int k = 1;

        while (fabs(z1 - z0) > epsilon) {
            double z2 = z1 - (J(z1) - f[i]) / (dJ(z1));
            z0 = z1;
            z1 = z2;
            k++;
        }
        file <<x[i] << '\t'<< f[i] << '\t' << z1 << '\t' << k << '\n';
    }
}

int main()
{
    std::ofstream file;
    file.open("table.txt");
    int n = 10, m = 10;

    double* x = Ravn(n);
    double* J = Tab(x, n);

    int *kol=Tab1(x,n);

    file << "Функция J0:\n";
    writeTab(J, x, n, file,kol);

    //узлы
    double* x1 = Ravn(m);
    double* L = polNewton(x, x1, n, m);
    file << "Полином, построенный по равномерно распределенным узлам\n";
    write(L, x, n, file);

    double* EPS = Eps(L, J, n);
    file << "Погрешность приближения полиномом, построенным по равномерно распределенным узлам\n";
    write(EPS, x, n, file);
    file << "Зависимость максимальной погрешности от количества узлов построения полинома при равномерно распределенных узлах\n";
    writeMaxEps(J, x, n, m, 'r', file);
    double* xCheb = Cheb(m);
    delete[] L;
    L = polNewton(x, xCheb, n, m);
    file << "Полином, построенный по узлам Чебышева\n";
    write(L, x, n, file);
    delete[] EPS;
    EPS = Eps(L, J, n);
    file << "Погрешность приближения полиномом, построенным по узлам Чебышева\n";
    write(EPS, x, n, file);
    file << "Зависимость максимальной погрешности от количества узлов построения полинома при узлах Чебышева\n";
    writeMaxEps(J, x, n, m, 'c', file);


    file << "Квадратурная формула центральных  прямоугольников\n";
    kvad(J, x, n, 'c', file);
    file << "Квадратурная формула трапеций\n";
    kvad(J, x, n, 't', file);
    file << "Квадратурная формула Гаусса\n";
    kvad(J, x, n, 'g', file);
    file << "Квадратурная формула Симпсона\n";
    kvad(J, x, n, 's', file);

    file << "\n";
    MethodChord(x, n, eps, file);
    file << "\n";
    MethodChord(x, n, dopEps, file);

    file << "\n";
    MethodTanget(x, n, eps, file);
    file << "\n";
    MethodTanget(x, n, dopEps, file);

    file.close();
    delete[] x;
    return 0;
}


