#define _USE_MATH_DEFINES
#include <format>
#include <iostream>
#include <vector>
#include <math.h>


class Solve
{
public:
    Solve(int n, int m) 
    {
        this->n = n, this->m = m; 

        h = (b - a) / (double)n;
        k = (d - c) / (double)m;

        h2 = 1. / (h * h);
        k2 = 1. / (k * k);

        A = -2. * (h2 + k2);
    }
    ~Solve() {}
    void mainTask(int per = 1, int iterMax = 100, double eps = 1e-8)
    {
        std::vector<double> resMain1 = chebMain(per, iterMax, eps, MAIN);
        std::cout << std::format("Для решения основной задачи использованы сетка с числом разбиений \n\
по x n = «{}» и числом разбиений по y m = «{}»,\n\
метод с чебышевским набором  параметров,\n\
параметры K = «{}»\n\
критерии остановки по точности epsмет = «{}» и по числу итераций Nmax = «{}»\n\
На решение схемы(СЛАУ) затрачено итераций N = «{}» и достигнута точность итерационного метода eps(N) = «{}»\n\
\
Схема(СЛАУ) решена с невязкой || R(N) || = «{}»\n\
для невязки СЛАУ использована норма «Евклидова»;\n\
В качестве начального приближения использовано «Нулевое»\n\
На основной сетке невязка СЛАУ на начальном приближении || R(0) || = «{}» \n", n, m, per, eps, iterMax, iter, eps_cur, Rlast, R0);


        n *= 2; m *= 2;
        std::vector<double> resMain2 = chebMain(per, iterMax, eps, MAIN);
        n /= 2; m /= 2;

        maxDiff = 0;

        for (int i = 0; i < resMain1.size(); i++)
        {
            int xInd = i % (n - 1);
            int yInd = i / (n - 1);

            int xIndRef = xInd * 2 + 1;
            int yIndRef = yInd * 2 + 1;
            //std::cout << xInd << ' ' << yInd << " -> " << xIndRef << ' ' << yIndRef << '\n';

            if (abs(resMain2[xIndRef + yIndRef*(2*n-1)] - resMain1[i]) > maxDiff)
            {
                maxDiff = abs(resMain2[xIndRef + yIndRef * (2 * n - 1)] - resMain1[i]);
                x = a + 2. * h * (xInd + 1);
                y = c + 2. * k * (yInd + 1);
            }
        }

        std::cout << std::format("Для контроля точности использована сетка с половинным шагом\n\
метод с чебышевским набором  параметров,\n\
параметры K = «{}»\n\
критерии остановки по точности epsмет2 = «{}» и по числу итераций Nmax2 = «{}»\n\
На решение схемы(СЛАУ) затрачено итераций N2 = «{}» и достигнута точность итерационного метода eps(N2) = «{}»\n\
\
Схема(СЛАУ) решена с невязкой || R(N2) || = «{}»\n\
для невязки СЛАУ использована норма «Евклидова»;\n\
Основная задача должна быть решена с точностью не хуже eps = 0.5*10 –6;\n\
задача решена с погрешностью eps2 = «{}»\n\
для погрешности использована норма «Бесконечность»;\n\
Максимальное отклонение точного и численного решений наблюдается в узле x = «{}»; y = «{}»\n\
В качестве начального приближения использовано «Нулевое»\
Невязка СЛАУ на начальном приближении || R(0) || = «{}» \n", per, eps, iterMax, iter, eps_cur, Rlast, maxDiff, x, y, R0);

    }
    void testTask(int per = 1, int iterMax = 100, double eps = 1e-10)
    {
        h = (b - a) / (double)n;
        k = (d - c) / (double)m;

        h2 = 1. / (h * h);
        k2 = 1. / (k * k);

        A = -2. * (h2 + k2);

        std::vector<double> U = createTestU();

        std::vector<double> V = chebMain(per, iterMax, eps, TEST);

        std::vector<double> Error((n - 1) * (m - 1));
        maxDiff = 0;

        for (int i = 0; i < (n - 1) * (m - 1); i++)
        {
            Error[i] = V[i] - U[i];
            if (abs(Error[i]) > maxDiff)
            {
                maxDiff = abs(Error[i]);
                x = a + h * (i % (n - 1) + 1);
                y = c + k * (i / (n - 1) + 1);
            }
        }

        std::cout << std::format("Для решения тестовой задачи использованы сетка с числом разбиений \n\
по x n = «{}» и числом разбиений по y m = «{}»,\n\
метод с чебышевским набором  параметров,\n\
параметры K = «{}»\n\
критерии остановки по точности epsмет = «{}» и по числу итераций Nmax = «{}»\n\
На решение схемы(СЛАУ) затрачено итераций N = «{}» и достигнута точность итерационного метода eps(N) = «{}»\n\
\
Схема(СЛАУ) решена с невязкой || R(N) || = «{}»\n\
для невязки СЛАУ использована норма «Евклидова»;\n\
Тестовая задача должна быть решена с погрешностью не более eps = 0.5*10 –6;\n\
задача решена с погрешностью eps1 = «{}»\n\
для погрешности использована норма «Бесконечность»;\n\
Максимальное отклонение точного и численного решений наблюдается в узле x = «{}»; y = «{}»\n\
В качестве начального приближения использовано «Нулевое»\
Невязка СЛАУ на начальном приближении || R(0) || = «{}» \n", n, m, per, eps, iterMax, iter, eps_cur, Rlast, maxDiff, x, y, R0);
    }

    double minLamda()
    {
        return  4. / (h * h) * sin(M_PI / (2. * n)) * sin(M_PI / (2. * n)) +
            4. / (k * k) * sin(M_PI / (2. * m)) * sin(M_PI / (2. * m));
    }
    double maxLamda()
    {
        return  4. / (h * h) * sin(M_PI * (n - 1.) / (2. * n)) * sin(M_PI * (n - 1.) / (2. * n)) +
            4. / (k * k) * sin(M_PI * (m - 1.) / (2. * m)) * sin(M_PI * (m - 1.) / (2. * m));
    }

private:
    int iter;
    double eps_cur, Rlast, maxDiff, x, y, R0;

    double a = 0, b = 1;
    double c = 0, d = 2;
    int n = 100, m = 100;
    double h, k;
    double h2, k2;
    double A;

    enum TASK
    {
        MAIN, TEST
    };
    std::vector<double> chebMain(int per = 1, int iterMax = 100, double eps = 1e-10, TASK task = TEST)
    {
        eps_cur = 10;

        h = (b - a) / (double)n;
        k = (d - c) / (double)m;

        h2 = 1. / (h * h);
        k2 = 1. / (k * k);

        A = -2. * (h2 + k2);

        std::vector<double> tau = createTau(per);

        std::vector<double> V((n - 1) * (m - 1));
        std::vector<double> R((n - 1) * (m - 1));
        std::vector<double>  F;
        if (task == TEST)
        {
            std::cout << "CREATE TEST";
            F = createTestRightPart();
        }
        else
            F = createMainRightPart();



        for (int i = 0; i < (n - 1) * (m - 1); i++)
            V[i] = 0;

        R = calcR(V, F);
        R0 = norm2(R);


        iter = 0;
        for (; iter < iterMax && eps_cur >= eps; iter++)
        {
            if (iterMax <= 1000 || iter % (iterMax / 100) == 0)
                std::cout << iter << "; eps = " << eps_cur << '\n';
            eps_cur = 0;

            R = calcR(V, F);
            //std::cout << normMax(R) << '\n';

            for (int i = 0; i < (n - 1) * (m - 1); i++)
            {
                V[i] = V[i] + tau[iter % per] * R[i];

                if (abs(tau[iter % per] * R[i]) > eps_cur)
                    eps_cur = abs(tau[iter % per] * R[i]);

            }
        }
        R = calcR(V, F);
        Rlast = norm2(R);
        
        return V;
    }


    double uxx(double x, double y)
    {
        return M_PI * M_PI * y * y * exp(pow(sin(M_PI * x * y), 2.)) * pow(sin(2. * M_PI * x * y), 2.) + 2 * M_PI * M_PI * y * y * exp(pow(sin(M_PI * x * y), 2.)) * cos(2. * M_PI * x * y);
    }
    double uyy(double x, double y)
    {
        return M_PI * M_PI * x * x * exp(pow(sin(M_PI * x * y), 2.)) * pow(sin(2. * M_PI * x * y), 2.) + 2 * M_PI * M_PI * x * x * exp(pow(sin(M_PI * x * y), 2.)) * cos(2. * M_PI * x * y);
    }

    double mainF(double x, double y)
    {
        return -abs(x - y);
    }
    double mainNu(double x, double y)
    {
        if (x == 0.)
            return sin(M_PI * y) * sin(M_PI * y);
        if (x == 1.)
            return abs(exp(sin(M_PI * y)) - 1);
        if (y == 0.)
            return x * (1. - x);
        if (y == 2.)
            return x * (1. - x) * exp(x);

        return 0;
    }

    double testU(double x, double y)
    {
        return exp(pow(sin(M_PI * x * y), 2.));
    }
    double testF(double x, double y)
    {
        return -uxx(x, y) - uyy(x, y);
    }
    double testNu(double x, double y)
    {
        return exp(pow(sin(M_PI * x * y), 2.));
    }



    std::vector<double> createTau(int K)
    {
        std::vector<double> res(K);
        double minL = minLamda();
        double maxL = maxLamda();
        for (int i = 0; i < K; i++)
        {
            res[i] = 1. / ((minL + maxL) / 2. + (maxL - minL) / 2. * cos(M_PI / (2. * K) * (1. + 2. * i)));
        }
        return res;
    }
    std::vector<double> createTestU()
    {
        std::vector<double> res((n - 1) * (m - 1));
        for (int i = 0; i < (n - 1) * (m - 1); i++)
        {
            int xInd = i % (n - 1) + 1;
            int yInd = i / (n - 1) + 1;

            double x = a + h * xInd;
            double y = c + k * yInd;

            res[i] = testU(x, y);
        }
        return res;
    }

    std::vector<double> createTestRightPart()
    {
        std::vector<double> res((n - 1) * (m - 1));
        for (int i = 0; i < (n - 1) * (m - 1); i++)
        {
            int xInd = i % (n - 1) + 1;
            int yInd = i / (n - 1) + 1;

            double x = a + h * xInd;
            double y = c + k * yInd;

            res[i] = -testF(x, y);

            if (xInd - 1 == 0)
            {
                res[i] -= testNu(x - h, y) * h2;
            }
            if (xInd + 1 == n)
            {
                res[i] -= testNu(x + h, y) * h2;
            }
            if (yInd - 1 == 0)
            {
                res[i] -= testNu(x, y - k) * k2;
            }
            if (yInd + 1 == m)
            {
                res[i] -= testNu(x, y + k) * k2;
            }
        }
        return res;
    }
    std::vector<double> createMainRightPart()
    {
        std::vector<double> res((n - 1) * (m - 1));
        for (int i = 0; i < (n - 1) * (m - 1); i++)
        {
            int xInd = i % (n - 1) + 1;
            int yInd = i / (n - 1) + 1;

            double x = a + h * xInd;
            double y = c + k * yInd;

            res[i] = -mainF(x, y);

            if (xInd - 1 == 0)
            {
                res[i] -= mainNu(x - h, y) * h2;
            }
            if (xInd + 1 == n)
            {
                res[i] -= mainNu(x + h, y) * h2;
            }
            if (yInd - 1 == 0)
            {
                res[i] -= mainNu(x, y - k) * k2;
            }
            if (yInd + 1 == m)
            {
                res[i] -= mainNu(x, y + k) * k2;
            }
        }
        return res;
    }

    std::vector<double> calcR(std::vector<double> V, std::vector<double> F)
    {
        std::vector<double> res((n - 1) * (m - 1));
        for (int i = 0; i < (n - 1) * (m - 1); i++)
        {
            int xInd = i % (n - 1) + 1;
            int yInd = i / (n - 1) + 1;


            res[i] = -F[i] + A * V[i];

            if (xInd - 1 != 0)
                res[i] += V[i - 1] * h2;
            if (xInd + 1 != n)
                res[i] += V[i + 1] * h2;
            if (yInd - 1 != 0)
                res[i] += V[i - (n - 1)] * k2;
            if (yInd + 1 != m)
                res[i] += V[i + (n - 1)] * k2;
        }
        return res;
    }
    double norm2(std::vector<double> v)
    {
        double res = 0;
        for (size_t i = 0; i < v.size(); i++)
        {
            res += v[i] * v[i];
        }
        return sqrt(res);
    }
    double normMax(std::vector<double> v)
    {
        double res = 0;
        for (size_t i = 0; i < v.size(); i++)
        {
            if (abs(v[i]) > res)
                res = abs(v[i]);
        }
        return res;
    }
    void printVec(std::vector<double> v)
    {
        for (size_t i = 0; i < v.size(); i++)
        {
            std::cout << v[i] << ' ';
        }
        std::cout << "\n";
    }
};



int main()
{
    setlocale(LC_ALL, "Russian");
    Solve s(4, 4);
    //s.mainTask(20, 500000, 1e-10);
    s.testTask(20, 100000, 1e-10);
}
