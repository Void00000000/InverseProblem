using MathNet.Numerics.LinearAlgebra;
using ScottPlot;
using System.Collections.ObjectModel;
namespace InverseProblem
{
    class Profile
    {
        private List<double> x; // X координаты приёмника
        public ReadOnlyCollection<double> X { get { return x.AsReadOnly(); } }
        public double Z { get;} // Расстояние между поверхностью и приёмниками
        public int N { get; } // Количество приёмников
        public double X0 { get; } // Начало профиля
        public double Xn { get; } // Конец профиля
        public Profile(double z, int n, double x0, double xn)
        {
            Z = z; N = n; X0 = x0; Xn = xn;
            x = new List<double>();
            x.Capacity = n;
            double h = (xn - x0) / (n - 1);
            for (int i = 0; i < x.Capacity; i++)
            {
                x.Add(x0 + i * h);
            }
        }
        public PointXZ GetPoint(int reciver_num)
        {
            return new PointXZ(x[reciver_num], Z);
        }
    }

    static class ForwardProblemSolver
    {
        static public (List<double>, List<double>) Solve(double I, Profile Profile, GridXZ Grid)
        {
            List<double> bx = new List<double>(Profile.N);
            List<double> bz = new List<double>(Profile.N);
            for (int i = 0; i < Profile.N; i++)
            {
                PointXZ m = Profile.GetPoint(i);
                double b_x = 0; double b_z = 0;
                foreach (CellXZ cell in Grid.Cells)
                {
                    PointXZ center = new PointXZ((cell.Xn + cell.X0) / 2, (cell.Zn + cell.Z0) / 2);
                    PointXZ newXZ = m + center;
                    double r = center.Distance(m);
                    double mes = (cell.Xn - cell.X0) * (cell.Zn - cell.Z0);
                    b_x += (mes * I) / (4 * Math.PI * r * r * r) *
                        (cell.P.X * (3 * newXZ.X * newXZ.X / (r * r) - 1) +
                        cell.P.Z * (3 * newXZ.X * newXZ.Z / (r * r)));
                    b_z += (mes * I) / (4 * Math.PI * r * r * r) *
                        (cell.P.X * (3 * newXZ.X * newXZ.Z / (r * r)) +
                        cell.P.Z * (3 * newXZ.Z * newXZ.Z / (r * r) - 1));
                }
                bx.Add(b_x);
                bz.Add(b_z);
            }
            return (bx, bz);
        }
    }

    // Решение обратной задачи
    static class InverseProblemSolver
    {
        // true если p1 отличается от p2 более чем в diffP раз
        static private bool Compare(double p1, double p2, double diffP)
        {
            double a, b;
            if (p1 > p2)
            {
                a = p1; b = p2;
            }
            else
            {
                a = p2; b = p1;
            }
            //if (b == 0 && a != 0) return true;
            if (b == 0) return true;
            if (a / b > diffP)
                return true;
            return false;
        }

        // diffP - Допустимая разница между параметрами в соседних ячейках
        // kGamma - коэффициент, на который увеличивается gamma
        static public List<PointXZ> Solve(double I, Profile Profile, GridXZ Grid, List<double> bx, List<double> bz, double alpha, List<double> gamma, double diffP, double kGamma = 2, int maxiter=200)
        {
            List<PointXZ> P = new List<PointXZ>();
            var S = Vector<double>.Build.Dense(Profile.N * 2);
            List<(int, int)> processed_cells = new List<(int, int)> ();

            for (int i = 0; i < Profile.N; i++)
            {
                S[2*i] = bx[i];
                S[2*i + 1] = bz[i];
            }
            var L = Matrix<double>.Build.Dense(Profile.N * 2, Grid.K * 2);
            var l = Matrix<double>.Build.Dense(2, 2);
            for (int i = 0; i < Profile.N; i++)
            {
                for (int j = 0; j < Grid.K; j++)
                {
                    CellXZ cell = Grid.Cells[j];
                    PointXZ m = Profile.GetPoint(i);
                    PointXZ center = new PointXZ((cell.Xn + cell.X0) / 2, (cell.Zn + cell.Z0) / 2);
                    PointXZ newXZ = m + center;
                    double r = center.Distance(m);
                    double mes = (cell.Xn - cell.X0) * (cell.Zn - cell.Z0);

                    l[0, 0] = 3 * newXZ.X * newXZ.X / (r*r) - 1;
                    l[0, 1] = 3 * newXZ.X * newXZ.Z / (r * r);
                    l[1, 0] = 3 * newXZ.X * newXZ.Z / (r * r);
                    l[1, 1] = 3 * newXZ.Z * newXZ.Z / (r * r) - 1;
                    l *= mes * I / (4 * Math.PI * r*r*r);
                    L[2 * i, 2 * j] = l[0, 0];
                    L[2 * i, 2 * j + 1] = l[0, 1];
                    L[2 * i + 1, 2 * j] = l[1, 0];
                    L[2 * i + 1, 2 * j + 1] = l[1, 1];
                }
            }
            // Создание матриц
            var A = Matrix<double>.Build.Dense(Grid.K * 2, Grid.K * 2);
            A = L.Transpose() * L + Matrix<double>.Build.DenseDiagonal(Grid.K * 2, Grid.K * 2, alpha);

            var b = Vector<double>.Build.Dense(Profile.N * 2);
            b = L.Transpose() * S;

            var C = Matrix<double>.Build.Dense(Grid.K * 2, Grid.K * 2);
            bool stop = false;
            int iter = 1;

            // Цикл подбора gamma --------------------------------
            while (!stop && iter <= maxiter)
            {
                // Создание матрицы С ------------------------------------
                for (int k = 0; k < Grid.K; k++)
                {
                    List<int> nearbyCells = Grid.GetNearbyCells(k);
                    for (int i = 0; i <= 1; i++)
                    {
                        int ck = 2 * k + i;
                        double diag = 0;
                        foreach (int m in nearbyCells)
                        {
                            int cm = 2 * m + i;
                            C[ck, cm] = -(gamma[ck] + gamma[cm]);

                            diag += 1 * gamma[ck];
                            diag += gamma[cm];
                        }
                        C[ck, ck] = diag;
                    }
                }
                // Конец создания матрицы C -----------------------------

                // Решение СЛАУ-----------------------
                A += C;
                var q = A.Solve(b);
                P = new List<PointXZ>(Grid.K);
                for (int i = 0; i < Grid.K; i++)
                {
                    P.Add(new PointXZ(q[2 * i], q[2 * i + 1]));
                }
                // Конец решения СЛАУ --------------------------------

                // Проверка ячеек---------------------------------
                stop = true;
                for (int k = 0; k < Grid.K; k++)
                {
                    List<int> nearbyCells = Grid.GetNearbyCells(k);
                    foreach (int m in nearbyCells)
                    {
                        if (processed_cells.Contains((k, m))) continue;
                        if (Compare(P[k].X, P[m].X, diffP)) {
                            stop = false;
                            gamma[2*m] *= kGamma;
                        }
                        if (Compare(P[k].Z, P[m].Z, diffP))
                        {
                            stop = false;
                            gamma[2 * m + 1] *= kGamma;
                        }
                        processed_cells.Add((k, m));
                        if (m != k)
                            processed_cells.Add((m, k));
                    }
                }
                // конец проверки ячеек -------------------------------------------
                processed_cells.Clear();
                A -= C;
                iter++;
            }
            // конец цикла --------------------------------------
            return P;
        }
    }
}
