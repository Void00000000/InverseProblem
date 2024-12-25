using System.Collections.ObjectModel;
using System.Windows;

namespace InverseProblem
{
    struct PointXZ
    {
        public double X { get; set; }
        public double Z { get; set; }
        public PointXZ(double x, double z)
        {
            X = x;
            Z = z;
        }
        public double Distance(PointXZ p) {
            return Math.Sqrt((p.X - X) * (p.X - X) + (p.Z - Z) * (p.Z - Z));
        }
        public static PointXZ operator +(PointXZ a, PointXZ b) => new PointXZ(a.X + b.X, a.Z + b.Z);
        public static PointXZ operator -(PointXZ a, PointXZ b) => new PointXZ(a.X - b.X, a.Z - b.Z);
        public static PointXZ operator *(PointXZ a, double alpha) => new PointXZ(a.X * alpha, a.Z * alpha);

    }
    abstract class RectangleXZ
    {
        public double X0 { get; init; }
        public double Xn { get; init; }
        public double Z0 { get; init; }
        public double Zn { get; init; }

        public RectangleXZ(double x0, double xn, double z0, double zn)
        {
            X0 = x0;
            Xn = xn;
            Z0 = z0;
            Zn = zn;
        }
    }

    class CellXZ
    {
        public double X0 { get; }
        public double Xn { get; }
        public double Z0 { get; }
        public double Zn { get; }
        public PointXZ P { get; set; }
        public CellXZ(double x0, double xn, double z0, double zn)
        {
            X0 = x0;
            Xn = xn;
            Z0 = z0;
            Zn = zn;
            P = new PointXZ(0, 0);
        }
        public CellXZ(double x0, double xn, double z0, double zn, PointXZ p)
        {
            X0 = x0;
            Xn = xn;
            Z0 = z0;
            Zn = zn;
            P = p;
        }
    }

    //// Ячейка
    //class CellXZ : RectangleXZ
    //{
    //    public PointXZ P { get; set; }

    //    public CellXZ(double x0, double xn, double z0, double zn, PointXZ p) : base(x0, xn, z0, zn) { P = p; }
    //    public CellXZ(double x0, double xn, double z0, double zn) : base(x0, xn, z0, zn) { P = new PointXZ(0, 0); }
    //}

    // Вся система координат--------------------------------------
    class AreaXZ : RectangleXZ
    {
        public AreaXZ(double x0, double xn, double z0, double zn) : base(x0, xn, z0, zn) { }
    }

    // Сетка из ячеек--------------------------------------
    class GridXZ : RectangleXZ
    {
        private List<double> x;
        private List<double> z;
        public ObservableCollection<CellXZ> Cells { get; set; }
        public int Nx { get { return x.Count; } }
        public int Nz { get { return z.Count; } }
        public int K { get { return Cells.Count; } }
        public ReadOnlyCollection<double> X { get { return x.AsReadOnly(); } }
        public ReadOnlyCollection<double> Z { get { return z.AsReadOnly(); } }

        public GridXZ(double x0, double xn, double z0, double zn, int nx, int nz, List<PointXZ> P) : 
            base(x0, xn, z0, zn)
        {
            x = new List<double>();
            z = new List<double>();
            MakeGrid(x, x0, xn, nx);
            MakeGrid(z, z0, zn, nz);
            Cells = new ObservableCollection<CellXZ>();
            for (int j = 0; j < Nz - 1; j++)
            {
                for (int i = 0; i < Nx - 1; i++)
                {
                    Cells.Add(new CellXZ(X[i], X[i + 1], Z[j], Z[j + 1], P[j * (Nx - 1) + i]));
                }
            }
        }

        public GridXZ(double x0, double xn, double z0, double zn, int nx, int nz) :
            base(x0, xn, z0, zn)
        {
            x = new List<double>();
            z = new List<double>();
            MakeGrid(x, x0, xn, nx);
            MakeGrid(z, z0, zn, nz);
            Cells = new ObservableCollection<CellXZ>();
            for (int j = 0; j < Nz - 1; j++)
            {
                for (int i = 0; i < Nx - 1; i++)
                {
                    Cells.Add(new CellXZ(X[i], X[i + 1], Z[j], Z[j + 1]));
                }
            }
        }

        private void MakeGrid(List<double> XZ, double xz0, double xzn, int n)
        {
            XZ.Capacity = n + 1;
            double h = (xzn - xz0) / n;
            for (int i = 0; i < XZ.Capacity; i++) {
                XZ.Add(xz0 + i * h);
            }
        }

        // Нумерация снизу вверх и слева направо
        public void SetP(List<PointXZ> P)
        {
            for (int i = 0; i < P.Count; i++) {
                Cells[i].P = P[i];
            }
        }

        private int global_node(int i, int j)
        {
            return j * (Nx - 1) + i;
        }

        // Возвращает номера соседних с k-о ячейкой
        public List<int> GetNearbyCells(int k) {
            int i = k % (Nx - 1);
            int j = k / (Nx - 1);
            List<int> nearbyNodes = new List<int>();
            bool left = false, right = false, top = false, bottom = false;
            if (i + 1 < Nx - 1)
            {
                nearbyNodes.Add(global_node(i + 1, j));
                right = true;
            }
            if (i - 1 >= 0)
            {
                nearbyNodes.Add(global_node(i - 1, j));
                left = true;
            }
            if (j + 1 < Nz - 1)
            {
                nearbyNodes.Add(global_node(i, j + 1));
                top = true;
            }
            if (j - 1 >= 0)
            {
                nearbyNodes.Add(global_node(i, j - 1));
                bottom = true;
            }
            if (left && bottom)
            {
                nearbyNodes.Add(global_node(i - 1, j - 1));
            }
            if (right && bottom)
            {
                nearbyNodes.Add(global_node(i + 1, j - 1));
            }
            if (left && top)
            {
                nearbyNodes.Add(global_node(i - 1, j + 1));
            }
            if (right && top)
            {
                nearbyNodes.Add(global_node(i + 1, j + 1));
            }
            return nearbyNodes;
        }
    }
}
