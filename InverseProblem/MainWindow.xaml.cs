using OpenTK.Audio.OpenAL.Extensions.Creative.EFX;
using OpenTK.Mathematics;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Media;
using System.Windows.Media.Media3D;
using System.Windows.Shapes;

namespace InverseProblem
{
    public partial class MainWindow : Window
    {
        AreaXZ areaXZ;
        GridXZ gridXZ;
        // Экранные границы окна с выводом расчётной области
        double bottomCanvasScreen;
        double topCanvasScreen;
        double leftCanvasScreen;
        double rightCanvasScreen;
        // Экранные границы координатной системы
        double bottomAreaScreen;
        double topAreaScreen;
        double leftAreaScreen;
        double rightAreaScreen;
        // Отсупы в процентах
        double bottomOffset = 0.1;
        double topOffset = 0.1;
        double leftOffset = 0.1;
        double rightOffset = 0.2;
        // Параметры
        Profile profile; // Профиль вдоль которого расположены установлены приёмники
        // Векторы магнитной индукции (решение прямой задачи с исходным P)
        List<double> bx;
        List<double> bz;
        // Векторы магнитной индукции (решение обратной задачи с вычисленным P)
        List<double> bx_inverse;
        List<double> bz_inverse;
        // Вектор намагниченности (истинный)
        List<PointXZ> P_true = new List<PointXZ>();
        // Векторы намагниченности (подсчитанный в результате прямой задачи)
        List<PointXZ> P = new List<PointXZ>();
        // Максимальное и минимальное значение Px или Py, нужны для цветовой шкалы
        double pxMin, pxMax;
        double pzMin, pzMax;
        // Невязка
        double PHI;
        // Для гамма регуляризации
        List<double> gamma;
        // Реальные границы координатной системы
        public double LeftAreaReal { get; set; }
        public double RightAreaReal { get; set; }
        public double BottomAreaReal { get; set; }
        public double TopAreaReal { get; set; }
        // Реальные границы сетки
        public double LeftGridReal { get; set; }
        public double RightGridReal { get; set; }
        public double BottomGridReal { get; set; }
        public double TopGridReal { get; set; }
        public int Nx { get; set; } // Количество ячеек по оси X
        public int Nz { get; set; } // Количество ячеек по оси Z
        // Для какой координаты строить P
        // true: X, false: Z
        public bool ShowX { get; set; } = true;
        public double I { get; set; } // Сила тока в амперах
        public double ProfileX0 { get; set; } // Левая граница профиля
        public double ProfileXn { get; set; } // Правая граница профиля
        public double ProfileZ { get; set; } // Высота профиля
        public int N { get; set; } // Количество приёмников
        public double Alpha { get; set; } // Используется для альфа регуляризации
        public double diffP { get; set; } // Используется для гамма регуляризации

        public MainWindow()
        {
            InitializeComponent();
            FunctionPlot.Plot.XLabel("X, Метр");
            FunctionPlot.Plot.YLabel("B, Тесла");
            FunctionPlot.Plot.Legend.FontSize = 18;
            InitDefaultValues();
            DataContext = this;
            ////Тест 16 ячеек---------------------------------- -
            //LeftGridReal = -200;
            //RightGridReal = 200;
            //BottomGridReal = -200;
            //TopGridReal = 0;
            //Nx = 4; Nz = 4;
            //PointXZ[] P_copy = new PointXZ[4];
            //P_true.CopyTo(P_copy);
            //P_true = new List<PointXZ>();
            //for (int i = 0; i < Nx * Nz; i++)
            //{
            //    if (i == 5)
            //        P_true.Add(P_copy[0]);
            //    else if (i == 6)
            //        P_true.Add(P_copy[1]);
            //    else if (i == 9)
            //        P_true.Add(P_copy[2]);
            //    else if (i == 10)
            //        P_true.Add(P_copy[3]);
            //    else
            //        P_true.Add(new PointXZ(0, 0));
            //}
            ////--------------------------------------------------

            // Тест 64 ячеек-----------------------------------
            LeftAreaReal = -600;
            RightAreaReal = 600;
            BottomAreaReal = -500;
            TopAreaReal = 300;

            LeftGridReal = -400;
            RightGridReal = 400;
            BottomGridReal = -450;
            TopGridReal = -50;
            Nx = 8; Nz = 8;
            PointXZ[] P_copy = new PointXZ[4];
            P_true.CopyTo(P_copy);
            P_true = new List<PointXZ>();
            for (int i = 0; i < Nx * Nz; i++)
            {
                if (i == 0 || i == 1)
                    P_true.Add(new PointXZ(1, 1));
                else
                    P_true.Add(new PointXZ(0, 0));
                //if (i == 27)
                //    P_true.Add(P_copy[0]);
                //else if (i == 28)
                //    P_true.Add(P_copy[1]);
                //else if (i == 35)
                //    P_true.Add(P_copy[2]);
                //else if (i == 36)
                //    P_true.Add(P_copy[3]);
                //else
                //    P_true.Add(new PointXZ(0, 0));
            }
            gamma = new List<double>(2 * Nx * Nz);
            for (int i = 0; i < 2 * Nx * Nz; i++) gamma.Add(0);
            //--------------------------------------------------
           
            CalculateForward();
            CalculateInverse();
        }

        public void InitDefaultValues()
        {
            ShowX = false;
            LeftAreaReal = -300;
            RightAreaReal = 300;
            BottomAreaReal = -400;
            TopAreaReal = 200;

            LeftGridReal = -100;
            RightGridReal = 100;
            BottomGridReal = -150;
            TopGridReal = -50;
            Nx = 2; Nz = 2;

            ProfileX0 = -600;
            ProfileXn = 600;
            ProfileZ = 50;
            N = 500;

            I = 1;
            P_true = new List<PointXZ>(Nx * Nz);
            for (int i = 0; i < Nx * Nz; i++){
                if (i == 0 && i == 1)
                    P_true.Add(new PointXZ(1, 1));
                else
                    P_true.Add(new PointXZ(0, 0));
            }
            Alpha = 10e-18;
            diffP = 100;
        }

        public void CalculateForward()
        {
            areaXZ = new AreaXZ(LeftAreaReal, RightAreaReal, BottomAreaReal, TopAreaReal);
            gridXZ = new GridXZ(LeftGridReal, RightGridReal, BottomGridReal, TopGridReal, Nx, Nz, P_true);
            profile = new Profile(ProfileZ, N, ProfileX0, ProfileXn);
            (bx, bz) = ForwardProblemSolver.Solve(I, profile, gridXZ);
            pxMin = P_true.Min(p => p.X);
            pxMax = P_true.Max(p => p.X);
            pzMin = P_true.Min(p => p.Z);
            pzMax = P_true.Max(p => p.Z);
            PlotGraphics(ShowX ? bx : bz, "Прямая" ,ShowX ? "Bx" : "Bz", ScottPlot.Colors.Blue);
            TextBlockStatusBar.Text = $"Количество ячеейк = {gridXZ.K}, Невязка = 0";

            DataGridCells.ItemsSource = gridXZ.Cells;
        }

        public void CalculateInverse()
        {
            areaXZ = new AreaXZ(LeftAreaReal, RightAreaReal, BottomAreaReal, TopAreaReal);
            gridXZ = new GridXZ(LeftGridReal, RightGridReal, BottomGridReal, TopGridReal, Nx, Nz);
            profile = new Profile(ProfileZ, N, ProfileX0, ProfileXn);
            P = InverseProblemSolver.Solve(I, profile, gridXZ, bx, bz, Alpha, gamma, diffP, 200, 1);
            gridXZ.SetP(P);
            (bx_inverse, bz_inverse) = ForwardProblemSolver.Solve(I, profile, gridXZ);
            PlotGraphics(ShowX ? bx_inverse : bz_inverse, "Обратная", ShowX ? "Bx" : "Bz", ScottPlot.Colors.Red);
            PHI = 0;
            for (int i = 0; i < profile.N; i++)
            {
                PHI += (bx[i] - bx_inverse[i]) * (bx[i] - bx_inverse[i]) + (bz[i] - bz_inverse[i]) * (bz[i] - bz_inverse[i]);
            }
            TextBlockStatusBar.Text = $"Количество ячеейк = {gridXZ.K}, Невязка = {PHI:E5}";

            DataGridCells.ItemsSource = gridXZ.Cells;
        }

        private void gridCanvas_Loaded(object sender, RoutedEventArgs e)
        {
            DrawGridCanvas();
        }

        private void PlotGraphics(List<double> b, string legengText, string title, ScottPlot.Color color)
        {
            var bPlot = FunctionPlot.Plot.Add.Scatter(profile.X.ToList(), b);
            bPlot.LegendText = legengText;
            bPlot.Color = color;
            FunctionPlot.Plot.Title(title);
            FunctionPlot.Refresh();
        }

        private void InitScreenCoordinates()
        {
            bottomCanvasScreen = gridCanvas.ActualHeight;
            topCanvasScreen = 0;
            leftCanvasScreen = 0;
            rightCanvasScreen = gridCanvas.ActualWidth;

            bottomAreaScreen = bottomCanvasScreen - gridCanvas.ActualHeight * bottomOffset;
            topAreaScreen = topCanvasScreen + gridCanvas.ActualHeight * topOffset;
            leftAreaScreen = leftCanvasScreen + gridCanvas.ActualWidth * leftOffset;
            rightAreaScreen = rightCanvasScreen - gridCanvas.ActualWidth * rightOffset;
        }

        private double ToScreenX(double x)
        {
            double res = leftAreaScreen * (areaXZ.Xn - x) / (areaXZ.Xn - areaXZ.X0) +
                       rightAreaScreen * (x - areaXZ.X0) / (areaXZ.Xn - areaXZ.X0);
            return res;
        }

        private double ToScreenZ(double z)
        {
            double res = bottomAreaScreen * (areaXZ.Zn - z) / (areaXZ.Zn - areaXZ.Z0) +
                       topAreaScreen * (z - areaXZ.Z0) / (areaXZ.Zn - areaXZ.Z0);
            return res;
        }

        private void DrawArea() {
            // Оси
            Line xAxis = new Line
            {
                X1 = leftAreaScreen,
                Y1 = bottomAreaScreen,
                X2 = rightAreaScreen,
                Y2 = bottomAreaScreen,
                Stroke = Brushes.Black,
                StrokeThickness = 3
            };

            Line zAxis = new Line
            {
                X1 = leftAreaScreen,
                Y1 = topAreaScreen,
                X2 = leftAreaScreen,
                Y2 = bottomAreaScreen,
                Stroke = Brushes.Black,
                StrokeThickness = 3
            };

            TextBlock xLabel = new TextBlock
            {
                Text = "X",
                Foreground = Brushes.Black
            };

            TextBlock yLabel = new TextBlock
            {
                Text = "Z",
                Foreground = Brushes.Black
            };
            Canvas.SetLeft(xLabel, rightAreaScreen + 3);
            Canvas.SetTop(xLabel, bottomAreaScreen);
            Canvas.SetLeft(yLabel, leftAreaScreen - 12);
            Canvas.SetTop(yLabel, topAreaScreen - 12);


            gridCanvas.Children.Add(xAxis);
            gridCanvas.Children.Add(zAxis);
            gridCanvas.Children.Add(xLabel);
            gridCanvas.Children.Add(yLabel);

            // Сеточка
            for (int i = 0; i < gridXZ.Nx; i++)
            {
                double x = gridXZ.X[i];
                Line xiLine = new Line
                {
                    X1 = ToScreenX(x),
                    X2 = ToScreenX(x),
                    Y1 = topAreaScreen,
                    Y2 = bottomAreaScreen,
                    Stroke = Brushes.Gray,
                    StrokeThickness = 1,
                    StrokeDashArray = new DoubleCollection() { 2, 2 }
                };
                gridCanvas.Children.Add(xiLine);
            }
            for (int i = 0; i < gridXZ.Nz; i++)
            {
                double z = gridXZ.Z[i];
                Line ziLine = new Line
                {
                    X1 = leftAreaScreen,
                    X2 = rightAreaScreen,
                    Y1 = ToScreenZ(z),
                    Y2 = ToScreenZ(z),
                    Stroke = Brushes.Gray,
                    StrokeThickness = 1,
                    StrokeDashArray = new DoubleCollection() { 2, 2 }
                };
                gridCanvas.Children.Add(ziLine);
            }

            Line topAreaScreenLine = new Line
            {
                X1 = leftAreaScreen,
                X2 = rightAreaScreen,
                Y1 = topAreaScreen,
                Y2 = topAreaScreen,
                Stroke = Brushes.Gray,
                StrokeThickness = 1
            };
            Line rightAreaScreenLine = new Line
            {
                X1 = rightAreaScreen,
                X2 = rightAreaScreen,
                Y1 = topAreaScreen,
                Y2 = bottomAreaScreen,
                Stroke = Brushes.Gray,
                StrokeThickness = 1
            };
            gridCanvas.Children.Add(topAreaScreenLine);
            gridCanvas.Children.Add(rightAreaScreenLine);
        }

        private void DrawGrid()
        {
            foreach (CellXZ cell in gridXZ.Cells)
            {
                double pi = ShowX ? cell.P.X : cell.P.Z;
                double minGradValue = ShowX ? pxMin : pzMin;
                double maxGradValue = ShowX ? pxMax : pzMax;
                double MinColorR = 0, MinColorG = 0, MinColorB = 0;
                double MaxColorR = 0, MaxColorG = 1, MaxColorB = 0;
                double r, g, b;

                if (pi < minGradValue) { r = MinColorR; g = MinColorG; b = MinColorB; }
                else if (pi > maxGradValue) { r = MaxColorR; g = MaxColorG; b = MaxColorB; }
                else
                {
                    double h = maxGradValue - minGradValue;
                    r = MinColorR * (maxGradValue - pi) / h + MaxColorR * (pi - minGradValue) / h;
                    g = MinColorG * (maxGradValue - pi) / h + MaxColorG * (pi - minGradValue) / h;
                    b = MinColorB * (maxGradValue - pi) / h + MaxColorB * (pi - minGradValue) / h;
                }
                Rectangle rect = new Rectangle
                {
                    Width = ToScreenX(cell.Xn) - ToScreenX(cell.X0),
                    Height = ToScreenZ(cell.Z0) - ToScreenZ(cell.Zn),
                    Stroke = Brushes.Black,
                    Fill = new SolidColorBrush(Color.FromRgb((byte)(r * 255), (byte)(g * 255), (byte)(b * 255)))
                };
                Canvas.SetLeft(rect, ToScreenX(cell.X0));
                Canvas.SetTop(rect, ToScreenZ(cell.Zn));
                gridCanvas.Children.Add(rect);
                // P of cell label
                TextBlock piLabel = new TextBlock
                {
                    Text = pi.ToString("0.###"),
                    Foreground = Brushes.White,
                };
                Canvas.SetLeft(piLabel, ToScreenX(cell.X0));
                Canvas.SetTop(piLabel, ToScreenZ(cell.Zn));
                gridCanvas.Children.Add(piLabel);
            }
            // x labels
            for (int i = 0; i < gridXZ.Nx; i++)
            {
                if (i % 2 != 0) continue;
                double x = gridXZ.X[i];
                TextBlock xiLabel = new TextBlock
                {
                    Text = x.ToString(),
                    Foreground = Brushes.Black
                };
                Canvas.SetLeft(xiLabel, ToScreenX(x) - 12);
                Canvas.SetTop(xiLabel, bottomAreaScreen);
                gridCanvas.Children.Add(xiLabel);
            }
            // z labels
            for (int i = 0; i < gridXZ.Nz; i++)
            {
                if (i % 2 != 0) continue;
                double z = gridXZ.Z[i];
                TextBlock ziLabel = new TextBlock
                {
                    Text = z.ToString(),
                    Foreground = Brushes.Black
                };
                Canvas.SetLeft(ziLabel, leftAreaScreen - 27);
                Canvas.SetTop(ziLabel, ToScreenZ(z) - 12);
                gridCanvas.Children.Add(ziLabel);
            }
        }

        private void DrawProfile()
        {
            Line profileLine = new Line
            {
                X1 = leftAreaScreen,
                Y1 = ToScreenZ(profile.Z),
                X2 = rightAreaScreen,
                Y2 = ToScreenZ(profile.Z),
                Stroke = Brushes.Red,
                StrokeThickness = 2,
                StrokeDashArray = new DoubleCollection() { 2, 2 }
            };
            gridCanvas.Children.Add(profileLine);
            TextBlock profileLabel = new TextBlock
            {
                Text = profile.Z.ToString(),
                Foreground = Brushes.Red
            };
            Canvas.SetLeft(profileLabel, leftAreaScreen - 27);
            Canvas.SetTop(profileLabel, ToScreenZ(profile.Z) - 12);
            gridCanvas.Children.Add(profileLabel);
        }

        private void DrawColorScale()
        {
            double colorScaleLeft = rightCanvasScreen - gridCanvas.ActualWidth * (leftOffset);
            double colorScaleRight = rightCanvasScreen - gridCanvas.ActualWidth * (leftOffset / 2);
            for (int i = 1; i < 10; i++)
            {
                Rectangle scaleRect = new Rectangle
                {
                    Width = colorScaleRight - colorScaleLeft,
                    Height = gridCanvas.ActualHeight / 11,
                    Fill = new SolidColorBrush(Color.FromRgb(0, (byte)((i - 1) * 255 / 8), 0))
                };
                Canvas.SetLeft(scaleRect, colorScaleLeft);
                Canvas.SetTop(scaleRect, gridCanvas.ActualHeight - (i + 1) * gridCanvas.ActualHeight / 11);
                gridCanvas.Children.Add(scaleRect);

                TextBlock piLabel = new TextBlock
                {
                    Text = ShowX ? ((9 - i) * pxMin / 8 + (i - 1)*pxMax / 8).ToString("0.####") : ((9 - i) * pzMin / 8 + (i - 1) * pzMax / 8).ToString("0.####"),
                    Foreground = Brushes.Black
                };

                Canvas.SetLeft(piLabel, colorScaleLeft - 38);
                Canvas.SetTop(piLabel, gridCanvas.ActualHeight - (i + 1) * gridCanvas.ActualHeight / 11);
                gridCanvas.Children.Add(piLabel);
            }
        }

        private void DrawGridCanvas()
        {
            gridCanvas.Children.Clear();
            InitScreenCoordinates();
            DrawArea();
            DrawGrid();
            DrawProfile();
            DrawColorScale();
        }

        private void gridCanvas_SizeChanged(object sender, SizeChangedEventArgs e)
        {
            if (!gridCanvas.IsLoaded) { return; }
            DrawGridCanvas();
        }
    }
}