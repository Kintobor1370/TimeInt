#define NOMINMAX

#include <iostream>
#include <string>
#include <cmath>
#include <fstream>
#include <sstream>
#include <cassert>
#include <time.h>
#include <windows.h>
#include <conio.h>

#include "types.h"

using namespace std;

#define KEY_ENT 13
#define KEY_ESC 27

const CalcType PI = 3.1415927;
const CalcType EPSILON = 1E-08;

enum KernelType { Normal, };

enum ExtrapolationType { Zero, SquaredN, Alpha, Exponential, };


/**************************************************************************************************************************
**                                            LOGISTIC EQUATION IMPLEMENTATION                                           **
***************************************************************************************************************************
**       The code below realizes the logistic equation for the first and second moments on plane (x, y), as well as      **
**                 properties of the spatial model, kernels generation and all the means of interaction                  **
**                 between parameters and results of the equation and further components of the program.                 **
**************************************************************************************************************************/

//___________________________________________________MODEL PROPERTIES______________________________________________________
struct Properties
{
    int SpeciesNumber;                                                  // Number of species                
        
    Vector b;                                                           // Birth rate oer unit time in species i  (b)
    Vector d;                                                           // Death rate per unit time in species i  (d)
    Matrix d_prime;                                                     // Linear dependence of death rate of species i on neighbours of species  (d')
        
    VectorOfMatrix M;                                                   // Rate of birth-movement of species i to distance (x,y)  (m(x))
    MatrixOfMatrix W;                                                   // Death-weighting for neighbours j of species i at distance (x,y)  (w(x))

    Matrix MaxNumForW;                                                  // Maximum number of non-zero bins for 'w' kernel
    Vector MaxNumForM;                                                  // Maximum number of non-zero bins for 'm' kernel

    Vector SigmaM;
    Matrix SigmaW;

    Vector MomentsGrid;                                                 // The non-uniform grid that corresponds to pair densities
    CalcType CutoffDistance;                                            // The largest distance between two individuals for which we compute the pair density
    Vector KernelsGrid;                                                 // The uniform grid that corresponds to kernels
    CalcType KernelsGridStep;                                           // The step of uniform grid for kernels

    //..............................GENERATE GRIDS FOR KERNELS AND PAIR DENSITIES
    void GenerateGrids(CalcType sigma_w, CalcType sigma_m)
    {
        CalcType cutoff_distance = std::max(sigma_w, sigma_m) * 3;
        CalcType max_radius = cutoff_distance;

        int num_for_moments = 20;
        auto moments_grid = Vector::setLinSpaced(num_for_moments, 0, cutoff_distance);
        
        int num_for_kernels = (num_for_moments - 1) * 2 - 1;
        if (std::min(sigma_w, sigma_m) < 0.015)
        {
            num_for_kernels = (num_for_kernels - 2) * 4 - 1;
            num_for_kernels = (num_for_kernels - 1) * 2 - 1;
        }
        cout << "number of mid-points for pair density grid: " << num_for_moments
             << "\nnumber of mid-points for kernels grid:      " << num_for_kernels << "\n\n";

        assert ((num_for_kernels % 2 == 1) && (num_for_kernels >= 5));
        auto kernels_grid = Vector::setLinSpaced(num_for_kernels, 0, max_radius);

        CutoffDistance  = moments_grid.last();
        MomentsGrid = moments_grid;
        KernelsGrid = kernels_grid;
        KernelsGridStep = kernels_grid[2] - kernels_grid[0];
    }

    //..............................GENERATE KERNEL
    pair<Matrix, CalcType> GenerateKernel(
        CalcType s,                                                     // Sigma 
        CalcType R,                                                     // Max radius 
        KernelType type
    ) {
        Vector grid = KernelsGrid;
        CalcType step = KernelsGridStep;

        Matrix kernel_values(grid.size());
        int max_non_zero_num = 0; 

        // Calculate
        for (int x=0; x<grid.size(); x++)
        {
            for (int y=0; y<grid.size(); y++)
            {
                const auto r = std::hypot(grid[x], grid[y]);
                if(r <= R)
                {
                    switch (type)
                    {
                        case Normal:
                            auto norm = s * std::sqrt(2 * PI) * (1 - std::exp((-R * R) / (2 * s * s)));
                            auto weight = std::exp((-r * r) / (2 * s * s));
                            kernel_values[x][y] = weight / norm;
                            break;
                    }
                    if (kernel_values[x][y] != 0)
                        max_non_zero_num = std::max({max_non_zero_num, x+1, y+1});
                }
                else kernel_values[x][y] = 0;
            }
        }
        if (grid.size() % 2 == 1 && max_non_zero_num % 2 == 0 && max_non_zero_num != 2) max_non_zero_num++;
        if (max_non_zero_num == 1) max_non_zero_num = 2;

        // Renormalize
        CalcType integral = 0;
        for (int x=1-max_non_zero_num; x<max_non_zero_num-1; x+=2)
            for (int y=1-max_non_zero_num; y<max_non_zero_num-1; y+=2)
                integral += kernel_values[abs(x + 1)][abs(y + 1)];
        
        assert (integral != 0);
        integral *= step * step;
        kernel_values = kernel_values / integral;

        return make_pair(kernel_values, max_non_zero_num);
    }

    void CheckKernels()
    {
        cout << "    Checking kernels...\n";
        for (int i=0; i<SpeciesNumber; i++)
        {
            for (int j=0; j<SpeciesNumber; j++)
            {
                cout << "      W_" << i + 1 << "_" << j + 1 << ": integral = ";
                
                CalcType integral = 0;
                for (int x=1-MaxNumForW[i][j]; x<MaxNumForW[i][j]-1; x+=2)
                    for (int y=1-MaxNumForW[i][j]; y<MaxNumForW[i][j]-1; y+=2)
                        integral += W(i, j)[abs(x + 1)][abs(y + 1)];
                
                integral *= std::pow(KernelsGridStep, 2);
                cout << integral;
                assert (std::abs(integral - 1) < EPSILON);
                cout << ".  Cool!!\n";
            }
        }

        for (int i=0; i<SpeciesNumber; i++)
        {
            cout << "      M_" << i + 1 << ":   integral = ";
                
            CalcType integral = 0;
            for (int x=1-MaxNumForM[i]; x<MaxNumForM[i]-1; x+=2)
                for (int y=1-MaxNumForM[i]; y<MaxNumForM[i]-1; y+=2)
                    integral += M(i)[abs(x + 1)][abs(y + 1)];
                
            integral *= std::pow(KernelsGridStep, 2);
            cout << integral;
            assert (std::abs(integral - 1) < EPSILON);
            cout << ".  Cool!!\n";
        }
        cout << "    Check complete!\n\n";
    }

    //..............................GENERATE MODEL PROPERTIES
    static Properties Generate(
        int num,                                                            // Number of species
        Matrix sW,                                                          // Sigma W
        Vector sM,                                                          // Sigma M
        Vector b,                                                           
        Vector d,
        Matrix dd,
        KernelType w_type = Normal,                                         // W kernel type
        KernelType m_type = Normal                                          // M kernel type
    ) {
        Properties props;
        props.SpeciesNumber = num;
        props.SigmaM = sM;
        props.SigmaW = sW;
        props.b = b;
        props.d = d;
        props.d_prime = dd;
        props.GenerateGrids(sW.max_value(), sM.max_value());

        MatrixOfMatrix w(num);
        Matrix max_num_w(num);
        
        Matrix max_radius_W;
        if (w_type == Normal)
            max_radius_W = sW * 3;
        
        for (int i=0; i<num; i++)
        {
            for (int j=0; j<num; j++)
            {
                auto kernel_W = props.GenerateKernel(sW[i][j], max_radius_W[i][j], w_type);
                w(i, j) = kernel_W.first;
                max_num_w[i][j] = kernel_W.second;
            }
        }    

        VectorOfMatrix m(num);
        Vector max_num_m(num);
        
        Vector max_radius_M;
        if (m_type == Normal)
            max_radius_M = sM * 3;

        for (int i=0; i<num; i++)
        {
            auto kernel_M = props.GenerateKernel(sM[i], max_radius_M[i], m_type);
            m(i) = kernel_M.first;
            max_num_m[i] = kernel_M.second;
        }

        props.M = m;
        props.W = w;
        props.MaxNumForW = max_num_w;
        props.MaxNumForM = max_num_m;
        return props;
    }
};


//_______________________________EQUATIONS FOR THE DYNAMICS OF SINGLET AND PAIR DENSITIES__________________________________
class Equation
{
    Properties props;                                                   // Model properties
    Vector N;                                                           // Vector of current singlet densities
    MatrixOfMatrix C;                                                   // Current pair densities
    Vector params;                                                      // Parameters for parametric closures
    ExtrapolationType C_extrapolation_type;

private:
    // Make Vector of N (Vector) and C (Matrix of Matrix) values
    Vector MakeVector(Vector N, MatrixOfMatrix C)
    {
        Vector result = N;
        result.add(MatrixOfMatrixToVector(C));

        assert (result.size() == props.SpeciesNumber + props.SpeciesNumber * props.SpeciesNumber * props.MomentsGrid.size() * props.MomentsGrid.size());
        return result;
    }

    MatrixOfMatrix CalculateIntegralOfWAndC(MatrixOfMatrix C, Vector N, bool square=false)
    {
        MatrixOfMatrix result(props.SpeciesNumber);

        for (int i=0; i<props.SpeciesNumber; i++) 
        {
            for (int j=0; j<props.SpeciesNumber; j++)
            {
                Matrix result_elem(props.MomentsGrid.size());
                for (int x_num=0; x_num<props.MomentsGrid.size(); x_num++)
                {
                    for (int y_num=0; y_num<props.MomentsGrid.size(); y_num++)
                    {
                        result_elem[x_num][y_num] = MakeIntegralWithPairDensity(x_num, y_num, props.MaxNumForW[i][j], props.W(i, j), C(i, j), N[i] * N[j], square);
                    }
                }
                result(i, j) = result_elem;
            }
        }
        return result;
    }

    // Calculte integral m(x)C(x)dx
    MatrixOfMatrix CalculateIntegralOfMAndC(MatrixOfMatrix C, Vector N)
    {
        MatrixOfMatrix result(props.SpeciesNumber);
        
        for (int i=0; i<props.SpeciesNumber; i++)
        {
            for (int j=0; j<props.SpeciesNumber; j++)
            {
                Matrix result_elem(props.MomentsGrid.size());
                for (int x_num=0; x_num<props.MomentsGrid.size(); x_num++)
                {
                    for (int y_num=0; y_num<props.MomentsGrid.size(); y_num++)
                    {
                        result_elem[x_num][y_num] = MakeIntegralWithPairDensity(x_num, y_num, props.MaxNumForM[i], props.M(i), C(i, j), N[i] * N[j]);
                    }
                }
                result(i, j) = result_elem;
            }
        }
        return result;
    }

    CalcType MakeIntegralWithPairDensity(int x_shift, int y_shift, int max_non_zero_num, Matrix kernel, Matrix C, CalcType asymptotic_value, bool square=false)
    {
        CalcType integral = 0;
        for (int i=1-max_non_zero_num; i<max_non_zero_num-1; i+=2)
        {
            for (int j=1-max_non_zero_num; j<max_non_zero_num-1; j+=2)
            {
                assert (i + 1 < max_non_zero_num);
                assert (j + 1 < max_non_zero_num);

                const CalcType sum_x = abs(props.MomentsGrid[x_shift] - (i >= 0 ? props.KernelsGrid[i+1] : -props.KernelsGrid[-i-1]));
                const CalcType sum_y = abs(props.MomentsGrid[y_shift] - (j >= 0 ? props.KernelsGrid[j+1] : -props.KernelsGrid[-j-1]));

                auto pair_density = BilinearInterpolation(props.MomentsGrid, C, sum_x, sum_y, props.CutoffDistance, C_extrapolation_type, asymptotic_value);
                CalcType pair_density2 = 1;
                if (square)
                {
                    assert (props.MomentsGrid[0] == 0);
                    CalcType x1 = abs(i >= 0 ? props.KernelsGrid[i+1] : -props.KernelsGrid[-i-1]);
                    CalcType y1 = abs(j >= 0 ? props.KernelsGrid[j+1] : -props.KernelsGrid[-j-1]);
                    pair_density2 = BilinearInterpolation(props.MomentsGrid, C, x1, y1, props.CutoffDistance, C_extrapolation_type, asymptotic_value);
                }
                integral += kernel[abs(i + 1)][abs(j + 1)] * pair_density * pair_density2;
            }
        }
        integral *= pow(props.KernelsGridStep, 2);
        return integral;
    }

    CalcType BilinearInterpolation(
        Vector grid,
        Matrix values,
        CalcType x,
        CalcType y,
        CalcType max_distance,
        ExtrapolationType type,
        CalcType asymptotic_value = 0
    ) {
        assert (grid.size() == values.size());

        CalcType distance = std::hypot(x, y);
        if (distance > max_distance)
            return Extrapolation(grid, values, distance, max_distance, type, asymptotic_value);
        
        // Find segments [a, b] and [c, d] so that x will be in [a, b] and y will be in [c, d]
        int i = 1, j = 1;
        while (x > grid[i]) i++;
        while (y > grid[j]) j++;
        
        CalcType a = grid[i-1];
        CalcType b = grid[i];
        CalcType c = grid[j-1];
        CalcType d = grid[j];

        // Compute weights
        CalcType wx = (x - a) / (b - a);
        CalcType wy = (y - c) / (d - c);

        // Compute weighted value
        CalcType value = values[i-1][j-1] * (1 - wx) * (1 - wy) + 
                         values[i][j-1] * wx * (1 - wy) +
                         values[i-1][j] * (1 - wx) * wy +
                         values[i][j] * wx * wy;
        return value;
    }

    CalcType Extrapolation(Vector grid, Matrix values, CalcType distance, CalcType max_distance, ExtrapolationType type, CalcType asymptotic_value)
    {
        if (std::abs(values[0].last() - asymptotic_value) < EPSILON)
            return asymptotic_value;
        
        CalcType extrapolated_value;
        switch (type)
        {
            case Zero:
                {
                    extrapolated_value = 0;
                    break;
                }

            case SquaredN:
                {
                    extrapolated_value = asymptotic_value;
                    break;
                }

            case Alpha:
                {
                    CalcType alpha = 2 / (grid.size() + 1);
                    Vector predicted_values;
                    
                    predicted_values.add(values[0][0]);
                    for (int k=1; k<grid.size(); k++)
                    {
                        CalcType new_predicted_value = alpha * values[0][k] + (1 - alpha) * predicted_values[k-1];
                        predicted_values.add(new_predicted_value);
                    }
                    extrapolated_value = alpha * asymptotic_value + (1 - alpha) * predicted_values.last();
                    break;
                }
            
            case Exponential:
                {
                    const int accuracy = 5;
                    CalcType decay_avg = 0;
                    for (int k=1; k<accuracy; k++)
                    {
                        CalcType pair_distance = max_distance - grid[grid.size() - 1 - k];
                        CalcType new_decay_value = std::abs((values[0][grid.size() - 1 - k] - asymptotic_value) / (values[0].last() - asymptotic_value));
                        
                        assert (pair_distance > 0);
                        assert (new_decay_value >= 0);

                        new_decay_value = std::log(new_decay_value) / pair_distance;
                        decay_avg += new_decay_value;
                    }
                    decay_avg = decay_avg / (accuracy - 1);
                    extrapolated_value = asymptotic_value + (values[0].last() - asymptotic_value) * std::exp(-1 * decay_avg * (distance - max_distance));
                }
        }
        return extrapolated_value;
    }

public:
    Equation(Properties my_props, Vector my_params) : props(my_props), params(my_params) {}

    // Make Vector of equilibrium values for N (Vector) and C (Matrix of Matrix). Here C=N*N
    Vector GenerateEquilibriumValuesForNAndC(CalcType value)
    {
        int N_size = props.SpeciesNumber;
        int C_size = props.SpeciesNumber * props.SpeciesNumber * props.MomentsGrid.size() * props.MomentsGrid.size();
        
        Vector result;  
        for (int i=0; i<N_size; i++)
            result.add(value);
        for (int i=0; i<C_size; i++)
            result.add(value * value);
        return result;
    }

    // Start Calculations
    Vector operator () (CalcType times, Vector vals, ExtrapolationType chosen_type)
    {
        C_extrapolation_type = chosen_type;

        N = vals.head(props.SpeciesNumber);
        vals.erase_first(props.SpeciesNumber);
        C = VectorToMatrixOfMatrix(vals, props.SpeciesNumber, props.MomentsGrid.size());

        Vector dN(props.SpeciesNumber);                                 // Vector of derivatives of current first moments
        MatrixOfMatrix dC(props.SpeciesNumber);                         // Derivatives of current second moments

        MatrixOfMatrix WC = CalculateIntegralOfWAndC(C, N);             // Integral of W and C
        MatrixOfMatrix WC2 = CalculateIntegralOfWAndC(C, N, true);      // Integral of W and C 2
        MatrixOfMatrix MC = CalculateIntegralOfMAndC(C, N);             // Integral of M and C

        //.......................................CALCULATE FIRST MOMENT
        dN = N * (props.b - props.d);                                   // Contribution of birth and death, density independent
        Vector d_prime_sum(props.SpeciesNumber, 0);                     // Death contribution, density dependent
        for (int i=0; i<props.SpeciesNumber; i++)
        {
            for (int j=0; j<props.SpeciesNumber; j++)
            {    
                d_prime_sum[i] = d_prime_sum[i] + props.d_prime[i][j] * WC(i, j)[0][0];
            }
        }
        dN = dN - d_prime_sum;

        //.......................................CALCULATE SECOND MOMENT
        for (int i=0; i<props.SpeciesNumber; i++)
        {
            for (int j=0; j<props.SpeciesNumber; j++)
            {
                Matrix dC_elem(props.MomentsGrid.size());

                for (int x_num=0; x_num<props.MomentsGrid.size(); x_num++)
                {
                    for (int y_num=0; y_num<props.MomentsGrid.size(); y_num++)
                    {
                        // Birth contribution, density independent
                        CalcType result = props.b[i] * MC(i, j)[x_num][y_num] + props.b[j] * MC(j, i)[x_num][y_num]; 

                        // Birth contribution, Kronecker symbols
                        CalcType temp;
                        CalcType interp;
                        if (i == j)
                        {
                            interp = BilinearInterpolation(
                                props.KernelsGrid,
                                props.M(i),
                                props.MomentsGrid[x_num],
                                props.MomentsGrid[y_num],
                                props.SigmaM[i] * 3,
                                Zero
                            );
                            temp = 2 * props.b[i] * N[i] * interp;
                            result += temp;
                        }

                        // Death contribution, density independent
                        result -= props.d[i] * C(i, j)[x_num][y_num];
                        result -= props.d[j] * C(j, i)[x_num][y_num];

                        //Parametric Power 2 closure
                        CalcType closure = params[0] * C(0, 0)[x_num][y_num] * WC(0, 0)[0][0] / N[0];
                        closure += params[1] * C(0, 0)[x_num][y_num] * WC(0, 0)[x_num][y_num] / N[0];
                        closure += params[2] * WC2(0, 0)[x_num][y_num] / N[0];
                        closure -= params[1] * N[0] * N[0] * N[0];
                        closure = closure * props.d_prime[0][0];
                        closure = closure / (params[0] + params[2]);
                        result -= 2 * closure;
                        
                        interp = BilinearInterpolation(
                            props.KernelsGrid,
                            props.W(i, j),
                            props.MomentsGrid[x_num],
                            props.MomentsGrid[y_num],
                            props.SigmaW[i][j] * 3,
                            Zero
                        );

                        temp = props.d_prime[i][j] * C(i, j)[x_num][y_num] * BilinearInterpolation(
                            props.KernelsGrid,
                            props.W(i, j),
                            props.MomentsGrid[x_num],
                            props.MomentsGrid[y_num],
                            props.SigmaW[i][j] * 3,
                            Zero
                        );
                        temp += props.d_prime[j][i] * C(j, i)[x_num][y_num] * BilinearInterpolation(
                            props.KernelsGrid,
                            props.W(j, i),
                            props.MomentsGrid[x_num],
                            props.MomentsGrid[y_num],
                            props.SigmaW[j][i] * 3,
                            Zero
                        );
                        result -= temp;
                        
                        dC_elem[x_num][y_num] = result;
                    }
                }
                dC(i, j) = dC_elem;
            }
        }
        return MakeVector(dN, dC);
    }
};


/**************************************************************************************************************************
**                                                         SOLVERS                                                       **
***************************************************************************************************************************
**         The code below realizes solvers for the logistic equation using vairous numerical methods, as well as         **
**                 a structure of all calculations' final results that will be used for creating plots.                  **
**************************************************************************************************************************/

//_____________________________________________________SOLVER RESULTS______________________________________________________
struct SolverResults
{
    Vector Times;                                                       // Time points
    vector<Vector> Values;                                              // Values of the solution at time
    vector<Vector> Derivatives;                                         // Derivatives

    void clear()                                                        // Clear solver results
    {
        Times.clear();
        Values.clear();
        Derivatives.clear();
    }
};


//___________________________________________________EULER METHOD SOLVER___________________________________________________
class EulerSolver
{
    SolverResults solver_results;

private:
    void Add(Vector value, CalcType time)
    {
        solver_results.Values.emplace_back(std::move(value));
        solver_results.Times.add(time);
    }

    void AddDerivative(Vector derivative)
    { solver_results.Derivatives.emplace_back(std::move(derivative)); }

public:
    // Solve given equation numericaly for fixed amount of iterations
    SolverResults Solve(
        Equation current_equation,
        ExtrapolationType type,
        CalcType time_max,
        CalcType step,
        Vector initial_value,
        bool debug_print_enabled = false
    ) {
        CalcType time = 0;
        Vector value = initial_value;
        Vector derivative = current_equation(time, value, type);

        Add(value, time);
        AddDerivative(derivative);
        if (debug_print_enabled)
                cout << "   t = " << time << ";   N = " << value.first() << '\n';
        
        time += step;
        while (time < time_max)
        {
            value = value + derivative * step;
            for (auto& val : value)
                if (val < 0)
                    val = 0;
            Add(value, time);
            derivative = current_equation(time, value, type);
            AddDerivative(derivative);
            
            if (debug_print_enabled)
                cout << "   t = " << time << ";   N = " << value.first() << '\n';

            time += step;
        }
    
        return solver_results;
    }
};


/**************************************************************************************************************************
**                                                    INPUT AND OUTPUT                                                   **
***************************************************************************************************************************
**          The code below realizes means of retreiving spatial model parameters from input data ".csv" file,            **
**                       saving the calculations' results in a file and visualising them in a plot.                      **
**************************************************************************************************************************/

//_____________________________________________________FILE READER_________________________________________________________
class FileReader
{
    fstream file;
    string filename;
    string new_line;
    string new_word;

private:
    CalcType GetCalcType()
    {
        getline(file, new_line);
        stringstream stream(new_line);

        getline(stream, new_word, ',');
        getline(stream, new_word, ',');
        return stod(new_word);
    }

    Vector GetVector(int vec_size)
    {
        Vector vector_from_file;
        getline(file, new_line);
        stringstream stream(new_line);

        getline(stream, new_word, ',');
        while (getline(stream, new_word, ','))
            vector_from_file.add(stod(new_word));
        assert (vector_from_file.size() == vec_size);
        return vector_from_file;
    }

    Matrix GetMatrix(int mat_size)
    {
        Matrix matrix_from_file(mat_size);
        for (int row_num=0; row_num<mat_size; row_num++)
        {
            getline(file, new_line);
            stringstream stream(new_line);

            int col_num = 0;
            getline(stream, new_word, ',');
            while (getline(stream, new_word, ','))
            {
                matrix_from_file[row_num][col_num] = stod(new_word);
                col_num++;
            }
        }
        for (int i=mat_size; i<3; i++)
            getline(file, new_line);
        
        return matrix_from_file;
    }

    void DebugPrint(Properties props)
    {
        cout << "MODEL PROPERTIES:\n"
             << "   Species number: " << props.SpeciesNumber << '\n'
             << "   Cutoff distance : " << props.CutoffDistance << '\n'
             << "   Sigma W:  " << props.SigmaW
             << "   Sigma M:  " << props.SigmaM
             << "         b:  " << props.b
             << "         d:  " << props.d
             << "         d': " << props.d_prime
             << "\n\n";
    }

public:
    FileReader(string str) : filename(str) {}

    pair<Properties, Vector> GetWithParameters(bool debug_print_enabled = false)
    {
        Vector params(0);
        file.open(filename);
        
        auto b = GetCalcType();
        auto d = GetCalcType();
        auto dd = GetCalcType();
        
        getline(file, new_line);
        auto sW = GetCalcType();

        getline(file, new_line);
        auto sM = GetCalcType();

        getline(file, new_line);
        CalcType alpha = GetCalcType();
        CalcType beta = GetCalcType();
        CalcType gamma = GetCalcType();
        params.add(alpha);
        params.add(beta);
        params.add(gamma);
        
        auto properties = Properties::Generate(
            1,
            Matrix(1, sW),
            Vector(1, sM),
            Vector(1, b),
            Vector(1, d),
            Matrix(1, dd)
        );
        if (debug_print_enabled)
        {
            DebugPrint(properties);
            cout << "Parameters:"
                 << "\n   Alpha: " << params[0]
                 << "\n   Beta: " << params[1]
                 << "\n   Gamma: " << params[2] << "\n\n\n";
        }
        properties.CheckKernels();

        file.close();
        return make_pair(properties, params);
    }
};


//_____________________________________________________FILE WRITER_________________________________________________________
class FileWriter
{
    fstream file;
    string filename;

public:
    FileWriter(string str) : filename(str) {}

    void WritePlotData(string number, string title, vector<string> labels, vector<Vector> points)
    {
        assert (points.size() >= 2);
        assert (points.size() == labels.size());

        int number_of_plots = points.size() - 1;
        int number_of_points = points.at(0).size();
        for (int i=1; i<points.size(); i++)
            assert (points.at(i).size() == number_of_points);

        file.open(filename, ios::out);
        file.clear();
        file << "Experiment" << ',' << number << '\n'
             << "Plot title" << ',' << title << '\n'
             << '\n'
             << "Number of" << ',' << number_of_plots << '\n'
             << "plots" << '\n'
             << '\n'
             << "Number of" << ',' << number_of_points << '\n'
             << "points" << '\n'
             << '\n'
             << labels.at(0);
        for (int i=1; i<labels.size(); i++)
            file << ',' << labels.at(i);
        file << '\n';

        for (int point_num=0; point_num<number_of_points; point_num++)
        {
            file << points.at(0)[point_num];
            for (int plot_num=1; plot_num<=number_of_plots; plot_num++)
                file << ',' << points.at(plot_num)[point_num];
            file << '\n';
        }
        
        file.close();
    }
};


//_____________________________________________________EXPERIMENTS_________________________________________________________
class Experiments
{
    Properties Props;
    string reader_file_path;
    string visdata_file_path;
    string visdatac_file_path;
    bool enable_debug_print;

private:
    // Sort pairs of values by the first value
    static bool SortByFirst(const pair<CalcType, CalcType> &a, const pair<CalcType, CalcType> &b)
    { return (a.first < b.first); }

    // Calculate first moment
    SolverResults CalculateFirstMoment(CalcType t_max, CalcType initN = 200, ExtrapolationType type=SquaredN)
    {
        FileReader Reader(reader_file_path);

        auto PropsAndParams = Reader.GetWithParameters(enable_debug_print=true);
        Properties Props = PropsAndParams.first;
        Vector Params = PropsAndParams.second;
        Equation equation(Props, Params);
        EulerSolver solver;

        cout << "\n==========================CALCULATIONS START==========================\n";
        clock_t Start = clock();
        
        auto solver_results = solver.Solve(equation, type, t_max, 0.1, equation.GenerateEquilibriumValuesForNAndC(initN), enable_debug_print=true);
        
        clock_t End = clock();
        cout << "\n===========================CALCULATIONS END===========================\n";
        
        Time calculation_time(End - Start);
        cout << "TIME:  " << calculation_time << '\n';

        return solver_results;
    }

    // Launch Visualiser.py     !! THIS FUNCTION WORKS ONLY ON WINDOWS !!
    void Visualise_Windows()
    {
        WinExec("python Visualiser.py", 1);
    }

public:
    Experiments() : reader_file_path("Data\\Properties.csv"), visdata_file_path("Data\\VisData.csv"), visdatac_file_path("Data\\VisDataC.csv") {}

    void Custom_Experiment()
    {
        CalcType max_time = 200;
        CalcType N_init = 200;
        auto solver_results = CalculateFirstMoment(max_time, N_init);
    
        Vector N;
        Vector CSquareRoot;
        for (auto Values : solver_results.Values)
        {
            CalcType N_value = Values.first();
            N.add(N_value);

            CalcType C_at_cutoff_distance = Values[Props.MomentsGrid.size()];
            CSquareRoot.add(std::sqrt(C_at_cutoff_distance));
        }


        vector<vector<CalcType>> C_result;
        bool unique;

        for (auto Values : solver_results.Values)
        {
            vector<pair<CalcType, CalcType>> C_vals_vec;
            auto C_vec = Values;
            C_vec.erase_first(1);
            MatrixOfMatrix C_mat = VectorToMatrixOfMatrix(C_vec, Props.SpeciesNumber, Props.MomentsGrid.size());
            vector<CalcType> C_vals;
            for (int x=0; x<Props.MomentsGrid.size(); x++)
            {
                for (int y=0; y<Props.MomentsGrid.size(); y++)
                {
                    unique = true;
                    CalcType r = std::hypot(Props.MomentsGrid[x], Props.MomentsGrid[y]);
                    for (auto item : C_vals_vec)
                    {
                        if (r == item.first)
                        {
                            unique = false;
                        }
                    }
                    if (unique)
                    {
                        C_vals_vec.push_back(make_pair(r, C_mat(0, 0)[x][y]));
                    }
                }
            }
            sort(C_vals_vec.begin(), C_vals_vec.end(), SortByFirst);
            C_vals.clear();
            for (auto item : C_vals_vec)
            {
                C_vals.push_back(item.second);
            }
            C_result.push_back(C_vals);
        }

        vector<Vector> C_plot_points;
        vector<string> C_labels;
        C_plot_points.push_back(solver_results.Times);
        C_labels.push_back("Time");
        int numberOfPlots = C_result.at(0).size();
        for (int i=0; i<numberOfPlots; i++)
        {
            Vector points;
            for (auto values : C_result)
            {
                points.add(values.at(i));
            }
            C_plot_points.push_back(points);
            C_labels.push_back(" ");
        }

        FileWriter CVisData(visdatac_file_path);
        CVisData.WritePlotData("1", "C", C_labels, C_plot_points);

        
        vector<Vector> plot_points;
        plot_points.push_back(solver_results.Times);
        plot_points.push_back(N);
        plot_points.push_back(CSquareRoot);

        vector<string> plot_labels;
        plot_labels.push_back("Time");
        plot_labels.push_back("N");
        plot_labels.push_back("Square root of C at cutoff distance");

        FileWriter VisData(visdata_file_path);
        VisData.WritePlotData("1", "N", plot_labels, plot_points);
        Visualise_Windows();
    }

    void KernelTest(char KernelSymbol, string case_str)
    {
        FileReader Reader(reader_file_path);

        auto PropsAndParams = Reader.GetWithParameters(enable_debug_print = true);
        Properties Props = PropsAndParams.first;
        Vector Params = PropsAndParams.second;
        
        Matrix kernel;
        string title;
        switch(KernelSymbol)
        {
            case 'W':
                kernel = Props.W(0, 0);
                title = "w(r) kernel test for experiment " + case_str;
                break;
        
            case 'M':
                kernel = Props.M(0);
                title = "m(r) kernel test for experiment " + case_str;
                break;
            
            default:
                cout << "ERROR: No such kernel\n";
                break;
        }
    
        vector<pair<CalcType, CalcType>> result;
        bool unique;

        for (int x=0; x<Props.KernelsGrid.size(); x++)
        {
            for (int y=0; y<Props.KernelsGrid.size(); y++)
            {
                unique = true;
                CalcType r = std::hypot(Props.KernelsGrid[x], Props.KernelsGrid[y]);
                for (auto item : result)
                {
                    //cout << "x = " << Props.KernelsGrid[x] << "    y = " << Props.KernelsGrid[y] << "     r = " << r << "    kernel = " << kernel[x][y] << '\n';
                    if (r == item.first)
                        unique = false;
                }
                if (unique)
                    result.push_back(make_pair(r, kernel[x][y]));
            }
        }
        sort(result.begin(), result.end(), SortByFirst);

        vector<CalcType> xValues;
        vector<CalcType> yValues;
        for (auto item : result)
        {
            xValues.push_back(item.first);
            yValues.push_back(item.second);
        }
        
        switch(KernelSymbol)
        {
            case 'W':
                cout << "MAX NON-ZERO NUMBERS: " << Props.MaxNumForW << '\n';
                break;
            
            case 'M':
                cout << "MAX NON-ZERO NUMBERS: " << Props.MaxNumForM << '\n';
                break;

            default:
                cout << "ERROR: No such kernel\n";
                break;
        }

        vector<Vector> kernel_plot_points;
        kernel_plot_points.push_back(xValues);
        kernel_plot_points.push_back(yValues);

        vector<string> kernel_plot_labels;
        kernel_plot_labels.push_back("r");
        kernel_plot_labels.push_back("distribution");

        FileWriter VisData(visdata_file_path);
        VisData.WritePlotData("1", title, kernel_plot_labels, kernel_plot_points);
        Visualise_Windows();
    }
};


/**************************************************************************************************************************
**                                                        MAIN BODY                                                      **
**************************************************************************************************************************/
int main()
{
    Experiments run;
    string MSG = "Welcome!\nPlease navigate to \'Data\' folder and fill the \'Properties.csv\' file.\nPress ENTER to confirm and start the experiment\nPress ESC to quit\n";
    
    system("Cls");
    cout << MSG;
    char input = getch();

        while (input != KEY_ESC)
        {
            if (input == KEY_ENT)
            {
                system("cls");
                run.Custom_Experiment();
                cout << MSG;
                input = getch();
            }
        }
        system("cls");   

    return 0;
}
