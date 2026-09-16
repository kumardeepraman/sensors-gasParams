// ============================================================
// C2F6 Rayleigh Scattering Length Extrapolation
//
// Input:
//   C2F6_Rayleigh_Absorption.csv
//
// Output:
//   C2F6_Rayleigh_1000nm.txt
//
// Model for extrapolation:
//   L_R(lambda) = A * lambda^n
//
// Measured values are preserved.
// Linear interpolation is used only between measured points.
// Extrapolation is performed only above the last measured point.
//
// Units:
//   wavelength : nm
//   length     : m
//   energy     : eV
// ============================================================

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iomanip>

using namespace std;

// ------------------------------------------------------------
// Constants
// ------------------------------------------------------------
const double HC = 1239.841984;   // eV*nm

const double OUTPUT_MIN_WL = 190.0;
const double OUTPUT_MAX_WL = 1000.0;
const double OUTPUT_STEP   = 1.0;

// Rayleigh fit range
const double FIT_MIN_WL = 400.0;
const double FIT_MAX_WL = 785.5;


// ------------------------------------------------------------
// Data structure
// ------------------------------------------------------------
struct OpticalPoint
{
    double wavelength;
    double energy;
    double absorption;
    double rayleigh;
};


// ------------------------------------------------------------
// Read CSV
// ------------------------------------------------------------
vector<OpticalPoint> ReadCSV(const char *filename)
{
    vector<OpticalPoint> data;

    ifstream file(filename);

    if (!file.is_open())
    {
        cerr << "ERROR: Cannot open file: "
             << filename << endl;
        return data;
    }

    string line;

    // --------------------------------------------------------
    // The CSV header occupies 4 physical lines because
    // the column names contain embedded newlines.
    // --------------------------------------------------------
    for (int i = 0; i < 4; i++)
        getline(file, line);

    while (getline(file, line))
    {
        if (line.empty())
            continue;

        stringstream ss(line);

        string s_wavelength;
        string s_energy;
        string s_absorption;
        string s_rayleigh;

        if (!getline(ss, s_wavelength, ','))
            continue;

        if (!getline(ss, s_energy, ','))
            continue;

        if (!getline(ss, s_absorption, ','))
            continue;

        if (!getline(ss, s_rayleigh, ','))
            continue;

        try
        {
            OpticalPoint p;

            p.wavelength = stod(s_wavelength);
            p.energy     = stod(s_energy);
            p.absorption = stod(s_absorption);
            p.rayleigh   = stod(s_rayleigh);

            data.push_back(p);
        }
        catch (...)
        {
            cerr << "WARNING: Could not read line:\n"
                 << line << endl;
        }
    }

    file.close();

    // --------------------------------------------------------
    // Always sort by wavelength
    // --------------------------------------------------------
    sort(data.begin(), data.end(),
         [](const OpticalPoint &a, const OpticalPoint &b)
         {
             return a.wavelength < b.wavelength;
         });

    return data;
}


// ------------------------------------------------------------
// Find exact measured point
// ------------------------------------------------------------
bool GetExactMeasuredValue(const vector<OpticalPoint> &data,
                            double wavelength,
                            double &rayleigh)
{
    const double tolerance = 1.0e-8;

    for (const auto &p : data)
    {
        if (fabs(p.wavelength - wavelength) < tolerance)
        {
            rayleigh = p.rayleigh;
            return true;
        }
    }

    return false;
}


// ------------------------------------------------------------
// Linear interpolation
//
// This is used ONLY inside the measured wavelength range.
// ------------------------------------------------------------
double InterpolateRayleigh(const vector<OpticalPoint> &data,
                           double wavelength)
{
    if (data.empty())
        return 0.0;

    // Below first measured point
    if (wavelength <= data.front().wavelength)
        return data.front().rayleigh;

    // Above last measured point
    if (wavelength >= data.back().wavelength)
        return data.back().rayleigh;

    for (size_t i = 0; i < data.size() - 1; i++)
    {
        double x1 = data[i].wavelength;
        double x2 = data[i + 1].wavelength;

        if (wavelength >= x1 && wavelength <= x2)
        {
            double y1 = data[i].rayleigh;
            double y2 = data[i + 1].rayleigh;

            double fraction =
                (wavelength - x1) / (x2 - x1);

            return y1 + fraction * (y2 - y1);
        }
    }

    return data.back().rayleigh;
}


// ------------------------------------------------------------
// Fit Rayleigh:
//
// log(L) = log(A) + n*log(lambda)
//
// Therefore:
// y = c + n*x
// ------------------------------------------------------------
bool FitRayleigh(const vector<OpticalPoint> &data,
                 double &A,
                 double &n)
{
    double Sx  = 0.0;
    double Sy  = 0.0;
    double Sxx = 0.0;
    double Sxy = 0.0;

    int N = 0;

    for (const auto &p : data)
    {
        if (p.wavelength < FIT_MIN_WL ||
            p.wavelength > FIT_MAX_WL)
            continue;

        if (p.rayleigh <= 0.0)
            continue;

        double x = log(p.wavelength);
        double y = log(p.rayleigh);

        Sx  += x;
        Sy  += y;
        Sxx += x*x;
        Sxy += x*y;

        N++;
    }

    if (N < 3)
    {
        cerr << "ERROR: Not enough points for Rayleigh fit."
             << endl;
        return false;
    }

    double denominator =
        N*Sxx - Sx*Sx;

    if (fabs(denominator) < 1e-20)
    {
        cerr << "ERROR: Singular Rayleigh fit."
             << endl;
        return false;
    }

    n = (N*Sxy - Sx*Sy) / denominator;

    double intercept =
        (Sy - n*Sx) / N;

    A = exp(intercept);

    cout << endl;
    cout << "============================================"
         << endl;
    cout << "Rayleigh Power-Law Fit"
         << endl;
    cout << "============================================"
         << endl;

    cout << "Fit range       : "
         << FIT_MIN_WL << " - "
         << FIT_MAX_WL << " nm" << endl;

    cout << "Number of points: "
         << N << endl;

    cout << "A               : "
         << scientific << setprecision(8)
         << A << endl;

    cout << "n               : "
         << fixed << setprecision(6)
         << n << endl;

    cout << endl;

    return true;
}


// ------------------------------------------------------------
// Main
// ------------------------------------------------------------
void extrapolate_C2F6_Rayleigh()
{
    const char *inputFile =
        "C2F6_Rayleigh_Absorption.csv";

    const char *outputFile =
        "C2F6_Rayleigh_1000nm.txt";

    // --------------------------------------------------------
    // Read data
    // --------------------------------------------------------
    vector<OpticalPoint> data =
        ReadCSV(inputFile);

    if (data.empty())
    {
        cerr << "ERROR: No data read." << endl;
        return;
    }

    cout << endl;
    cout << "============================================"
         << endl;
    cout << "C2F6 Rayleigh Extrapolation"
         << endl;
    cout << "============================================"
         << endl;

    cout << "Number of data points : "
         << data.size() << endl;

    cout << "Minimum wavelength    : "
         << data.front().wavelength
         << " nm" << endl;

    cout << "Maximum wavelength    : "
         << data.back().wavelength
         << " nm" << endl;


    // --------------------------------------------------------
    // Fit
    // --------------------------------------------------------
    double A = 0.0;
    double n = 0.0;

    if (!FitRayleigh(data, A, n))
        return;


    // --------------------------------------------------------
    // Open output file
    // --------------------------------------------------------
    ofstream out(outputFile);

    if (!out.is_open())
    {
        cerr << "ERROR: Cannot open output file."
             << endl;
        return;
    }

    out << "# C2F6 Rayleigh scattering length\n";
    out << "# Energy[eV]   RayleighLength[m]\n";
    out << "# Rayleigh fit: L = A * lambda^n\n";
    out << "# A = " << setprecision(12)
        << A << "\n";
    out << "# n = " << setprecision(12)
        << n << "\n";
    out << "#\n";


    // --------------------------------------------------------
    // Generate output
    //
    // IMPORTANT:
    // Exact measured values are used whenever the wavelength
    // exists in the CSV.
    //
    // For wavelengths inside the measured range but not
    // explicitly measured, linear interpolation is used.
    //
    // Above the last measured wavelength, the power-law fit
    // is used.
    // --------------------------------------------------------
    for (double wavelength = OUTPUT_MAX_WL;
         wavelength >= OUTPUT_MIN_WL - 1e-9;
         wavelength -= OUTPUT_STEP)
    {
        double energy =
            HC / wavelength;

        double rayleighLength = 0.0;

        double exactValue = 0.0;

        // ----------------------------------------------------
        // Exact measured point
        // ----------------------------------------------------
        if (GetExactMeasuredValue(data,
                                  wavelength,
                                  exactValue))
        {
            rayleighLength = exactValue;
        }

        // ----------------------------------------------------
        // Extrapolation above measured range
        // ----------------------------------------------------
        else if (wavelength > data.back().wavelength)
        {
            rayleighLength =
                A * pow(wavelength, n);
        }

        // ----------------------------------------------------
        // Interpolation inside measured range
        // ----------------------------------------------------
        else
        {
            rayleighLength =
                InterpolateRayleigh(data, wavelength);
        }


        // ----------------------------------------------------
        // Geant4 style output
        // ----------------------------------------------------
        out << fixed << setprecision(6)
            << energy << "*eV   "
            << rayleighLength << "*m"
            << endl;
    }

    out.close();


    // --------------------------------------------------------
    // 1000 nm result
    // --------------------------------------------------------
    double rayleigh1000 =
        A * pow(1000.0, n);

    cout << "============================================"
         << endl;

    cout << "Rayleigh length at 1000 nm = "
         << scientific << setprecision(8)
         << rayleigh1000
         << " m" << endl;

    cout << "                     = "
         << rayleigh1000 / 1000.0
         << " km" << endl;

    cout << endl;

    cout << "Output file:"
         << endl;
    cout << "  " << outputFile
         << endl;

    cout << "============================================"
         << endl;
}
