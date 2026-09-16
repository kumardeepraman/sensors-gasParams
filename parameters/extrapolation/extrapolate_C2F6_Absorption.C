// ============================================================
// C2F6 Absorption Length
//
// Input:
//   C2F6_Rayleigh_Absorption.csv
//
// Output:
//   C2F6_Absorption_1000nm.txt
//
// IMPORTANT:
//
//   No absorption extrapolation fit is used.
//
//   For wavelengths ABOVE the last measured wavelength:
//       Absorption Length = 1350.0 m
//
//   Inside the measured range:
//       - exact measured values are preserved
//       - linear interpolation is used only if necessary
//
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

// High-wavelength absorption plateau
const double HIGH_WAVELENGTH_ABSORPTION =
    1350.0;   // m


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
    // Header occupies 4 physical lines
    // because the CSV column names contain newlines.
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
    // Sort by wavelength
    // --------------------------------------------------------
    sort(data.begin(), data.end(),
         [](const OpticalPoint &a,
            const OpticalPoint &b)
         {
             return a.wavelength < b.wavelength;
         });

    return data;
}


// ------------------------------------------------------------
// Exact measured value
// ------------------------------------------------------------
bool GetExactMeasuredValue(
        const vector<OpticalPoint> &data,
        double wavelength,
        double &absorption)
{
    const double tolerance = 1.0e-8;

    for (const auto &p : data)
    {
        if (fabs(p.wavelength - wavelength)
            < tolerance)
        {
            absorption = p.absorption;
            return true;
        }
    }

    return false;
}


// ------------------------------------------------------------
// Linear interpolation
//
// Used only inside the measured range.
// ------------------------------------------------------------
double InterpolateAbsorption(
        const vector<OpticalPoint> &data,
        double wavelength)
{
    if (data.empty())
        return 0.0;


    // Below first measured wavelength
    if (wavelength <= data.front().wavelength)
        return data.front().absorption;


    // Above last measured wavelength
    if (wavelength >= data.back().wavelength)
        return data.back().absorption;


    for (size_t i = 0;
         i < data.size() - 1;
         i++)
    {
        double x1 = data[i].wavelength;
        double x2 = data[i + 1].wavelength;


        if (wavelength >= x1 &&
            wavelength <= x2)
        {
            double y1 = data[i].absorption;
            double y2 = data[i + 1].absorption;


            double fraction =
                (wavelength - x1)
                / (x2 - x1);


            return y1 +
                   fraction * (y2 - y1);
        }
    }


    return data.back().absorption;
}


// ------------------------------------------------------------
// Main
// ------------------------------------------------------------
void extrapolate_C2F6_Absorption()
{
    const char *inputFile =
        "C2F6_Rayleigh_Absorption.csv";

    const char *outputFile =
        "C2F6_Absorption_1000nm.txt";


    // --------------------------------------------------------
    // Read data
    // --------------------------------------------------------
    vector<OpticalPoint> data =
        ReadCSV(inputFile);


    if (data.empty())
    {
        cerr << "ERROR: No data read."
             << endl;
        return;
    }


    cout << endl;
    cout << "============================================"
         << endl;
    cout << "C2F6 Absorption Length"
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

    cout << "High-wavelength value : "
         << HIGH_WAVELENGTH_ABSORPTION
         << " m" << endl;


    // --------------------------------------------------------
    // Open output
    // --------------------------------------------------------
    ofstream out(outputFile);


    if (!out.is_open())
    {
        cerr << "ERROR: Cannot open output file."
             << endl;
        return;
    }


    out << "# C2F6 absorption length\n";
    out << "# Energy[eV]   AbsorptionLength[m]\n";
    out << "#\n";
    out << "# Above "
        << data.back().wavelength
        << " nm: L_abs = "
        << HIGH_WAVELENGTH_ABSORPTION
        << " m\n";
    out << "#\n";


    // --------------------------------------------------------
    // Generate output
    //
    // 1000 nm --> 190 nm
    // --------------------------------------------------------
    for (double wavelength = OUTPUT_MAX_WL;
         wavelength >= OUTPUT_MIN_WL - 1e-9;
         wavelength -= OUTPUT_STEP)
    {
        double energy =
            HC / wavelength;


        double absorptionLength = 0.0;

        double exactValue = 0.0;


        // ----------------------------------------------------
        // 1. EXACT measured value
        //
        // This guarantees that the CSV value is not modified.
        // ----------------------------------------------------
        if (GetExactMeasuredValue(
                data,
                wavelength,
                exactValue))
        {
            absorptionLength = exactValue;
        }


        // ----------------------------------------------------
        // 2. HIGH-WAVELENGTH REGION
        //
        // Above the last measured point:
        //
        //       L_abs = 1350 m
        //
        // exactly.
        // ----------------------------------------------------
        else if (wavelength >
                 data.back().wavelength)
        {
            absorptionLength =
                HIGH_WAVELENGTH_ABSORPTION;
        }


        // ----------------------------------------------------
        // 3. Inside measured range
        //
        // Linear interpolation only where the requested
        // wavelength is not directly present in the CSV.
        // ----------------------------------------------------
        else
        {
            absorptionLength =
                InterpolateAbsorption(
                    data,
                    wavelength);
        }


        // ----------------------------------------------------
        // Geant4 style
        // ----------------------------------------------------
        out << fixed << setprecision(6)
            << energy << "*eV   "
            << absorptionLength << "*m"
            << endl;
    }


    out.close();


    // --------------------------------------------------------
    // Explicit 1000 nm result
    // --------------------------------------------------------
    cout << endl;
    cout << "============================================"
         << endl;

    cout << "Absorption length at 1000 nm = "
         << fixed << setprecision(3)
         << HIGH_WAVELENGTH_ABSORPTION
         << " m"
         << endl;

    cout << endl;

    cout << "Therefore:"
         << endl;

    cout << "  785.5 nm < lambda <= 1000 nm"
         << endl;

    cout << "  L_abs = 1350.0 m"
         << endl;

    cout << endl;

    cout << "Output file:"
         << endl;

    cout << "  " << outputFile
         << endl;

    cout << "============================================"
         << endl;
}
