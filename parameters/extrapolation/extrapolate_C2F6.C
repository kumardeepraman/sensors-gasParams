
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <string>

#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TAxis.h"

using namespace std;


// ============================================================
// DATA STRUCTURE
// ============================================================

struct C2F6Data
{
    double wavelength;   // nm
    double energy;       // eV
    double absorption;   // m
    double rayleigh;     // m
};


// ============================================================
// CONSTANTS
// ============================================================

// h*c = 1239.841984 eV*nm
const double HC = 1239.841984;

// Output wavelength range
const double MIN_WAVELENGTH = 190.0;     // nm
const double MAX_WAVELENGTH = 1000.0;    // nm

// Output step
const double WAVELENGTH_STEP = 1.0;      // nm


// ============================================================
// RAYLEIGH FIT SETTINGS
// ============================================================
//
// L_R(lambda) = A_R * lambda^n_R
//
// We fit:
//
// ln(L_R) = ln(A_R) + n_R * ln(lambda)
//
// The measured C2F6 data show an excellent lambda^4 behavior,
// therefore the fit is performed over the long-wavelength
// measured region.
//

const double RAYLEIGH_FIT_MIN = 400.0;   // nm
const double RAYLEIGH_FIT_MAX = 785.5;   // nm


// ============================================================
// ABSORPTION FIT SETTINGS
// ============================================================
//
// For extrapolation we use:
//
// L_abs(lambda) = A_A * lambda^n_A
//
// or:
//
// ln(L_abs) = ln(A_A) + n_A * ln(lambda)
//
// Saturated values are excluded from the fit.
//
// The original dataset contains values reaching 1350 m.
// These are treated as upper-limit/capped values rather than
// genuine measurements.
//

const double ABSORPTION_SATURATION = 1300.0;


// ============================================================
// READ CSV
// ============================================================

vector<C2F6Data> ReadCSV(const char* filename)
{
    vector<C2F6Data> data;

    ifstream file(filename);

    if (!file.is_open())
    {
        cerr << "\nERROR: Cannot open input file:\n"
             << filename << endl;

        return data;
    }

    string line;


    // --------------------------------------------------------
    // The CSV header occupies four physical lines because
    // the column names contain embedded newlines.
    //
    // "Wavelength
    // (nm)","Energy
    // (eV)","Absorption
    // Length (m)","Rayleigh
    // Length (m)"
    // --------------------------------------------------------

    for (int i = 0; i < 4; ++i)
        getline(file, line);


    // --------------------------------------------------------
    // Read data
    // --------------------------------------------------------

    while (getline(file, line))
    {
        if (line.empty())
            continue;

        stringstream ss(line);

        string s1, s2, s3, s4;

        if (!getline(ss, s1, ','))
            continue;

        if (!getline(ss, s2, ','))
            continue;

        if (!getline(ss, s3, ','))
            continue;

        if (!getline(ss, s4, ','))
            continue;


        try
        {
            C2F6Data p;

            p.wavelength = stod(s1);
            p.energy     = stod(s2);
            p.absorption = stod(s3);
            p.rayleigh   = stod(s4);

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
         [](const C2F6Data& a, const C2F6Data& b)
         {
             return a.wavelength < b.wavelength;
         });


    // --------------------------------------------------------
    // Print information
    // --------------------------------------------------------

    cout << "\n====================================================\n";
    cout << "C2F6 DATA\n";
    cout << "====================================================\n";

    cout << "Number of points = "
         << data.size() << endl;

    if (!data.empty())
    {
        cout << "Minimum wavelength = "
             << data.front().wavelength
             << " nm\n";

        cout << "Maximum wavelength = "
             << data.back().wavelength
             << " nm\n";
    }

    cout << "====================================================\n";


    return data;
}


// ============================================================
// RAYLEIGH POWER-LAW FIT
// ============================================================
//
// L_R(lambda) = A_R * lambda^n_R
//
// ln(L_R) = ln(A_R) + n_R * ln(lambda)
// ============================================================

bool FitRayleigh(
    const vector<C2F6Data>& data,
    double& A,
    double& n)
{
    double Sx  = 0.0;
    double Sy  = 0.0;
    double Sxx = 0.0;
    double Sxy = 0.0;

    int N = 0;


    for (const auto& p : data)
    {
        if (p.wavelength < RAYLEIGH_FIT_MIN)
            continue;

        if (p.wavelength > RAYLEIGH_FIT_MAX)
            continue;

        if (p.rayleigh <= 0.0)
            continue;


        double x = log(p.wavelength);
        double y = log(p.rayleigh);

        Sx  += x;
        Sy  += y;
        Sxx += x*x;
        Sxy += x*y;

        ++N;
    }


    if (N < 3)
    {
        cerr << "\nERROR: Not enough points for Rayleigh fit."
             << endl;

        A = 0.0;
        n = 0.0;

        return false;
    }


    double denominator =
        N*Sxx - Sx*Sx;


    if (fabs(denominator) < 1e-20)
    {
        cerr << "\nERROR: Singular Rayleigh fit."
             << endl;

        A = 0.0;
        n = 0.0;

        return false;
    }


    n =
        (N*Sxy - Sx*Sy)
        / denominator;


    double intercept =
        (Sy - n*Sx) / N;


    A = exp(intercept);


    cout << "\n====================================================\n";
    cout << "RAYLEIGH FIT\n";
    cout << "====================================================\n";

    cout << "Fit range = "
         << RAYLEIGH_FIT_MIN
         << " - "
         << RAYLEIGH_FIT_MAX
         << " nm\n";

    cout << "Number of points = "
         << N << endl;

    cout << fixed << setprecision(6);

    cout << "Power n = "
         << n << endl;

    cout << scientific;

    cout << "A = "
         << A << endl;

    cout << "\nModel:\n";

    cout << "L_R(lambda) = A * lambda^n\n";

    cout << "====================================================\n";


    return true;
}


// ============================================================
// ABSORPTION POWER-LAW FIT
// ============================================================
//
// L_abs(lambda) = A_A * lambda^n_A
//
// ln(L_abs) = ln(A_A) + n_A * ln(lambda)
//
// Saturated values are excluded.
// ============================================================

bool FitAbsorption(
    const vector<C2F6Data>& data,
    double& A,
    double& n)
{
    double Sx  = 0.0;
    double Sy  = 0.0;
    double Sxx = 0.0;
    double Sxy = 0.0;

    int N = 0;


    // --------------------------------------------------------
    // Use all valid non-saturated absorption measurements
    // --------------------------------------------------------

    for (const auto& p : data)
    {
        if (p.absorption <= 0.0)
            continue;

        if (p.absorption >= ABSORPTION_SATURATION)
            continue;


        double x = log(p.wavelength);
        double y = log(p.absorption);


        Sx  += x;
        Sy  += y;
        Sxx += x*x;
        Sxy += x*y;

        ++N;
    }


    if (N < 3)
    {
        cerr << "\nERROR: Not enough valid points for "
             << "absorption fit.\n";

        cerr << "Usable absorption points = "
             << N << endl;

        A = 0.0;
        n = 0.0;

        return false;
    }


    double denominator =
        N*Sxx - Sx*Sx;


    if (fabs(denominator) < 1e-20)
    {
        cerr << "\nERROR: Singular absorption fit."
             << endl;

        A = 0.0;
        n = 0.0;

        return false;
    }


    n =
        (N*Sxy - Sx*Sy)
        / denominator;


    double intercept =
        (Sy - n*Sx) / N;


    A = exp(intercept);


    cout << "\n====================================================\n";
    cout << "ABSORPTION FIT\n";
    cout << "====================================================\n";

    cout << "Fit type = Power law\n";

    cout << "Number of usable points = "
         << N << endl;

    cout << fixed << setprecision(6);

    cout << "Power n = "
         << n << endl;

    cout << scientific;

    cout << "A = "
         << A << endl;

    cout << "\nModel:\n";

    cout << "L_abs(lambda) = A * lambda^n\n";

    cout << "====================================================\n";


    return true;
}


// ============================================================
// FIND EXACT MEASURED POINT
// ============================================================
//
// Returns true only if the wavelength actually exists in the
// original dataset.
//
// This is important because:
//
//   measured point -> use original value exactly
//
// rather than interpolating it.
//

bool FindMeasuredPoint(
    double wavelength,
    const vector<C2F6Data>& data,
    C2F6Data& point)
{
    const double tolerance = 1e-6;


    for (const auto& p : data)
    {
        if (fabs(p.wavelength - wavelength) < tolerance)
        {
            point = p;
            return true;
        }
    }


    return false;
}


// ============================================================
// LINEAR INTERPOLATION
// ============================================================
//
// Used only for wavelengths that are NOT exact measured points
// but lie between two measured wavelengths.
//
// This is particularly useful because the last measured point
// is 785.5 nm while the output grid is 1 nm.
//

double Interpolate(
    double wavelength,
    const vector<C2F6Data>& data,
    bool rayleigh)
{
    if (data.empty())
        return -1.0;


    if (wavelength < data.front().wavelength)
        return -1.0;


    if (wavelength > data.back().wavelength)
        return -1.0;


    for (size_t i = 0; i < data.size()-1; ++i)
    {
        double x1 = data[i].wavelength;
        double x2 = data[i+1].wavelength;


        if (wavelength >= x1 &&
            wavelength <= x2)
        {
            double y1 =
                rayleigh
                ? data[i].rayleigh
                : data[i].absorption;


            double y2 =
                rayleigh
                ? data[i+1].rayleigh
                : data[i+1].absorption;


            if (fabs(x2-x1) < 1e-12)
                return y1;


            double fraction =
                (wavelength-x1)
                / (x2-x1);


            return y1 +
                   fraction*(y2-y1);
        }
    }


    return -1.0;
}


// ============================================================
// GET OPTICAL PROPERTY
// ============================================================
//
// Priority:
//
// 1. Exact measured value
// 2. Interpolation between measured values
// 3. Extrapolation above measured range
//
// ============================================================

double GetRayleigh(
    double wavelength,
    const vector<C2F6Data>& data,
    double A_R,
    double n_R)
{
    C2F6Data p;


    // --------------------------------------------------------
    // Exact measured point
    // --------------------------------------------------------

    if (FindMeasuredPoint(wavelength, data, p))
        return p.rayleigh;


    // --------------------------------------------------------
    // Inside measured range -> interpolation
    // --------------------------------------------------------

    if (wavelength <= data.back().wavelength)
    {
        double value =
            Interpolate(
                wavelength,
                data,
                true);

        if (value > 0.0)
            return value;
    }


    // --------------------------------------------------------
    // Above measured range -> extrapolation
    // --------------------------------------------------------

    return A_R * pow(wavelength, n_R);
}


// ============================================================

double GetAbsorption(
    double wavelength,
    const vector<C2F6Data>& data,
    double A_A,
    double n_A)
{
    C2F6Data p;


    // --------------------------------------------------------
    // Exact measured point
    // --------------------------------------------------------

    if (FindMeasuredPoint(wavelength, data, p))
        return p.absorption;


    // --------------------------------------------------------
    // Inside measured range -> interpolation
    // --------------------------------------------------------

    if (wavelength <= data.back().wavelength)
    {
        double value =
            Interpolate(
                wavelength,
                data,
                false);

        if (value > 0.0)
            return value;
    }


    // --------------------------------------------------------
    // Above measured range -> extrapolation
    // --------------------------------------------------------

    return A_A * pow(wavelength, n_A);
}


// ============================================================
// MAIN FUNCTION
// ============================================================

void extrapolate_C2F6()
{
    // ========================================================
    // INPUT FILE
    // ========================================================

    const char* inputFile =
        "C2F6_Rayleigh_Absorption.csv";


    // ========================================================
    // READ DATA
    // ========================================================

    vector<C2F6Data> data =
        ReadCSV(inputFile);


    if (data.empty())
    {
        cerr << "\nERROR: No data found."
             << endl;

        return;
    }


    // ========================================================
    // FIT RAYLEIGH
    // ========================================================

    double A_R = 0.0;
    double n_R = 0.0;


    if (!FitRayleigh(
            data,
            A_R,
            n_R))
    {
        return;
    }


    // ========================================================
    // FIT ABSORPTION
    // ========================================================

    double A_A = 0.0;
    double n_A = 0.0;


    if (!FitAbsorption(
            data,
            A_A,
            n_A))
    {
        return;
    }


    // ========================================================
    // OUTPUT FILE
    // ========================================================

    const char* outputFile =
        "C2F6_OpticalProperties_1000nm.txt";


    ofstream output(outputFile);


    if (!output.is_open())
    {
        cerr << "\nERROR: Cannot create output file:\n"
             << outputFile << endl;

        return;
    }


    // ========================================================
    // HEADER
    // ========================================================

    output << "# C2F6 optical properties\n";

    output << "#\n";

    output << "# Format:\n";

    output << "# Energy*eV   RayleighLength*m   AbsorptionLength*m\n";

    output << "#\n";

    output << "# Original measured values are preserved exactly.\n";

    output << "# Interpolation is used only between measured points.\n";

    output << "# Extrapolation starts above "
           << data.back().wavelength
           << " nm.\n";

    output << "#\n";


    // ========================================================
    // GENERATE TABLE
    // ========================================================
    //
    // 1000 nm -> 190 nm
    //
    // Therefore photon energy increases monotonically.
    //
    // ========================================================

    cout << "\n";
    cout << "====================================================\n";
    cout << "GENERATING TABLE\n";
    cout << "====================================================\n";


    for (double wavelength = MAX_WAVELENGTH;
         wavelength >= MIN_WAVELENGTH - 1e-9;
         wavelength -= WAVELENGTH_STEP)
    {
        // ----------------------------------------------------
        // Photon energy
        // ----------------------------------------------------

        double energy =
            HC / wavelength;


        // ----------------------------------------------------
        // Rayleigh length
        // ----------------------------------------------------

        double rayleigh =
            GetRayleigh(
                wavelength,
                data,
                A_R,
                n_R);


        // ----------------------------------------------------
        // Absorption length
        // ----------------------------------------------------

        double absorption =
            GetAbsorption(
                wavelength,
                data,
                A_A,
                n_A);


        // ----------------------------------------------------
        // Check values
        // ----------------------------------------------------

        if (!isfinite(rayleigh) ||
            rayleigh <= 0.0)
        {
            cerr << "WARNING: Invalid Rayleigh value at "
                 << wavelength
                 << " nm"
                 << endl;

            continue;
        }


        if (!isfinite(absorption) ||
            absorption <= 0.0)
        {
            cerr << "WARNING: Invalid absorption value at "
                 << wavelength
                 << " nm"
                 << endl;

            continue;
        }


        // ----------------------------------------------------
        // Write Geant4-style table
        // ----------------------------------------------------

        output
            << fixed
            << setprecision(8)
            << energy
            << "*eV   "
            << rayleigh
            << "*m   "
            << absorption
            << "*m"
            << "\n";
    }


    output.close();


    // ========================================================
    // 1000 nm VALUES
    // ========================================================

    double rayleigh_1000 =
        A_R * pow(1000.0, n_R);


    double absorption_1000 =
        A_A * pow(1000.0, n_A);


    double energy_1000 =
        HC / 1000.0;


    cout << "\n";
    cout << "====================================================\n";
    cout << "EXTRAPOLATED VALUES AT 1000 nm\n";
    cout << "====================================================\n";

    cout << fixed << setprecision(6);

    cout << "Wavelength        = 1000 nm\n";

    cout << "Energy            = "
         << energy_1000
         << " eV\n";

    cout << "Rayleigh length   = "
         << rayleigh_1000
         << " m\n";

    cout << "Absorption length = "
         << absorption_1000
         << " m\n";

    cout << "====================================================\n";


    // ========================================================
    // PRINT LAST MEASURED VALUE
    // ========================================================

    const C2F6Data& last =
        data.back();


    cout << "\n";
    cout << "LAST MEASURED POINT\n";
    cout << "====================================================\n";

    cout << "Wavelength        = "
         << last.wavelength
         << " nm\n";

    cout << "Energy            = "
         << last.energy
         << " eV\n";

    cout << "Rayleigh length   = "
         << last.rayleigh
         << " m\n";

    cout << "Absorption length = "
         << last.absorption
         << " m\n";

    cout << "====================================================\n";


    cout << "\nOutput file:\n";
    cout << outputFile << "\n";


    // ========================================================
    // DIAGNOSTIC GRAPHS
    // ========================================================

    TGraph* gRayleigh =
        new TGraph();

    TGraph* gAbsorption =
        new TGraph();


    TGraph* gRayleighExtra =
        new TGraph();

    TGraph* gAbsorptionExtra =
        new TGraph();


    int ir  = 0;
    int ia  = 0;
    int ire = 0;
    int iae = 0;


    // --------------------------------------------------------
    // Original measured data
    // --------------------------------------------------------

    for (const auto& p : data)
    {
        gRayleigh->SetPoint(
            ir++,
            p.wavelength,
            p.rayleigh);


        gAbsorption->SetPoint(
            ia++,
            p.wavelength,
            p.absorption);
    }


    // --------------------------------------------------------
    // Extrapolated region
    // --------------------------------------------------------

    for (double wavelength =
             data.back().wavelength + 0.1;

         wavelength <= MAX_WAVELENGTH;

         wavelength += 1.0)
    {
        double rayleigh =
            A_R *
            pow(wavelength, n_R);


        double absorption =
            A_A *
            pow(wavelength, n_A);


        gRayleighExtra->SetPoint(
            ire++,
            wavelength,
            rayleigh);


        gAbsorptionExtra->SetPoint(
            iae++,
            wavelength,
            absorption);
    }


    // ========================================================
    // RAYLEIGH PLOT
    // ========================================================

    TCanvas* c1 =
        new TCanvas(
            "c1",
            "C2F6 Rayleigh Length",
            900,
            700);


    gRayleigh->SetTitle(
        "C2F6 Rayleigh Length;"
        "Wavelength (nm);"
        "Rayleigh Length (m)");


    gRayleigh->SetMarkerStyle(20);

    gRayleigh->Draw("AP");


    gRayleighExtra->SetLineStyle(2);

    gRayleighExtra->SetLineWidth(2);

    gRayleighExtra->Draw("L SAME");


    TLegend* leg1 =
        new TLegend(
            0.55,
            0.75,
            0.88,
            0.88);


    leg1->AddEntry(
        gRayleigh,
        "Measured",
        "p");


    leg1->AddEntry(
        gRayleighExtra,
        "Extrapolated",
        "l");


    leg1->Draw();


    c1->SetLogy();


    c1->SaveAs(
        "C2F6_Rayleigh_Extrapolation.png");


    // ========================================================
    // ABSORPTION PLOT
    // ========================================================

    TCanvas* c2 =
        new TCanvas(
            "c2",
            "C2F6 Absorption Length",
            900,
            700);


    gAbsorption->SetTitle(
        "C2F6 Absorption Length;"
        "Wavelength (nm);"
        "Absorption Length (m)");


    gAbsorption->SetMarkerStyle(20);

    gAbsorption->Draw("AP");


    gAbsorptionExtra->SetLineStyle(2);

    gAbsorptionExtra->SetLineWidth(2);

    gAbsorptionExtra->Draw("L SAME");


    TLegend* leg2 =
        new TLegend(
            0.55,
            0.75,
            0.88,
            0.88);


    leg2->AddEntry(
        gAbsorption,
        "Measured",
        "p");


    leg2->AddEntry(
        gAbsorptionExtra,
        "Extrapolated",
        "l");


    leg2->Draw();


    c2->SetLogy();


    c2->SaveAs(
        "C2F6_Absorption_Extrapolation.png");


    // ========================================================
    // DONE
    // ========================================================

    cout << "\n";
    cout << "Diagnostic plots created:\n";

    cout << "  C2F6_Rayleigh_Extrapolation.png\n";

    cout << "  C2F6_Absorption_Extrapolation.png\n";

    cout << "\nDONE.\n";
}
