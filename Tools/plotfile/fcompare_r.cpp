// SPDX-FileCopyrightText: 2023 - 2025 Yadong Zeng<zdsjtu@gmail.com> & ZhuXu Li<1246206018@qq.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include <AMReX.H>
#include <AMReX_Geometry.H>
#include <AMReX_Print.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_MultiFabUtil.H>
#include <algorithm>
#include <limits>
#include <cmath>
#include <cstdlib>

using namespace amrex;

struct ErrZone {
    Real max_abs_err = std::numeric_limits<Real>::lowest();
    int level;
    int grid_index;
    IntVect cell;
};

void PrintUsage()
{
    amrex::Print()
        << "\n"
        << " Compare two plotfiles, zone by zone, to machine precision\n"
        << " and report the maximum absolute and relative errors for each\n"
        << " variable.\n"
        << "\n"
        << " usage:\n"
        << "    fcompare [-n|--norm num] [-d|--diffvar var] [-z|--zone_info var] [-a|--allow_diff_grids] [-l|--allow_diff_num_levels] [-r|rel_tol] [--abs_tol] [--abort_if_not_all_found] file1 file2\n"
        << "\n"
        << " optional arguments:\n"
        << "    --infile1                : plotfile 1"
        << "    --infile2                : plotfile 2"
        << "    --l1to2                  : plotfile 1 level compare to plotfile2"
        << "    -n|--norm num            : what norm to use (default is 0 for inf norm)\n"
        << "    --abort_if_not_all_found : abort if not all variables are present in both files\n"
        << '\n';
}

int main_main()
{
    const int narg = amrex::command_argument_count();

    Real global_error = 0.0;
    bool any_nans = false;
    ErrZone err_zone;
    bool all_variables_found = true;
    bool all_variables_passed = true;
    int plot1_to_plot2_level = 0;

    // defaults
    int norm = 0;
    std::string plotfile_a;
    std::string plotfile_b;
    std::string diffvar;
    int zone_info = false;
    std::string zone_info_var_name;
    Vector<std::string> plot_names(1);
    bool abort_if_not_all_found = false;

    int farg = 1;
    while (farg <= narg) {
        const std::string fname = amrex::get_command_argument(farg);
        if (fname == "-h" || fname == "--help"){
            PrintUsage();
            return EXIT_SUCCESS;
        } else if (fname == "--infile1") {
            plotfile_a = amrex::get_command_argument(++farg);
        } else if (fname == "--infile2") {
            plotfile_b = amrex::get_command_argument(++farg);
        } else if (fname == "-l1to2"){
            plot1_to_plot2_level = std::stoi(amrex::get_command_argument(++farg));
        } else if (fname == "-n" || fname == "--norm") {
            norm = std::stoi(amrex::get_command_argument(++farg));
        } else if (fname == "-z" || fname == "--zone_info") {
            zone_info_var_name = amrex::get_command_argument(++farg);
            zone_info = true;
        } else if (fname == "-d" || fname == "--diffvar") {
            diffvar = amrex::get_command_argument(++farg);
            plot_names[0] = diffvar;
        } else if (fname == "--abort_if_not_all_found") {
            abort_if_not_all_found = true;
        } else {
            break;
        }
        ++farg;
    };

    if (plotfile_a.empty() || plotfile_b.empty()) {
        PrintUsage();
        return EXIT_FAILURE;
    }

    PlotFileData pf_a(plotfile_a);
    PlotFileData pf_b(plotfile_b);
    pf_b.syncDistributionMap(pf_a);

    const int dm = pf_a.spaceDim();
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(pf_a.spaceDim() == pf_b.spaceDim(),
                                     "ERROR: plotfiles have different numbers of spatial dimensions");

    const int finest_level = std::min(pf_a.finestLevel(), pf_b.finestLevel());
    const int nlevels = finest_level+1;

    const int ncomp_a = pf_a.nComp();
    const int ncomp_b = pf_b.nComp();

    if (ncomp_a != ncomp_b) {
        amrex::Print() << "\n WARNING: number of variables do not match\n";
    }

    int save_var_a = -1;
    int zone_info_var_a = -1;

    // get fab component name
    const Vector<std::string>& names_a = pf_a.varNames();
    const Vector<std::string>& names_b = pf_b.varNames();

    // compare the components are same order
    Vector<int> ivar_b(ncomp_a,-1); // in case the variables are not in the same order
    for (int n_a = 0; n_a < ncomp_a; ++n_a) {
        auto r = std::find(std::begin(names_b), std::end(names_b), names_a[n_a]);
        if (r == std::end(names_b)) {
            amrex::Print() << " WARNING: variable " << names_a[n_a] << " not found in plotfile 2\n";
            all_variables_found = false;
        } else {
            ivar_b[n_a] = static_cast<int>(std::distance(std::begin(names_b), r));
        }

        if (names_a[n_a] == diffvar) {
            save_var_a = n_a;
        }

        if (names_a[n_a] == zone_info_var_name) {
            zone_info_var_a = n_a;
        }
    }

    // also print out, as a diagnostic, those variables in plotfile 1 that
    // are not in plotfile 2
    for (int n_b = 0; n_b < ncomp_b; ++n_b) {
        auto r = std::find(std::begin(names_a),std::end(names_a),names_b[n_b]);
        if (r == std::end(names_a)) {
            amrex::Print() << " WARNING: variable " << names_b[n_b] << " not found in plotfile 1\n";
            all_variables_found = false;
        }
    }

    // create a multifab to store the difference for output, if desired
    Vector<MultiFab> mf_array(nlevels);
    if (save_var_a >= 0) {
        for (int ilev = 0; ilev < nlevels; ++ilev) {
            mf_array[ilev].define(pf_a.boxArray(ilev),
                                  pf_a.DistributionMap(ilev),
                                  1, 0);
        }
    }

    amrex::Print() << "\n"
                   << "  " << std::setw(24) << std::right << "variable name"
                   << "  " << std::setw(24) << "absolute error"
                   << "  " << std::setw(24) << "relative error" << "\n"
                   << "  " << std::setw(24) << " "
                   << "  " << std::setw(24) << "(||A - B||)"
                   << "  " << std::setw(24) << "(||A - B||/||A||)" << "\n"
                   << "  " << std::string(76,'-') << "\n";

    // go level-by-level and patch-by-patch and compare the data
    // only coarse level

    #define ilev 0
    if (pf_a.boxArray(ilev).empty() && pf_b.boxArray(ilev).empty()) {
        amrex::Abort("do not have any fab in plotfile");
    }

    Vector<Real> aerror(ncomp_a, 0.0);
    Vector<Real> rerror(ncomp_a, 0.0);
    Vector<Real> rerror_denom(ncomp_a, 0.0);
    Vector<int> has_nan_a(ncomp_a, false);
    Vector<int> has_nan_b(ncomp_a, false);
    for (int icomp_a = 0; icomp_a < ncomp_a; ++icomp_a) {
        if (ivar_b[icomp_a] >= 0) {
            // get mf_d from pf_a
            const MultiFab& mf_d = pf_a.get(ilev, names_a[icomp_a]);
            // get mf_b from pf_b
            MultiFab mf_b = pf_b.get(ilev, names_b[ivar_b[icomp_a]]);
            // down mf_d to mf_a
            MultiFab mf_a;
            Geometry geom_a(pf_a.probDomain(0), RealBox(pf_a.probLo(), pf_a.probHi()),pf_a.coordSys(), {0,0,0});
            Geometry geom_b(pf_b.probDomain(0), RealBox(pf_b.probLo(), pf_b.probHi()),pf_b.coordSys(), {0,0,0});
            amrex::average_down(mf_d, mf_a, geom_a, geom_b, 0, 1, IntVect(plot1_to_plot2_level));

            has_nan_a[icomp_a] = mf_a.contains_nan();
            has_nan_b[icomp_a] = mf_b.contains_nan();
            MultiFab::Subtract(mf_b,mf_a,0,0,1,0); // b = b - a
            Real max_err = mf_b.norm0();
            if (norm == 1) {
                aerror[icomp_a] = mf_b.norm1();
                rerror[icomp_a] = aerror[icomp_a];
                rerror_denom[icomp_a] = mf_a.norm1();
            } else if (norm == 2) {
                aerror[icomp_a] = mf_b.norm2();
                rerror[icomp_a] = aerror[icomp_a];
                rerror_denom[icomp_a] = mf_a.norm2();
            } else {
                aerror[icomp_a] = max_err;
                rerror[icomp_a] = aerror[icomp_a];
                rerror_denom[icomp_a] = mf_a.norm0();
            }

            if (norm == 0) {
                rerror[icomp_a] /= rerror_denom[icomp_a];
            } else {
                const auto& dx = pf_a.cellSize(ilev);
                Real dv = 1.0;
                for (int idim = 0; idim < dm; ++idim) {
                    dv *= dx[idim];
                }
                aerror[icomp_a] *= std::pow(dv,Real(1.)/static_cast<Real>(norm));
                rerror[icomp_a] = rerror[icomp_a]/rerror_denom[icomp_a];
            }

            if (icomp_a == save_var_a || icomp_a == zone_info_var_a) {
                mf_b.abs(0,1);
            }

            if (icomp_a == save_var_a) {
                MultiFab::Copy(mf_array[ilev], mf_b, 0, 0, 1, 0);
            }

            if (icomp_a == zone_info_var_a) {
                if (max_err > err_zone.max_abs_err) {
                    err_zone.max_abs_err = max_err;
                    err_zone.level = ilev;
                    err_zone.cell = mf_b.maxIndex(0);
                    auto isects = pf_a.boxArray(ilev).intersections
                        (Box(err_zone.cell,err_zone.cell), true, 0);
                    err_zone.grid_index = isects[0].first;
                }
            }
        }
    }

    amrex::Print() << " level = " << ilev << "\n";
    for (int icomp_a = 0; icomp_a < ncomp_a; ++icomp_a) {
        if (ivar_b[icomp_a] < 0) {
            amrex::Print() << " " << std::setw(24) << std::left << names_a[icomp_a]
                            << "  " << std::setw(50)
                            << "< variable not present in both files > "
                            << "\n";
        } else if (has_nan_a[icomp_a] && has_nan_b[icomp_a]) {
            amrex::Print() << " " << std::setw(24) << std::left << names_a[icomp_a]
                            << "  " << std::setw(50)
                            << "< NaN present in both A and B > "
                            << "\n";
        } else if (has_nan_a[icomp_a]) {
            amrex::Print() << " " << std::setw(24) << std::left << names_a[icomp_a]
                            << "  " << std::setw(50)
                            << "< NaN present in A > "
                            << "\n";
        } else if (has_nan_b[icomp_a]) {
            amrex::Print() << " " << std::setw(24) << std::left << names_b[icomp_a]
                            << "  " << std::setw(50)
                            << "< NaN present in B > "
                            << "\n";
        } else {
            Real aerr = 0., rerr = 0.;
            if (aerror[icomp_a] > 0.) {
                aerr = std::min(
                    std::max(aerror[icomp_a], std::numeric_limits<Real>::min()),
                    std::numeric_limits<Real>::max());
            }
            if (rerror[icomp_a] > 0.) {
                rerr = std::min(
                    std::max(rerror[icomp_a], std::numeric_limits<Real>::min()),
                    std::numeric_limits<Real>::max());
            }
            amrex::Print() << " " << std::setw(24) << std::left << names_a[icomp_a]
                            << std::right
                            << "  " << std::setw(24) << std::setprecision(10) << aerr
                            << "  " << std::setw(24) << std::setprecision(10) << rerr
                            << "\n";
        }
    }

    global_error = std::max(global_error,
                            *(std::max_element(aerror.begin(),
                                                aerror.end())));

    for (int icomp_a = 0; icomp_a < ncomp_a; ++icomp_a) {
        any_nans = any_nans || has_nan_a[icomp_a] || has_nan_b[icomp_a];
    }

    if (zone_info) {
        if (err_zone.max_abs_err > 0.) {
            ParallelDescriptor::Barrier();
            const DistributionMapping& dmap = pf_a.DistributionMap(err_zone.level);
            bool owner_proc = ParallelDescriptor::MyProc() == dmap[err_zone.grid_index];

            if (owner_proc) {
                amrex::AllPrint() << '\n'
                                  << " maximum error in " << zone_info_var_name << "\n"
                                  << "   level = " << err_zone.level << " (i,j,k) = " << err_zone.cell << "\n";
            }

            for (int icomp_a = 0; icomp_a < ncomp_a; ++icomp_a) {
                const MultiFab& mf = pf_a.get(err_zone.level,names_a[icomp_a]);
                if (owner_proc) {
                    Real v = mf[err_zone.grid_index](err_zone.cell);
                    amrex::AllPrint() << " " << std::setw(24)
                                      << names_a[icomp_a] << "  "
                                      << std::setw(24) << std::right
                                      << v << "\n";
                }
            }
        }
    }

    if (! all_variables_found) {
        amrex::Print() << " WARNING: not all variables present in both files\n";
        if (abort_if_not_all_found) { return EXIT_FAILURE; }
    }

    if (any_nans) {
        return EXIT_FAILURE;
    } else if (global_error == 0.0) {
        amrex::Print() << " PLOTFILE AGREE" << '\n';
        return EXIT_SUCCESS;
    } else {
        return EXIT_FAILURE;
    }
}

int main (int argc, char* argv[])
{
    amrex::SetVerbose(0);
    amrex::Initialize(argc, argv, false);
    int r = main_main();
    amrex::Finalize();
    return r;
}
