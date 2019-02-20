
#include <iostream>
#include <regex>
#include "AMReX_DataServices.H"
#include <Radiation_F.H>
#include <Radiation.H>

using namespace amrex;

std::string inputs_name = "";

void Print_Help() {
	Print() << "\nPrint out the radiation quantities at a specified distance from"
	        << "\nthe origin.  This is written for the 1-d radiating sphere problem."
	        << "\n"
	        << "\n./fradsphere -p plotfile -r radius -g groupfile"
	        << "\n"
	        << "\nHere groupfile is the file containing the group structure information"
	        << "\nas output by Castro (usually group_structure.dat)."
	        << "\n\n" << std::endl;
}

int main(int argc, char* argv[])
{

	amrex::Initialize(argc, argv, false);

	// timer for profiling
	BL_PROFILE_VAR("main()", pmain);

	{
		// Input arguments
		string pltfile, groupfile;
		Real radius = 0.;
		int j = 1; // skip program name

		while ( j < argc)
		{

			if ( !strcmp(argv[j], "-p") || !strcmp(argv[j],"--pltfile") )
			{
				pltfile = argv[++j];
			}
			else if ( !strcmp(argv[j], "-g") || !strcmp(argv[j],"--groupfile") )
			{
				groupfile = argv[++j];
			}
			else if ( !strcmp(argv[j], "-r") || !strcmp(argv[j],"--radius") )
			{
				radius = std::atof(argv[++j]);
			}
			else
			{
				std::cout << "\n\nOption " << argv[j] << " not recognized" << std::endl;
				Print_Help();
				exit ( EXIT_FAILURE );
			}

			// Go to the next parameter name
			++j;
		}

		if (pltfile.empty() || groupfile.empty())
		{
			Print_Help();
			Abort("Missing input file");
		}

		Print() << "\nplotfile  = \"" << pltfile << "\"" << std::endl;
		Print() <<   "groupfile = \"" << groupfile << "\"" << std::endl;
		Print() << std::endl;

		// Start dataservices (no clue why we need to do this)
		DataServices::SetBatchMode();

		// Define the type of file
		Amrvis::FileType fileType(Amrvis::NEWPLT);
		DataServices dataServices (pltfile, fileType);

		if (!dataServices.AmrDataOk())
			DataServices::Dispatch(DataServices::ExitRequest, NULL);

		// get data from plot file
		AmrData& data = dataServices.AmrDataRef();

		int finestLevel = data.FinestLevel();

		// get variable names
		const Vector<string>& varNames = data.PlotVarNames();

		// get the index bounds and dx.
		Box domain = data.ProbDomain()[finestLevel];

		Vector<Real> dx(AMREX_SPACEDIM);
		for (int i = 0; i < AMREX_SPACEDIM; i++)
			dx[i] = data.ProbSize()[i] / domain.length(i);

		const Vector<Real>& problo = data.ProbLo();
		const Vector<Real>& probhi = data.ProbHi();
		Vector<int> rr = data.RefRatio();

		if (radius < problo[0] || radius > probhi[0])
			Abort("ERROR: specified observer radius outside of domain");

		std::cout.setf(std::ios::scientific);
		std::cout.precision(12);
		std::cout << "rmin = " << problo[0] << std::endl;
		std::cout << "rmax = " << probhi[0] << std::endl << std::endl;

		int nbins = domain.length(0);

		// find variable indices
		Vector<std::string> compVarNames = {"rad0"};

		auto varComps = GetComponents(data, compVarNames);
		auto rad_comp = varComps[0];

		const int nvars = data.NComp();

		// allocate storage for data
		Vector<Real> vars_bin(nbins * (nvars + 1), 0.);
		auto r1 = 1.0;

		// fill a multifab with the data
		Vector<int> fill_comps(data.NComp());
		for (auto i = 0; i < data.NComp(); i++)
			fill_comps[i] = i;

		// ! imask will be set to false if we've already output the data.
		// ! Note, imask is defined in terms of the finest level.  As we loop
		// ! over levels, we will compare to the finest level index space to
		// ! determine if we've already output here
		int mask_size = nbins;
		for (auto i = 0; i < finestLevel - 1; i++)
			mask_size *= rr[i];


#if (AMREX_SPACEDIM == 1)
		Vector<int> imask(mask_size);
#elif (AMREX_SPACEDIM >=2)
		Vector<int> imask(pow(mask_size, AMREX_SPACEDIM));
#endif

		for (auto it=imask.begin(); it!=imask.end(); ++it)
			*it = 1;

		int cnt = 0;

		// extract the 1d data
		for (int l = finestLevel; l >= 0; l--) {

			Vector<Real> level_dx = data.CellSize(l);

			const BoxArray& ba = data.boxArray(l);
			const DistributionMapping& dm = data.DistributionMap(l);

			MultiFab lev_data_mf(ba, dm, data.NComp(), data.NGrow());
			data.FillVar(lev_data_mf, l, varNames, fill_comps);

			for (MFIter mfi(lev_data_mf, true); mfi.isValid(); ++mfi) {
				const Box& bx = mfi.tilebox();

				fradsphere(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
				           ZFILL(problo), ZFILL(probhi),
				           BL_TO_FORTRAN_FAB(lev_data_mf[mfi]),
				           nbins, vars_bin.dataPtr(),
				           imask.dataPtr(), mask_size, r1,
				           ZFILL(dx), &cnt);
			}

			// adjust r1 for the next lowest level
			if (l != 0) r1 *= rr[l-1];
		}

		// sort the data based on the coordinates

		Vector<Real> coords(cnt);
		for (auto i = 0; i < cnt; i++)
			coords[i] = vars_bin[i];

		auto isv = sort_indexes(coords);

		// open the group file and read in the group information
		std::ifstream group_file;
		group_file.open(groupfile);

		string header_line;

		string line;

		Vector<Real> nu_groups;
		Vector<Real> dnu_groups;

		bool first_line = true;

		int ngroups = 0;

		while(std::getline(group_file, line)) {
			if (first_line) {

				const std::regex re("=\\s*([0-9]*)");
				std::smatch m;
				std::regex_search(line, m, re);

				ngroups = stoi(m[1].str());

				first_line = false;

				// skip next line
				std::getline(group_file, line);

			} else {

				// read in the group centers and weights
				Real nu, dnu;
				std::istringstream iss(line);
				iss >> nu >> dnu;

				nu_groups.push_back(nu);
				dnu_groups.push_back(dnu);
			}
		}

		// find the index corresponding to the desired observer radius
		auto idx_obs = -1;

		for (auto i = 0; i < cnt; i++) {
			if (radius >= vars_bin[isv[i]*(data.NComp()+1)] && radius < vars_bin[isv[i+1]*(data.NComp()+1)])
				idx_obs = i;
			break;
		}

		if (idx_obs == -1) Abort("ERROR: radius not found in domain");

		// output all the radiation energies

		const auto w = 28;

		Print() << std::setw(15) << "group name"
		        << std::setw(w) << "group center energy"
		        << std::setw(w) << "E_rad(nu)*dnu (erg/cm^3)"
		        << std::setw(w) << "E_rad(nu) (erg/cm^3/Hz)" << std::endl;

		for (auto i = 0; i < ngroups; i++) {
			Print() << std::setw(15) << varNames[rad_comp+i]
			        << std::setw(w) << nu_groups[i]
			        << std::setw(w) << vars_bin[isv[idx_obs] + (rad_comp+i+1)*nbins]
			        << std::setw(w) << vars_bin[isv[idx_obs] + (rad_comp+i+1)*nbins] / dnu_groups[i] << std::endl;
		}

	}

	// destroy timer for profiling
	BL_PROFILE_VAR_STOP(pmain);

	amrex::Finalize();
}
