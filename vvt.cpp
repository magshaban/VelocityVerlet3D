/*
 *	vtt.cpp
 *	V5.0, (Optmized VERSION).
 *	
 * 	
 *	Info: Implementation of Velocity Verlet Algorithm from Lennard Jonnes System.
 *
 * -To run: after the compiling, from cmd or the terminal enter <vtt.exe atoms.xyz>
 *      or: From Windows OS after compiling, you will find A Batch file Run_vtt.bat this 
 * 	     	will run the program directlly 
 *    		where input parameters are listed in the file "atoms.xyz", 
 *    		the first number in the file is the number of atoms or particles
 *		and each line is the initial position vectors xyz of each atom.
 *							  	
 *	you can un the produced trajectory file by using OVITTO MD program.
 * 
 *  Maged Shaban
 *	maged.shaban[at]mathmods.eu
 *  Gdansk University of Technology, February 2,2017.
 */

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <sstream>
using namespace std;

typedef struct
{
	string name;
	double r[3];	//the particle position vector.
	double f[3];	//force vector.
	double v[3];	//the paricle velocity vector.
	double a[3]; 	//accelaeration vector.
	double pa[3]; 	//previous acceleration vector.
	double pot;		//the potential
} particle;
	
int nthermo , step;
string config_file_name;
double h,number_of_steps ,mass , sigma , epsilon ; //the parameters , h is the time step.
double sigma2,sigma6 , sigma12 ,massinverse;
double  Utot , Ktot, Etot;

string trajectory_filename;
string thermo_filename;
ofstream thermofile;
ofstream trajectory_file;

const double acceleration_converter = 0.0096485356;	 
const double Energy_converter = 103.642670915; 

///variables to optimize the program speed. 
double MassInversPerAccelerationConvert; 
double half_h,half_h_square;

particle* p; 	//p will define the particles or atoms. 

int NumberOfParticles ;// number of particles;

/* Function decleration */
void checkCMDinput(int argc);
void readInputFile(string filename);
void printRunParameters ();
void readsystemFile (string filename);
void Intialization ();
void printsystem();
void open_thermofile();
void write_thermoFile ();
void close_thermoFile ();
double cal_vec_rij2(double* rij_v, int i, int j) ;
void intialForceAndPotential();
void LJ_ForcesAndEnergies();
void print_ForcesAndEnergies ();
void calculate_Acceleratoins();
void calculate_Positions ();
void store_PreviousAccelerations();
void calculate_ThermoParameters();
void calculate_Velocities();
void patchfile();  
void close_TrajectoryFile();
void write_TrajectoryFile();
void open_TrajectoryFile();
void Sumulate();
void freeMemory();

int main(int argc, char* argv[]) 
{
    patchfile();  
    
	checkCMDinput(argc);
    string input_file_name = argv[1]; 
	
	readInputFile (input_file_name);
    printRunParameters();
    
    Intialization ();
    
    readsystemFile (config_file_name);
	printsystem();
	open_thermofile();
	open_TrajectoryFile();
	
	Sumulate();

	close_thermoFile();
    close_TrajectoryFile();
    
    cout<<endl<<" >> DONE!"<<endl;
			
	freeMemory();		
	return 0 ; 
}

void checkCMDinput(int argc)
{
	if (argc!=2)
	{
		cout << endl<< "> ERROR : Wrong Invocation!" << endl;
		cout <<endl<< "Please try with:" ;cout << endl<<" >>vvt.exe <string input_filename.txt>" << endl;
		exit(1);
	}
}

void readInputFile(string filename)
{
	ifstream in;
	string line;
	istringstream iss;
	string keyword;
	
	in.open(filename.c_str());
	if ( in.is_open() == 0 )
	{
		cout << " > Error: could not open file \"" << filename << "\"!" << endl;
		exit(2);
	}
	
	while ( getline(in, line) )
	{	
		iss.str(line);
		iss >> keyword;
		
		  if ( keyword == "epsilon" )
			iss >> epsilon;								
		else if ( keyword == "atomic_sys" )
			iss >> config_file_name;
		else if ( keyword == "sigma" )
			iss >> sigma;
		else if ( keyword == "mass" )
			iss >> mass;
		else if ( keyword == "trajectory" )
			iss >> trajectory_filename;
		else if ( keyword == "timestep" )
			iss >> h;
		else if ( keyword == "nsteps" )
			iss >> number_of_steps;
		else if ( keyword == "thermo" )
			iss >>  nthermo;	
		else if ( keyword == "log" )
			iss >> thermo_filename;
		
		iss.clear();
	}
	
	in.close();
}

void printRunParameters () 
{
    cout << "Run Parameters:" << endl;
    cout << "==============" << endl<<endl;
    cout << "atomic system       " << config_file_name << endl;
    cout << "epsilon             " << epsilon << endl;
    cout << "sigma               " << sigma << endl;
    cout << "mass                " << mass << endl;
    cout << "timestep            " << h << endl;
    cout << "nsteps              " << number_of_steps <<endl;
    cout << "thermo              " << nthermo <<endl;
    cout << "thermo_filename     " << thermo_filename << endl<<endl;
    cout << "===============================" <<endl<<endl;
}

void Intialization ()
{
	massinverse = 1.0 / mass;
    sigma2 = sigma * sigma;
    sigma6 = sigma2 * sigma2 * sigma2 ;
    sigma12 = sigma6 * sigma6;
    MassInversPerAccelerationConvert = massinverse * acceleration_converter;
    half_h = 0.5 * h;
    half_h_square = 0.5 * h * h;
}

void readsystemFile (string filename)
{
    ifstream myfile;
    string line;
    istringstream iss;
    
    myfile.open(filename.c_str());
    
    if (myfile.is_open() == 0) 
	{
        cout << " > Error: could not locate atomic system file" << endl;
        exit(3);
    }
    
    getline(myfile, line);
    iss.clear();
    iss.str(line);
    iss >> NumberOfParticles; 
    
	p = new particle[NumberOfParticles];
    
    getline(myfile, line);
    
    for (int i = 0; i < NumberOfParticles; i++) 
	{
        getline(myfile , line);
        iss.clear();
        iss.str(line);
        iss >> p[i].name >>p[i].r[0] >> p[i].r[1]>> p[i].r[2]; 		 // input the intial position 
        
		p[i].v[0] = 0.0; p[i].v[1] = 0.0; p[i].v[2] = 0.0; //zero intial velocities
    }
    
    myfile.close();
}

void printsystem()
{	
	cout<<"Total Number Of Atoms = "<< NumberOfParticles <<endl<<endl;

	for (int i = 0; i < NumberOfParticles ; i++)
	{
    cout<< "r = ["<<p[i].r[0] << ", "<<p[i].r[1]<<", " << p[i].r[2]<<"]"<<endl;
	}	
    cout<<endl<< "===============================" <<endl<<endl;
    cout<<" > Please wait, The simulation is running..."<<endl<<endl;
}

void open_thermofile() 
{
    thermofile.open(thermo_filename.c_str());
    
    if (thermofile.is_open() == 0) {
        cout << " > Error: Could Not Open Thermo File. " << endl;
        exit(4);
    }
    
	thermofile << "#timestep,			 #V,				#K.E,	 			#Total Energy" << endl;
	thermofile << "#=========================================================================" << endl;
}

void write_thermoFile() 
{
    thermofile << step  << "			 " << Utot << "		" <<  Ktot<< "		" << Etot << endl;
}

void close_thermoFile() 
{
    thermofile.close();
}

double cal_vec_rij2(double *rij_v, int i, int j) 
{
	double  r_ij2;
	
    rij_v[0] = p[j].r[0] - p[i].r[0];
    rij_v[1] = p[j].r[1] - p[i].r[1];
    rij_v[2] = p[j].r[2] - p[i].r[2];
    
    r_ij2 = rij_v[0] * rij_v[0] + rij_v[1] * rij_v[1] + rij_v[2] * rij_v[2];
    
    return r_ij2 ;
}

void intialForceAndPotential()
 {
    
    for (int i = 0 ; i < NumberOfParticles ; i++)
	{
        p[i].f[0] = 0.0;
        p[i].f[1] = 0.0;
        p[i].f[2] = 0.0;
        p[i].pot  = 0.0;
    }
        
}

void LJ_ForcesAndEnergies()  //calculateForcesAndEnergies
{
    double rij_v[3]; // r_ij vector ;
    double r_ij2;	//the distance squared 
    double inv_rij2;	
    double r6i; 
    double v_ij, f_ij;
    double fij_v[3];
    
    intialForceAndPotential();
 
    for (int i = 0; i < NumberOfParticles ; i++)
	{
    	for (int j=0; j < NumberOfParticles; j++) 
		{
  		if(i==j)
  			continue ;
  		
        r_ij2 = cal_vec_rij2(rij_v , i , j) ;
         
        inv_rij2 = 1.0 / r_ij2;
        r6i = inv_rij2 * inv_rij2 * inv_rij2;
   
        v_ij = 4*epsilon * r6i * (sigma12 * r6i - sigma6);
        f_ij = - 24 * epsilon * inv_rij2 * r6i * (2 * sigma12 * r6i - sigma6) ;
            
        fij_v[0] = f_ij * rij_v[0];
        fij_v[1] = f_ij * rij_v[1];
        fij_v[2] = f_ij * rij_v[2];
            
        p[i].f[0] += fij_v[0];
        p[i].f[1] += fij_v[1];
        p[i].f[2] += fij_v[2];

        p[i].pot += 0.5 * v_ij;
        }
    }
}

void print_ForcesAndEnergies () 
{
    for (int i = 0; i < NumberOfParticles ; i++ )
	{
    cout << "Atom " << i+1 << " f = [" << p[i].f[0] << ", " << p[i].f[1] << ", " << p[i].f[2] << "] "<< " v = " << p[i].pot << endl;
    }
}

void calculate_Acceleratoins() 
{
    
    for (int i = 0; i < NumberOfParticles ; i++)
	{
        p[i].a[0] = p[i].f[0] * MassInversPerAccelerationConvert;
        p[i].a[1] = p[i].f[1] * MassInversPerAccelerationConvert;
        p[i].a[2] = p[i].f[2] * MassInversPerAccelerationConvert;
    }
}

void calculate_Positions () 
{
    for (int i = 0 ; i < NumberOfParticles; i++) 
	{
        p[i].r[0] = p[i].r[0] + p[i].v[0] * h + half_h_square * p[i].a[0];
        p[i].r[1] = p[i].r[1] + p[i].v[1] * h + half_h_square * p[i].a[1];
        p[i].r[2] = p[i].r[2] + p[i].v[2] * h + half_h_square * p[i].a[2];
    }
}

void store_PreviousAccelerations() 
{
    for (int i = 0 ; i < NumberOfParticles ; i++)
	{
        p[i].pa[0] = p[i].a[0];
        p[i].pa[1] = p[i].a[1];
        p[i].pa[2] = p[i].a[2];
    }
}

void calculate_Velocities() 
{
    for (int i = 0; i < NumberOfParticles; i++) 
	{
        p[i].v[0] = p[i].v[0] +  half_h * (p[i].pa[0] + p[i].a[0]);
        p[i].v[1] = p[i].v[1] +  half_h * (p[i].pa[1] + p[i].a[1]);
        p[i].v[2] = p[i].v[2] +  half_h * (p[i].pa[2] + p[i].a[2]);
    }
}
 
void calculate_ThermoParameters() 
{
    double v2 ; 
    
    Ktot = 0;  Utot = 0; Etot = 0;
    
    for ( int i = 0 ; i < NumberOfParticles ; i++)
	{
        v2 = p[i].v[0] * p[i].v[0] + p[i].v[1] * p[i].v[1] + p[i].v[2] * p[i].v[2];
        
        Ktot = Ktot + v2 ;
        
        Utot = Utot + p[i].pot;
    }
    
    Etot = Ktot * 0.5 * mass * Energy_converter + Utot;
}

void patchfile()
{
	ofstream patch_file;
	
	patch_file.open("Run_vvt.bat");
			
    patch_file<< "@ECHO OFF "<<endl
	   	<<"cls"<<endl
		<<"title  vtt V5.0  Maged Shaban, "<<endl
		<<"color 1f"<<endl
		<<"ECHO. "<<endl 
		<<"cmd.exe /K \"vvt.exe input.txt  &&  vvt.exe input.txt %> out_vtt.dat \""<<endl
		<<"pause;"<<endl;
			
	patch_file.close();
}

void open_TrajectoryFile() 
{
    trajectory_file.open(trajectory_filename.c_str());
    
    if (trajectory_file.is_open() == 0) 
	{
        cout << " >> Error: Could not create trajectory file!" << endl;
        exit(5);
    }
    
    
}

void write_TrajectoryFile() 
{    
    trajectory_file << NumberOfParticles << endl;
    trajectory_file << "Timestep " << step << endl;

    
    for (int i = 0; i < NumberOfParticles; i++) 
	{
        trajectory_file	<< p[i].name <<" " << p[i].r[0] << " " << p[i].r[1] << " " << p[i].r[2] << endl;
    }
  

    trajectory_file.flush();    
    
}

void close_TrajectoryFile() 
{
    trajectory_file.close();
}

void Sumulate()
{
	LJ_ForcesAndEnergies();
    calculate_Acceleratoins();
	for (step = 1 ; step <= number_of_steps ; step ++) 
	{ 
    	calculate_Positions();
    	LJ_ForcesAndEnergies();
    	store_PreviousAccelerations();
    	calculate_Acceleratoins();
    	calculate_Velocities();
        calculate_ThermoParameters();
		write_thermoFile();  
        write_TrajectoryFile();
	}
}

void freeMemory () 
{
    delete [] p;
}

