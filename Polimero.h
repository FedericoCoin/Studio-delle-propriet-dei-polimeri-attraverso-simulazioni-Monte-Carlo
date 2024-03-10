#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <random>
#include <string>
#include <algorithm> 
#include <fstream>
#include <string>
#include <chrono>

//Costanti
using namespace std;
const double pi = 3.14159265358979323846;
const double LJ_shift = 4 * (1 /pow(2.5,12) - 1 /(2.5,6));
const double cut=pow(2,1./6);
const double c1=7.5745810510344;
//Funzioni globali
inline double LJ(double x2) {
    double x6 = x2 * x2 * x2;
    double x12 = x6 * x6;
	if(x2<6.25)
	{
       return 4 * (1 / x12 - 1 / x6)-LJ_shift;
	}
    else
	{
		return 0;
	}
}
inline double Fene(double x, double k, double r0) {
	if(x<cut)
	{
     return -k * r0 * r0 * log(1 - (x * x) / (r0 * r0)) / 2+LJ(x*x)-c1;
	}
	else
	{
		return -k * r0 * r0 * log(1 - (x * x) / (r0 * r0))/2-c1-1;
	}
    
}
inline double gaussian(double r, double sigma)
{
	return exp(-(r * r) / (2 * sigma * sigma));
}
//----------------- CLASSE POSIZIONI

class Posizioni {
    public:
		double** r;
		Posizioni(int N)
		{
			r = new double* [N];
			for (int i = 0; i < N; i++)
			{
				r[i] = new double[3];
			}
	    }
};
std::mt19937 gen(0); // Use a fixed seed
std::uniform_real_distribution<> dis(0.0, 1.0);

// Function to generate a random double between -0.5 and 0.5
double random_half_range() {
    return dis(gen) - 0.5;
}
//----------------- CLASSE DEL SISTEMA 
class Polimero {
	public:
	Polimero(double* costanti_d, int* costanti_i)
	{
		//Dichiaro le costanti
		this->T = *(costanti_d);
		this->k = *(costanti_d + 1);
		this->r0 = *(costanti_d + 2);
		this->Delta[0] = *(costanti_d + 3);

		this->N = *costanti_i;
		this->N_eq = *(costanti_i + 1);
		this->N_stat = *(costanti_i + 2);
		this->N_bin = *(costanti_i + 3);
	
		N_move=N;
		Delta[1] = 0.05;
		Delta[2] = 0.05;
		Delta[3] = 0.2;
		Delta[4] = 0.02;
		dDelta[0] = Delta[0] * 0.5;
		dDelta[1] = Delta[1] * 0.5;
		dDelta[2] = Delta[2] * 0.5;
		dDelta[3] = Delta[3] * 0.5;
		dDelta[4] = Delta[4] * 0.5;
		dR = 5.0/ ((double)N_bin);
		dR2 = 1.5 / ((double)N_bin);
		EM = new double[N_stat];
		RGM = new double[N_stat];
		RMAXM = new double[N_stat];
		R2M = new double[N_stat];

		gr = new double[N_bin];
		gr2 = new double[N_bin];
		norm = new double[N_bin];
		r1 = new Posizioni(N);
		r2 = new Posizioni(N);
		temp = new Posizioni(N);
		backup = new Posizioni(N);

		neighbours = new int* [N];
		neighboursprec = new int* [N];
		for (int i = 0; i < N; i++)
		{
			neighbours[i] = new int[N];
			neighboursprec[i] = new int[N];
		}
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				neighbours[i][j] = -1;
			}
		}
	}
	void CreaSistema(bool load=false,string* positions=NULL)
	{
		if (load)
		{
			load_r(positions);
		}
		else
		{
			init_r();
		}
		find_neighbours(r1);
		calcE(r1);
		init_Gr();
		E = E_new;
		test[0][0] = 2;
		test[1][0] = 2;
		test[2][0] = 2;
		test[3][0] = 2;
		test[4][0] = 2;
	    stepposidx=0;
	}
	void print_condizioni_iniziali()
	{
		cout << " Simulazione Polimero\n" << endl;
		cout << " Temperatura: " << T << endl;
		cout << " Energia potenziale iniziale: " << E << endl;
		cout << " N: " << N << endl;
	}
	void fase_equilibratura()//
	{
		acc_tot[0] = 0;
		acc_tot[1] = 0;
		acc_tot[2] = 0;
		acc_tot[3] = 0;
		acc_tot[4] = 0;
		ofstream filepos;
		cout << "\nInizio fase di equilibratura\n" << endl;
		for (int i = 0; i < N_eq; i++)
		{		
			if (stepposidx==Step_pos)
			{
				find_neighbours(r1);
			}
			perform_moves();
			if (i % Step_control == 0 && i < N_eq / 2 && i != 0) {
				update_deltas();
			}
			if (i % Step_print == 0) {
				cout << "i: " << i << " U: " << E << " Rmax:  " << Rmax << " Rg: " << Rg << " R2: " << R2 << " " << (double)acc_tot[0]/(i)<<" "<<Delta[0]<< endl;
			}
			stepposidx++;
		}
		cout << "\nFine fase di equilibratura\n" << endl;
	}
	void fase_statistica(double *Gr,double *Dist,double *Gr2,double *Dist2)
	{
		E = E_new;
		fasestat = true;
		cout << "\nInizio fase di statistica\n" << endl;
		idx = 0;
		acc_tot[0] = 0;
		acc_tot[1] = 0;
		acc_tot[2] = 0;
		for (int i = 0; i < N_stat; i++)
		{
			if (i % Step_pos == 0)
			{
				find_neighbours(r1);
			}
			perform_moves();
			EM[i] = E;
			RMAXM[i] = Rmax;
			RGM[i] = Rg;
			R2M[i] = R2;
			Esum += E;
			E2sum += E * E;
			Rmaxsum += Rmax;
			Rgsum += Rg;
			R2sum += R2;
			Rmax2sum += Rmax * Rmax;
			Rg2sum += Rg * Rg;
			R22sum += R2 * R2;

			if (i % Step_print == 0) {
				cout << "i: " << i << " U: " << E << " Rmax:  " << Rmax << " Rg: " << Rg << " R2: " << R2 << endl;
			}
		}

		cout << "\nFine fase di statistica\n" << endl;
		norm_gr();
		for (int i = 0; i < N_bin; i++)
		{
			Dist[i] = ((double)i + 0.5) * dR;
			Dist2[i] = ((double)i + 0.5) * dR2;
			Gr[i]=gr[i];
			Gr2[i]=gr2[i];
		}
		
	}
	void calcola_osservabili(double* oss,string* filenamepos = NULL) {

		double Emean = Esum / N_stat;
		double Rgmean = Rgsum / N_stat;
		double Rmaxmean = Rmaxsum / N_stat;
		double R2mean = R2sum / N_stat;
		double Evar = E2sum / N_stat - pow(Emean, 2);
		double Estd = sqrt(Evar / N_stat);
		double Rgvar = Rg2sum / N_stat - pow(Rgmean, 2);
		double Rgstd = sqrt(Rgvar / N_stat);
		double R2var = R22sum / N_stat - pow(R2mean, 2);
		double R2std = sqrt(R2var / N_stat);
		double Rmaxvar = Rmax2sum / N_stat - pow(Rmaxmean, 2);
		double Rmaxstd = sqrt(Rmaxvar / N_stat);
        
		double tauE=calculateTau(EM,N_stat,25,Emean);
		double tauRmax=calculateTau(RMAXM,N_stat,25,Rmaxmean);
		double tauRg=calculateTau(RGM,N_stat,25,Rgmean);
		double tauR2=calculateTau(R2M,N_stat,25,R2mean);

		double Estdcorr =    sqrt(Evar*(1+tauE)/N_stat);//block(EM, N_stat);
		double Rmaxstdcorr = sqrt(Rmaxvar*(1+tauRmax)/N_stat);//block(RMAXM, N_stat);
        double Rgstdcorr =    sqrt(Rgvar*(1+tauRg)/N_stat);//block(RGM, N_stat);
		double R2stdcorr =   sqrt(R2var*(1+tauR2)/N_stat);//block(R2M, N_stat);
		

		cout << "\nEnergia potenziale: " << Emean << " +- " << Estd << " - " << Estdcorr <<" - "<<tauE<< endl;
		cout << "\nRaggio massimo: " << Rmaxmean << " +- " << Rmaxstd << " - " << Rmaxstdcorr <<" - "<<tauRmax<< endl;
		cout << "\nRaggio G: " << Rgmean << " +- " << Rgstd << " - " << Rgstdcorr <<" - "<<tauRg<< endl;
		cout << "\nRaggio G: " << R2mean << " +- " << R2std << " - " << R2stdcorr <<" - "<<tauR2<< endl;

		oss[0] = N;
		oss[1] = Emean;
		oss[2] = Estdcorr;
		oss[3] = Rmaxmean;
		oss[4] = Rmaxstdcorr;
		oss[5] = Rgmean;
		oss[6] = Rgstdcorr;
		oss[7] = R2mean;
		oss[8] = R2stdcorr;
		oss[9] = tauE;
		oss[10] = tauRmax;
		oss[11] = tauRg;
		oss[12] = tauR2;
		ofstream filepos;
		if (filenamepos != NULL)
		{
			filepos.open(*filenamepos);
		}
		for (int i = 0; i < N; i++)
			{
				filepos << r1->r[i][0] << ";" << r1->r[i][1] << ";" << r1->r[i][2]<< endl;
			}
		filepos.close();
	}
	private:
	int N,//Numero di monomeri
	 N_eq,//Numero di iterazioni per fase equilibratura
	  N_stat,//Numero di iterazioni per fase statistica
	 //N_cv,
	 N_bin,
	 N_move,
	 grcount = 0,
	Step_pos = 50,//Numero di iterazioni per la matrice dei primi vicini
	acc[5] = {0,0,0,0,0},//Numero di mosse accettate temporanee
	acc_tot[5] = { 0,0,0,0,0 },///Numero di mosse accettate totali
	idx,//Indice per funzione di distribuzione di coppia
	idxgr = 0,//Altro indice
	Step_print = 1000, 
	stepposidx,
	Step_control = 200,//Numero di iterazioni di controllo dei Delta
	Step_save = 10,//Numero di iterazioni per salvataggio dati
    goodmoves=0,
	test[5][2];
	 double rcm[3];
	double k, T, E, Esum = 0, E2sum = 0, Rmaxsum = 0, Rmax2sum = 0, Rgsum = 0, Rg2sum = 0, R2sum = 0, R22sum = 0, dDelta[5], Delta[5], fraz_acc, Rmax, Rg, R2, r0;
	double acc_max = 0.55, acc_min = 0.45, E_new,angle_max=0.2;
	double* gr, * gr2, * norm,* EM, * RGM, * RMAXM, * R2M;
	int** neighbours,** neighboursprec;
	double shift = LJ(2.5),Ossbackup[4];
	double dR, dR2;
	bool fasestat = false;
	Posizioni *r1,*r2,*temp,*backup;
	  //---INIZIALIZZA LE POSIZIONI INIZIALI DEI POLIMERI
   void init_r()
	{
		int i = 0;
		for (int i = 0; i < N; i++)
		{
			r1->r[i][0] = (double)i - ((double)N) / 2;
			r1->r[i][1] = 0;
			r1->r[i][2] = 0;
			backup->r[i][0]=r1->r[i][0];
            backup->r[i][1]=r1->r[i][1];
            backup->r[i][2]=r1->r[i][2];
		}
	}
   void load_r(string* positions)
   {
	   ifstream file(*positions);
	   string line;
	   int row = 0;
	   int i = 0;
	   while (getline(file, line) && row < N) {
		   istringstream iss(line);
		   string value;
		   int col = 0;
		   if (row < N)
		   {
			   while (getline(iss, value, ';') && col < 3) {
				   r1->r[row][col] = stod(value);
				   backup->r[row][col] = r1->r[row][col];
				   col++;
			   }
		   }
		   row++;
	   }
	   file.close();
   }
	//----------------- TRASLA IL CENTRO DI MASSA NELL'ORIGINE
	void Norm_r()
	{
		rcm[0] = 0;
		rcm[1] = 0;
		rcm[2] = 0;
		for (int i = 0; i < N; i++) {

			rcm[0] += r1->r[i][0];
			rcm[1] += r1->r[i][1];
			rcm[2] += r1->r[i][2];
		}
		rcm[0] /= N;
		rcm[1] /= N;
		rcm[2] /= N;
		for (int i = 0; i < N; i++) {

			r1->r[i][0] -= rcm[0];
			r1->r[i][1] -= rcm[1];
			r1->r[i][2] -= rcm[2];
		}
		Rmax = pow(r1->r[N - 1][0] - r1->r[0][0], 2) + pow(r1->r[N - 1][1] - r1->r[0][1], 2) + pow(r1->r[N - 1][2] - r1->r[0][2], 2);
		Rg = 0;
		R2 = 0;
		for (int i = 0; i < N; i++)
		{
			Rg += pow(r1->r[i][0], 2) + pow(r1->r[i][1], 2) + pow(r1->r[i][2], 2);
			R2 += pow(r1->r[i][0]- r1->r[0][0], 2) + pow(r1->r[i][1]-r1->r[0][1], 2) + pow(r1->r[i][2]-r1->r[0][2], 2);

		}
		Rg /= (N - 1);
		R2 /= (N - 1);
	}

	void norm_gr()
	{
		for (int i = 0; i < N_bin; i++)
		{
			gr[i] /= norm[i] * N_stat;
			gr2[i] /= (double)grcount;
		}
	}
    //----------------- AGGIORNA I DELTA DELLE MOSSE MONTECARLO
	void update_deltas()
	{
		for (int i = 0; i < 5; i++)
		{

			fraz_acc = (double)acc[i] / Step_control;
			
			if (fraz_acc > acc_max) {
				Delta[i] += dDelta[i];
				test[i][1] = 1;
			}
			else if (fraz_acc < acc_min) {
				Delta[i] -= dDelta[i];
				test[i][1] = 0;
			}

			if (Delta[i] <= 0) {
				Delta[i] += dDelta[i];
				dDelta[i] *= 0.5;

			}
			else if (test[i][0] + test[i][1] == 1)
			{
				dDelta[i] *= 0.5;
			}
			test[i][0] = test[i][1];
			acc[i] = 0;
		}
	}
	//-------AGGIORNA LA MATRICE DEI PRIMI VICINI
	void find_neighbours(Posizioni* p)
	{
		double dr2, dr_mod,Eprec=E;
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				neighbours[i][j] = -1;
			}
		}
		for (int i = 1; i < N; i++)
		{
			for (int j = 0; j < i; j++)
			{
				dr2 = 0;
				for (int k = 0; k < 3; k++)
				{
					dr2 += pow((p->r[i][k] - p->r[j][k]), 2);
				}
				dr_mod = sqrt(dr2);
				if (dr_mod < 2.5) {
					neighbours[i][j] = j;
					neighbours[j][i] = -1;
				}
			}
           
		}
		for (int i = 0; i < N; i++)
		{
			// Sort the line to move all -1 values to the end
			std::sort(neighbours[i], neighbours[i] + N, std::greater<int>());
		}
		 calcE(r1);
           if(E_new>Eprec)
			{
              double prob = exp(-(E_new - Eprec) / T);
			  double si=random_half_range()+0.5;
			if(prob<si)
			{
              if(Step_pos>20)
			{
				Step_pos--;
			}
			else if(Step_pos==20&&angle_max>0.1)
			{
              angle_max-=0.005;
			}
			 for (int i = 0; i < N; i++)
			{
                p->r[i][0] = backup->r[i][0];
				p->r[i][1] = backup->r[i][1];
				p->r[i][2] = backup->r[i][2];
				for (int j = 0; j < N; j++)
				{
					neighbours[i][j] = neighboursprec[i][j];
				}
				
			}
	            E=Ossbackup[0];
				Rmax=Ossbackup[1];
				Rg=Ossbackup[2];
				R2=Ossbackup[3];
                goodmoves=0;

			}
			}
			else
			{
                goodmoves+=1;
			    for (int i = 0; i < N; i++)
			    {
				backup->r[i][0] = p->r[i][0];
				backup->r[i][1] = p->r[i][1];
				backup->r[i][2] = p->r[i][2];
				for (int j = 0; j < N; j++)
				{
					neighboursprec[i][j] = neighbours[i][j];
				}
			    }
			  Norm_r();
			  Ossbackup[0]=E_new;
			  Ossbackup[1]=Rmax;
			  Ossbackup[2]=Rg;
			  Ossbackup[3]=R2;
			  E=E_new;
                if(goodmoves==10&&Step_pos<500)
                {
                    Step_pos++;
                    goodmoves=0;
                }
			}
		 stepposidx=0;
			
	}
	void perform_moves()
	{
		if (passo_metropolis())
		{
			acc[0]++;
			acc_tot[0]++;
		}
		if (passo_metropolis2())
		{
			acc[1]++;
			acc_tot[1]++;
		}
		if (gaussian_move())
		{
			acc[2]++;
			acc_tot[2]++;
		}
		if (angle_move())
		{
			acc[3]++;
			acc_tot[3]++;
		}
		if (angle_move2())
		{
			acc[4]++;
			acc_tot[4]++;
		}
		Norm_r();
	}
	//-----INIZIALIZZA LE DISTRIBUZIONI DI COPPIA
	void init_Gr()
	{
		for (int i = 0; i < N_bin; i++)
		{
			norm[i] = 4. / 3 * pi * (pow((i + 1) * dR, 3) - pow(i * dR, 3)) * N;
			gr[i] = 0;
			gr2[i] = 0;
		}
	}
	//---CALCOLA ENERGIA INTERNA DEL POLIMERO
	void calcE(Posizioni* p) {
    double dr2, dr_mod,diff;
    int j, idxgr;
    E_new = 0;
    for (int i = 1; i < N; i++) {
        j = 0;
        while (neighbours[i][j] != -1) {
            dr2 = 0;
            for (int k = 0; k < 3; k++) {
                diff = p->r[i][k] - p->r[neighbours[i][j]][k];
                dr2 += diff * diff;  
            }
            dr_mod = sqrt(dr2);

            if (abs(i - neighbours[i][j]) == 1) {
               
                E_new += Fene(dr_mod, k, r0);  
                if (fasestat && dr_mod < 1.5) {
                    idxgr = (int)(dr_mod / dR2);
                    gr2[idxgr] += 1;
                    grcount++;
                }
            } else {
                E_new += LJ(dr2); 
            }

            if (fasestat && dr_mod <=5) {
                idxgr = (int)(dr_mod / dR);
                gr[idxgr] += 2;
            }
            j++;
        }
    }
}
	void copy_r()
	{
		for (int i = 0; i < N; i++)
		{
			for (int k = 0; k < 3; k++)
			{
				r2->r[i][k] = r1->r[i][k];
			}
		}
	}
	bool testmove()
	{
		double prob;
		calcE(r2);
		if ((random_half_range() + 0.5) < exp(-(E_new - E) / T))
		{
			temp = r1;
			r1 = r2;
			r2 = temp;
			E = E_new;
			return true;
		}
		return false;
	}
    bool passo_metropolis()
	{
		copy_r();
		idx=rand()%N;
        for (int k = 0; k < 3; k++)
		{
			
			r2->r[idx][k] = r1->r[idx][k] + Delta[0] * random_half_range();
		}
		return testmove();
	}
	bool passo_metropolis2()
	{
		for (int i = 0; i < N; i++)
		{
			for (int k = 0; k < 3; k++)
			{
				r2->r[i][k] = r1->r[i][k] + Delta[1] * random_half_range();
			}
		}
		return testmove();
	}
	bool gaussian_move()
	{

		copy_r();
		int particle_idx = rand() % N;
		double move1 = Delta[2] * random_half_range();
		double move2 = Delta[2] * random_half_range();
		double move3 = Delta[2] * random_half_range();
		double distance, dr2 = 0;
		for (int i = 0; i < N; i++)
		{
			for (int k = 0; k < 3; k++)
			{
				dr2 += pow((r2->r[i][k] - r2->r[particle_idx][k]), 2);
			}
			distance = sqrt(dr2);
			r2->r[i][0] = r2->r[i][0] + gaussian(distance, 10) * move1;
			r2->r[i][1] = r2->r[i][1] + gaussian(distance, 10) * move2;
			r2->r[i][2] = r2->r[i][2] + gaussian(distance, 10) * move3;
			dr2 = 0;
		}
		return testmove();
	}
	void fold(double angle)
	{
		int particle_idx = rand() % N;
		bool direction = (rand() % 2);
		double axis_x = (random_half_range() + 0.5) * 2.0 - 1.0;
		double axis_y = (random_half_range() + 0.5) * 2.0 - 1.0;
		double axis_z = (random_half_range() + 0.5) * 2.0 - 1.0;
		double axis_norm = sqrt(axis_x * axis_x + axis_y * axis_y + axis_z * axis_z);
		axis_x /= axis_norm;
		axis_y /= axis_norm;
		axis_z /= axis_norm;
		
		int start, end;
		if (direction)
		{
			start = particle_idx;
			end = N;
		}
		else
		{
			start = 0;
			end = particle_idx;
		}

		for (int i = start; i < end; i++)
		{
			double x = r2->r[i][0] - r2->r[particle_idx][0];
			double y = r2->r[i][1] - r2->r[particle_idx][1];
			double z = r2->r[i][2] - r2->r[particle_idx][2];
			double dot_product = x * axis_x + y * axis_y + z * axis_z;

			double cross_x = y * axis_z - z * axis_y;
			double cross_y = z * axis_x - x * axis_z;
			double cross_z = x * axis_y - y * axis_x;

			x = x * cos(angle) + cross_x * sin(angle) + axis_x * dot_product * (1 - cos(angle));
			y = y * cos(angle) + cross_y * sin(angle) + axis_y * dot_product * (1 - cos(angle));
			z = z * cos(angle) + cross_z * sin(angle) + axis_z * dot_product * (1 - cos(angle));
			r2->r[i][0] = x + r2->r[particle_idx][0];
			r2->r[i][1] = y + r2->r[particle_idx][1];
			r2->r[i][2] = z + r2->r[particle_idx][2];
		}
	}
   bool angle_move()
	{
		if (Delta[3] > angle_max)
		{
			Delta[3] = angle_max;
		}
		copy_r();
		double angle = Delta[3] * random_half_range();
		fold(angle);
		return testmove();
	}
   bool angle_move2()
   {
	   if (Delta[4] > (angle_max/10))
	   {
		   Delta[4] = (angle_max/10);
	   }
	   copy_r();
	   double angle = Delta[4] * random_half_range();
	   for (int j = 0; j < 10; j++)
	   {
		   fold(angle);
	   }
	   return testmove();
   }
  double calculateTau(double* Data, int N, int len,double mean) {
    double *C, *M,tau=0,inv_e = exp(-1);
    int idx,idx2;
    while(tau==0&&len<=N)
    {
    C = new double[len];
    M = new double[len];
    for (int i = 0; i < N; i++) {
        idx = i % len;
        M[idx] = Data[i];
        if (i >= len) {
            for (int j = 0; j < len; j++) {
                idx2 = (idx - j + len) % len;
                C[j] += M[idx] * M[idx2];
            }
        }
    }

    C[0] = C[0] / (double)(N - len + 1) - mean * mean;
    for (int i = 1; i < len; i++) {
        C[i] = ((C[i] / (double)(N - len + 1)) - mean * mean) / C[0];
    }
    C[0] = 1;
    for (int i = 1; i < len; i++) {
				if (C[i] <= inv_e && tau == 0) {

					if (i != 1)
					{
						tau = (inv_e - C[i]) / (C[i] - C[i - 1]) + (double)i;
					}
					else
					{
						tau = (inv_e - C[i]) / (C[i] - 1) + (double)i;
					}
				}
			}
    if(tau==0)
    {
        len*=2;
    }
    delete[] C;
    delete[] M;
    }
   if(tau==0)
   {
    return -1;
   }
   return tau;
}
};