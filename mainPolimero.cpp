#include <iostream>
#include <sstream> 
#include <iomanip>
#include "Polimero - update2.0.h"
using namespace std;

int main()
{
	//Nomi dei file su cui salvare i dati
	string path;
	string dati_pos;
	string dati_load_pos;
	string dati_oss;
	string dati_gr;
    
	double T= 4.5; //Temperatura adimensionale
	double k = 10;//Densità adimensionale
	double x0 = 1.5;
	double Delta = 0.01;//Delta per le mosse sulle coordinate
	double costanti_d[4] = {T,k,x0,Delta};//Vettore costanti double
	int N = 1;//Numero di monomeri
	int N_eq = 500000; //Numero di passi per la fase di equilibratura
	int N_stat = 500000; //Numero di passi per la fase statistica
	int N_bin = 1000;

	int costanti_i[4] = {N,N_eq,N_stat,N_bin}; //Vettore costanti int
	double oss[13];
   double Gr[1000];
	double Dist[1000];
	double Gr2[1000];
	double Dist2[1000];

    path="C:/Users/trent/Desktop/Tesi2.0/PoliDef/10/";
    dati_oss = path + string("ossT45(-1).csv");
    dati_gr=path + string("grT45(-1).csv");
    ofstream filesave,filegr;
	filesave.open(dati_oss);
	filegr.open(dati_gr);
	int numb[16]={3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90};
	int i=0;
    Polimero* s[16];
    while (i<16)
	{
		N=numb[i];
		dati_pos = path + to_string(N) + string("posT45.csv");
		//dati_load_pos=path+to_string(N)+string("pos")+"T05.csv";
		costanti_i[0] = N;
		//Polimero* s = new Polimero(costanti_d, costanti_i);
        s[i] = new Polimero(costanti_d, costanti_i);
		s[i]->CreaSistema();
		s[i]->print_condizioni_iniziali();
		s[i]->fase_equilibratura();
		s[i]->fase_statistica(Gr, Dist, Gr2, Dist2);
		s[i]->calcola_osservabili(oss,&dati_pos);
		filesave << oss[0] << ';';
		filesave << oss[1] << ';';
		filesave << oss[2] << ';';
		filesave << oss[3] << ';';
		filesave << oss[4] << ';';
		filesave << oss[5] << ';';
		filesave << oss[6] << ';';
		filesave << oss[7] << ';';
		filesave << oss[8] << ';';
		filesave << oss[9] << ';';
		filesave << oss[10] << ';';
		filesave << oss[11] << ';';
		filesave << oss[12] << ';'<<endl; 
		
		for (int i = 0; i < N_bin; i++)
		{
			filegr << Dist[i] << ";" << Gr[i] << ";" << Dist2[i] << ";" << Gr2[i] << endl;
		}
		//delete s;
		i++;
	}
   	filegr.close();
	filesave.close();
   

	return 0;
}

