#include <list>
#include <string>
#include <cstring>
#include "TFile.h"
#include <iostream>
#include <fstream>
using namespace std;

void check_combine()
{
    list<char *> daynum;
    list<char *>::iterator day;
    int max = 0;

    char command[200];
    char filename[200];
    char outputfilename[200];
    char resubcommand[200];

    daynum.push_back("09");
    daynum.push_back("10");
    daynum.push_back("11");
    //daynum.push_back("12");
    //daynum.push_back("13");
    //daynum.push_back("14");
    //daynum.push_back("15");
    //daynum.push_back("16");
    
    int count = 0;

    for(day = daynum.begin(); day != daynum.end(); ++day)
    {
        if(strcmp(*day, "09") == 0)
        {
            max = 247;
        }
        else if(strcmp(*day, "10") == 0)
        {
            max = 525;
        }
        else if(strcmp(*day, "11") == 0)
        {
            max = 250;
        }
	else if(strcmp(*day, "12") == 0)
        {
            max = 130;
        }
	else if(strcmp(*day, "13") == 0)
        {
            max = 75;
        }
	else if(strcmp(*day, "14") == 0)
        {
            max = 93;
        }
	else if(strcmp(*day, "15") == 0)
        {
            max = 155;
        }
	else if(strcmp(*day, "16") == 0)
        {
            max = 103;
        }

        for(int cent = 5; cent <= 8; cent++)
        {
            sprintf(outputfilename, "./missing_%s_%d.sh", *day, cent);
            ofstream outputfile;
            outputfile.open (outputfilename);
            
        	for(int i = 0; i <= max; i++){
                
                if((count%10) == 1) cout << *day << " " << cent << " " << i << endl;

                //sprintf(outputfilename, "./temp/Data_%s_%d/Data_%d/missing.sh", *day, cent, i);
                //TFile *outputfile = new TFile(outputfilename);

                sprintf(filename, "./temp/Data_%s_%d/Data_%d/t_%dcen%d.weight_112_module_new.root", *day, cent, i, i, cent);
                TFile *w_file = new TFile(filename); 
                if(!w_file->Get("TPCmeanPhi_FF_2")){
                    cout << filename << " didn't work" << endl;
                    sprintf(resubcommand, "rm ./temp/Data_%s_%d/Data_%d/t_* \n", *day, cent, i);
                    outputfile << resubcommand;
                    sprintf(resubcommand, "./while_loop.sh %s %d %d %d %d \n", *day, max, (max - i), (max - i), cent);
                    outputfile << resubcommand;
                    continue;
                }

                sprintf(filename, "./temp/Data_%s_%d/Data_%d/t_%dcen%d.v2_fullEP_eff_pT02_module.root", *day, cent, i, i, cent);
                TFile *v_file = new TFile(filename);
                if(!v_file->Get("Hist_cos")){
                    cout << filename << " didn't work" << endl;
                    sprintf(resubcommand, "rm ./temp/Data_%s_%d/Data_%d/t_* \n", *day, cent, i);
                    outputfile << resubcommand;
                    sprintf(resubcommand, "./while_loop.sh %s %d %d %d %d \n", *day, max, (max - i), (max - i), cent);
                    outputfile << resubcommand;
                    continue;
                }

		/*sprintf(filename, "./temp/Data_%s_%d/Data_%d/t_%dcen%d.gamma112_fullEP_eff_pT02_module.root", *day, cent, i, i, cent);
                TFile *v_file = new TFile(filename);
                if(!v_file->Get("Hist_cos")){
                    cout << filename << " didn't work" << endl;
                    sprintf(resubcommand, "rm ./temp/Data_%s_%d/Data_%d/t_* \n", *day, cent, i);
                    outputfile << resubcommand;
                    sprintf(resubcommand, "./while_gamma.sh %s %d %d %d %d \n", *day, max, (max - i), (max - i), cent);
                    outputfile << resubcommand;
                    continue;
                }

                sprintf(filename, "./temp/Data_%s_%d/Data_%d/t_%dplam.root", *day, cent, i, i, cent);
                TFile *h_file = new TFile(filename);
                if(!h_file->Get("V0Mass")){
                    cout << filename << " didn't work" << endl;
                    sprintf(resubcommand, "rm ./temp/Data_%s_%d/Data_%d/t_* \n", *day, cent, i);
                    outputfile << resubcommand;
                    sprintf(resubcommand, "./while_gamma.sh %s %d %d %d %d \n", *day, max, (max - i), (max - i), cent);
                    outputfile << resubcommand;
                    continue;
                }*/

		sprintf(filename, "./temp/Data_%s_%d/Data_%d/t_%dplam.root", *day, cent, i, i, cent);
                TFile *h_file = new TFile(filename);
                if(!h_file->Get("V0Mass")){
			cout << filename << " didn't work" << endl;
                    sprintf(resubcommand, "rm ./temp/Data_%s_%d/Data_%d/t_* \n", *day, cent, i);
                    outputfile << resubcommand;
                    sprintf(resubcommand, "./while_loop.sh %s %d %d %d %d \n", *day, max, (max - i), (max - i), cent);
                    outputfile << resubcommand;
                    continue;
                }

                
                w_file->Close();
                v_file->Close();
                h_file->Close();
                
                count++;

            }
        }

    }
    outputfile.close();
}
