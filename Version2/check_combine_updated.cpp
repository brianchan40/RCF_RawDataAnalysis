#include <list>
#include <string>
#include <cstring>
#include "TFile.h"
#include <iostream>
#include <fstream>
using namespace std;

void check_combine_updated(char *name_num)
{
    int max = 0;

    char command[200];
    char filename[200];
    char outputfilename[200];
    char resubcommand[200];
    
    int count = 0;
    max = 10398;
    
    TFile *h_file;
    TFile *t_file;

        for(int cent = 0; cent <= 0; cent++)
        {
            sprintf(outputfilename, "./missing_%d.sh", cent);
            ofstream outputfile;
            outputfile.open (outputfilename);
            
        	for(int i = 0; i <= max; i++){
                
                if((count%10) == 1) cout << " " << cent << " " << i << endl;

                sprintf(filename, "./output/lam_low/Data_14_4/sched%s_%d_plam.root", name_num, i);
                h_file = new TFile(filename);
                if(!h_file->Get("V0Mass")){
                    cout << filename << " didn't work" << endl;
                    sprintf(resubcommand, "rm ./output/lam_low/Data_14_4/sched%s_%d_* \n", name_num, i);
                    outputfile << resubcommand;
                    h_file->Close();
	            continue;
                }
                
                sprintf(filename, "./output/lam_low/Data_14_4/sched%s_%d_tree.root", name_num, i);
                t_file = new TFile(filename);
                if(!t_file->IsOpen()){
                    cout << filename << " didn't work" << endl;
                    sprintf(resubcommand, "rm ./output/lam_low/Data_14_4/sched%s_%d_* \n", name_num, i);
                    outputfile << resubcommand;
                    h_file->Close();
		    t_file->Close();
		    continue;
                }
                
                t_file->Close();
                h_file->Close();
                
                count++;

            }
        }

    //}
    outputfile.close();
}
