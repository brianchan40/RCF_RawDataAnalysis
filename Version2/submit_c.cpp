#include <list>
#include <string>
#include <cstring>
using namespace std;

void submit_c()
{
    list<char *> daynum;
    list<char *>::iterator day;
    int max = 0;

    char command[200];

    daynum.push_back("09");
    daynum.push_back("10");
    daynum.push_back("11");
    //daynum.push_back("12");
    //daynum.push_back("13");
    //daynum.push_back("14");
    //daynum.push_back("15");
    //daynum.push_back("16");

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
        	int temp_max = 0;

            for(int first = 0; first < (int)max/20; first++)
            {
                sprintf(command, "./submit_combine.sh %s %d %d %d %d", *day, cent, (first*20), ((first+1)*20-1), max);
                //sprintf(command, "./submit_gamma.sh %s %d %d %d %d", *day, cent, (first*20), ((first+1)*20-1), max);
		//cout << command << endl;
                system(command);
                temp_max = (first*20-1);
            }

            sprintf(command, "./submit_combine.sh %s %d %d %d %d", *day, cent, (temp_max+1), max, max);
            //sprintf(command, "./submit_gamma.sh %s %d %d %d %d", *day, cent, (first*20), ((first+1)*20-1), max);
	    cout << command << endl;
            
            system(command);
        }

    }
}
