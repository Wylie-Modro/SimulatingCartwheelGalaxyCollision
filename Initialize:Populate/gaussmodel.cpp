// Midterm
// Pappas

// initialization code

#include <fstream>
#include <iostream>
#include <cmath>
#include <ctime>
#include <cstdlib>

using namespace std;

struct initializer
{
    double rstar[4];
};

struct radiushold
{
    double rstartheo, theopopdens;
};

int main()
{
    //Define Consts. (for)
    const int N = 500;
    const int N_STEPS = 50000;
    const int N_REGIONS = 5000;
    double random;
    
    //Define Structures
    initializer rgen[N];
    radiushold rhold[N_STEPS]; 
    
    // Define Constants and Changed Variables
    double Rmax = 0.999999;
    double Rmin = 0.000001;
    double amillion = pow(10,6); 
    double ahundredthousand = pow(10,5);
    double tenthousand = pow(10,4);
    double thousand = pow(10,3);
    double onehundred = 100;
    double ten = 10;
    int qi = 0;
    int qf = 5;
    int ri = 0;
    int qc = 0;
    double r0 = Rmax/2;
    double X1, X2; 
    bool swap = true;
    double temp;
    int treat;
    double count = 0;
    double riti, ritf;
    
    //Define Files
    fstream theoFile;
    fstream stardistFile;
    fstream randFile;
    fstream stardistorganizedFile;
    fstream stardistcheckFile;
    
    //Open Files
    theoFile.open("theopopdens.txt", ios::out);
    randFile.open("rands.txt", ios::out);
    
    // Populate the Theoretical Star Density 
    for(int i = 0; i < N_STEPS; i++)
    {
        rhold[i].rstartheo = ((double)(i) * 1/N_STEPS);
        rhold[i].theopopdens = (rhold[i].rstartheo/r0) * (exp(( - (pow(rhold[i].rstartheo,2)))/(2 * pow(r0,2))));
        theoFile << rhold[i].rstartheo << " " << rhold[i].theopopdens << " " << endl;
    }      
    
    // Create a seed
    unsigned int seed = time(NULL);
    srand(seed); 
    
    // initialize stars by setting 0 < X1 < g(Rmax) & 0 < g(X2) < X1 and if so then a point of the dist. is X2
    do 
    {
        X1 = (((double)rand())/(RAND_MAX)); // generate a random value between 0 and 1
        X2 = ((((double)rand())/(RAND_MAX))*100); // ditto 
        cout << X1 << " " << X2 << " " << endl;
        if(X1 < ((1/r0) * (exp((-1 / (2*pow(r0,2)))))))
        {
            if((X2 * exp((- (pow(X2,2)/(2))))) < X1)
            {
                rgen[ri].rstar[0] = X2 * r0; // set the ri-th star equal to the randomly selected radius
                ri++; // increase the count of the radius for the next star
            }
        }
    }while(ri < N);
        
    theoFile.close();
    randFile.close();
    
    stardistFile.open("stardist.txt",ios::out);

    //output dist.
    for(int i = 0; i < N; i++)
    {
        
        stardistFile << rgen[i].rstar[0] << " " << endl;
    
    }
    
    stardistorganizedFile.open("stardistorg.txt", ios::out);
    
    //sort distribution
    do
    {
        treat = 0;
        for(int i = 0; i < N; i++)
        {
            if(rgen[i].rstar[0] > rgen[i+1].rstar[0])
            {
                temp = rgen[i+1].rstar[0];
                rgen[i+1].rstar[0] = rgen[i].rstar[0];
                rgen[i].rstar[0] = temp;
                swap = true;
                treat = 1;
            }
        }
        if(treat == 0)
        {
            swap = false;
        }
        
    }while(swap == true);
    
    //output sorted dist. 
    for(int i = 0; i < N; i++)
    {
        
        stardistorganizedFile << rgen[i].rstar[0] << " " << endl;
    
    }
    
    //close stardist files
    stardistorganizedFile.close();
    stardistFile.close();
    
    //open stardistcheckfile
    stardistcheckFile.open("stardistcheck.txt", ios::out);
    
    // Check distribution
    for(int j = 0; j + 1 < N_STEPS; j = j + 10)
    {
        count = 0;
        for(int i = 0; i < N; i++)
        {
            riti = double(j)/N_STEPS;
            ritf = (double(j)+10)/N_STEPS;
            if(rgen[i].rstar[0] >= riti && rgen[i].rstar[0] <= ritf)
            { 
               count = count + 1/((double)N);
            }
        }
        stardistcheckFile << riti << " " << ritf << " " << count << " " << endl;
    }
        
    stardistcheckFile.close();
    
    return 0;
}
