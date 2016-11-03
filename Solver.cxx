#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;


class programming
{
   private:

   public:
      int nx;

      void input_value()
      {
         cout << "In function input_value, Enter an integer\n";
         cin >> nx;
      }

      void output_value()
      {
         cout << "Variable entered is ";
         cout << nx << "\n";
      }

      void initializePhi(double* Phi)
      {
        for(int i=0; i<nx; i++ ){
            Phi[i] = 0;
        }
      }

      void showList(double* Phi)
      {
        cout << "\n";
        for(int i=0; i<nx; i++ ){
            cout << Phi[i] << " ";
        }
      }

      /*void solveForPhi(double* Phi, double* Rho){
        double Phi2[nx-2];
        for(int i =0; i<nx-2; i++)
        {
         Phi2[i] = Phi[i]+Phi[i+2];
        }

        for(int i =0 ; i<nx-2; i++){
            Phi[i+1] = (Phi2[i] - Rho[i])/2;
        }
      }*/

      void solveForPhiLU(double* X, double* Rho, double* Phi, int n){
        const double dx = X[1] - X[0];
        // Saving the last value :
        const double Phi_end = Phi[n-1];
        // Initial Guess:
        Phi[1] = Phi[0]; //+ Rho[1]*dx;
        // Calculation of Potential Profile
        for(int j = 0; j < n-2; j++){
            Phi[j+2] = 2*Phi[j+1] - Phi[j] + Rho[j+1]*dx*dx;
        }
        // Termination of linear degree of freedom :
        const double Phi_1 = Phi_end - Phi[n-1];
        const double X_range = X[n-1] - X[0];
        for(int i=0; i<n; i++){
            Phi[i] += Phi_1*(X[i]-X[0])/X_range;
        }
      }

      void linspace(double* List, double min, double max, int n){
        for(int i = 0; i< n; i++){
            List[i] = min+(max-min)*i/(n-1);
        }
      }

      void writeFile(string Name, double* List1, double* List2, double* List3, int n)
      {
        ofstream myFile;
        myFile.open("Example.txt");
        for(int i =0; i <n; i ++){
            myFile << List1[i] << "\t" << List2[i] << "\t" << List3[i] << "\n";
        }
        myFile.close();
      }

};

main()
{
    programming a;

    // Input of Grid frequency
    a.input_value();
    a.output_value();

    // Cration of X-array
    double X[a.nx];
    a.linspace(X, -1.0, 2.0, a.nx);

    // Creation of Potential-Array
    double Phi[a.nx];
    a.initializePhi(Phi);

    double Rho[a.nx];
    // Definition of inital Conditions
    // Definition of Charge Density
    /*
    Example #1
    double Rho0 = 6.0;
    a.linspace(Rho, -Rho0 , 2.0*Rho0, a.nx);
    // Initial Condition for the Potential
    Phi[0] = -1; Phi[a.nx-1] = 8;

    // Example #2
    for(int i = 0; i<a.nx ; i++){
        Rho[i] = X[i]*X[i]*12.0;
    }
    // Initial Condition for the Potential
    Phi[0] = 1; Phi[a.nx-1] = 16;

    // Example #3
    for(int i = 0; i<a.nx ; i++){
        Rho[i] = sin(X[i]*2.0*M_PI);
    }
    // Initial Condition for the Potential
    Phi[0] = 0; Phi[a.nx-1] = 0;
    */
    // Example #4
    for(int i = 0; i<a.nx ; i++){
        Rho[i] = exp(X[i]);
    }
    // Initial Condition for the Potential
    Phi[0] = exp(X[0]); Phi[a.nx-1] = exp(X[a.nx-1]);



    // Applying the Solver
    a.solveForPhiLU(X, Rho, Phi, a.nx);
    a.writeFile("File1", X, Rho, Phi, a.nx);

    return 0;
}
//object.variable;  Will produce an error because variable is private
