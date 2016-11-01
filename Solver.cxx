#include<iostream>
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
        //cout << "Initializing" << "\n";
        for(int i=0; i<nx; i++ ){
            Phi[i] = i;
            //cout << Phi[i];
        }
      }

      void showList(double* Phi)
      {
        cout << "\n";
        for(int i=0; i<nx; i++ ){
            cout << Phi[i] << " ";
        }
      }

      void solveForPhi(double* Phi, double* Rho){
        double Phi2[nx-2];
        for(int i =0; i<nx-2; i++)
        {
         Phi2[i] = Phi[i]+Phi[i+2];
        }

        for(int i =0 ; i<nx-2; i++){
            Phi[i+1] = (Phi2[i] - Rho[i])/2;
        }
      }

      void linspace(double* List, double min, double max, int n){
          //double Result[n];
          for(int i = 0; i< n; i++){
            List[i] = min+(max-min)*i/(n-1);
          }
      }

      writeFile(String Name, double* List1, double* List2, int n)
      {
        ofstream myFile;
        for(int i =0; i <n; i ++){
            myFile.open()
        }
      }

int main () {
  ofstream myfile;
  myfile.open ("example.txt");
  myfile << "Writing this to a file.\n";
  myfile.close();
  return 0;
}

};

main()
{
   programming a;

   a.input_value();
   a.output_value();
   double Phi[a.nx];
   a.initializePhi(Phi);
   double Rho[a.nx];
   a.linspace(Rho, -1.0, 1.0, a.nx);
   //a.showList(Rho);
   a.solveForPhi(Phi, Rho);
   a.showList(Phi);

   //cout << a.nx ;
   //a.Poisson();

   //object.variable;  Will produce an error because variable is private

   return 0;
}
