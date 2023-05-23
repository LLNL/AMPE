#include <vector>
#include <string>
#include <iostream>

class TemperatureHistory
{
 public:
   TemperatureHistory(std::ostream& os) : d_os(os){};

   void readCSV(const std::string& filename);

   // get T and grad(T) at a specific time and position by interpolating
   // in between points in 2D array d_temperatures
   int getTandGradT(const double time, const double position,
                    double& temperature, double& gradT);

 private:
   std::ostream& d_os;

   std::vector<double> d_times;
   std::vector<double> d_positions;

   /// temperatures at times and positions
   /// Each row corresponds to spatial data at a given time
   std::vector<std::vector<double>> d_temperatures;
};
