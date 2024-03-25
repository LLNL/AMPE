#include "TemperatureHistory.h"

#include "SAMRAI/tbox/PIO.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <regex>
#include <cassert>
#include <algorithm>

using namespace SAMRAI;

void TemperatureHistory::readCSV(const std::string& filename)
{
   d_os << "Read temperatures from CSV file..." << std::endl;
   std::ifstream input;
   input.open(filename);
   const std::regex delimiter(",");

   // following rows store position followed by temperatures at all times
   std::vector<std::vector<std::string>> csvData;
   for (std::string line{}; std::getline(input, line);) {
      // std::cout<<line<<std::endl;
      // skip lines starting with '#' (comments)
      if (line[0] == '#') continue;
      csvData.emplace_back(std::vector<std::string>(
          std::sregex_token_iterator(line.begin(), line.end(), delimiter, -1),
          {}));
   }

   // first row stores positions
   std::vector<std::string> v(csvData[0]);
   for (auto e : v) {
      d_positions.push_back(std::stod(e));
   }
   std::reverse(d_positions.begin(), d_positions.end());
   d_os << "Number of locations: " << d_positions.size() << std::endl;

   bool first_row = true;
   for (auto v : csvData) {
      if (!first_row) {  // skip first row already read
         std::vector<double> row;
         bool first_col = true;
         for (auto s : v) {
            if (first_col) {
               double time = std::stod(s);
               d_os << "Time: " << time << std::endl;
               d_times.push_back(time);
            } else
               row.push_back(std::stod(s));
            first_col = false;
         }
         std::reverse(row.begin(), row.end());
         assert(row.size() == d_positions.size());
         for (auto t : row)
            d_os << t << std::endl;
         d_temperatures.push_back(row);
      }
      first_row = false;
   }
   d_os << "Number of times: " << d_times.size() << std::endl;

   assert(d_temperatures.size() == d_times.size());

   input.close();
}

int TemperatureHistory::getTandGradT(const double time, const double position,
                                     double& temperature, double& gradT)
{
   assert(d_times.size() > 1);
   assert(d_positions.size() > 1);
   // d_os << "TemperatureHistory::getTandGradT()..."<<std::endl;
   // find times to interpolate in between
   if (time < d_times[0]) {
      d_os << "Time " << time
           << " is beyond range defined in TemperatureHistory" << std::endl;
      d_os << "Minimum time is " << d_times[0] << std::endl;
      return 1;
   }
   unsigned itime = 0;
   for (auto t : d_times) {
      if (time >= t)
         itime++;
      else
         break;
   }
   itime--;
   assert(itime >= 0);
   if ((itime + 1) >= d_times.size()) {
      d_os << "Time " << time
           << " is beyond range defined in TemperatureHistory" << std::endl;
      d_os << "Maximum time is " << d_times[d_times.size() - 1] << std::endl;
      return 1;
   }
   const double tm = d_times[itime];
   const double tp = d_times[itime + 1];
   const double tfrac = (time - tm) / (tp - tm);

   // find positions to interpolate in between
   if (position < d_positions[0]) {
      d_os << "Position " << position
           << " is beyond range defined in TemperatureHistory" << std::endl;
      d_os << "Minimum position is " << d_positions[0] << std::endl;
      return 1;
   }
   unsigned ipos = 0;
   for (auto x : d_positions) {
      if (position > x)
         ipos++;
      else
         break;
   }
   ipos--;
   if ((ipos + 1) >= d_positions.size()) {
      d_os << "Position " << position
           << " is beyond range defined in TemperatureHistory" << std::endl;
      d_os << "Maximum position is " << d_positions[d_positions.size() - 1]
           << std::endl;
      return 1;
   }
   const double xm = d_positions[ipos];
   const double xp = d_positions[ipos + 1];
   const double xfrac = (position - xm) / (xp - xm);

   // get temperatures at four surrounding corners
   const double temp_tm_xm = d_temperatures[itime][ipos];
   const double temp_tm_xp = d_temperatures[itime][ipos + 1];
   const double temp_tp_xm = d_temperatures[itime + 1][ipos];
   const double temp_tp_xp = d_temperatures[itime + 1][ipos + 1];

   // evaluate temperature at "position" by interpolation
   temperature =
       (1. - tfrac) * ((1. - xfrac) * temp_tm_xm + xfrac * temp_tm_xp) +
       tfrac * ((1. - xfrac) * temp_tp_xm + xfrac * temp_tp_xp);

   // now evaluate gradient of T at "position"
   double grad_tm = (temp_tm_xp - temp_tm_xm) / (xp - xm);
   double grad_tp = (temp_tp_xp - temp_tp_xm) / (xp - xm);

   // add weighted contributions from interval on left or right
   if (ipos > 0 && xfrac < 0.5) {
      const double xmm = d_positions[ipos - 1];

      const double temp_tm_xmm = d_temperatures[itime][ipos - 1];
      double grad_tm_xm = (temp_tm_xm - temp_tm_xmm) / (xm - xmm);

      const double temp_tp_xmm = d_temperatures[itime + 1][ipos - 1];
      double grad_tp_xm = (temp_tp_xm - temp_tp_xmm) / (xm - xmm);

      grad_tm = (0.5 + xfrac) * grad_tm + (0.5 - xfrac) * grad_tm_xm;
      grad_tp = (0.5 + xfrac) * grad_tp + (0.5 - xfrac) * grad_tp_xm;
   } else if ((ipos + 2 < d_positions.size()) && xfrac > 0.5) {
      const double xpp = d_positions[ipos + 2];

      const double temp_tm_xpp = d_temperatures[itime][ipos + 2];
      double grad_tm_xp = (temp_tm_xpp - temp_tm_xp) / (xpp - xp);

      const double temp_tp_xpp = d_temperatures[itime + 1][ipos + 2];
      double grad_tp_xp = (temp_tp_xpp - temp_tp_xp) / (xpp - xp);

      grad_tm = (1.5 - xfrac) * grad_tm + (xfrac - 0.5) * grad_tm_xp;
      grad_tp = (1.5 - xfrac) * grad_tp + (xfrac - 0.5) * grad_tp_xp;
   }

   // finally interpolate in time
   gradT = (1. - tfrac) * grad_tm + tfrac * grad_tp;

   return 0;
}
