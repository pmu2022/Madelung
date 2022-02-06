//
// Created by F.Moitzi on 20.01.2022.
//

#include "MultipoleMoments.hpp"

#include "integer_factors.hpp"
#include "integrator.hpp"

std::vector<std::complex<double>> lsms::multipole_moms::construct_multi_mom(const std::vector<double> &rmesh,
                                                                            const lsms::NDArray<std::complex<double>, 2> &density_tilda_l,
                                                                            std::size_t nr_ps,
                                                                            std::size_t lmax,
                                                                            double tolerance) {

  std::vector<double> sqrt_r(nr_ps);

  for(auto ir = 0; ir < nr_ps; ir++) {
    sqrt_r[ir] = std::sqrt(rmesh[ir]);
  }

  auto jmax = lsms::get_jmax(lmax);

  std::vector<std::complex<double>> multi_moms(jmax, 1.0);

  std::vector<std::complex<double>> func(nr_ps, 0.0);

  auto jl = 0;

  for (auto l = 0; l <= lmax; l++) {

    for (auto m = 0; m <= l; m++) {

      // Check tolerance
      bool is_small = true;

      for (auto ir = 0; ir < nr_ps; ir++) {
        if (std::abs(density_tilda_l(ir, jl)) >= tolerance) {
          is_small = false;
        }
      }

      if (!is_small) {

        // Create the right weight
        for (auto ir = 0; ir < nr_ps; ir++) {
          func[ir] = std::pow(rmesh[ir], 2 + l) * 4 * M_PI * density_tilda_l(ir, jl);
        }

        auto res = lsms::simpson_nonuniform(rmesh, func, nr_ps);

        multi_moms[jl] = res;
      }


      jl++;

    }

  }

  return multi_moms;
}
