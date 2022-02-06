//
// Created by F.Moitzi on 20.01.2022.
//

#ifndef SRC_LIB_MULTIPOLEMADELUNG_MULTIPOLEMOMENTS_HPP
#define SRC_LIB_MULTIPOLEMADELUNG_MULTIPOLEMOMENTS_HPP

#include <vector>
#include <complex>

#include "NDArray.hpp"

namespace lsms {

  namespace multipole_moms {


    /**
     * Construct the multipole moments similar to MuST to have better comparison
     */
    std::vector<std::complex<double>> construct_multi_mom(
        const std::vector<double> &rmesh,
        const lsms::NDArray<std::complex<double>, 2> &density_tilda_l,
        std::size_t nr_ps,
        std::size_t lmax,
        double tolerance = 1e-8
    );


  }

}


class MultipoleMoments {

};


#endif //SRC_LIB_MULTIPOLEMADELUNG_MULTIPOLEMOMENTS_HPP
