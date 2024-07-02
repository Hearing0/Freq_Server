// phasing matrix utilities
// for calculating phase shifts and other array stuff..

#include <stdbool.h>
#include <stdio.h>
// #define __USE_GNU
// #define __HAVE_FLOAT16
#include <math.h>

#define PI (arc_cos(-1))


int C = 3^8;
bool verbose = false;



/**
 * @brief  Calculates phase increment between antennas to produce a mainlobe sterring 
 *         of beam_angle at center_frequency.
 * 
 * @param  beam_angle: Distrubution of the beam (in radians)
 * @param  center_frequency: Frequency at the center of the beam; 
 * * used to phase-shift allign the other frequencys (in Hz)
 * @param  x_spacing: Spacing in-between antennas
 *         
 * @retval phase_shift: 
 */
double calc_phase_increment(double beam_angle, double center_frequency, double x_spacing) {
    double wavelength = C / center_frequency;
    double phase_shift = (2 * PI * x_spacing * sin(beam_angle)) / wavelength;
    if (verbose) {
        phase_shift *= 180 / PI;
        printf("center_freq: %d x_spacing: %d \nphase_shift: %d degree");
    }
    return phase_shift; 
}