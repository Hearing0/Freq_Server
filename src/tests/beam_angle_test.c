#include <stdio.h>
#include <math.h>

// Define PI if not already defined
#ifndef PI
#define PI 3.14159265358979323846
#endif

int verbose = 1; // Set to 1 for verbose output

double calc_beam_angle(int n_beams, int beam_num, double beam_sep) {
    // Calculate Beamforming shift
    double center_beam = (n_beams - 1) / 2;

    // Calculate Beam Azimuth
    double b_azi = ((beam_num - center_beam) * beam_sep) * (PI / 180);
    if (verbose) {
        printf("n_beams: %d, beam_num: %d, beam_sep: %lf\n", n_beams, beam_num, beam_sep);
        printf("    beam = %lf degree\n", (b_azi * 180 / PI));
    }
    return b_azi;
}

int main() {
    // Example usage of calc_beam_angle
    double angle = calc_beam_angle(16, 3, 3.24);
    printf("Calculated beam angle: %lf radians\n", angle);

    printf("\n\nTesting\n");

    int test = -1;
    if (test > -1) printf("test exists!\n");
    else if (test == -1) printf("test doesn't exist.\n");

    return 0;
}
