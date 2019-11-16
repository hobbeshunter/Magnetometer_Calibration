# Magnetometer_Calibration
This repository contains a platformio library and an explanatory jupyter notebook for calculating the hard and soft iron offsets of an magnetometer as descriped in [Calibrating an eCompass in the Presence of Hard- and Soft-Iron Interference](https://www.nxp.com/docs/en/application-note/AN4246.pdf) by NXP.

## Jupyter notebook

Simply run `jupyter notebook` int his folder to start jupyter. Then open the mag_calibration notebook. There you all informations needed.

## C++ implementation

This implementation depends on Eigen. For platform.io projects you can use [Eigen_Platformio_Header](https://github.com/hobbeshunter/Eigen_Platformio_Header).

Usage:
```C++
#define NO_CALIBRATION_MEASUREMENTS 5000
Eigen::Array<float, NO_CALIBRATION_MEASUREMENTS, 1> mxs;
Eigen::Array<float, NO_CALIBRATION_MEASUREMENTS, 1> mys;
Eigen::Array<float, NO_CALIBRATION_MEASUREMENTS, 1> mzs;
Eigen::Array<float, NO_CALIBRATION_MEASUREMENTS, 1> rolls;
Eigen::Array<float, NO_CALIBRATION_MEASUREMENTS, 1> pitchs;

// Fill measurement arrays

Magnetometer_Calibration mag_cal;
mag_cal.calibrate(mxs, mys, mzs, rolls, pitchs);
Winv = mag_cal.getWinv();
W = mag_cal.getW();
V = mag_cal.getV();
B = mag_cal.getB();
incl = mag_cal.getInclination();
```