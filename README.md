# cxd5602pwbimu_logger
CXD5602PWBIMU_LOCALIZER : Sample code to localize the IMU raw data

This code is a modified version of the code in the following repository example.
https://github.com/sonydevworld/spresense

In addition to the angular velocity and acceleration data from the sensor, the program outputs quaternions, speed, and position via uart.
After writing the program to the spresense, you can use test.py to write the measurement data to log.csv.
```
python3 test.py
```