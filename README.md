#Speaker Driver Parameters Calculator

This is a Python-based tool for calculating the Thiele-Small loudspeaker parameters from an input wave file.

We feed a white noise signal from a low-impedance source into a voltage divider circuit that includes a test resistor and the target speaker driver, as follows:

```
AMP OUT -----+--------> LINE In Right
             |
           Ztest
             |
             +--------> LINE In Left
             |   
          Rseries
             |
GND----------+--------> LINE In GND
```
The wave file used as an input to the program should be the left and right channels as measured above.

The paramaters passed to the program should be the filename and the resistance of the series resistor (in Ohms):

```
python tsparam.py input.wav 10
```

The software should return the following parameters:

Fs (resonant frequency)
Qms (mechanical damping factor)
Qes (electrical damping factor)
Qts (total damping factor)
Zmax (peak impedance)
Re (voice coil impedance)

Under development (as of August 2018) is the ability to calculate Vas (equivalent volume). Enter the changed resonant frequency Fs as mass is added to the cone, or as the driver is confined to a sealed box.
