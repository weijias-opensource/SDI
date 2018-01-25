SEISMIC DAYLIGHT IMAGING
========================

This novel receiver-based approach is based on the priciple by Claerbout (1968), 
which constructes the P reflectivity via autocorrelation. 
For detailed descriptions, please refer to Sun and Kennett (2016).

Citations
=========
If you use the package, the author would appreciate it if you could cite the listed references below.
1. Weijia Sun and B. L. N. Kennett, 2016, Receiver structure from teleseisms: Autocorrelation and cross correlation, Geophys Res Lett, 43, 6234–6242.
2. Weijia Sun and B. L. N. Kennett, 2017, Mid–lithosphere discontinuities beneath the western and central North China Craton, Geophys Res Lett, 44, 1302–1310.

Installation
============
Put the python source files in any folder you select and add the path to the PYTHONPATH enviroment.

For bashrc,

	export PYTHONPATH=$PYTHONPATH:/where/you/put/the/python/files
	
The way is very convenient for users without root access.

Before you run the example, make sure you have installed the package of OBSPY https://github.com/obspy/obspy/wiki.

How to use
==========
First you should import SDI::
	
	import sdi
	
Then apply these procedures.

	sdi.pre-processing()
	sdi.autocorrelation()
	sdi.moveout()
	sdi.stack()

For detailed, please go to the examples folder and check the script of run_sdi.py and plot.py.	

Applications
============
If you use the package and publish your work, Please send a copy to me via swj [at] mail.iggcas.ac.cn. I would list your work here.

Contacts
========
If you have any problems, please contact with me. I would response your problems as quick.

Email: swj<at>mail.iggcas.ac.cn; weijia_sun<at>163.com