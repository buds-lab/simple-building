#SimpleBuilding

SimpleBuilding is a lightweight building performance simulation engine designed to predict the heating and cooling energy in multi-zone commercial buildings. It was designed based on the simplified, dynamic hourly method from [ISO 13790:2009](http://www.iso.org/iso/catalogue_detail.htm?csnumber=41974) and from [Stephane Bertagnolio's PhD thesis](http://ulg.academia.edu/StephaneBertagnolio).

Users should be cautioned that this engine is not as comprehensive as a complex whole building simulation program like [EnergyPlus](http://apps1.eere.energy.gov/buildings/energyplus/), as it is meant to simulate system level metrics of common types of buildings with common systems. It is an attempt to apply the [Pareto Principle](http://en.wikipedia.org/wiki/Pareto_principle) to building performance simulation through a transparent, open-source, scalable simulation engine written in the high-level, easy-to-learn [Python](http://www.python.org/) scripting language.

##Getting Started
The Engine is designed to be utilized through an IPython console. The easiest way to get started is to download the free [Enthought Python Distribution (EPD)](http://www.enthought.com/products/epd_free.php) which includes all of the libraries required by the engine. 

After installation of EPD, the `pandas` library will need to be upgraded to at least [Version 0.90](http://pandas.pydata.org/).

From the IPython console, navigate to the source code folder and run:
    run -t SimpleBuildingEngine.py

##Beta Version
This engine is in the beta phase and much of the input/output functionality common to building simulation engines is very crude and code-integrated.