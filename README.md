# roses2020

File Repository for ROSES 2020. You can view and run the course materials online by clicking on the badge below:

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/roseseismo/roses2020/HEAD)

ROSES = Remote Online Sessions for Emerging Seismologists. ROSES was conceived by the leadership of the Seismology Section of the American Geophysical Union, and ROSES 2020 was organized by Fransiska Dannemann Dugick and Suzan van der Lee. The school is targeted towards advanced Ph.D. students, who have used Python before and are familiar navigating in Linux/Unix. Lectures cover topics at the intersection of Seismology and Data Science.

Recorded lectures are hosted on the IRIS website [here](https://www.iris.edu/hq/inclass/course/roses), while this repository serves to maintain data exercises and lectures for the seismology community as a whole.

![ROSES logo](color_full.png)

## Summer School Topics and Instructors
6/25 (Th) ObsPy — Sydney Dybing, U. of Oregon

7/2 (Th) Data and Metadata — Emily Wolin, USGS

7/9 (Th) Time Series Analysis — German Prieto, U. Nacional de Colombia

7/14 (T) Waveform Cross Correlation — Elizabeth Berg, U. of Utah

7/21 (T) Array Seismology/Network Analysis — Stephen Arrowsmith, S. Methodist U.

7/28 (T) Polarization Analysis — Tolulope Olugboji, U. of Rochester

8/4 (T) Machine Learning — Zachary Ross, CalTech

8/11 (T) PyGMT — [Liam Toney](https://liam.earth), U. of Alaska Fairbanks

8/18 (T) Inversion, Bayesian — Steve Myers, L. Livermore National Lab

8/25 (T) Inversion, kriging — Tony Lowry, Utah State U.

9/1 (T) Inversion, gridding (and/or tomography) — Suzan van der Lee, Northwestern U.

## Authorship
Individual unit lectures and data exercises are produced by the respective instructor.  This GitHub repository is maintained by Fransiska Dannemann Dugick with assistance from the ROSES 2020 instructor team.

[![DOI](https://zenodo.org/badge/279129879.svg)](https://zenodo.org/badge/latestdoi/279129879)

## Installation
A `.yml` file is provided in order to create a Python environment using conda for running data exercises.  This environment can be created by running:
```
conda env create -f environment.yml
```

If this command executes correctly and finishes without errors, it should print out instructions on how to activate and deactivate the new environment. To activate the environment, use:
```
conda activate roses
```
To deactivate the environment, use:
```
conda deactivate
```

## Additional Resources
* [Python tutorial](https://docs.python.org/3/tutorial/index.html)
* [Getting Started with Anaconda](https://docs.anaconda.com/anaconda/user-guide/getting-started/)
* [Jupyter Notebook Overview](https://jupyter-notebook.readthedocs.io/en/stable/)
* [ObsPy Tutorial](https://docs.obspy.org/tutorial/)
