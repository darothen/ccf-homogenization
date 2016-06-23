# ccf-homogenization

This repository is an archive of a project conducted as part of Google's Summer of Code in 2011 under the supervision of the Climate Code Foundation. We aimed to re-implement the ["Pairwise Homogenization Algorithm"][PHA] ([Menne and Williams, 2009][MenneWilliams2009]) used to quality-control the [US Historical Climatology Network (version 2)][USHCN] by detecting and correcting undocumented shifts and breaks in station temperature and precipitation data.

The code herein *works*; that is, it will (for the most part) reproduce the original Fortran implementation of the PHA on small networks of real or synthetic station data. It's written entirely in Python 2.7, without leaning on any numerical libraries like NumPy, SciPy, or pandas - in fact, it was written before pandas was debuted to the broader scientific Python community!

I don't maintain this code and it's unlicensed - I'm not certain what license is suitable here, but if you are interested in using this code for any reason, please [contact me](mailto:darothen@mit.edu) and we can look into the appropriate details.

## Links for More Information

- Climate Code Foundation blog articles from GSoC 2011
  - [Welcome post](http://climatecode.org/blog/2011/05/welcome-daniel-rothenberg/)
  - [Homogenization Project Progress](http://climatecode.org/blog/2011/07/homogenization-project-progress/)
  - [Homogenization Report](http://climatecode.org/blog/2011/09/homogenization-report/)
- [GHCN-M v3.1.0 - Showing the value of engaging with software engineers](http://surfacetemperatures.blogspot.com/2011/11/ghcn-m-v310-showing-value-of-engaging.html)
- [Lessons from Deploying the USHCN Pairwise Homogenization Algorithm in Python](https://ams.confex.com/ams/92Annual/webprogram/Paper198219.html) - Talk at the 92nd AMS Annual Meeting (New Orleans, 2012)

[MenneWilliams2009]: ftp://ftp.ncdc.noaa.gov/pub/data/ushcn/v2/monthly/menne-williams2009.pdf
[PHA]: http://www.ncdc.noaa.gov/oa/climate/research/ushcn/#phas
[USHCN]: http://www.ncdc.noaa.gov/oa/climate/research/ushcn/
