# Global River Widths from Landsat (GRWL) Database

Allen and Pavelsky (2018) paper: https://doi.org/10.1126/science.aat0636
Data: https://zenodo.org/record/1297434

Raster mask: https://zenodo.org/record/1297434/files/GRWL_mask_V01.01.zip?download=1
Vector mask: https://zenodo.org/record/1297434/files/GRWL_vector_V01.01.zip?download=1
Stats: https://zenodo.org/record/1297434/files/GRWL_summaryStats_V01.01.zip?download=1

Tiles in lonlat: https://drive.google.com/file/d/1K6x1E0mmLc0k7er4NCIeaZsTfHi2wxxI/view?usp=sharing

The raster mask was clipped using lonlat tiles but without line densification when reprojecting the tile
boundaries. This led to small gaps between the tiles (personal comm. by G. Allen). The unclipped/buffered
tiles were provided by G. Allen on Feb 3, 2022 and are in `RW_in`.

River mask pixel classifications: 
DN = 256 : No Data
DN = 255 : River
DN = 180 : Lake/reservoir 
DN = 126 : Tidal rivers/delta 
DN = 86 : Canal
DN = 0 : Land/water not connected to the GRWL river network