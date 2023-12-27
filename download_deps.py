# This simple script downloads content from https://caps.gsfc.nasa.gov/simpson/software/m33inv_f90.txt
# and https://caps.gsfc.nasa.gov/simpson/software/m33inv_f90.txt 

import os

if not os.path.exists("./downloaded_f90_deps"):
    print("creating deps directory")
    os.system("mkdir downloaded_f90_deps")

os.system("wget https://caps.gsfc.nasa.gov/simpson/software/m44det_f90.txt -O downloaded_f90_deps/m44det.f90")
os.system("wget https://caps.gsfc.nasa.gov/simpson/software/m33inv_f90.txt -O downloaded_f90_deps/m33inv.f90")

