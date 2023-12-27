# This simple script downloads content from:
# - https://caps.gsfc.nasa.gov/simpson/software/m33inv_f90.txt
# - https://caps.gsfc.nasa.gov/simpson/software/m44det_f90.txt
# - https://caps.gsfc.nasa.gov/simpson/software/m33det_f90.txt

import os

if not os.path.exists("./downloaded_f90_deps"):
    print("creating deps directory")
    os.system("mkdir downloaded_f90_deps")
    os.system("wget https://caps.gsfc.nasa.gov/simpson/software/m33det_f90.txt -O downloaded_f90_deps/m33det.f90")
    os.system("wget https://caps.gsfc.nasa.gov/simpson/software/m44det_f90.txt -O downloaded_f90_deps/m44det.f90")
    os.system("wget https://caps.gsfc.nasa.gov/simpson/software/m33inv_f90.txt -O downloaded_f90_deps/m33inv.f90")
    os.system("wget https://caps.gsfc.nasa.gov/simpson/software/m44inv_f90.txt -O downloaded_f90_deps/m44inv.f90")

# ================================================================================
# ================================================================================

def append_lines( base_file, input_file, keyword_start_read, keyword_end_read, complex_bool=False ):
    print("hello")

    fbase = open(base_file,"a")
    fread = open(input_file,"r")

    started_copying = False
    finished_copying = False

    while not finished_copying:
        line = fread.readline()
        # print("line: "+line)

        if keyword_start_read in line:
            started_copying = True

        if started_copying:
            if (complex_bool and ("M33INV" in line) ):
               print("M33INV!!")
               line = line.replace("M33INV","M33INVComp")

            if (complex_bool and ("M33DET" in line) ):
               print("M33DET!!")
               line = line.replace("M33DET","M33DETComp")

            if (complex_bool and ("DOUBLE PRECISION" in line) and (not "PARAMETER" in line) ):
               print("DOUBLE PRECISION")
               line = line.replace("DOUBLE PRECISION","DOUBLE complex")

            fbase.write(line)
            
        if keyword_end_read in line:
            finished_copying = True

# ================================================================================
# ================================================================================

fdeps = open("src/depsmodule.f90","w")
fdeps.write("module depsmodule\n  implicit none\n  contains\n")
fdeps.close()

append_lines("src/depsmodule.f90",
             "downloaded_f90_deps/m33inv.f90",
             "SUBROUTINE M33INV (A, AINV, OK_FLAG)",
             "END SUBROUTINE M33INV")

append_lines("src/depsmodule.f90",
             "downloaded_f90_deps/m44inv.f90",
             "SUBROUTINE M44INV (A, AINV, OK_FLAG)",               
             "END SUBROUTINE M44INV")

append_lines("src/depsmodule.f90",
             "downloaded_f90_deps/m33det.f90",
             "FUNCTION M33DET (A) RESULT (DET)",
             "END FUNCTION M33DET")

append_lines("src/depsmodule.f90",
             "downloaded_f90_deps/m44det.f90",
             "FUNCTION M44DET (A) RESULT (DET)",               
             "END FUNCTION M44DET")

append_lines("src/depsmodule.f90",
             "downloaded_f90_deps/m33inv.f90",
             "SUBROUTINE M33INV (A, AINV, OK_FLAG)",
             "END SUBROUTINE M33INV",
             True)

append_lines("src/depsmodule.f90",
             "downloaded_f90_deps/m33det.f90",
             "FUNCTION M33DET (A) RESULT (DET)",
             "END FUNCTION M33DET",
             True)

fdeps = open("src/depsmodule.f90","a")
fdeps.write("end module\n")
fdeps.close()
